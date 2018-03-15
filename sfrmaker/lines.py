import time
import collections
import os
import numpy as np
import pandas as pd
import fiona
from fiona.crs import to_string
from shapely.ops import unary_union
from shapely.geometry import box
from .gis import shp2df, get_proj4, crs, read_feature, project
from . import grid
from .utils import make_graph, find_path, unit_conversion, \
    width_from_arbolate_sum, arbolate_sum, consolidate_reach_conductances
from .nhdplus_utils import load_NHDPlus_v2
from .sfrdata import sfrdata

class lines:
    """Class for working with linestring feature input.

    Parameters
    ----------
    df : DataFrame
        Dataframe with linestring features and attribute information.
        Must have the following columns:
        id: sequential integers for identifying each feature
        toid: integers representing routing connections
        ...
        geometry: shapely LineString objects for each feature
    length_units : 'm' or 'ft'
        Length units for feature attributes (e.g. width, arbolate sum, etc.)
        (default 'm')
    height_units : 'm' or 'ft'
        Length units for elevation attributes
        (default 'm')
    epsg: int
        EPSG code identifying Coordinate Reference System (CRS)
        for features in df.geometry
        (optional)
    proj4: str
        PROJ4 string identifying CRS for features in df.geometry
        (optional)
    prjfile: str
        File path to projection (.prj) file identifying CRS
        for features in df.geometry
        (optional)
    """

    def __init__(self, df=None,
                 length_units='m', height_units='m',
                 epsg=None, proj4=None, prjfile=None):

        self.df = df
        self.length_units = length_units
        self.height_units = height_units
        self.crs = crs(epsg=epsg, proj4=proj4, prjfile=prjfile)

        self._routing = None # dictionary of routing connections
        self._paths = None # routing sequence from each segment to outlet

    @property
    def to_m(self):
        if self.length_units == 'ft':
            return 0.3048
        return 1.0

    @property
    def to_ft(self):
        if self.length_units == 'ft':
            return 1.
        return 1/0.3048

    @property
    def routing(self):
        if self._routing is None or self._routing_changed():
            self._routing = make_graph(self.df.id.values, self.df.toid.values)
        return self._routing

    @property
    def paths(self):
        if self._paths is None:
            self._set_paths()
            return self._paths
        if self._routing_changed():
            self._routing = None
            self._set_paths()
        return self._paths

    def _set_paths(self):
        routing = self.routing
        self._paths = {seg: find_path(routing, seg) for seg in routing.keys()}

    def _routing_changed(self):
        # check to see if routing in segment data was changed
        ids = np.array(sorted(self._paths.keys()), dtype=int)
        ids = ids[ids > 0].copy()
        toids = np.array([self._paths[k][1] for k in ids])

        return not np.array_equal(ids, self.df.id.values) or \
               not np.array_equal(toids, self.df.toid.values)

    def cull(self, feature, inplace=False):
        """Cull linework; retaining only the
        lines that intersect a polygon feature.

        Parameters
        ----------
        feature : shapely Polygon, list of Polygons, or shapefile path
            Polygons must be in same CRS as linework; shapefile
            features will be reprojected if their crs is different.
        inplace : bool
            If True, the attribute .df is modified;
            if False, a copy of .df is returned.
        """
        feature = read_feature(feature, self.crs)

        intersects = np.array([g.intersects(feature) for g in self.df.geometry])
        if inplace:
            self.df = self.df.loc[intersects]
        else:
            return self.df.loc[intersects].copy()

    def intersect(self, grid, size_thresh=1e5):
        """Intersect linework with a model grid.

        Parameters
        ----------
        grid : instance of sfrmaker.grid
        size_thresh : int
            Determines whether spatial indexing will be used. If
            number of grid cells x number of flowlines > size_thresh,
            a spatial index (rtree package) will be used to speed
            up the intersections.

        Returns
        -------
        reach_data : DataFrame
            DataFrame containing intersected reaches with grid cell information
            and original linework IDs.
        """
        from .gis import intersect, intersect_rtree
        from .utils import setup_reach_data

        # reproject the flowlines in they aren't in same CRS as grid
        if self.crs != grid.crs:
            self.reproject(grid.crs.proj4)

        grid_polygons = grid.df.geometry.tolist()
        stream_linework = self.df.geometry.tolist()
        id_list = self.df.id.tolist()

        print("intersecting flowlines with grid cells...")
        # building the spatial index takes a while
        # only use spatial index if number of tests exceeds size_thresh
        size = len(grid_polygons) * len(stream_linework)
        if size < size_thresh:
            grid_intersections = intersect(grid_polygons, stream_linework)
        else:
            grid_intersections = intersect_rtree(grid_polygons, stream_linework)

        print("setting up reach data... (may take a few minutes for large grids)")
        reach_data = setup_reach_data(stream_linework, id_list,
                                      grid_intersections, grid_polygons, tol=.001)

        column_order = ['node', 'k', 'i', 'j', 'reachID',
                        'reach', 'segment', 'line_id', 'geometry']
        if grid.structured:
            reach_data['k'] = 0
            reach_data['i'] = np.floor(reach_data.node / grid.ncol).astype(int)
            reach_data['j'] = reach_data.node.values % grid.ncol
        else:
            column_order.remove('k')
            column_order.remove('i')
            column_order.remove('j')
        return reach_data[column_order].copy()

    def reproject(self, dest_proj4):
        assert self.crs.proj4 is not None, "No proj4 string for flowlines"
        assert dest_proj4 is not None, "No destination CRS."

        geoms = project(self.df.geometry, self.crs.proj4, dest_proj4)
        assert np.isfinite(np.max(geoms[0].xy[0])), \
            "Invalid reprojection; check CRS for lines and grid."
        self.df['geometry'] = geoms

    @staticmethod
    def from_shapefile(shapefile,
                       id_column='id',
                       routing_column='toid',
                       arbolate_sum_column='asum',
                       width_column='width',
                       up_elevation_column='elevup',
                       dn_elevation_column='elevdn',
                       length_units='m', height_units='m',
                       epsg=None, proj4=None, prjfile=None):

        if prjfile is None:
            prjfile = shapefile.replace('.shp', '.prj')
            prjfile = prjfile if os.path.exists(prjfile) else None

        df = shp2df(shapefile)
        assert 'geometry' in df.columns, "No feature geometries found in {}.".format(shapefile)

        return lines.from_dataframe(df,
                                    id_column=id_column,
                                    routing_column=routing_column,
                                    arbolate_sum_column=arbolate_sum_column,
                                    width_column=width_column,
                                    up_elevation_column=up_elevation_column,
                                    dn_elevation_column=dn_elevation_column,
                                    length_units=length_units, height_units=height_units,
                                    epsg=epsg, proj4=proj4, prjfile=prjfile)

    @staticmethod
    def from_dataframe(df,
                       id_column='id',
                       routing_column='toid',
                       arbolate_sum_column='asum',
                       width_column='width',
                       up_elevation_column='elevup',
                       dn_elevation_column='elevdn',
                       geometry_column='geometry',
                       length_units='m', height_units='m',
                       epsg=None, proj4=None, prjfile=None):

        assert geometry_column in df.columns, \
            "No feature geometries found: dataframe column '{}' doesn't exist.".format(geometry_column)
        assert routing_column in df.columns, \
            "No routing information found; dataframe column '{}' doesn't exist.".format(routing_column)

        # rename the columns for consistency
        rename_cols = {id_column: 'id',
                       routing_column: 'toid',
                       arbolate_sum_column: 'asum',
                       width_column: 'width',
                       up_elevation_column: 'elevup',
                       dn_elevation_column: 'elevdn'}

        rename_cols = {k:v for k, v in rename_cols.items() if k in df.columns}
        df.rename(columns=rename_cols, inplace=True)
        column_order = ['id', 'toid', 'asum', 'width', 'elevup', 'elevdn', 'geometry']
        for c in column_order:
            if c not in df.columns:
                df[c] = 0
        df = df[column_order].copy()

        return lines(df, length_units=length_units, height_units=height_units,
                     epsg=epsg, proj4=proj4, prjfile=prjfile)

    @staticmethod
    def from_NHDPlus_v2(NHDFlowlines=None, PlusFlowlineVAA=None, PlusFlow=None, elevslope=None,
             filter=None,
             epsg=None, proj4=None, prjfile=None):
        """
        Parameters
        ==========
        NHDFlowlines : str or list of strings.
            Shapefile or list of NHDFlowline shapefiles containing
            feature geometries (line arcs) for stream network. Must contain
            the following attribute fields:
            COMID : common identifier number
        PlusFlowlineVAA : str or list of strings.
            DBF file or list of DBF files with NHDPlus attribute information.
            Must contain the following attribute fields:
            COMID : common identifier number
        PlusFlow : str or list of strings.
            DBF file or list of DBF files with NHDPlus routing information.
            Must contain the following attribute fields:
            COMID : common identifier number
        elevslope : str or list of strings.
            DBF file or list of DBF files with end elevations for each
            line arc in NHDFlowlines. Must contain the following attribute fields:
            COMID : common identifier number
        filter : tuple or str
            Bounding box (tuple) or shapefile of model stream network area.
        """
        df = load_NHDPlus_v2(NHDFlowlines=NHDFlowlines, PlusFlowlineVAA=PlusFlowlineVAA,
                             PlusFlow=PlusFlow, elevslope=elevslope,
                             filter=filter,
                             epsg=epsg, proj4=proj4, prjfile=prjfile)

        if prjfile is None:
            prjfile = NHDFlowlines if not isinstance(NHDFlowlines, list) else NHDFlowlines[0]

        # convert arbolate sums from km to m
        df['asum'] = df.ArbolateSu * 1000

        # convert comid end elevations from cm to m
        for c in ['MAXELEVSMO', 'MINELEVSMO']:
            if c in df.columns:
                df[c] = df[c] / 100.
        if 'MAXELEVSMO' in df.columns:
            df['elevup'] = df.MAXELEVSMO / 100.
        if 'MINELEVSMO' in df.columns:
            df['elevdn'] = df.MINELEVSMO / 100.

        return lines.from_dataframe(df, id_column='COMID',
                                    routing_column='tocomid',
                                    length_units='m', height_units='m',
                                    epsg=epsg, proj4=proj4, prjfile=prjfile)

    def to_sfr(self, grid,
               minimum_reach_length=None,
               consolidate_conductance=False, one_reach_per_cell=False):

        if grid.active_area is not None:
            self.cull(grid.active_area)
        elif grid._bounds is not None: # cull to grid bounding box if already computed
            self.cull(box(*grid._bounds))

        mult = unit_conversion[self.length_units+grid.model_units]

        # intersect lines with model grid to get preliminary reaches
        rd = self.intersect(grid)

        print("computing widths...")
        # length in linework units (not model units)
        rd['length'] = np.array([g.length for g in rd.geometry])
        lengths = rd[['line_id', 'length']].copy() # length of each intersected line fragment
        groups = lengths.groupby('line_id') # fragments grouped by parent line

        # compute arbolate sums for lines if they weren't provided
        if 'asum' not in self.df.columns:
            asums = arbolate_sum(self.df.id,
                                 dict(zip(self.df.id, [g.length for g in self.df.geometry])),
                                 self.routing)
        else:
            asums = dict(zip(self.df.id, self.df.asum))

        # compute arbolate sum at reach midpoints
        segment_asums = [asums[id] for id in lengths.line_id]
        reach_cumsums = np.concatenate([np.cumsum(grp.length.values[::-1])[::-1] - 0.5*grp.length.values
                                      for s, grp in groups])
        reach_asums = -1 * reach_cumsums + segment_asums
        # maintain positive asums; lengths in NHD often aren't exactly equal to feature lengths
        reach_asums[reach_asums < 0.] = 0
        rd['asum'] = reach_asums
        # width estimation formula expects meters
        width = width_from_arbolate_sum(reach_asums * self.to_m)
        width = width / self.to_m # convert back to original units

        # convert length and width to model units
        rd['length'] *= mult
        rd['widht'] *= mult

        # discard very small reaches; redo numbering
        # set minimum reach length based on cell size
        thresh = 0.05 # fraction of cell length (based on square root of area)
        if minimum_reach_length is None:
            cellgeoms = grid.df.node.loc[rd.node.values, 'geometry']
            mean_area = np.mean([g.area for g in cellgeoms])
            minimum_reach_length = np.sqrt(mean_area) * thresh

        inds = rd.length > minimum_reach_length
        print('Dropping {} reaches with length < {:.2f} ...'.format(np.sum(~inds), minimum_reach_length))
        rd = rd.loc[inds].copy()

        # handle co-located reaches
        if consolidate_conductance or one_reach_per_cell:
            rd = consolidate_reach_conductances(rd, keep_only_dominant=one_reach_per_cell)

        # patch the routing
        # 1) reduce one to many routing to one-to-one routing
        # 2) create new graph with just one-to-one segments
        # 3) list new paths
        # 4) code below will update new graph to only include remaining segments
        # 5) create sfrdata instance; numbering will be messed up
        # 6) run methods on sfrdata instance to fix numbering and route reaches with unique numbers
        print('\nRepairing routing connections...')
        remaining_segments = rd.segment.unique()

        old_graph = self.routing
        paths = self.paths
        new_graph = {}
        for k in remaining_segments:
            for s in paths[k][1:]:
                if s in remaining_segments:
                    new_graph[k] = s
                    break
                else:
                    new_graph[k] = 0


        # set up segment data
        sd = pd.DataFrame()




        sfrd = sfrdata(reach_data=rd, segment_data=sd)
