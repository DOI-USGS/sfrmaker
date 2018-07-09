import time
import collections
import os
import numpy as np
import pandas as pd
import fiona
from fiona.crs import to_string
from shapely.ops import unary_union
from shapely.geometry import box
from .gis import shp2df, get_proj4, crs, read_polygon_feature, project, get_bbox
from . import grid
from .utils import make_graph, find_path, unit_conversion, \
    pick_toids, width_from_arbolate_sum, arbolate_sum, \
    consolidate_reach_conductances, renumber_segments
from .nhdplus_utils import load_NHDPlus_v2, get_prj_file
from .grid import grid as gridclass
from .sfrdata import sfrdata
from .checks import routing_is_circular

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
    length_units : 'meters' or 'feet'
        Length units for feature attributes (e.g. width, arbolate sum, etc.)
        (default 'meters')
    height_units : 'meters' or 'feet'
        Length units for elevation attributes
        (default 'meters')
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
                 length_units='meters', height_units='meters',
                 epsg=None, proj4=None, prjfile=None):

        self.df = df
        self.length_units = length_units
        self.height_units = height_units
        self.crs = crs(epsg=epsg, proj4=proj4, prjfile=prjfile)

        self._routing = None # dictionary of routing connections
        self._paths = None # routing sequence from each segment to outlet

    @property
    def to_m(self):
        if self.length_units == 'feet':
            return 0.3048
        return 1.0

    @property
    def to_ft(self):
        if self.length_units == 'feet':
            return 1.
        return 1/0.3048

    @property
    def routing(self):
        if self._routing is None or self._routing_changed():
            toid = self.df.toid.values
            # check whether or not routing is
            # many-to-one or one-to-one (no diversions)
            # squeeze it down
            to_one = False
            # if below == True, all toids are scalar or length 1 lists
            to_one = np.isscalar(np.squeeze(toid)[0])
            # if not, try converting any scalars to lists
            if not to_one:
                toid = [[l] if np.isscalar(l) else l for l in toid]
                to_one = np.isscalar(np.squeeze(toid)[0])
            toid = np.squeeze(toid)
            self._routing = make_graph(self.df.id.values, toid,
                                       one_to_many=not to_one)
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
        # compare the private routing attribute
        # to current values in reach data
        df_routing = dict(zip(self.df.id, self.df.toid))
        return df_routing != self._routing

    def cull(self, feature, simplify=False, tol=None,
             feature_crs=None, inplace=False):
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
        print('\nCulling hydrography to active area...')
        ta = time.time()
        df = self.df.copy()
        feature = read_polygon_feature(feature, self.crs,
                                       feature_crs=feature_crs)
        if simplify:
            print('simplification tolerance: {:.2f}'.format(tol))
            feature_s = feature.simplify(tol).buffer(tol)
        else:
            feature_s = feature
        lines = df.geometry.tolist()
        print('starting lines: {:,d}'.format(len(lines)))
        #isn = np.array([g.intersection(feature_s) for g in lines])
        #df['geometry'] = isn
        #drop = np.array([g.is_empty for g in isn])
        #df = df.loc[~drop]
        intersects = [g.intersects(feature_s) for g in lines]
        if not np.any(intersects):
            print('No lines in active area. Check CRS.')
            quit()

        df = df.loc[intersects]
        df['geometry'] = [g.intersection(feature) for g in df.geometry]
        drop = np.array([g.is_empty for g in df.geometry.tolist()])
        if len(drop) > 0:
            df = df.loc[~drop]
        print('remaining lines: {:,d}'.format(len(df)))
        if inplace:
            self.df = df
        else:
            return df
        print("finished in {:.2f}s\n".format(time.time() - ta))

    def intersect(self, grid, size_thresh=1e5):
        """Intersect linework with a model grid.

        Parameters
        ----------
        grid : instance of sfrmaker.grid
            Must have a valid Coordinate Reference System (CRS).
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

        # reproject the flowlines if they aren't in same CRS as grid
        if self.crs != grid.crs:
            self.reproject(grid.crs.proj4)

        grid_polygons = grid.df.geometry.tolist()
        stream_linework = self.df.geometry.tolist()
        id_list = self.df.id.tolist()

        ncells, nlines = len(grid_polygons), len(stream_linework)
        print("\nIntersecting {:,d} flowlines with {:,d} grid cells...".format(nlines, ncells))
        # building the spatial index takes a while
        # only use spatial index if number of tests exceeds size_thresh
        size = ncells * nlines
        # don't spend time on a spatial index if it isn't created and problem is small
        if size < size_thresh and grid._idx is None:
            grid_intersections = intersect(grid_polygons, stream_linework)
        else:
            idx = grid.spatial_index
            grid_intersections = intersect_rtree(grid_polygons, stream_linework, index=idx)

        # create preliminary reaches
        reach_data = setup_reach_data(stream_linework, id_list,
                                      grid_intersections, grid_polygons, tol=.001)

        column_order = ['node', 'k', 'i', 'j', 'rno',
                        'ireach', 'iseg', 'line_id', 'name', 'geometry']

        # transfer names
        names = dict(zip(self.df.id, self.df.name))
        reach_data['name'] = [names[line_id] for line_id in reach_data.line_id]

        # assign rows and columns from node
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

        print('\nreprojecting hydrography from\n{}\nto\n{}\n'.format(self.crs.proj4,
                                                                     dest_proj4))
        geoms = project(self.df.geometry, self.crs.proj4, dest_proj4)
        assert np.isfinite(np.max(geoms[0].xy[0])), \
            "Invalid reprojection; check CRS for lines and grid."
        self.df['geometry'] = geoms
        self.crs.proj4 = dest_proj4

    @staticmethod
    def from_shapefile(shapefile,
                       id_column='id',
                       routing_column='toid',
                       arbolate_sum_column='asum',
                       width_column='width',
                       up_elevation_column='elevup',
                       dn_elevation_column='elevdn',
                       name_column='name',
                       length_units='meters', height_units='meters',
                       filter=None,
                       epsg=None, proj4=None, prjfile=None):
        """
        Parameters
        ----------

        filter : tuple or str
            Bounding box (tuple) or shapefile of model stream network area.
        """

        if prjfile is None:
            prjfile = shapefile.replace('.shp', '.prj')
            prjfile = prjfile if os.path.exists(prjfile) else None

        shpfile_crs = crs(epsg=epsg, proj4=proj4, prjfile=prjfile)

        # ensure that filter bbox is in same crs as flowlines
        if filter is not None and not isinstance(filter, tuple):
            filter = get_bbox(filter, shpfile_crs)

        df = shp2df(shapefile, filter)
        assert 'geometry' in df.columns, "No feature geometries found in {}.".format(shapefile)

        return lines.from_dataframe(df,
                                    id_column=id_column,
                                    routing_column=routing_column,
                                    arbolate_sum_column=arbolate_sum_column,
                                    width_column=width_column,
                                    up_elevation_column=up_elevation_column,
                                    dn_elevation_column=dn_elevation_column,
                                    name_column=name_column,
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
                       name_column='name',
                       length_units='meters', height_units='meters',
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
                       dn_elevation_column: 'elevdn',
                       name_column: 'name'}

        rename_cols = {k:v for k, v in rename_cols.items() if k in df.columns}
        df.rename(columns=rename_cols, inplace=True)
        column_order = ['id', 'toid', 'asum', 'width', 'elevup', 'elevdn', 'name', 'geometry']
        for c in column_order:
            if c not in df.columns:
                df[c] = 0
        df = df[column_order].copy()

        return lines(df, length_units=length_units, height_units=height_units,
                     epsg=epsg, proj4=proj4, prjfile=prjfile)

    @staticmethod
    def from_NHDPlus_v2(NHDPlus_paths=None,
                        NHDFlowlines=None, PlusFlowlineVAA=None, PlusFlow=None, elevslope=None,
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
        df = load_NHDPlus_v2(NHDPlus_paths=NHDPlus_paths,
                             NHDFlowlines=NHDFlowlines, PlusFlowlineVAA=PlusFlowlineVAA,
                             PlusFlow=PlusFlow, elevslope=elevslope,
                             filter=filter,
                             epsg=epsg, proj4=proj4, prjfile=prjfile)

        if prjfile is None:
            prjfile = get_prj_file(NHDPlus_paths, NHDFlowlines)

        # convert arbolate sums from km to m
        df['asum'] = df.ArbolateSu * 1000

        # convert comid end elevations from cm to m
        if 'MAXELEVSMO' in df.columns:
            df['elevup'] = df.MAXELEVSMO / 100.
        if 'MINELEVSMO' in df.columns:
            df['elevdn'] = df.MINELEVSMO / 100.

        return lines.from_dataframe(df, id_column='COMID',
                                    routing_column='tocomid',
                                    name_column='GNIS_NAME',
                                    length_units='meters', height_units='meters',
                                    epsg=epsg, proj4=proj4, prjfile=prjfile)

    def to_sfr(self, grid=None, sr=None,
               active_area=None, isfr=None,
               minimum_reach_length=None,
               consolidate_conductance=False, one_reach_per_cell=False,
               model_name=None,
               **kwargs):
        """Create a streamflow routing dataset from the information
        in sfrmaker.lines class instance and a supplied sfrmaker.grid class instance.

        Parameters
        ----------
        grid : sfrmaker.grid instance
            Describes the numerical model grid. Grid must have
        minimum_reach_length : float
            Minimum reach length to retain. Default is to compute
            an effective mean model cell length by taking the square root
            of the average cell area, and then set minimum_reach_length
            to 5% of effective mean cell length.
        consolidate_conductance : bool
            If True, total reach conductance each cell is computed, and
            assigned to the most downstream reach via the hydraulic conductivity
            parameter.
        one_reach_per_cell : bool
            If True, streambed conductance in each reach is consolidated
            (consolidate_conductance = True), and additional reaches besides
            the most downstream reach are dropped.
        model_name : str
            Base name for writing sfr output.
        kwargs : keyword arguments to sfrmaker.sfrdata

        Returns
        -------
        sfrdata : sfrmaker.sfrdata instance
        """
        print("\nCreating sfr dataset...")
        totim = time.time()

        if grid is None and sr is not None:
            print('\nCreating grid class instance from flopy SpatialReference...')
            ta = time.time()
            grid = gridclass.from_sr(sr, active_area=active_area, isfr=isfr)
            print("grid class created in {:.2f}s\n".format(time.time() - ta))

        # print grid information to screen
        print(grid)

        # reproject the flowlines if they aren't in same CRS as grid
        if self.crs != grid.crs:
            self.reproject(grid.crs.proj4)
        if grid.active_area is not None:
            self.cull(grid.active_area, inplace=True, simplify=True, tol=2000)
        elif grid._bounds is not None: # cull to grid bounding box if already computed
            self.cull(box(*grid._bounds), inplace=True)
        if model_name is None:
            model_name = 'model'

        mult = unit_conversion.get(self.length_units+grid.model_units, 1.)
        mult_h = unit_conversion.get(self.height_units+grid.model_units, 1.)

        # convert routing connections (toid column) from lists (one-to-many)
        # to ints (one-to-one or many-to-one)
        elevup = dict(zip(self.df.id, self.df.elevup))
        routing = self.routing.copy()
        to_one = np.isscalar(np.squeeze(list(routing.values()))[0])
        if not to_one:
            routing = pick_toids(routing, elevup)
        valid_ids = routing.keys()
        # df.toid column is basis for routing attributes
        # all paths terminating in invalid toids (outside of the model)
        # will be none; set invalid toids = 0
        self.df.toid = [routing[i] if routing[i] in valid_ids else 0
                        for i in self.df.id.tolist()]

        # intersect lines with model grid to get preliminary reaches
        rd = self.intersect(grid)

        print("Computing widths...")
        # length in linework units (not model units)
        rd['rchlen'] = np.array([g.length for g in rd.geometry])
        lengths = rd[['line_id', 'rchlen']].copy() # length of each intersected line fragment
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
        reach_cumsums = np.concatenate([np.cumsum(grp.rchlen.values[::-1])[::-1] - 0.5*grp.rchlen.values
                                      for s, grp in groups])
        reach_asums = -1 * reach_cumsums + segment_asums
        # maintain positive asums; lengths in NHD often aren't exactly equal to feature lengths
        reach_asums[reach_asums < 0.] = 0
        rd['asum'] = reach_asums
        # width estimation formula expects meters
        width = width_from_arbolate_sum(reach_asums * self.to_m)
        width = width / self.to_m # convert back to original units

        # convert length and width to model units
        rd['rchlen'] *= mult
        rd['width'] = width * mult

        # discard very small reaches; redo numbering
        # set minimum reach length based on cell size
        thresh = 0.05 # fraction of cell length (based on square root of area)
        if minimum_reach_length is None:
            cellgeoms = grid.df.loc[rd.node.values, 'geometry']
            mean_area = np.mean([g.area for g in cellgeoms])
            minimum_reach_length = np.sqrt(mean_area) * thresh

        inds = rd.rchlen > minimum_reach_length
        print('\nDropping {} reaches with length < {:.2f} {}...'.format(np.sum(~inds),
                                                                        minimum_reach_length,
                                                                        self.length_units))
        rd = rd.loc[inds].copy()

        # handle co-located reaches
        if consolidate_conductance or one_reach_per_cell:
            rd = consolidate_reach_conductances(rd, keep_only_dominant=one_reach_per_cell)

        # patch the routing
        # 1) reduce one to many routing to one-to-one routing (pick_toids() above)
        # 2) create new graph with just one-to-one segments
        # 3) list new paths;  2) and 3) should be automatic following 1)
        # 4) code below will update new graph to only include remaining segments
        # 5) create sfrdata instance; numbering will be messed up
        # 6) run methods on sfrdata instance to fix numbering and route reaches with unique numbers

        print('\nRepairing routing connections...')
        remaining_ids = rd.line_id.unique()
        # routing and paths properties should update automatically
        # when id and toid columns are changed in self.df
        # but only rd (reach_data) has been changed
        new_routing = {}
        # for each segment
        for k in remaining_ids:
            # interate through successive downstream segments
            for s in self.paths[k][1:]:
                # assign the first segment that still exists as the outseg
                if s in remaining_ids:
                    new_routing[k] = s
                    break
            # if no segments are left downstream, assign outlet
            if k not in new_routing.keys():
                new_routing[k] = 0

        # map remaining_ids to segment numbers
        segment = dict(zip(rd.line_id, rd.iseg))
        line_id = {s: lid for lid, s in segment.items()}

        # get the segment associated with each line id
        nseg = [segment[rid] for rid in remaining_ids]
        # get the segment associated with new connection for each line id
        outseg = [segment.get(new_routing[line_id[s]], 0) for s in nseg]

        # renumber the segments to be consecutive,
        # starting at 1 and only increasing downstream
        r = renumber_segments(nseg, outseg)
        # map new segment numbers to line_ids
        line_id = {r[s]: lid for s, lid in line_id.items()}
        #segment2 = {lid: r[s] for lid, s in segment.items()}
        #line_id2 = {s: lid for lid, s in segment2.items()}

        # update reach_data
        rd['iseg'] = [r[s] for s in rd.iseg]

        print('\nSetting up segment data...')
        sd = pd.DataFrame()
        #sd['nseg'] = [segment2[lid] for lid in remaining_ids]
        #sd['outseg'] = [segment2.get(new_routing[line_id2[s]], 0) for s in nseg]
        sd['nseg'] = [r[s] for s in nseg]
        sd['outseg'] = [r[s] for s in outseg]
        sd.sort_values(by='nseg', inplace=True)

        # verify that no segments route to themselves
        assert not routing_is_circular(sd.nseg, sd.outseg)

        # (elevup dict was created above)
        elevdn = dict(zip(self.df.id, self.df.elevdn))
        sd['elevup'] = [elevup[line_id[s]] for s in sd.nseg]
        sd['elevdn'] = [elevdn[line_id[s]] for s in sd.nseg]
        # convert elevation units
        sd['elevup'] *= mult_h
        sd['elevdn'] *= mult_h

        # create sfrdata instance
        # this class has methods for fix segment and reach numbering,
        # assigning elevations and other properties by reach,
        # smoothing elevations, writing sfr package files
        # and other output
        rd = rd[[c for c in sfrdata.rdcols if c in rd.columns]].copy()
        sfrd = sfrdata(reach_data=rd, segment_data=sd, grid=grid,
                       model_name=model_name, **kwargs)
        print("\nTime to create sfr dataset: {:.2f}s\n".format(time.time() - totim))
        return sfrd
