import os
import time

import numpy as np
import pandas as pd
from shapely.geometry import box
import flopy
from gisutils import shp2df, df2shp, project
import sfrmaker
from sfrmaker.routing import pick_toids, find_path, make_graph, renumber_segments
from .checks import routing_is_circular, is_to_one
from .gis import crs, read_polygon_feature, get_bbox
from .grid import StructuredGrid
from .nhdplus_utils import load_nhdplus_v2, get_prj_file
from .sfrdata import SFRData
from .units import convert_length_units, get_length_units
from .utils import (width_from_arbolate_sum, arbolate_sum, \
                    consolidate_reach_conductances,
                    interpolate_to_reaches)


class Lines:
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
    attr_length_units : 'meters' or 'feet'
        Length units for feature attributes (e.g. width, arbolate sum, etc.)
        (default 'meters')
    attr_height_units : 'meters' or 'feet'
        Length units for elevation attributes
        (default 'meters')
    epsg: int
        EPSG code identifying Coordinate Reference System (CRS)
        for features in df.geometry
        (optional)
    proj_str: str
        proj_str string identifying CRS for features in df.geometry
        (optional)
    prjfile: str
        File path to projection (.prj) file identifying CRS
        for features in df.geometry
        (optional)
    """

    def __init__(self, df=None,
                 attr_length_units='meters',
                 attr_height_units='meters',
                 epsg=None, proj_str=None, prjfile=None):

        self.df = df
        self.attr_length_units = attr_length_units
        self.attr_height_units = attr_height_units
        self.crs = crs(epsg=epsg, proj_str=proj_str, prjfile=prjfile)
        self._geometry_length_units = None

        self._routing = None  # dictionary of routing connections
        self._paths = None  # routing sequence from each segment to outlet

    @property
    def attr_to_m(self):
        if self.attr_length_units == 'feet':
            return 0.3048
        return 1.0

    @property
    def attr_to_ft(self):
        if self.attr_length_units == 'feet':
            return 1.
        return 1 / 0.3048

    @property
    def geometry_to_m(self):
        if self.geometry_length_units == 'feet':
            return 0.3048
        return 1.0

    @property
    def geometry_to_ft(self):
        if self.geometry_length_units == 'feet':
            return 1.
        return 1 / 0.3048

    @property
    def geometry_length_units(self):
        self._geometry_length_units = self.crs.length_units
        if self._geometry_length_units is None:
            print("Warning: No length units specified in CRS for input LineStrings; "
                  "defaulting to meters.")
            self._geometry_length_units = 'meters'
        return self._geometry_length_units

    @property
    def routing(self):
        if self._routing is None or self._routing_changed():
            toid = self.df.toid.values
            # check whether or not routing is
            # many-to-one or one-to-one (no diversions)
            # squeeze it down
            to_one = False
            # if below == True, all toids are scalar or length 1 lists
            if len(toid) > 1:
                to_one = is_to_one(toid)
                # if not, try converting any scalars to lists
                if not to_one:
                    toid = [[l] if np.isscalar(l) else l for l in toid]
                    to_one = is_to_one(toid)
                toid = np.squeeze(list(toid))
                routing = make_graph(self.df.id.values, toid,
                                     one_to_many=not to_one)
                if not to_one:
                    routing = pick_toids(routing, self.elevup)
            else:
                routing = {self.df.id.values[0]: 0}
            self._routing = routing
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
            feature_s = feature.simplify(tol).buffer(tol).buffer(0)
        else:
            feature_s = feature.buffer(0)  # in case feature is invalid, might fix

        lines = df.geometry.tolist()
        print('starting lines: {:,d}'.format(len(lines)))
        # isn = np.array([g.intersection(feature_s) for g in lines])
        # df['geometry'] = isn
        # drop = np.array([g.is_empty for g in isn])
        # df = df.loc[~drop]
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
            self.reproject(grid.crs.proj_str)

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
        if isinstance(grid, StructuredGrid):
            reach_data['k'] = 0
            reach_data['i'] = np.floor(reach_data.node / grid.ncol).astype(int)
            reach_data['j'] = reach_data.node.values % grid.ncol
        else:
            column_order.remove('k')
            column_order.remove('i')
            column_order.remove('j')
        return reach_data[column_order].copy()

    def reproject(self, dest_proj_str):
        assert self.crs.proj_str is not None, "No proj_str string for flowlines"
        assert dest_proj_str is not None, "No destination CRS."

        print('\nreprojecting hydrography from\n{}\nto\n{}\n'.format(self.crs.proj_str,
                                                                     dest_proj_str))
        geoms = project(self.df.geometry, self.crs.proj_str, dest_proj_str)
        assert np.isfinite(np.max(geoms[0].xy[0])), \
            "Invalid reprojection; check CRS for lines and grid."
        self.df['geometry'] = geoms
        self.crs.proj_str = dest_proj_str

    def write_shapefile(self, outshp='flowlines.shp'):
        df2shp(self.df, outshp, epsg=self.crs.epsg, prj=self.crs.prjfile)

    @classmethod
    def from_shapefile(cls, shapefile,
                       id_column='id',
                       routing_column='toid',
                       arbolate_sum_column2='asum2',
                       width1_column='width1',
                       width2_column='width2',
                       up_elevation_column='elevup',
                       dn_elevation_column='elevdn',
                       name_column='name',
                       attr_length_units='meters', attr_height_units='meters',
                       filter=None,
                       epsg=None, proj_str=None, prjfile=None):
        """
        Parameters
        ----------

        filter : tuple or str
            Bounding box (tuple) or shapefile of model stream network area.
        """

        if prjfile is None:
            prjfile = shapefile.replace('.shp', '.prj')
            prjfile = prjfile if os.path.exists(prjfile) else None

        shpfile_crs = crs(epsg=epsg, proj_str=proj_str, prjfile=prjfile)

        # ensure that filter bbox is in same crs as flowlines
        if filter is not None and not isinstance(filter, tuple):
            filter = get_bbox(filter, shpfile_crs)

        df = shp2df(shapefile, filter=filter)
        assert 'geometry' in df.columns, "No feature geometries found in {}.".format(shapefile)

        return cls.from_dataframe(df,
                                  id_column=id_column,
                                  routing_column=routing_column,
                                  arbolate_sum_column2=arbolate_sum_column2,
                                  width1_column=width1_column,
                                  width2_column=width2_column,
                                  up_elevation_column=up_elevation_column,
                                  dn_elevation_column=dn_elevation_column,
                                  name_column=name_column,
                                  attr_length_units=attr_length_units,
                                  attr_height_units=attr_height_units,
                                  epsg=epsg, proj_str=proj_str, prjfile=prjfile)

    @classmethod
    def from_dataframe(cls, df,
                       id_column='id',
                       routing_column='toid',
                       arbolate_sum_column2='asum2',
                       width1_column='width1',
                       width2_column='width2',
                       up_elevation_column='elevup',
                       dn_elevation_column='elevdn',
                       geometry_column='geometry',
                       name_column='name',
                       attr_length_units='meters',
                       attr_height_units='meters',
                       epsg=None, proj_str=None, prjfile=None):

        assert geometry_column in df.columns, \
            "No feature geometries found: dataframe column '{}' doesn't exist.".format(geometry_column)
        assert routing_column in df.columns, \
            "No routing information found; dataframe column '{}' doesn't exist.".format(routing_column)

        # rename the columns for consistency
        rename_cols = {id_column: 'id',
                       routing_column: 'toid',
                       arbolate_sum_column2: 'asum2',
                       width1_column: 'width1',
                       width2_column: 'width2',
                       up_elevation_column: 'elevup',
                       dn_elevation_column: 'elevdn',
                       name_column: 'name'}

        # dictionary for assigning new column names
        rename_cols = {k: v for k, v in rename_cols.items() if k in df.columns and v != k}
        # drop any existing columns that have one of the new names
        # (otherwise pandas will create a DataFrame
        # instead of a Series under that column name)
        to_drop = set(rename_cols.values()).intersection(df.columns)
        df.drop(to_drop, axis=1, inplace=True)
        df.rename(columns=rename_cols, inplace=True)

        column_order = ['id', 'toid',
                        'asum1', 'asum2',
                        'width1', 'width2',
                        'elevup', 'elevdn',
                        'name', 'geometry']
        for c in column_order:
            if c not in df.columns:
                df[c] = 0
            else:
                assert isinstance(df[c], pd.Series)
        df = df[column_order].copy()

        return cls(df, attr_length_units=attr_length_units,
                   attr_height_units=attr_height_units,
                   epsg=epsg, proj_str=proj_str, prjfile=prjfile)


    @classmethod
    def from_nhdplus_v2(cls, NHDPlus_paths=None,
                        NHDFlowlines=None, PlusFlowlineVAA=None, PlusFlow=None, elevslope=None,
                        filter=None,
                        epsg=None, proj_str=None, prjfile=None):
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
        df = load_nhdplus_v2(NHDPlus_paths=NHDPlus_paths,
                             NHDFlowlines=NHDFlowlines, PlusFlowlineVAA=PlusFlowlineVAA,
                             PlusFlow=PlusFlow, elevslope=elevslope,
                             filter=filter,
                             epsg=epsg, proj_str=proj_str, prjfile=prjfile)

        if prjfile is None:
            prjfile = get_prj_file(NHDPlus_paths, NHDFlowlines)

        # convert arbolate sums from km to m
        df['asum2'] = df.ArbolateSu * 1000

        # convert comid end elevations from cm to m
        if 'MAXELEVSMO' in df.columns:
            df['elevup'] = df.MAXELEVSMO / 100.
        if 'MINELEVSMO' in df.columns:
            df['elevdn'] = df.MINELEVSMO / 100.

        return cls.from_dataframe(df, id_column='COMID',
                                  routing_column='tocomid',
                                  name_column='GNIS_NAME',
                                  attr_length_units='meters',
                                  attr_height_units='meters',
                                  epsg=epsg, proj_str=proj_str, prjfile=prjfile)


    def to_sfr(self, grid=None, sr=None,
               active_area=None, isfr=None,
               model=None,
               model_length_units='undefined',
               minimum_reach_length=None,
               minimum_reach_width=1.,
               cull_flowlines_to_active_area=True,
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
                minimum_reach_length : float
        minimum_reach_width : float
            Minimum reach width to specify (in model units), if computing widths from
            arbolate sum values. (default = 1)
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
        print("\nSFRmaker version {}".format(sfrmaker.__version__))
        print("\nCreating sfr dataset...")
        totim = time.time()

        if grid is None and sr is not None:
            print('\nCreating grid class instance from flopy SpatialReference...')
            ta = time.time()
            grid = StructuredGrid.from_sr(sr, active_area=active_area, isfr=isfr)
            print("grid class created in {:.2f}s\n".format(time.time() - ta))
        elif isinstance(grid, flopy.discretization.StructuredGrid):
            print('\nCreating grid class instance from flopy Grid instance...')
            ta = time.time()
            grid = StructuredGrid.from_modelgrid(grid, active_area=active_area, isfr=isfr)
            print("grid class created in {:.2f}s\n".format(time.time() - ta))
        elif not isinstance(grid, sfrmaker.grid.Grid):
            raise TypeError('Unrecognized input for grid: {}'.format(grid))

        # print grid information to screen
        print(grid)

        # print model information to screen
        print(model)

        model_length_units = get_length_units(model_length_units, grid, model)
        mult = convert_length_units(self.attr_length_units, model_length_units)
        mult_h = convert_length_units(self.attr_height_units, model_length_units)
        gis_mult = convert_length_units(self.geometry_length_units, model_length_units)

        # reproject the flowlines if they aren't in same CRS as grid
        if self.crs != grid.crs:
            self.reproject(grid.crs.proj_str)
        if cull_flowlines_to_active_area:
            if grid.active_area is not None:
                self.cull(grid.active_area, inplace=True, simplify=True, tol=2000)
            elif grid._bounds is not None:  # cull to grid bounding box if already computed
                self.cull(box(*grid._bounds), inplace=True)
        if model_name is None:
            model_name = 'model'

        # convert routing connections (toid column) from lists (one-to-many)
        # to ints (one-to-one or many-to-one)
        self.elevup = dict(zip(self.df.id, self.df.elevup))
        routing = self.routing.copy()

        # one to many routing is not supported
        to_one = is_to_one(routing.values())
        assert to_one, "routing is still one-to-many"
        # if not to_one:
        #    routing = pick_toids(routing, elevup)
        valid_ids = routing.keys()
        # df.toid column is basis for routing attributes
        # all paths terminating in invalid toids (outside of the model)
        # will be none; set invalid toids = 0
        # TODO: write a test for pick_toids if some IDs route to more than one connection
        assert not np.any([isinstance(r, list) for r in routing.items()]), "one to many routing not supported"
        self.df.toid = [routing[i] if routing[i] in valid_ids else 0
                        for i in self.df.id.tolist()]

        # intersect lines with model grid to get preliminary reaches
        rd = self.intersect(grid)

        # length of intersected line fragments (in model units)
        rd['rchlen'] = np.array([g.length for g in rd.geometry]) * gis_mult

        # estimate widths if they aren't supplied
        if self.df.width1.sum() == 0:
            print("Computing widths...")

            # compute arbolate sums for original LineStrings if they weren't provided
            if 'asum2' not in self.df.columns:
                raise NotImplementedError('Check length unit conversions before using this option.')
                asums = arbolate_sum(self.df.id,
                                     dict(zip(self.df.id,
                                              np.array([g.length for g in self.df.geometry]) * self.geometry_to_m
                                              )),
                                     self.routing)
            else:
                asums = dict(zip(self.df.id, self.df.asum2 * self.attr_to_m))

            # populate starting asums (asum1)
            routing_r = {v: k for k, v in self.routing.items() if v != 0}
            self.df['asum1'] = [asums.get(routing_r.get(id, 0), 0) for id in self.df.id.values]
            asum1s = dict(zip(self.df.id, self.df.asum1))

            # compute arbolate sum at reach midpoints (in meters)
            lengths = rd[['line_id', 'ireach', 'geometry']].copy()
            lengths['rchlen'] = np.array([g.length for g in lengths.geometry]) * self.geometry_to_m
            groups = lengths.groupby('line_id')  # fragments grouped by parent line

            # segment_asums = [asums[id] for id in lengths.line_id]
            # reach_cumsums = np.concatenate([np.cumsum(grp.rchlen.values[::-1])[::-1] - 0.5*grp.rchlen.values
            #                              for s, grp in groups])
            # reach_asums = -1 * reach_cumsums + segment_asums
            reach_cumsums = []
            ordered_ids = rd.line_id.loc[rd.line_id.diff() != 0].values
            for id in ordered_ids:
                grp = groups.get_group(id).sort_values(by='ireach')
                dist = np.cumsum(grp.rchlen.values) - 0.5 * grp.rchlen.values
                reach_cumsums.append(dist)
            reach_cumsums = np.concatenate(reach_cumsums)
            segment_asums = [asum1s[id] for id in lengths.line_id]
            reach_asums = segment_asums + reach_cumsums
            # maintain positive asums; lengths in NHD often aren't exactly equal to feature lengths
            # reach_asums[reach_asums < 0.] = 0
            rd['asum'] = reach_asums
            # width estimation formula expects meters
            width = width_from_arbolate_sum(reach_asums)
            rd['width'] = width * convert_length_units('meters', model_length_units) # convert back to model units
            rd.loc[rd.width < minimum_reach_width, 'width'] = minimum_reach_width

            # assign width1 and width2 back to segment data
            self.df['width1'] = width_from_arbolate_sum(self.df.asum1.values) * mult
            self.df['width2'] = width_from_arbolate_sum(self.df.asum2.values) * mult

        # interpolate linestring end widths to intersected reaches
        else:
            # verify that each linestring has only 1 segment associated with it
            # (interpolation might be wrong for multiple segments otherwise)
            assert rd.groupby('line_id').iseg.nunique().max() == 1
            # sort the linestring and reach data so that they are aligned
            self.df.sort_values(by='id', inplace=True)
            rd.sort_values(by=['line_id', 'ireach'], inplace=True)
            rd['width'] = interpolate_to_reaches(reach_data=rd,
                                                 segment_data=self.df,
                                                 segvar1='width1', segvar2='width2',
                                                 reach_data_group_col='line_id',
                                                 segment_data_group_col='id'
                                                 ) * mult

        # discard very small reaches; redo numbering
        # set minimum reach length based on cell size
        thresh = 0.05  # fraction of cell length (based on square root of area)
        if minimum_reach_length is None:
            cellgeoms = grid.df.loc[rd.node.values, 'geometry']
            mean_area = np.mean([g.area for g in cellgeoms])
            minimum_reach_length = np.sqrt(mean_area) * thresh * gis_mult

        inds = rd.rchlen > minimum_reach_length
        print('\nDropping {} reaches with length < {:.2f} {}...'.format(np.sum(~inds),
                                                                        minimum_reach_length,
                                                                        model_length_units))
        rd = rd.loc[inds].copy()
        rd['strhc1'] = 1.  # default value of streambed Kv for now
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
        paths = self.paths.copy()
        # for each segment
        for k in remaining_ids:
            # interate through successive downstream segments
            for s in paths[k][1:]:
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
        # segment2 = {lid: r[s] for lid, s in segment.items()}
        # line_id2 = {s: lid for lid, s in segment2.items()}

        # update reach_data
        rd['iseg'] = [r[s] for s in rd.iseg]

        print('\nSetting up segment data...')
        sd = pd.DataFrame()
        # sd['nseg'] = [segment2[lid] for lid in remaining_ids]
        # sd['outseg'] = [segment2.get(new_routing[line_id2[s]], 0) for s in nseg]
        sd['nseg'] = [r[s] for s in nseg]
        sd['outseg'] = [r[s] for s in outseg]
        sd.sort_values(by='nseg', inplace=True)

        # verify that no segments route to themselves
        assert not routing_is_circular(sd.nseg, sd.outseg)

        # (elevup dict was created above)
        elevup = self.elevup
        elevdn = dict(zip(self.df.id, self.df.elevdn))
        sd['elevup'] = [elevup[line_id[s]] for s in sd.nseg]
        sd['elevdn'] = [elevdn[line_id[s]] for s in sd.nseg]
        # convert elevation units
        sd['elevup'] *= mult_h
        sd['elevdn'] *= mult_h

        # apply widths if they were included
        if self.df[['width1', 'width2']].sum().sum() > 0:
            width1 = dict(zip(self.df.id, self.df.width1))
            width2 = dict(zip(self.df.id, self.df.width2))
            sd['width1'] = [width1[line_id[s]] for s in sd.nseg]
            sd['width2'] = [width2[line_id[s]] for s in sd.nseg]
            sd['width1'] *= mult
            sd['width2'] *= mult  # convert length units from source data to model
        elif self.df.width2.sum() == 0:
            raise NotImplementedError('Need to supply width1 and width2 or use arbolate sum.')

        # create sfrdata instance
        # this class has methods for fix segment and reach numbering,
        # assigning elevations and other properties by reach,
        # smoothing elevations, writing sfr package files
        # and other output
        rd = rd[[c for c in SFRData.rdcols if c in rd.columns]].copy()
        sfrd = SFRData(reach_data=rd, segment_data=sd, grid=grid,
                       model=model, model_length_units=model_length_units,
                       model_name=model_name, **kwargs)
        print("\nTime to create sfr dataset: {:.2f}s\n".format(time.time() - totim))
        return sfrd
