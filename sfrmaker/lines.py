import os
import time

import numpy as np
import pandas as pd
from shapely.geometry import box
import flopy
from gisutils import shp2df, df2shp, project, get_authority_crs
import sfrmaker
from sfrmaker.routing import pick_toids, find_path, make_graph, renumber_segments
from sfrmaker.checks import routing_is_circular, is_to_one
from sfrmaker.gis import read_polygon_feature, get_bbox, get_crs
from sfrmaker.grid import StructuredGrid
from sfrmaker.nhdplus_utils import load_nhdplus_v2, get_prj_file
from sfrmaker.sfrdata import SFRData
from sfrmaker.units import convert_length_units, get_length_units
from sfrmaker.utils import (width_from_arbolate_sum, arbolate_sum)
from sfrmaker.reaches import consolidate_reach_conductances, interpolate_to_reaches, setup_reach_data
from sfrmaker.routing import get_previous_ids_in_subset


class Lines:
    """Class for working with linestring feature input.

    Parameters
    ----------
    df : DataFrame
        Dataframe with linestring features and attribute information.
        Must have the following columns:
        
        ============ =====================================================
        **id**       sequential integers for identifying each feature
        **toid**     integers representing routing connections
        **geometry** shapely :class:`LineString` objects for each feature
        ============ =====================================================
        
    attr_length_units : str, {'meters', 'feet', ..}
        Length units for feature attributes (e.g. width, arbolate sum, etc.)
        (default 'meters')
    attr_height_units : str, {'meters', 'feet', ..}
        Length units for elevation attributes
        (default 'meters')
    crs : obj
        Coordinate reference system for features on df.geometry.
        A Python int, dict, str, or :class:`pyproj.crs.CRS` instance
        passed to :meth:`pyproj.crs.CRS.from_user_input`

        Can be any of:
          - PROJ string
          - Dictionary of PROJ parameters
          - PROJ keyword arguments for parameters
          - JSON string with PROJ parameters
          - CRS WKT string
          - An authority string [i.e. 'epsg:4326']
          - An EPSG integer code [i.e. 4326]
          - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
          - An object with a `to_wkt` method.
          - A :class:`pyproj.crs.CRS` class

        By default, None
        
    Attributes
    ----------
    crs : :class:`sfrmaker.gis.CRS` instance
    """

    def __init__(self, df=None,
                 attr_length_units='meters',
                 attr_height_units='meters',
                 crs=None, epsg=None, proj_str=None, prjfile=None):

        self.df = df
        self.attr_length_units = attr_length_units
        self.attr_height_units = attr_height_units
        self.crs = get_crs(prjfile=prjfile, epsg=epsg, proj_str=proj_str, crs=crs)
        self._geometry_length_units = None

        self._routing = None  # dictionary of routing connections
        self._paths = None  # routing sequence from each segment to outlet

        # dictionary of elevations at the upstream ends of flowlines
        self.elevup = dict(zip(self.df.id, self.df.elevup))
        # static dictionary of original flowline routing
        self._original_routing = self.routing.copy()


    @property
    def geometry_length_units(self):
        """Length units of reach LineString geometries.
        """
        valid_units = {'feet': 'feet',
                       'foot': 'feet',
                       'meters': 'meters',
                       'metre': 'meters'}
        self._geometry_length_units = valid_units.get(self.crs.axis_info[0].unit_name)
        if self._geometry_length_units is None:
            print("Warning: No length units specified in CRS for input LineStrings "
                  "or length units not recognized"
                  "defaulting to meters.")
            self._geometry_length_units = 'meters'
        return self._geometry_length_units

    @property
    def routing(self):
        """Dictionary of routing connections from ids (keys)
        to to_ids (values).
        """
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
        """Dictionary of paths, where each value is a list
        of downstream lines constituting a flow path to an outlet
        for a given line (key).
        """        
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
        simplify : bool
            Option to simplify the polygon, which can speed intersection 
            with the lines.
        tol: float
            Simplification tolerance (distance), in the units of the LineStrings
            (usually meters).
        feature_crs : obj
            A Python int, dict, str, or :py:class:`pyproj.crs.CRS` instance
            passed to the :py:meth:`pyproj.crs.CRS.from_user_input`
            See http://pyproj4.github.io/pyproj/stable/api/crs/crs.html#pyproj.crs.CRS.from_user_input.
            Can be any of:

              - PROJ string
              - Dictionary of PROJ parameters
              - PROJ keyword arguments for parameters
              - JSON string with PROJ parameters
              - CRS WKT string
              - An authority string [i.e. 'epsg:4326']
              - An EPSG integer code [i.e. 4326]
              - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
              - An object with a `to_wkt` method.
              - A :class:`pyproj.crs.CRS` class
        inplace : bool
            If True, the attribute .df is modified;
            if False, a copy of .df is returned.
            
        Returns
        -------
        df : DataFrame
            Version of the :py:attr:`Lines.df` DataFrame
            containing only the lines that intersect the ``feature``.
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

        # to_crs the flowlines if they aren't in same CRS as grid
        if self.crs != grid.crs:
            self.to_crs(grid.crs)

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

    def to_crs(self, dest_crs):
        """Reproject the LineStrings in :py:attr:`Lines.df` to
        a different Coordinate Reference System.

        Parameters
        ----------
        dest_crs : obj
            A Python int, dict, str, or :py:class:`pyproj.crs.CRS` instance
            passed to the :py:meth:`pyproj.crs.CRS.from_user_input`
            See http://pyproj4.github.io/pyproj/stable/api/crs/crs.html#pyproj.crs.CRS.from_user_input.
            Can be any of:

              - PROJ string
              - Dictionary of PROJ parameters
              - PROJ keyword arguments for parameters
              - JSON string with PROJ parameters
              - CRS WKT string
              - An authority string [i.e. 'epsg:4326']
              - An EPSG integer code [i.e. 4326]
              - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
              - An object with a `to_wkt` method.
              - A :class:`pyproj.crs.CRS` class
        """        
        assert self.crs is not None, "No crs for flowlines"
        assert dest_crs is not None, "No destination CRS."

        dest_crs = get_authority_crs(dest_crs)
        print('\nreprojecting hydrography from\n{}\nto\n{}\n'.format(self.crs,
                                                                     dest_crs))
        geoms = project(self.df.geometry, self.crs, dest_crs)
        assert np.isfinite(np.max(geoms[0].xy[0])), \
            "Invalid reprojection; check CRS for lines and grid."
        self.df['geometry'] = geoms
        self.crs = dest_crs

    def write_shapefile(self, outshp='flowlines.shp'):
        """Write a shapefile of :py:attr:`Lines.df`.

        Parameters
        ----------
        outshp : str, optional
            Shapefile name, by default 'flowlines.shp'
        """        
        df2shp(self.df, outshp, crs=self.crs)

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
                       crs=None, epsg=None, proj_str=None, prjfile=None):
        """Create a Lines instance from a shapefile.

        Parameters
        ----------
        shapefile : str
            Input shapefile
        id_column : str, optional
            Attribute field with line identifiers, 
            by default 'id'
        routing_column : str, optional
            Attribute field with downstream routing connections,
            by default 'toid'
        arbolate_sum_column2 : str, optional
            Attribute field with arbolate sums at downstream ends of lines, 
            by default 'asum2'
        width1_column : str, optional
            Attribute field with channel widths at upstream ends of lines,
            by default 'width1'
        width2_column : str, optional
            Attribute field with channel widths at downstream ends of lines, 
            by default 'width2'
        up_elevation_column : str, optional
            Attribute field with elevations at upstream ends of lines, 
            by default 'elevup'
        dn_elevation_column : str, optional
            Attribute field with elevations at downstream ends of lines,
            by default 'elevdn'
        name_column : str, optional
            Attribute field with feature names, 
            by default 'name'
        attr_length_units : str, optional
            Length units for feature attributes (e.g. width, arbolate sum, etc.)
            By default, meters.
        attr_height_units : str, optional
            Length units for elevation attributes
            By default, 'meters'.
        filter : tuple, optional
            (xmin, ymin, xmax, ymax) bounding box to filter which records 
            are read from the shapefile. By default None.
        epsg: int, optional
            EPSG code identifying Coordinate Reference System (CRS)
            for features in the input shapefile.
        proj_str: str, optional
            proj_str string identifying CRS for features in the input shapefile.
        prjfile: str, optional
            File path to projection (.prj) file identifying CRS
            for features in the input shapefile. By default,
            the projection file included with the input shapefile
            will be used.

        Returns
        -------
        lines : :class:`Lines` instance
        """        
        

        if prjfile is None:
            prjfile = shapefile.replace('.shp', '.prj')
            prjfile = prjfile if os.path.exists(prjfile) else None

        shpfile_crs = get_crs(prjfile=prjfile, epsg=epsg, proj_str=proj_str, crs=crs)

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
        """[summary]

        Parameters
        ----------
        df : DataFrame
            Pandas DataFrame with flowline information, including
            shapely :class:`LineStrings <LineString>` in a `'geometry'` column.
        id_column : str, optional
            Attribute field with line identifiers, 
            by default 'id'
        routing_column : str, optional
            Attribute field with downstream routing connections,
            by default 'toid'
        arbolate_sum_column2 : str, optional
            Attribute field with arbolate sums at downstream ends of lines, 
            by default 'asum2'
        width1_column : str, optional
            Attribute field with channel widths at upstream ends of lines,
            by default 'width1'
        width2_column : str, optional
            Attribute field with channel widths at downstream ends of lines, 
            by default 'width2'
        up_elevation_column : str, optional
            Attribute field with elevations at upstream ends of lines, 
            by default 'elevup'
        dn_elevation_column : str, optional
            Attribute field with elevations at downstream ends of lines,
            by default 'elevdn'
        name_column : str, optional
            Attribute field with feature names, 
            by default 'name'
        attr_length_units : str, optional
            Length units for feature attributes (e.g. width, arbolate sum, etc.)
            By default, meters.
        attr_height_units : str, optional
            Length units for elevation attributes
            By default, 'meters'.
        filter : tuple, optional
            (xmin, ymin, xmax, ymax) bounding box to filter which records 
            are read from the shapefile. By default None.
        epsg: int, optional
            EPSG code identifying Coordinate Reference System (CRS)
            for features in the input shapefile.
        proj_str: str, optional
            proj_str string identifying CRS for features in the input shapefile.
        prjfile: str, optional
            File path to projection (.prj) file identifying CRS
            for features in the input shapefile. By default,
            the projection file included with the input shapefile
            will be used.

        Returns
        -------
        lines : :class:`Lines` instance
        """
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
        NHDPlus_paths : str or list of strings
            List of paths to the root folders of NHDPlus drainage basins
            to include, assuming the file structure is the same as
            downloaded from the NHDPlus version 2 website. For example::
            
                NHDPlus_paths=['/NHDPlusGL/NHDPlus04/',
                               '/NHDPlusMS/NHDPlus07/']    
                                     
            for the Great Lakes (04) and Upper Mississippi (07) basins.      
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
        epsg: int, optional
            EPSG code identifying Coordinate Reference System (CRS)
            for features in the input shapefile.
        proj_str: str, optional
            proj_str string identifying CRS for features in the input shapefile.
        prjfile: str, optional
            File path to projection (.prj) file identifying CRS
            for features in the input shapefile. By default,
            the projection file included with the input shapefile
            will be used.

        Returns
        -------
        lines : :class:`Lines` instance
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

    def to_sfr(self, grid=None,
               active_area=None, isfr=None,
               model=None,
               model_length_units='undefined',
               model_time_units='days',
               minimum_reach_length=None,
               width_from_asum_a_param=0.1193,
               width_from_asum_b_param=0.5032,
               minimum_reach_width=1.,
               consolidate_conductance=False, one_reach_per_cell=False,
               add_outlets=None,
               package_name=None,
               **kwargs):
        """Create a streamflow routing dataset from the information
        in sfrmaker.lines class instance and a supplied sfrmaker.grid class instance.

        Parameters
        ----------
        grid : sfrmaker.grid or flopy.discretization.StructuredGrid
            Numerica model grid instance. Required unless an attached model
            has a valid modelgrid attribute.
        active_area : shapely Polygon, list of shapely Polygons, or shapefile path; optional
            Shapely Polygons must be in same CRS as input flowlines; shapefile
            features will be reprojected if their crs is different.
        isfr : ndarray, optional
            Numpy integer array of the same size as the model grid, designating area that will
            be populated with SFR reaches (0=no SFR; 1=SFR). An isfr array of shape
            nrow x ncol will be broadcast to all layers. Only required if a model is not
            supplied, or if SFR is only desired in a subset of active model cells.
            By default, None, in which case the model ibound or idomain array will be used.
        model : flopy.modflow.Modflow or flopy.mf6.ModflowGwf, optional
            Flopy model instance
        model_length_units : str; e.g. {'ft', 'feet', 'meters', etc.}, optional
            Length units of the model. While SFRmaker will try to read these
            from a supplied grid (first) and then a supplied model (second),
            it is good practice to specify them explicitly here.
        model_time_units : str; e.g. {'d', 'days'}, optional
            Time units for model. By default, days.
        minimum_reach_length : float, optional
            Minimum reach length to retain. Default is to compute
            an effective mean model cell length by taking the square root
            of the average cell area, and then set minimum_reach_length
            to 5% of effective mean cell length.
        width_from_asum_a_param : float, optional
            :math:`a` parameter used for estimating channel width from arbolate sum.
            Only needed if input flowlines are lacking width information.
            See :func:`~sfrmaker.utils.width_from_arbolate`. By default, 0.1193.
        width_from_asum_b_param : float, optional
            :math:`b` parameter used for estimating channel width from arbolate sum.
            Only needed if input flowlines are lacking width information.
            See :func:`~sfrmaker.utils.width_from_arbolate`. By default, 0.5032.
        minimum_reach_width : float, optional
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
        add_outlets : sequence of ints
            Option to add breaks in routing at specified line ids. For example
            if controlled flows out of a reservoir are specified as inflows
            to the SFR network, an outlet can be added above to the dam to
            prevent double-counting of flow. By default, None
        package_name : str
            Base name for writing sfr output.
        kwargs : keyword arguments to :class:`SFRData`

        Returns
        -------
        sfrdata : sfrmaker.SFRData instance

        """
        print("\nSFRmaker version {}".format(sfrmaker.__version__))
        print("\nCreating sfr dataset...")
        totim = time.time()

        if flopy and active_area is None and isfr is None and model is not None:
            if model.version == 'mf6':
                isfr = np.sum(model.dis.idomain.array == 1, axis=0) > 0
            else:
                isfr = np.sum(model.bas6.ibound.array == 1, axis=0) > 0
        if flopy and isinstance(grid, flopy.discretization.StructuredGrid):
            print('\nCreating grid class instance from flopy Grid instance...')
            ta = time.time()
            grid = StructuredGrid.from_modelgrid(grid, active_area=active_area, isfr=isfr)
            print("grid class created in {:.2f}s\n".format(time.time() - ta))
        elif flopy and model is not None:
            grid = StructuredGrid.from_modelgrid(model.modelgrid, active_area=active_area, isfr=isfr)
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

        # to_crs the flowlines if they aren't in same CRS as grid
        if self.crs != grid.crs:
            self.to_crs(grid.crs)
        # cull the flowlines to the active part of the model grid
        if grid.active_area is not None:
            self.cull(grid.active_area, inplace=True, simplify=True, tol=2000)
        elif grid._bounds is not None:  # cull to grid bounding box if already computed
            self.cull(box(*grid._bounds), inplace=True)
        if package_name is None:
            if model is not None:
                package_name = model.name
            else:
                package_name = 'model'

        # convert routing connections (toid column) from lists (one-to-many)
        # to ints (one-to-one or many-to-one)
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
                                              np.array([g.length for g in self.df.geometry]) * convert_length_units(self.geometry_length_units, 'meters')
                                              )),
                                     self.routing)
            else:
                asums = dict(zip(self.df.id, 
                                 self.df.asum2 * convert_length_units(self.attr_length_units, 
                                                                      'meters')))

            # populate starting asums (asum1)
            routing_r = {v: k for k, v in self.routing.items() if v != 0}
            self.df['asum1'] = [asums.get(routing_r.get(id, 0), 0) for id in self.df.id.values]
            asum1s = dict(zip(self.df.id, self.df.asum1))

            # compute arbolate sum at reach midpoints (in meters)
            lengths = rd[['line_id', 'ireach', 'geometry']].copy()
            lengths['rchlen'] = np.array([g.length for g in lengths.geometry]) * convert_length_units(self.geometry_length_units, 'meters')
            groups = lengths.groupby('line_id')  # fragments grouped by parent line

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
            width = width_from_arbolate_sum(reach_asums,
                                            a=width_from_asum_a_param,
                                            b=width_from_asum_b_param,
                                            minimum_width=minimum_reach_width,
                                            input_units='meters', output_units=model_length_units)
            rd['width'] = width
            rd.loc[rd.width < minimum_reach_width, 'width'] = minimum_reach_width

            # assign width1 and width2 back to segment data
            self.df['width1'] = width_from_arbolate_sum(self.df.asum1.values,
                                            a=width_from_asum_a_param,
                                            b=width_from_asum_b_param,
                                            minimum_width=minimum_reach_width,
                                            input_units=self.attr_length_units,
                                            output_units=model_length_units)
            self.df['width2'] = width_from_arbolate_sum(self.df.asum2.values,
                                                        a=width_from_asum_a_param,
                                                        b=width_from_asum_b_param,
                                                        minimum_width=minimum_reach_width,
                                                        input_units=self.attr_length_units,
                                                        output_units=model_length_units)

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

        # add any outlets to the stream network
        # for now handle int or str ids
        if add_outlets is not None:
            # get the
            if isinstance(add_outlets, str) or isinstance(add_outlets, int):
                add_outlets = [add_outlets]
            for outlet_id in add_outlets:
                if rd.line_id.dtype == np.object:
                    outlet_id = str(outlet_id)
                    outlet_toid = '0'
                else:
                    outlet_id = int(outlet_id)
                    outlet_toid = 0
                valid_outlet_ids = get_previous_ids_in_subset(rd.line_id, self.routing, outlet_id)
                loc = rd.line_id.isin(valid_outlet_ids)
                rd.loc[loc, 'toid'] = outlet_toid
                for valid_outlet_id in valid_outlet_ids:
                    new_routing[valid_outlet_id] = outlet_toid

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
                       model_time_units=model_time_units,
                       package_name=package_name, **kwargs)
        print("\nTime to create sfr dataset: {:.2f}s\n".format(time.time() - totim))
        return sfrd
