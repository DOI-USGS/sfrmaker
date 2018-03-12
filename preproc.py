__author__ = 'aleaf'
import warnings
import time
import operator
import numpy as np
import pandas as pd
import fiona
from shapely.geometry import Point, LineString, Polygon, shape, box
from shapely.ops import unary_union
from GISio import shp2df, df2shp, get_proj4
from GISops import project, projectdf, build_rtree_index, intersect_rtree
import GISops

class linesBase(object):

    output_cols = ['node', 'k', 'i', 'j', 'segment', 'reach', 'strtop', 'width', 'length', 'strthick',
                   'strhc1', 'roughch', 'asum', 'reachID', 'outreach', 'comid']

    def __init__(self, lines=None,
                 mf_grid=None, mf_grid_node_col=None,
                 nrows=None, ncols=None,
                 mfdis=None, xul=None, yul=None, rot=0,
                 model_domain=None,
                 lines_proj4=None, mfgrid_proj4=None, domain_proj4=None,
                 mf_units='feet'):
        """Class for working with information from NHDPlus v2.
        See the user's guide for more information:
        <http://www.horizon-systems.com/NHDPlus/NHDPlusV2_documentation.php#NHDPlusV2 User Guide>

        Parameters
        ==========
        lines : str, list of strings or dataframe
            Shapefile, list of shapefiles, or dataframe with linework defining SFR network;
            assigned to the Flowline attribute.
        mf_grid : str or dataframe
            Shapefile or dataframe containing MODFLOW grid
        mf_grid_node_col : str
            Column in grid shapefile or dataframe with unique node numbers.
            In case the grid isn't sorted!
            (which will result in mixup if rows and columns are assigned later using the node numbers)
        nrows : int
            (structured grids) Number of model rows
        ncols : int
            (structured grids) Number of model columns
        mfdis : str
            MODFLOW discretization file (not yet supported for this class)
        xul : float, optional
            x offset of upper left corner of grid. Only needed if using mfdis instead of shapefile
        yul : float, optional
            y offset of upper left corner of grid. Only needed if using mfdis instead of shapefile
        rot : float, optional (default 0)
            Grid rotation; only needed if using mfdis instead of shapefile.
        model_domain : str (shapefile) or shapely polygon, optional
            Polygon defining area in which to create SFR cells.
            Default is to create SFR at all intersections between the model grid and NHD flowlines.
        lines_proj4 : str, optional
            Proj4 string for coordinate system of NHDFlowlines.
            Only needed if flowlines are supplied in a dataframe.
        domain_proj4 : str, optional
            Proj4 string for coordinate system of model_domain.
            Only needed if model_domain is supplied as a polygon.
        mf_units : str, 'feet' or 'meters'
            Length units of MODFLOW model
        """
        self.df = lines
        self.mf_grid = mf_grid
        self.model_domain = model_domain
        self.nrow = nrows
        self.ncol = ncols
        self.mfdis = mfdis
        self.xul = xul
        self.yul = yul
        self.rot = rot

        # unit conversions (set below after grid projection is verified)
        self.mf_units = mf_units
        self.mf_units_mult = 1.0 # go from GIS units to model units
        self.GISunits = None #
        self.to_km = None # converts GIS units to km for arbolate sum

        self.proj4 = lines_proj4
        self.mf_grid_proj4 = mfgrid_proj4
        self.domain_proj4 = domain_proj4

        print("Reading input...")
        # handle dataframes or shapefiles as arguments
        # get proj4 for any shapefiles that are submitted
        for attr, input in {'df': lines,
                            'grid': mf_grid}.items():
            if isinstance(input, pd.DataFrame):
                self.__dict__[attr] = input
            else:
                self.__dict__[attr] = shp2df(input)
        if isinstance(model_domain, Polygon):
            self.domain = model_domain
        elif isinstance(model_domain, str):
            self.domain = shape(fiona.open(model_domain).next()['geometry'])
            self.domain_proj4 = get_proj4(model_domain)
        else:
            print('setting model domain to extent of grid ' \
                  'by performing unary union of grid cell geometries...\n' \
                  '(may take a few minutes for large grids)')
            # add tiny buffer to overcome floating point errors in gridcell geometries
            # (otherwise a multipolygon feature may be returned)
            geoms = [g.buffer(0.001) for g in self.grid.geometry.tolist()]
            self.domain = unary_union(geoms)

        # sort and pair down the grid
        if mf_grid_node_col is not None:
            self.grid.sort_values(by=mf_grid_node_col, inplace=True)
            if len(self.grid) != self.grid[mf_grid_node_col].max():
                warnings.warn(NodeIndexWarning(mf_grid, mf_grid_node_col))
            self.grid.index = self.grid[mf_grid_node_col].values
        else:
            warnings.warn(NodeIndexWarning(mf_grid))
        self.grid = self.grid[['geometry']]

        # get projections
        if self.mf_grid_proj4 is None and not isinstance(mf_grid, pd.DataFrame):
            self.mf_grid_proj4 = get_proj4(mf_grid)
        if self.proj4 is None:
            if isinstance(lines, list):
                self.proj4 = get_proj4(lines[0])
            elif not isinstance(lines, pd.DataFrame):
                self.proj4 = get_proj4(lines)

        # first check that grid is in projected units
        if self.mf_grid_proj4.split('proj=')[1].split()[0].strip() == 'longlat':
            raise ProjectionError(self.mf_grid)

        # reproject the NHD Flowlines and model domain to model grid if they aren't
        # (prob a better way to check for same projection)

        # set GIS units from modflow grid projection (used for arbolate sum computation)
        # assumes either m or ft!
        self.GISunits = parse_proj4_units(self.mf_grid_proj4)
        self.mf_units_mult = 1/0.3048 if self.GISunits == 'm' and self.mf_units == 'feet' \
                             else 0.3048 if not self.GISunits == 'm' and self.mf_units == 'meters' \
                             else 1.0
        self.to_km = 0.001 if self.GISunits == 'm' else 0.001/0.3048


        if different_projections(self.proj4, self.mf_grid_proj4):
            print("reprojecting NHDFlowlines from\n{}\nto\n{}...".format(self.proj4, self.mf_grid_proj4))
            self.df['geometry'] = projectdf(self.df, self.proj4, self.mf_grid_proj4)

        if model_domain is not None \
                and different_projections(self.domain_proj4, self.mf_grid_proj4):
            print("reprojecting model domain from\n{}\nto\n{}...".format(self.domain_proj4, self.mf_grid_proj4))
            self.domain = project(self.domain, self.domain_proj4, self.mf_grid_proj4)

    def renumber_segments(self):
        """Renumber segments so that segment numbering is continuous and always increases
        in the downstream direction. Experience suggests that this can substantially speed
        convergence for some models using the NWT solver."""
        r = renumber_segments(self.m2.segment.values, self.m2.outseg.values)

        self.m2['segment'] = [r.get(s, s) for s in self.m2.segment]
        self.m2['outseg'] = [r.get(s, s) for s in self.m2.outseg]
        self.m2.sort_values(by=['segment'], inplace=True)
        self.m2.index = self.m2.segment.values # reset the index to new segment numbers
        assert _in_order(self.m2.segment.values, self.m2.outseg.values)
        assert len(self.m2.segment) == self.m2.segment.max()
        self.m1['segment'] = [r.get(s, s) for s in self.m1.segment]
        self.m1['outseg'] = [r.get(s, s) for s in self.m1.outseg]

    def write_tables(self, basename='SFR', m1=None, m2=None):
        """Write tables with SFR reach (Mat1) and segment (Mat2) information out to csv files.

        Parameters
        ----------
        basename: string
            e.g. Mat1 is written to <basename>Mat1.csv
        """
        #m2_cols = ['segment', 'icalc', 'outseg', 'elevup', 'elevdn']
        if m1 is None:
            m1 = self.m1
        if m2 is None:
            m2 = self.m2
        m1_cols = self.output_cols
        m2_cols = ['segment', 'icalc', 'outseg', 'elevup', 'elevdn', 'width1', 'width2', 'comid']
        for ind in ['i', 'j']:
            if ind not in m1.columns:
                m1_cols.remove(ind)
        print("writing Mat1 to {0}{1}, Mat2 to {0}{2}".format(basename, 'Mat1.csv', 'Mat2.csv'))
        m1[m1_cols].to_csv(basename + 'Mat1.csv', index=False)
        m2[m2_cols].to_csv(basename + 'Mat2.csv', index=False)

    def write_linework_shapefile(self, basename='SFR'):
        """Write a shapefile containing linework for each SFR reach,
        with segment, reach, model node number, and NHDPlus COMID attribute information

        Parameters
        ----------
        basename: string
            Output will be written to <basename>.shp
        """
        outfile = basename.split('.')[0] + '.shp'
        print("writing reach geometries to {}".format(outfile))
        df2shp(self.m1[['reachID', 'node', 'segment', 'reach', 'outseg', 'comid', 'asum', 'width', 'geometry']],
               outfile, proj4=self.mf_grid_proj4)


class NHDdata(object):

    convert_elevslope_to_model_units = {'feet': 0.0328084,
                                        'meters': 0.01}

    output_cols = ['node', 'k', 'i', 'j', 'segment', 'reach', 'strtop', 'width', 'length', 'strthick',
                   'strhc1', 'roughch', 'asum', 'reachID', 'outreach', 'comid']

    def __init__(self, NHDFlowline=None, PlusFlowlineVAA=None, PlusFlow=None, NHDFcode=None,
                 elevslope=None,
                 elevup_col=None, elevdn_col=None, routing_col=None, arbolate_sum_col=None,
                 mf_grid=None, mf_grid_node_col=None,
                 nrows=None, ncols=None,
                 mfdis=None, xul=None, yul=None, rot=0, sr=None,
                 model_domain=None, filter=True,
                 lines_proj4=None, mfgrid_proj4=None, domain_proj4=None,
                 mf_units='feet'):
        """Class for working with information from NHDPlus v2.
        See the user's guide for more information:
        <http://www.horizon-systems.com/NHDPlus/NHDPlusV2_documentation.php#NHDPlusV2 User Guide>

        Parameters
        ==========
        NHDFlowline : str, list of strings or dataframe
            Shapefile, list of shapefiles, or dataframe defining SFR network;
            assigned to the Flowline attribute.
        PlusFlowlineVAA : str, list of strings or dataframe
            DBF file, list of DBF files with NHDPlus attribute information;
            assigned to PlusFlowlineVAA attribute.
        PlusFlow : str, list of strings or dataframe
            DBF file, list of DBF files with routing information;
            assigned to PlusFlow attribute.
        mf_grid : str or dataframe
            Shapefile or dataframe containing MODFLOW grid
        mf_grid_node_col : str
            Column in grid shapefile or dataframe with unique node numbers.
            In case the grid isn't sorted!
            (which will result in mixup if rows and columns are assigned later using the node numbers)
        nrows : int
            (structured grids) Number of model rows
        ncols : int
            (structured grids) Number of model columns
        mfdis : str
            MODFLOW discretization file (not yet supported for this class)
        xul : float, optional
            x offset of upper left corner of grid. Only needed if using mfdis instead of shapefile
        yul : float, optional
            y offset of upper left corner of grid. Only needed if using mfdis instead of shapefile
        rot : float, optional (default 0)
            Grid rotation; only needed if using mfdis instead of shapefile.
        model_domain : str (shapefile) or shapely polygon, optional
            Polygon defining area in which to create SFR cells.
            Default is to create SFR at all intersections between the model grid and NHD flowlines.
        lines_proj4 : str, optional
            Proj4 string for coordinate system of NHDFlowlines.
            Only needed if flowlines are supplied in a dataframe.
        domain_proj4 : str, optional
            Proj4 string for coordinate system of model_domain.
            Only needed if model_domain is supplied as a polygon.
        mf_units : str, 'feet' or 'meters'
            Length units of MODFLOW model
        """
        self.Flowline = NHDFlowline
        self.PlusFlowlineVAA = PlusFlowlineVAA

        self.PlusFlow = PlusFlow
        self.elevslope = elevslope
        self.fl_cols = ['COMID', 'FCODE', 'FDATE', 'FLOWDIR',
                          'FTYPE', 'GNIS_ID', 'GNIS_NAME', 'LENGTHKM',
                          'REACHCODE', 'RESOLUTION', 'WBAREACOMI', 'geometry']
        self.pfvaa_cols = ['ArbolateSu', 'Hydroseq', 'DnHydroseq',
                      'LevelPathI', 'StreamOrde']

        self.mf_grid = mf_grid
        self.model_domain = model_domain
        self.nrow = nrows
        self.ncol = ncols
        self.mfdis = mfdis
        self.xul = xul
        self.yul = yul
        self.rot = rot

        # unit conversions (set below after grid projection is verified)
        self.mf_units = mf_units
        self.mf_units_mult = 1.0 # go from GIS units to model units
        self.GISunits = None #
        self.to_km = None # converts GIS units to km for arbolate sum

        self.fl_proj4 = lines_proj4
        self.mf_grid_proj4 = mfgrid_proj4
        self.domain_proj4 = domain_proj4

        print("Reading input...")

        # get projections
        if self.mf_grid_proj4 is None and not isinstance(mf_grid, pd.DataFrame):
            self.mf_grid_proj4 = get_proj4(mf_grid)
        if self.fl_proj4 is None:
            if isinstance(NHDFlowline, list):
                self.fl_proj4 = get_proj4(NHDFlowline[0])
            elif not isinstance(NHDFlowline, pd.DataFrame):
                self.fl_proj4 = get_proj4(NHDFlowline)

        # model domain
        if isinstance(model_domain, Polygon):
            self.domain = model_domain
        elif isinstance(model_domain, str):
            self.domain = shape(fiona.open(model_domain).next()['geometry'])
            self.domain_proj4 = get_proj4(model_domain)
        else:
            print('setting model domain to extent of grid ' \
                  'by performing unary union of grid cell geometries...\n' \
                  '(may take a few minutes for large grids)')
            # add tiny buffer to overcome floating point errors in gridcell geometries
            # (otherwise a multipolygon feature may be returned)
            geoms = [g.buffer(0.001) for g in self.grid.geometry.tolist()]
            self.domain = unary_union(geoms)

        # model grid
        if sr is not None:
            print('reading grid from flopy SpatialReference...')
            d = {'geometry': [Polygon(vrts) for vrts in sr.vertices],
                 'j': np.arange(sr.ncol).tolist() * sr.nrow,
                 'i': sorted(np.arange(sr.nrow).tolist() * sr.ncol),
                 'node': np.arange(sr.nrow*sr.ncol, dtype=int)}
            self.grid = pd.DataFrame(d)
            self.grid['i'] = self.grid.i.astype(int)
            self.grid['j'] = self.grid.j.astype(int)
            self.grid = self.grid[['node', 'i', 'j', 'geometry']]
            self.nrow = sr.nrow
            self.ncol = sr.ncol
            mf_grid_node_col = 'node'

        # handle dataframes or shapefiles as arguments
        # get proj4 for any shapefiles that are submitted
        for attr, input in {'fl': NHDFlowline,
                            'pf': PlusFlow,
                            'pfvaa': PlusFlowlineVAA,
                            'elevs': elevslope,
                            'grid': mf_grid}.items():
            if isinstance(input, pd.DataFrame):
                self.__dict__[attr] = input
            else:
                # filter flowlines to speed reading them in
                if input is None:
                    continue
                if attr == 'fl' and filter:
                    if model_domain is not None \
                            and different_projections(self.domain_proj4, self.fl_proj4):
                        print(self.domain_proj4)
                        print(self.fl_proj4)
                        self.domain_filter = project(self.domain, self.domain_proj4, self.fl_proj4).bounds
                        print(self.domain_filter)
                    else:
                        self.domain_filter = self.domain.bounds
                else:
                    self.domain_filter = None
                self.__dict__[attr] = shp2df(input, filter=self.domain_filter)

        # set the indices
        for attr in ['fl', 'pfvaa', 'elevs']:
            if attr == 'elevs' and elevslope is None:
                continue
            comid_col = [c for c in self.__dict__[attr].columns
                         if c.lower() == 'comid'][0]
            self.__dict__[attr].dropna(subset=[comid_col], axis=0, inplace=True)
            self.__dict__[attr][comid_col] = self.__dict__[attr][comid_col].astype(int)
            if not self.__dict__[attr].index.name == comid_col:
                self.__dict__[attr].index = self.__dict__[attr][comid_col]

        # create a working dataframe
        self.df = self.fl[self.fl_cols].join(self.pfvaa[self.pfvaa_cols], how='inner')

        # elevations already supplied in flowlines file
        if elevslope is None:
            if elevdn_col in self.fl.columns and elevup_col in self.fl.columns:
                self.elevs = self.fl[[elevup_col, elevdn_col]]
                self.elevs.rename(columns={elevdn_col: 'Min',
                                           elevup_col: 'Max'}, inplace=True)
            else:
                print('No elevslope.dbf file(s) or up/down elevations columns supplied!')
        else: # bring in elevations from elevslope.dbf
            # convert the elevations from elevslope table
            self.elevs['Max'] = self.elevs.MAXELEVSMO * self.convert_elevslope_to_model_units[self.mf_units]
            self.elevs['Min'] = self.elevs.MINELEVSMO * self.convert_elevslope_to_model_units[self.mf_units]
        self.df = self.df.join(self.elevs[['Max', 'Min']], how='inner')
        self.df.rename(columns={'Max': 'elevup', 'Min': 'elevdn'}, inplace=True)

        # routing information provided in flowlines file
        if PlusFlow is None:
            if routing_col in self.fl.columns:
                graph = dict(zip(self.fl.index, self.fl[routing_col]))
                graph_r = {}
                allcomids = set(self.fl[routing_col].tolist() + self.fl.index.tolist())
                for toid in allcomids:
                    frmids = [k for k, v in graph.items() if v == toid]
                    graph_r[toid] = frmids
                self.df['dncomids'] = [[graph[c]] for c in self.df.index]
                self.df['upcomids'] = [graph_r[c] for c in self.df.index]
            else:
                print('No PlusFlow.dbf file(s) or routing column supplied!')

        # arbolate sum information provided in flowlines file
        if arbolate_sum_col is not None:
            self.df['ArbolateSu'] = self.fl.loc[self.df.index, arbolate_sum_col]

        # sort and pair down the grid
        if mf_grid_node_col is not None:
            self.grid.sort_values(by=mf_grid_node_col, inplace=True)
            if len(self.grid) != self.grid[mf_grid_node_col].max() +1:
                warnings.warn(NodeIndexWarning(mf_grid, mf_grid_node_col))
            self.grid.index = self.grid[mf_grid_node_col].values
        else:
            warnings.warn(NodeIndexWarning(mf_grid))

        # first check that grid is in projected units
        if self.mf_grid_proj4.split('proj=')[1].split()[0].strip() == 'longlat':
            raise ProjectionError(self.mf_grid)

        # reproject the NHD Flowlines and model domain to model grid if they aren't
        # (prob a better way to check for same projection)

        # set GIS units from modflow grid projection (used for arbolate sum computation)
        # assumes either m or ft!
        self.GISunits = parse_proj4_units(self.mf_grid_proj4)
        self.mf_units_mult = 1/0.3048 if self.GISunits == 'm' and self.mf_units == 'feet' \
                             else 0.3048 if not self.GISunits == 'm' and self.mf_units == 'meters' \
                             else 1.0
        self.to_km = 0.001 if self.GISunits == 'm' else 0.001 * 0.3048



        if different_projections(self.fl_proj4, self.mf_grid_proj4):
            print("reprojecting NHDFlowlines from\n{}\nto\n{}...".format(self.fl_proj4, self.mf_grid_proj4))
            self.df['geometry'] = project(self.df.geometry, self.fl_proj4, self.mf_grid_proj4)

        if model_domain is not None \
                and different_projections(self.domain_proj4, self.mf_grid_proj4):
            print("reprojecting model domain from\n{}\nto\n{}...".format(self.domain_proj4, self.mf_grid_proj4))
            self.domain = project(self.domain, self.domain_proj4, self.mf_grid_proj4)



    def list_updown_comids(self):
        print('getting routing information from NHDPlus Plusflow table...')
        # setup local variables and cull plusflow table to comids in model
        comids = self.df.index.tolist()
        pf = self.pf.loc[(self.pf.FROMCOMID.isin(comids)) |
                        (self.pf.TOCOMID.isin(comids))].copy()

        # subset PlusFlow entries for comids that are not in flowlines dataset
        # comids may be missing because they are outside of the model
        # or if the flowlines dataset was edited (resulting in breaks in the routing)
        missing_tocomids = ~pf.TOCOMID.isin(comids) & (pf.TOCOMID != 0)
        missing = pf.loc[missing_tocomids, ['FROMCOMID', 'TOCOMID']].copy()
        # recursively crawl the PlusFlow table
        # to try to find a downstream comid in the flowlines dataest
        missing['nextCOMID'] = [find_next(tc, self.pf, comids) for tc in missing.TOCOMID]
        pf.loc[missing_tocomids, 'TOCOMID'] = missing.nextCOMID

        # set any remaining comids not in model to zero
        # (outlets or inlets from outside model)
        #pf.loc[~pf.TOCOMID.isin(comids), 'TOCOMID'] = 0 (these should all be handled above)
        pf.loc[~pf.FROMCOMID.isin(comids), 'FROMCOMID'] = 0
        tocomid = pf.TOCOMID.values
        fromcomid = pf.FROMCOMID.values
        self.df['dncomids'] = [tocomid[fromcomid == c].tolist() for c in comids]
        self.df['upcomids'] = [fromcomid[tocomid == c].tolist() for c in comids]

    def assign_segments(self):
        print('assigning segment numbers...')
        # create segment numbers
        self.df['segment'] = np.arange(len(self.df)) + 1

        # dncomids must all be lists
        assert np.all([isinstance(c, list) for c in self.df.dncomids])

        # reduce dncomids to 1 per segment
        # divergences are where there is more than 1 dncomid
        diverges = np.array([len(d) for d in self.df.dncomids]) > 1
        diverging_dncomids = self.df.loc[diverges, 'dncomids'].tolist()
        all_comids = set(self.df.index)
        dncomids_at_divergences = []
        #for i, r in divergences.iterrows():
        # iterate through segments with a divergence and pick
        # dncomid with lowest elevation
        for dncomids in diverging_dncomids:
            # select the dncomid that has a matching levelpath
            #in_index = set(r.dncomids).intersection(set(self.df.index))
            in_index = set(dncomids).intersection(all_comids)
            if len(in_index) > 0:
                #df = self.df.loc[r.dncomids].sort_values(by='elevdn')
                # sort dncomids list by elevation (ascending)
                df = self.df.loc[dncomids].sort_values(by='elevdn')
                dncomid = [df.COMID.values[0]]
            else:
                dncomid = [0]
            dncomids_at_divergences.append(dncomid)
            '''
            if len(in_index) > 0:
                levelpath_matches = self.df.loc[r.dncomids, 'LevelPathI'].values == r.LevelPathI
                in_same_levelpath = np.array(r.dncomids)[levelpath_matches]
                
            else:
                in_same_levelpath = []
            # if none match, select the first dncomid
            if len(in_same_levelpath) == 0:
                
                dncomid = [r.dncomids[0]]
            else:
                dncomid = np.unique(in_same_levelpath).tolist()
            '''
            #self.df.set_value(i, 'dncomids', dncomid)
        assert np.all([isinstance(c, list) for c in dncomids_at_divergences])
        all_dncomids = self.df.dncomids.tolist()
        for i, dncomid in enumerate(all_dncomids):
            if diverges[i]:
                all_dncomids[i] = dncomids_at_divergences.pop(0)
        assert len(dncomids_at_divergences) == 0
        #all_dncomids[diverges] = dncomids_at_divergences
        self.df['dncomids'] = all_dncomids
        #self.df.loc[diverges, 'dncomids'] = dncomids_at_divergences
        assert np.all([isinstance(c, list) for c in self.df.dncomids.values])

        # assign upsegs and outsegs based on NHDPlus routing
        segments = dict(zip(self.df.index, self.df.segment))
        #self.df['upsegs'] = [[self.df.segment[c] if c != 0 else 0 for c in comids] for comids in self.df.upcomids]
        #self.df['dnsegs'] = [[self.df.segment[c] if c != 0 else 0 for c in comids] for comids in self.df.dncomids]
        self.df['upsegs'] = [[segments.get(c, 0) for c in comids] for comids in self.df.upcomids]
        for comids in self.df.dncomids:
            if not isinstance(comids, list):
                break
        self.df['dnsegs'] = [[segments.get(c, 0) for c in comids] for comids in self.df.dncomids]

        # make a column of outseg integers
        self.df['outseg'] = [d[0] for d in self.df.dnsegs]

    def to_sfr(self, roughch=0.037, strthick=1, strhc1=1,
               icalc=1, minimum_length=1.,
               consolidate_conductance=False, one_reach_per_cell=False,
               iupseg=0, iprior=0, nstrpts=0, flow=0, runoff=0, etsw=0, pptsw=0,
               roughbk=0, cdepth=0, fdepth=0, awdth=0, bwdth=0):

        print('\nclipping flowlines to active area... (may take awhile for active area polygons with many vertices)')
        # this step is slow for domains with many vertices
        # (i.e. those created from an ibound array)
        # would be ideal to discard this step and intersect directly with grid using rtree
        ta = time.time()
        inside = np.array([g.intersects(self.domain) for g in self.df.geometry])
        self.df = self.df.loc[inside].copy()
        self.df.sort_values(by='COMID', inplace=True)
        flowline_geoms = [g.intersection(self.domain) for g in self.df.geometry]
        grid_geoms = self.grid.geometry.tolist()
        print("finished in {:.2f}s\n".format(time.time() - ta))

        print("setting up segments... (may take a few minutes for large networks)")
        ta = time.time()
        if 'upcomids' not in self.df.columns or 'dncomids' not in self.df.columns:
            self.list_updown_comids()
        self.assign_segments()
        fl_segments = self.df.segment.tolist()
        fl_comids = self.df.COMID.tolist()
        print("finished in {:.2f}s\n".format(time.time() - ta))

        print("intersecting flowlines with grid cells...") # this part crawls in debug mode
        grid_intersections = GISops.intersect_rtree(grid_geoms, flowline_geoms)

        print("setting up reaches and Mat1... (may take a few minutes for large grids)")
        ta = time.time()
        m1 = make_mat1(flowline_geoms, fl_segments, fl_comids, grid_intersections, grid_geoms, tol=.001)
        print("finished in {:.2f}s\n".format(time.time() - ta))

        print("computing widths...")
        m1['length'] = np.array([g.length for g in m1.geometry])
        lengths = m1[['segment', 'length']].copy()
        groups = lengths.groupby('segment')
        # compute arbolate sum at reach midpoints
        # verify that self.df and m1 are still in the same order
        assert np.sum(groups.segment.max().values - self.df.segment.values) == 0
        reach_cumsums = np.concatenate([np.cumsum(grp.length.values[::-1])[::-1] - 0.5*grp.length.values
                                      for s, grp in groups])
        segment_asums = np.array([self.df.ArbolateSu.values[s-1] for s in m1.segment.values])
        reach_asums = -1 * self.to_km * reach_cumsums + segment_asums # arbolate sums are computed in km
        # maintain positive asums; lengths in NHD often aren't exactly equal to feature lengths
        reach_asums[reach_asums < 0.] = 0
        m1['asum'] = reach_asums
        width = width_from_arbolate(reach_asums) # widths are returned in m
        if self.GISunits != 'm':
            width = width / 0.3048

        print("multiplying length units by {} to convert from GIS to MODFLOW...".format(self.mf_units_mult))
        m1['width'] = width * self.mf_units_mult
        m1['length'] = m1.length * self.mf_units_mult

        m1['roughch'] = float(roughch)
        m1['strthick'] = float(strthick)
        m1['strhc1'] = float(strhc1)

        if self.nrow is not None:
            m1['i'] = np.floor(m1.node / self.ncol).astype(int)
        if self.ncol is not None:
            j = m1.node.values % self.ncol
            m1['j'] = j
        m1['k'] = 0

        self.m1 = m1

        print("\nsetting up Mat2...")
        ta = time.time()
        self.m2 = self.df[['segment', 'outseg', 'elevup', 'elevdn']].copy()
        self.m2['comid'] = self.df.COMID
        self.m2['icalc'] = icalc
        self.m2.index = self.m2.segment
        self.m2.sort_values(by='segment', inplace=True)
        self.m2['width2'] = self.m1.groupby('segment')['width'].max().values
        self.m2['width1'] = self.m1.groupby('segment')['width'].min().values
        print("finished in {:.2f}s\n".format(time.time() - ta))

        self.m1['strtop'] = self.interpolate_to_reaches('elevup', 'elevdn')

        # add outseg information to Mat1;
        # make copy of original state prior to dropping any reaches
        self.m1.sort_values(by=['segment', 'reach'], inplace=True)
        self.m1['ReachID'] = np.arange(1, len(self.m1) + 1)
        self.set_outreaches()  # fixes any reach numbering gaps first
        graph = self.make_graph()
        self.m1['outseg'] = [graph[s] for s in self.m1.segment]
        self.m1_original_reaches = self.m1.copy()

        # discard very small reaches; redo numbering
        inds = self.m1.length > minimum_length
        print('Dropping {} reaches with length < {:.2f} ...'.format(np.sum(~inds), minimum_length))
        self.m1 = self.m1.loc[inds].copy()

        # handle co-located reaches
        if consolidate_conductance or one_reach_per_cell:
            self.m1 = consolidate_reach_conductances(self.m1, keep_only_dominant=one_reach_per_cell)

        # patch the routing
        print('\nRepairing routing connections...')
        remaining_segments = self.m1.segment.unique()

        old_graph = graph#self.make_graph()
        paths = self.paths(old_graph)
        new_graph = {}
        for k in remaining_segments:
            for s in paths[k][1:]:
                if s in remaining_segments:
                    new_graph[k] = s
                    break
                else:
                    new_graph[k] = 0

        # sync m2 with m1; update the routing
        remaining = np.array([True if s in remaining_segments else False
                              for s in self.m2.segment.values])
        self.m2 = self.m2.loc[remaining].copy()
        self.m2['outseg'] = [new_graph[s] for s in self.m2.segment]

        # renumber the segments to be consecutive
        self.renumber_segments()  # enforce best segment numbering
        self.m1.sort_values(by=['segment', 'reach'], inplace=True)
        self.m1['ReachID'] = np.arange(1, len(self.m1) + 1)
        self.set_outreaches() # fixes any reach numbering gaps first
        graph = self.make_graph()

        # add outseg information to Mat1
        self.m1['outseg'] = [graph[s] for s in self.m1.segment]

        print('\nDone creating SFR dataset.')

    def make_graph(self):
        graph = dict(zip(self.m2.segment.values, self.m2.outseg.values))
        outlets = set(graph.values()).difference(set(graph.keys()))  # including lakes
        graph.update({o: 0 for o in outlets})
        return graph

    def paths(self, graph):
        return {seg: find_path(graph, seg) for seg in graph.keys()}

    def interpolate_to_reaches(self, segvar1, segvar2):
        """Interpolate values in mat2 (segment data) to reaches in mat1 (reach data)

        Parameters
        ----------
        segvar1 : str
            Column/variable name in segment_data array for representing start of segment
            (e.g. elevup)
        segvar2 : str
            Column/variable name in segment_data array for representing start of segment
            (e.g. elevdn)

        Returns
        -------
        reach_values : 1D array
            One dimmensional array of interpolated values of same length as mat1.
            For example, elevup and elevdn could be entered as inputs to get values for the
            strtop column in reach_data.

        """
        reach_data = self.m1.sort_values(by=['segment', 'reach'])
        segment_data = self.m2.sort_values(by='segment')
        var1 = dict(zip(segment_data.segment, segment_data[segvar1]))
        var2 = dict(zip(segment_data.segment, segment_data[segvar2]))
        reach_values = []
        for seg in segment_data.segment:
            reaches = reach_data.loc[reach_data.segment == seg]
            segment_length = np.cumsum(reaches.length.values)[-1]
            # compute cumulative distance at midpoint of each reach
            dist = (np.cumsum(reaches.length) - 0.5 * reaches.length).values
            fp = [var1[seg], var2[seg]]
            xp = [0, segment_length]
            reach_values += np.interp(dist, xp, fp).tolist()
        return np.array(reach_values)

    def renumber_segments(self):
        """Renumber segments so that segment numbering is continuous and always increases
        in the downstream direction. This may speed convergence of the NWT solver 
        in some situations.
        """

        self.m2.sort_values(by='segment', inplace=True)

        nseg = self.m2.segment
        outseg = self.m2.outseg

        # explicitly fix any gaps in the numbering
        # (i.e. from removing segments)
        nseg2 = np.arange(1, len(nseg) + 1)
        # intermediate mapping that
        r1 = dict(zip(nseg, nseg2))
        r1[0] = 0
        outseg2 = np.array([r1[s] for s in outseg])

        # function re-assigning upseg numbers consecutively at one level relative to outlet(s)
        # counts down from the number of segments
        def reassign_upsegs(r, nexts, upsegs):
            nextupsegs = []
            for u in upsegs:
                r[u] = nexts if u > 0 else u  # handle lakes
                nexts -= 1
                nextupsegs += list(nseg2[outseg2 == u])
            return r, nexts, nextupsegs

        ns = len(nseg)

        # start at outlets with nss;
        # renumber upsegs consecutively at each level
        # until all headwaters have been reached
        nexts = ns
        r2 = {0: 0}
        nextupsegs = nseg2[outseg2 == 0]
        for i in range(ns):
            r2, nexts, nextupsegs = reassign_upsegs(r2, nexts, nextupsegs)
            if len(nextupsegs) == 0:
                break
        # map original segment numbers to new numbers
        r = {k: r2.get(v, v) for k, v in r1.items()}

        # renumber segments
        self.m2['segment'] = [r.get(s, s) for s in self.m2.segment]
        self.m2['outseg'] = [r.get(s, s) for s in self.m2.outseg]
        self.m2.sort_values(by='segment', inplace=True)
        inds = (self.m2.outseg > 0) & (self.m2.segment > self.m2.outseg)
        assert not np.any(inds)
        assert len(self.m2.segment) == self.m2.segment.max()

        # renumber segments in reach_data
        self.m1['segment'] = [r.get(s, s) for s in self.m1.segment]
        self.m1.sort_values(by=['segment', 'reach'], inplace=True)
        self.m1['reachID'] = np.arange(1, len(self.m1) + 1)
        self.set_outreaches()  # reset the outreaches to ensure continuity

    def reset_reaches(self):
        self.m1.sort_values(by=['segment', 'reach'], inplace=True)
        reach_data = self.m1
        segment_data = self.m2
        ireach = []
        for iseg in segment_data.segment.values:
            nreaches = np.sum(reach_data.segment == iseg)
            ireach += list(range(1, nreaches+1))
        self.m1['reach'] = ireach

    def set_outreaches(self):
        """Determine the outreach for each SFR reach (requires a reachID column in reach_data).
        Uses the segment routing specified for the first stress period to route reaches between segments.
        """
        self.m1.sort_values(by=['segment', 'reach'], inplace=True)
        self.reset_reaches() # ensure that each segment starts with reach 1
        rd = self.m1
        graph = self.make_graph()
        outseg = graph
        reach1IDs = dict(zip(rd.loc[rd.reach == 1].segment.values,
                             rd.loc[rd.reach == 1].reachID.values))
        outreach = []
        for i in range(len(rd)):
            # if at the end of reach data or current segment
            if i + 1 == len(rd) or rd.reach.values[i + 1] == 1:
                nextseg = outseg[rd.segment.values[i]]  # get next segment
                if nextseg > 0:  # current reach is not an outlet
                    nextrchid = reach1IDs[nextseg]  # get reach 1 of next segment
                else:
                    nextrchid = 0
            else:  # otherwise, it's the next reachID
                nextrchid = rd.reachID.values[i + 1]
            outreach.append(nextrchid)
        self.m1['outreach'] = outreach

    def write_tables(self, basename='SFR', m1=None, m2=None):
        """Write tables with SFR reach (Mat1) and segment (Mat2) information out to csv files.

        Parameters
        ----------
        basename: string
            e.g. Mat1 is written to <basename>Mat1.csv
        """
        if m1 is None:
            m1 = self.m1
        if m2 is None:
            m2 = self.m2
        m1_cols = self.output_cols
        m2_cols = ['segment', 'icalc', 'outseg', 'elevup', 'elevdn', 'width1', 'width2', 'comid']
        for ind in ['i', 'j']:
            if ind not in m1.columns:
                m1_cols.remove(ind)
        print("writing Mat1 to {0}{1}, Mat2 to {0}{2}".format(basename, 'Mat1.csv', 'Mat2.csv'))
        m1[m1_cols].to_csv(basename + 'Mat1.csv', index=False)
        m2[m2_cols].to_csv(basename + 'Mat2.csv', index=False)

    def write_linework_shapefile(self, basename='SFR', m1=None,
                                 output_cols=None):
        """Write a shapefile containing linework for each SFR reach,
        with segment, reach, model node number, and NHDPlus COMID attribute information

        Parameters
        ----------
        basename: string
            Output will be written to <basename>.shp
        """
        outfile = basename.replace('.shp', '') + '.shp'
        print("writing reach geometries to {}".format(outfile))
        if m1 is None:
            m1 = self.m1
        if output_cols is None:
            output_cols = self.output_cols
        if 'geometry' not in output_cols:
            output_cols += ['geometry']
        for ind in ['i', 'j']:
            if ind not in m1.columns:
                output_cols.remove(ind)
        df2shp(m1[output_cols],
               outfile, proj4=self.mf_grid_proj4)


class lines(linesBase):
    """Class for building SFR from generic GIS linework."""

    def __init__(self, lines, minElev_field=None, maxElev_field=None,
                 model_domain=None,
                 mf_grid=None, mf_grid_node_col=None, mf_units='feet',
                 lines_proj4=None,
                 routing_tol=200):

        linesBase.__init__(self, lines=lines, model_domain=model_domain,
                           mf_grid=mf_grid, mf_grid_node_col=mf_grid_node_col,
                           mf_units=mf_units,
                           lines_proj4=lines_proj4)

        self.routing_tol = routing_tol
        self.df['elevup'] = self.df[maxElev_field] if maxElev_field is not None else 0
        self.df['elevdn'] = self.df[minElev_field] if minElev_field is not None else 0
        self.allupsegs = {} # dict containing set of all upstream segments for each segment

    def append2sfr(self, sfrlinework, route2reach1=True,
                   trim_buffer=20, routing_tol=None,
                   roughch=0.037, strthick=1, strhc1=1,
                   icalc=1,
                   iupseg=0, iprior=0, nstrpts=0, flow=0, runoff=0, etsw=0, pptsw=0,
                   roughbk=0, cdepth=0, fdepth=0, awdth=0, bwdth=0):
        """Convert linework to input that can be appended to an existing SFR dataset.

        Creates Mat1 (m1) and Mat2 (m2) attributes.
        """

        if not isinstance(sfrlinework, pd.DataFrame):
            self.sfr = shp2df(sfrlinework)
        else:
            self.sfr = sfrlinework.copy()

        self.sfr.sort_values(by=['segment', 'reach'], inplace=True)
        # assign reachIDs if they don't exist
        if 'reachID' not in self.sfr.columns:
            self.sfr['reachID'] = np.arange(1, len(self.sfr) + 1)

        # set the segment numbers to being after other dataset
        starting_segment = self.sfr.segment.max() + 1
        self.df['segment'] = np.arange(len(self.df)) + starting_segment
        self.route_lines_by_proximity()
        self.route_lines_to_sfr(sfrlinework=self.sfr, route2reach1=route2reach1,
                   trim_buffer=trim_buffer, routing_tol=routing_tol)
        return self.to_sfr(starting_reachID=self.sfr.reachID.max()+1,
                    roughch=roughch, strthick=strthick, strhc1=strhc1,
                    icalc=icalc,
                    iupseg=iupseg, iprior=iprior, nstrpts=nstrpts, flow=flow, runoff=runoff, etsw=etsw, pptsw=pptsw,
                    roughbk=roughbk, cdepth=cdepth, fdepth=fdepth, awdth=awdth, bwdth=bwdth)

    def get_end_elevs_from_dem(self, dem):

        from GISio import get_values_at_points
        self.df['elevup'] = get_values_at_points(dem, self.start_cds)
        self.df['elevdn'] = get_values_at_points(dem, self.end_cds)

    def get_segment_asums(self):
        """Using allupsegs dictionary, sum lengths of all upstream segments for each segment

        Returns
        -------
        segment_asums : dict
            Dictionary of arbolate sums for each segment
        """
        if 'length' not in self.df.columns:
            self.df['length'] = [g.length for g in self.df.geometry]
        return {s: self.df.length[self.df.segment.isin(self.allupsegs[s])].sum() * self.to_km
                if s in self.allupsegs.keys() else 0
                for s in self.df.segment.tolist()}

    @property
    def start_cds(self):
        return [(g.xy[0][0], g.xy[1][0]) for g in self.df.geometry]

    @property
    def end_cds(self):
        return [(g.xy[0][-1], g.xy[1][-1]) for g in self.df.geometry]

    def renumber_segments(self):
        """Renumber segments so that segment numbering is continuous and always increases
        in the downstream direction. Experience suggests that this can substantially speed
        convergence for some models using the NWT solver."""
        r = renumber_segments(self.m2.segment.values, self.m2.outseg.values)

        self.m2['segment'] = [r.get(s, s) for s in self.m2.segment]
        self.m2['outseg'] = [r.get(s, s) for s in self.m2.outseg]
        self.m2.sort_values(by='segment', inplace=True)
        assert _in_order(self.m2.segment.values, self.m2.outseg.values)
        assert len(self.m2.segment) == self.m2.segment.max()
        self.m1['segment'] = [r.get(s, s) for s in self.m1.segment]
        self.m1['outseg'] = [r.get(s, s) for s in self.m1.outseg]

    def route_lines_by_proximity(self, tol=50):
        """Route lines based on proximity of starts and ends.

        Parameters
        ----------
        tol : numeric
            Only consider starting coordinates within tol of each end coordinate. This
            number should be fairly small, otherwise circular routing may occur.
        """
        '''
        nearest_start = get_nearest(self.start_cds, self.end_cds)

        # record the preliminary seg. number of nearest start if within tol
        segments = self.df.segment.values
        self.df['outseg'] = [segments[n] if Point(*self.start_cds[n]).distance(\
                                            Point(*self.end_cds[i])) < tol
                             else 0 for i, n in enumerate(nearest_start)]

        self.df['upsegs'] = [self.df.segment[self.df.outseg == s].tolist() for s in self.df.segment]

        self.allupsegs = get_upsegs(self.df.segment.values, self.df.outseg.values)
        '''
        nearest_start = get_nearest(self.start_cds, self.end_cds)
        next_starts = np.array(self.start_cds)[np.array(nearest_start)]
        end_crds = np.array(self.end_cds)

        # record the preliminary seg. number of nearest start if within tol
        segments = self.df.segment.values
        distances = distance(end_crds, next_starts)
        outsegs = np.zeros(len(segments))
        outsegs[distances < tol] = nearest_start[distances < tol] + 1
        self.df['outseg'] = outsegs
        self.df['upsegs'] = [self.df.segment[self.df.outseg == s].tolist() for s in self.df.segment]

        self.allupsegs = get_upsegs(self.df.segment.values, self.df.outseg.values)

        #check for circular routing (a segment shouldn't be upstream of itself)
        for k, v in self.allupsegs.items():
            assert len({k}.intersection(v)) == 0

    def route_lines_to_sfr(self, sfrlinework, route2reach1=False,
                           trim_buffer=20, routing_tol=None):
        """Route the linework in the lines class to an existing
        set of lines representing an SFR package.

        Assigns SFR outsegs and out reaches
        Modifies connecting lines so that they end at point on SFR network
        Renumbers segments in added lines so SFR segment numbers aren't duplicated

        Parameters
        ----------
        sfrlinework : str (shapefile path) or dataframe
            Contains linework representing SFR package (e.g. already broken by grid)
        route2reach1 : boolean
            If true, linework is routed to closest starting coordinate of an SFR segment.
            Otherwise, routing is to closest starting coordinate of an SFR reach
            (existing SFR segment in that location will have to be subdivided).
        """
        if not isinstance(sfrlinework, pd.DataFrame):
            self.sfr = shp2df(sfrlinework)
        else:
            self.sfr = sfrlinework.copy()

        tol = self.routing_tol if routing_tol is None else routing_tol

        if route2reach1:
            segments = self.sfr.loc[self.sfr.reach == 1, 'segment'].tolist()
            geoms = self.sfr.loc[self.sfr.reach == 1, 'geometry'].tolist()
        else:
            segments = self.sfr.segment.tolist()
            geoms = [g for g in self.sfr.geometry]

        # update segment numbering so that it starts after highest seg in sfr dataset
        if len(set(self.df.segment).intersection(self.sfr.segment)) != 0:
            maxseg = self.sfr.segment.max()
            self.df['segment'] += maxseg
            self.df.loc[self.df.outseg > 0, 'outseg'] += maxseg
            self.df['upsegs'] = [[u + maxseg for u in us] for us in self.df.upsegs.tolist()]
            self.allupsegs = get_upsegs(self.df.segment.values, self.df.outseg.values)
            # self.allupsegs = {s + maxseg: {ss + maxseg if ss > 0 else ss for ss in v}
            #                  for s, v in self.allupsegs.items()}

        sfr_start_cds = [(g.xy[0][0], g.xy[1][0]) for g in geoms]

        print('routing new lines within {} to SFR...'.format(tol))

        is_outlet = self.df.outseg.values == 0
        new_lines_outlet_cds = list(map(tuple, np.array(self.end_cds)[is_outlet]))

        # get index of nearest start to each end

        nearest_sfr = get_nearest(sfr_start_cds, new_lines_outlet_cds)
        nearest_sfr = np.array(nearest_sfr)
        next_sfr = np.array(sfr_start_cds)[nearest_sfr]
        new_lines_outlet_cds = np.array(new_lines_outlet_cds)

        # record the preliminary seg. number of nearest start if within tol
        outlets = np.zeros(len(self.df.loc[is_outlet, 'outseg']))
        distances = distance(new_lines_outlet_cds, next_sfr)
        within_tol = distances < tol
        outlets[within_tol] = nearest_sfr[within_tol]
        self.df.loc[is_outlet, 'outseg'] = outlets
        '''
        #self.df.loc[is_outlet, 'outseg'] = [segments[n]
                                       if Point(*sfr_start_cds[n]).distance(\
                                          Point(*new_lines_outlet_cds[i])) < tol
                                       else 0 for i, n in enumerate(nearest_sfr)]
        #if not route2reach1:
        #reaches = self.sfr.reach.tolist()
        '''
        reaches = self.sfr.reach.values
        outreaches = np.zeros(len(self.df.loc[is_outlet]))
        outreaches[within_tol] = reaches[within_tol]
        self.df.loc[is_outlet, 'outreachID'] = outreaches

        '''
        self.df.loc[is_outlet, 'outreachID'] = [reaches[n]
                                   if Point(*sfr_start_cds[n]).distance(\
                                      Point(*new_lines_outlet_cds[i])) < tol
                                   else 0 for i, n in enumerate(nearest_sfr)]
        '''

        def fix_newline_end(line, sfr_start_coord, trim_buffer=20):
            # trim the ends of the new lines (in case of overlap)
            # and connect them to start of closest segment or reach on SFR network
            diff = line.difference(Point(line.coords[-1]).buffer(trim_buffer))
            if diff.length > 0:
                return LineString(list(diff.coords) + [sfr_start_coord]) \
                    if sfr_start_cds not in diff.coords else LineString(list(diff.coords))
            return LineString([line.coords[0], sfr_start_coord])

        # only connect the lines that are within the routing tolerance
        geoms = self.df.geometry.values

        '''
        new_outlet_geoms = []
        for i, l in enumerate(self.df.geometry[is_outlet]):
            if within_tol[i]:
                new_geom = fix_newline_end(l, next_sfr[i], trim_buffer=trim_buffer)
                new_outlet_geoms.append(new_geom)
            else:
                new_outlet_geoms.append(i)
        geoms[is_outlet] = new_outlet_geoms

        geoms[is_outlet] = [fix_newline_end(l, sfr_start_cds[nearest_sfr[i]], trim_buffer=trim_buffer)
                            if Point(l.coords[-1]).distance(\
                               Point(sfr_start_cds[nearest_sfr[i]])) < tol
                            else l
                            for i, l in enumerate(self.df.geometry[is_outlet])]
        '''
        self.df['geometry'] = geoms


    def to_sfr(self, starting_reachID=1,
               roughch=0.037, strthick=1., strhc1=1,
               icalc=1,
               iupseg=0, iprior=0, nstrpts=0, flow=0, runoff=0, etsw=0, pptsw=0,
               roughbk=0, cdepth=0, fdepth=0, awdth=0, bwdth=0,
               tol=0.01):
        """Convert linework to SFR input.

        Creates Mat1 (m1) and Mat2 (m2) attributes.
        """

        print('\nclipping lines to active area...')
        inside = np.array([g.intersects(self.domain) for g in self.df.geometry])
        self.df = self.df.loc[inside].copy()
        if 'segment' not in self.df.columns:
            self.df['segment'] = np.arange(1, len(self.df) + 1)
        if 'outseg' not in self.df.columns:
            self.df['outseg'] = 0
        if 'upsegs' not in self.df.columns:
            self.df['upsegs'] = [[]] * len(self.df)

        line_geoms = [g.intersection(self.domain) for g in self.df.geometry]
        grid_geoms = self.grid.geometry.tolist()

        # segments may already be routed if appending to SFR
        if self.df.outseg.sum() == 0:
            print("establishing routing...")
            self.route_lines_by_proximity(tol=tol)

        print("intersecting lines with grid cells...") # this part crawls in debug mode
        grid_intersections = GISops.intersect_rtree(grid_geoms, line_geoms)

        print("setting up reaches and Mat1... (may take a few minutes for large grids)")
        ta = time.time()
        segments = self.df.segment.tolist()
        m1 = make_mat1(line_geoms, segments, segments, grid_intersections, grid_geoms, tol=tol)
        m1.sort_values(by=['segment', 'reach'], inplace=True)
        m1['reachID'] = np.arange(starting_reachID, len(m1) + starting_reachID)
        print("finished in {:.2f}s\n".format(time.time() - ta))

        print("computing lengths...")
        m1['length'] = np.array([g.length for g in m1.geometry])
        lengths = m1[['segment', 'length']].copy()
        groups = lengths.groupby('segment')

        print("computing arbolate sums at reach midpoints...")
        ta = time.time()
        reach_asums = np.concatenate([np.cumsum(grp.length.values[::-1])[::-1] - 0.5*grp.length.values
                                      for s, grp in groups])
        segment_asums_d = self.get_segment_asums()
        segment_asums = np.array([segment_asums_d[s] for s in m1.segment.values])
        reach_asums = -1 * self.to_km * reach_asums + segment_asums # arbolate sums are computed in km
        m1['asum'] = reach_asums
        print("finished in {:.2f}s\n".format(time.time() - ta))

        print("computing widths...")
        width = width_from_arbolate(reach_asums) # widths are returned in m
        if self.GISunits != 'm':
            width = width / 0.3048

        print("multiplying length units by {} to convert from GIS to MODFLOW...".format(self.mf_units_mult))
        m1['width'] = width * self.mf_units_mult
        m1['length'] = m1.length * self.mf_units_mult

        m1['roughch'] = float(roughch)
        m1['strthick'] = float(strthick)
        m1['strhc1'] = float(strhc1)
        m1['strtop'] = 0.

        if self.nrow is not None:
            m1['i'] = np.floor(m1.node / self.ncol).astype(int)
        if self.ncol is not None:
            j = m1.node.values % self.ncol
            m1['j'] = column
        m1['k'] = 0

        print("\nsetting up Mat2...")
        ta = time.time()
        write_columns = list(set(['segment', 'outseg', 'outreachID', 'elevup', 'elevdn'])\
                             .intersection(set(self.df.columns)))
        self.m2 = self.df[write_columns].copy()
        self.m2['icalc'] = icalc
        self.m2.index = self.m2.segment
        print("finished in {:.2f}s\n".format(time.time() - ta))

        self.m1 = m1
        self.m1.sort_values(by=['segment', 'reach'], inplace=True)
        self.m1['ReachID'] = np.arange(1, len(self.m1) + 1)
        self.set_outreaches()  # fixes any reach numbering gaps first
        graph = self.make_graph()
        self.m1['outseg'] = [graph[s] for s in self.m1.segment]
        self.m1_original_reaches = self.m1.copy()

        # discard very small reaches; redo numbering
        inds = self.m1.length > minimum_length
        print('Dropping {} reaches with length < {:.2f} ...'.format(np.sum(~inds), minimum_length))
        self.m1 = self.m1.loc[inds].copy()

        # handle co-located reaches
        if consolidate_conductance or one_reach_per_cell:
            self.m1 = consolidate_reach_conductances(self.m1, keep_only_dominant=one_reach_per_cell)

        # patch the routing
        print('\nRepairing routing connections...')
        remaining_segments = self.m1.segment.unique()

        old_graph = graph  # self.make_graph()
        paths = self.paths(old_graph)
        new_graph = {}
        for k in remaining_segments:
            for s in paths[k][1:]:
                if s in remaining_segments:
                    new_graph[k] = s
                    break
                else:
                    new_graph[k] = 0

        # sync m2 with m1; update the routing
        remaining = np.array([True if s in remaining_segments else False
                              for s in self.m2.segment.values])
        self.m2 = self.m2.loc[remaining].copy()
        self.m2['outseg'] = [new_graph[s] for s in self.m2.segment]

        # renumber the segments to be consecutive
        self.renumber_segments()  # enforce best segment numbering
        self.m1.sort_values(by=['segment', 'reach'], inplace=True)
        self.m1['ReachID'] = np.arange(1, len(self.m1) + 1)
        self.set_outreaches()  # fixes any reach numbering gaps first
        graph = self.make_graph()

        # add outseg information to Mat1
        self.m1['outseg'] = [graph[s] for s in self.m1.segment]

        print('\nDone creating SFR dataset.')
        '''
        # add outseg information to Mat1
        m1['outseg'] = [m2.outseg[s] for s in m1.segment]
        m1.sort_values(by=['segment', 'reach'], inplace=True)
        m1['ReachID'] = np.arange(1, len(m1) + 1)
        self.m1 = m1
        self.m2 = m2
        self.renumber_segments() # enforce best segment numbering
        print('\nDone creating SFR dataset.')
        '''
        return self.m1, self.m2


        '''
        # at this point a new 'SFR' object should be created (with its own methods such as "renumber segments")
        self.m1 = m1

        print("\nsetting up Mat2...")
        ta = time.time()
        self.m2 = self.df[['segment', 'outseg', 'elevMax', 'elevMin']].copy()
        self.m2['icalc'] = icalc
        self.renumber_segments() # enforce best segment numbering
        self.m2.index = self.m2.segment
        print("finished in {:.2f}s\n".format(time.time() - ta))

        # add outseg information to Mat1
        self.m1['outseg'] = [self.m2.outseg[s] for s in self.m1.segment]
        self.m1.sort_values(by=['segment', 'reach'], inplace=True)
        print('\nDone creating SFR dataset.')
        '''


class SFRdata(object):

    m1_columns2flopy = {'layer': 'k',
                        'row': 'i',
                        'column': 'j',
                        'segment': 'iseg',
                        'reach': 'ireach',
                        'sbtop': 'strtop',
                        'length': 'rchlen',
                        'sbthick': 'strthick',
                        'sbK': 'strhc1'
                        }
    m2_columns2flopy = {'segment': 'nseg'
                        }

    def __init__(self, reach_data, segment_data):
        pass

def _in_order(nseg, outseg):
    """Check that segment numbering increases in downstream direction.

    Parameters
    ----------
    nseg : 1-D array of segment numbers
    outseg : 1-D array of outseg numbers for segments in nseg.

    Returns
    -------
    True if there are no decreases in segment number in downstream direction, False otherwise.
    """
    inds = (outseg > 0) & (nseg > outseg)
    if not np.any(inds):
        return True
    return False

def _get_headwaters(segments, outsegs):
        """List all segments that are not outsegs (that do not have any segments upstream).

        Parameters
        ----------
        nseg : 1-D array of segment numbers
        outseg : 1-D array of outseg numbers for segments in nseg.

        Returns
        -------
        headwaters : np.ndarray (1-D)
            One dimmensional array listing all headwater segments.
        """
        upsegs = [segments[outsegs == s].tolist()
                  for s in segments]
        return segments[np.array([i for i, u in enumerate(upsegs) if len(u) == 0])]

def _get_outlets(segment_seguences_array):
    """Create a dictionary listing outlets associated with each segment
    outlet is the last value in each row of segment sequences array that is != 0 or 999999

    Parameters
    ----------
    segment_sequences_array : 2-D array produced by map_segment_sequences()

    Returns
    -------
    outlets : dict
        Dictionary of outlet number (values) for each segment (keys).
    """
    return {i + 1: r[(r != 0) & (r != 999999)][-1]
            if len(r[(r != 0) & (r != 999999)]) > 0
            else i + 1
            for i, r in enumerate(segment_seguences_array.T)}


def consolidate_reach_conductances(m1, keep_only_dominant=False):
    """For model cells with multiple SFR reaches, shift all conductance to widest reach,
    by adjusting the length, and setting the lengths in all smaller collocated reaches to 1,
    and the K-values in these to bedKmin

    Parameters
    ----------

    strhc1_min : float
        Hydraulic conductivity value to use for collocated SFR reaches in a model cell that are not the dominant
        (widest) reach. This is used to effectively set conductance in these reaches to 0, to avoid circulation
        of water between the collocated reaches.

    Returns
    -------
        Modifies the SFR Mat1 table dataframe attribute of the SFRdata object. A new SFR package file can be
        written from this table.

    Notes
    -----
        See the ConsolidateConductance notebook in the Notebooks folder.
    """
    print('\nAssigning total SFR conductance to dominant reach in cells with multiple reaches...')
    #m1['line_len'] = m1.length
    m1['cond'] = m1.width * m1.length * m1.strhc1 / m1.strthick

    # make a new column that designates whether a reach is dominant in each cell
    # dominant reaches include those not collocated with other reaches, and the longest collocated reach
    m1['Dominant'] = [True] * len(m1)

    for c, nreaches in m1.node.value_counts().iteritems():
        # this is apparently nearly twice as fast as
        # a vectorized approach on all of m1.
        # apparently because it only operates on collocated reaches
        if nreaches > 1:
            # select the collocated reaches for this cell
            df = m1[m1.node == c].sort_values(by='width', ascending=False)

            # set all of these reaches except the largest to not Dominant
            m1.loc[df.index[1:], 'Dominant'] = False

    # Sum up the conductances for all of the collocated reaches
    # returns a series of conductance sums by model cell, put these into a new column in Mat1
    Cond_sums = m1[['node', 'cond']].groupby('node').agg('sum').cond
    m1['Cond_sum'] = [Cond_sums[c] for c in m1.node]

    # Calculate a new length for widest reaches, set length in secondary collocated reaches to 1
    # also set the K values in the secondary cells to bedKmin
    #m1['length'] = 1.
    m1d = m1.loc[m1.Dominant]
    #m1.loc[m1.Dominant, 'length'] = m1d.Cond_sum * m1d.strthick / (m1d.strhc1 * m1d.width)
    #m1.loc[~m1.Dominant, 'strhc1'] = strhc1_min
    m1.loc[m1.Dominant, 'strhc1'] = m1d.Cond_sum * m1d.strthick / (m1d.length * m1d.width)
    m1.loc[~m1.Dominant, 'strhc1'] = 0.
    if keep_only_dominant:
        print('Dropping {} non-dominant reaches...'.format(np.sum(~m1.Dominant)))
        return m1.loc[m1.Dominant].copy()
    return m1

def create_reaches(part, segment_nodes, grid_geoms, tol=0.01):
    """Creates SFR reaches for a segment by ordering model cells intersected by a LineString

    Parameters
    ----------
    part: LineString
        shapely LineString object (or a part of a MultiLineString)

    segment_nodes: list of integers
        Model node numbers intersected by *part*

    grid_geoms: list of Polygons
        List of shapely Polygon objects for the model grid cells, sorted by node number

    Returns
    -------
    ordered_reach_geoms: list of LineStrings
        List of LineString objects representing the SFR reaches for the segment

    ordered_node_numbers: list of model cells containing the SFR reaches for the segment
    """
    reach_nodes = {}
    reach_geoms = {}
    # interesct flowline part with grid nodes
    reach_intersections = {c: part.intersection(grid_geoms[c]) for c in segment_nodes}
    reach_intersections = {c:g for c, g in reach_intersections.items() if g.length > 0} #drops points and empty geometries
    # empty geometries are created when segment_nodes variable includes nodes intersected by
    # other parts of a multipart line. Not sure what causes points besides duplicate vertices.

    # "flatten" all grid cell intersections to single part geometries
    n = 1
    for node, g in reach_intersections.items():
        if g.type == 'LineString':
            reach_nodes[n] = node
            reach_geoms[n] = g
            n += 1
        else:
            geoms = [gg for gg in g.geoms if gg.type == 'LineString']
            reach_nodes.update({n+nn: node for nn, gg in enumerate(geoms)})
            reach_geoms.update({n+nn: gg for nn, gg in enumerate(geoms)})
            n += len(geoms)

    # make point features for start and end of flowline part
    start = Point(part.coords[0])
    end = Point(part.coords[-1])

    ordered_reach_geoms = []
    ordered_node_numbers = []
    current_reach = start
    nreaches = len(reach_geoms) # length before entries are popped

    # for each flowline part (reach)
    for i in range(nreaches):
        '''
        # find the next flowline part (reach) that touches the current reach
        try:
            r = [j for j, g in reach_geoms.items() if current_reach.intersects(g.buffer(tol))
                ][0]
        except:
            pass
        '''
        dist = {j: g.distance(current_reach) for j, g in reach_geoms.items()}
        dist_sorted = sorted(dist.items(), key=operator.itemgetter(1))
        r = dist_sorted[0][0]
        next_reach = reach_geoms.pop(r)
        ordered_reach_geoms.append(next_reach)
        ordered_node_numbers.append(reach_nodes[r])
        current_reach = next_reach

        if current_reach.touches(end.buffer(tol)) and len(ordered_node_numbers) == nreaches:
            break
    assert len(ordered_node_numbers) == nreaches # new list of ordered node numbers must include all flowline parts
    return ordered_reach_geoms, ordered_node_numbers

def distance(c1, c2):
    return np.sqrt(np.sum((c2 - c1)**2, axis=1))

def different_projections(proj4, common_proj4):
    if not proj4 == common_proj4 \
        and not proj4 is None \
        and not common_proj4 is None:
        return True
    else:
        for prj in proj4, common_proj4:
            if prj is None:
                print("Warning, no projection information for {}!".format(prj))
        return False

def find_next(comid, pftable, comids, max_levels=10):
    """Crawls the PlusFlow table to find the next downstream comid that
    is in the set comids. Looks up subsequent downstream comids to a
    maximum number of iterations, specified by max_levels (default 10).
    """
    pftable = pftable.copy()
    nextocomid = [comid]
    comids = set(comids)
    for i in range(max_levels):
        nextocomid = pftable.loc[pftable.FROMCOMID.isin(nextocomid), 'TOCOMID'].tolist()
        if len(set(nextocomid).intersection(comids)) > 0:
            # if more than one comid is found, simply take the first
            # (often these will be in different levelpaths,
            # so there is no way to determine a preferred routing path)
            return list(set(nextocomid).intersection(comids))[0]
    return 0

def get_nearest(starts, ends):
    """Returns index of nearest start coordinate to each coordinate in ends.

    Parameters
    ----------
    starts : list of tuples
        Could be starting coordinates of each LineString.
    ends : list of tuples
        Could be ending coordinates of each LineString.
    """
    try:
        from rtree import index
    except:
        raise ImportError("This method requires the rtree package.")

    idx = index.Index()
    for i, start in enumerate(starts):
        idx.insert(i, start)

    nearest = []
    for i, end in enumerate(ends):
        results = list(idx.nearest(end, 2))
        results = [r for r in results if r != i]
        nearest.append(results[0])
    # assert that the nearest segment can't be itself
    nearest = np.array(nearest)
    # assert that the nearest segment can't be itself
    assert np.all(nearest != np.arange(len(nearest)))

    #return [list(idx.nearest(end, 2))[0] if list(idx.nearest(end, 1))[0] != i
    #        else list(idx.nearest(end, 2))[1]
    #        for i, end in enumerate(ends)]
    return nearest

def get_upsegs(nseg, outseg):
    """From segment_data, returns nested dict of containing sets of all
    segments upstream of each segment.

    Parameters
    ----------
    nseg : 1-D array of segment numbers
    outseg : 1-D array of outseg numbers for each segment in nseg.

    Returns
    -------
    upsegs : dict
        Nested dictionary of form {stress period: {segment: [set of upsegs]}}
    """
    '''
    # make a list of adjacent upsegments keyed to outseg list in Mat2
    upsegs = {o: set(nseg[outseg == o]) for o in np.unique(outseg)}

    outsegs = outseg[outseg != 0].tolist() # exclude 0, the outlet designator

    # for each outseg key, for each upseg, check for more upsegs, append until headwaters has been reached
    for s in outsegs:
        upsegslist = upsegs[s]
        for i in range(len(nseg)): # limit iterations to number of segments
            added_upsegs = set()
            [added_upsegs.update(set(upsegs[us])) for us in upsegslist if us in outsegs]
            if len(added_upsegs) == 0:
                break
            else:
                upsegslist = added_upsegs
                upsegs[s].update(added_upsegs)
    '''
    def get_nextupsegs(upsegs):
        nextupsegs = []
        for s in upsegs:
            nextupsegs += nseg[outseg == s].tolist()
        return nextupsegs

    def get_upsegs(seg):
        upsegs = nseg[outseg == seg].tolist()
        all_upsegs = upsegs
        for i in range(len(nseg)):
            upsegs = get_nextupsegs(upsegs)
            if len(upsegs) > 0:
                all_upsegs.extend(upsegs)
            else:
                break
        return all_upsegs

    return {s: get_upsegs(s) for s in nseg}

def map_segment_sequences(segments, outsegs, verbose=True):
    """Generate array containing all segment routing sequences from each headwater
    to the respective outlet.

    Parameters
    ----------
    nseg : 1-D array of segment numbers
    outseg : 1-D array of outseg numbers for segments in nseg.

    Returns
    -------
    """
    all_outsegs = np.vstack([segments, outsegs])
    nseg = len(segments)
    max_outseg = all_outsegs[-1].max()
    knt = 1
    txt = '' # text recording circular routing instances, if encountered
    while max_outseg > 0:

        nextlevel = np.array([outsegs[s - 1] if s > 0 and s < 999999 else 0
                              for s in all_outsegs[-1]])

        all_outsegs = np.vstack([all_outsegs, nextlevel])
        max_outseg = nextlevel.max()
        if max_outseg == 0:
            break
        knt += 1
        if knt > nseg:
            # subset outsegs map to only include rows with outseg number > 0 in last column
            circular_segs = all_outsegs.T[all_outsegs[-1] > 0]

            # only retain one instance of each outseg number at iteration=nss
            vals = []  # append outseg values to vals after they've appeared once
            mask = [(True, vals.append(v))[0]
                    if v not in vals
                    else False for v in circular_segs[-1]]
            circular_segs = circular_segs[:, np.array(mask)]

            # cull the circular segments array to remove duplicate instances of routing circles
            circles = []
            duplicates = []
            for i in range(np.shape(circular_segs)[0]):
                # find where values in the row equal the last value;
                # record the index of the second to last instance of last value
                repeat_start_ind = np.where(circular_segs[i] == circular_segs[i, -1])[0][-2:][0]
                # use that index to slice out the repeated segment sequence
                circular_seq = circular_segs[i, repeat_start_ind:].tolist()
                # keep track of unique sequences of repeated segments
                if set(circular_seq) not in circles:
                    circles.append(set(circular_seq))
                    duplicates.append(False)
                else:
                    duplicates.append(True)
            circular_segs = circular_segs[~np.array(duplicates), :]

            txt += '{0} instances where an outlet was not found after {1} consecutive segments!\n' \
                .format(len(circular_segs), nseg)
            txt += '\n'.join([' '.join(map(str, row)) for row in circular_segs]) + '\n'

            f = 'circular_routing.csv'
            np.savetxt(f, circular_segs, fmt='%d', delimiter=',', header=txt)
            txt += 'See {} for details.'.format(f)
            if verbose:
                print(txt)
            break

    # the array of segment sequence is useful for other other operations,
    # such as plotting elevation profiles
    return all_outsegs

def make_mat1(flowline_geoms, fl_segments, fl_comids, grid_intersections, grid_geoms, tol=0.01):

    reach = []
    segment = []
    node = []
    geometry = []
    comids = []

    for i in range(len(flowline_geoms)):
        segment_geom = flowline_geoms[i]
        segment_nodes = grid_intersections[i]
        if segment_geom.type != 'MultiLineString' and segment_geom.type != 'GeometryCollection':
            ordered_reach_geoms, ordered_node_numbers = create_reaches(segment_geom, segment_nodes, grid_geoms, tol=tol)
            reach += list(np.arange(len(ordered_reach_geoms)) + 1)
            geometry += ordered_reach_geoms
            node += ordered_node_numbers
            segment += [fl_segments[i]] * len(ordered_reach_geoms)
            comids += [fl_comids[i]] * len(ordered_reach_geoms)
        else:
            start_reach = 0
            for j, part in enumerate(list(segment_geom.geoms)):
                geoms, node_numbers = create_reaches(part, segment_nodes, grid_geoms)
                if j > 0:
                    start_reach = reach[-1]
                reach += list(np.arange(start_reach, start_reach+len(geoms)) + 1)
                geometry += geoms
                node += node_numbers
                segment += [fl_segments[i]] * len(geoms)
                comids += [fl_comids[i]] * len(geoms)
        if len(reach) != len(segment):
            print('bad reach assignment!')
            break

    m1 = pd.DataFrame({'reach': reach, 'segment': segment, 'node': node,
                            'geometry': geometry, 'comid': comids})
    m1.sort_values(by=['segment', 'reach'], inplace=True)
    m1['reachID'] = np.arange(len(m1)) + 1
    return m1

def renumber_segments(nseg, outseg):
    """Renumber segments so that segment numbering is continuous, starts at 1, and always increases
        in the downstream direction. Experience suggests that this can substantially speed
        convergence for some models using the NWT solver.

    Parameters
    ----------
    nseg : 1-D array
        Array of segment numbers
    outseg : 1-D array
        Array of outsegs for segments in nseg.

    Returns
    -------
    r : dict
        Dictionary mapping old segment numbers (keys) to new segment numbers (values). r only
        contains entries for number that were remapped.
    """
    def reassign_upsegs(r, nexts, upsegs):
        nextupsegs = []
        for u in upsegs:
            r[u] = nexts if u > 0 else u # handle lakes
            nexts -= 1
            nextupsegs += list(nseg[outseg == u])
        return r, nexts, nextupsegs

    print('enforcing best segment numbering...')
    # enforce that all outsegs not listed in nseg are converted to 0
    # but leave lakes alone
    r = {0: 0}
    r.update({o: 0 for o in outseg if o > 0 and o not in nseg})
    outseg = np.array([o if o in nseg or o < 0 else 0 for o in outseg])


    # if reach data are supplied, segment/outseg pairs may be listed more than once
    if len(nseg) != len(np.unique(nseg)):
        d = dict(zip(nseg, outseg))
        nseg, outseg = np.array(d.keys()), np.array(d.values())
    ns = len(nseg)

    nexts = ns
    nextupsegs = nseg[outseg == 0]
    for i in range(ns):
        r, nexts, nextupsegs = reassign_upsegs(r, nexts, nextupsegs)
        if len(nextupsegs) == 0:
            break
    return r

def parse_proj4_units(proj4string):
    """Determine units from proj4 string. Not tested extensively.
    """
    l = proj4string.split('+')
    units = [i for i in l if 'units' in i]
    if len(units) == 0:
        return 'latlon'
    elif 'm' in units[0]:
        return 'm'
    else:
        return 'ft'

def width_from_arbolate(arbolate, minimum_width=1):
    """Estimate stream width in feet from arbolate sum in meters, using relationship
    described by Feinstein et al (2010), Appendix 2, p 266.

    Parameters
    ----------
    arbolate: float
        Arbolate sum in km.

    Returns
    -------
    width: float
        Estimated width in meters (original formula returned width in ft)
    """
    w = 0.3048 * 0.1193 * (1000 * arbolate) ** 0.5032
    if isinstance(arbolate, np.ndarray):
        w[np.isnan(w)] = float(minimum_width)
    elif np.isnan(w):
        w = minimum_width
    else:
        pass
    return w

def find_path(graph, start, end=0, path=[]):
    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    if not isinstance(graph[start], list):
        graph[start] = [graph[start]]
    for node in graph[start]:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath: return newpath
    return None

class ProjectionError(Exception):
    def __init__(self, infile):
        self.infile = infile
    def __str__(self):
        return('\n\nModel grid shapefile is in lat-lon. Please use a projected coordinate system.')

class NodeIndexWarning(Warning):
    def __init__(self, grid_shapefile, node_field=None):
        self.grid_shapefile = grid_shapefile
        self.node_field = node_field
    def __str__(self):
        if self.node_field is None:
            return("\nWarning: Node field for {} not supplied."
                   "Node numbers will be assigned sequentially to grid cells starting at 0."
                   "This will result in incorrect location of SFR reaches if the grid shapefile"
                   "is not sorted.".format(self.node_field))
        else:
            return("\nWarning: Node numbers in field {} of grid are not consecutive intergers"
                    "starting at 0. SFRmaker will sort the grid cells by node, but then use"
                    "base-0 consecutive integers for indexing. These will be"
                    "the node numbers reported in the output (sfr package, tables, and shapefile). This"
                    "may lead to confusion if the SFR output is joined back to {}."
                   .format(self.node_field, self.grid_shapefile))
