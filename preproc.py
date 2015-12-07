__author__ = 'aleaf'

import time
import numpy as np
import pandas as pd
import fiona
from shapely.geometry import Point, Polygon, shape
from shapely.ops import unary_union
from GISio import shp2df, df2shp, get_proj4
from GISops import project, projectdf, build_rtree_index, intersect_rtree
import GISops

class NHDdata(object):

    convert_elevslope_to_model_units = {'feet': 0.0328084,
                                        'meters': 0.01}

    def __init__(self, NHDFlowline=None, PlusFlowlineVAA=None, PlusFlow=None, NHDFcode=None,
                 elevslope=None,
                 mf_grid=None, mf_grid_node_col=None,
                 nrows=None, ncols=None,
                 mfdis=None, xul=None, yul=None, rot=0,
                 model_domain=None,
                 flowlines_proj4=None, mfgrid_proj4=None, domain_proj4=None,
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
        flowlines_proj4 : str, optional
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
        self.nrows = nrows
        self.ncols = ncols
        self.mfdis = mfdis
        self.xul = xul
        self.yul = yul
        self.rot = rot

        # unit conversions (set below after grid projection is verified)
        self.mf_units = mf_units
        self.mf_units_mult = 1.0 # go from GIS units to model units
        self.GISunits = None #
        self.to_km = None # converts GIS units to km for arbolate sum

        self.fl_proj4 = flowlines_proj4
        self.mf_grid_proj4 = mfgrid_proj4
        self.domain_proj4 = domain_proj4

        print("Reading input...")
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
            self.grid.index = self.grid[mf_grid_node_col].values
        else:
            print('Warning: Node field for grid shape file not supplied. \
                  Node numbers will be assigned using index. \
                  This may result in incorrect location of SFR reaches.')
        self.grid = self.grid[['geometry']]

        # get projections
        if self.mf_grid_proj4 is None and not isinstance(mf_grid, pd.DataFrame):
            self.mf_grid_proj4 = get_proj4(mf_grid)
        if self.fl_proj4 is None:
            if isinstance(NHDFlowline, list):
                self.fl_proj4 = get_proj4(NHDFlowline[0])
            elif not isinstance(NHDFlowline, pd.DataFrame):
                self.fl_proj4 = get_proj4(NHDFlowline)

        # set the indices
        for attr, index in {'fl': 'COMID',
                            'pfvaa': 'ComID',
                            'elevs': 'COMID'}.items():
            if not self.__dict__[attr].index.name == index:
                self.__dict__[attr].index = self.__dict__[attr][index]

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

        # convert the elevations from elevslope table
        self.elevs['Max'] = self.elevs.MAXELEVSMO * self.convert_elevslope_to_model_units[self.mf_units]
        self.elevs['Min'] = self.elevs.MINELEVSMO * self.convert_elevslope_to_model_units[self.mf_units]

        if different_projections(self.fl_proj4, self.mf_grid_proj4):
            print("reprojecting NHDFlowlines from\n{}\nto\n{}...".format(self.fl_proj4, self.mf_grid_proj4))
            self.fl['geometry'] = projectdf(self.fl, self.fl_proj4, self.mf_grid_proj4)

        if model_domain is not None \
                and different_projections(self.domain_proj4, self.mf_grid_proj4):
            print("reprojecting model domain from\n{}\nto\n{}...".format(self.domain_proj4, self.mf_grid_proj4))
            self.domain = project(self.domain, self.domain_proj4, self.mf_grid_proj4)

    def list_updown_comids(self):
        print('getting routing information from NHDPlus Plusflow table...')
        # setup local variables and cull plusflow table to comids in model
        comids = self.df.index.tolist()
        pf = self.pf.ix[(self.pf.FROMCOMID.isin(comids)) |
                        (self.pf.TOCOMID.isin(comids))].copy()

        # subset PlusFlow entries for comids that are not in flowlines dataset
        # comids may be missing because they are outside of the model
        # or if the flowlines dataset was edited (resulting in breaks in the routing)
        missing_tocomids = ~pf.TOCOMID.isin(comids) & (pf.TOCOMID != 0)
        missing = pf.ix[missing_tocomids, ['FROMCOMID', 'TOCOMID']].copy()
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

        # reduce dncomids to 1 per segment
        braids = self.df[np.array([len(d) for d in self.df.dncomids]) > 1]
        for i, r in braids.iterrows():
            # select the dncomid that has a matching levelpath
            levelpath_matches = self.df.ix[r.dncomids, 'LevelPathI'].values == r.LevelPathI
            in_same_levelpath = np.array(r.dncomids)[levelpath_matches]
            # if none match, select the first dncomid
            if len(in_same_levelpath) == 0:
                dncomid = [r.dncomids[0]]
            else:
                dncomid = np.unique(in_same_levelpath).tolist()
            self.df.set_value(i, 'dncomids', dncomid)

        # assign upsegs and outsegs based on NHDPlus routing
        self.df['upsegs'] = [[self.df.segment[c] if c != 0 else 0 for c in comids] for comids in self.df.upcomids]
        self.df['dnsegs'] = [[self.df.segment[c] if c != 0 else 0 for c in comids] for comids in self.df.dncomids]

        # make a column of outseg integers
        self.df['outseg'] = [d[0] for d in self.df.dnsegs]

    def to_sfr(self, roughness=0.037, streambed_thickness=1, streambedK=1,
               icalc=1,
               iupseg=0, iprior=0, nstrpts=0, flow=0, runoff=0, etsw=0, pptsw=0,
               roughch=0, roughbk=0, cdepth=0, fdepth=0, awdth=0, bwdth=0):


        # create a working dataframe
        self.df = self.fl[self.fl_cols].join(self.pfvaa[self.pfvaa_cols], how='inner')

        # bring in elevations from elevslope table
        self.df = self.df.join(self.elevs[['Max', 'Min']], how='inner')

        print('\nclipping flowlines to active area...')
        inside = np.array([g.intersects(self.domain) for g in self.df.geometry])
        self.df = self.df.ix[inside].copy()
        self.df.sort_values(by='COMID', inplace=True)
        flowline_geoms = [g.intersection(self.domain) for g in self.df.geometry]
        grid_geoms = self.grid.geometry.tolist()

        print("intersecting flowlines with grid cells...") # this part crawls in debug mode
        grid_intersections = GISops.intersect_rtree(grid_geoms, flowline_geoms)

        print("setting up segments... (may take a few minutes for large networks)")
        ta = time.time()
        self.list_updown_comids()
        self.assign_segments()
        fl_segments = self.df.segment.tolist()
        fl_comids = self.df.COMID.tolist()
        print("finished in {:.2f}s\n".format(time.time() - ta))

        print("setting up reaches and Mat1... (may take a few minutes for large grids)")
        ta = time.time()
        m1 = make_mat1(flowline_geoms, fl_segments, fl_comids, grid_intersections, grid_geoms)
        print("finished in {:.2f}s\n".format(time.time() - ta))

        print("computing widths...")
        m1['length'] = np.array([g.length for g in m1.geometry])
        lengths = m1[['segment', 'length']].copy()
        groups = lengths.groupby('segment')
        # compute arbolate sum at reach midpoints
        reach_asums = np.concatenate([np.cumsum(grp.length.values[::-1])[::-1] - 0.5*grp.length.values
                                      for s, grp in groups])
        segment_asums = np.array([self.df.ArbolateSu.values[s-1] for s in m1.segment.values])
        reach_asums = -1 * self.to_km * reach_asums + segment_asums # arbolate sums are computed in km
        width = width_from_arbolate(reach_asums) # widths are returned in m
        if self.GISunits != 'm':
            width = width / 0.3048

        print("multiplying length units by {} to convert from GIS to MODFLOW...".format(self.mf_units_mult))
        m1['width'] = width * self.mf_units_mult
        m1['length'] = m1.length * self.mf_units_mult

        m1['roughness'] = roughness
        m1['sbthick'] = streambed_thickness
        m1['sbK'] = streambedK
        m1['sbtop'] = 0

        if self.nrows is not None:
            m1['row'] = np.floor(m1.node / self.ncols) + 1
        if self.ncols is not None:
            column = m1.node.values % self.ncols
            column[column == 0] = self.ncols # last column has remainder of 0
            m1['column'] = column
        m1['layer'] = 1

        self.m1 = m1

        print("\nsetting up Mat2...")
        ta = time.time()
        self.m2 = self.df[['segment', 'outseg', 'Max', 'Min']].copy()
        self.m2['icalc'] = icalc
        self.renumber_segments() # enforce best segment numbering
        self.m2.index = self.m2.segment
        print("finished in {:.2f}s\n".format(time.time() - ta))

        # add outseg information to Mat1
        self.m1['outseg'] = [self.m2.outseg[s] for s in self.m1.segment]
        print('\nDone creating SFR dataset.')

    def renumber_segments(self):
        """Renumber segments so that segment numbering is continuous and always increases
        in the downstream direction. Experience suggests that this can substantially speed
        convergence for some models using the NWT solver."""
        r = renumber_segments(self.m2.segment.values, self.m2.outseg.values)

        self.m2.segment = [r[s] for s in self.m2.segment]
        self.m2.outseg = [r[s] for s in self.m2.outseg]
        self.m2.sort_values(by='segment', inplace=True)
        assert _in_order(self.m2.segment.values, self.m2.outseg.values)
        assert len(self.m2.segment) == self.m2.segment.max()
        self.m1.segment = [r[s] for s in self.m1.segment]

    def write_tables(self, basename='SFR'):
        """Write tables with SFR reach (Mat1) and segment (Mat2) information out to csv files.

        Parameters
        ----------
        basename: string
            e.g. Mat1 is written to <basename>Mat1.csv
        """
        m1_cols = ['node', 'layer', 'segment', 'reach', 'sbtop', 'width', 'length', 'sbthick', 'sbK', 'roughness', 'reachID']
        m2_cols = ['segment', 'icalc', 'outseg', 'Max', 'Min']
        if self.nrows is not None:
            m1_cols.insert(1, 'row')

        if self.ncols is not None:
            m1_cols.insert(2, 'column')
        print("writing Mat1 to {0}{1}, Mat2 to {0}{2}".format(basename, 'Mat1.csv', 'Mat2.csv'))
        self.m1[m1_cols].to_csv(basename + 'Mat1.csv', index=False)
        self.m2[m2_cols].to_csv(basename + 'Mat2.csv', index=False)

    def write_linework_shapefile(self, basename='SFR'):
        """Write a shapefile containing linework for each SFR reach,
        with segment, reach, model node number, and NHDPlus COMID attribute information

        Parameters
        ----------
        basename: string
            Output will be written to <basename>.shp
        """
        print("writing reach geometries to {}".format(basename+'.shp'))
        df2shp(self.m1[['reachID', 'node', 'segment', 'reach', 'outseg', 'comid', 'geometry']],
               basename+'.shp', proj4=self.mf_grid_proj4)


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

def create_reaches(part, segment_nodes, grid_geoms):
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
    reach_intersections = [part.intersection(grid_geoms[c]) for c in segment_nodes]

    # "flatten" all grid cell intersections to single part geometries
    n = 1
    for i, g in enumerate(reach_intersections):
        if g.type == 'LineString':
            reach_nodes[n] = segment_nodes[i]
            reach_geoms[n] = g
            n += 1
        else:
            reach_nodes.update({n+nn: segment_nodes[i] for nn, gg in enumerate(g.geoms)})
            reach_geoms.update({n+nn: gg for nn, gg in enumerate(g.geoms)})
            n += len(g.geoms)

    # make point features for start and end of flowline part
    partlist = list(zip(part.xy[0], part.xy[1]))
    start = Point(partlist[0])
    end = Point(partlist[-1])

    ordered_reach_geoms = []
    ordered_node_numbers = []
    nreaches = len(reach_geoms)
    current_reach = start

    # for each flowline part (reach)
    for i in range(nreaches):
        # find the next flowline part (reach) that touches the current reach
        r = [j for j, g in reach_geoms.items() if current_reach.intersects(g.buffer(0.001))
             ][0]
        next_reach = reach_geoms.pop(r)
        ordered_reach_geoms.append(next_reach)
        ordered_node_numbers.append(reach_nodes[r] + 1)
        current_reach = next_reach

        if current_reach.touches(end.buffer(0.001)) and len(ordered_node_numbers) == nreaches:
            break
    assert len(ordered_node_numbers) == nreaches # new list of ordered node numbers must include all flowline parts
    return ordered_reach_geoms, ordered_node_numbers

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
        nextocomid = pftable.ix[pftable.FROMCOMID.isin(nextocomid), 'TOCOMID'].tolist()
        if len(set(nextocomid).intersection(comids)) > 0:
            # if more than one comid is found, simply take the first
            # (often these will be in different levelpaths,
            # so there is no way to determine a preferred routing path)
            return list(set(nextocomid).intersection(comids))[0]
    return 0

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

def make_mat1(flowline_geoms, fl_segments, fl_comids, grid_intersections, grid_geoms):

    reach = []
    segment = []
    node = []
    geometry = []
    comids = []

    for i in range(len(flowline_geoms)):
        segment_geom = flowline_geoms[i]
        segment_nodes = grid_intersections[i]

        if segment_geom.type != 'MultiLineString':
            ordered_reach_geoms, ordered_node_numbers = create_reaches(segment_geom, segment_nodes, grid_geoms)
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
    """Renumber segments so that segment numbering is continuous and always increases
        in the downstream direction. Experience suggests that this can substantially speed
        convergence for some models using the NWT solver.

    Parameters
    ----------
    nseg : 1-D array
        Array of segment numbers
    outseg : 1-D array
        Array of outsegs for segments in neg.

    Returns
    -------
    r : dict
        Dictionary mapping old segment numbers (keys) to new segment numbers (values)
    """
    def reassign_upsegs(r, nexts, upsegs):
        nextupsegs = []
        for u in upsegs:
            r[u] = nexts if u > 0 else u # handle lakes
            nexts -= 1
            nextupsegs += list(nseg[outseg == u])
        return r, nexts, nextupsegs

    print('enforcing best segment numbering...')
    ns = len(nseg)

    nexts = ns
    r = {0: 0}
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


class ProjectionError(Exception):
    def __init__(self, infile):
        self.infile = infile
    def __str__(self):
        return('\n\nModel grid shapefile is in lat-lon. Please use a projected coordinate system.')