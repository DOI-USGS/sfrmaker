__author__ = 'aleaf'

import os
import numpy as np
from postproc import SFRdata



class diagnostics(SFRdata):
    """
    Daniel Feinstein's top 10 SFR problems:
    1) cell gaps btw adjacent reaches in a single segment
    2) cell gaps btw routed segments. possibly because of re-entry problems at domain edge
    3) adjacent reaches with STOP sloping the wrong way
    4) routed segments with end/start sloping the wrong way
    5) STOP>TOP1 violations, i.e.,floaters
    6) STOP<<TOP1 violations, i.e., exaggerated incisions
    7) segments that end within one diagonal cell distance from another segment, inviting linkage
    8) circular routing of segments
    9) multiple reaches with non-zero conductance in a single cell
    10) reaches in inactive cells
    """

    def __init__(self, sfrobject=None, Mat1=None, Mat2=None, sfr=None, node_column=None,
                     mfpath=None, mfnam=None, mfdis=None,
                     xll=0, yll=0, outpath=os.getcwd()):

        SFRdata.__init__(self, sfrobject=sfrobject, Mat1=Mat1, Mat2=Mat2, sfr=sfr, node_column=node_column,
                         mfpath=mfpath, mfnam=mfnam, mfdis=mfdis, xll=xll, yll=yll)

    def check_numbering(self):
        """checks for continuity in segment and reach numbering
        """
        print 'Checking for continuity in segment and reach numbering...'
        self.m1.sort(['segment', 'reach'], inplace=True)
        self.m2.sort('segment', inplace=True)
        uniquesegs = self.m2.segment.values
        if np.max(np.diff(uniquesegs)) != 1 or np.min(np.diff(uniquesegs)) != 1:
            print 'gaps in segment numbering!, check Mat2'
            return
        segments, reaches = self.m1.segment.values, self.m1.reach.values

        for seg in uniquesegs:
            r = reaches[segments == seg]
            if len(r) > 1:
                if np.max(np.diff(r)) != 1 or np.min(np.diff(r)) != 1:
                    print 'Invalid reach numbering at segment {:.0f}!'.format(seg)
                    print self.m1.reach[self.m1.segment == seg]
                    break
        print 'passed.'

    def check_routing(self, max_levels=1000):
        """checks for breaks in routing and does comprehensive check for circular routing
        """
        print '\nChecking for circular routing...'
        circular_routing = self.map_outsegs(max_levels=max_levels)
        if circular_routing is not None:
            print circular_routing
        else:
            print 'passed.'

    def check_overlapping(self, tol=1e-6):
        """checks for multiple SFR reaches in one cell; and whether more than one reach has Cond > 0
        """
        print '\nChecking for model cells with multiple non-zero SFR conductances...'

        m1 = self.m1.copy()
        m1nodes = m1.node.values
        # Calculate SFR conductance for each reach
        def cond(X):
            c = X['sbK'] * X['width'] * X['SFRlength'] / X['sbthick']
            return c

        m1['Cond'] = m1.apply(cond, axis=1)

        nodes_with_multiple_conductance = []
        for c in self.shared_cells:

            # select the collocated reaches for this cell
            df = m1[m1nodes == c].sort('Cond', ascending=False)
            cnd = df.Cond.values

            if cnd[1] / cnd[0] > tol:
                nodes_with_multiple_conductance.append(c)

        if len(nodes_with_multiple_conductance) > 0:
            print '{} model cells with multiple non-zero SFR conductances found. ' \
                  'This may lead to circular routing between collocated reaches. ' \
                  '\nSee collocated_reaches.csv.' \
                  '\nRun consolidate_conductance() to establish a dominant reach in each SFR cell ' \
                  '(sets remaining reaches to zero-conductance).'.format(len(nodes_with_multiple_conductance))

            df = m1.ix[np.in1d(m1nodes, nodes_with_multiple_conductance),
                       ['node', 'segment', 'reach', 'length', 'SFRlength', 'width', 'sbK', 'Cond']]
            df.sort('node', inplace=True)
            df.to_csv('collocated_reaches.csv', index=False)
        else:
            print 'passed.'

    def check_elevations(self):
        """checks that elevations decrease monotonically downstream
        """
        outfile = 'elevation_increases.csv'
        print '\nChecking Mat1 for downstream rises in streambed elevation...'
        diffs = self.m1.sbtop.diff()
        pos_diffs = self.m1.ix[diffs > 0, ['segment', 'reach', 'sbtop']]
        pos_diffs['sb_rise'] = diffs[diffs > 0]

        # reach 1s may have positive differences due to other segment in previous index
        if pos_diffs.reach.max() > 1:
            bad_elevs = pos_diffs[pos_diffs.reach > 1]
            print 'Found {} instances of bad elevations, see {}'.format(len(bad_elevs), outfile)
            bad_elevs.sort('sb_rise', ascending=False, inplace=True)
            bad_elevs.to_csv(os.path.join(self.outpath, outfile))
        else:
            print 'passed.'

        print '\nChecking Mat1 for segments with reach 1 higher than last reach ...'
        sbtop = self.m1.sbtop.values
        m1segments = self.m1.segment.values
        m1reaches = self.m1.reach.values
        m2segments = self.m2.segment.values
        r1elev = np.array([sbtop[(m1segments == s) & (m1reaches == 1)][0] for s in m2segments])
        lr = [np.max(m1reaches[(m1segments == s)]) for s in m2segments]
        lrelev = np.array([sbtop[(m1segments == s) & (m1reaches == lr[i])][0] for i, s in enumerate(m2segments)])
        diffs = m2segments[(lrelev - r1elev) > 0]
        if len(diffs) > 0:
            print 'Found {} segments with lower reach 1, see Mat1_backwards_segments.csv'.format(len(diffs))
            df = self.m1.ix[self.m1.segment.isin(diffs)]
            df.sort(['segment', 'reach'], inplace=True)
            df.to_csv('Mat1_backwards_segments.csv', index=False)
        else:
            print 'passed.'

        '''
        Need to add check for sb elevations above cell tops! (floating reaches)
        '''
        print '\nChecking for MODFLOW altitude errors...'
        if self.elevs is None:
            print 'Model elevations not found. Please run read_dis2() with the model DIS file.'
            return
        rows = self.m1.row.values - 1
        columns = self.m1.column.values - 1
        layers = self.m1.layer.values - 1
        bots = self.elevs[layers+1, rows, columns]
        sb = self.m1.sbtop.values - self.m1.sbthick.values
        altitude_errors = ~np.array([sb[i] > bots[i] for i in np.arange(len(sb))])

        outfile = os.path.join(self.outpath, 'altitude_errors.csv')
        if len(self.m1[altitude_errors]) > 0:
            print 'Altitude errors found. Run either' \
                  '\n1) reset_model_top_2streambed(), to put all SFR cells in layer 1, or' \
                  '\n2) assign_layers() to put each SFR cell in the model layer bracketing streambed bottom elevations' \
                  '\nsee {} for more info'.format(outfile)
            df = self.m1[altitude_errors].copy()
            df['layer_bottom'] = bots[altitude_errors]
            df.to_csv(outfile)
        else:
            print 'passed.'

    def check_Mat2_min_max(self):

        print '\nChecking for segment minimum elevations > segment maximum elevations in Mat2...'
        # compute differences in min and max elevations listed in Mat2
        diffs = np.ravel(np.diff(self.m2[['Min', 'Max']].values))
        df = self.m2.ix[diffs < 0, ['segment', 'Min', 'Max']]

        if len(df) > 0:
            df['rise'] = df.Min - df.Max
            df.sort('rise', ascending=False, inplace=True)
            df.to_csv('Mat2_elevation_violations.csv', index=False)
            print '{} segments found with Min elevations > Max. See Mat2_elevation_violations.csv'.format(len(df))
            return
        print 'passed.'

    def check_minimum_slope(self):

        print '\nChecking to make sure that minimum slope was honored...'
        if self.m1.slope.min() != self.minimum_slope:
            loc = self.m1.iloc[np.argmin(self.m1.slope.values)][['segment', 'reach']]
            print '{:.1e} was specified as a minimum slope, but slope of {:.1e} found at ' \
                  'segment {}, reach {}'.format(loc.segment, loc.reach)
            return
        else:
            print 'passed.'

    def check_outlets(self, model_domain=None, buffer=100):

        print "\nChecking for outlets in the model interior..."
        if model_domain is None:
            print 'Need a shapefile of the model domain edge to check for interior outlets.\n' \
                  '(rotated coordinate systems not supported; ' \
                  'for rotated grids submit shapefile of unrotated domain in coordinates consistent with model cells.)'
            return
        else:
            import GISio

            df = GISio.shp2df(model_domain)
            self.domain = df.iloc[0].geometry

        self._get_sfr_cell_geoms()

        print '\nChecking for breaks in routing (outlets) within the SFR network'
        if 'Outlet' not in self.m1.columns:
            circular_routing = self.map_outsegs()
            if circular_routing is not None:
                print circular_routing
                print '\nCannot evaluate interior outlets until circular routing is fixed.'
                return
            else:
                pass

        outlets = np.unique(self.m1.Outlet.values)
        self.m1.sort(['segment', 'reach'], inplace=True)
        outlet_nodes = [self.m1.ix[(self.m1.segment == o), 'node'].values[-1] for o in outlets]
        outlet_geoms = [self.m1.ix[(self.m1.segment == o), 'geometry'].values[-1] for o in outlets]

        interior_outlets = [n for i, n in enumerate(outlet_nodes)
                            if outlet_geoms[i].buffer(buffer).within(self.domain)]

        if len(interior_outlets) > 0:
            print 'Interior outlets found at the following nodes:\n'
            for i in interior_outlets:
                print '{} '.format(i),
            print '\n'
        else:
            # make sure that at least 1 SFR cell is inside the domain
            for g in self.m1.geometry:
                if g.within(self.domain):
                    print 'passed.'
                    break
                else:
                    continue
            print "No SFR cells were inside of the supplied domain! Check domain shapefile coordinates,\n" \
                  "and that the correct model origin was supplied."

    def check_4gaps_in_routing(self, model_domain=None, tol=0):

        print "\nChecking for gaps in routing between segments..."
        if model_domain is None:
            print 'No model_domain supplied. ' \
                  'Routing gaps for segments intersecting model domain boundary will not be considered.\n' \
                  '(rotated coordinate systems not supported; ' \
                  'for rotated grids submit shapefile of unrotated domain in coordinates consistent with model cells.)'
        else:
            import GISio

            df = GISio.shp2df(model_domain)
            self.domain = df.iloc[0].geometry

        self._get_sfr_cell_geoms()

        m1 = self.m1.copy()
        m2 = self.m2.copy()

        if tol is None:
            try:
                tol = np.min([np.min(self.dis.delr),
                              np.min(self.dis.delc)])
            except:
                tol = 0
        # get number of reaches for each segment
        max_reach_numbers = [np.max(m1.reach[m1.segment == s]) for s in m2.segment]

        # get cell centroids for last reach in each segment
        end_centroids = [m1.geometry[(m1.segment == s) &
                                         (m1.reach == max_reach_numbers[i])].values[0].centroid
                         for i, s in enumerate(m2.segment.tolist())]

        # get cell centroids for first reach in each segment
        start_centroids = [m1.geometry[(m1.segment == s) &
                                           (m1.reach == 1)].values[0].centroid
                           for s in m2.segment]

        # compute distances between end reach cell centroid, and cell centroid of outseg reach 1
        distances = [end_centroids[i].distance(start_centroids[os-1]) if 0 < os < 999999 else 0
             for i, os in enumerate(m2.outseg.astype(int))]

        m2['routing_distance'] = distances
        m2['end_reach_geom'] = [m1.geometry[(m1.segment == s) &
                                            (m1.reach == max_reach_numbers[i])].values[0]
                                for i, s in enumerate(m2.segment.tolist())]

        routing = m2.ix[m2.routing_distance > tol, ['segment', 'outseg', 'routing_distance', 'end_reach_geom']]\
            .sort('routing_distance', ascending=False)

        reportfile = 'routing_gaps.csv'
        if model_domain is not None:

            # identify routing gaps that do not occur along the model boundary
            # (HWR code in sfr_classes may route segments along the model boundary to each other)
            interior_gaps = [g.within(self.domain) for g in routing.end_reach_geom]
            routing = routing[interior_gaps]

            if len(routing) > 0:
                print '{:.0f} gaps in routing greater than {} (default is minimum model cell size)' \
                      '\nfound that do not coincide with {}. See {}'.format(len(routing),
                                                                                                                 tol,
                                                                                                                 model_domain,
                                                                                                                 reportfile)
                routing.drop('end_reach_geom', axis=1).to_csv(reportfile, index=False)
                return

        if len(routing) > 0:
            print '{:.0f} gaps in routing greater than {} found, ' \
                  'but these may coincide with model domain boundary. See {}'.format(len(routing), tol, reportfile)
            routing.drop('end_reach_geom', axis=1).to_csv(reportfile, index=False)
            return

        print 'passed.'

    def check_grid_intersection(self, sfr_linework_shapefile=None):

        print '\Checking SFR cell geometries for consistancy with linework...'
        if sfr_linework_shapefile is None:
            print 'Need linework to check SFR cell geometries.'
            return
        else:
            import GISio

            df = GISio.shp2df(sfr_linework_shapefile)
            try:
                self.lines = df[['segment', 'reach', 'geometry']].copy()
                self.lines.sort(['segment', 'reach'], inplace=True)
            except IndexError:
                print 'Linework shapfile must have segment and reach information!'
        self._get_sfr_cell_geoms()
        self.m1.sort(['segment', 'reach'], inplace=True)
        lines, cells = self.lines.geometry.tolist(), self.m1.geometry.tolist()
        intersections = np.array([lines[i].intersects(cell) for i, cell in enumerate(cells)])

        nmismatch = len(self.m1[~intersections])
        if nmismatch > 0:
            reportfile = 'bad_intersections.csv'
            print "{} SFR reaches don't coincide with linework in {}!" \
                  "see {}.".format(nmismatch, sfr_linework_shapefile, reportfile)
            self.m1[~intersections].to_csv(reportfile, index=False)
        else:
            print 'passed.'

    def plot_segment_linkages(self, linkshp='segment_linkages.shp', outletshp='outlets.shp'):

        from shapely.geometry import LineString
        from GISio import df2shp

        print '\nMaking shapefiles of SFR segment linkages and outlet locations...'
        self._get_sfr_cell_geoms()

        segments = self.m1[['segment', 'geometry']].groupby('segment')
        geoms = [(df.geometry.tolist()[0].centroid,
                  df.geometry.tolist()[-1].centroid) for s, df in segments]

        linksdf = self.m2[['segment', 'outseg']].copy()
        linksdf['geometry'] = [LineString([geoms[s][1], geoms[o-1][0]])
                               for s, o in enumerate(linksdf.outseg)]
        linksdf = linksdf[linksdf.outseg != 0]
        df2shp(linksdf, linkshp, prj=self.prj)
        df2shp(linksdf, outletshp, prj=self.prj)

    def _get_sfr_cell_geoms(self):

        if 'geometry' not in self.m1.columns:
            try:
                self.get_cell_geometries()
            except:
                print "No geometry column in attribute m1; couldn't generate geometries from dis attribute.\n" \
                      "Read in the DIS file by running read_dis2()."

        if self.xll == 0 and self.yll == 0:
            print 'Warning, model origin for SFR object is 0, 0.'