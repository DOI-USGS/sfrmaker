__author__ = 'aleaf'

import os
import numpy as np
import pandas as pd
from postproc import SFRdata



class diagnostics(SFRdata):

    def __init__(self, sfrobject=None, Mat1=None, Mat2=None, sfr=None, node_column=None,
                     mfpath=None, mfnam=None, mfdis=None, outpath=os.getcwd()):

        SFRdata.__init__(self, sfrobject=sfrobject, Mat1=Mat1, Mat2=Mat2, sfr=sfr, node_column=node_column,
                         mfpath=mfpath, mfnam=mfnam, mfdis=mfdis)

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

    def check_routing(self):
        """checks for breaks in routing and does comprehensive check for circular routing
        """
        print '\nChecking for circular routing...'
        circular_routing = self.map_outsegs()
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
            c = X['bed_K'] * X['width_in_cell'] * X['length_in_cell'] / X['bed_thickness']
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
                       ['node', 'segment', 'reach', 'length_in_cell', 'width_in_cell', 'bed_K', 'Cond']]
            df.sort('node', inplace=True)
            df.to_csv('collocated_reaches.csv', index=False)
        else:
            print 'passed.'

    def check_elevations(self):
        """checks that elevations decrease monotonically downstream
        """
        outfile = 'elevation_increases.csv'
        print '\nChecking Mat1 for downstream rises in streambed elevation...'
        diffs = self.m1.top_streambed.diff()
        pos_diffs = self.m1.ix[diffs > 0, ['segment', 'reach', 'top_streambed']]
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
        top_streambed = self.m1.top_streambed.values
        m1segments = self.m1.segment.values
        m1reaches = self.m1.reach.values
        m2segments = self.m2.segment.values
        r1elev = np.array([top_streambed[(m1segments == s) & (m1reaches == 1)][0] for s in m2segments])
        lr = [np.max(m1reaches[(m1segments == s)]) for s in m2segments]
        lrelev = np.array([top_streambed[(m1segments == s) & (m1reaches == lr[i])][0] for i, s in enumerate(m2segments)])
        diffs = m2segments[(lrelev - r1elev) > 0]
        if len(diffs) > 0:
            print 'Found {} segments with lower reach 1, see Mat1_backwards_segments.csv'.format(len(diffs))
            df = self.m1.ix[self.m1.segment.isin(diffs)]
            df.sort(['segment', 'reach'], inplace=True)
            df.to_csv('Mat1_backwards_segments.csv', index=False)
        else:
            print 'passed.'


        print '\nChecking for MODFLOW altitude errors...'
        if self.elevs is None:
            print 'Model elevations not found. Please run read_dis2() with the model DIS file.'
            return
        rows = self.m1.row.values - 1
        columns = self.m1.column.values - 1
        layers = self.m1.layer.values - 1
        bots = self.elevs[layers+1, rows, columns]
        sb = self.m1.top_streambed.values - self.m1.bed_thickness.values
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
        if self.m1.bed_slope.min() != self.minimum_slope:
            loc = self.m1.iloc[np.argmin(self.m1.bed_slope.values)][['segment', 'reach']]
            print '{:.1e} was specified as a minimum slope, but slope of {:.1e} found at ' \
                  'segment {}, reach {}'.format(loc.segment, loc.reach)
            return
        else:
            print 'passed.'

    def check_outlets(self, model_domain=None, buffer=100):

        if model_domain is None:
            print 'Need a shapefile of the model domain edge to check for interior outlets.'
            return
        else:
            import GISio

            df = GISio.shp2df(model_domain)
            self.domain = df.iloc[0].geometry

        if 'geometry' not in self.m1.columns:
            try:
                self.get_cell_geometries()
            except:
                print "No geometry column in attribute m1; couldn't generate geometries from dis attribute.\n" \
                      "Read in the DIS file by running read_dis2()."

        if self.xll == 0 and self.yll == 0:
            print 'Warning, model origin for SFR object is 0, 0.'

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





