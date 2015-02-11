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
        print 'Checking for circular routing...'
        self.map_outsegs()
        print 'passed.'

    def check_overlapping(self):
        """checks for multiple SFR reaches in one cell; and whether more than one reach has Cond > 0
        """
        print 'Overlapping reaches check not implemented yet.'

    def check_elevations(self):
        """checks that elevations decrease monotonically downstream
        """
        outfile = 'elevation_increases.csv'
        print 'Checking for downstream rises in streambed elevation...'
        diffs = self.m1.top_streambed.diff()
        pos_diffs = self.m1.ix[diffs > 0, ['segment', 'reach', 'top_streambed']]
        pos_diffs['sb_rise'] = diffs[diffs > 0]

        if pos_diffs.reach.max() > 1:
            print 'Bad elevations found, see {}'.format(outfile)
            bad_elevs = pos_diffs[pos_diffs.reach > 1]
            bad_elevs.sort('sb_rise', ascending=False, inplace=True)
            bad_elevs.to_csv(os.path.join(self.outpath, outfile))
        else:
            print 'passed.'

        print 'Checking for MODFLOW altitude errors...'
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
            return
        print 'passed.'