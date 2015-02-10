__author__ = 'aleaf'

import os
import numpy as np
import pandas as pd
from postproc import SFRdata



class diagnostics(SFRdata):

    def __init__(self, Mat1=None, Mat2=None, sfr=None, node_column=None,
                     mfpath=None, mfnam=None, mfdis=None, outpath=os.getcwd()):

        SFRdata.__init__(self, Mat1=Mat1, Mat2=Mat2, sfr=sfr, node_column=node_column,
                         mfpath=mfpath, mfnam=mfnam, mfdis=mfdis)

    def check_numbering(self):
        """checks for continuity in segment and reach numbering
        """
        self.m1.sort(['segment', 'reach'], inplace=True)
        self.m2.sort('segment', inplace=True)
        uniquesegs = self.m2.segment.values
        if np.max(np.diff(uniquesegs)) != 1 or np.min(np.diff(uniquesegs)) != 1:
            print 'gaps in segment numbering!, check Mat2'
        segments, reaches = self.m1.segment.values, self.m1.reach.values

        for seg in uniquesegs:
            r = reaches[segments == seg]
            if len(r) > 1:
                if np.max(np.diff(r)) != 1 or np.min(np.diff(r)) != 1:
                    print 'Invalid reach numbering at segment {:.0f}!'.format(seg)
                    print self.m1.reach[self.m1.segment == seg]
                    break

    def check_routing(self):
        """checks for breaks in routing and does comprehensive check for circular routing
        """
        pass

    def check_overlapping(self):
        """checks for multiple SFR reaches in one cell; and whether more than one reach has Cond > 0
        """
        pass

    def check_elevations(self):
        """checks that elevations decrease monotonically downstream
        """
        outfile = 'elevation_increases.csv'
        print 'Checking for streambed elevation violations...'
        pos_diffs = self.m1.ix[self.m1.top_streambed.diff() > 0, ['segment', 'reach', 'top_streambed']]

        if pos_diffs.reach.max() > 1:
            print 'Bad elevations found, see {}'.format(outfile)
            bad_elevs = pos_diffs[pos_diffs.reach > 1]
            bad_elevs.to_csv(os.path.join(self.outpath, outfile))
        else:
            print 'Elevation check passed.'