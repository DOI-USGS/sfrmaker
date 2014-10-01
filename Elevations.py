__author__ = 'aleaf'
'''
Smooth streambed elevations outside of the context of the objects in SFR classes
(works off of information in Mat1 and Mat2; generates updated versions of these files
'''
import numpy as np
import pandas as pd
import discomb_utilities as disutil

class SFRdata(object):

    def __init__(self, Mat1, Mat2, landsurface, nrow=2, ncol=2, node_column=None):

        self.m1 = pd.read_csv(Mat1).sort(['segment', 'reach'])
        self.m2 = pd.read_csv(Mat2).sort('segment')
        self.m2.index = self.m2.segment
        self.landsurface = np.fromfile(landsurface, sep=' ') # array of elevations to use (sorted by cellnumber)
        self.nrow = nrow
        self.ncol = ncol

        self.segments = sorted(np.unique(self.m1.segment))

        # node numbers for cells with SFR
        if not node_column:
            self.node_column = 'node'
            self.m1[node_column] = (self.ncol * (self.m1['row'] - 1) + self.m1['column']).astype('int')
        else:
            self.node_column = node_column

        # assign land surface elevations based on node number
        self.m1['landsurface'] = [self.landsurface[c] for c in self.m1[node_column]]

        # assign upstream segments to Mat2
        self.m2['upsegs'] = [self.m2.segment[self.m2.outseg == s].tolist() for s in self.segments]

        # check for circular routing
        c = [s for s in self.m2.segment if s == self.m2.outseg[s]]
        if len(c) > 0:
            raise ValueError('Warning! Circular routing in segments {}.\n'
                             'Fix manually in Mat2 before continuing'.format(', '.join(map(str, c))))


class smooth(SFRdata):

    def check4backwards(self):
        '''
        check for segments where the constraining elevations (in up/downstream segments) are against routing dir
        tests up to max elevation of upstream segment and min elevation of downstream segment
        a more general approach would expand the search until a higher elevation was encountered upstream,
        or until a headwater or outlet was encountered.
        '''

        # assign constraining (upstream min, downstream max) elevations to each segment in Mat2
        self.m2['upstreamMin'] = [np.min([self.seg_maxmin[useg-1][1] for useg in self.m2.ix[s, 'upsegs']])
                                  if len(self.m2.ix[s, 'upsegs']) > 0
                                  else self.seg_maxmin[s-1][0] for s in self.m2.index]
        self.m2['upstreamMax'] = [np.min([self.seg_maxmin[useg-1][0] for useg in self.m2.ix[s, 'upsegs']])
                                  if len(self.m2.ix[s, 'upsegs']) > 0
                                  else self.seg_maxmin[s-1][0] for s in self.m2.index]

        self.m2['downstreamMax'] = [self.seg_maxmin[s-1][0] for s in self.m2.outseg]
        self.m2['downstreamMin'] = [self.seg_maxmin[s-1][1] for s in self.m2.outseg]






        # set the maximum (starting) elevation in all segments to the upstream minimum
        self.seg_maxmin = np.array([(self.m2.upstreamMin.values[s-1], self.seg_maxmin[s-1][1]) for s in self.segments])

        # check for higher minimum elevations in each segment
        bw = dict([(i, maxmin) for i, maxmin in enumerate(self.seg_maxmin) if maxmin[0] < maxmin[1]])

        # assign any violating minimum elevations to downstream max, check again
        if len(bw) > 0:
            bw_inds = bw.keys()
            #self.lower_downstream(bw_inds, )


        # identify segments where the upstream min is less than the downstream max
        bw_inds1 = np.where((self.m2.upstreamMin.values - self.m2.downstreamMax.values) < 0)[0]

        # set the maximum (starting) elevation in these segments to the upstream minimum
        self.seg_maxmin[bw_inds1] = [(self.m2.upstreamMin.values[b], self.seg_maxmin[b][1]) for b in bw_inds1]

        # if first level (immediately adjoining ends up up/down segments) is backwards, check next level
        # (if min in downstream segment is > min in upstream segment)
        if len(bw_inds1) > 0:
            bw_inds2 = np.where((self.m2.upstreamMin.values - self.m2.downstreamMin.values) < 0)[0]



    def get_outsegs(self):
        '''
        returns dataframe of all downstream segments (will not work with circular routing!)
        '''
        outsegs = pd.DataFrame(self.m2.outseg)
        max_outseg = np.max(outsegs[outsegs.columns[-1]])
        knt = 2
        nsegs = len(self.m2)
        while max_outseg > 0:
            outsegs['outseg{}'.format(knt)] = [self.m2.outseg[s] if s > 0 else 0 for s in outsegs[outsegs.columns[-1]]]
            max_outseg = np.max(outsegs[outsegs.columns[-1]])
            print max_outseg
            if max_outseg == 0:
                break
            knt +=1
            if knt > nsegs:
                print 'Circular routing encountered in segment {}'.format(max_outseg)
                break

        return outsegs


    def lower_downstream(self, bw_inds, downstream):
        self.seg_maxmin[bw_inds] = np.array([(self.seg_maxmin[i][0], downstream[i]) for i in bw_inds])


    def fix_backwards(self):

        # determine min and max elevations for each segment
        self.seg_maxmin = np.array([(np.max(self.m1.landsurface[self.m1.segment == s]),
                      np.min(self.m1.landsurface[self.m1.segment == s])) for s in self.segments])

        # check if there are any min elevations in upstream segments that are below max elevations in downstream segments
        bw = self.check4backwards()

        #if len(bw) > 0:



    def smoothDownstream(self):

        for seg in self.segments:

            sbtop = self.m1.top_streambed[self.m1.segment == seg]
            len = self.m1.length_in_cell[self.m1.segment == seg]
            dist = np.cumsum(len) - 0.5 * len # cumulative distance at cell centers
