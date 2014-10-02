__author__ = 'aleaf'
'''
Smooth streambed elevations outside of the context of the objects in SFR classes
(works off of information in Mat1 and Mat2; generates updated versions of these files
'''
import numpy as np
import pandas as pd
import discomb_utilities as disutil







class SFRdata(object):

    def __init__(self, Mat1, Mat2, landsurface=None, nrow=2, ncol=2, node_column=None):

        # read Mats 1 and 2
        self.m1 = pd.read_csv(Mat1).sort(['segment', 'reach'])
        self.m2 = pd.read_csv(Mat2).sort('segment')
        self.m2.index = self.m2.segment

        self.nrow = nrow
        self.ncol = ncol

        self.segments = sorted(np.unique(self.m1.segment))

        # node numbers for cells with SFR
        if not node_column:
            self.node_column = 'node'
            self.m1[node_column] = (self.ncol * (self.m1['row'] - 1) + self.m1['column']).astype('int')
        else:
            self.node_column = node_column

        if landsurface:
            self.landsurface = np.fromfile(landsurface, sep=' ') # array of elevations to use (sorted by cellnumber)

            # assign land surface elevations based on node number
            self.m1['landsurface'] = [self.landsurface[c] for c in self.m1[node_column]]

        # assign upstream segments to Mat2
        self.m2['upsegs'] = [self.m2.segment[self.m2.outseg == s].tolist() for s in self.segments]

        # check for circular routing
        c = [s for s in self.m2.segment if s == self.m2.outseg[s]]
        if len(c) > 0:
            raise ValueError('Warning! Circular routing in segments {}.\n'
                             'Fix manually in Mat2 before continuing'.format(', '.join(map(str, c))))


    def map_outsegs(self):
        '''
        from Mat2, returns dataframe of all downstream segments (will not work with circular routing!)
        '''
        self.outsegs = pd.DataFrame(self.m2.outseg)
        max_outseg = np.max(self.outsegs[self.outsegs.columns[-1]])
        knt = 2
        nsegs = len(self.m2)
        while max_outseg > 0:
            self.outsegs['outseg{}'.format(knt)] = [self.m2.outseg[s] if s > 0 else 0
                                                    for s in self.outsegs[self.outsegs.columns[-1]]]
            max_outseg = np.max(self.outsegs[self.outsegs.columns[-1]])
            if max_outseg == 0:
                break
            knt +=1
            if knt > nsegs:
                print 'Circular routing encountered in segment {}'.format(max_outseg)
                break



class smooth(SFRdata):

    def segment_ends(self):
        '''
        check for segments where the constraining elevations (in up/downstream segments) are against routing dir
        tests up to max elevation of upstream segment and min elevation of downstream segment
        a more general approach would expand the search until a higher elevation was encountered upstream,
        or until a headwater or outlet was encountered.
        '''

        # set initial max / min elevations for segments based on max and min land surface elevations along each segment
        self.seg_maxmin = np.array([(np.max(self.m1.landsurface[self.m1.segment == s]),
                      np.min(self.m1.landsurface[self.m1.segment == s])) for s in self.segments])

        # determine constraining upstream min elevations for each segment in Mat2
        upstreamMin = np.array([np.min([self.seg_maxmin[useg-1][1] for useg in self.m2.ix[s, 'upsegs']])
                                  if len(self.m2.ix[s, 'upsegs']) > 0
                                  else self.seg_maxmin[s-1][0] for s in self.m2.index])

        # if the upstream minimum elevation is lower than the max elevation in the segment,
        # reset the max to the upstream minimum
        self.seg_maxmin = np.array([(upstreamMin[s-1], self.seg_maxmin[s-1][1])
                                    if upstreamMin[s-1] < self.seg_maxmin[s-1][0]
                                    else self.seg_maxmin[s-1] for s in self.segments])

        # check for higher minimum elevations in each segment
        bw = dict([(i, maxmin) for i, maxmin in enumerate(self.seg_maxmin) if maxmin[0] < maxmin[1]])


        if len(bw) > 0:

            self.map_outsegs() # creates dataframe of all outsegs for each segment

            self.fix_backwards_ends(bw) # iterate through the mapped outsegs, replacing minimum elevations until there are no violations


        # populate Mat2 dataframe with bounding elevations for upstream and downstream segments, so they can be checked
        self.m2['upstreamMin'] = [np.min([self.seg_maxmin[useg-1][1] for useg in self.m2.ix[s, 'upsegs']])
                                  if len(self.m2.ix[s, 'upsegs']) > 0
                                  else self.seg_maxmin[s-1][0] for s in self.m2.index]
        self.m2['upstreamMax'] = [np.min([self.seg_maxmin[useg-1][0] for useg in self.m2.ix[s, 'upsegs']])
                                  if len(self.m2.ix[s, 'upsegs']) > 0
                                  else self.seg_maxmin[s-1][0] for s in self.m2.index]
        self.m2['Max'] = self.seg_maxmin[:, 0]
        self.m2['Min'] = self.seg_maxmin[:, 1]
        self.m2['downstreamMax'] = [self.seg_maxmin[s-1][0] for s in self.m2.outseg]
        self.m2['downstreamMin'] = [self.seg_maxmin[s-1][1] for s in self.m2.outseg]


    def replace_downstream(self, bw, level, ind):
        '''
        replace minimum elevations in segments with either the max or min elevation in downstream segment
        bw = dict with keys that are indices of segements to modify
        level = column of outsegs to reference in outsegs table
        ind = 0 (replace with downstream max elev) or 1 (min elev)
        get downstream max elevations for level from outsegs at that level (for segments with backwards elevs)
        '''
        # make list of downstream elevations (max or min, depending on ind)
        # if segment is an outlet, use minimum elevation of segment
        downstream_elevs = [self.seg_maxmin[self.outsegs.ix[s, level] - 1][ind]
                            if self.outsegs.ix[s, level] > 0
                            else np.min(self.seg_maxmin[s-1]) for s in self.segments]

        # assign any violating minimum elevations to downstream elevs from above
        bw_inds = bw.keys()
        #print 'segment {}: {}, downseg: {}'.format(bw_inds[0]+1, self.seg_maxmin[bw_inds[0]], downstream_elevs[bw_inds[0]])
        self.seg_maxmin[bw_inds] = np.array([(self.seg_maxmin[i][0], downstream_elevs[i]) for i in bw_inds])
        for i in bw_inds:
            print 'segment {}: {}, downseg: {}'.format(i+1, self.seg_maxmin[i], downstream_elevs[i])


    def replace_upstream(self):
        '''
        replace maximum elevations in segments with minimum elevations in upstream segments,
        in segments where the max is above the upstream minimum
        '''
        # make list of min upstream elevations
        # if segment is a headwater, use maximum elevation of segment
        upstreamMin = np.array([np.min([self.seg_maxmin[useg-1] for useg in self.m2.ix[s, 'upsegs']])
                       if len(self.m2.ix[s, 'upsegs']) > 0
                       else self.seg_maxmin[s-1][0] for s in self.segments])

        # if the upstream minimum elevation is lower than the max elevation in the segment,
        # reset the max to the upstream minimum
        self.seg_maxmin = np.array([(upstreamMin[s-1], self.seg_maxmin[s-1][1])
                                    if upstreamMin[s-1] < self.seg_maxmin[s-1][0]
                                    else self.seg_maxmin[s-1] for s in self.segments])


    def fix_backwards_ends(self, bw):

        for level in self.outsegs.columns:
            print 'segment 48: {}, {}'.format(self.seg_maxmin[47][0], self.seg_maxmin[47][1])
            print 'segment 50: {}, {}'.format(self.seg_maxmin[49][0], self.seg_maxmin[49][1])
            # replace minimum elevations in backwards segments with downstream maximums
            self.replace_downstream(bw, level, 0)

            # check for higher minimum elevations in each segment
            bw = dict([(i, maxmin) for i, maxmin in enumerate(self.seg_maxmin) if maxmin[0] < maxmin[1]])

            # if still backwards elevations, try using downstream min elevations
            if len(bw) > 0:

                # replace minimum elevations in backwards segments with downstream minimums
                self.replace_downstream(bw, level, 1)

            # in segments where the max elevation is higher than the upstream minimum,
            # replace max elev with upstream minimum
            self.replace_upstream()

            # check again for higher minimum elevations in each segment
            bw = dict([(i, maxmin) for i, maxmin in enumerate(self.seg_maxmin) if maxmin[0] < maxmin[1]])

            # stop if there are no longer any backwards segments, otherwise continue
            if len(bw) == 0:
                break


    def segment_interiors(self):

        for seg in self.segments:

            sbtop = self.m1.top_streambed[self.m1.segment == seg]
            len = self.m1.length_in_cell[self.m1.segment == seg]
            dist = np.cumsum(len) - 0.5 * len # cumulative distance at cell centers
