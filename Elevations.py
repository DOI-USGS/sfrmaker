__author__ = 'aleaf'
'''
Smooth streambed elevations outside of the context of the objects in SFR classes
(works off of information in Mat1 and Mat2; generates updated versions of these files
'''
import numpy as np
import pandas as pd
import discomb_utilities as disutil


class SFRdata(object):

    def __init__(self, Mat1, Mat2, landsurface=None, nrow=2, ncol=2, node_column=False, MFdis=None):

        # read Mats 1 and 2
        self.MAT1 = Mat1
        self.MAT2 = Mat2
        self.m1 = pd.read_csv(Mat1).sort(['segment', 'reach'])
        self.m2 = pd.read_csv(Mat2).sort('segment')
        self.m2.index = self.m2.segment

        self.nrow = nrow
        self.ncol = ncol
        self.MFdis = MFdis

        self.segments = sorted(np.unique(self.m1.segment))

        # node numbers for cells with SFR
        if not node_column:
            self.node_column = 'node'
            self.gridtype = 'structured'
            self.m1[self.node_column] = (self.ncol * (self.m1['row'] - 1) + self.m1['column']).astype('int')
        else:
            self.node_column = node_column

        self.node_attribute = self.node_column # compatibility with sfr_plots and sfr_classes

        if landsurface:
            self.landsurface = np.fromfile(landsurface, sep=' ') # array of elevations to use (sorted by cellnumber)

            # assign land surface elevations based on node number
            self.m1['landsurface'] = [self.landsurface[n] for n in self.m1[self.node_column]]

        # assign upstream segments to Mat2
        self.m2['upsegs'] = [self.m2.segment[self.m2.outseg == s].tolist() for s in self.segments]

        # check for circular routing
        c = [s for s in self.m2.segment if s == self.m2.outseg[s]]
        if len(c) > 0:
            raise ValueError('Warning! Circular routing in segments {}.\n'
                             'Fix manually in Mat2 before continuing'.format(', '.join(map(str, c))))


    def read_DIS(self, MFdis):


        DX, DY, NLAY, self.NROW, self.NCOL, i = disutil.read_meta_data(MFdis)

        # get layer tops/bottoms
        self.layer_elevs = np.zeros((NLAY+1, self.NROW, self.NCOL))
        for c in range(NLAY + 1):
            tmp, i = disutil.read_nrow_ncol_vals(MFdis, self.NROW, self.NCOL, 'float', i)
            self.layer_elevs[c, :, :] = tmp

        # make dictionary of model top elevations by cellnum
        self.elevs_by_cellnum = {}
        for c in range(self.NCOL):
            for r in range(self.NROW):
                cellnum = r*self.NCOL + c + 1
                self.elevs_by_cellnum[cellnum] = self.layer_elevs[0, r, c]


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

    def smooth_segment_ends(self, report_file='smooth_segment_ends.txt'):
        '''
        smooth segment end elevations so that they decrease monotonically down the stream network
        '''

        print '\nSmoothing segment ends...\n'
        # open a file to report a summary of the elevation adjustments
        self.ofp = open(report_file, 'w')
        self.ofp.write('Segment end smoothing report:\n'
                  'segment,max_elev,min_elev,downstream_min_elev\n')

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

        # smooth any backwards elevations
        if len(bw) > 0:

            # creates dataframe of all outsegs for each segment
            self.map_outsegs()

            # iterate through the mapped outsegs, replacing minimum elevations until there are no violations
            self.fix_backwards_ends(bw)

        print 'segment ends smoothing finished in {} iterations.\n' \
              'segment ends saved to {}\n' \
              'See {} for report.'\
            .format(self.smoothing_iterations, self.MAT2[:-4] + '_elevs.csv', report_file)

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

        elevations_summary = self.m2[['segment','upstreamMax','upstreamMin','Max','Min','downstreamMax','downstreamMin']]
        elevations_summary.to_csv(self.ofp)
        self.ofp.close()

        # save new copy of Mat2 with segment end elevations
        self.m2 = self.m2.drop(['upstreamMax', 'upstreamMin', 'downstreamMax', 'downstreamMin', 'upsegs'], axis=1)
        self.m2.to_csv(self.MAT2[:-4] + '_elevs.csv', index=False)


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
            self.ofp.write('{},{},{}\n'.format(i+1, self.seg_maxmin[i], downstream_elevs[i]))


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

        knt = 1
        for level in self.outsegs.columns:

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

            knt += 1

        self.smoothing_iterations = knt


    def interpolate(self, istart, istop):
        '''
        istart = index of previous minimum elevation
        istop = index of current minimum elevation
        dx = interpolation distance
        dS = change in streambed top over interpolation distance
        '''
        dx = self.cdist[istart] - self.cdist[istop]
        dS = self.minelev - self.elevs[istop]

        dist = 0

        if dS == 0:
            slope = 0
        else:
            slope = dS/dx

        for i in np.arange(istart+1, istop+1):

            dist += self.cdist[i-1] - self.cdist[i]
            self.sm.append(self.minelev - dist * slope)

            self.ofp.write('{:.0f},{:.0f},{:.2f},{:.2f},{:.2f},{:.2e},{:.2f}\n'
                           .format(self.seg, i, self.elevs[i-1], self.minelev, dist, slope, self.sm[-1]))

        # reset the minimum elevation to current
        self.minelev = self.sm[-1]


    def smooth_segment_interiors(self, report_file='smooth_segment_interiors.txt'):

        print '\nSmoothing segment interiors...\n'

        # open a file to report on interpolations
        self.ofp = open(report_file, 'w')
        self.ofp.write('segment, reach, land_surface, minelev, dist, slope, sb_elev\n')

        for seg in self.segments:
            self.seg = seg

            try:
                max, min = self.m2['Max'], self.m2['Min']
            except:
                raise IndexError("Max, Min elevation columns not found in Mat2. Call smooth_segment_ends method first")

            # make a view of the Mat 1 dataframe that only includes the segment
            df = self.m1[self.m1.segment == seg].sort('reach')

            # start with land surface elevations along segment
            self.elevs = df.landsurface.values


            # get start and end elevations from Mat2; if they are equal, continue
            start, end = self.m2.ix[seg, ['Max', 'Min']]
            if start == end:
                self.sm = [start] * len(self.elevs)
                self.m1.loc[self.m1.segment == seg, 'top_streambed'] = self.sm
                continue
            else:
                self.sm = [start]
                self.elevs[-1] = end

            # calculate cumulative distances at cell centers
            lengths = df.length_in_cell.values
            self.cdist = np.cumsum(lengths) - 0.5 * lengths

            self.minelev = start
            minloc = 0
            nreaches = len(self.elevs)

            for i in range(nreaches)[1:]:

                # if the current elevation is equal to or below the minimum
                if self.elevs[i] <= self.minelev:

                    # if the current elevation is above the end, interpolate from previous minimum to current
                    if self.elevs[i] >= end:
                        self.interpolate(minloc, i)
                        minloc = i

                    # otherwise if it is below the end, interpolate from previous minimum to end
                    elif self.elevs[i] <= end:

                        self.interpolate(minloc, nreaches-1)
                        break
                else:
                    continue

            # update Mat1 with smoothed streambed tops
            self.m1.loc[self.m1.segment == seg, 'top_streambed'] = self.sm

        # save updated Mat1
        self.m1.to_csv(self.MAT1[:-4] + '_elevs.csv', index=False)
        self.ofp.close()
        print 'Done, updated Mat1 saved to {}.\nsee {} for report.'.format(self.MAT1[:-4] + '_elevs.csv', report_file)
