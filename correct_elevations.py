import numpy as np
import SFR_classes as SFRc
import shutil

class segs:
    def __init__(self, segnum):
        self.seg = segnum
        self.reaches = None
        self.strtop = None
        self.stage = None
        self.slope = None
        self.altered = False

infile = 'SFR_input_NACP.xml'

SFRdata = SFRc.SFRInput(infile)

inreachdat = np.genfromtxt(SFRdata.MAT1, names=True, dtype=None, delimiter=',')
fileheader = open(SFRdata.MAT1).readline()
allsegs = np.unique(inreachdat['segment'])

segdata = dict()

# read in the segment data from the input file
for cseg in allsegs:
    segdata[cseg] = segs(cseg)
    inds = np.squeeze(np.where(inreachdat['segment'] == cseg))
    segdata[cseg].reaches = np.atleast_1d(inreachdat['reach'][inds])
    segdata[cseg].strtop = np.atleast_1d(inreachdat['top_streambed'][inds])
    segdata[cseg].stage = np.atleast_1d(inreachdat['stage'][inds])
    segdata[cseg].slope = np.atleast_1d(inreachdat['bed_slope'][inds])
    segdata[cseg].lengths = np.atleast_1d(inreachdat['length_in_cell'][inds])
    segdata[cseg].row = np.atleast_1d(inreachdat['row'][inds])
    segdata[cseg].col = np.atleast_1d(inreachdat['column'][inds])
    segdata[cseg].layer = np.atleast_1d(inreachdat['layer'][inds])
    segdata[cseg].width = np.atleast_1d(inreachdat['width_in_cell'][inds])
    segdata[cseg].bedK = np.atleast_1d(inreachdat['bed_K'][inds])
    segdata[cseg].bedthick = np.atleast_1d(inreachdat['bed_thickness'][inds])
    segdata[cseg].bedrough = np.atleast_1d(inreachdat['bed_roughness'][inds])



    if np.min(segdata[cseg].strtop) < SFRdata.min_elev or np.min(segdata[cseg].slope) < SFRdata.minimum_slope:
        segdata[cseg].altered = True
        if len(segdata[cseg].reaches) > 1:
            for i in reversed(segdata[cseg].reaches):
                # now fix the elevations due to lowest stage being below the min elev allowable
                if i == len(segdata[cseg].reaches):
                    if segdata[cseg].strtop[i-1] < SFRdata.min_elev:
                        segdata[cseg].strtop[i-1] = SFRdata.min_elev
                else:
                    if segdata[cseg].strtop[i-1] <= segdata[cseg].strtop[i]:
                        cslope = np.max(SFRdata.minimum_slope, segdata[cseg].slope[i-1])
                        segdata[cseg].strtop[i-1] = (segdata[cseg].strtop[i] +
                                                    cslope*segdata[cseg].lengths[i-1] +
                                                    cslope*segdata[cseg].lengths[i]
                        )
                # now fix the slopes if they are too low
                if segdata[cseg].slope[i-1] < SFRdata.minimum_slope:
                    segdata[cseg].slope[i-1] = SFRdata.minimum_slope
        else:
            if segdata[cseg].strtop < SFRdata.min_elev:
                segdata[cseg].strtop = SFRdata.min_elev
                segdata[cseg].strtop = np.atleast_1d(segdata[cseg].strtop)
            if segdata[cseg].slope < SFRdata.minimum_slope:
                segdata[cseg].slope = SFRdata.minimum_slope
                segdata[cseg].slope = np.atleast_1d(segdata[cseg].slope)

# now write out the entire file again with updated values
print 'saving backup of file {0:s} to --> {1:s}'.format(SFRdata.MAT1, SFRdata.MAT1 + '.backup')
shutil.copyfile(SFRdata.MAT1, SFRdata.MAT1 + '.backup')
ofp = open(SFRdata.MAT1, 'w')

# write out the header
ofp.write(fileheader.strip() + '\n')
for cseg in allsegs:
    currsegdata = segdata[cseg]
    if segdata[cseg].altered:
        print 'Made changes to segment --> {0:d}'.format(cseg)
    for i in np.arange(len(currsegdata.reaches)):
        ofp.write('{0:d},'.format(currsegdata.row[i]))
        ofp.write('{0:d},'.format(currsegdata.col[i]))
        ofp.write('{0:d},'.format(currsegdata.layer[i]))
        ofp.write('{0:f},'.format(currsegdata.strtop[i] + SFRdata.stream_depth))
        ofp.write('{0:f},'.format(currsegdata.strtop[i]))
        ofp.write('{0:d},'.format(currsegdata.reaches[i]))
        ofp.write('{0:d},'.format(cseg))
        ofp.write('{0:f},'.format(currsegdata.width[i]))
        ofp.write('{0:f},'.format(currsegdata.lengths[i]))
        ofp.write('{0:f},'.format(currsegdata.bedK[i]))
        ofp.write('{0:f},'.format(currsegdata.bedthick[i]))
        ofp.write('{0:f},'.format(currsegdata.slope[i]))
        ofp.write('{0:f}\n'.format(currsegdata.bedrough[i]))


ofp.close()

