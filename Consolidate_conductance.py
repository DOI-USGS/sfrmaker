# script for running Assign Layers outside of SFR_main
__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc
import numpy as np

class cprops:
    def __init__(self, leng, wid, k, thick):
        self.leng = leng
        self.wid = wid
        self.k = k
        self.thick = thick
        self.cond = leng*wid*k/thick
        self.newk = k
        self.newlen = leng
        self.newthick = thick
        self.newleng = leng
infile = 'SFR_input_NACP.xml'

SFRdata = SFRc.SFRInput(infile)


SFRops = SFRc.SFROperations(SFRdata)

# now read in the reach information
SFRinfo = np.genfromtxt(SFRdata.MAT1, delimiter=',', names=True, dtype=None)
rcl = list()
condprops = dict()
# find each individual cell (identified by RCL)
for i in np.arange(len(SFRinfo)):
    indy1 = (SFRinfo['row'][i], SFRinfo['column'][i], SFRinfo['layer'][i])
    rcl.append(indy1)
rcl = list(set(rcl))
rcl_segreach_lookup = {ind: list() for ind in rcl}

# link up segments and reaches with cells and read in conductance parameters
for i in np.arange(len(SFRinfo)):
    indy1 = (SFRinfo['row'][i], SFRinfo['column'][i], SFRinfo['layer'][i])
    indy2 = (SFRinfo['segment'][i], SFRinfo['reach'][i])
    rcl_segreach_lookup[indy1].append(indy2)
    condprops[indy2] = cprops(SFRinfo['length_in_cell'][i],
                             SFRinfo['width_in_cell'][i],
                             SFRinfo['bed_K'][i],
                             SFRinfo['bed_thickness'][i])
maxCloc = dict()
maxWloc = dict()
for ccell in rcl:
    if ccell == (92,82,3):
        print 'ginger'
    if len(rcl_segreach_lookup[ccell]) > 1:
        #  more than one segment/reach in this cell. let's consolidate!
        #  first read in the important information and then
        maxwid = -999
        totcond = 0
        totreaches = len(rcl_segreach_lookup[ccell])

        for csr in rcl_segreach_lookup[ccell]:
            # total up the conductivity
            totcond += condprops[csr].cond
            # find the stream reach with maximum width - we will assign COND to that cell
            if condprops[csr].wid > maxwid:
                maxWloc[ccell] = csr
                maxwid = condprops[csr].wid
        for csr in rcl_segreach_lookup[ccell]:
            # make len=1.0 for the cells so COND will be calculated as (1.0 * K * W) /t
            condprops[maxWloc[ccell]].newlen = 1.0
            # now assign (COND * t)/maxwid to be newk for the cell and bedKmin to the others
            if csr == maxWloc[ccell]:
                condprops[maxWloc[ccell]].newk = (condprops[maxWloc[ccell]].thick * totcond)/maxwid
                condprops[maxWloc[ccell]].newthick = condprops[maxWloc[ccell]].thick
                condprops[maxWloc[ccell]].newleng = 1.0

            else:
                condprops[csr].newk = SFRdata.bedKmin
                condprops[csr].newthick = 1.0

# now, rewrite the file in SFRdata.MAT1
header = open(SFRdata.MAT1, 'r').readline().strip()
ofp = open(SFRdata.MAT1, 'w')
ofp.write(header + '\n')

for i in np.arange(len(SFRinfo)):
    indy1 = (SFRinfo['row'][i], SFRinfo['column'][i], SFRinfo['layer'][i])
    indy2 = (SFRinfo['segment'][i], SFRinfo['reach'][i])
    #  if more than one reach in a cell, consolidate
    if len(rcl_segreach_lookup[indy1]) > 1:
        outline = '{0:d},{1:d},{2:d},{3:e},{4:e},{5:d},{6:d},{7:e},{8:e},{9:e},{10:e},{11:e},{12:e}\n'.format(
            SFRinfo['row'][i],
            SFRinfo['column'][i],
            SFRinfo['layer'][i],
            SFRinfo['stage'][i],
            SFRinfo['top_streambed'][i],
            SFRinfo['reach'][i],
            SFRinfo['segment'][i],
            SFRinfo['width_in_cell'][i],
            condprops[indy2].newleng,
            condprops[indy2].newk,
            condprops[indy2].newthick,
            SFRinfo['bed_slope'][i],
            SFRinfo['bed_roughness'][i]
        )
    #  otherwise, just write it out
    else:
        outline = '{0:d},{1:d},{2:d},{3:e},{4:e},{5:d},{6:d},{7:e},{8:e},{9:e},{10:e},{11:e},{12:e}\n'.format(
            SFRinfo['row'][i],
            SFRinfo['column'][i],
            SFRinfo['layer'][i],
            SFRinfo['stage'][i],
            SFRinfo['top_streambed'][i],
            SFRinfo['reach'][i],
            SFRinfo['segment'][i],
            SFRinfo['width_in_cell'][i],
            SFRinfo['length_in_cell'][i],
            SFRinfo['bed_K'][i],
            SFRinfo['bed_thickness'][i],
            SFRinfo['bed_slope'][i],
            SFRinfo['bed_roughness'][i]
        )
    # finally, write the line to the file
    ofp.write(outline)
ofp.close()
SFRoutput = SFRc.SFRoutput(SFRdata)
SFRoutput.build_SFR_package()
