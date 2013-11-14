import numpy as np
import discomb_utilities as disutil # this is a utility I made to read in a dis file
import matplotlib.pyplot as plt

#
indis = 'NACP_SWI_Regional.dis'
inmat1 = 'SFR_GWVmat1.txt'
inmat2 = 'SFR_GWVmat2.txt'
output_seg_data = 'All_Segment_Routing.out'
output_reach_data = 'All_Reach_Routing.out'

class reach:
    def __init__(self,row,col,lay,stage,top_bed,bedK,width,length,slope,roughness):
        self.row = row
        self.col = col
        self.lay = lay
        self.stage = stage
        self.top_bed = top_bed
        self.bedK = bedK
        self.width = width
        self.length = length
        self.slope = slope
        self.roughness = roughness

class segment:
    def __init__(self,iupseg,outseg):
        self.reaches = dict()
        self.iupseg = iupseg
        self.outseg = outseg
        self.numreach = -999
        self.dist_outseg = -999
        self.outreach_row = -999
        self.outreach_col = -999
        self.inreach_row = -999
        self.inreach_col = -999
        self.dist_reaches = []
        self.circular = 0
        
# first read in the model top and all layer elevations
DX,DY,NLAY,NROW,NCOL,i = disutil.read_meta_data(indis)
deltaX = np.max(np.diff(DX))
deltaY = np.max(np.diff(DY))

print '>>>>>>>>>>>>> All Layer elevations read in <<<<<<<<<<<<<'

#  read in the SFR information from inmat1 --> reaches
print '>>> read in SFR information from %s' %inmat1
SFRdata1 = np.genfromtxt(inmat1,dtype=None, names=True, delimiter = ',')

#  read in the SFR information from inmat1 --> segments
print '>>> read in SFR information from %s' %inmat2
SFRdata2 = np.genfromtxt(inmat2,dtype=None, names=True, delimiter = ',')

numreaches = len(SFRdata1)
numsegs = len(SFRdata2)

# set up the segment data
SegData = dict()
for i,cseg in enumerate(SFRdata2['segment']):
    SegData[cseg] = segment(SFRdata2['iupseg'][i],SFRdata2['outseg'][i])

# now set up the reaches data
for i,cseg in enumerate(SFRdata1['segment']):
    SegData[cseg].reaches[SFRdata1['reach'][i]] = reach(SFRdata1['row'][i],
                                                       SFRdata1['column'][i],
                                                       SFRdata1['layer'][i],
                                                       SFRdata1['stage'][i],
                                                       SFRdata1['top_streambed'][i],
                                                       SFRdata1['bed_K'][i],
                                                       SFRdata1['width_in_cell'][i],
                                                       SFRdata1['length_in_cell'][i],
                                                       SFRdata1['bed_slope'][i],
                                                       SFRdata1['bed_roughness'][i])

# determine the location of the furthest downstream and furthest upstream reach per segment
for cseg in SegData:
    currseg = SegData[cseg]
    lastreach = np.max(currseg.reaches.keys())
    firstreach = np.min(currseg.reaches.keys())
    currseg.outreach_row = currseg.reaches[lastreach].row
    currseg.outreach_col = currseg.reaches[lastreach].col
    currseg.inreach_row = currseg.reaches[firstreach].row
    currseg.inreach_col = currseg.reaches[firstreach].col
    

# now, find the distance between the furthest downstream reach of each segment
# and the furthest upstream segment of its outseg

# make an output file for the reach data
ofp = open(output_reach_data,'w')
ofp.write('%12s'*8 %('Segment','Reach','row_from','col_from','row_to','col_to', 'Dist_Units','Dist_cells') + '\n')

for cseg in SegData:
    currseg = SegData[cseg]
    coutseg = currseg.outseg

    if coutseg == 0:
        currseg.dist_outseg = 999999
    else:
        dx = np.abs(currseg.outreach_col-SegData[coutseg].inreach_col)*deltaX
        dy = np.abs(currseg.outreach_row-SegData[coutseg].inreach_row)*deltaY
        currseg.dist_outseg = np.sqrt(dx**2 + dy**2)

    # finally, also check the distance between reaches within each segment and look for spirals
    if len(currseg.reaches.keys()) > 1:
        allrowcols = []
        for creach in currseg.reaches:
            r1 = currseg.reaches[creach].row
            c1 = currseg.reaches[creach].col
            allrowcols.append('%d_%d' %(r1,c1))
            if creach>1:
                r2 = currseg.reaches[creach-1].row
                c2 = currseg.reaches[creach-1].col
                cdist = np.sqrt(((r1-r2)*deltaY)**2+((c1-c2)*deltaX)**2)
                currseg.dist_reaches.append(cdist)
                ofp.write('%12d' * 6 %(cseg,creach,r2,c2,r1,c1) + 
                              '%12.2f'*2 %(cdist,cdist/deltaX) + '\n')
            else:
                ofp.write('%12d' * 6 %(cseg,creach,0,0,r1,c1) + 
                              '%12.2f'*2 %(0,0) + '\n')
                
        currseg.allrowcols = np.array(allrowcols)
        if len(np.unique(currseg.allrowcols)) < currseg.numreach:
            currseg.circular = 1
    
ofp.close()
            


# make an output file for the segment data
ofp = open(output_seg_data,'w')
ofp.write('%12s'*10 %('Segment','Numreach','outseg','row_from','col_from','row_to','col_to', 'Dist_Units','Dist_cells','circular') + '\n')
for cseg in SegData:
    currseg = SegData[cseg]
    coutseg = currseg.outseg
    if coutseg == 0:
        row_to = -9999
        col_to = -9999
    else:
        row_to = SegData[coutseg].inreach_row
        col_to = SegData[coutseg].inreach_col
    row_from = currseg.outreach_row
    col_from = currseg.outreach_col
    ofp.write('%12d' * 7 %(cseg,len(currseg.reaches),coutseg,row_from,col_from,row_to,col_to) + 
              '%12.2f'*2 %(currseg.dist_outseg,(currseg.dist_outseg/deltaX)) + 
              '%12d' %currseg.circular + '\n')
ofp.close()


