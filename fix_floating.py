import os
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


infile='br_stop.dif' # from Daniel's CORR-SFR script
GWVmat1='BR_Lay_minSlope.COR' # GWV SFR input mat1, from COFF-SFR
outfile='BR_GVW_SFRmat1.dat' # revised GWV SFR input mat1
float_cutoff=-200 # process all segments with reaches that are float_cutoff above land surface
increment=0.01 # if a minimum in the dem is lower than the segment end, new elev. will be end + increment


# bring in segment, reach, streamtop, land surface, and residuals
infile_data=np.genfromtxt(infile, names=True, dtype=None)
float_inds=np.where(infile_data['DIF']>float_cutoff)[0]
floaters=infile_data[float_inds]
segments=list(np.unique(floaters['MSEG']))

# bring in reach lengths
GWVmat1_data=np.genfromtxt(GWVmat1, skiprows=1,dtype=None) # can't use names because of commas in header

seg_distdict=defaultdict(list)
TOPNEWdict=defaultdict(list)
STOP1dict=defaultdict(list)
STOP2dict=defaultdict(list)
for segnum in segments:

    seg_inds=list(np.where(infile_data['MSEG']==segnum)[0])
    segment=infile_data[seg_inds]
    
    if len(segment)<3: # no interior points; keep ending elevations
        continue
    
    seg_indsGWVmat=list(np.where(GWVmat1_data['f6']==segnum)[0])
    seg_distances=list(np.cumsum(GWVmat1_data['f8'][seg_indsGWVmat])) # cumulative distances along segment
    
    
    # for now, start with existing STOP at RCH1. Could improve by first having a routine that corrects all segment ends to eliminate floaters.
    STOP2=[]
    STOP2.append(segment['STOP'][0]) 
    end=segment['STOP'][-1]
    
    descending=False # If true, means that current reach is lower than all previous
    flattening=False # means that STOP2 elevation is getting close to segment end elevation
    if STOP2[0]>segment['TOPNEW'][1]:
        descending=True
    
    for i in range(len(segment['TOPNEW']))[1:-1]:
        dem_el=segment['TOPNEW'][i]
        
        # in case a minimum in DEM is below segment end, set elev. to increment above, interpolate to end
        if segment['TOPNEW'][i]<end:
            dem_el=end+increment
            nextlower=end
            slope=(nextlower-STOP2[-1])/(seg_distances[len(segment)-1]-seg_distances[len(STOP2)-1])
            sub_inds=range(len(segment)-1)[len(STOP2):(len(segment)-1)]
            for s in sub_inds:
                dist=seg_distances[s]-seg_distances[s-1]
                STOPint=STOP2[-1]+slope*dist
                STOP2.append(STOPint)            
            break
        
        # if last rch is the next, need to work with end STOP elevation
        if i==np.max(range(len(segment['TOPNEW']))[1:-1]):
            nextlower=end
            
            # current reach is below all previous, and above end
            if descending:
                STOP2.append(dem_el)
                
            # current is not below all previous
            # perform linear interpolation between previous low and end, recording points            
            else:
                slope=(nextlower-STOP2[-1])/(seg_distances[i+1]-seg_distances[len(STOP2)-1])
                sub_inds=range(i+1)[len(STOP2):i+1]
                for s in sub_inds:
                    dist=seg_distances[s]-seg_distances[s-1]
                    STOPint=STOP2[-1]+slope*dist
                    STOP2.append(STOPint)
                    
        # next rch is lower but also lower than the end 
        elif segment['TOPNEW'][i+1]<end:
            flattening=True
            continue
        
        # next rch elevation is lower        
        elif segment['TOPNEW'][i+1]<dem_el:
            nextlower=segment['TOPNEW'][i+1]
            if descending and not flattening:
                STOP2.append(dem_el)
            
            # next reach is lower than current but higher than a previous reach, keep going    
            elif nextlower>STOP2[-1]:
                continue
            
            # next reach is lower than current, and all previous reaches
            # perform linear interpolation between previous low and next, recording points
            else:   
                slope=(nextlower-STOP2[-1])/(seg_distances[i+1]-seg_distances[len(STOP2)-1])
                sub_inds=range(i+1)[len(STOP2):i+1]
                for s in sub_inds:
                    dist=seg_distances[s]-seg_distances[s-1]
                    STOPint=STOP2[-1]+slope*dist
                    STOP2.append(STOPint)
                descending=True
                continue
        
        # next rch elevation is equal or higher
        elif segment['TOPNEW'][i+1]>=dem_el:
            
            # current reach is still lower than all previous, record
            if descending and not flattening:
                STOP2.append(dem_el)
                descending=False
            else:
                continue
    
    STOP2.append(segment['STOP'][-1])
    
    # if some condition isn't handled properly by one of the above statements, chances are this will occur:
    if len(STOP2)!=len(seg_distances):
        raise NameError('unequal lengths in seg_distances and STOPs!')    
    STOP2dict[segnum] = STOP2
    seg_distdict[segnum] = seg_distances
    TOPNEWdict[segnum] = list(segment['TOPNEW'])
    STOP1dict[segnum] = list(segment['STOP'])
     
# Pull everything together into new GWV SFR mat. 1
input_file=open(GWVmat1).readlines()
ofp=open(outfile,'w')
ofp.write(input_file[0])
for line in input_file[1:]:
    segment=int(line.split()[6])
    if segment in STOP2dict:
        reach=int(line.split()[5])
        linestart=','.join(map(str,line.split()[:3]))
        lineend=','.join(map(str,line.split()[5:-1]))
        STOP=STOP2dict[segment][reach-1]
        SBOT=STOP2dict[segment][reach-1]-1
        ofp.write('%s,%s,%s,%s\n' %(linestart,STOP,SBOT,lineend))
    else:
        ofp.write(','.join(map(str,line.split()))+'\n')
ofp.close()

segs2plot=[1081,324,1188,1060,1201,1076,341,1113,1131,949,523,408,469,1064,200,300,400,500,600,700,800,900,1000,1100,1200]
segs2plot=sorted(segs2plot)
for seg in segs2plot:
    if seg not in STOP2dict:
        print seg
pdf=PdfPages('selected_segments.pdf')
for seg in segs2plot:
    fig=plt.figure(seg)
    plt.plot(seg_distdict[seg],TOPNEWdict[seg],'b',seg_distdict[seg],STOP1dict[seg],'r',seg_distdict[seg],STOP2dict[seg],'g')
    plt.title('segment ' + str(seg))
    plt.xlabel('distance along segment (ft.)')
    plt.ylabel('Elevation (ft)')
    pdf.savefig(fig)
pdf.close()
plt.close('all')




















'''
#look at some of the segments that were initially some of the worst floaters
plt.close('all')
for seg in segs2plot:
    fig=plt.figure(seg)   plt.plot(seg_distdict[seg],TOPNEWdict[seg],'b',seg_distdict[seg],STOP1dict[seg],'r',seg_distdict[seg],STOP2dict[seg],'g')
    plt.draw()

    pdf.savefig(fig)
pdf.close()
'''

'''
# look at results in comparison to landsurface and previous stream elevation for one segment

plt.close()
plt.figure(4)
plt.plot(seg_distances,segment['TOPNEW'],'b',seg_distances,segment['STOP'],'r',seg_distances,STOP2,'g')
plt.draw()

'''













'''
Another approach is to find local minima in the DEM and then interpolate between those, but this is problematic in that it depends on their being minima to begin with (what if the DEM is convex?). And then at least a second round of corrections must be made to get rid of any resulting floating segments.

# identify minima in DEM stream trace
diffs=[]
for i in range(len(segment))[1:]:
    diff=segment['TOPNEW'][i]-segment['TOPNEW'][i-1]
    diffs.append(diff)

min_inds=[]
for i in range(len(diffs))[:-1]:
    if diffs[i]<0 and diffs[i+1]>0:
        minimum_ind=i+1
        min_inds.append(minimum_ind)


# get list of TOPNEW values at minima; insert at front and append first and last STOP values in segment
# (for now, maintain segment ends at corrected values)
minima=list(segment['TOPNEW'][min_inds])
min_inds.insert(0,0)
min_inds.append(range(len(segment))[-1])
minima.insert(0,segment['STOP'][0])
minima.append(segment['STOP'][-1])

# add in routine that enforce monotonicity in minima list




# Linear interpolation between TOPNEW minima using segment stream distances
STOP2=[]       
for i in range(len(min_inds))[1:]:
    STOP2.append(minima[i-1])
    slope=(minima[i]-minima[i-1])/(minima_dist[i]-minima_dist[i-1])
    sub_inds=range(min_inds[i])[min_inds[i-1]:min_inds[i]]
    for s in sub_inds[1:]:
        dist=seg_distances[s]-seg_distances[s-1]
        STOPint=STOP2[-1]+slope*dist
        STOP2.append(STOPint)
STOP2.append(minima[-1])
        
# Second loop that forces all floating STOP2 elevations to land surface (TOPNEW)

# The check again for minima; if minima found, drop next segment to next lowest STOP elevation
# on the first round, instead of interpolating between minima, could instead simply always look for lower elevations going downstream.


    
    
    
    
    
    
    
    
    
    
    
    
    
    

plt.plot(segment['MRCH'],segment['TOPNEW'],'b',segment['MRCH'],segment['STOP'],'r',min_inds,minima,'g')
plt.hold = True
plt.plot(segment['MRCH'],segment['STOP'])
plt.draw()
plt.close()

'''
