# Program to correct SFR reach streamtop elevations using land surface elevations, while enforcing downstream monotonicity

import os
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import STOP_compare
import pdb

# Global Input file for SFR utilities
infile="SFR_setup.in"

# Get input files locations
infile=open(infile,'r').readlines()
inputs=defaultdict()
for line in infile:
    if "#" not in line.split()[0]:
        varname=line.split("=")[0]
        var=line.split("=")[1].split()[0]
        inputs[varname]=var.replace("'","") # strip extra quotes

# Inputs
infile='STOP_comparison.csv' # from STOP_compare.py
GWVmat1=inputs["MAT1"] # GWV SFR input mat1, from Finalize_SFR.py
GWVmat2=inputs["MAT2"] # GWV SFR input mat2
L1top=inputs["Layer1top"] # Layer 1 (model) TOP elevations

# Outputs
outfile=inputs["MAT1"] # revised GWV SFR input mat1
pdffile='selected_segments.pdf' # plots of selected segments
end_interp_report='fix_w_DEM_interps.txt' # file recording adjustments made using up/dn_increments
error_report='fix_w_DEM_errors.txt' # file for reporting 0 slope errors, etc.

# Settings
float_cutoff=-500 # process all segments with reaches that are float_cutoff above land surface
up_increment=1.0 # if a minimum in the dem is lower than the segment end, new elev. will be end + increment
dn_increment=0.1 # if end+increment is greater than upstream STOP, decrease by dn_increment until current STOP is less
slope_min=1.0e-6 # replace 0 slope values (resulting from NHD) with an aritrary value
enforce_slope_min=True # True or False, if True, replaces all slopes <slope_min with slope_min. Affected slopes are listed in error file.
plot_slope=False # flag to include plots of streambed slope in output PDF
num_segs2plot=50 # number of segments to plot (will choose using regular interval)

print "bringing in segment, reach, streamtop, land surface, and residuals..."
diffs=STOP_compare.stopcomp(L1top,GWVmat1,'STOP_compare_SFR_utilities.csv') # generates column text file from MAT1 and Layer1 top elevs
STOP_compare.plot_diffstats(diffs,'Land-surface - SFR comparison after SFR_utilities.py')
infile_data=np.genfromtxt(infile, delimiter=',',names=True, dtype=None)
float_inds=np.where(infile_data['DIF']>float_cutoff)[0]
floaters=infile_data['MSEG'][float_inds]
segments=list(np.unique(floaters))

print "bringing in reach lengths, slopes and routing info..."
GWVmat1_data=np.genfromtxt(GWVmat1, delimiter=',',names=True,dtype=None)
GWVmat2_data=np.genfromtxt(GWVmat2, delimiter=',',names=True,dtype=None)

# make backup of previous MAT1
GWVmat1old=inputs["MAT1"]+"_old_SFR_utilities"
os.rename(inputs["MAT1"],GWVmat1old)

seg_distdict=defaultdict(list)
reach_lengthsdict=defaultdict(list)
TOPNEWdict=defaultdict(list)
STOP1dict=defaultdict(list)
STOP2dict=defaultdict(list)
slopesdict=defaultdict(list)

# open output file to record segments where STOP was interpolated to end, to avoid going below end
ofp=open(end_interp_report,'w')
ofp2=open(error_report,'w')
print "fixing reaches in segments..."
for segnum in segments:
    
    print str(segnum)
    seg_inds=list(np.where(infile_data['MSEG']==segnum)[0])
    segment=infile_data[seg_inds]
    
    seg_indsGWVmat=list(np.where(GWVmat1_data['segment']==segnum)[0])
    seg_distances=list(np.cumsum(GWVmat1_data['length_in_cell'][seg_indsGWVmat])) # cumulative distances along segment
    reach_lengths=list(GWVmat1_data['length_in_cell'][seg_indsGWVmat])
    NHD_slopes=list(GWVmat1_data['bed_slope'][seg_indsGWVmat])
    
    if len(segment)<3: # no interior points; keep ending elevations
        STOP1dict[segnum] = list(segment['STOP'])
        seg_distdict[segnum] = seg_distances
        reach_lengthsdict[segnum] = reach_lengths
        TOPNEWdict[segnum] = list(segment['TOPNEW'])
        slopesdict[segnum]=NHD_slopes
        continue
    
    # for now, start with existing STOP at RCH1. Could improve by first having a routine that corrects all segment ends to eliminate floaters.
    STOP2=[]
    STOP2.append(segment['STOP'][0]) 
    end=segment['STOP'][-1]
    
    descending=False # If true, means that current reach is lower than all previous
    flattening=False # means that a reach has been encountered with a STOP below the segment end elev.
    if STOP2[0]>segment['TOPNEW'][1]:
        descending=True
    
    for i in range(len(segment['TOPNEW']))[1:-1]:
        dem_el=segment['TOPNEW'][i]
        
        # in case a minimum in DEM is below segment end, set elev. to increment above, interpolate to end
        if segment['TOPNEW'][i]<=end:
            dem_el=end+up_increment
            print "streamtop below segment end elevation,see output file %s " %(end_interp_report)
            ofp.write("STOP at segment %s reach %s is %s, segment end elevation is %s\nReadjusting to %s\n" %(segnum,i,segment['TOPNEW'][i],end,dem_el))
            # check to make sure current STOP elevation is still less than previous
            if dem_el>STOP2[-1]:
                ofp.write("readjusting again to ")
                while dem_el>STOP2[-1]:
                    dem_el-=dn_increment
                    ofp.write(str(dem_el)+'\n')
            ofp.write("interpolating to segment end elevation...\n")   
            nextlower=end
            STOP2.append(dem_el) # first append adjusted elevation, then interp to end
            #slope=(nextlower-STOP2[-1])/(seg_distances[len(segment)-1]-seg_distances[len(STOP2)-1])
            slope=(nextlower-dem_el)/(seg_distances[len(segment)-1]-seg_distances[len(STOP2)])        
            sub_inds=range(len(segment)-1)[len(STOP2):(len(segment)-1)]
            for s in sub_inds:
                dist=seg_distances[s]-seg_distances[s-1]
                STOPint=STOP2[-1]+slope*dist
                STOP2.append(STOPint)
                print '%s %s' %(s,STOPint)
                ofp.write(str(STOPint)+'\n')
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
        elif segment['TOPNEW'][i+1]<=end:
            flattening=True
            continue
        
        # next rch elevation is lower        
        elif segment['TOPNEW'][i+1]<dem_el:
            nextlower=segment['TOPNEW'][i+1]
            if descending and not flattening:
                STOP2.append(dem_el)
            
            # next reach is lower than current but equal to or higher than a previous reach, keep going    
            elif nextlower>=STOP2[-1]:
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
    reach_lengthsdict[segnum] = reach_lengths
    TOPNEWdict[segnum] = list(segment['TOPNEW'])
    STOP1dict[segnum] = list(segment['STOP'])
    
    # update stream bottom slopes for each reach, by interpolating streambed elevations at reach boundaries
    
    # Slope in reach 1 is NHD slope
    # For reaches 2+, interpolate STOP at upstream reach boundary (STOPup), using reach 1 slope
    # Calculate reach slope: STOPup-STOPreach / 0.5* reachlength = reach slope
    # If slope is <0, calculate slope by interpolation between previous reach STOP, and current STOP (because of the above STOP2 algorithm, this value will be positive.
    
    slopes=[]
    STOP2ups=[]
    slopes.append(NHD_slopes[0]) # start with NHD slope for segment
    for reach in range(len(STOP2))[1:]:
        STOP2_up=STOP2[reach-1]-slopes[reach-1]*0.5*reach_lengths[reach-1]
        slope=(STOP2_up-STOP2[reach])/(0.5*reach_lengths[reach])
        
        if enforce_slope_min and slope<slope_min:
            if slope<=0:
                slope=(STOP2[reach-1]-STOP2[reach])/(0.5*reach_lengths[reach-1]+0.5*reach_lengths[reach])
            if slope<slope_min:
                ofp2.write("slope less than slope mininum at segment "+str(segnum)+' reach '+str(reach)+'! Reassigned a value of ' +str(slope_min)+'\n')
                slope=slope_min
        if slope<=0:
            slope=(STOP2[reach-1]-STOP2[reach])/(0.5*reach_lengths[reach-1]+0.5*reach_lengths[reach])
            if slope<=0:
                ofp2.write("zero or negative slope at segment "+str(segnum)+' reach '+str(reach)+'! Reassigned a value of ' +str(slope_min)+'\n')
                # slope is likely 0 because NHD segment end points are equal; assign arbitrary slope
                slope=slope_min
            STOP2_up=STOP2[reach-1]-0.5*reach_lengths[reach-1]*slope
        slopes.append(slope)
        STOP2ups.append(STOP2_up)
    slopesdict[segnum]=slopes

ofp.close()
ofp2.close()
   
print "saving new streamtop elevations and slopes to GWV SFR mat. 1 file"
input_file=open(GWVmat1old).readlines()
ofp=open(outfile,'w')
ofp.write(input_file[0])
for line in input_file[1:]:
    segment=int(line.split(',')[6])
    if segment in STOP2dict:
        reach=int(line.split(',')[5])
        linestart=','.join(map(str,line.split(',')[:3]))
        linemid=','.join(map(str,line.split(',')[5:-2]))
        lineend=line.strip().split(',')[-1]
        STAGE=STOP2dict[segment][reach-1]+1
        STOP=STOP2dict[segment][reach-1]
        slope=slopesdict[segment][reach-1]
        ofp.write('%s,%s,%s,%s,%s,%s\n' %(linestart,STAGE,STOP,linemid,slope,lineend))
    else:
        ofp.write(line)
ofp.close()

# run STOP_compare again, to get results post fix_w_DEM
diffs=STOP_compare.stopcomp(L1top,outfile,'STOP_compare_fix_w_DEM.csv')
STOP_compare.plot_diffstats(diffs,'Land-surface - SFR comparison after fix_w_DEM.py')

# list of segments to plot
#segs2plot=[2,1081,324,1188,1060,1201,1076,341,1113,1131,949,523,408,469,1064,200,300,400,500,600,700,800,900,1000,1100,1200]
segs2plot=segments[::int(np.floor(len(segments)/num_segs2plot))]
segs2plot=sorted(segs2plot)

pdf=PdfPages(pdffile)

print "saving plots of selected segments to " + pdffile
for seg in segs2plot:
    if plot_slope:
        fig, ((ax1, ax2)) = plt.subplots(2,1,sharex=True,sharey=False)
    else:
        fig,ax1=plt.subplots(1,1,1)
    ax1.grid(True)
    p1=ax1.plot(seg_distdict[seg],TOPNEWdict[seg],'b',label='land surface')
    p2=ax1.plot(seg_distdict[seg],STOP1dict[seg],'r',label='sfr_utilities')
    if len(seg_distdict[seg])==len(STOP2dict[seg]):
        p3=ax1.plot(seg_distdict[seg],STOP2dict[seg],'g',label='fix_w_DEM')
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles,labels)
    ax1.set_title('segment ' + str(seg))
    plt.xlabel('distance along segment (ft.)')
    ax1.set_ylabel('Elevation (ft)')
    if plot_slope:
        ax2.grid(True)
        p4=ax2.plot(seg_distdict[seg],slopesdict[seg],color='0.75',label='streambed slopes')
        ax2.set_ylabel('Streambed slope')
        ax3=ax2.twinx()
        p5=ax3.plot(seg_distdict[seg],reach_lengthsdict[seg],'b',label='reach length')
        ax3.set_ylabel('reach length (ft)')
        handles, labels = ax2.get_legend_handles_labels()
        ax2.legend(handles,labels)
        ax3.legend(loc=0) 

    pdf.savefig(fig)
pdf.close()
plt.close('all')

print "finished OK"


