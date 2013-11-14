# Program to correct SFR reach streamtop elevations using land surface elevations, while enforcing downstream monotonicity

import os, shutil
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import STOP_compare
import Fix_segment_ends
import Plot_segment_profiles
import discomb_utilities as dis_util
import pandas as pd # this is just used to sort segments for plotting, prob. can be done in np
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
GWVmat1=inputs["MAT1"] # GWV SFR input mat1, from Finalize_SFR.py
GWVmat2=inputs["MAT2"] # GWV SFR input mat2
L1top=inputs["Layer1top"] # Layer 1 (model) TOP elevations
bots=inputs["Layer_bots"] # Bottom elevations for all layers
MFdis = inputs["MFdis"]

# Outputs
outfile=inputs["MAT1"] # revised GWV SFR input mat1
Mat2out=inputs["MAT2"]
pdffile='selected_segments.pdf' # plots of selected segments
fix_routing_report='Fix_routing_report.txt' # records which segments were given new outsegs by Fix_routing
fix_ends_report='Fixed_segment_ends.csv' # records changes made to segment end elevations
fix_ends_errors='Fix_segment_ends_errors.txt' # records instances where segments were left with backwards routing
end_interp_report='fix_w_DEM_interps.txt' # file recording adjustments made using up/dn_increments
error_report='fix_w_DEM_errors.txt' # file for reporting 0 slope errors, etc.
STOP_comp_SFR_utilities='STOP_compare_SFR_utilities.csv'# from STOP_compare.py
STOP_comp_fixwDEM='STOP_compare_fix_w_DEM.csv' # from STOP_compare.py

# Settings
bed_thickness=3 # SFR streambed thickness
headwater_incise=3 # amount to drop headwater segments below land surface
float_cutoff=-500 # process all segments with reaches that are float_cutoff above land surface
up_increment=1.0 # if a minimum in the dem is lower than the segment end, new elev. will be end + increment
dn_increment=0.1 # if end+increment is greater than upstream STOP, decrease by dn_increment until current STOP is less
slope_min=1.0e-6 # replace 0 slope values (resulting from NHD) with an aritrary value
enforce_slope_min=True # True or False, if True, replaces all slopes <slope_min with slope_min. Affected slopes are listed in error file.
plot_slope=False # flag to include plots of streambed slope in output PDF
num_segs2plot=50 # number of segments to plot (will choose using regular interval)
Fix_ends=False # attempt to adjust segment end elevations to reduce floating and incising
Fix_routing=False # Fixes instance of incorrect routing where two or more streams come together; 
# (all lower-order reaches at junction are routed to segment with lowest maximum (i.e. reach 1) elevation
maxsearch=1 # number of cells to search in, in each direction, to find downstream segments when Fix_routing
backwards_routing_tol=0.1 # amount of backwards routing that will be tolerated (model units) needed because of floating point precision issues with small or flat gradients

print "bringing in segment, reach, streamtop, land surface, and residuals..."
Diffs=STOP_compare.stopcomp(L1top,GWVmat1,STOP_comp_SFR_utilities) # generates column text file from MAT1 and Layer1 top elevs
STOP_compare.plot_diffstats(Diffs['DIF'],'Land-surface - SFR comparison after SFR_utilities.py')
infile_data=np.genfromtxt(STOP_comp_SFR_utilities, delimiter=',',names=True, dtype=None)
float_inds=np.where(infile_data['DIF']>float_cutoff)[0]
floaters=infile_data['MSEG'][float_inds]
segments=list(np.unique(floaters))

print "bringing in reach lengths, slopes and routing info..."
GWVmat1_data=np.genfromtxt(GWVmat1, delimiter=',',names=True,dtype=None)
GWVmat2_data=np.genfromtxt(GWVmat2, delimiter=',',names=True,dtype=None)

# make backups of previous MAT1 and MAT2
GWVmat1old=inputs["MAT1"]+"_old_SFR_utilities"
GWVmat2old=inputs["MAT2"]+"_old"
shutil.move(inputs["MAT1"],GWVmat1old)
shutil.move(inputs["MAT2"],GWVmat2old)


if Fix_routing:
    print "fixing routing at confluences..."
    ofp=open(Mat2out,'w')
    ofp2=open(fix_routing_report,'w')
    ofp.write(','.join(GWVmat2_data.dtype.names)+'\n')
    ofp2.write('segment,old_outseg,new_outseg\n')
    for segment in GWVmat2_data:
        segnum=segment[0]
        outseg_old=segment[2]
        if outseg_old>0:
            # find other SFR cells within maxsearch cells of segment end
            nearbySFR=Fix_segment_ends.find_nearby_SFR(GWVmat1_data,segnum,'last',maxsearch)
    
            # return nearbySFR segment with lowest maximum elevaton
            outseg=Fix_segment_ends.find_outSEG(segnum,GWVmat1_data,nearbySFR,outseg_old)
    
            if outseg<>outseg_old:
                print '%s: changed outseg %s -> %s' %(segnum,outseg_old,outseg)
                newline=map(str,segment.tolist())
                newline=','.join(newline[0:2]+[str(outseg)]+newline[3:])
                ofp.write('%s\n' %(newline))
                ofp2.write('%s,%s,%s\n' %(segnum,outseg_old,outseg))
            else:
                ofp.write(','.join(map(str,segment))+'\n')
        else:
            ofp.write(','.join(map(str,segment))+'\n')
    ofp.close()
    ofp2.close()
    # refresh MAT2 array using new MAT2 file
    GWVmat2_data=np.genfromtxt(Mat2out, delimiter=',',names=True,dtype=None)

upSEGs=Fix_segment_ends.get_upSEGs(GWVmat2_data)
Seg_ends=defaultdict(list)
if Fix_ends:
    print "fixing segment end elevations to reduce floating and incising..."
    nbackwards=0 # number of segments that are routed backwards
    ofp3=open(fix_ends_report,'w') 
    ofp3.write('segment,old_start,old_end,upSEG,dnSEG,up_elev,dn_elev,new_start,new_end\n')
    ofp4=open(fix_ends_errors,'w')
    ofp4.write('These segments were left with backwards routing!\n')
    ofp4.write('segment,upSTOP,dnSTOP,reason\n')

    for segnum in sorted(segments):
        seg_inds=list(np.where(infile_data['MSEG']==segnum)[0])
        segment=infile_data[seg_inds]
        
        if len(segment)<2: # no interior points; keep ending elevations
            
            # run Fix_segment_ends just to get the start elevation, add to Seg_ends dict and move to next
            start,end,nbackwards=Fix_segment_ends.fixends(segnum,segment['STOP'],segment['TOPNEW'],upSEGs,Seg_ends,GWVmat1_data,GWVmat2_data,headwater_incise,ofp3,ofp4,nbackwards,backwards_routing_tol)
            Seg_ends[segnum]=[start]
            continue
        
        if len(segment)>1:
            start,end,nbackwards=Fix_segment_ends.fixends(segnum,segment['STOP'],segment['TOPNEW'],upSEGs,Seg_ends,GWVmat1_data,GWVmat2_data,headwater_incise,ofp3,ofp4,nbackwards,backwards_routing_tol)    
            Seg_ends[segnum]=[start,end]
        
    ofp3.close()
    
    # final check for backwards routing in segment ends
    print "Checking for backwards routing..."
    nbackwards_final,nbackwards_final_file=Fix_segment_ends.check4backwards_routing(Seg_ends,GWVmat2_data,backwards_routing_tol)
    if nbackwards>0:
        print "%s segments had end elevations that were initially backward during correction..." %(nbackwards)
        if nbackwards_final==0:
            print "these were all fixed."
        else:
            print "Warning! %s segments had final end elevations that were backwards!!" %(nbackwards_final)
            print "See %s" %(nbackwards_final_file)

    # Now that segment ends have gone through an initial adjustment and checked for downhill monotonicity;
    # Go through and reset last reach elevations that are a specified amount higher than their downstream segments
    # reopen fix ends report and record any adjustments
    fix_ends_reported=np.genfromtxt(fix_ends_report,delimiter=',',names=True,dtype=None)
    previouslyfixed=fix_ends_reported['segment']
    ofp3=open(fix_ends_report,'w')
    ofp3.write('segment,old_start,old_end,upSEG,dnSEG,up_elev,dn_elev,new_start,new_end_initial,new_end_final\n')
    for segnum in sorted(segments):
        lastSTOP2=Fix_segment_ends.fix_last_reaches(segnum,Seg_ends,GWVmat2_data)
        if segnum in previouslyfixed:
            inds=np.where(fix_ends_reported['segment']==segnum)[0][0]
            ofp3.write(','.join(map(str,fix_ends_reported[inds]))+',%s\n' %(lastSTOP2))
        else:
            ofp3.write('Not previously adjusted'+len(fix_ends_reported[0])*','+'\n')
            
        if len(Seg_ends[segnum])>1:
            Seg_ends[segnum]=[Seg_ends[segnum][0],lastSTOP2]
        else:
            Seg_ends[segnum]=[lastSTOP2]
        

Bottomsdict=STOP_compare.getbottoms(L1top,bots,GWVmat1old)
seg_distdict=defaultdict(list)
reach_lengthsdict=defaultdict(list)
TOPNEWdict=defaultdict(list)
STOP1dict=defaultdict(list)
STOP2dict=defaultdict(list)
slopesdict=defaultdict(list)

# open output file to record segments where STOP was interpolated to end, to avoid going below end
ofp=open(end_interp_report,'w')
ofp2=open(error_report,'w')

print "\nfixing reaches in segments..."
for segnum in sorted(segments):
    
    #print str(segnum)
    seg_inds=list(np.where(infile_data['MSEG']==segnum)[0])
    segment=infile_data[seg_inds]
    
    seg_indsGWVmat=list(np.where(GWVmat1_data['segment']==segnum)[0])
    seg_distances=list(np.cumsum(GWVmat1_data['length_in_cell'][seg_indsGWVmat])) # cumulative distances along segment
    reach_lengths=list(GWVmat1_data['length_in_cell'][seg_indsGWVmat])
    NHD_slopes=list(GWVmat1_data['bed_slope'][seg_indsGWVmat])  
    
    # for segments with one reach
    if len(segment)<2: # no interior points; keep ending elevations
        STOP1dict[segnum] = list(segment['STOP'])
        seg_distdict[segnum] = seg_distances
        reach_lengthsdict[segnum] = reach_lengths
        TOPNEWdict[segnum] = list(segment['TOPNEW'])
        slopesdict[segnum]=NHD_slopes
        
        if not Fix_ends:
            STOP2dict[segnum]=STOP1dict[segnum]
        else:
            STOP2dict[segnum]=Seg_ends[segnum]
        continue
    
    # if segment ends weren't adjusted above, set end elevs at previous values (post-SFR Utilities)
    STOP2=[]
    if not Fix_ends:
        STOP2.append(segment['STOP'][0])
        end=segment['STOP'][-1]
    else:
        STOP2.append(Seg_ends[segnum][0])
        end=Seg_ends[segnum][-1]
    
    
    descending=False # If true, means that current reach is lower than all previous
    flattening=False # means that a reach has been encountered with a STOP below the segment end elev.
    belowbot=False # means that a reach has been encounted with a SBOT below the model bottom
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
            slope=(nextlower-dem_el)/(seg_distances[len(segment)-1]-seg_distances[len(STOP2)-1])        
            sub_inds=range(len(segment)-1)[len(STOP2):(len(segment)-1)]
            for s in sub_inds:
                dist=seg_distances[s]-seg_distances[s-1]
                STOPint=STOP2[-1]+slope*dist
                STOP2.append(STOPint)
                #print '%s %s' %(s,STOPint)
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
            
    #STOP2.append(Seg_ends[segnum][-1])
    STOP2.append(end)
    #STOP2.append(segment['STOP'][-1])
    
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


print "saving new streamtop elevations and slopes to GWV_SFRmat1.dat file"
print "Appending cellnum onto lines"
print "reading in NROW, NCOL, etc from dis file: %s" %MFdis
DX,DY,NLAY,NROW,NCOL,foo = dis_util.read_meta_data(MFdis)
input_file=open(GWVmat1old).readlines()
ofp=open(outfile,'w')
ofp.write('%s,%s\n' %(input_file[0].strip(),'cellnum'))
for line in input_file[1:]:
    line = line.strip().split(',')
    # add cellnum for later evaluation in Arc
    line.append((int(line[1])-1)*NCOL + int(line[0]))
    segment=int(line[6])
    if segment in STOP2dict.keys():
        reach=int(line[5])
        linestart=','.join(map(str,line[:3]))
        linemid=','.join(map(str,line[5:-3]))
        lineend=','.join(map(str,line[-2:]))
        STAGE=STOP2dict[segment][reach-1]+1
        STOP=STOP2dict[segment][reach-1]
        slope=slopesdict[segment][reach-1]
        ofp.write('%s,%s,%s,%s,%s,%s\n' %(linestart,STAGE,STOP,linemid,slope,lineend))
    else:
        ofp.write(line+'\n')
ofp.close()

# run STOP_compare again, to get results post fix_w_DEM
Diffs=STOP_compare.stopcomp(L1top,outfile,STOP_comp_fixwDEM)
STOP_compare.plot_diffstats(Diffs['DIF'],'Land-surface - SFR comparison after fix_w_DEM.py')

# profile plots showing where SFR elev goes below model bottom
#sinkers=Plot_segment_profiles.get_sinkers('below_bot.csv','segment')
#Plot_segment_profiles.plot_profiles(sinkers,seg_distdict,[TOPNEWdict,STOP2dict],['model top','streambed top post-fix_w_DEM'],'Below_bottom.pdf',Bottoms=Bottomsdict)

# profiles of 50 worst floaters
df=pd.DataFrame(Diffs['DIF'],index=Diffs['MSEG'])
Floats=df.sort(columns=0,ascending=True)
floats=list(np.unique(list(Floats[0][0:50].index)))
Plot_segment_profiles.plot_profiles(floats,seg_distdict,[TOPNEWdict,STOP2dict],['model top','streambed top post-fix_w_DEM'],'50_worst_floating.pdf')

# profiles of 50 worst incised
Incised=df.sort(columns=0,ascending=False)
incised=list(np.unique(list(Incised[0][0:50].index)))
Plot_segment_profiles.plot_profiles(incised,seg_distdict,[TOPNEWdict,STOP2dict],['model top','streambed top post-fix_w_DEM'],'50_worst_incised.pdf')

# list of selected segments to plot
segs2plot=segments[::int(np.floor(len(segments)/num_segs2plot))]
segs2plot=sorted(segs2plot)
Plot_segment_profiles.plot_profiles(segs2plot,seg_distdict,[TOPNEWdict,STOP2dict,STOP1dict],['model top','streambed top post-fix_w_DEM','streambed top post-SFR_utilities'],'selected_segments.pdf',plot_slopes=plot_slope,slopes=slopesdict,reach_lengths=reach_lengthsdict)

print "finished OK"

