<<<<<<< HEAD
# SFR_utilities.py
'''
 Utilities to work with SFR output files generated using NHDPlus
 and MODFLOW grid shapefile and correct problems

 SFR Segments are routed from original routine, most of them have
 an elevation consistent with DEMs used to assign elevation to the
 MODFLOW grid cells.  Some segments, however, have elevations at the
 ends that either float above land surface or are incised a large
 amount.  Utility to identify these segments.

 The biggest problem noted in the SFR package was that reach elevations
 were derived by linear interpolation along the stream length from the
 upstream end to the downstream end.  If the land surface does not
 vary linearly, then individual reaches could either float above
 land surface or be too incised.  Utility to correct in-segment
 floaters and in-segment incision.
 
 version for ArcGIS 10 and arcpy module

 Requirements: Python and ArcGIS 10, modules- os, arcpy, numpy, re, defaultdict
 Author: H.W. Reeves; USGS Michigan Water Science Center
 Date: 9/14/2012
'''
import pdb
import arcpy
import os
import numpy as np
import re
from collections import defaultdict

def check_segment(seglist,segment,gwv,sfrtop,invcellnum,maxincise,maxfloating):
    seenflags=dict()
    reachlist=sorted(seglist[segment]) # changed "seg" to segment
    segdiff=[]
    for rch in reachlist:
        flag=0
        cellnum=invcellnum[segment][rch]
        diff=gwv[cellnum]-sfrtop[cellnum]
        segdiff.append(diff)
        if diff < maxfloating:
            flag=1
        if diff > maxincise:
            flag=2
        seenflags[rch]=flag
        
    return (seenflags, segdiff)

def fix_routing(sortedreaches, invcellnum, cortop, corslope, segment):
    #checks routing within each segment, a few slipped by initial processing
    #first take care of ends
    if len(sortedreaches)>3:
        cellnum=invcellnum[segment][sortedreaches[0]]
        cellnext=invcellnum[segment][sortedreaches[1]]
        cellafter=invcellnum[segment][sortedreaches[2]]
        change=cortop[cellnum]-cortop[cellnext]
        if change<0:
            tempdist=(length[cellnum]+length[cellnext])/2.
            cortop[cellnext]=0.5*(cortop[cellafter]+cortop[cellnum])
            if segment==86:
                print cellnum, cellnext, cellafter, '******'
                print cortop[cellnum], cortop[cellnext], cortop[cellafter]
            corslope[cellnum]=(cortop[cellnum]-cortop[cellnext])/tempdist
            tempdist=(length[cellnext]+length[cellafter])/2.
            corslope[cellnext]=(cortop[cellnext]-cortop[cellafter])/tempdist
        cellnum=invcellnum[segment][sortedreaches[-1]]
        cellprev=invcellnum[segment][sortedreaches[-2]]
        cellbefore=invcellnum[segment][sortedreaches[-3]]
        change=cortop[cellprev]-cortop[cellnum]
        if change<0:
            tempdist=(length[cellnum]+length[cellprev])/2.
            cortop[cellprev]=0.5*(cortop[cellbefore]+cortop[cellnum])
            corslope[cellnum]=(cortop[cellprev]-cortop[cellnum])/tempdist
            tempdist=(length[cellprev]+length[cellbefore])/2.
            corslope[cellprev]=(cortop[cellbefore]-cortop[cellprev])/tempdist    
    for i in range(2,len(sortedreaches)-1):
        cellnumprev=invcellnum[segment][sortedreaches[i-1]]
        cellnum=invcellnum[segment][sortedreaches[i]]
        cellnext=invcellnum[segment][sortedreaches[i+1]]
        change=cortop[cellnumprev]-cortop[cellnum]
        if change<0:
            tempdist=(length[cellnum]+length[cellnext])/2.
            cortop[cellnum]=0.5*(cortop[cellnext]+cortop[cellnumprev])
            if segment==86:
                print i, cellnum
                print cortop[cellnumprev], cortop[cellnum], cortop[cellnext]
            corslope[cellnum]=(cortop[cellnum]-cortop[cellnext])/tempdist
        
    return cortop,corslope

def check_routing(sortedreaches, invcellnum, cortop, segment):
    routingflag=0    
    for i in range(0,len(sortedreaches)-1):
        cellnum=invcellnum[segment][sortedreaches[i]]
        cellnext=invcellnum[segment][sortedreaches[i+1]]
        change=cortop[cellnum]-cortop[cellnext]
        if change<0:
            routingflag=1
        
    return routingflag

def fix_incise_interior(segment,sglist,invcellnum,gwv,cortop,corslope,length,maxincise):
    #fix incised reaches within the interior of a segment, this subroutine is
    #not called if one of the two ends has a problem.

    reachlist=sorted(sglist)
    fixedreach=[]
    #build list of node-to-node stream distances
    nodedist=[]
    #first 'node' is the beginning of the segment set distance to zero.
    nodedist.append(0.)
    for i in range(0,len(reachlist)-1):
        cell1=invcellnum[segment][reachlist[i]]
        cell2=invcellnum[segment][reachlist[i+1]]
        tempdist=(length[cell1]+length[cell2])/2.
        nodedist.append(tempdist)
    
    diff=dict()
    totlength=0.
    badreaches=[]
    for reach in reachlist:
        cellnum=invcellnum[segment][reach]
        diff[reach]=gwv[cellnum]-cortop[cellnum]
        totlength+=length[cellnum]
        if diff[reach]>maxincise:     # 'greater than' because incised difference is positive
            badreaches.append(reach)
        
    #start with first badreach, fix SFR elevation and work to last badreach
    #checks all reaches between first and last even if some of the intermediate
    #ones were OK
    if len(badreaches)>0:
        #checks if the problem was already fixed in an earlier step
        for rch in range(badreaches[0],badreaches[-1]+1):
            cellnum=invcellnum[segment][rch]
            prevcell=invcellnum[segment][rch-1]
            nextcell=invcellnum[segment][rch+1]
            if diff[rch]>maxincise:
                target=gwv[cellnum]-maxincise+1.
                interval=cortop[prevcell]-target
##                if segment==8:
##                    print 'before',rch,gwv[cellnum],cortop[cellnum]
                if interval<0:
                    #the correction will wreck routing
                    cortop[cellnum]=cortop[prevcell]
                else:
                    cortop[cellnum]=target
            corslope[cellnum]=(cortop[cellnum]-cortop[nextcell])/nodedist[rch]
##            if segment==8:
##                print 'after',rch,gwv[cellnum],cortop[cellnum]
        
            
    return cortop, corslope

def linear_interp(segment,sglist,invcellnum,cortop,corslope,length):
    #script to linerally interpolate between first and last reach using
    #the stream distance between nodes.  Called when either the first
    #or last reach is moved, so that the starting interpolation is
    #OK (although it may still be floating or incised too much)

    reachlist=sorted(sglist)
    #build list of node-to-node stream distances
    nodedist=[]
    streamdist=[]
    #first 'node' is the beginning of the segment set distance to zero.
    nodedist.append(0.)
    streamdist.append(0.)
    dist=0.
    for i in range(0,len(reachlist)-1):
        cell1=invcellnum[segment][reachlist[i]]
        cell2=invcellnum[segment][reachlist[i+1]]
        tempdist=(length[cell1]+length[cell2])/2.
        dist=dist+tempdist
        nodedist.append(tempdist)
        streamdist.append(dist)
        
    streamlength=streamdist[-1]
    firstcell=invcellnum[segment][reachlist[0]]
    lastcell=invcellnum[segment][reachlist[-1]]
    slp=(cortop[firstcell]-cortop[lastcell])/streamlength
    if slp < 0:
        slp=0.0001
        cortop[firstcell]=cortop[lastcell]+slp*streamlength
    corslope[firstcell]=slp
    corslope[lastcell]=slp
    cells=[]
    cells.append(firstcell)
    for i in range(1,len(reachlist)-1):
        cellnum=invcellnum[segment][reachlist[i]]
        cortop[cellnum]=cortop[firstcell]-slp*streamdist[i]
        cells.append(cellnum)

    cells.append(lastcell)
    plist=[]
    for cell in cells:
        plist.append(cortop[cell])

    return cortop, corslope

def fix_floating_interior(segment,sglist,invcellnum,gwv,cortop,corslope,length,maxfloating):
    #fix floating reaches within the interior of a segment, this subroutine is
    #not called if one of the two ends has a problem.

    reachlist=sorted(sglist)
    fixedreach=[]
    #build list of node-to-node stream distances
    nodedist=[]
    #first 'node' is the beginning of the segment set distance to zero.
    nodedist.append(0.)
    for i in range(0,len(reachlist)-1):
        cell1=invcellnum[segment][reachlist[i]]
        cell2=invcellnum[segment][reachlist[i+1]]
        tempdist=(length[cell1]+length[cell2])/2.
        nodedist.append(tempdist)
    
    diff=dict()
    totlength=0.
    badreaches=[]
    for reach in reachlist:
        cellnum=invcellnum[segment][reach]
        diff[reach]=gwv[cellnum]-cortop[cellnum]
        totlength+=length[cellnum]
        if diff[reach]<maxfloating:     # 'less than' because floating difference is negative
            badreaches.append(reach)
        
    #start with last badreach, fix SFR elevation and work backwards to first badreach
    #updates all reaches between first and last even if some of the intermediate
    #ones were OK
    if segment==1159:
        print 'here I am'
        
    if len(badreaches)>0:
        if badreaches[0]==1:
            badreaches[0]=2
        if badreaches[-1]==reachlist[-1]:
            badreaches[-1]=reachlist[-2]
        #check if the problem was already fixed in an earlier step
        for rch in range(badreaches[-1],badreaches[0]-1,-1):
            cellnum=invcellnum[segment][rch]
            nextcell=invcellnum[segment][rch+1]
            interval=gwv[cellnum]-cortop[nextcell]
            if interval<0:
                #make this reach elevation equal to the next one and live with the floating,
                #this happens if gwv elevation falls below the elevation of the next cell
                cortop[cellnum]=cortop[nextcell]
            else:
                cortop[cellnum]=gwv[cellnum]-0.5*interval
            #check if the new cortop is greater than the previous cortop, if so, split the difference
            #between previous and next...keep the routing
            prevcell=invcellnum[segment][rch-1]
            check=cortop[prevcell]-cortop[cellnum]
            if check<0:
                cortop[cellnum]=0.5*(cortop[prevcell]+cortop[nextcell])
            corslope[cellnum]=(cortop[cellnum]-cortop[nextcell])/nodedist[rch]
            
    return cortop, corslope
#
#  MAIN PROGRAM
#

if __name__ == '__main__':
    
    # Global Input file for SFR utilities
    infile="SFR_setup.in"
    
    #set maximum incision
    maxincise=500
    #set maximum floating, use negative to indicate above land surface
    maxfloating=0.0

    # Get input files locations
    infile=open(infile,'r').readlines()
    inputs=defaultdict()
    for line in infile:
        if "#" not in line.split()[0]:
            varname=line.split("=")[0]
            var=line.split("=")[1].split()[0]
            inputs[varname]=var.replace("'","") # strip extra quotes
    
    MFgrid=inputs["MFgrid"]
    sfrin=inputs["MAT1"]
    ABOVE=open('reachstatus_meanelev.txt','w')
    ABOVE.write('segment,reach,cellnum,gwvelev,sfrtop,flag\n')
    TOPNEW=open(inputs["Layer1top"],'r') # Text file exported from GWV in Matrix format (with all columns on one line selected)
    MAT1corr=sfrin[:-4]+'corr.txt'
    COR=open(MAT1corr,'w')
    CHK2=open('routing_check.txt','w')
    CHK2.write('cellnum,segment,reach,gridelev,sfrelev\n')
    SFR2=open(inputs["MAT2"],'r')
    
    #set rows and columns
    print "getting grid information.."
    table=arcpy.SearchCursor(MFgrid)
    numrow=[]
    numcols=[]
    delx=[]
    dely=[]
    for row in table:
        numrow.append(row.getValue("row"))
        numcols.append(row.getValue("column"))
        delx.append(row.getValue("delx"))
        dely.append(row.getValue("dely"))

    nrow=np.max(numrow)
    ncol=np.max(numcols)

    print "\t%s rows x %s columns" %(nrow,ncol)

    
    # read SFR input from text file
    SFR=open(sfrin,'r')
    header=SFR.readline()
    COR.write(header)
    sfrtop=dict()
    segment=dict()
    reach=dict()
    rchlist=defaultdict(list)
    invcellnum=defaultdict(dict)
    row=dict()
    column=dict()
    width=dict()
    length=dict()
    slope=dict()
    for line in SFR:
        vals=re.split(",",line.strip())
        rowtemp=int(vals[0])
        columntemp=int(vals[1])
        cellnum=(rowtemp-1)*ncol+columntemp
        row[cellnum]=rowtemp
        column[cellnum]=columntemp
        sfrtop[cellnum]=float(vals[4])
        segment[cellnum]=int(vals[6])
        reach[cellnum]=int(vals[5])
        width[cellnum]=float(vals[7])
        length[cellnum]=float(vals[8])
        slope[cellnum]=float(vals[11])
        rchlist[segment[cellnum]].append(reach[cellnum])
        invcellnum[segment[cellnum]][reach[cellnum]]=cellnum
        
    SFR.close()
    ### reach cell elevation from shapefile
    ###set workspace
    ##path=os.getcwd()
    ##arcpy.env.workspace=path
    ##arcpy.env.overwriteOutput=True
    ##gwv=dict()
    ##cells=arcpy.SearchCursor(GWVSHP)
    ##for cell in cells:
    ##    row=int(cell.row)
    ##    column=int(cell.column)
    ##    cellnum=(row-1)*740+column
    ##    gwv[cellnum]=float(cell.top)
    ##del cell, cells
    ##
    ##print 'have shapefile'

    #read cell elevation from text file exported by GWV
    gwv=dict()
    for i in range(0,nrow):
        line=TOPNEW.readline()
        vals=re.split('\s+',line.strip())
        for j in range(0,ncol):
            cellnum=(i*ncol)+(j+1)
            gwv[cellnum]=float(vals[j])
                        
    firstrch=0
    lastrch=0
    bothrch=0
    firstincise=0
    lastincise=0
    bothincise=0
    difference=[]
    segmentstatus=defaultdict(dict)
    numrch=0
    inciserch=0
    numseg=0
    numinciseseg=0
    for seg in sorted(rchlist.keys()):
        (tempseen, segdiff)=check_segment(rchlist,seg,gwv,sfrtop,invcellnum,maxincise,maxfloating)
        #put results in summary dictionary and count up offenders
        difference.extend(segdiff)
        bothflag=0
        bothinciseflag=0
        segflag=0
        seginciseflag=0
        sortedreaches=sorted(tempseen.keys())
        segmentstatus[seg]=tempseen
        for rch in sortedreaches:
            cellnum=invcellnum[seg][rch]
            if segmentstatus[seg][rch]==1:
                segflag+=1
                if rch==1:
                    #first reach is floating
                    firstrch+=1
                    bothflag+=1
                    numrch+=1
                if rch==sortedreaches[-1]:
                    #last reach is floating
                    lastrch+=1
                    bothflag+=1
                    numrch+=1
                if rch>1 and rch<sortedreaches[-1]:
                    #interior reach
                    numrch+=1
            if segmentstatus[seg][rch]==2:
                seginciseflag+=1
                if rch==1: 
                    #first reach is incised
                    firstincise+=1
                    bothinciseflag+=1
                    inciserch+=1
                if rch==sortedreaches[-1]:
                    #last reach is incised
                    lastincise+=1
                    bothinciseflag+=1
                    inciserch+=1
                if rch>1 and rch<sortedreaches[-1]:
                    #interior reach
                    inciserch+=1
            ABOVE.write('{0:d},{1:d},{2:d},{3:.2f},{4:.2f},{5:d}\n'.format(seg,rch,cellnum,gwv[cellnum],sfrtop[cellnum],segmentstatus[seg][rch]))

        if bothflag==2:
            bothrch+=1
            firstrch-=1
            lastrch-=1
        if bothinciseflag==2:
            bothincise+=1
            firstincise-=1
            lastincise-=1
        if segflag>0:
            numseg+=1
        if seginciseflag>0:
            numinciseseg+=1
    ABOVE.close()
    CHK3=open('exteriorcheck.txt','w')
    CHK3.write('cellnum,segment,reach,flag\n')
    for segment in sorted(rchlist.keys()):
        sortedreaches=sorted(segmentstatus[segment].keys())
        if segmentstatus[segment][1]>0 or segmentstatus[segment][sortedreaches[-1]]>0:
            for reach in sortedreaches:
                cellnum=invcellnum[segment][reach]
                if segmentstatus[segment][reach]>0:
                    CHK3.write('{0:d},{1:d},{2:d},{3:d}\n'.format(cellnum,segment,reach,segmentstatus[segment][reach]))
    CHK3.close()
    
    #Now have all the counts and vector of difference
    print 'number of segments with diff < {1:.2f}: {0:d}'.format(numseg,maxfloating)
    print 'number of reaches with diff < {1:.2f}: {0:.0f}'.format(numrch,maxfloating)
    print 'number of first reaches with flag=1: {0:d}'.format(firstrch)
    print 'number of last reaches with flag=1: {0:d}'.format(lastrch)
    print 'number of both first and last reaches with flag=1: {0:d}'.format(bothrch)
    print 'number of segments with diff > {0:.2f}: {1:d}'.format(maxincise, numinciseseg)
    print 'number of reaches with diff > {0:.2f}: {1:.0f}'.format(maxincise,inciserch)
    print 'number of first reaches with flag=2: {0:d}'.format(firstincise)
    print 'number of last reaches with flag=2: {0:d}'.format(lastincise)
    print 'number of both first and last reaches with flag=2: {0:d}'.format(bothincise)

    diffnp=np.array(difference,dtype=float)

    print 'minimum difference {0:f}'.format(np.min(diffnp))
    print 'mean difference {0:f}'.format(np.mean(diffnp))
    print 'maximum difference {0:f}'.format(np.max(diffnp))
    print 'median difference {0:f}'.format(np.median(diffnp))
    print 'standard deviation {0:f}'.format(np.std(diffnp))

    #read in GWV matrix2 file (SEGMENT information) to get output segments
    #for each segment -> inverting identifies headwaters, segments that
    #don't appear as an output segment are headwaters.  Inletsegments
    #dictionary has upstream as key and next segment downstream as
    #value

    outletsegments=dict()
    inletsegments=dict()
    header=SFR2.readline()
    for line in SFR2:
        vals=re.split(',',line)
        segment=int(vals[0])
        outsegment=int(vals[2])
        outletsegments[outsegment]=segment
        inletsegments[segment]=outsegment   
    SFR2.close()

    #loop over segments
    cortop=sfrtop
    corslope=slope
    for segment in sorted(rchlist.keys()):
        sortedreaches=sorted(segmentstatus[segment].keys())
        
        fixSegmentFloat=False
        fixSegmentIncise=False
            
        if segmentstatus[segment][1]>0:
            #first reach has a problem, check if its a headwater, if so
            #fix the first reach, set flag to fix segment
            if not segment in outletsegments:
                if segmentstatus[segment][1]==1:
                    #floating
                    if not segmentstatus[segment][sortedreaches[-1]]>0:
                        fixSegmentFloat=True
                        cellnum=invcellnum[segment][1]
                        cortop[cellnum]=gwv[cellnum]
                        (cortop, corslope)=linear_interp(segment,rchlist[segment],invcellnum,cortop,corslope,length)
                        segmentstatus[segment][1]=0
                if segmentstatus[segment][1]==2:
                    #incised
                    if not segmentstatus[segment][sortedreaches[-1]]>0:
                        fixSegmentIncise=True
                        cellnum=invcellnum[segment][1]
                        cortop[cellnum]=gwv[cellnum]-maxincise
                        (cortop, corslope)=linear_interp(segment,rchlist[segment],invcellnum,cortop,corslope,length)
                        segmentstatus[segment][1]=0

        if segmentstatus[segment][sortedreaches[-1]]==1 and segment in inletsegments:
            if not inletsegments[segment]==0:
                #last reach is floating, check first reach of the next
                #setment and see if there is room to move this reach down
                nextcell=invcellnum[inletsegments[segment]][1]
                cellnum=invcellnum[segment][sortedreaches[-1]]
                if segment==224 or segment==407:
                    print 'segment, next segment',segment,inletsegments[segment]
                    print 'cell, nextcell',cellnum, nextcell
                    print 'gwv[cell],cortop[cell],gwv[next],cortop[next]',gwv[cellnum],cortop[cellnum],gwv[nextcell],cortop[nextcell]
                if cortop[nextcell] < cortop[cellnum]:
                    cortop[cellnum]=cortop[nextcell]+0.1*(cortop[cellnum]-cortop[nextcell])
                    if cortop[cellnum] > gwv[cellnum]:
                         cortop[cellnum]=cortop[nextcell]+0.001*(cortop[cellnum]-cortop[nextcell])
                    if cortop[cellnum] > gwv[cellnum]:
                         segmentstatus[segment][sortedreaches[-1]]=1
                         #fixSegmentFloat=True #might need to move fixSegment Float switch up from line 499
                    else:
                        segmentstatus[segment][sortedreaches[-1]]=0
                    fixSegmentFloat=True
                if segment==224 or segment==407:
                    print 'segment, next segment',segment,inletsegments[segment]
                    print 'cell, nextcell',cellnum, nextcell
                    print 'gwv[cell],cortop[cell],gwv[next],cortop[next]',gwv[cellnum],cortop[cellnum],gwv[nextcell],cortop[nextcell]
        
        #loop over rest of reaches: if any have flag==1 then call fix floating
        #if equal to 2 then call fix incise
        
        if not segmentstatus[segment][1]>0 and not segmentstatus[segment][sortedreaches[-1]]>0:
            for reach in sortedreaches:
                if segmentstatus[segment][reach]==1:
                    fixSegmentFloat=True
                if segmentstatus[segment][reach]==2:
                    fixSegmentIncise=True
       
        if fixSegmentFloat and len(rchlist[segment])>3:
            (cortop,corslope)=fix_floating_interior(segment,rchlist[segment],
                                                        invcellnum,gwv,cortop,corslope,
                                                        length,maxfloating)                                                        
        if fixSegmentIncise and len(rchlist[segment])>3:
            (cortop,corslope)=fix_incise_interior(segment,rchlist[segment],
                                                    invcellnum,gwv,cortop,corslope,
                                                    length,maxincise)    

        (cortop,corslope)=fix_routing(sortedreaches, invcellnum, cortop, corslope, segment)
        routing_flag=check_routing(sortedreaches, invcellnum, cortop, segment)
        if routing_flag>0:
            for reach in sortedreaches:
                cellnum=invcellnum[segment][reach]
                CHK2.write('{0:d},{1:d},{2:d},{3:.2f},{4:.2f}\n'.format(cellnum,segment,reach,gwv[cellnum],cortop[cellnum]))
            
            
        for reach in sortedreaches:
            cellnum=invcellnum[segment][reach]
            printlist=(row[cellnum],
                       column[cellnum],
                       1,
                       cortop[cellnum]+1.,
                       cortop[cellnum],
                       reach,
                       segment,
                       width[cellnum],
                       length[cellnum],
                       5.,
                       1.,
                       corslope[cellnum],
                       0.0370,
                       )
            COR.write(",".join(map(str,printlist))+'\n')
    
    COR.close()
    os.rename(sfrin,sfrin+'_old')
    os.rename(MAT1corr,inputs["MAT1"])
    CHK2.close()
    #end of main



=======
# SFR_utilities.py
'''
 Utilities to work with SFR output files generated using NHDPlus
 and MODFLOW grid shapefile and correct problems

 SFR Segments are routed from original routine, most of them have
 an elevation consistent with DEMs used to assign elevation to the
 MODFLOW grid cells.  Some segments, however, have elevations at the
 ends that either float above land surface or are incised a large
 amount.  Utility to identify these segments.

 The biggest problem noted in the SFR package was that reach elevations
 were derived by linear interpolation along the stream length from the
 upstream end to the downstream end.  If the land surface does not
 vary linearly, then individual reaches could either float above
 land surface or be too incised.  Utility to correct in-segment
 floaters and in-segment incision.
 
 version for ArcGIS 10 and arcpy module

 Requirements: Python and ArcGIS 10, modules- os, arcpy, numpy, re, defaultdict
 Author: H.W. Reeves; USGS Michigan Water Science Center
 Date: 9/14/2012
'''
import pdb
import arcpy
import os
import numpy as np
import re
from collections import defaultdict

def check_segment(seglist,segment,gwv,sfrtop,invcellnum,maxincise,maxfloating):
    seenflags=dict()
    reachlist=sorted(seglist[segment]) # changed "seg" to segment
    segdiff=[]
    for rch in reachlist:
        flag=0
        cellnum=invcellnum[segment][rch]
        diff=gwv[cellnum]-sfrtop[cellnum]
        segdiff.append(diff)
        if diff < maxfloating:
            flag=1
        if diff > maxincise:
            flag=2
        seenflags[rch]=flag
        
    return (seenflags, segdiff)

def fix_routing(sortedreaches, invcellnum, cortop, corslope, segment):
    #checks routing within each segment, a few slipped by initial processing
    #first take care of ends
    if len(sortedreaches)>3:
        cellnum=invcellnum[segment][sortedreaches[0]]
        cellnext=invcellnum[segment][sortedreaches[1]]
        cellafter=invcellnum[segment][sortedreaches[2]]
        change=cortop[cellnum]-cortop[cellnext]
        if change<0:
            tempdist=(length[cellnum]+length[cellnext])/2.
            cortop[cellnext]=0.5*(cortop[cellafter]+cortop[cellnum])
            if segment==86:
                print cellnum, cellnext, cellafter, '******'
                print cortop[cellnum], cortop[cellnext], cortop[cellafter]
            corslope[cellnum]=(cortop[cellnum]-cortop[cellnext])/tempdist
            tempdist=(length[cellnext]+length[cellafter])/2.
            corslope[cellnext]=(cortop[cellnext]-cortop[cellafter])/tempdist
        cellnum=invcellnum[segment][sortedreaches[-1]]
        cellprev=invcellnum[segment][sortedreaches[-2]]
        cellbefore=invcellnum[segment][sortedreaches[-3]]
        change=cortop[cellprev]-cortop[cellnum]
        if change<0:
            tempdist=(length[cellnum]+length[cellprev])/2.
            cortop[cellprev]=0.5*(cortop[cellbefore]+cortop[cellnum])
            corslope[cellnum]=(cortop[cellprev]-cortop[cellnum])/tempdist
            tempdist=(length[cellprev]+length[cellbefore])/2.
            corslope[cellprev]=(cortop[cellbefore]-cortop[cellprev])/tempdist    
    for i in range(2,len(sortedreaches)-1):
        cellnumprev=invcellnum[segment][sortedreaches[i-1]]
        cellnum=invcellnum[segment][sortedreaches[i]]
        cellnext=invcellnum[segment][sortedreaches[i+1]]
        change=cortop[cellnumprev]-cortop[cellnum]
        if change<0:
            tempdist=(length[cellnum]+length[cellnext])/2.
            cortop[cellnum]=0.5*(cortop[cellnext]+cortop[cellnumprev])
            if segment==86:
                print i, cellnum
                print cortop[cellnumprev], cortop[cellnum], cortop[cellnext]
            corslope[cellnum]=(cortop[cellnum]-cortop[cellnext])/tempdist
        
    return cortop,corslope

def check_routing(sortedreaches, invcellnum, cortop, segment):
    routingflag=0    
    for i in range(0,len(sortedreaches)-1):
        cellnum=invcellnum[segment][sortedreaches[i]]
        cellnext=invcellnum[segment][sortedreaches[i+1]]
        change=cortop[cellnum]-cortop[cellnext]
        if change<0:
            routingflag=1
        
    return routingflag

def fix_incise_interior(segment,sglist,invcellnum,gwv,cortop,corslope,length,maxincise):
    #fix incised reaches within the interior of a segment, this subroutine is
    #not called if one of the two ends has a problem.

    reachlist=sorted(sglist)
    fixedreach=[]
    #build list of node-to-node stream distances
    nodedist=[]
    #first 'node' is the beginning of the segment set distance to zero.
    nodedist.append(0.)
    for i in range(0,len(reachlist)-1):
        cell1=invcellnum[segment][reachlist[i]]
        cell2=invcellnum[segment][reachlist[i+1]]
        tempdist=(length[cell1]+length[cell2])/2.
        nodedist.append(tempdist)
    
    diff=dict()
    totlength=0.
    badreaches=[]
    for reach in reachlist:
        cellnum=invcellnum[segment][reach]
        diff[reach]=gwv[cellnum]-cortop[cellnum]
        totlength+=length[cellnum]
        if diff[reach]>maxincise:     # 'greater than' because incised difference is positive
            badreaches.append(reach)
        
    #start with first badreach, fix SFR elevation and work to last badreach
    #checks all reaches between first and last even if some of the intermediate
    #ones were OK
    if len(badreaches)>0:
        #checks if the problem was already fixed in an earlier step
        for rch in range(badreaches[0],badreaches[-1]+1):
            cellnum=invcellnum[segment][rch]
            prevcell=invcellnum[segment][rch-1]
            nextcell=invcellnum[segment][rch+1]
            if diff[rch]>maxincise:
                target=gwv[cellnum]-maxincise+1.
                interval=cortop[prevcell]-target
##                if segment==8:
##                    print 'before',rch,gwv[cellnum],cortop[cellnum]
                if interval<0:
                    #the correction will wreck routing
                    cortop[cellnum]=cortop[prevcell]
                else:
                    cortop[cellnum]=target
            corslope[cellnum]=(cortop[cellnum]-cortop[nextcell])/nodedist[rch]
##            if segment==8:
##                print 'after',rch,gwv[cellnum],cortop[cellnum]
        
            
    return cortop, corslope

def linear_interp(segment,sglist,invcellnum,cortop,corslope,length):
    #script to linerally interpolate between first and last reach using
    #the stream distance between nodes.  Called when either the first
    #or last reach is moved, so that the starting interpolation is
    #OK (although it may still be floating or incised too much)

    reachlist=sorted(sglist)
    #build list of node-to-node stream distances
    nodedist=[]
    streamdist=[]
    #first 'node' is the beginning of the segment set distance to zero.
    nodedist.append(0.)
    streamdist.append(0.)
    dist=0.
    for i in range(0,len(reachlist)-1):
        cell1=invcellnum[segment][reachlist[i]]
        cell2=invcellnum[segment][reachlist[i+1]]
        tempdist=(length[cell1]+length[cell2])/2.
        dist=dist+tempdist
        nodedist.append(tempdist)
        streamdist.append(dist)
        
    streamlength=streamdist[-1]
    firstcell=invcellnum[segment][reachlist[0]]
    lastcell=invcellnum[segment][reachlist[-1]]
    slp=(cortop[firstcell]-cortop[lastcell])/streamlength
    if slp < 0:
        slp=0.0001
        cortop[firstcell]=cortop[lastcell]+slp*streamlength
    corslope[firstcell]=slp
    corslope[lastcell]=slp
    cells=[]
    cells.append(firstcell)
    for i in range(1,len(reachlist)-1):
        cellnum=invcellnum[segment][reachlist[i]]
        cortop[cellnum]=cortop[firstcell]-slp*streamdist[i]
        cells.append(cellnum)

    cells.append(lastcell)
    plist=[]
    for cell in cells:
        plist.append(cortop[cell])

    return cortop, corslope

def fix_floating_interior(segment,sglist,invcellnum,gwv,cortop,corslope,length,maxfloating):
    #fix floating reaches within the interior of a segment, this subroutine is
    #not called if one of the two ends has a problem.

    reachlist=sorted(sglist)
    fixedreach=[]
    #build list of node-to-node stream distances
    nodedist=[]
    #first 'node' is the beginning of the segment set distance to zero.
    nodedist.append(0.)
    for i in range(0,len(reachlist)-1):
        cell1=invcellnum[segment][reachlist[i]]
        cell2=invcellnum[segment][reachlist[i+1]]
        tempdist=(length[cell1]+length[cell2])/2.
        nodedist.append(tempdist)
    
    diff=dict()
    totlength=0.
    badreaches=[]
    for reach in reachlist:
        cellnum=invcellnum[segment][reach]
        diff[reach]=gwv[cellnum]-cortop[cellnum]
        totlength+=length[cellnum]
        if diff[reach]<maxfloating:     # 'less than' because floating difference is negative
            badreaches.append(reach)
        
    #start with last badreach, fix SFR elevation and work backwards to first badreach
    #updates all reaches between first and last even if some of the intermediate
    #ones were OK
    if segment==1159:
        print 'here I am'
        
    if len(badreaches)>0:
        if badreaches[0]==1:
            badreaches[0]=2
        if badreaches[-1]==reachlist[-1]:
            badreaches[-1]=reachlist[-2]
        #check if the problem was already fixed in an earlier step
        for rch in range(badreaches[-1],badreaches[0]-1,-1):
            cellnum=invcellnum[segment][rch]
            nextcell=invcellnum[segment][rch+1]
            interval=gwv[cellnum]-cortop[nextcell]
            if interval<0:
                #make this reach elevation equal to the next one and live with the floating,
                #this happens if gwv elevation falls below the elevation of the next cell
                cortop[cellnum]=cortop[nextcell]
            else:
                cortop[cellnum]=gwv[cellnum]-0.5*interval
            #check if the new cortop is greater than the previous cortop, if so, split the difference
            #between previous and next...keep the routing
            prevcell=invcellnum[segment][rch-1]
            check=cortop[prevcell]-cortop[cellnum]
            if check<0:
                cortop[cellnum]=0.5*(cortop[prevcell]+cortop[nextcell])
            corslope[cellnum]=(cortop[cellnum]-cortop[nextcell])/nodedist[rch]
            
    return cortop, corslope
#
#  MAIN PROGRAM
#

if __name__ == '__main__':
    
    # Global Input file for SFR utilities
    infile="SFR_setup.in"
    
    #set maximum incision
    maxincise=500
    #set maximum floating, use negative to indicate above land surface
    maxfloating=0.0

    # Get input files locations
    infile=open(infile,'r').readlines()
    inputs=defaultdict()
    for line in infile:
        if "#" not in line.split()[0]:
            varname=line.split("=")[0]
            var=line.split("=")[1].split()[0]
            inputs[varname]=var.replace("'","") # strip extra quotes
    
    MFgrid=inputs["MFgrid"]
    sfrin=inputs["MAT1"]
    ABOVE=open('reachstatus_meanelev.txt','w')
    ABOVE.write('segment,reach,cellnum,gwvelev,sfrtop,flag\n')
    TOPNEW=open(inputs["Layer1top"],'r') # Text file exported from GWV in Matrix format (with all columns on one line selected)
    MAT1corr=sfrin[:-4]+'corr.txt'
    COR=open(MAT1corr,'w')
    CHK2=open('routing_check.txt','w')
    CHK2.write('cellnum,segment,reach,gridelev,sfrelev\n')
    SFR2=open(inputs["MAT2"],'r')
    
    #set rows and columns
    print "getting grid information.."
    table=arcpy.SearchCursor(MFgrid)
    numrow=[]
    numcols=[]
    delx=[]
    dely=[]
    for row in table:
        numrow.append(row.getValue("row"))
        numcols.append(row.getValue("column"))
        delx.append(row.getValue("delx"))
        dely.append(row.getValue("dely"))

    nrow=np.max(numrow)
    ncol=np.max(numcols)

    print "\t%s rows x %s columns" %(nrow,ncol)

    
    # read SFR input from text file
    SFR=open(sfrin,'r')
    header=SFR.readline()
    COR.write(header)
    sfrtop=dict()
    segment=dict()
    reach=dict()
    rchlist=defaultdict(list)
    invcellnum=defaultdict(dict)
    row=dict()
    column=dict()
    width=dict()
    length=dict()
    slope=dict()
    for line in SFR:
        vals=re.split(",",line.strip())
        rowtemp=int(vals[0])
        columntemp=int(vals[1])
        cellnum=(rowtemp-1)*ncol+columntemp
        row[cellnum]=rowtemp
        column[cellnum]=columntemp
        sfrtop[cellnum]=float(vals[4])
        segment[cellnum]=int(vals[6])
        reach[cellnum]=int(vals[5])
        width[cellnum]=float(vals[7])
        length[cellnum]=float(vals[8])
        slope[cellnum]=float(vals[11])
        rchlist[segment[cellnum]].append(reach[cellnum])
        invcellnum[segment[cellnum]][reach[cellnum]]=cellnum
        
    SFR.close()
    ### reach cell elevation from shapefile
    ###set workspace
    ##path=os.getcwd()
    ##arcpy.env.workspace=path
    ##arcpy.env.overwriteOutput=True
    ##gwv=dict()
    ##cells=arcpy.SearchCursor(GWVSHP)
    ##for cell in cells:
    ##    row=int(cell.row)
    ##    column=int(cell.column)
    ##    cellnum=(row-1)*740+column
    ##    gwv[cellnum]=float(cell.top)
    ##del cell, cells
    ##
    ##print 'have shapefile'

    #read cell elevation from text file exported by GWV
    gwv=dict()
    for i in range(0,nrow):
        line=TOPNEW.readline()
        vals=re.split('\s+',line.strip())
        for j in range(0,ncol):
            cellnum=(i*ncol)+(j+1)
            gwv[cellnum]=float(vals[j])
                        
    firstrch=0
    lastrch=0
    bothrch=0
    firstincise=0
    lastincise=0
    bothincise=0
    difference=[]
    segmentstatus=defaultdict(dict)
    numrch=0
    inciserch=0
    numseg=0
    numinciseseg=0
    for seg in sorted(rchlist.keys()):
        (tempseen, segdiff)=check_segment(rchlist,seg,gwv,sfrtop,invcellnum,maxincise,maxfloating)
        #put results in summary dictionary and count up offenders
        difference.extend(segdiff)
        bothflag=0
        bothinciseflag=0
        segflag=0
        seginciseflag=0
        sortedreaches=sorted(tempseen.keys())
        segmentstatus[seg]=tempseen
        for rch in sortedreaches:
            cellnum=invcellnum[seg][rch]
            if segmentstatus[seg][rch]==1:
                segflag+=1
                if rch==1:
                    #first reach is floating
                    firstrch+=1
                    bothflag+=1
                    numrch+=1
                if rch==sortedreaches[-1]:
                    #last reach is floating
                    lastrch+=1
                    bothflag+=1
                    numrch+=1
                if rch>1 and rch<sortedreaches[-1]:
                    #interior reach
                    numrch+=1
            if segmentstatus[seg][rch]==2:
                seginciseflag+=1
                if rch==1: 
                    #first reach is incised
                    firstincise+=1
                    bothinciseflag+=1
                    inciserch+=1
                if rch==sortedreaches[-1]:
                    #last reach is incised
                    lastincise+=1
                    bothinciseflag+=1
                    inciserch+=1
                if rch>1 and rch<sortedreaches[-1]:
                    #interior reach
                    inciserch+=1
            ABOVE.write('{0:d},{1:d},{2:d},{3:.2f},{4:.2f},{5:d}\n'.format(seg,rch,cellnum,gwv[cellnum],sfrtop[cellnum],segmentstatus[seg][rch]))

        if bothflag==2:
            bothrch+=1
            firstrch-=1
            lastrch-=1
        if bothinciseflag==2:
            bothincise+=1
            firstincise-=1
            lastincise-=1
        if segflag>0:
            numseg+=1
        if seginciseflag>0:
            numinciseseg+=1
    ABOVE.close()
    CHK3=open('exteriorcheck.txt','w')
    CHK3.write('cellnum,segment,reach,flag\n')
    for segment in sorted(rchlist.keys()):
        sortedreaches=sorted(segmentstatus[segment].keys())
        if segmentstatus[segment][1]>0 or segmentstatus[segment][sortedreaches[-1]]>0:
            for reach in sortedreaches:
                cellnum=invcellnum[segment][reach]
                if segmentstatus[segment][reach]>0:
                    CHK3.write('{0:d},{1:d},{2:d},{3:d}\n'.format(cellnum,segment,reach,segmentstatus[segment][reach]))
    CHK3.close()
    
    #Now have all the counts and vector of difference
    print 'number of segments with diff < {1:.2f}: {0:d}'.format(numseg,maxfloating)
    print 'number of reaches with diff < {1:.2f}: {0:.0f}'.format(numrch,maxfloating)
    print 'number of first reaches with flag=1: {0:d}'.format(firstrch)
    print 'number of last reaches with flag=1: {0:d}'.format(lastrch)
    print 'number of both first and last reaches with flag=1: {0:d}'.format(bothrch)
    print 'number of segments with diff > {0:.2f}: {1:d}'.format(maxincise, numinciseseg)
    print 'number of reaches with diff > {0:.2f}: {1:.0f}'.format(maxincise,inciserch)
    print 'number of first reaches with flag=2: {0:d}'.format(firstincise)
    print 'number of last reaches with flag=2: {0:d}'.format(lastincise)
    print 'number of both first and last reaches with flag=2: {0:d}'.format(bothincise)

    diffnp=np.array(difference,dtype=float)

    print 'minimum difference {0:f}'.format(np.min(diffnp))
    print 'mean difference {0:f}'.format(np.mean(diffnp))
    print 'maximum difference {0:f}'.format(np.max(diffnp))
    print 'median difference {0:f}'.format(np.median(diffnp))
    print 'standard deviation {0:f}'.format(np.std(diffnp))

    #read in GWV matrix2 file (SEGMENT information) to get output segments
    #for each segment -> inverting identifies headwaters, segments that
    #don't appear as an output segment are headwaters.  Inletsegments
    #dictionary has upstream as key and next segment downstream as
    #value

    outletsegments=dict()
    inletsegments=dict()
    header=SFR2.readline()
    for line in SFR2:
        vals=re.split(',',line)
        segment=int(vals[0])
        outsegment=int(vals[2])
        outletsegments[outsegment]=segment
        inletsegments[segment]=outsegment   
    SFR2.close()

    #loop over segments
    cortop=sfrtop
    corslope=slope
    for segment in sorted(rchlist.keys()):
        sortedreaches=sorted(segmentstatus[segment].keys())
        
        fixSegmentFloat=False
        fixSegmentIncise=False
            
        if segmentstatus[segment][1]>0:
            #first reach has a problem, check if its a headwater, if so
            #fix the first reach, set flag to fix segment
            if not segment in outletsegments:
                if segmentstatus[segment][1]==1:
                    #floating
                    if not segmentstatus[segment][sortedreaches[-1]]>0:
                        fixSegmentFloat=True
                        cellnum=invcellnum[segment][1]
                        cortop[cellnum]=gwv[cellnum]
                        (cortop, corslope)=linear_interp(segment,rchlist[segment],invcellnum,cortop,corslope,length)
                        segmentstatus[segment][1]=0
                if segmentstatus[segment][1]==2:
                    #incised
                    if not segmentstatus[segment][sortedreaches[-1]]>0:
                        fixSegmentIncise=True
                        cellnum=invcellnum[segment][1]
                        cortop[cellnum]=gwv[cellnum]-maxincise
                        (cortop, corslope)=linear_interp(segment,rchlist[segment],invcellnum,cortop,corslope,length)
                        segmentstatus[segment][1]=0

        if segmentstatus[segment][sortedreaches[-1]]==1 and segment in inletsegments:
            if not inletsegments[segment]==0:
                #last reach is floating, check first reach of the next
                #setment and see if there is room to move this reach down
                nextcell=invcellnum[inletsegments[segment]][1]
                cellnum=invcellnum[segment][sortedreaches[-1]]
                if segment==224 or segment==407:
                    print 'segment, next segment',segment,inletsegments[segment]
                    print 'cell, nextcell',cellnum, nextcell
                    print 'gwv[cell],cortop[cell],gwv[next],cortop[next]',gwv[cellnum],cortop[cellnum],gwv[nextcell],cortop[nextcell]
                if cortop[nextcell] < cortop[cellnum]:
                    cortop[cellnum]=cortop[nextcell]+0.1*(cortop[cellnum]-cortop[nextcell])
                    if cortop[cellnum] > gwv[cellnum]:
                         cortop[cellnum]=cortop[nextcell]+0.001*(cortop[cellnum]-cortop[nextcell])
                    if cortop[cellnum] > gwv[cellnum]:
                         segmentstatus[segment][sortedreaches[-1]]=1
                         #fixSegmentFloat=True #might need to move fixSegment Float switch up from line 499
                    else:
                        segmentstatus[segment][sortedreaches[-1]]=0
                    fixSegmentFloat=True
                if segment==224 or segment==407:
                    print 'segment, next segment',segment,inletsegments[segment]
                    print 'cell, nextcell',cellnum, nextcell
                    print 'gwv[cell],cortop[cell],gwv[next],cortop[next]',gwv[cellnum],cortop[cellnum],gwv[nextcell],cortop[nextcell]
        
        #loop over rest of reaches: if any have flag==1 then call fix floating
        #if equal to 2 then call fix incise
        
        if not segmentstatus[segment][1]>0 and not segmentstatus[segment][sortedreaches[-1]]>0:
            for reach in sortedreaches:
                if segmentstatus[segment][reach]==1:
                    fixSegmentFloat=True
                if segmentstatus[segment][reach]==2:
                    fixSegmentIncise=True
       
        if fixSegmentFloat and len(rchlist[segment])>3:
            (cortop,corslope)=fix_floating_interior(segment,rchlist[segment],
                                                        invcellnum,gwv,cortop,corslope,
                                                        length,maxfloating)                                                        
        if fixSegmentIncise and len(rchlist[segment])>3:
            (cortop,corslope)=fix_incise_interior(segment,rchlist[segment],
                                                    invcellnum,gwv,cortop,corslope,
                                                    length,maxincise)    

        (cortop,corslope)=fix_routing(sortedreaches, invcellnum, cortop, corslope, segment)
        routing_flag=check_routing(sortedreaches, invcellnum, cortop, segment)
        if routing_flag>0:
            for reach in sortedreaches:
                cellnum=invcellnum[segment][reach]
                CHK2.write('{0:d},{1:d},{2:d},{3:.2f},{4:.2f}\n'.format(cellnum,segment,reach,gwv[cellnum],cortop[cellnum]))
            
            
        for reach in sortedreaches:
            cellnum=invcellnum[segment][reach]
            printlist=(row[cellnum],
                       column[cellnum],
                       1,
                       cortop[cellnum]+1.,
                       cortop[cellnum],
                       reach,
                       segment,
                       width[cellnum],
                       length[cellnum],
                       5.,
                       1.,
                       corslope[cellnum],
                       0.0370,
                       )
            COR.write(",".join(map(str,printlist))+'\n')
    
    COR.close()
    os.rename(sfrin,sfrin+'_old')
    os.rename(MAT1corr,inputs["MAT1"])
    CHK2.close()
    #end of main



>>>>>>> parent of e9ab5ff... Minor bugfixes
