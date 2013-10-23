# RouteRiverCells.py
# Description: Makes a to/from table of river cells by cellnumber
# to help build and check river routing
#
# also uses hydroseq from NHDPlus to build the SFR Segment/reach
# ordering table which is printed to a file
#
# Requirements: os, sys, re, arcgisscripting, defaultdict
#
# Author: H.W. Reeves; USGS Michigan Water Science Center
# Date 2/3/2010
#

import os
import sys
import re
import arcgisscripting
from collections import defaultdict
import math
from operator import itemgetter
import pdb

# Global Input file for SFR utilities
infile="SFR_setup.in"

# Settings
#set elevation difference to assign same value to all cells
#so that check does not get messed up by rounding floats
eps=1.0e-02

# Get input files locations
infile=open(infile,'r').readlines()
inputs=defaultdict()
for line in infile:
    if "#" not in line.split()[0]:
        varname=line.split("=")[0]
        var=line.split("=")[1].split()[0]
        inputs[varname]=var.replace("'","") # strip extra quotes

CELLS=inputs["CELLS"]
INTERSECT=inputs["INTERSECT"]
NHD=inputs["NHD"]  # from intersect.py
VAA=inputs["PlusflowVAA"]

try:
    OUT=open("routed_cells.txt",'w')
except:
    print "could not open output file"
    exit()

OUT.write("fromcell, tocell\n")



path=os.getcwd()
gp=arcgisscripting.create()
gp.Workspace=path

#make a dictionary of hydroseq, uphydroseq, dnhydroseq numbers by COMID
hydroseq=dict()
uphydroseq=dict()
dnhydroseq=dict()
segments=gp.SearchCursor(VAA)
segment=segments.Next()
while segment:
    comid=int(segment.comid)
    hydroseq[comid]=int(segment.hydroseq)
    uphydroseq[comid]=int(segment.uphydroseq)
    dnhydroseq[comid]=int(segment.dnhydroseq)
    segment=segments.Next()

del segment, segments

print 'have dictionary'

#delete any working layers if found
if gp.Exists("temp_lyr"):
    gp.Delete_management("temp_lyr")
if gp.Exists("temp_nhd_lyr"):
    gp.Delete_management("temp_nhd_lyr")
gp.overwriteoutput = 1

#read in from->to table as a dictionary of lists
#generated in RouteStreamNetwork.py
    
fromCOMIDlist=defaultdict(list)
toCOMIDlist=defaultdict(list)

#fromCOMIDlist - key is the 'from', value is the 'to'
#toCOMIDlist - key is the 'to', value is the 'from'

#need to correct routing lists for segments
#that intersected the boundary
#read in the file with information for clipped segments at the boundary
print "Processing clipped segments along domain boundary..."
CLIP=open("boundaryClipsRouting.txt",'r').readlines()
header = CLIP.pop(0) #header
clipto=dict()
clipfrom=dict()
for line in CLIP:
    vals=line.strip().split(',')
    fromin=int(float(vals[0]))
    toin=int(float(vals[1]))
    if fromin==99999:
        clipto[toin]=1
    if toin==99999:
        clipfrom[fromin]=1


infile=open("check_network.txt",'r')
line=infile.readline()
while line:
    vals=re.split(",",line)
    fromin=int(vals[0])
    toin=int(vals[1])
    if not toin in clipto:
        toCOMIDlist[toin].append(fromin)
        if not fromin in clipfrom:
            fromCOMIDlist[fromin].append(toin)
        else:
            fromCOMIDlist[fromin].append(99999)
    line=infile.readline()

infile.close()



#now read in segment level information to order cells for each comid
#using logic from AssignRiverElev.py
count=0
RCH=open('reach_ordering.txt','w')
RCH.write('CELLNUM,COMID, hydroseq, uphydroseq, dnhydroseq, SFRseq, localseq\n')

SFRseq=0
hydroseqused=dict()
#get a list of the comids/hydroseq in this model
for comid in fromCOMIDlist.iterkeys():
    if comid in hydroseq:
        hydroseqused[comid]=hydroseq[comid]
    else:
        exit('check VAA table for COMID %d' % comid)

del hydroseq    
#sort from high hydroseq to low hydroseq and then
#use the sorted list to loop over comids

sortedtuples=hydroseqused.items()
sortedtuples.sort(key=itemgetter(1), reverse=True)
sortedlist=map(itemgetter(0),sortedtuples)

icount=0
needfix = 0 # flag regarding whether or not manual fixes will be needed
outfixfilename = 'COMIDS_to_fixRRC.dat'
ofpfix = open(outfixfilename,'w')
ofpfix.write('The following COMIDs require manual attention.\nFix them in check_network.txt.\n')
ofpfix.write('Then rerun RouteRiverCells.py.\n' + '#'*42 + '\n')
for comid in sortedlist:
    icount+=1
    SFRseq=SFRseq+1
    print SFRseq
    count=count+1
    #find all the intersected cells with the same comid
    theQuery="COMID= %d" %(comid)
    print theQuery
    gp.MakeFeatureLayer(INTERSECT,"temp_lyr",theQuery)
    #make hashes of starting and ending x,y and
    #then find the start without a corresponding end
    #and the end without the corresponding start:  order
    #the stream from start to end..
    shortsegs=gp.SearchCursor("temp_lyr")
    seg=shortsegs.Next()
    count=0
    startx=defaultdict(list)
    starty=defaultdict(list)
    endx=defaultdict(list)
    endy=defaultdict(list)
    top=dict()
    cellnumber=dict()
    while seg:
        count=count+1
        startx[seg.FID].append(seg.X_start)
        starty[seg.FID].append(seg.Y_start)
        endx[seg.FID].append(seg.X_end)
        endy[seg.FID].append(seg.Y_end)
        top_elev=float(seg.MAXELEVSMO)*3.28/100.
        top[seg.FID]=top_elev
        cellnumber[seg.FID]=seg.CELLNUM
        seg=shortsegs.Next()
    del shortsegs
    orderedsegs=[]
    start_has_end=dict()
    end_has_start=dict()
    fact=0.95  # because small segments were dropped
    for test in startx.iterkeys():
        j=0
        while j < len(startx[test]):
            for othertest in startx.iterkeys():
                if(othertest == test):
                    continue
                i=0
                diffx=startx[test][j]-endx[othertest][i]
                diffy=starty[test][j]-endy[othertest][i]
                while i < len(startx[othertest]):
                    if(math.fabs(diffx)<fact and math.fabs(diffy)<fact):
                        start_has_end[test]=othertest
                    i=i+1
                i=0
                diffx=endx[test][j]-startx[othertest][i]
                diffy=endy[test][j]-starty[othertest][i]
                while i < len(startx[othertest]):
                    if(math.fabs(diffx)<fact and math.fabs(diffy)<fact):
                        end_has_start[test]=othertest
                    i=i+1
            j=j+1
    
    #find key in start that didn't match an end and key in end that didn't match a start
    numstart=0
    numend=0
    startingFID=[]
    endingFID=[]
    for test in startx.iterkeys():
        if(test not in start_has_end):
            startingFID.append(test)
            numstart=numstart+1
        if(test not in end_has_start):
            endingFID.append(test)
            numend=numend+1

    if (numstart != 1 or numend != 1):
        print "starting FIDs: "+",".join(map(str,startingFID))
        print "ending FIDs: " + ",".join(map(str,endingFID))
        gp.refreshcatalog(path)
        print "manually fix COMID =" + str(comid)
        ofpfix.write("manually fix COMID = %d\n"  %int(comid))
        needfix +=1
    else:
        # MNF DEBUG ----- BIG ASSUMPTION! 
        # Seems to me that this loop needed indentation to avoid being trapped by error above. Maybe be way off!
        orderedFID=[]
        orderedFID.append(startingFID[0])
        i=1
        while i<len(end_has_start):
            orderedFID.append(end_has_start[orderedFID[i-1]])
            i=i+1
        orderedFID.append(endingFID[0])
    
        #now the COMID is ordered
        if orderedFID[0]==orderedFID[-1]:
            #comid starts and ends in same cell, doesn't help with cell routing, skip
            print "starts and ends in same cell %d" % comid
            continue    
        if(top[orderedFID[0]] < top[orderedFID[-1]]):
            orderedFID.reverse()
        elif(top[orderedFID[0]] == top[orderedFID[-1]]):
            bigtest=0
            testcomid=comid
            seen=dict()
            checkdirflag=0
            endsegment=0
            startsegment=0
            i=0
            while i < len(fromCOMIDlist[testcomid]):
                if fromCOMIDlist[testcomid][i] == 99999:
                    endsegment=1
                i=i+1
            i=0
            if testcomid not in toCOMIDlist:
                startsegment=1
    ##        if (startsegment == 1 and endsegment==1):
    ##            print "check "+str(testcomid)+" start/end"
    ##            exit()
            while bigtest==0:
                theQuery="COMID= "+str(testcomid)
                gp.MakeFeatureLayer(NHD,"temp_nhd_lyr",theQuery)
                nhds=gp.SearchCursor("temp_nhd_lyr")
                nhd=nhds.Next()
                print theQuery
                checkelev=nhd.MINELEVSMO
                del nhds
                flag=0
                if startsegment != 1:
                #there is an upstream comid, keep looking upstream until
                #a non-equal elevation is found and see if it is more or less
                    for prevCOMID in toCOMIDlist[testcomid]:
                        print "comid = "+str(testcomid)+" previous comid= "+str(prevCOMID)
                        theQuery="COMID = "+str(prevCOMID)
                        gp.MakeFeatureLayer(NHD,"temp_nhd_lyr",theQuery)
                        nhds=gp.SearchCursor("temp_nhd_lyr")
                        nhd=nhds.Next()
                        prevelev=nhd.MINELEVSMO
                        if prevelev < checkelev:
                            flag=1
                        if prevelev == checkelev and flag == 0:
                            flag=2
                        if prevelev > checkelev and flag ==2:
                            flag=0  
                elif endsegment != 1:
                    #there isn't an upstream COMID, look downstream and check
                    for nextCOMID in fromCOMIDlist[testcomid]:
                        print "comid = "+str(testcomid)+" next comid= "+str(nextCOMID)
                        if nextCOMID==99999:
                            flag=0
                            continue
                        theQuery="COMID = "+str(nextCOMID)
                        gp.MakeFeatureLayer(NHD,"temp_nhd_lyr",theQuery)
                        nhds=gp.SearchCursor("temp_nhd_lyr")
                        nhd=nhds.Next()
                        nextelev=nhd.MINELEVSMO
                        if nextelev > checkelev:
                            flag=1
                        if nextelev == checkelev and flag == 0:
                            flag=2
                        if nextelev < checkelev and flag ==2:
                            flag=0
                else:
                    #a one-segment long piece, no upstream or downstream
                    flag=0
                if flag==0:
                    bigtest=1
                elif flag==1:
                    orderedFID.reverse()
                    bigtest=1
                elif flag==2:
                    print "flag = 2, COMID="+str(testcomid)
                    if startsegment != 1:
                        k=0
                        while k < len(toCOMIDlist[testcomid]):
                            if not toCOMIDlist[testcomid][k] in seen:
                                seen[toCOMIDlist[testcomid][k]]=1
                                testcomid=toCOMIDlist[testcomid][k]
                                checkdirflag=1
                            k=k+1
                    else:
                        k=0
                        while k < len(fromCOMIDlist[testcomid]):
                            if not fromCOMIDlist[testcomid][k] in seen:
                                seen[fromCOMIDlist[testcomid][k]]=1
                                testcomid=fromCOMIDlist[testcomid][0]
                            k=k+1
                    bigtest=0
    
        
        for i in range(0,len(orderedFID)-1):
            fromcell=cellnumber[orderedFID[i]]
            tocell=cellnumber[orderedFID[i+1]]
            OUT.write(str(fromcell)+","+str(tocell)+"\n")
    
        for i in range(0,len(orderedFID)):
            RCH.write(",".join(map(str,[cellnumber[orderedFID[i]],comid,hydroseqused[comid],uphydroseq[comid],dnhydroseq[comid],SFRseq,i+1]))+"\n")
        percentdone=round(100*icount/len(sortedlist),2)
        print "%s %% done" %(percentdone)
ofpfix.close()
if needfix:
    print 'Manual fixes of some COMIDs required -->'
    print 'see %s for COMIDs and instruction.' %(outfixfilename)
    print 'then rern RouteRiverCells.py'
OUT.close()
RCH.close()
gp.refreshcatalog(path)
del gp

    
