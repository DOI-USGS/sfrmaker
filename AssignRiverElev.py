# AssignRiverElev.py
# Description: Takes shapefile generated from NHDplus with start (X,Y), end(X,Y),
# and line length added and a cell-by-cell intersection of the NHDplus
# and assigns elevations to the cell-by-sell intersection river segments
# based on the distance travelled along the line and the maximum and
# minimum hydro-smoothed elevations from NHDplus.
#
# input shapefiles must have line starts and ends added
#
# edited version for Red Cliff refinement work
#
# Requirements: os, sys, arcpy, defaultdict, math
#
# Author: H.W. Reeves; USGS Michigan Water Science Center
# Date 7/9/12
#
import os
import sys
import arcpy
from collections import defaultdict
import math

# HARD CODE INPUT SHAPEFILES

INTERSECT="river_explode.shp"
OUTTAB="river_elevs.dbf"
outfile=open("fix_comids.txt",'w')

# Global Input file for SFR utilities (see also for input instructions)
infile="SFR_setup.in"

# Get input files locations
infile=open(infile,'r').readlines()
inputs=defaultdict()
for line in infile:
    if "#" not in line.split()[0]:
        varname=line.split("=")[0]
        var=line.split("=")[1].split()[0]
        inputs[varname]=var.replace("'","") # strip extra quotes

ELEV=inputs["ELEV"]

#set elevation difference to assign same value to all cells
#so that check does not get messed up by rounding floats
eps=1.0e-02

#set workspace
path=os.getcwd()
arcpy.env.workspace=path

#create output DBF table
if arcpy.Exists(OUTTAB):
    arcpy.Delete_management(OUTTAB)
arcpy.CreateTable_management(path,OUTTAB)
arcpy.AddField_management(OUTTAB,"OLDFID","LONG")
arcpy.AddField_management(OUTTAB,"CELLNUM","LONG")
arcpy.AddField_management(OUTTAB,"ELEVMAX","DOUBLE")
arcpy.AddField_management(OUTTAB,"ELEVAVE","DOUBLE")
arcpy.AddField_management(OUTTAB,"ELEVMIN","DOUBLE")

#loop over COMIDs explode shapefile and read information into
#hashes
startx=defaultdict(list)
starty=defaultdict(list)
endx=defaultdict(list)
endy=defaultdict(list)
maxsmoothelev=dict()
minsmoothelev=dict()
totallength=dict()
FID=defaultdict(list)
lengthft=defaultdict(defaultdict)
cellnum=defaultdict(defaultdict)

segments=arcpy.SearchCursor(INTERSECT)
for seg in segments:
    comid=int(seg.COMID)
    startx[comid].append(float(seg.X_start))
    starty[comid].append(float(seg.Y_start))
    endx[comid].append(float(seg.X_end))
    endy[comid].append(float(seg.Y_end))
    FID[comid].append(int(seg.FID))
    lengthft[comid][int(seg.FID)]=float(seg.LengthFt)
    cellnum[comid][int(seg.FID)]=int(seg.node)
    lengthkmft=float(seg.LENGTHKM)*1000.0*3.28
    if(lengthkmft<0.01):
        print "check LENGTHKM of comid= ",comid
        exit()
    #the 04 smoothed elevations are in cm
    #convert to ft.
    maxsmoothelev[comid]=float(seg.MAXELEVSMO)*3.2808/100.
    minsmoothelev[comid]=float(seg.MINELEVSMO)*3.2808/100.
del seg
del segments

print "shapefile read in"

OUT2=open('check.txt','w')
#pstring=['comid','cellnum','orderedFID','maxcellrivelev','avecellrivelev','mincellrivelev']
#OUT2.write(",".join(pstring)+"\n")
#have all the information, loop over COMIDs
#if min smooth elevation = max smooth elevation, set the elevation
#the same for all the cells with the same COMID and go to next one
#otherwise use distance along the line to estimate elevation for
#each cell
fact=0.95
segcounter=0
rows=arcpy.InsertCursor(OUTTAB)
for comid in FID.iterkeys():
    segcounter=segcounter+1
    print "COMID = ",comid," segment is ",segcounter
    if (maxsmoothelev[comid]<=minsmoothelev[comid]+eps and
    maxsmoothelev[comid]>=minsmoothelev[comid]-eps):
        for i in range(0,len(FID[comid])):
            #pstring=[cellnum[comid][FID[comid][i]],maxsmoothelev[comid],maxsmoothelev[comid],maxsmoothelev[comid]]
            #OUT2.write(",".join(map(str,pstring))+"\n")
            row=rows.newRow()
            row.OLDFID=FID[comid][i]
            row.CELLNUM=cellnum[comid][FID[comid][i]]
            row.ELEVMAX=maxsmoothelev[comid]
            row.ELEVAVE=maxsmoothelev[comid]
            row.ELEVMIN=maxsmoothelev[comid]
            rows.insertRow(row)
        continue
    #order the entries for each COMID-match starts and ends
    #in developing the shapefile, segments < 1.0 ft were deleted, so
    #check if start(x,y) is within a fact (set above) of end(x,y)
    fidlist=FID[comid]
    start_has_end=dict()
    end_has_start=dict()

    for i in range(0,len(fidlist)):
        for j in range(0,len(fidlist)):
            if j==i:
                continue
            diffx=startx[comid][i]-endx[comid][j]
            diffy=starty[comid][i]-endy[comid][j]
            if comid==1799009:
                pstring=(i, j, startx[comid][i], endx[comid][j],starty[comid][i], endy[comid][j], diffx, diffy)
                OUT2.write(",".join(map(str,pstring))+"\n")
            if(math.fabs(diffx)<fact and math.fabs(diffy)<fact):
                start_has_end[fidlist[i]]=fidlist[j]
                if comid==1799009:
                    pstring=('start has end', fidlist[i], fidlist[j])
                    OUT2.write(",".join(map(str,pstring))+"\n")
                break
        for j in range(0,len(fidlist)):
            if j==i:
                continue
            diffx=endx[comid][i]-startx[comid][j]
            diffy=endy[comid][i]-starty[comid][j]
            if comid==1799009:
                pstring=(i, j, startx[comid][i], endx[comid][j],starty[comid][i], endy[comid][j], diffx, diffy)
                OUT2.write(",".join(map(str,pstring))+"\n")
            if(math.fabs(diffx)<fact and math.fabs(diffy)<fact):
                end_has_start[fidlist[i]]=fidlist[j]
                if comid==1799009:
                    pstring=('end has start', fidlist[i], fidlist[j])
                    OUT2.write(",".join(map(str,pstring))+"\n")
                break
    #find key in start_has_end that didn't match and end and
    #key in end_has_start that didn't match a start
        
    numstart=0
    numend=0
    startingFID=[]
    endingFID=[]
    for test in fidlist:
        if test not in start_has_end:
            startingFID.append(test)
            numstart=numstart+1
        if test not in end_has_start:
            endingFID.append(test)
            numend=numend+1
    if (numstart!=1 or numend !=1):
        outfile.write("numstart ="+ str(numstart)+" \n")
        outfile.write("numend = "+ str(numend)+" \n")
        outfile.write("starting FIDs: " + ",".join(map(str,startingFID))+"\n")
        outfile.write("ending FIDs: " + ",".join(map(str,endingFID))+"\n")
        outfile.write("manually fix COMID = " + str(comid)+" \n")
        continue

    orderedFID=[]
    orderedFID.append(startingFID[0])
    for i in range(1,len(end_has_start)):
        orderedFID.append(end_has_start[orderedFID[i-1]])
    orderedFID.append(endingFID[0])
    #total length read through lengthkm didn't always match up
    #to the sum of the lengths of the segments (exactly), sumup total length
    totallength=0
    for i in range(0,len(orderedFID)):
        totallength=totallength+lengthft[comid][orderedFID[i]]
    if totallength==0:
        exit('check length ft for FIDs in COMID= %d' % comid)
    slope=(maxsmoothelev[comid]-minsmoothelev[comid])/totallength
    distance=0.
    for i in range(0,len(orderedFID)):
        maxcellrivelev=maxsmoothelev[comid]-slope*distance
        distance=distance+lengthft[comid][orderedFID[i]]
        mincellrivelev=maxsmoothelev[comid]-slope*distance
        avecellrivelev=0.5*(maxcellrivelev+mincellrivelev)
        #pstring=[comid,cellnum[comid][orderedFID[i]],orderedFID[i],maxcellrivelev,avecellrivelev,mincellrivelev]
        #OUT2.write(",".join(map(str,pstring))+"\n")
        row=rows.newRow()
        row.OLDFID=orderedFID[i]
        row.CELLNUM=cellnum[comid][orderedFID[i]]
        row.ELEVMAX=maxcellrivelev
        row.ELEVAVE=avecellrivelev
        row.ELEVMIN=mincellrivelev
        rows.insertRow(row)
        

del row
del rows
OUT2.close()
outfile.close()
arcpy.RefreshCatalog(path)

print "\nJoining %s to %s; saving as %s (might take a minute)" %(OUTTAB,INTERSECT,ELEV)
arcpy.CopyFeatures_management(INTERSECT,ELEV)
arcpy.JoinField_management(ELEV,"FID",OUTTAB,"OLDFID")

print "Done! Check fix_comids.txt for segments with multiple starts/ends that must be manually deleted."