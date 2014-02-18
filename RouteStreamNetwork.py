# RouteStreamNetwork.py
# Description: Takes river_elev shapefile from AssignRiverElev.py
# and uses PlusFlow.dbf table (fromCOMID -> toCOMID) information
# to check routing of the whole network
#
# Requirements: os, sys, re, arcgisscripting, defaultdict
#
# Author: H.W. Reeves; USGS Michigan Water Science Center
# Date 1/29/2010
#
import os
import sys
import re
import arcgisscripting
from collections import defaultdict

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

ELEV=inputs["ELEV"] # 'river_w_elevations.shp'
FLOW=inputs["FLOW"] # 'PlusFlow_merge.dbf'

INTERSECT="river_explode.shp"
RIVER_ELEVS="river_elevs.dbf"
'''
print "\nJoining %s to %s; saving as %s (might take a minute)" %(RIVER_ELEVS,INTERSECT,ELEV)
arcpy.CopyFeatures_management(INTERSECT,ELEV)
arcpy.JoinField_management(ELEV,"FID",RIVER_ELEVS,"OLDFID")
'''
path=os.getcwd()
gp=arcgisscripting.create(9.3)
gp.Workspace=path

# sets over-write output option to true for output shapefile (river_elev.shp)
gp.overwriteoutput = 1

#delete any working layers if found
if gp.Exists("temp_lyr"):
    gp.Delete_management("temp_lyr")
if gp.Exists("flow_lyr"):
    gp.Delete_management("flow_lyr")
if gp.Exists("temp_table"):
    gp.Delete_management("temp_table")

#read through river_elev shapefile and get list of all COMIDs; then
#read through NHDFlow and pull the from/to information for those COMIDs
comids=dict()
segments=gp.SearchCursor(ELEV)
for seg in segments:
    #print str(seg.COMID)
    comids[seg.COMID]=1

del seg
del segments

print "have all the comids"
#if check_network.txt (from and to comid list)
#does not exist, create on by setting maketable=1
lcount=0
maketable=1
lencomids = len(comids)
print "Creating to/from routing list..."
if maketable==1:
    outfile=open("check_network.txt",'w')
    outfile_FAILS = open("check_network_COM_FAILS.txt",'w')
    for com in comids.iterkeys():
        lcount+=1
        percentdone=100.0*lcount/lencomids
        print "%4.2f percent done" %(percentdone)
        theQuery="FROMCOMID= %d" %(com)
        print theQuery
        gp.MakeTableView(FLOW,"temp_table",theQuery)
        result=gp.GetCount("temp_table")
        number=int(result.GetOutput(0))
        print number
        if number> 0:
            selected=gp.SearchCursor("temp_table")
            sel=selected.Next()
            while sel:
                if sel.TOCOMID in comids:
                    outfile.write("%d,%d\n" %(com,sel.TOCOMID))
                else:
                    outfile.write("%d,99999\n" %(com))
                sel=selected.Next()
        else:
            theQuery="TOCOMID= "+str(com)
            print theQuery
            gp.MakeTableView(FLOW,"temp_table",theQuery)
            result=gp.GetCount("temp_table")
            number=int(result.GetOutput(0))
            print number
            if number > 0:
                selected=gp.SearchCursor("temp_table")
                sel=selected.Next()
                while sel:
                    if sel.FROMCOMID in comids:
                        outfile.write("%d,%d\n" %(sel.FROMCOMID,com))
                    else:
                        outfile.write("99999,%d\n" %(com))
                    sel=selected.Next()
            else:
                outfile_FAILS.write("%d,check routing, COMID in neither FROM or TO\n" %(com))
    try:
        outfile.close()
    except:
        print "problem closing outfile\n"
        exit()
    try:
        outfile_FAILS.close()
    except:
        print "problem closing outfile_FAILS\n"
        exit()

fromCOMIDlist=defaultdict(list)

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

infile=open("check_network.txt",'r').readlines()

for line in infile:
    vals=line.strip().split(",")
    fromin=int(float(vals[0]))
    toin=int(float(vals[1]))
    if not toin in clipto:
        if not fromin in clipfrom:
            fromCOMIDlist[fromin].append(toin)
        else:
            fromCOMIDlist[fromin].append(99999)
#now have full from->to table as a dictionary of lists
#go through river_elev file COMID by COMID and check
#that next downstream COMID has lower river elevation

outfile=open("Flagged_comids.txt",'w')
outfile.write("COMID,DWNCOMID,MINELEV,DWNELEV\n")
gp.MakeFeatureLayer(ELEV,"flow_lyr")
total=len(comids)
icount=1
print "Checking for downstream routing..."
for com in comids.iterkeys():
    print icount,"out of ",total
    icount=icount+1
    flag=0
    theQuery="COMID= "+ str(com)
    #print theQuery
    if com==1798409:
        print 'here I am'
    gp.SelectLayerByAttribute("flow_lyr","NEW_SELECTION",theQuery)
    gp.MakeFeatureLayer("flow_lyr","temp_lyr")
    selected=gp.SearchCursor("temp_lyr")
    sel=selected.Next()
    maxelev=0.0
    minelev=999999.
    while sel:
        if sel.ELEVAVE > maxelev:
            maxelev=sel.ELEVAVE
        if sel.ELEVAVE < minelev:
            minelev=sel.ELEVAVE
        sel=selected.Next()
    del selected
    i=0
    while i<len(fromCOMIDlist[com]):
        dwnstrm=fromCOMIDlist[com][i]
        theQuery="COMID= %d" %(dwnstrm)
        gp.SelectLayerByAttribute("flow_lyr","NEW_SELECTION",theQuery)
        gp.MakeFeatureLayer("flow_lyr","temp_lyr")
        selected=gp.SearchCursor("temp_lyr")
        sel=selected.Next()
        dwdnmax=0.0
        dwdnmin=999999.
        while sel:
            if sel.ELEVAVE > dwdnmax:
                dwdnmax=sel.ELEVAVE
            if sel.ELEVAVE < dwdnmin:
                dwdnmin=sel.ELEVAVE
            sel=selected.Next()
        del selected
        if(dwdnmax > (minelev+0.01)):
            outfile.write("%d,%d,%f,%f\n" %(com,dwnstrm,minelev,dwdnmax))
            flag=1
        i=i+1
    if(flag == 0):
        sys.stdout.write(str(com) + " is OK\n")
        
outfile.close()
gp.refreshcatalog(path)
del gp

