# Driver script for steps 1-9 in Howard Reeves' (USGS MI Water Science Center) SFR Workflow. Includes code from Howard's explode_add_geometry.py
# Authored by Andrew Leaf, USGS WI Water Water Science Center


import math
import numpy as np
from collections import defaultdict
import os
import re
import arcpy
from arcpy.sa import *
import datetime, time
import SFR_arcpy

ts = time.time()
st_time = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

# Global Input file for SFR utilities (see also for input instructions)
infile="SFR_setup.in"

# set up an output file for logging
ofp = open('SFR_preproc.log','w')
ofp.write('SFR_preproc log.')
ofp.write('\n' + '#' * 25 + '\nStart Time: %s\n' %st_time +'#' * 25 + '\n')
# Get input files locations
infile=open(infile,'r').readlines()
inputs=defaultdict()
for line in infile:
    if "#" not in line.split()[0]:
        varname=line.split("=")[0]
        var=line.split("=")[1].split()[0]
        inputs[varname]=var.replace("'","") # strip extra quotes

reach_cutoff=float(inputs["reach_cutoff"])
MFgrid=inputs["MFgrid"]
MFdomain=inputs["MFdomain"]
DEM=inputs["DEM"]
PlusflowVAA=inputs["PlusflowVAA"]
Elevslope=inputs["Elevslope"]
Flowlines=inputs["Flowlines"]
Flowlines_unclipped=inputs["Flowlines_unclipped"]
arcpy_path=inputs["arcpy_path"]

# environmental settings
path=os.getcwd()
arcpy.env.workspace = path
arcpy.env.overwriteOutput = True 
arcpy.env.qualifiedFieldNames = False

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# get name of field (useful for case issues and appended field names in joined tables, etc)
def getfield(table,joinname):
    Fields=arcpy.ListFields(table)
    joinname=joinname.lower()
    for field in Fields:
        if joinname in field.name.lower():
            joinname=field.name
            break
    return(joinname)


print "Clip original NHD flowlines to model domain..."
# this creates a new file so original dataset untouched
arcpy.Clip_analysis(Flowlines_unclipped,MFdomain,Flowlines)

print "joining PlusflowVAA and Elevslope tables to NHD Flowlines..."
# copy original tables and import copies
print "importing tables..."
arcpy.MakeTableView_management(Elevslope,"Elevslope")
arcpy.MakeTableView_management(PlusflowVAA,"PlusflowVAA")

# delete all unneeded fields
fields2keep=["comid","divergence","lengthkm","thinnercod","maxelevsmo","minelevsmo","hydroseq","uphydroseq","dnhydroseq","reachcode","streamorde","arbolatesu","fcode"]
fields2keep=[x.lower() for x in fields2keep]
ofp.write('Joining %s with %s: fields kept:\n' %(Elevslope,Flowlines))
ofp.write('%s\n' %(','.join(fields2keep)))
print "\nkeeping: %s fields; deleting the rest" %(','.join(fields2keep))
Join=True # whether or not PlusflowVAA and Elevslope need to be joined
for table in ["PlusflowVAA","Elevslope"]:
    fields2delete=[]
    Fields=arcpy.ListFields(table)

    for field in Fields:
        name=field.name
        namel=name.lower()
        if namel in fields2keep:
            if table=="PlusflowVAA" and namel=="maxelevsmo":
                Join=False
            continue
        elif namel=="oid":
            continue
        else:
            fields2delete.append(field.name)
            print field.name
    if len(fields2delete)>0:
        arcpy.DeleteField_management(table,fields2delete)

# Join PlusflowVAA and Elevslope to Flowlines
arcpy.MakeFeatureLayer_management(Flowlines,"Flowlines")

# If not already, permanently join PlusflowVAA and Elevslope, then to Flowlines    
if Join:
    print "\nJoining Elevslope to PlusflowVAA...\n"
    comid1=getfield("PlusflowVAA","comid")
    comid2=getfield("Elevslope","comid")
    arcpy.JoinField_management("PlusflowVAA",comid1,"Elevslope",comid2)
else:
    print "PlusflowVAA and Elevslope already joined from previous run..."

print "Joining PlusflowVAA to NHDFlowlines...\n"    
comid1=getfield("Flowlines","comid")

# join to Flowlines, keeping only common
SFR_arcpy.general_join(Flowlines, "Flowlines",comid1,PlusflowVAA,"comid",True)

# reopen flowlines as "Flowlines" --> clunky a bit to save and reopen, but must do so
arcpy.MakeFeatureLayer_management(Flowlines,"Flowlines")

print "\n"
ofp.write('\n' + 25*'#' + '\nRemoving segments with no elevation information, and with ThinnerCod = -9..\n')
print "Removing segments with no elevation information, and with ThinnerCod = -9..."
ThinnerCod=getfield("Flowlines","thinnercod")
MaxEl=getfield("Flowlines","maxelevsmo")
comid=getfield("Flowlines","comid")
FLtable=arcpy.UpdateCursor("Flowlines")
zerocount=0
tcount=0
for segments in FLtable:
    if segments.getValue(MaxEl)==0:
        print "%d no elevation data" %segments.getValue(comid)
        ofp.write("%d no elevation data\n" %segments.getValue(comid))
        FLtable.deleteRow(segments)
        zerocount+=1
    elif segments.getValue(ThinnerCod)==-9:
        print "%d ThinnerCod=-9" %segments.getValue(comid)
        ofp.write("%d ThinnerCod=-9\n" %segments.getValue(comid))
        FLtable.deleteRow(segments)
        tcount+=1
            
print "...removed %s segments with no elevation data" %(zerocount)
ofp.write("...removed %s segments with no elevation data\n" %(zerocount))
print "...removed %s segments with ThinnerCod = -9\n" %(tcount)
ofp.write("...removed %s segments with ThinnerCod = -9\n" %(tcount))
print "Performing spatial join (one-to-many) of NHD flowlines to model grid to get river cells...(this step may take several minutes or more)\n"
arcpy.SpatialJoin_analysis(MFgrid,"Flowlines","river_cells.shp","JOIN_ONE_TO_MANY","KEEP_COMMON")

# add in cellnum field for backwards compatibility
arcpy.AddField_management("river_cells.shp","CELLNUM","LONG")
arcpy.CalculateField_management("river_cells.shp","CELLNUM","!node!","PYTHON")

print "Dissolving river cells on cell number to isolate unique cells...\n"
node=getfield("river_cells.shp","node")
arcpy.Dissolve_management("river_cells.shp","river_cells_dissolve.shp",node)

print "Exploding NHD segments to grid cells using Intersect and Multipart to Singlepart..."
arcpy.Intersect_analysis(["river_cells_dissolve.shp","Flowlines"],"river_intersect.shp")
arcpy.MultipartToSinglepart_management("river_intersect.shp","river_explode.shp")
print "\n"
print "Adding in stream geometry"
#set up list and dictionary for fields, types, and associated commands
fields=('X_start','Y_start','X_end','Y_end','LengthFt')
types={'X_start':'DOUBLE',
       'Y_start':'DOUBLE',
       'X_end':'DOUBLE',
       'Y_end':'DOUBLE',
       'LengthFt':'DOUBLE'}
commands={'X_start':"float(!SHAPE.firstpoint!.split()[0])",
       'Y_start':"float(!SHAPE.firstpoint!.split()[1])",
       'X_end':"float(!SHAPE.lastpoint!.split()[0])",
       'Y_end':"float(!SHAPE.lastpoint!.split()[1])",
       'LengthFt':"float(!SHAPE.length!)"}
    
#add fields for start, end, and length
for fld in fields:
    arcpy.AddField_management("river_explode.shp",fld,types[fld])

#calculate the fields
for fld in fields:
    print "\tcalculating %s(s)..." %(fld)
    arcpy.CalculateField_management("river_explode.shp",fld,commands[fld],"PYTHON")
ofp.write('\n' + 25*'#' + '\nRemoving reaches with lengths less than or equal to %s...\n' %(reach_cutoff))
print "\nRemoving reaches with lengths less than or equal to %s..." %(reach_cutoff)
comid=getfield("river_explode.shp","comid")
node=getfield("river_explode.shp","node")
Length=getfield("river_explode.shp","lengthft")
table=arcpy.UpdateCursor("river_explode.shp")
count=0
for reaches in table:
    if reaches.getValue(Length)<=reach_cutoff:
        print "segment: %d, cell: %s, length: %s" %(reaches.getValue(comid),reaches.getValue(node),reaches.getValue(Length))
        ofp.write("segment: %d, cell: %s, length: %s\n" 
                  %(reaches.getValue(comid),reaches.getValue(node),reaches.getValue(Length)))
        table.deleteRow(reaches)
        count+=1
print "removed %s reaches with lengths <= %s\n" %(count,reach_cutoff)
ofp.write("removed %s reaches with lengths <= %s\n" %(count,reach_cutoff))
print "removing cells corresponding to those reaches..."

# temporarily join river_cells_dissolve to river explode; record nodes with no elevation information
arcpy.MakeFeatureLayer_management("river_cells_dissolve.shp","river_cells_dissolve")
arcpy.MakeFeatureLayer_management("river_explode.shp","river_explode")
arcpy.AddJoin_management("river_cells_dissolve","node","river_explode","node","KEEP_ALL") # this might not work as-in in stand-alone mode
node=getfield("river_cells_dissolve","node")
maxelevsmo=getfield("river_cells_dissolve","maxelevsmo")
table=arcpy.SearchCursor("river_cells_dissolve")
nodes2delete=[]
for row in table:
    if row.isNull(maxelevsmo):
        nodes2delete.append(row.getValue(node))
arcpy.RemoveJoin_management("river_cells_dissolve","river_explode")

# remove nodes with no elevation information from river_explode
ofp.write('\n' + 25*'#' + '\nRemoving nodes with no elevation information from river_explode\n')
print 'Removing nodes with no elevation information from river_explode'
node=getfield("river_cells_dissolve.shp","node")
table=arcpy.UpdateCursor("river_cells_dissolve.shp")
count=0
for cells in table:
    if cells.getValue(node) in nodes2delete:
        print "%d" %(cells.getValue(node))
        ofp.write('%d\n' %(cells.getValue(node)))
        table.deleteRow(cells)
        count+=1
print "removed %s cells\n" %(count)
ofp.write("removed %s cells" %(count))

print "removing any remaining disconnected reaches..."
ofp.write('\n' + 25*'#' + '\nremoving any remaining disconnected reaches...\n')
node=getfield("river_explode.shp","node")
comid=getfield("river_explode.shp","comid")
table=arcpy.UpdateCursor("river_explode.shp")
count=0
for reaches in table:
    if reaches.getValue(node) in nodes2delete:
        print "%d" %(reaches.getValue(node))
        ofp.write('%d\n' %(reaches.getValue(node)))
        table.deleteRow(reaches)
        count+=1
if count>0:
    print "removed %s disconnected reaches" %(count)
    ofp.write("removed %s disconnected reaches\n" %(count))
else:
    print "no disconnected reaches found!"
    ofp.write("no disconnected reaches found!\n")
print "\n"
print "Done with pre-processing, ready to run intersect.py!"
ofp.write('\n' + '#' * 25 + '\nDone with pre-processing, ready to run intersect.py!')
ts = time.time()
end_time = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

ofp.write('\n' + '#' * 25 + '\nEnd Time: %s\n' %end_time +'#' * 25)
ofp.close()