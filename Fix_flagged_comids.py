# Fix flagged COMIDs from RouteStreamNetwork

import arcpy
import os
from collections import defaultdict
import pdb

# Global Input file for SFR utilities
infile="SFR_setup.in"
outfile="fix_flagged.out"

# Settings
increment=1 # amount to raise headwater elevation if gradient=0

path=os.getcwd()
arcpy.env.workspace = path
arcpy.env.overwriteOutput = True

if arcpy.Exists("temp_lyr"):
    arcpy.Delete_management("temp_lyr")


# Get input files locations
infile=open(infile,'r').readlines()
inputs=defaultdict()
for line in infile:
    if "#" not in line.split()[0]:
        varname=line.split("=")[0]
        var=line.split("=")[1].split()[0]
        inputs[varname]=var.replace("'","") # strip extra quotes

ELEV=inputs["ELEV"]
PlusflowVAA=inputs["PlusflowVAA"]

# bring in to/from routing list and list of upstream routing
check_network=open('check_network.txt','r').readlines()
flaggedcomids=open('Flagged_comids.txt','r').readlines()
names=flaggedcomids[0].split()

checknw=defaultdict()
for lines in check_network:
    line=lines.strip().split(',')
    checknw[line[0]]=line[1]

flagged=defaultdict()

for lines in flaggedcomids[1:]:
    line=lines.strip().split(',')
    line=[l.strip(' ') for l in line]
    flagged[line[0]]=line[1:]
    
# bring in river_w_elevations table
arcpy.MakeFeatureLayer_management(ELEV,"river_w_elevations")

# Make list of COMIDs to delete
# Delete COMIDs that are artificial paths and with with Divergence=2)
print "Removing segments that are artificial paths and with with Divergence=2 (these are secondary connections for braids in channel)..."

# just in case, wipe out any existing selections
arcpy.SelectLayerByAttribute_management("river_w_elevations","CLEAR_SELECTION")

for comids in flagged.iterkeys():
    theQuery="COMID= "+ str(comids)
    arcpy.SelectLayerByAttribute_management("river_w_elevations","ADD_TO_SELECTION",theQuery)
    
arcpy.MakeFeatureLayer_management("river_w_elevations","temp_lyr")
selected=arcpy.UpdateCursor("temp_lyr")
count=0

deletedcomids=[]

for rows in selected:
    ArtificialPath=False
    if rows.getValue("FTYPE")=="ArtificialPath":
        ArtificialPath=True
        if rows.getValue("Divergence")==2:
            comid=rows.getValue("COMID")
            count+=1
            if comid not in deletedcomids:
                deletedcomids.append(comid)
                print "removing " + str(comid)
                # delete references to comid in check_network.txt
                for line in check_network:
                    if str(comid) in line:
                        print "removing lines from check_network.txt: " + line.strip()
                        check_network.remove(line)
            selected.deleteRow(rows)
            
print "deleted %s comids (%s reaches)" %(len(deletedcomids),count) 
del selected

# Make list of COMIDs to fix routing
comidsremaining=[c for c in flagged.iterkeys() if int(c) not in deletedcomids]
ELEVtable=arcpy.SearchCursor(ELEV)
reaches=defaultdict(list)
ofp=open(outfile,'w')
ofp.write('comid,maxelev,minelev,endnode_maxelev,endnode_minelev,Upcomid,min_elev_DOWNcomid,DWNELEV_from_Flagged,min-max,initial_reach_elevations...\n')
print "fixing segment elevations..."
for comid in comidsremaining[1:]:
    print "fixing comid " + str(comid)
    theQuery="COMID= "+ str(comid)
    arcpy.SelectLayerByAttribute_management("river_w_elevations","NEW_SELECTION",theQuery)
    arcpy.MakeFeatureLayer_management("river_w_elevations","temp_lyr")
    
    # In case comidsremaining contains comids that have already been deleted
    # table will have no rows, skip to next
    # check that these have been removed from check_network
    rows=int(arcpy.GetCount_management("temp_lyr").getOutput(0))
    if rows==0:
        print "COMID %s already deleted" %(comid)
        for line in check_network:
            if str(comid) in line:
                print "removing lines from check_network.txt: " + line.strip()
                check_network.remove(line)
        continue
    
    selected=arcpy.SearchCursor("temp_lyr")
    # assuming that segment end reaches can be identified by min/max reach elevations
    # identify nodes (model cells) correspoding to segment end reaches
    maxelev=0.0
    minelev=999999.
    endnodes=[]
    for row in selected:
        if row.getValue("ELEVMAX") > maxelev:
            maxelev=row.getValue("ELEVMAX")
            endnode_maxelev=row.getValue("node")
        if row.getValue("ELEVMIN")< minelev:
            minelev=row.getValue("ELEVMIN")
            endnode_minelev=row.getValue("node")
    endnodes.append(endnode_maxelev)
    endnodes.append(endnode_minelev)
    del selected
    # end nodes are identified
    # identify for other segments in same end nodes (i.e. upstream or downstream segments)
    #  find which node has segment with lowest ELEVMAX (this should be downstream), and highest ELEVMIN (should be upstream)
    # if no other segments, than this is a headwater (upstream end), flows into/out of the grid, or the segment terminates at a grid cell boundary
    maxelev=0.0
    minelev=999999.
    headwater=False
    for node in endnodes:
        theQuery="node= "+ str(node)
        arcpy.SelectLayerByAttribute_management("river_w_elevations","NEW_SELECTION",theQuery)
        arcpy.MakeFeatureLayer_management("river_w_elevations","temp_lyr")
        selected=arcpy.SearchCursor("temp_lyr")
        rows=int(arcpy.GetCount_management("temp_lyr").getOutput(0))
        # if there aren't any other segments sharing the end node, check routing list to see whether it flows out of the grid
        if rows==1:
            '''TOcomid=checknw[comid]
            if int(TOcomid)==99999:
                headwater=True''' # This is wrong, 99999 in this position means the segment is not routed to any downstream segments
        for row in selected:
            if row.getValue("ELEVMIN") > maxelev:
                if headwater:
                    UPcomid=99999
                    maxelev=row.getValue("ELEVMAX")
                elif rows==1:
                    maxelev=row.getValue("ELEVMAX")
                else:
                    maxelev=row.getValue("ELEVMIN")
                    UPcomid=row.getValue("COMID") # this will be the upstream COMID
            if row.getValue("ELEVMAX") < minelev:
                minelev=row.getValue("ELEVMAX")
                DOWNcomid=row.getValue("COMID") # this will be the downstream COMID; only need this if we're going to change the routing.
    
    DOWNELEV_original=float(flagged[comid][2]) # original downstream elevation based on NHD routing
    
    # if headwater and min/max elevs are the same, bump up maxelev by increment
    # could extend this to fixing zero gradients for non-headwaters, but would need additional info on up and downstream segments
    dh=DOWNELEV_original-maxelev
    if headwater:
        while dh >=0:
            print "fixing headwater elevation" %(comid)
            maxelev+=increment
            dh=DOWNELEV_original-maxelev
            print "maxelev=%s, gradient=%s" %(maxelev,dh)
    del selected
    
    
    ofp.write(','.join(map(str,[comid,maxelev,minelev,endnode_maxelev,endnode_minelev,UPcomid,DOWNcomid,DOWNELEV_original,dh]))+',')  
      
    # have segment max/min elevations and corresponding nodes, now interpolate reaches in between
    # query river_w_elevations for all reaches with comid again
    theQuery="COMID= "+ str(comid)
    arcpy.SelectLayerByAttribute_management("river_w_elevations","NEW_SELECTION",theQuery)
    arcpy.MakeFeatureLayer_management("river_w_elevations","temp_lyr")
    selected=arcpy.SearchCursor("temp_lyr","","","","ELEVAVE A") # sort ascending by ELEVAVE field
    
    # reset segment end elevations, and fill in with linear interpolation
    # first get total segment length and ave slope
    dist=0
    for rows in selected:
        clength=rows.getValue("LengthFt")
        dist+=clength
    slope=dh/dist
    print "dh = %s" %(dh)
    print "slope = %s" %(slope)
    del selected
    selected=arcpy.UpdateCursor("temp_lyr","","","","ELEVAVE D")
    
    dist=0
    print "updating river_w_elevations..."
    for rows in selected:
        if dist==0:
            rows.ELEVMAX=maxelev
        else:
            rows.ELEVMAX=prev_min    
        length=rows.LengthFt
        dist+=length
        rows.ELEVAVE=rows.ELEVMAX+slope*0.5*length
        ofp.write(str(rows.ELEVAVE)+',')
        rows.ELEVMIN=rows.ELEVMAX+slope*length
        prev_min=rows.ELEVMIN
        print "%s,%s,%s" %(rows.ELEVMAX,rows.ELEVAVE,rows.ELEVMIN)
        selected.updateRow(rows)
    ofp.write('\n')
    del selected
ofp.close()

# write out new check network file
print "updating check_network.txt; previous version renamed to 'check_network.old'"
os.rename('check_network.txt','check_network.old')

ofp=open('check_network.txt','w')
for line in check_network:
    ofp.write(line)
ofp.close()
print "Done!"
    


