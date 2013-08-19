# intersect.py
# authoed by Howard Reeves, MI Water Science Center
#
import os
import re
import arcpy
import math
from collections import defaultdict

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
        
#hard code filenames and grid information
NHD=inputs["Flowlines_unclipped"] # Unclipped NHDFlowlines, in same coordinates as everything else
GRID=inputs["MFdomain"] # polygon of grid boundary, in same coords
MODLNHD=inputs["Flowlines"] # from step 2)
NEWMOD=inputs["NHD"]
outfile="boundaryClipsRouting.txt"

#set workspace
path=os.getcwd()
arcpy.env.workspace=path
arcpy.env.overwriteOutput=True

#convert to layers
arcpy.MakeFeatureLayer_management(NHD, 'nhd_lyr')
arcpy.MakeFeatureLayer_management(GRID, 'grid_lyr')

print "selecting streams that cross model grid boundary..."
arcpy.SelectLayerByLocation_management('nhd_lyr', 'CROSSED_BY_THE_OUTLINE_OF','grid_lyr')
arcpy.CopyFeatures_management('nhd_lyr','intersect.shp')

#copy the model NHD streams to a temp shapefile, find the ends that match
#the streams that were clipped by the boundary and update the lengthKm,
#minsmoothelev, maxsmoothelev in the MODLNHD
arcpy.CopyFeatures_management(MODLNHD,'temp.shp')

#add fields for start, end, and length to the temp and intersect
#shapefiles  (use LENKM as field name because temp already has LENGHTKM)
print "adding fields for start, end and length..."
shapelist=('temp.shp','intersect.shp')
for shp in shapelist:
    arcpy.AddField_management(shp,'STARTX','DOUBLE')
    arcpy.AddField_management(shp,'STARTY','DOUBLE')
    arcpy.AddField_management(shp,'ENDX','DOUBLE')
    arcpy.AddField_management(shp,'ENDY','DOUBLE')
    arcpy.AddField_management(shp,'LENKM','DOUBLE')

print "calculating new info for fields..."    
#calculate the fields, convert length to km - projection is in feet
for shp in shapelist:
    arcpy.CalculateField_management(shp,'STARTX',"float(!SHAPE.firstpoint!.split()[0])","PYTHON")
    arcpy.CalculateField_management(shp,'STARTY',"float(!SHAPE.firstpoint!.split()[1])","PYTHON")
    arcpy.CalculateField_management(shp,'ENDX',"float(!SHAPE.lastpoint!.split()[0])","PYTHON")
    arcpy.CalculateField_management(shp,'ENDY',"float(!SHAPE.lastpoint!.split()[1])","PYTHON")
    arcpy.CalculateField_management(shp,'LENKM',"float(!SHAPE.length@kilometers!)","PYTHON")
    
#go through intersect, identify which end matches COMID in temp
#find new length in temp and use linear interpolation to get new elev
#finally put COMID out to intersect_comid.txt file to indicate
#if the cut-off is flowing out of grid or into grid so the
#routing tables in the final datasets do not route out of grid
#and back in. Also identifies to user ends of streams that
#could accept inflow conditions for SFR
print "fixing routing for streams that cross the grid boundary..."

#set up some dictionaries to hold changes
newstartx=dict()
newstarty=dict()
newendx=dict()
newendy=dict()
newlength=dict()
inout=dict()
newmaxel=dict()
newminel=dict()

intersects=arcpy.SearchCursor('intersect.shp')
manual_intervention = 0 
for stream in intersects:
    comid=stream.COMID
    query="COMID ="+str(comid)
    stx=stream.STARTX
    sty=stream.STARTY
    endx=stream.ENDX
    endy=stream.ENDY
    lenkm=stream.LENKM
    clippedstream=arcpy.SearchCursor('temp.shp',query)
    for clip in clippedstream:
        clcomid=clip.COMID
        clstx=clip.STARTX
        clsty=clip.STARTY
        clendx=clip.ENDX
        clendy=clip.ENDY
        cllen=clip.LENKM
        clmaxel=clip.MAXELEVSMO
        clminel=clip.MINELEVSMO
        #find the end that matches
        # ATL: indented this paragraph and below if/elif/else statements from original version
        eps=0.01
        stdiffx=stx-clstx
        stdiffy=sty-clsty
        enddiffx=endx-clendx
        enddiffy=endy-clendy
        if math.fabs(stdiffx)<eps and math.fabs(stdiffy)<eps:
            #beginning of stream is kept, it flows out, and
            #maximum elevation is the same
            inout[comid]='OUT'
            newstartx[comid]=clstx
            newstarty[comid]=clsty
            newendx[comid]=clendx
            newendy[comid]=clendy
            newlength[comid]=cllen
            newmaxel[comid]=round(clmaxel)
            slope=clmaxel-clminel
            newminel[comid]=round(clmaxel-slope*cllen/lenkm)
        elif math.fabs(enddiffx)<eps and math.fabs(enddiffy)<eps:
            #end of stream is kept, it flows in, and
            #minimum elevation is the same
            inout[comid]='IN'
            newstartx[comid]=clstx
            newstarty[comid]=clsty
            newendx[comid]=clendx
            newendy[comid]=clendy
            newlength[comid]=cllen
            slope=clmaxel-clminel
            clippedlength=lenkm-cllen
            newmaxel[comid]=round(clmaxel-slope*clippedlength/lenkm)
            newminel[comid]=round(clminel)
        else:
            manual_intervention +=1
            if manual_intervention == 1:
                ofp = open('boundary_manual_fix_issues.txt','w')
                ofp.write('The following COMIDs identify streams that need manual attention.\n')
                ofp.write('Fix in the Flowlines variable from the input file. Then rerun intersect.py\n')
                ofp.write('#' * 16 + '\n')
            print 'both ends are cut off for comid ', comid
            ofp.write('both ends are cut off for comid %d\n' %(comid))
            print 'need to fix it manually'
del intersects, stream, clip, clippedstream
if manual_intervention:
    ofp.write('#' * 16 + '\n')    
    ofp.close()

#create a new working NHD shapefile incorporating new values just found

print "Creating new shapefile " + NEWMOD
arcpy.CopyFeatures_management(MODLNHD,NEWMOD)
intersects=arcpy.UpdateCursor(NEWMOD)
for stream in intersects:
    comid=stream.COMID
    if comid in newstartx:
        stream.LENGTHKM=newlength[comid]
        stream.MAXELEVSMO=newmaxel[comid]
        stream.MINELEVSMO=newminel[comid]
        intersects.updateRow(stream)

del intersects, stream

#put out the routing information
print "Saving routing information to " + outfile
OUT=open(outfile,'w')
OUT.write("FROMCOMID,TOCOMID\n")
for comid in inout:
    if re.match('IN',inout[comid]):
        OUT.write(",".join(map(str,(99999,comid)))+'\n')
    else:
        OUT.write(",".join(map(str,(comid,99999)))+'\n')

OUT.close()
        
    
if manual_intervention:
    print 'Some manual intervention required:\nSee boundary_manual_fix_issues.txt for details'
    
