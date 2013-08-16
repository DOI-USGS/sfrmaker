import arcpy
from collections import defaultdict
import os, shutil

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

Flowlines=inputs["Flowlines"]
Flowlines_unclipped=inputs["Flowlines_unclipped"]

#set workspace
path=os.getcwd()
arcpy.env.workspace=path
arcpy.env.overwriteOutput = True     # allow files to overwrite

# open up the necessary layers
#convert to layers
arcpy.MakeFeatureLayer_management( 'river_cells.shp','rivcells')
arcpy.MakeFeatureLayer_management('river_cells_dissolve.shp','rivcellsdiss')
arcpy.MakeFeatureLayer_management('river_explode.shp','rivexplode')

# clip the river cells to the remaining streams and re-save
print 'interesecting river_cells with river_explode...\n\n'
rivtrim = arcpy.SelectLayerByLocation_management('rivcells',"INTERSECT",'rivexplode')
print '...Saving down results...\n\n'
if arcpy.Exists(os.path.join(path,'river_cells.shp')):
    print 'first removing old version of river_cells.shp'
    arcpy.Delete_management(os.path.join(path,'river_cells.shp'))
arcpy.CopyFeatures_management(rivtrim,os.path.join(path,'river_cells'))
# do the same for the dissolved river cells layer
print 'interesecting river_cells_dissolve with river_explode...\n\n'
rivdisstrim = arcpy.SelectLayerByLocation_management('rivcellsdiss',"INTERSECT",'rivexplode')
print '...Saving down results...\n\n'
if arcpy.Exists(os.path.join(path,'river_cells_dissolve.shp')):
    arcpy.Delete_management(os.path.join(path,'river_cells_dissolve.shp'))
    print 'first removing old version of river_cells_dissolve.shp'
arcpy.CopyFeatures_management(rivdisstrim,os.path.join(path,'river_cells_dissolve'))

# make a backup copy of fix_comids.txt
shutil.copy('fix_comids.txt','fix_comids_backup.txt')