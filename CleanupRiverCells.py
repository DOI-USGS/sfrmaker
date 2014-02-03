# Utility to merge river cells with river_explode

import arcpy
from collections import defaultdict
import os, shutil

# Global Input file for SFR utilities (see also for input instructions)
import SFR_classes as SFRc

infile = 'SFR_input.xml'

SFRdata = SFRc.SFRInput(infile)

#set workspace
path=os.getcwd()
arcpy.env.workspace=path

# open up the necessary layers
#convert to layers
# NB --> for now, making a scratch copy to open to avoid arcpy from not being able to save to an open file
if arcpy.Exists(os.path.join(path, 'river_cells_TMP.shp')):
    arcpy.Delete_management(os.path.join(path, 'river_cells_TMP.shp'))
arcpy.Copy_management(os.path.join(path,SFRdata.CELLS),os.path.join(path,'river_cells_TMP.shp'))
arcpy.MakeFeatureLayer_management( 'river_cells_TMP.shp','rivcells')
if arcpy.Exists(os.path.join(path,'river_cells_dissolve_TMP.shp')):
    arcpy.Delete_management(os.path.join(path,'river_cells_dissolve_TMP.shp'))
arcpy.Copy_management(os.path.join(path,SFRdata.CELLS_DISS),os.path.join(path,'river_cells_dissolve_TMP.shp'))
arcpy.MakeFeatureLayer_management('river_cells_dissolve_TMP.shp','rivcellsdiss')
# opening river_explode.shp as read only so no need for scratch file
arcpy.MakeFeatureLayer_management(SFRdata.intersect,'rivexplode')

# clip the river cells to the remaining streams and re-save
print '\ninteresecting river_cells with river_explode...'
rivtrim = arcpy.SelectLayerByLocation_management('rivcells',"INTERSECT",'rivexplode')
print '...Saving down results...'
if arcpy.Exists(os.path.join(path,SFRdata.CELLS)):
    print 'first removing old version of {0:s}'.format(SFRdata.CELLS)
    arcpy.Delete_management(os.path.join(path,SFRdata.CELLS))
arcpy.CopyFeatures_management(rivtrim,os.path.join(path,'river_cells'))
arcpy.Delete_management('rivcells')
arcpy.Delete_management(os.path.join(path,'river_cells_TMP.shp'))
# do the same for the dissolved river cells layer
print 'interesecting river_cells_dissolve with river_explode...\n\n'
rivdisstrim = arcpy.SelectLayerByLocation_management('rivcellsdiss',"INTERSECT",'rivexplode')
print '...Saving down results...'
if arcpy.Exists(os.path.join(path,SFRdata.CELLS_DISS)):
    print 'first removing old version of {0:s}'.format(SFRdata.CELLS_DISS)
arcpy.Delete_management(os.path.join(path,SFRdata.CELLS_DISS))
arcpy.CopyFeatures_management(rivdisstrim,os.path.join(path,'river_cells_dissolve'))
arcpy.Delete_management('rivcellsdiss')
arcpy.Delete_management(os.path.join(path,'river_cells_dissolve_TMP.shp'))

# make a backup copy of fix_comids.txt
print 'making a copy of "boundary_manual_fix_issues.txt" --> "boundary_manual_fix_issues_backup.txt"'
if os.path.exists('boundary_manual_fix_issues.txt'):
    shutil.copy('boundary_manual_fix_issues.txt','boundary_manual_fix_issues_backup.txt')
# make a backup copy of fix_comids.txt
print 'making a copy of "fix_comids.txt" --> "fix_comids_backup.txt"'
if os.path.exists('fix_comids.txt'):
    shutil.copy('fix_comids.txt','fix_comids_backup.txt')
