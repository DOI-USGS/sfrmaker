# SFR-related utilities for arcpy
# --> a set of macros to make more efficient operations
#     and the summarize some common sets of operations
#
#  a m!ke@usgs joint --- mnfienen@usgs.gov
#

import arcpy


def general_join(target_name, target_lay, joinfield1, joined, joinfield2, keep_common=True):
    # input:
    # target_name: Name of the shapefile to save the joined data to
    # target_lay: this is a layer as made with arcpy.MakeFeatureLayer_management 
    # joinfield1: the name of the field in target_lay to join on
    # joined: this is a filename for the table/shapefile to join to
    # joinfield2: the name of the field in joined to join on
    # keep_common (optional): True means use KEEP_COMMON, False means KEEP_ALL
    
    # need the next line to make sure that column names join properly
    arcpy.env.qualifiedFieldNames = False
    
    if keep_common:
        join_type = "KEEP_COMMON"
    else:
        join_type = "KEEP_ALL"
    print "\nJoining field {0:s} of {1:s} to field {2:s} of {3:s}".format(
              joinfield2, target_lay, joinfield2, joined)
    arcpy.AddJoin_management(target_lay,joinfield1,joined,joinfield2,join_type)
    # save back down the results
    if arcpy.Exists('tmpjunkus.shp'):
        print 'Removing old version of tmpjunkus.shp'
        print 'This is a holding temporary file to save down %s' %target_name
        print 'tmpjunkus.shp will be deleted'
        arcpy.Delete_management('tmpjunkus.shp')
    arcpy.CopyFeatures_management(target_lay,'tmpjunkus.shp')
    if arcpy.Exists(target_name):
        print 'Removing old version of %s' %target_name
        arcpy.Delete_management(target_name)
    arcpy.Rename_management('tmpjunkus.shp',target_name)
