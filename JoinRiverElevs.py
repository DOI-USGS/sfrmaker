# file to join rivers with elevations
# Fienen --> 8/2013


import arcpy


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

# MNF added the joining of elevations to river_explode.shp as indicated in the notes. Saves as ELEV file
arcpy.MakeTableView_management("river_elevs.dbf",'elevations') # this is like bringing something into the table of contents in ArcMap
arcpy.MakeFeatureLayer_management('river_explode.shp',"rivexp") 
oldfid = getfield('elevations','OLDFID')
fid = getfield('rivexp','fid')
print 'joining new elevations (river_elevs.dbf) to trimmed river lines (river_explode.shp)'
arcpy.JoinField_management('rivexp',fid,'elevations',oldfid)
print 'Saving joined results to --> %s' %(ELEV)
if arcpy.Exists(ELEV):
    arcpy.Delete_management(ELEV)
arcpy.CopyFeatures_management('rivexp',ELEV)