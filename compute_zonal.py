# Program to sample a raster using the zonal statistics tool in Arc Toolbox

import arcpy
from collections import defaultdict
import numpy as np

# Input file
'''
infile="SFR_setup.in"

# Get input parameters
infile=open(infile,'r').readlines()
inputs=defaultdict()
for line in infile:
    if "#" not in line.split()[0]:
        varname=line.split("=")[0]
        var=line.split("=")[1].split()[0]
        inputs[varname]=var.replace("'","")
'''
nrows,ncolumns=800,800
delx=250
MFgrid='BadRiver_cells.shp'
path='D:\\ATLData\\BadRiver\\Grid\\'
DEM='ned_ft_utm'

# Settings
arcpy.env.workspace = path
arcpy.env.overwriteOutput = True 

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")


#arcpy.AddField_management(MFgrid,"cellnum","LONG") #isn't needed if there is a node column

print "extracting model top elevations from DEM..."
# create raster of gridcells based on node (cellnum)
# run zonal statistics on DEM for each grid cell
print "\tbuilding raster of model grid cells for sampling DEM..."
DEMres=arcpy.GetRasterProperties_management(DEM,"CELLSIZEX")
DEMres=float(DEMres.getOutput(0)) # cell size of DEM in project units
arcpy.PolygonToRaster_conversion(MFgrid,"node",path+"MFgrid_raster","MAXIMUM_AREA","node",DEMres)
print "\tbuilding raster attribute table..."
arcpy.BuildRasterAttributeTable_management(path+"MFgrid_raster")
print "\trunning zonal statistics... (this step will likely take several minutes or more)"
arcpy.sa.ZonalStatisticsAsTable(path+"MFgrid_raster","VALUE",DEM,path+"demzstats.dbf","DATA","ALL")