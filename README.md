SFR
===
Set of programs for automating the construction of the MODFLOW Streamflow-Routing Package, using NHD Plus v2 and a digital elevation model (DEM).
NHDPlus datasets are available at: http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php


Dependencies:
In addition to standard Python modules, ESRI Arcpy is needed.
Some input and output are designed for the Groundwater Vistas GUI, but this is not necessarily required

Input requirements:

Polygon shapefile export of model grid (e.g. as exported by Groundwater Vistas)
rows x columns ascii matrix of model TOP elevations (e.g. as exported by Groundwater Vistas)
Shapefile polygon of model domain (merged polygon of gridcells)
A DEM for model area, doesn't need to be clipped
PlusflowVAA database from NHDPlus v2
Elevslope database from NHDPlus v2
NHDFlowline shapefile from NHDPlus v2
Flowlines='Flowlines.shp' # from NHDPlus
PlusFlow database from NHDPlus v2
NHDFcode database from NHDPlus v2


Outputs:
SFR package file*
Text file with reach information (r,c,l,elevation, width, slope, etc.) for importing SFR package into Groundwater Vistas
Text file with segment/routing data (e.g. icalc, outseg, iupseg, etc.), for copying and pasting into Groundwater Vistas
