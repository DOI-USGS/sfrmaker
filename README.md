SFR
===
Set of programs for automating the construction of the MODFLOW Streamflow-Routing Package, using NHD Plus v2


#### Dependencies:

In addition to standard Python modules, ESRI Arcpy is needed.
Some input and output are designed for the Groundwater Vistas GUI, but this is not necessarily required

#### Input requirements:

1) NHD Plus v2 hydrography datasets:
	* Available at http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php
	* NHDPlusV21_SR_**_NHDSnapshot_**.7z
		* NHDFcode.dbf
		* NHDFlowline.dbf, .prj, .shp, .shx
	* NHDPlusV21_SR_**_NHDPlusAttributes_**.7z
		* elevslope.dbf
		* PlusFlow.dbf
		* PlusFlowlineVAA.dbf

	* If model domain contains a major divide (i.e. multiple HUC2 hydrologic regions), need to merge relevant NHD datasets (e.g. 04 and 07) prior to running this script
 

2) DEM and/or topographic contours for area of model domain:
	* Available from the National Map Viewer and Download Platform: (http://viewer.nationalmap.gov/viewer/)
	* In "Overlays" tab, select "Elevation Availability" to view available DEM resolutions
	* click "Download Data" link in upper right to download DEM(s)
	* selected "elevation" and/or "contours" to download
	
	* Internal WIWSC source for DEMs in Wisconsin: 		\\igsarmewfsapa\GISData\BaseData\Wisconsin\Statewide\Elevation\NED_2013_10m\NED_10m_DEM_WI_UP.gdb

Note: If model domain area has multiple DEMs or elevation contour shapefiles, they need to be merged prior to setting up SFR.

3) Model grid information:
	* polygon shapefile export of model grid
	* discretization file for model

Download elevation contour data from the National Map
		
* Polygon shapefile export of model grid (e.g. as exported by Groundwater Vistas)
* rows x columns ascii matrix of model TOP elevations (e.g. as exported by Groundwater Vistas)
* Shapefile polygon of model domain (merged polygon of gridcells)
* A DEM for model area, doesn't need to be clipped
* PlusflowVAA database from NHDPlus v2 --> PlusFlowlineVAA.dbf from NHDPlusV21\_XX\_YY\_NHDPlusAttributes\_03.7z
* Elevslope database from NHDPlus v2 --> elevslope.dbf from NHDPlusV21\_XX\_YY\_NHDPlusAttributes\_03.7z
* NHDFlowline shapefile from NHDPlus v2
* Flowlines='Flowlines.shp' # from NHDPlus
* PlusFlow database from NHDPlus v2 --> PlusFlow.dbf from NHDPlusV21\_XX\_YY\_NHDPlusAttributes\_03.7z
* NHDFcode database from NHDPlus v2

NOTE: XX is Drainage Area ID and YY is VPU (vector processing unit) in the above (see NHDPlus website for details).

#### Outputs:

* SFR package file*
* Text file with reach information (r,c,l,elevation, width, slope, etc.) for importing SFR package into Groundwater Vistas
* Text file with segment/routing data (e.g. icalc, outseg, iupseg, etc.), for copying and pasting into Groundwater Vistas

#### Notes:

* All shps should be in (same) model coordinate units (e.g. ft.)
* If model domain contains a major divide, need to merge relevant NHD datasets (e.g. 04 and 07) prior to running this script


### Workflow for building SFR input:
