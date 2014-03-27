SFR
==
Set of programs for automating the construction of the MODFLOW Streamflow-Routing Package, using NHD Plus v2


## Dependencies:

In addition to standard Python modules, ESRI Arcpy is required, and **must be added to the path of your main python distribution**. To do that:  
  
1) Make a file called: Desktop10.pth, with the following lines:
```
C:\ArcGIS\Desktop10.1\arcpy  
C:\ArcGIS\Desktop10.1\bin64  
C:\ArcGIS\Desktop10.1\ArcToolbox\Script
```
Notes:  
The second line may be "bin" for 32 bit or "bin64" for 64 bit.  
If you are using ArcMap 10.0 or 10.2, "Desktop10.1" in the above path needs to be modified accordingly.

2) Place this file in your python path where all your site-packages are installed. For example, for users of the Enthought Canopy Distribution, the file would need to be placed at:
`C:\Users\<username>\AppData\Local\Enthought\Canopy\User\Lib\site-packages\`



## Input requirements:

#####1) NHDPlus v2 hydrography datasets  
 * Available at <http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php>
 * Archives needed/relevant files:
 	* **NHDPlusV21_XX_YY_NHDSnapshot_**.7z**   
 		* NHDFcode.dbf  
 		* NHDFlowline.dbf, .prj, .shp, .shx  
 	* **NHDPlusV21\_XX\_YY\_NHDPlusAttributes\_**.7z**  
 		* elevslope.dbf  
		* PlusFlow.dbf  
		* PlusFlowlineVAA.dbf
	* If model domain contains a major divide (i.e. multiple HUC2 hydrologic regions), need to merge relevant NHD datasets (e.g. 04 and 07) prior to running this script

**Notes:**  

* XX is Drainage Area ID (e.g., GL for Great Lakes) and YY is the Vector Processing Unit (VPU; e.g. 04) in the  above (see NHDPlus website for details).  
* If model domain encompases multiple Drainage Areas, need to merge NHD datasets prior to running this script (e.g. in ArcMap, using the Merge tool under *ArcToolbox>Data Management>General>Merge*)  
* The NHDFlowline shapefile should be reprojected into the same coordinate system of the model, if necessary. In ArcMap this can be done under *ArcToolbox>Data Management>Projections and Transformations*.
* **Edit the NHD datasets at your own risk!!** Editing may corrupt routing or other information, greatly complicating automated setup of an SFR package. NHD elevation units should be left as is (cm).
 
#####2) DEM and/or topographic contours for area of model domain  
* Available from the National Map Viewer and Download Platform: <http://viewer.nationalmap.gov/viewer/>  
	* In "Overlays" tab, select "Elevation Availability" to view available DEM resolutions  
	* click "Download Data" link in upper right to download DEM(s)  
	* select "elevation" and/or "contours" to download  
* Internal WIWSC source for DEMs in Wisconsin:  
  `\\igsarmewfsapa\GISData\BaseData\Wisconsin\Statewide\Elevation\NED_2013_10m\NED_10m_DEM_WI_UP.gdb`  
  
**Notes:**  

* If model domain area has multiple DEMs or elevation contour shapefiles, they need to be merged prior to setting up SFR. The merged DEM should be in the same coordinate system as the model.  

#####3) Model grid information
* polygon shapefile export of model grid (in same coordinate system as Flowlines and DEM)
* polygon shapefile of model domain (can be created by dissolving the grid shapefile using *Geoprocessing>dissolve* in the ArcToolbox)
* discretization file for model
  
## Outputs:
* SFR package file
* Text file table with reach information (r,c,l,elevation, width, slope, etc.) 
* Text file table with segment/routing data (e.g. icalc, outseg, iupseg, etc.)  
  
  
## Workflow for building SFR input:
  
1) **Setup XML input file** (see EXAMPLE.XML in \<InputFiles> section) to point to the above input datasets  
  
* check settings in \<GlobalSettings> section; make sure that \<preproc> is set to **True**  

2) Make sure that the "infile" variable in **SFR_main.py** (see SFR_main_EXAMPLE.py) points to the XML input file  

3) **Run SFR_main.py** by typing *python SFR_main.py* at the command prompt  

4) If a "manual intervention" message is encountered in screen output,  

* **check the following files:  **
	* **fix_com_IDs.txt:** Lists stream segments that have multiple ends and starts; usually these are 		streams that are broken into mulitple parts by the grid boundary. 
	* **boundary_manual_fix_issues.txt:** Lists stream segments that don't have any connections to other 		segments.  
* **edit the offending segments** (COMIDs) in the shapefile specified under \<IntermediateFiles>\<intersect> in the XML input file (usually this means deleting the parts of multi-part COMIDs that are isolated by the grid boundary, and possibly deleting any unconnected segments).  
  
5) set \<preproc> in the XML input file to False (meaning the existing \<intersect> shapefile will be read in lieu of the preprocessing operations). Then **rerun SFR_main.py**.  
6) Once the reach and segment information tables have been written, they can be subsequently edited, and an SFR package file rebuilt by running the **Assign_layers.py** script. Since Assign_layers also reads the MODFLOW Discretizaiton file, it can aslo re-assign layering for SFR cells following any modifications to the grid elevations.
 
## Visualizating the results:
#####Plotting streambed profiles
SFR_main_EXAMPLE.py has examples for producing comparison plots of streambed profiles for the different elevation methods (NHDPlus, DEM, and topographic contours), and also for plotting profiles of the final streambed elevations (segments and reaches have been created). Both of these methods need to be run in the MAIN program, in an order similar to that shown in the example MAIN file.
#####Visualizing routing
SFR_routing_checker.py can be run independently of the MAIN program to visualize routing. Simply edit the SFR_routing_checker.XML input file, and run by typing *python SFR_routing_checker.py SFR_routing_checker.XML* at the command prompt. Requires an SFR package file, and a grid spec. (SPC) file (written by GWV). Works best with \<all_layers> set to False. To visualize shapefile output in ArcMap, after importing, under Properties>Symbology choose categories and click "Add all values".

#####Visualizing streamflow and aquifer interactions  
Edit the plot_SFR_flows.py example. Requires a MODFLOW DIS file, an "exploded" stream linework file that has stream fragments by model cell number (i.e. \<intersect> in the XML input file, or similar), and an SFR package output file (i.e. "*streamflow.dat"). Produces a shapefile of the same name as the SFR output file.
To view in Arc, after importing, under Properties>Symbology, click Import and choose:  

* **SFR_flow_symbology.lyr** to plot flow by line thickness
* **SFR_interactions_symbology.lyr** to plot gaining, loosing, and dry segments by color
* **SFR_interactions_graduated_symbology.lyr** to plot stream/aquifer interactions as graduated colors (i.e. dark blue for largest gains, gray for dry, red for largest losses)
