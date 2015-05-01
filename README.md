SFRmaker
==
Python package to facilitate use of the streamflow-routing (SFR) package in MODFLOW. Currently, SFRmaker consists of three primary modules:  
  
####sfr_classes.py
Intersects linework with model grid and develops routed segments and reaches from information in NHDPlus. Produces two tables, **Mat1** and **Mat2**, which contain SFR package information by reach and by segment, respectively. Also produces a shapefile of linework broken by the grid, which contains the stream geometries represented by the SFR reaches.

####postproc.py
Postproc operates on **Mat1** and **Mat2** as pandas dataframe attributes, to produce and visualize an SFR package. Includes methods to:  

*  sample streambed top elevations from a DEM, and smooth them so they decrease monotonically downstream
*  incorporate field measurements of streambed elevation
*  visualize routing by outlet
*  visualize stream profiles (sequences of SFR segments) from headwaters to outlet
*  adjust SFR streambed conductance so that only one reach per cell has conductance
*  adjust the MODFLOW discretization file so that the geometry of layer 1 is consistent with the SFR elevations
*  read in SFR results from a MODFLOW run, allowing for interactive plotting of SFR stages and mass balance components.
*  write shapefiles of the SFR package input and results
*  write an SFR package file and updated versions of Mat1 and Mat2

####diagnostics.py
Contains methods to check for common SFR problems such as  

* continuity in segment/reach numbering
* circular routing, spatial gaps in routing, breaks in routing (outlets in the interior of the model)
* multiple SFR conductances in a single cell (can lead to circular routing of water between SFR reaches)
* elevations that rise in the downstream direction (either within a segment, or segment ends that rise)
* slope values below the a user-defined minimum (can lead to artifically high stages)



## Dependencies:

SFRmaker runs in Python 2.7. The Anaconda Scientific Python Distribution (<https://store.continuum.io/cshop/anaconda/>) is available for free, and provides an easy way for managing python packages through its package manager **conda**. Additional dependencies are listed by module below:

####sfr_classes  
* **ESRI Arcpy**, which **must be added to the path of your main python distribution**. See instructions below.
* **flopy**, available via **pip**, or at <https://github.com/modflowpy/flopy>

####postproc  
* **fiona**
* **shapely** (fiona and shapely are available via **conda** or **pip**; instructions on installing pip and these packages are available here: <https://github.com/aleaf/LinesinkMaker/blob/master/README.md>, under **"Installing the required Python packages""**)
* **rasterstats** (also available via **pip**, see <https://github.com/perrygeo/python-raster-stats> for more info)
* **flopy**
* **GIS_utils** (see <https://github.com/aleaf/GIS_utils>)




  

####Adding Arcpy to the python path:  
1) **Make a file called: Desktop10.pth**, with the following lines:
```
C:\ArcGIS\Desktop10.1\arcpy  
C:\ArcGIS\Desktop10.1\bin64  
C:\ArcGIS\Desktop10.1\ArcToolbox\Script
```
Notes:  
The second line may be "bin" for 32 bit or "bin64" for 64 bit.  
If you are using ArcMap 10.0 or 10.2, "Desktop10.1" in the above path needs to be modified accordingly.

2) **Place this file in your python path where all your site-packages are installed**. For example, for users of the Enthought Canopy Distribution, the file would need to be placed at:
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
* row and column information must be supplied as attributes to the shapefile, named "row" and "column" respectively  
* polygon shapefile of model domain (can be created by dissolving the grid shapefile using *Geoprocessing>dissolve* in the ArcToolbox). **Note:** this polygon (it must be a shapefile polygon, not a line) defines the area where streams are represented by SFR within the model. If SFR in only a subset (i.e. nearfield) area of the model is desired, then this polygon is not the same as the model boundary.
* discretization file for model
  
## Outputs:
* SFR package file
* Text file table with reach information (r,c,l,elevation, width, slope, etc.) 
* Text file table with segment/routing data (e.g. icalc, outseg, iupseg, etc.)  
  
  
## Workflow for building SFR input:
  
1) **Setup XML input file** (see EXAMPLE.XML in \<InputFiles> section) to point to the above input datasets  
  
* check settings in \<GlobalSettings> section; make sure that \<preproc> is set to **True**  
* select an elevation method under \<elevflag

2) Make sure that the "infile" variable in **SFR_main.py** (see SFR_main_EXAMPLE.py) points to the XML input file. Also make sure that the calls to classes are entered correctly.

3) **Run SFR_main.py** by typing *python SFR_main.py* at the command prompt  

4) If a "manual intervention" message is encountered in screen output,  

* **check the following files:  **
	* **fix_com_IDs.txt:** Lists stream segments that have multiple ends and starts; usually these are 		streams that are broken into mulitple parts by the grid boundary. 
	* **boundary_manual_fix_issues.txt:** Lists stream segments that don't have any connections to other 		segments.  
* **edit the offending segments** (COMIDs) in the shapefile specified under \<IntermediateFiles>\<intersect> in the XML input file (usually this means deleting the parts of multi-part COMIDs that are isolated by the grid boundary, and possibly deleting any unconnected segments).  
  
5) set \<preproc> in the XML input file to False (meaning the existing \<intersect> shapefile will be read in lieu of the preprocessing operations). Then **rerun SFR_main.py**.  
6) Once the reach and segment information tables have been written, they can be subsequently edited, and an SFR package file rebuilt by running the **Assign_layers.py** script. Since Assign_layers also reads the MODFLOW Discretizaiton file, it can aslo re-assign layering for SFR cells following any modifications to the grid elevations.
 
## Visualizating the results:  
#####Obtaining SFR output  
Prior to running your model, make sure that variable ISTCB2 in the SFR package file (see SFR manual) is set to a positive unit number, which is also listed in the MODFLOW NAM file, in association with an ouput text file name. The value of this variable can be specified in the XML input file to these python scripts, but may be overwritten if a pre-processor GUI (i.e. GW Vistas) is used to generate MODFLOW input. The text file generated by this setting has information on streamflow and steam/aquifer exchange which is needed for the visualization tools below.  

#####Plotting streambed profiles
SFR_main_EXAMPLE.py has examples for producing comparison plots of streambed profiles for the different elevation methods (NHDPlus, DEM, and topographic contours), and also for plotting profiles of the final streambed elevations (segments and reaches have been created). Both of these methods need to be run in the MAIN program, in an order similar to that shown in the example MAIN file.  

#####Visualizing routing
SFR_routing_checker.py can be run independently of the MAIN program to visualize routing. Simply edit the SFR_routing_checker.XML input file, and run by typing *python SFR_routing_checker.py SFR_routing_checker.XML* at the command prompt. Requires an SFR package file, and a grid spec. (SPC) file (written by GWV). To create a grid spec. file, in Groundwater Vistas, go to Model \> PEST \> Create Grid Spec. File. SFR_routing_checker works best with \<all_layers> set to False. To visualize shapefile output in ArcMap, after importing, under Properties>Symbology choose categories and click "Add all values".

#####Visualizing streamflow and aquifer interactions  
Edit the plot_SFR_flows.py example. Requires a MODFLOW DIS file, an "exploded" stream linework file that has stream fragments by model cell number (i.e. \<intersect> in the XML input file, or similar), and an SFR package output textfile (i.e. "*streamflow.dat"; this is the same file discussed in the "Obtaining SFR output" section above). Produces a shapefile of the same name as the SFR output file.
To view in Arc, after importing, under Properties>Symbology, click Import and choose:  

* **SFR_flow_symbology.lyr** to plot flow by line thickness
* **SFR_interactions_symbology.lyr** to plot gaining, loosing, and dry segments by color
* **SFR_interactions_graduated_symbology.lyr** to plot stream/aquifer interactions as graduated colors (i.e. dark blue for largest gains, gray for dry, red for largest losses)
