SFRmaker
=
Python package to facilitate use of the streamflow-routing (SFR) package in MODFLOW.  
###Disclaimer:
This package is provided as-is, and is not an official USGS product. While the authors are continuing to refine and improve the code as time allows, the project remains unfunded and for now, largerly undocumented and somewhat held together by "chewing gum and duct tape." If used properly, the code can in most cases produce a robust, reproducable SFR package in a fraction of the time that would be required for manual construction in a GUI. Some diagnostic and visualization tools are included; however, **it is up to the modeler to check their own results for correctness, using their own tests if necessary.** Some methods in the code are better than others, some are dead-ends, and other methods were written for specific situations encountered in individual projects. In lieu of comprehensive documentation, the best approach is to follow the instructions outlined in this readme file, as well as the IPython Notebook examples linked to below.

###Code structure
Currently, SFRmaker consists of three primary modules:  
  
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



### Dependencies:

SFRmaker runs in Python 2.7. The Anaconda Scientific Python Distribution (<https://store.continuum.io/cshop/anaconda/>) is available for free, and provides an easy way for managing python packages through its package manager **conda**. Additional dependencies are listed by module below:

####sfr_classes  
* #####ESRI Arcpy
which **must be added to the path of your main python distribution**.

	**Adding Arcpy to the python path**:  	
	1) **Make a file called: Desktop10.pth**, with the following lines:  

	```
	C:\ArcGIS\Desktop10.2\arcpy  
	C:\ArcGIS\Desktop10.2\bin64  
	C:\ArcGIS\Desktop10.2\ArcToolbox\Scripts
	```
Notes:  
These lines tell python where to find arcpy and its associated packages/libraries. The second line may be "bin" for 32 bit or "bin64" for 64 bit. If you are using ArcMap 10.0 or 10.1, "Desktop10.2" in the above path needs to be modified accordingly.

	2) **Place this file where python can find it.** For Anaconda on Windows 7, this path should work (replacing ```aleaf``` with your username):

		C:\Users\aleaf\AppData\Local\Continuum\Anaconda
		
	The **Lib** subfolder in this folder or the **site-packages** folder within Lib may also work. Python checks these folders by default when looking for modules; within these folders, it looks for files with the extension **.pth**, and checks for additional paths within those files (i.e., the arcpy paths listed above).
	
* ####flopy
available via **pip**, or at <https://github.com/modflowpy/flopy> 

####postproc  
* **fiona**
* **shapely** (fiona and shapely are available via **conda** or **pip**; instructions on installing pip and these packages are available here: <https://github.com/aleaf/LinesinkMaker/blob/master/README.md>, under **"Installing the required Python packages""**)
* **rasterstats** (also available via **pip**, see <https://github.com/perrygeo/python-raster-stats> for more info)
* **flopy**
* **GIS_utils** (see <https://github.com/aleaf/GIS_utils>)




  





## Input data requirements

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

	* XX is Drainage Area ID (e.g., GL for Great Lakes) and YY is the Vector Processing Unit (VPU; e.g. 04) in the  above (see NHDPlus 	website for details).  
	* If model domain encompases multiple Drainage Areas, need to merge NHD datasets prior to running this script (e.g. in ArcMap, using the Merge tool under *ArcToolbox>Data Management>General>Merge*)  
	* The NHDFlowline shapefile should be reprojected into the same coordinate system of the model, if necessary. In ArcMap this can be done under *ArcToolbox>Data Management>Projections and Transformations*.
	* **Edit the NHD datasets at your own risk!!** Editing may corrupt routing or other information, greatly complicating automated setup of an SFR package. NHD elevation units should be left as is (cm).
 
#####2) DEM  
Used for deriving streambed top elevations for individual reaches (below the 100,000k scale of NHDPlus)  
**Note:** SFRmaker doesn’t edit the DEM. The general approach is to sample elevations from the DEM and then smooth them so that they decrease in the downstream direction. See workflow instructions below.

* Available from the National Map Viewer and Download Platform: <http://viewer.nationalmap.gov/viewer/>  
* In "Overlays" tab, select "Elevation Availability" to view available DEM resolutions  
* click "Download Data" link in upper right to download DEM(s)  	
* select "elevation" and/or "contours" to download  

  
	**Notes:**  
	
	If model domain area has multiple DEMs or elevation contour shapefiles, they need to be merged prior to setting up SFR. The merged 	DEM should be in the same coordinate system as the model.

#####3) Model grid information
* polygon shapefile export of model grid (in same coordinate system as Flowlines and DEM)  
* row and column information must be supplied as attributes to the shapefile, named "row" and "column" respectively  
* polygon shapefile of model domain (can be created by dissolving the grid shapefile using *Geoprocessing>dissolve* in the ArcToolbox). **Note:** this polygon (it must be a shapefile polygon, not a line) defines the area where streams are represented by SFR within the model. If SFR in only a subset (i.e. nearfield) area of the model is desired, then this polygon is not the same as the model boundary.
* discretization file for model
  
## Outputs:
* SFR package file
* **Mat1** (Text file table with reach information (r,c,l,elevation, width, slope, etc.))
* **Mat2** (Text file table with segment/routing data (e.g. icalc, outseg, iupseg, etc.) ) 
  
  
## Workflow for building SFR package:

####Run sfr_classes
With the XML input file, the SFR_main.py script (calling the methods in sfr_classes) will produce two tables, **Mat1** and **Mat2**, with SFR reach and segment information.
  
1) **Setup XML input file** (see EXAMPLE.XML in \<InputFiles> section) to point to the above input datasets  
  
* check settings in \<GlobalSettings> section; make sure that \<preproc> is set to **True**  
* Set \<elevflag> to NHDPlus (this will set segment end elevations from NHDPlus COMID min/max elevations in the **elevslope** table)

2) Make sure that the "infile" variable in **SFR_main.py** (see SFR_main_EXAMPLE.py) points to the XML input file. Also make sure that the calls to classes are entered correctly.

3) **Run SFR_main.py** by typing *python SFR_main.py* at the command prompt  

4) If a "manual intervention" message is encountered in screen output,  

* **check the following files:  **
	* **fix_com_IDs.txt:** Lists stream segments that have multiple ends and starts; usually these are 		streams that are broken into mulitple parts by the grid boundary. 
	* **boundary_manual_fix_issues.txt:** Lists stream segments that don't have any connections to other 		segments.  
* **edit the offending segments** (COMIDs) in the shapefile specified under \<IntermediateFiles>\<intersect> in the XML input file (usually this means deleting the parts of multi-part COMIDs that are isolated by the grid boundary, and possibly deleting any unconnected segments).  
  
5) set \<preproc> in the XML input file to False (meaning the existing \<intersect> shapefile will be read in lieu of the preprocessing operations). Then **rerun SFR_main.py**.

####Run postproc and diagnostics 
1) Once the reach and segment information tables (**Mat1** and **Mat2**) have been written, the methods in **postproc.py** can be run to address all of the issues listed under the **postproc** and **diagnostics** modules above. An example workflow for **postproc** is given here:  
  
<http://nbviewer.ipython.org/github/aleaf/SFRmaker/blob/master/Examples/Example_postproc_workflow.ipynb>  

Streambed top elevations are set by two methods in postproc.py:  
* **reset_m1_streambed_top_from_dem()** will sample the minimum DEM elevation for each SFR reach, using zonal statistics.  
* **smooth_segment_interiors()** will take the sampled minimum elevations, and enforce that they decrease in the downstream direction. The segment end elevations from NHDPlus are retained, ensuring that segment elevations do not rise from start to end. 

As shown in the above example, it is a good idea to run the diagnostics module. This will perform a lot of checks to help avoid common problems with SFR setup.

## Visualizating the results:  


#####Obtaining SFR output  
Prior to running your model, make sure that variable ISTCB2 in the SFR package file (see SFR manual) is set to a positive unit number, which is also listed in the MODFLOW NAM file, in association with an ouput text file name. The text file generated by this setting has information on streamflow and steam/aquifer exchange which is needed for the visualization tools below.  

#####Plotting streambed profiles
An example of reading in SFR results from the output text file is given here:
<http://nbviewer.ipython.org/github/aleaf/SFRmaker/blob/master/Examples/SFR_results_visualization_example.ipynb>
The example includes:  

* using **SFRmaker** to read the SFR output, and write it to a shapefile
* plotting profiles of streamflow, stream-aquifer interactions, and stage for a selected segment
* reading in the MODFLOW dis file using **flopy**, and comparing grid cell elevations to profiles of stage  
* using **flopy** to read SFR information from the MODFLOW cell-by-cell budget (**cbb**) file

#####Visualizing routing
Routing can be visualized in a shapefile, as produced in the first (postproc workflow) example above. SFR reaches in the shapefile are attributed by segment, reach, outseg, and outlet.

#####Visualizing streamflow and aquifer interactions  
The shapefile of SFR output (including streamflow, stream-aquifer interactions, etc.) can be visualized in Arcmap with the symbology contained in the layer files below (available in the **Symbology** folder).

To view in Arc, after importing the shapefile, under Properties>Symbology, click Import and choose:  

* **SFR_flow_symbology.lyr** to plot flow by line thickness
* **SFR_interactions_symbology.lyr** to plot gaining, loosing, and dry segments by color
* **SFR_interactions_graduated_symbology.lyr** to plot stream/aquifer interactions as graduated colors (i.e. dark blue for largest gains, gray for dry, red for largest losses)
