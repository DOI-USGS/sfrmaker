SFRmaker
=
Python package to facilitate use of the streamflow-routing (SFR) package in MODFLOW.  
###Disclaimer:
This package is provided as-is, and is not an official USGS product. While the authors are continuing to refine and improve the code as time allows, the project remains unfunded and for now, largerly undocumented and somewhat held together by "chewing gum and duct tape." If used properly, the code can in most cases produce a robust, reproducable SFR package in a fraction of the time that would be required for manual construction in a GUI. Some diagnostic and visualization tools are included; however, **it is up to the modeler to check their own results for correctness, using their own tests if necessary.** Some methods in the code are better than others, some are dead-ends, and other methods were written for specific situations encountered in individual projects. In lieu of comprehensive documentation, the best approach is to follow the instructions outlined in this readme file, as well as the IPython Notebook examples linked to below.

###Code structure
Currently, SFRmaker consists of four primary modules. Either sfr_classes or preproc can be used for the preprocessing. See the **Workflow** section further down.
  
####sfr_classes.py
This is essentially the original SFRmaker written by Reeves, Fienen and Leaf in 2014.
Intersects linework with model grid and develops routed segments and reaches from information in NHDPlus. Produces two tables, **Mat1** and **Mat2**, which contain SFR package information by reach and by segment, respectively. Also produces a shapefile of linework broken by the grid, which contains the stream geometries represented by the SFR reaches. **Requires Arcpy**.

####preproc.py
This is written as a replacement to sfr_classes, using the [python GIS toolchain](https://github.com/aleaf/SFRmaker/blob/master/pythonGIS_install_readme.md "python GIS toolchain") instead of ArcPy. Accepts NHDPlus datasets and a shapefile of the model grid as input, and produces Mat1 (reach) and Mat2 (segment) tables as output, as well as a shapefile of the linework represented by each SFR reach, with segment and reach information. The speed, flexibility and simplicity of the python GIS packages allow for several improvements over sfr_classes, including  
   
* better performance (minutes vs. hours for very large networks)
* no intermediate files
* less code
* a one-to-one correspondance between NHDPlus COMIDs and SFR segments, which is documented in the shapefile output, allowing for easier referencing to the original linestring geometries and other auxillary information on stream reahces/segments
* automatic reprojection of input datasets
* no Arcpy requirement


Preproc should also fix several known issues with sfr_classes (although it is still being tested!):  

* routing skipping segments near confluences (unless editing of the input linework has resulted in the removal of COMIDs in the stream sequence)
* segment numbering decreasing in the downstream direction, which can slow convergance with MODFLOW-NWT in some cases.

Preproc can be used instead of sfr_classes to set up the SFR segments and reaches. The **example_preproc_postproc.py** script in the examples folder shows a complete workflow from input datasets to SFR package input file, using the preproc and postproc modules.


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

#####SFRmaker runs in Python 2 or 3
The Anaconda Scientific Python Distribution (<https://store.continuum.io/cshop/anaconda/>) is available for free, and provides an easy way for managing python packages through its package manager **conda**. Additional dependencies are listed by module below:

#####flopy
available via **pip** (see instructions below), or at <https://github.com/modflowpy/flopy> 

#####GIS packages
 The postproc module depends on a collection of packages (**fiona, shapely, gdal, pyproj, rasterio, rasterstats**, and **GIS_utils**) that provide python bindings to open source GIS libraries. Instructions for installing these packages on Windows can be found here:
 <https://github.com/aleaf/SFRmaker/blob/master/pythonGIS_install_readme.md>



**GIS_utils** is just a collection of macros to facilitate basic GIS tasks. If SFRmaker ever gets made into a coherent package, 	GIS_utils will likely be integrated. In the meantime, it can be installed by running  

	```
	>pip install https://github.com/aleaf/GIS_utils/archive/master.zip
	```

#### Note
for running the **sfr_classes / ArcPy** version of SFRmaker see
<https://github.com/aleaf/SFRmaker/blob/master/SFRmaker_arcpy_workflow.md>


  





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

#### For a complete workflow see the examples:
<http://nbviewer.ipython.org/github/aleaf/SFRmaker/blob/master/Examples/example_SFRmaker_script.py>  

<https://github.com/aleaf/SFRmaker/blob/master/Examples/example_Notebook.ipynb>


As shown in the above examples, it is a good idea to run diagnostics using the Flopy checker. This will perform a lot of checks to help avoid common problems with SFR setup.

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
