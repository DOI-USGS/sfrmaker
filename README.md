SFR
===
Set of programs for automating the construction of the MODFLOW Streamflow-Routing Package, using NHD Plus v2 and a digital elevation model (DEM).
NHDPlus datasets are available at: http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php


#### Dependencies:

In addition to standard Python modules, ESRI Arcpy is needed.
Some input and output are designed for the Groundwater Vistas GUI, but this is not necessarily required

#### Input requirements:

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

* currently the SFR package file is create prior to final elevation corrections with fix_w_dem.py. Until this is fixed, SFR package input should be manually created from the Groundwater Vistas input tables.
Note:
* All shps should be in (same) model coordinate units (e.g. ft.)
* If model domain contains a major divide, need to merge relevant NHD datasets (e.g. 04 and 07) prior to running this script


### Workflow for building SFR input:

##### 0) run computeZonal.py
	
     Inputs:
     Polygon shapefile export of model grid (e.g. as exported by Groundwater Vistas)
     A DEM for model area, doesn't need to be clipped
     
##### 1) run SFR_preproc.py
(see Steps 1-9 in Howard Reeves' SFR notes)

     Inputs: 
     Shapefile of model grid cells with model top elevations (from computeZonal.py above)
     Shapefile polygon of model domain (merged polygon of gridcells)
     PlusflowVAA database from NHDPlus v2
     Elevslope database from NHDPlus v2
     original NHDFlowline shapefile from NHDPlus v2
     NHDFlowline shapefile clipped to model grid * this can be taken out if code is added to perform the clipping
     
     Outputs:
     river_explode.shp (lines)
	 river_cells.shp
     river_cells_dissolve.shp (grid cells (SFR reaches) that contain NHD segment elevations and channel location information (start and end xy; reach length in cell))



##### 2) run intersect.py

     Inputs:
          - <casename>.in
          - <Flowlines_unclipped>.shp
          - <MFdomain>.shp
          - <Flowlines>.shp
          
     Outputs:
          - <NHD>.shp
          - boundaryClipsRouting.txt
          - boundary_manual_fix_issues.txt (if manual intervention required)

Look at boundary\_manual\_fix\_issues.txt output file and find COMIDs that leave and reenter the domain. Fix them in Flowlines (from input file) and in river\_explode.py (if necessary).

Once runs through with no manual fix messages, move on.

##### 2a) run CleanupRiverCells.py,
which trims down the river cells shapefiles to reflect the deletiosn that were made in the previous step deleting COMIDs from river_explode.shp
This also makes a backup copy of fix\_comids.txt --> fix\_comids\_backup.txt which can be used to inspect the results of rerunning AssignRiverElev.py

##### 2b) rerun intersect.py
Only if there were COMID repair issues on the previous run. 

##### 3) run Assign_and_Route.py
This works through the COMID definitions along with model cell boundaries to both assign COMIDs to model cells and assign elevations based on those COMID segments. Also handles the routing from NHDplus as much as possible.

Once this runs with no manual fixes, move to step 4. Otherwise handle COMIDs listed in fix\_comids.txt

##### 3a) run CleanupRiverCells.py
If fix\_comids.txt indicated necessary manual fixes, follow direction in that file then run CleanupRiverCells.py

##### 3b) rerun Assign_and_Route.py


##### 4) run Finalize\_SFR\_revised\_OCT.py
Creates the actual SFR file with elevations and routing
     
##### 5) run Fix\_w\_DEM.py
	- runs Fix\_segment\_ends.py to adjust segment end elevations to match model TOP as closely as possible (minus an adjustable incising parameter)
	- checks for backward routing in segment ends
	- adjust segment interior reaches to model TOP, as long as they are monotonically downhill and don't go below segment end
	- uses linear interpolation otherwise
	- generates PDF files showing comparative statistics for SFR vs. model top before and after
	- generates PDF files with comparative plots of SFR before and after, with model top, also 50 worst floating and incised segments
	- can also generate a plot of segments that go below the model bottom, but need to run Assign_Layers.py first to get below_bot.csv
	
	
	Dependencies:
		-STOP_compare.py (generates columns of model top elev, SFR streambed top, and their differece, indexed by segment)
		-Fix_segment_ends.py  (Adjusts segment end elevations within the constraints of up/down segment min/max elevations)
		-Plot_segment_profiles.py  (Has generic code to plot SFR segments with land surface or any other elevation)
		-pandas (non-standard module of Python- used to quickly sort 50 biggest floating and incised reaches; this could probably be done by numpy)
		
	Inputs:
		- GWVmat1 (from SFR_utilities.py)
		- GWVmat2 (from SFR_utilties.py)
		- ascii array of model top elevations
		- ascii multi-layer array of model bottom elevations
	
	Outputs:
		- GWVmat1
		- GWVmat2
		- selected_segments.pdf (profile plots of selected segments compared to land surface)
		- PDF files showing profiles of most floating and most incised segments
		- PDF files showing cumulative distribution of landsurface/streambed elevation differences and summary statistics before and after
		- fix_routing_report='Fix_routing_report.txt' # records which segments were given new outsegs by Fix_routing
		- fix_ends_report='Fixed_segment_ends.csv' # records changes made to segment end elevations
		- fix_ends_errors='Fix_segment_ends_errors.txt' # records instances where segments were left with backwards routing
		- end_interp_report='fix_w_DEM_interps.txt' # file recording adjustments made using up/dn_increments
		- error_report='fix_w_DEM_errors.txt' # file for reporting 0 slope errors, etc.
		- STOP_comp_SFR_utilities='STOP_compare_SFR_utilities.csv'# from STOP_compare.py
		- STOP_comp_fixwDEM='STOP_compare_fix_w_DEM.csv' # from STOP_compare.py

	###### NOTE: there are also a number of hard-coded settings for fix_w_DEM that are listed with descriptions under the #Settings comment in the source code
		- in particular,
			- when Fix_ends is True, functions in Fix_segment_ends.py will be run, which try to reduce floating and incision of segment end elevations relative to the land surface
			- when Fix_routing is True, functions in Fix_segment_ends.py will search within a specified radius for nearby SFR cells, and reroute to the nearby cell with the lowest starting elevation. This can help with problems at confluences, where two streams come together and route into eachother before routing downstream.
			- These two features are run before the main body of fix_w_DEM, which attemps to correct interior reaches within the segment
			- It's not clear if Fix_ends and Fix_routing are always worth using. Best to run without, and then look over the results (saving the summary PDFs to different file names). And then run with these features turned on to see if there is any improvement.
		- At the start, Fix_w_DEM saves original copies of Mat1 and Mat2, with "_old" appended to the name. If rerunning, need to restore the old versions first (won't save output otherwise)
		
##### 6) run Assign_Layers.py
	Using the stream bottom elevations, assigns layers to all of the SFR cells
	If for some reason stream bottoms go below the model bottom, will post a warning, and write all of the violations to output.
	If lowering the model bottom correct the violations is reasonable, then can rerun the program with Lowerbot=True; this will adjust the model bottom elevations to accomodate the SFR cells
	Also, there may have been some segments after fix_w_DEM that weren't handled properly. In the Columbia Co. model, many of these were headwaters that were in NHD incorrectly (i.e., connecting across a saddle to the wrong valley, etc.)
	If these segments are deleted, the segment numbering is no long continuous, and MODFLOW won't run.
	Assign_Layers will adjust the segment numbering accordingly, to maintain continuity.

	Inputs:
		- rows x columns text matrix of bottom elevations for all layers (either exported from GW Vistas, or generated seperate from DIS file)
		- rows x columns text matrix of model top elevations
		- GWVmat1 from above
		- GWVmat2 from above
	Outputs:
		- MODFLOW SFR package file
		- revised GWVmat1 file with layers assigned
		- revised GWVmat2 file with updated segment numbers (if segment renumbering was necessary)
		- 'SFR_layer_assignments.txt' (summarizes number of reaches assigned to each layer)
		- 'Model_bottom_adjustments.pdf' PDF images of model bottom before and after adjustments
		- 'below_bot.csv' - records all stream bottoms that go below the model bottom, with corresponding layer elevation info


	


