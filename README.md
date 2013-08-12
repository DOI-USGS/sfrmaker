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

* currently the SFR package file is create prior to final elevation corrections with fix_w_dem.py. Until this is fixed, SFR package input should be manually created from the Groundwater Vistas input tables.
Note:
# All shps should be in (same) model coordinate units (e.g. ft.)
If model domain contains a major divide, need to merge relevant NHD datasets (e.g. 04 and 07) prior to running this script


Workflow for building SFR input:

1) run SFR_preproc.py (see Steps 1-9 in Howard Reeves' SFR notes)

     Inputs: 
     Polygon shapefile export of model grid (e.g. as exported by Groundwater Vistas)
     Shapefile polygon of model domain (merged polygon of gridcells)
     A DEM for model area, doesn't need to be clipped
     PlusflowVAA database from NHDPlus v2
     Elevslope database from NHDPlus v2
     original NHDFlowline shapefile from NHDPlus v2
     NHDFlowline shapefile clipped to model grid * this can be taken out if code is added to perform the clipping
     
     Outputs:
     river_explode.shp (lines)
     river_cells_dissolve.shp (grid cells (SFR reaches) that contain NHD segment elevations and channel location information (start and end xy; reach length in cell))



2) run AssignRiverElev.py, which produces a table of reach midpoint elevations via linear interpolation from NHD segment endpoint elevations

     Inputs: river_explode.shp
     Outputs: river_elevs.dbf
                    fix_comids.txt      # list of comids with two or more sets of endpoints.  These are usually segments that meander out of and then back into the grid.

3) manually delete flowlines and corresponding gridcell polygons for segments that have multiple starts/ends (e.g. those that meander out of the grid and back in)

4) run intersect.py
     Inputs:
          - flowlines clipped from previous step
          - unclipped (original) flowlines from NHDPlus
          - polygon of grid boundary

     Outputs:
          - boundaryclipsrouting.txt (for segments intersecting model boundary)
          - 'NHD_intersect_edited.shp'


5) run RouteStreamNetwork.py

      Inputs:
          - boundaryclipsrouting.txt
          - river_w_elevations.shp    # "ELEV" in input file; contains line segments for each SFR reach, with preliminary elevation, length and width information
          - PlusFlow database from NHDPlus v2
          
      Outputs:
          - check_network.txt     List of to/from routing for segments
          - flagged_comids.txt    List instances where downstream segment has higher start elevation

     Check Network lines are created as follows:
          - list of COMIDs is obtained from river_w_elevations
          - FROMCOMID field in PlusFlow is queried for each comid
          - If the corresponding TOCOMID is in the list of COMIDs:
               - FROMCOMID, TOCOMID is written
          - Else:
               - FROMCOMID, 99999 is written (this will happen if stream flows outside of grid)
          - If there are no corresponding FROMCOMIDs:
               - TOCOMID field is queried for comid in list
          - If the corresponding FROMCOMID is in the list of comids:
               - FROMCOMID, TOCOMID is written
          - Else:
               - 99999, TOCOMID is written (this will happen if stream flows into grid)

5.1) (if necessary) run Fix_flagged_comids.py (this program still needs some work):

      Handles the segments in flagged_comids.txt (e.g., making elevations consistent with routing) by:
          - deleting segments that have an FTYPE of "ArtificialPath" and a "Divergence" value of 2. These segments appeared mostly to be small parts of braids in the channel. References to these segments are removed from check_network.txt
          - for all other segments in flagged_comids.txt:
               - identifies segment end-reaches, and corresponding cells
               - looks for other segments in those cells, and identifies global max/min elevations (these could be used to edit routing)
               - takes segment max elevation, and max elevation of next downstream segment (from check_network.txt)
                    - re-calculates elevations of all reaches in segment by linear interpolation
                    - replaces original elevations in river_w_elevations.txt
                    - 0 gradients on headwater segments are corrected by increasing first reach elevation incrementally
                    - other zero gradients are passed to next stages of SFR setup process.
     
     Note: the part of this script that identifies headwaters is incorrect, also, for some reason not all of the ArtificialPaths from above were removed from check_network.txt; the logic in this section of Fix_flagged likely needs to be fixed. Also the resetting of elevations in river_w_elevations could be greatly streamlined by instead using the fixsegelevs function (see below)
     18 segments were deleted, 22 were corrected
     
     fix_segment_elevs.py has a function that enforces downhill monotonicity within a segment, if first and last cellnum (node number), and a desired minimum and maximum elevation are known
     - easy to execute manually in the Arcpy window (while looking at the flagged segments in ArcMap)
     - would be better if integrated into a correct version of Fix_flagged_comids.py
          
6) run RouteRiverCells.py
     Inputs:
          - 'river_cells_dissolve.shp' # from SFR_preproc.py
          - 'river_w_elevations.shp' # uses this to route reaches within segments
          - 'NHD_intersect_edited.shp'  # from intersect.py
          - PlusFlowlineVAA database from NHDPlus
          - check_network.txt  # from RouteStreamNetwork.py
          - boundaryclipsrouting.txt  # from intersect.py

     Outputs:
          - reach_ordering.txt
          - routed_cells.txt

7) run Finalize_SFR.py
added manual step of joining river_w_elevations to river_cells to this script
also removed hard-coding for grid information (obtains from shapefile of grid exported from GWV) 
Note: getting grid info this way is slow, should be done ONCE at beginning, in SFR_preproc. Could be appended to input file for rapid re-referencing in following scripts.

had to fix one segment using fixsegelevs function (stored in fix_segment_elevs.py). This function, given the COMID, cellnums of first and last reaches, and desired min/max elevations, resets elevations of all reaches using linear interpolation. Fix_flagged_comids.py (and maybe other scripts) could be cleaned up and made more robust by importing this function, instead of relying on built-in code


8) SFR_utilities.py, to improve any segment start/end elevations where there is room for improvement.
