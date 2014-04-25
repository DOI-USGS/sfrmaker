__author__ = 'Fienen, Reeves, Leaf - USGS'

import os
import SFR_classes as SFRc
import sfr_plots

# Specify input XML file
infile = 'Wbasin.xml'

# Instantiate class for storing input data
SFRdata = SFRc.SFRInput(infile)

# Instantiate class for GIS pre-processing of inputs
SFRpre = SFRc.SFRpreproc(SFRdata)

# Instantiate class for performing various SFR-building operations
SFRops = SFRc.SFROperations(SFRdata)

# Run preprocessing routine if it hasn't been run yet
if SFRdata.preproc:
    print 'Running preprocessing routine'
    SFRpre.clip_and_join_attributes(SFRops)

# Calculate length of stream fragments in each grid cell
# Identify COMIDs that cross the SFR area boundary
# Identify COMIDs that appear to be isolated (not routed and not near other COMIDs)
# Identify COMIDs with multiple starts and ends (usually these are streams with multiple meanders across the boundary
SFRops.intersect()

# Instantiate class for storing information on stream "fragments"
FragIDdata = SFRc.FragIDPropsAll()

# Populate the "fragment" class with information stored in SFRdata.intersect (created by the preprocessing routine)
FragIDdata.populate(SFRdata)

# Cross reference fragments with their respective COMIDs
FragIDdata.return_FragID_comid_list()

# Build a table of NHD streambed elevation information; also add this information to fragment properties
SFRops.make_rivers_table(FragIDdata)
FragIDdata.populate_elevations(SFRdata)

# Instantiate class for storing information on COMIDs
COMIDdata = SFRc.COMIDPropsAll()

# Prune cells with only short lengths of stream (debris left over from intersection of NHD linework with grid)
CELLdata = SFRc.CellPropsAll()
CELLdata.populate_cells(SFRdata)

# Instantiate class for storing NHD "LevelPath" information
# use LevelPath informaiton to route COMIDs
LevelPathdata = SFRc.LevelPathIDpropsAll()
COMIDdata.populate_routing(SFRdata, FragIDdata, LevelPathdata, CELLdata)
COMIDdata.return_hydrosequence_comid()
FragIDdata.return_cellnum_LevelPathID(LevelPathdata)
LevelPathdata.return_cutoffs(FragIDdata, CELLdata, SFRdata)

# Get streambed elevation information from topographic contours
if SFRdata.elevflag == 'elevation_contours' or SFRdata.elev_comp:
    SFRpre.intersect_contours(SFRdata) # this needs to be in the example Main!
    ContourElevs = SFRc.ElevsFromContours(SFRdata)
    ContourElevs.get_contour_intersections(FragIDdata, COMIDdata)
    ContourElevs.assign_elevations_to_FragID(FragIDdata, COMIDdata)

# Get streambed elevation information from DEM
if SFRdata.elevflag == 'smoothed_DEM' or SFRdata.elev_comp:
    SFRpre.intersect_DEM(SFRdata) # this needs to be in the example Main!
    DEMelevs = SFRc.ElevsFromDEM()
    DEMelevs.DEM_elevs_by_FragID(SFRdata, SFRops)
    DEMelevs.connect_downhill(FragIDdata)

# Initialize plotting class
SFRp = sfr_plots.plot_elevation_profiles(SFRdata)
SFRp.read_DIS()

# Comparison plots of streambed elevations (by COMID) for different elevation methods
if SFRdata.elev_comp:
    SFRp.get_comid_plotting_info(FragIDdata, COMIDdata, SFRdata)
    SFRp.plot_profiles(os.path.join(SFRdata.working_dir, 'Elevation_method_comparison.pdf'))


saveme = {'COMIDdata': COMIDdata,
         'FragIDdata': FragIDdata,
         'LevelPathdata': LevelPathdata,
         'CELLdata': CELLdata,
         'SFRdata':  SFRdata}
SFRc.savetmp(saveme)

#a = SFRc.loadtmp['SFRdata', 'COMIDdata', 'CELLdata', 'FragIDdata', 'LevelPathdata']

#COMIDdata = a['COMIDdata']
#FragIDdata = a['FragIDdata']
#LevelPathdata = a['LevelPathdata']
#CELLdata = a['CELLdata']
#SFRdata = a['SFRdata']


# Build SFR package segments
SFRops.reach_ordering(COMIDdata, FragIDdata, LevelPathdata)
Segmentdata = SFRc.SFRSegmentsAll()
Segmentdata.divide_at_confluences(LevelPathdata, FragIDdata, COMIDdata, CELLdata, SFRdata)
Segmentdata.accumulate_same_levelpathID(LevelPathdata, COMIDdata, FragIDdata, SFRdata, CELLdata)


# Plot SFR segment profiles with final elevations, output csv file comparing elevations to land surface
#SFRp.get_segment_plotting_info(Segmentdata)
#SFRp.plot_profiles(os.path.join(SFRdata.working_dir, 'Segment_profiles.pdf'))
SFRp.check4elevation_issues(FragIDdata, COMIDdata, SFRdata, Segmentdata)

# Write Output:
SFRoutput = SFRc.SFRoutput(SFRdata)

# Write the "MAT1" and "MAT2" tables, which contain the reach and segment level information needed for SFR input
SFRoutput.write_SFR_tables(Segmentdata)

# Assign model layers to SFR cells based on streambed elevations and cell bottoms in DIS file
SFRops.assign_layers(SFRdata)

# Write the SFR package input (.sfr) file
SFRoutput.build_SFR_package()

# Make a shapefile of SFR cells with attribute information for visualisation in GIS
#SFRoutput.build_SFR_shapefile(Segmentdata)
SFRoutput.buildSFRshapefile2()
print "\n\nDone!"

i = 2