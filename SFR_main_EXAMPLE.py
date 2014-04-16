__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc
import sfr_plots

infile = 'Example.xml'

SFRdata = SFRc.SFRInput(infile)

SFRpre = SFRc.SFRpreproc(SFRdata)

SFRops = SFRc.SFROperations(SFRdata)

if SFRdata.preproc:
    print 'Running preprocessing routine'

    SFRpre.clip_and_join_attributes(SFRops)

SFRops.intersect()

FragIDdata = SFRc.FragIDPropsAll()

FragIDdata.populate(SFRdata)

FragIDdata.return_FragID_comid_list()

SFRops.make_rivers_table(FragIDdata)

FragIDdata.populate_elevations(SFRdata)

COMIDdata = SFRc.COMIDPropsAll()

CELLdata = SFRc.CellPropsAll()

CELLdata.populate_cells(SFRdata)

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

# Comparison plots of streambed elevations (by COMID) for different elevation methods
SFRp = sfr_plots.plot_elevation_profiles(SFRdata)
SFRp.read_DIS()
SFRp.get_comid_plotting_info(FragIDdata, COMIDdata, SFRdata)
SFRp.plot_profiles('Elevation_method_comparison.pdf')

# COMMENT
SFRops.reach_ordering(COMIDdata,
                      FragIDdata,
                      LevelPathdata)

Segmentdata = SFRc.SFRSegmentsAll()

Segmentdata.divide_at_confluences(LevelPathdata, FragIDdata,
                                  COMIDdata, CELLdata, SFRdata)
Segmentdata.accumulate_same_levelpathID(LevelPathdata, COMIDdata,
                                        FragIDdata, SFRdata, CELLdata)

# Plot SFR segment profiles with final elevations
SFRp.get_segment_plotting_info(Segmentdata)
SFRp.plot_profiles('Segment_profiles.pdf')

# Write Ouput
SFRoutput = SFRc.SFRoutput(SFRdata)
SFRoutput.write_SFR_tables(Segmentdata)
SFRops.assign_layers(SFRdata)
SFRoutput.build_SFR_package()
#SFRoutput.build_SFR_shapefile(Segmentdata)
print "Done"

i = 2