__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc
import sfr_plots

infile = 'SFR_input.xml'

SFRdata = SFRc.SFRInput(infile)

SFRops = SFRc.SFROperations(SFRdata)

SFRpre = SFRc.SFRpreproc(SFRdata)


if SFRdata.preproc:
    print 'Running preprocessing routine'

    SFRpre.clip_and_join_attributes(SFRops)


#SFRops.assign_layers(SFRdata)

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

SFRops.reach_ordering(COMIDdata, FragIDdata, LevelPathdata)
'''
SFRpre.intersect_contours(SFRdata) # this needs to be in the example Main!
ContourElevs = SFRc.ElevsFromContours(SFRdata)
ContourElevs.get_contour_intersections(FragIDdata, COMIDdata)
ContourElevs.assign_elevations_to_FragID(FragIDdata, COMIDdata)
'''
SFRpre.intersect_DEM(SFRdata) # this needs to be in the example Main!
DEMelevs = SFRc.ElevsFromDEM()
DEMelevs.DEM_elevs_by_FragID(SFRdata, SFRops)
DEMelevs.connect_downhill(FragIDdata)

SFRp = sfr_plots.plot_elevation_profiles(SFRdata, COMIDdata)
SFRp.read_DIS()
SFRp.get_comid_plotting_info(FragIDdata)
SFRp.plot_profiles('elevs_from_contours.pdf')
'''
saveme ={'COMIDdata' : COMIDdata,
         'FragIDdata' : FragIDdata,
         'LevelPathdata' : LevelPathdata,
         'CELLdata' : CELLdata}
SFRc.savetmp(saveme)


loadme = ['COMIDdata', 'FragIDdata', 'LevelPathdata']
instuff = SFRc.loadtmp(loadme)

LevelPathdata = instuff['LevelPathdata']
FragIDdata = instuff['FragIDdata']
COMIDdata = instuff['COMIDdata']
'''



Segmentdata = SFRc.SFRSegmentsAll()
Segmentdata.divide_at_confluences(LevelPathdata, FragIDdata, COMIDdata, CELLdata, SFRdata)
Segmentdata.accumulate_same_levelpathID(LevelPathdata, COMIDdata, FragIDdata, SFRdata, CELLdata)

#make some output

SFRoutput = SFRc.SFRoutput(SFRdata)
SFRoutput.write_SFR_tables(Segmentdata)
SFRoutput.build_SFR_package()
#SFRoutput.build_SFR_shapefile(Segmentdata)



i = 2