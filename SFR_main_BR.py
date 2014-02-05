__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc
import sfr_plots

infile = 'SFR_input_BR.xml'

SFRdata = SFRc.SFRInput(infile)

SFRpre = SFRc.SFRpreproc(SFRdata)

SFRops = SFRc.SFROperations(SFRdata)

if SFRdata.preproc:
    print 'Running preprocessing routine'

    #SFRpre.clip_and_join_attributes(SFRops)



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



'''
saveme ={'SFRdata':  SFRdata, 'COMIDdata': COMIDdata, 'FragIDdata': FragIDdata,
              'SFRops': SFRops, 'LevelPathdata' : LevelPathdata, 'CELLdata' : CELLdata}

SFRc.savetmp(saveme)

a = SFRc.loadtmp(['SFRdata', 'CELLdata', 'COMIDdata', 'FragIDdata', 'SFRops', 'LevelPathdata'])

SFRdata = a['SFRdata']
COMIDdata = a['COMIDdata']
FragIDdata = a['FragIDdata']
SFRops = a['SFRops']
LevelPathdata = a['LevelPathdata']
CELLdata = a['CELLdata']
'''
SFRpre.intersect_contours(SFRdata) # this needs to be in the example Main!
ContourElevs = SFRc.ElevsFromContours(SFRdata)
ContourElevs.get_contour_intersections(FragIDdata, COMIDdata)
ContourElevs.assign_elevations_to_FragID(FragIDdata, COMIDdata)

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
         'CELLdata' : CELLdata, 'SFRdata':  SFRdata,}
SFRc.savetmp(saveme)


a = SFRc.loadtmp['SFRdata', 'COMIDdata', 'CELLdata', 'FragIDdata', 'LevelPathdata']

COMIDdata = a['COMIDdata']
FragIDdata = a['FragIDdata']
LevelPathdata = a['LevelPathdata']
CELLdata = a['CELLdata']
SFRdata = a['SFRdata']
'''
SFRops.reach_ordering(COMIDdata,
                      FragIDdata,
                      LevelPathdata)

Segmentdata = SFRc.SFRSegmentsAll()

Segmentdata.divide_at_confluences(LevelPathdata, FragIDdata,
                                  COMIDdata, CELLdata, SFRdata)
Segmentdata.accumulate_same_levelpathID(LevelPathdata, COMIDdata,
                                        FragIDdata, SFRdata, CELLdata)

#make some output

SFRoutput = SFRc.SFRoutput(SFRdata)
SFRoutput.write_SFR_tables(Segmentdata)
SFRoutput.build_SFR_package()
SFRoutput.build_SFR_shapefile(Segmentdata)


i = 2