__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc
import sfr_plots

infile = 'SFR_input.xml'

SFRdata = SFRc.SFRInput(infile)

'''
if SFRdata.preproc:
    print 'Running preprocessing routine'
    SFRpre = SFRc.SFRpreproc(SFRdata)

    #SFRpre.clip_and_join_attributes()
'''
SFRops = SFRc.SFROperations(SFRdata)

SFRops.assign_layers(SFRdata)

SFRops.intersect()
'''
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

LevelPathdata.return_cutoffs(FragIDdata, CELLdata, SFRdata)

SFRops.reach_ordering(COMIDdata, FragIDdata, LevelPathdata)

SFRpre.intersect_contours(SFRdata)
ContourElevs = SFRc.ElevsFromContours(SFRdata)
ContourElevs.get_contour_intersections(FragIDdata, COMIDdata)

SFRpre.intersect_DEM(SFRdata)
DEMelevs = SFRc.ElevsFromDEM()
DEMelevs.DEM_elevs_by_FragID(SFRdata, SFRops)


saveme ={'SFRdata':  SFRdata, 'COMIDdata': COMIDdata, 'FragIDdata': FragIDdata,
              'SFRops': SFRops, 'ContourElevs': ContourElevs, 'DEMelevs': DEMelevs, 'LevelPathdata' : LevelPathdata}

SFRc.savetmp(saveme)
'''
a = SFRc.loadtmp(['SFRdata', 'COMIDdata', 'FIDdata', 'SFRops', 'ContourElevs', 'LevelPathdata'])

SFRdata = a['SFRdata']
COMIDdata = a['COMIDdata']
FragIDdata = a['FIDdata']
SFRops = a['SFRops']
Contour_elevs = a['ContourElevs']
LevelPathdata = a['LevelPathdata']

DEMelevs = SFRc.ElevsFromDEM()
DEMelevs.DEM_elevs_by_FragID(SFRdata, SFRops)

Contour_elevs.assign_elevations_to_FragID(FragIDdata, COMIDdata)
DEMelevs.connect_downhill(FragIDdata)

SFRp = sfr_plots.plot_segments(SFRdata, COMIDdata)
SFRp.read_DIS()
SFRp.get_segment_plotting_info(FragIDdata)
SFRp.plot_profiles('elevs_from_contours.pdf')

saveme ={'COMIDdata' : COMIDdata,
         'FragIDdata' : FragIDdata,
         'LevelPathdata' : LevelPathdata}
SFRc.savetmp(saveme)

'''
loadme = ['COMIDdata', 'FragIDdata', 'LevelPathdata']

instuff = SFRc.loadtmp(loadme)


SFRops.reach_ordering(instuff['COMIDdata'],
                      instuff['FragIDdata'],
                      instuff['LevelPathdata'])


Segmentdata = SFRc.SFRSegmentsAll()
Segmentdata.divide_at_confluences(LevelPathdata, FragIDdata, COMIDdata, CELLdata)
Segmentdata.accumulate_same_levelpathID(LevelPathdata, COMIDdata, FragIDdata, SFRdata)


saveme ={}


SFRc.savetmp(saveme)

a = SFRc.loadtmp(saveme)

SFRdata = a['SFRdata']
'''

i = 2