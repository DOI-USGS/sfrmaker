__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc

infile = 'SFR_input.xml'

SFRdata = SFRc.SFRInput(infile)


if SFRdata.preproc:
    print 'Running preprocessing routine'
    SFRpre = SFRc.SFRpreproc(SFRdata)

    SFRpre.clip_and_join_attributes()

SFRops = SFRc.SFROperations(SFRdata)

SFRops.assign_layers(SFRdata)

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

LevelPathdata.return_cutoffs(FragIDdata, CELLdata, SFRdata)

SFRops.reach_ordering(COMIDdata, FragIDdata, LevelPathdata)

'''
saveme ={'COMIDdata' : COMIDdata,
         'FragIDdata' : FragIDdata,
         'LevelPathdata' : LevelPathdata}
SFRc.savetmp(saveme)


loadme = ['COMIDdata', 'FragIDdata', 'LevelPathdata']

instuff = SFRc.loadtmp(loadme)


SFRops.reach_ordering(instuff['COMIDdata'],
                      instuff['FragIDdata'],
                      instuff['LevelPathdata'])

'''
Segmentdata = SFRc.SFRSegmentsAll()
Segmentdata.divide_at_confluences(LevelPathdata, FragIDdata, COMIDdata, CELLdata)

'''
saveme ={}


SFRc.savetmp(saveme)

a = SFRc.loadtmp(saveme)

SFRdata = a['SFRdata']
'''

i = 2