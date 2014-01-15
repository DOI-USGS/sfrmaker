__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc

infile = 'SFR_input.xml'

SFRdata = SFRc.SFRInput(infile)

if SFRdata.preproc:
    print 'Running preprocessing routine'
    SFRpre = SFRc.SFRpreproc(SFRdata)

    SFRpre.clip_and_join_attributes()


SFRops = SFRc.SFROperations(SFRdata)

SFRops.intersect()

FIDdata = SFRc.FIDPropsAll()

FIDdata.populate(SFRdata)

FIDdata.return_fid_comid_list()

SFRops.make_rivers_table(FIDdata)

FIDdata.populate_elevations(SFRdata)

COMIDdata = SFRc.COMIDPropsAll()

LevelPathdata = SFRc.LevelPathIDpropsAll()

COMIDdata.populate_routing(SFRdata, FIDdata, LevelPathdata)

LevelPathdata.return_cutoffs(FIDdata)

i = 2