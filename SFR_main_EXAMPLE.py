__author__ = 'Fienen, Reeves, Leaf - USGS'

import sys
sys.path.append('D:\\ATLData\\Documents\\GitHub\\SFR\\') # path to SFRmaker files
import SFR_classes as SFRc

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


# COMMENT
SFRops.reach_ordering(COMIDdata,
                      FragIDdata,
                      LevelPathdata)

Segmentdata = SFRc.SFRSegmentsAll()

Segmentdata.divide_at_confluences(LevelPathdata, FragIDdata,
                                  COMIDdata, CELLdata, SFRdata)
Segmentdata.accumulate_same_levelpathID(LevelPathdata, COMIDdata,
                                        FragIDdata, SFRdata, CELLdata)


# Write Ouput
SFRoutput = SFRc.SFRoutput(SFRdata)
SFRoutput.write_SFR_tables(Segmentdata)