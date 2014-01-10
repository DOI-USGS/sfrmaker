__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc

infile = 'SFR_input.xml'

SFRdata = SFRc.SFRInput(infile)

SFRpre = SFRc.SFRpreproc(SFRdata)

SFRpre.clip_and_join_attributes()


SFRops = SFRc.SFROperations(SFRdata)

SFRops.intersect()


FIDdata = SFRc.FIDPropsAll()

FIDdata.populate(SFRdata)

FIDdata.return_fid_comid_list()


'''
FIDdata.populate_elevations(SFRdata)
'''
FIDdata.populate_routing(SFRdata)

i = 2