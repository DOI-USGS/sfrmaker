__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc

infile = 'SFR_input.xml'

SFRdata = SFRc.SFRInput(infile)
'''
SFRpre = SFRc.SFRpreproc(SFRdata)

SFRpre.clip_and_join_attributes()
'''

SFRops = SFRc.SFROperations(SFRdata)

SFRops.intersect()

COMIDdata = SFRc.COMIDPropsAll()

COMIDdata.populate(SFRdata)



COMIDdata.populate_elevations(SFRdata)

COMIDdata.populate_routing(SFRdata)


i = 2