__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc

infile = 'SFR_setup.xml'

SFRdata = SFRc.SFRInput(infile)
"""
SFRops = SFRc.SFROperations(SFRdata)

SFRops.intersect()
"""
COMIDdata = SFRc.COMIDPropsAll()

COMIDdata.populate(SFRdata)

COMIDdata.populate_elevations(SFRdata)

COMIDdata.populate_routing(SFRdata)


i = 2