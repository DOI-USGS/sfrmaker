__author__ = 'mnfienen'

import SFR_classes as SFRc

infile = 'SFR_setup.xml'

SFRdata = SFRc.SFRInput(infile)

COMIDdata = SFRc.COMIDPropsAll()

COMIDdata.populate(SFRdata)

COMIDdata.populate_elevations(SFRdata)

COMIDdata.populate_routing(SFRdata)


i = 2