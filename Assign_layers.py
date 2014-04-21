# script for running Assign Layers outside of SFR_main
__author__ = 'Fienen, Reeves, Leaf - USGS'

import SFR_classes as SFRc
#import sfr_plots

infile = 'SFR_input_RPU.xml'

SFRdata = SFRc.SFRInput(infile)


SFRops = SFRc.SFROperations(SFRdata)

SFRops.assign_layers(SFRdata)

SFRoutput = SFRc.SFRoutput(SFRdata)
SFRoutput.build_SFR_package()
