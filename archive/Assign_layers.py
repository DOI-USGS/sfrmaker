# script for running Assign Layers outside of SFR_main
__author__ = 'Fienen, Reeves, Leaf - USGS'

import sys
sys.path.append('D:/ATLData/Documents/GitHub/SFR')
import SFR_classes as SFRc
import sfr_plots

infile = 'SFR_input_BR.xml'

SFRdata = SFRc.SFRInput(infile)
'''
SFRp = sfr_plots.plot_elevation_profiles(SFRdata)
SFRp.plot_Mat1_profiles('Segment_profiles_from_Mat1.pdf')
'''
SFRops = SFRc.SFROperations(SFRdata)

SFRops.assign_layers(SFRdata)


SFRoutput = SFRc.SFRoutput(SFRdata)
SFRoutput.build_SFR_package()
#SFRoutput.buildSFRshapefile2()
