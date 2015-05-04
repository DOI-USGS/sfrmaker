# adds streamflow output from SFR package to river_explode linework
import sfr_plots as SFRp
import shutil
import os

# copy SFR output file to folder with post-processed GIS results
shutil.copy(<...streamflow.dat>, os.path.join(<path>, <...streamflow.dat>))

# plot_streamflows arguments are (<MODFLOW DIS file>,<intersect> from SFR XML input file (or similar, <SFR package output file>)
stuff = SFRp.plot_streamflows(<MFdis>, <intersect>, <...streamflow.dat>, <node_attribute> )
stuff.join_SFR_out2streams()

# can use symbology in SFR_flow_symbology.lyr to visualize flow in in ArcMap as line thickness