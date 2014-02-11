# adds streamflow output from SFR package to river_explode linework
import sfr_plots as SFRp

# plot_streamflows arguments are (<MODFLOW DIS file>,<intersect> from SFR XML input file (or similar, <SFR package output file>)
stuff = SFRp.plot_streamflows('BadRiver.dis', 'river_explode.shp', 'BadRiver_streamflow.dat')
stuff.join_SFR_out2streams()
# can use symbology in SFR_flow_symbology.lyr to visualize flow in in ArcMap as line thickness