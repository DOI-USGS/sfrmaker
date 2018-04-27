import sys
import glob
import os
sys.path.append('..')
sys.path += glob.glob('/Users/aleaf/Documents/GitHub/*')
import numpy as np
import flopy
fm = flopy.modflow
from sfrmaker import lines, grid, sfrdata

# NHDPlus input files (see the input requirements in the SFRmaker readme file
# (Note that multiple datasets can be supplied as lists;
# when the SFR area covers multiple drainage basins)
NHDPlus_paths = ['/Users/aleaf/Documents/NHDPlus/NHDPlusMS/NHDPlus08/']

# output folder
outdir = 'temp/'
if not os.path.isdir(outdir):
    os.mkdir(outdir)

lns = lines.from_NHDPlus_v2(NHDPlus_paths=NHDPlus_paths,
                            filter='quiver_outlet_area.shp')
nrow, ncol = 112, 160
xll, yll = 497812, 1160724
dxy = 1000 /.3048
sr = flopy.utils.SpatialReference(delr=np.ones(ncol)*dxy,
                                  delc=np.ones(nrow)*dxy,
                                  lenuni=1,
                                  xll=xll, yll=yll, rotation=0,
                                  proj4_str='+init=epsg:5070')

sfr = lns.to_sfr(sr=sr, active_area='quiver_outlet_area.shp')

sfr.reach_data['strtop'] = sfr.interpolate_to_reaches('elevup', 'elevdn')
sfr.get_slopes()

sfr.write_package(outdir + 'quiver_outlet.sfr')
sfr.export_sfrlines(outdir + 'quiver_outlet_sfrlines.shp')
sfr.export_routing(outdir + 'quiver_outlet_routing.shp')
j=2