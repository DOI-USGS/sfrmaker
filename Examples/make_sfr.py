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
pfvaa_files = ['data/PlusFlowlineVAA.dbf']
plusflow_files = ['data/PlusFlow.dbf']
elevslope_files = ['data/elevslope.dbf']
flowlines = ['data/NHDFlowlines.shp']

# output folder
outdir = 'temp/'
if not os.path.isdir(outdir):
    os.mkdir(outdir)

lns = lines.from_NHDPlus_v2(NHDFlowlines=flowlines,
                            PlusFlowlineVAA=pfvaa_files,
                            PlusFlow=plusflow_files,
                            elevslope=elevslope_files,
                            filter='data/grid.shp')

nrow, ncol = 112, 160
dxy = 250
sr = flopy.utils.SpatialReference(delr=np.ones(ncol)*dxy,
                                  delc=np.ones(nrow)*dxy,
                                  lenuni=1,
                                  xll=682688, yll=5139052, rotation=0,
                                  proj4_str='+init=epsg:26715')

m = fm.Modflow('example', model_ws=outdir)
dis = fm.ModflowDis(m, nlay=1, nrow=nrow, ncol=nrow,
                    top=1000, botm=0)
m.sr = sr

grd = grid.from_sr(sr, active_area='data/active_area.shp')
#grd = grid.from_shapefile('data/grid.shp', icol='row', jcol='column')

sfr = lns.to_sfr(grd)

sfr.reach_data['strtop'] = sfr.interpolate_to_reaches('elevup', 'elevdn')
sfr.get_slopes()

sfr.write_package(outdir + 'example.sfr')
sfr.export_sfrlines(outdir + 'example.shp')
j=2