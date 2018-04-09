"""Example workflow for generating an SFR package
from NHDPlus data, a DEM, and MODFLOW grid.
"""
import sys
import os
import glob
sys.path += glob.glob('/Users/aleaf/Documents/GitHub/*')
# if preproc.py and postproc.py are in a different folder from this script,
# add that folder to the path so python can find it
sys.path.append('../')
import time
import numpy as np
import flopy
from preproc import NHDdata
from postproc import SFRdata
fm = flopy.modflow

# modflow files
mfdis ='example.dis'
mfnam ='example.nam'
model_ws = 'data'

outdir = 'temp' # folder for writing output
basename = 'example' # base filename for output
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# output options
write_new_layer_bottoms = True # write new layer bottoms that are compatible with streambed bottoms
new_layer_bottoms_path = outdir
one_reach_per_cell = True # whether or not to allow multiple reaches in cells with meandering, confluences, etc.
# whether or not to consolidate conductances of SFR reaches
# (by adjusting hydraulic conductivity)
# so that only one reach in a cell has a non-zero conductance
# (True by defaul if one_reach_per_cell=True)
consolidate_conductance = False

# unit numbers
sfr_unit_number = 29,
sfr_leakage_unit_number = 53
sfr_output_unit_number = 66

# shapefiles defining the model grid, and the area that will include SFR
mf_grid = 'data/grid.shp'
mf_domain = 'data/domain.shp'

# NHDPlus input files (see the input requirements in the SFRmaker readme file
# (Note that multiple datasets can be supplied as lists;
# when the SFR area covers multiple drainage basins)
pfvaa_files = ['data/PlusFlowlineVAA.dbf']
plusflow_files = ['data/PlusFlow.dbf']
elevslope_files = ['data/elevslope.dbf']
flowlines = ['data/NHDFlowlines.shp']

# dem used for streambed elevations
dem = 'data/dem.tif'
dem_units_mult = 1. # convert dem elevation units to modflow units

# Flopy SpatialReference
xul, yul = 682668.5, 5147586.7 # UTM zone 15 N, nad27 meters (epsg code 26715)
nrow, ncol = 100, 160
dxdy = 250 # ft
delr = np.ones(nrow) * dxdy
delc = np.ones(ncol) * dxdy

# load model into Flopy, just the DIS and BAS packages
m = fm.Modflow.load(mfnam, model_ws=model_ws, forgive=False, check=False,
                    load_only=['DIS', 'BAS6'])
m.sr = flopy.utils.SpatialReference(delr=m.dis.delr,
                                    delc=m.dis.delc,
                                    xul=682668.5, yul=5147586.7,
                                    lenuni=1, # feet
                                    proj4_str='EPSG:26715', rotation=0)


ta = time.time()
# Read in the NHD datasets
nhd = NHDdata(NHDFlowline=flowlines, PlusFlowlineVAA=pfvaa_files, PlusFlow=plusflow_files,
              elevslope=elevslope_files,
              mf_grid=mf_grid, mf_units='feet',
              nrows=100, ncols=160,
              model_domain=mf_domain)

# Setup the segments, reaches, and other basic SFR parameters
nhd.to_sfr(minimum_length=10., one_reach_per_cell=one_reach_per_cell,
           consolidate_conductance=consolidate_conductance)

# Write out this information to Mat1 and Mat2 tables
nhd.write_tables(basename='{}/{}'.format(outdir, basename))

# Write out a shapefile that has the SFR linework
nhd.write_linework_shapefile(basename='{}/{}'.format(outdir, basename))
print("preproc finished in {:.2f}s\n".format(time.time() - ta))

# Mat 1 and Mat2 files generated from preproc.py above
m1 = 'exampleMat1.csv'
m2 = 'exampleMat2.csv'

# Read in Mat1 and Mat2 into an SFRdata object (postproc module)
# also include MODFLOW DIS file, NAM file, and path to model datasets
sfr = SFRdata(Mat1=m1, Mat2=m2, mfgridshp=mf_grid,
              mfdis=mfdis, mfpath=model_ws,
              mfnam=mfnam)

# For interior stream reaches (not at the ends of segments),
# assign streambed tops from the minimum DEM elevations in the model cell
sfr.reset_m1_streambed_top_from_dem(dem, dem_units_mult=dem_units_mult)

# Often the NHDPlus elevations don't match DEM at scales below 100k.
# reset the segment end elevations to the minimum dem elevation
# encountered in the end reach and all reaches upstream
# (to use the NHDPlus elevations for the segment ends, comment this out)
sfr.reset_segment_ends_from_dem()

# Remove any bumps in the DEM elevations,
# so that interior elevations in the segments always decrease in the downstream direction
sfr.smooth_interior_elevations()

# convert the m1 and m2 dataframes to input for Flopy
m1 = sfr.m1
m2 = sfr.m2
rd = fm.ModflowSfr2.get_empty_reach_data(len(m1))
sd = fm.ModflowSfr2.get_empty_segment_data(len(m2))

# dictionary to translate column names
sfrname = {
    'reach': 'ireach',
    'SFRlength': 'rchlen',
    'length': 'rchlen',
    'sbtop': 'strtop',
    'sbthick': 'strthick',
    'sbK': 'strhc1',
    'reachID': 'reachID',
    'outseg': 'outseg',
    'elevup': 'elevup',
    'elevdn': 'elevdn',
    'width1': 'width1',
    'width2': 'width2'
}

# populate the reach_data rec array
rd['i'] = m1.i.values
rd['j'] = m1.j.values
rd['k'] = 0
rd['iseg'] = m1.segment.values

for k, v in sfrname.items():
    if k in m1.columns and v in rd.dtype.names:
        rd[v] = m1[k].values

# populate the segment_data recarray
sd['nseg'] = m2.segment.values
sd['icalc'] = 1
sd['roughch'] = .037
for k, v in sfrname.items():
    if k in m2.columns:
        sd[v] = m2[k].values

# create the Flopy SFR object
sfr = fm.ModflowSfr2(m, nstrm=-len(m1), reach_data=rd, segment_data={0: sd},
                     isfropt=1, unit_number=sfr_unit_number,
                     ipakcb=sfr_leakage_unit_number,
                     istcb2=sfr_output_unit_number)
# assign K values from vertical
# compute the slopes
m.sfr.get_slopes()
# assign layers to the SFR
m.sfr.assign_layers(adjust_botms=True)

# write new external files for layer bottoms
if write_new_layer_bottoms:
    for i, l in enumerate(m.dis.botm.array):
        if not os.path.isdir(new_layer_bottoms_path):
            os.mkdir(new_layer_bottoms_path)
        ofn = os.path.join(new_layer_bottoms_path, 'botm{}.dat'.format(i + 1))
        np.savetxt(ofn, l, fmt='%.2f')
        print('wrote {}'.format(ofn))

# run diagnostics on the SFR and DIS packages
m.sfr.check('SFR.chk')
m.dis.check('DIS.chk')

# export shapefiles of SFR cells, routing connections and outlets
m.sfr.export('{}/{}_SFRcells.shp'.format(outdir, basename))
m.sfr.export_linkages('{}/{}_SFR_routing.shp'.format(outdir, basename))
m.sfr.export_outlets('{}/{}_SFR_outlets.shp'.format(outdir, basename))

# write the SFR package file
m.sfr.write_file('{}/{}.sfr'.format(outdir, basename))
