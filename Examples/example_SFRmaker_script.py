"""Example workflow for generating an SFR package
from NHDPlus data, a DEM, and MODFLOW grid.
"""
import sys

# if preproc.py and postproc.py are in a different folder from this script,
# add that folder to the path so python can find it
sys.path.append('../')
import time
from GISio import shp2df, get_proj4
from preproc import NHDdata
from postproc import SFRdata

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

ta = time.time()
# Read in the NHD datasets
nhd = NHDdata(NHDFlowline=flowlines, PlusFlowlineVAA=pfvaa_files, PlusFlow=plusflow_files,
              elevslope=elevslope_files,
              mf_grid=mf_grid, mf_units='feet',
              nrows=100, ncols=160,
              model_domain=mf_domain)

# Setup the segments, reaches, and other basic SFR parameters
nhd.to_sfr()

# Write out this information to Mat1 and Mat2 tables
nhd.write_tables(basename='example')

# Write out a shapefile that has the SFR linework
nhd.write_linework_shapefile(basename='SFRlines')
print("preproc finished in {:.2f}s\n".format(time.time() - ta))

# Mat 1 and Mat2 files generated from preproc.py above
m1 = 'exampleMat1.csv'
m2 = 'exampleMat2.csv'

# Read in Mat1 and Mat2 into an SFRdata object (postproc module)
# also include MODFLOW DIS file, NAM file, and path to model datasets
sfr = SFRdata(Mat1=m1, Mat2=m2, mfgridshp=mf_grid,
              mfdis='example.dis', mfpath='data',
              mfnam='example.nam')

# For interior stream reaches (not at the ends of segments),
# assign streambed tops from the minimum DEM elevations in the model cell
sfr.reset_m1_streambed_top_from_dem(dem, dem_units_mult=dem_units_mult)

# Often the NHDPlus elevations don't match DEM at scales below 100k.
# reset the segment end elevations to the minimum dem elevation
# encountered in the end reach and all reaches upstream
# (to use the NHDPlus elevations for the segment ends, comment this out)
sfr.reset_segment_ends_from_dem()

# Create array listing all unique segment sequences from headwaters to outlet
# used by other methods, and ensures that there is no circular routing
sfr.map_outsegs()

# Remove any bumps in the DEM elevations,
# so that interior elevations in the segments always decrease in the downstream direction
sfr.smooth_interior_elevations()

# Create a PDF showing the elevation profile of each unique segment sequence,
# in comparison to the mean DEM elevation in each cell along the profile
# Note: this is optional but allows for inspection to verify that the elevations are reasonable
# on large networks (with thousands of sequences) this step may take an hour or two.
sfr.plot_stream_profiles()

# In cells with multiple SFR reaches (at confluences), put all conductance in the dominant reach
# (to avoid circular routing)
sfr.consolidate_conductance()

# Put all SFR reaches in layer 1, and adjust layer bottoms in SFR cells accordingly,
# so that no streambed bottoms are below the bottom of layer 1, and
# so that there are no cell thicknesses less than 1
# Note: This produces a new MODFLOW DIS file with the suffix "_adjusted to streambed.dis"
# (unless another suffix is specified)
sfr.reset_model_top_2streambed(minimum_thickness=1)

# write out a shapefile of the SFR dataset
sfr.write_shapefile('SFR_package.shp')

# write out updated Mat1 and Mat2 tables
sfr.write_tables(basename='example_pp_')

# write the SFR package file
sfr.write_sfr_package(basename='data/example')

# flopy now includes an SFR module, with a checker
# run the flopy suite of diagnostics
import flopy
m = flopy.modflow.Modflow.load('example.nam', model_ws='data')
#dis = flopy.modflow.ModflowDis.load('ozarkgwmod/3-Input/A-calibration/ozark_adjusted_to_streambed.dis', m)
sfr = flopy.modflow.ModflowSfr2.load('data/example.sfr', m)
sfr.check(f='flopy_SFR_check')
