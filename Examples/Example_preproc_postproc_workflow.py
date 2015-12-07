"""Example workflow for generating an SFR package 
from NHDPlus data, a DEM, and MODFLOW grid.
"""
import sys

# if preproc.py and postproc.py are in a different folder from this script,
# add that folder to the path so python can find it
sys.path.append('/Users/aleaf/Documents/GitHub/SFR')
import time
from GISio import shp2df, get_proj4
from preproc import NHDdata
from postproc import SFRdata

# shapefiles defining the model grid, and the area that will include SFR
mf_grid = 'SFR/modelgrid.shp'
mf_domain = 'SFR/ozarkBoundary.shp'

# NHDPlus input files (see the input requirements in the SFRmaker readme file
# (Note that multiple datasets can be supplied as lists; 
# when the SFR area covers multiple drainage basins)
pfvaa_files = 'SFR/nhdozarkAtts.dbf'
plusflow_path = '/Users/aleaf/Documents/NHDPlus/NHDPlusMS/NHDPlus{}/NHDPlusAttributes/PlusFlow.dbf'
regions = ['07', '08', '10L', '11']
pf_files = [plusflow_path.format(r) for r in regions]
elevslope_files = 'SFR/ozarkselevslope.dbf'
flowlines = 'SFR/NHDmrgBufClpAlbGT1.shp'

ta = time.time()
# Read in the NHD datasets
nhd = NHDdata(NHDFlowline=flowlines, PlusFlowlineVAA=pfvaa_files, PlusFlow=pf_files,
              elevslope=elevslope_files,
              mf_grid=mf_grid, mf_grid_node_col='cellnum', mf_units='feet',
              model_domain=mf_domain)
# Setup the segments, reaches, and other basic SFR parameters
nhd.to_sfr()

# Write out this information to Mat1 and Mat2 tables
nhd.write_tables(basename='ozark')

# Write out a shapefile that has the SFR linework
nhd.write_linework_shapefile(basename='SFRlines')
print("preproc finished in {:.2f}s\n".format(time.time() - ta))

# Mat 1 and Mat2 files generated from preproc.py above
m1 = 'ozarkMat1.csv'
m2 = 'ozarkMat2.csv'

# Read in Mat1 and Mat2 into an SFRdata object (postproc module)
# also include MODFLOW DIS file, NAM file, and path to model datasets
sfr = SFRdata(Mat1=m1, Mat2=m2, mfgridshp=mf_grid,
              mfdis='ozark.dis', mfpath='ozarkgwmod/3-Input/A-calibration/',
              mfnam='ozark.nam',
              mfgridshp_node_field='cellnum',
              mfgridshp_row_field='irow',
              mfgridshp_column_field='icol')

# For interior stream reaches (not at the ends of segments), 
# assign streambed tops from the minimum DEM elevations in the model cell
sfr.reset_m1_streambed_top_from_dem('dem30m1.tif', dem_units_mult=1/0.3048)

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
sfr.write_shapefile('SFR_postproc.shp')

# run some diagnostics to look for common problems
sfr.run_diagnostics(model_domain=mf_domain, sfr_linework_shapefile='SFRlines.shp')

# write out updated Mat1 and Mat2 tables
sfr.write_tables(basename='ozark_pp_')

# write the SFR package file
sfr.write_sfr_package(basename='ozark')

# the development branch of flopy now includes an SFR module, with a checker
# run the flopy suite of diagnostics as well
import flopy
m = flopy.modflow.Modflow(model_ws='ozarkgwmod/3-Input/A-calibration/')
dis = flopy.modflow.ModflowDis.load('ozarkgwmod/3-Input/A-calibration/ozark_adjusted_to_streambed.dis', m)
sfr = flopy.modflow.ModflowSfr2.load('ozark.sfr', m)
sfr.check(f='flopy_SFR_check')
