"""Example script showing use of postproc.py module.
* Input tables (Mat1 and Mat2) containing routing information 
    have already been created using sfr_classes.py  
* A shapefile is supplied for model grid input instead of a MODFLOW DIS file,  
    as the model grid is rotated in this case 
    (the DIS file reader in SFRmaker does not support rotated grids)
"""

import sys
sys.path.insert(0, 'D:/ATLData/Documents/GitHub/SFR')
#sys.path.append('D:\JointBaseModel\SFRMakerData\TestModel')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rasterstats import zonal_stats
from postproc import SFRdata

path = 'D:/ATLData/SFR_testing/TestModel/'

# intantiate SFRdata object
sfr = SFRdata(Mat1=path + 'SFR_GWVmat1.txt',
              Mat2=path + 'SFR_GWVmat2.txt',
              mfgridshp='ForsyGrd.shp',
              mfgridshp_node_field='CellNum', 
              mfgridshp_row_field='ROW', mfgridshp_column_field='COLUMN')

# create columns in Mat2 of Min and Max elevation for each segment
# (these columns are not created by sfr_classes.py)
sfr.update_Mat2_elevations()

# update the reach elevations in Mat1 with minimum elevations from DEM
sfr.reset_m1_streambed_top_from_dem(dem=path + 'forsy_lid')

# trace routing from headwater segments to outlets; assign outlet to each segment
sfr.map_outsegs()

# creates table of segment confluences
sfr.map_confluences()

# smooth DEM elevations in segment interiors so that they decrease in downstream direction
sfr.smooth_interior_elevations()

# read in the DIS file (this is needed for some of the methods below;
# e.g. model top elevations are added to the stream profiles by default)
sfr.read_dis2(mfdis='Forsy.DIS', mfnam='ForsySFRMaker.nam')

# plot profiles of streambed elevations in comparison to model top and DEM minimum
sfr.plot_stream_profiles(add_profiles={'Minimum DEM elevation': 'landsurface'})

# enforce only one SFR conductance for each model cell 
# (other reaches in cell assigned near-zero conductance)
sfr.consolidate_conductance()

# adjust model grid so that all SFR reaches are in layer 1
# outputs a new DIS file for model
sfr.reset_model_top_2streambed(outdisfile='Forsy_adjusted_to_streambed.dis')

# run suite of diagnostics to test for common problems with SFR package
sfr.run_diagnostics()

# create shapefile for visualizing SFR package
sfr.write_shapefile(outshp='Forsy.shp', prj='ForsyGrd.prj')

# write updated tables
sfr.write_tables(basename='Forsy')

# write an SFR package file
sfr.write_sfr_package(basename='Forsy')
