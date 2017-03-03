"""Test for adding linework to sfr"""

import sys
sys.path.append('..')
sys.path.append('/Users/aleaf/Documents/GitHub/flopy3')
sys.path.append('D:/ATLData/Documents/GitHub/flopy')
import numpy as np
import pandas as pd
import fiona
from shapely.geometry import Point, LineString
from rtree import index
import flopy
from GISio import shp2df, df2shp, get_values_at_points, get_proj4
from preproc import NHDdata, renumber_segments, get_upsegs, lines

outpath = 'temp/'
nhd_lines = '../Examples/data/NHDflowlines.shp'
newlines = '../Examples/data/added_lines2.shp'
sfrlines = '../Examples/data/SFRlines2.shp'
grid = '../Examples/data/grid2.shp'
domain = '../Examples/data/domain2.shp'
#dem = 'D:/ATLData/BR/BadRiver/grid/dem/dem_utm_ft/w001001.adf'
#dem = '/Users/aleaf/Documents/BR/BadRiver/grid/dem/dem_utm_ft/w001001.adf'
dem = '../Examples/data/dem.tif'

#nhdpathGL = '/Users/aleaf/Documents/NHDPlus/NHDPlusGL/NHDPlus04/NHDPlusAttributes/'
nhdpathGL = '../Examples/data/'

pfvaa_files = nhdpathGL + 'PlusFlowlineVAA.dbf'
pf_files = nhdpathGL + 'PlusFlow.dbf'
elevslope_files = nhdpathGL + 'elevslope.dbf'
nhd_lines_proj4 = get_proj4(nhd_lines)

nhd = NHDdata(NHDFlowline=nhd_lines, PlusFlowlineVAA=pfvaa_files, PlusFlow=pf_files,
              elevslope=elevslope_files,
              mf_grid=grid, mf_units='feet',
              model_domain=domain,
              lines_proj4=nhd_lines_proj4)

nhd.to_sfr()
nhd.write_linework_shapefile(basename=outpath + 'SFRlines2')

lns = lines(newlines, mf_grid=grid, mf_grid_node_col='node', model_domain=domain)
lns.get_end_elevs_from_dem(dem)
m1, m2 = lns.append2sfr(nhd.m1, route2reach1=False, routing_tol=200)

nhd.m1 = nhd.m1.append(m1)
nhd.m2 = nhd.m2.append(m2)

nhd.write_linework_shapefile(basname=outpath + 'SFRlines_combined')
nhd.write_tables(basename=outpath + 'combined')