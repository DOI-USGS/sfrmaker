import os
from GISio import shp2df, get_proj4
from preproc import lines, NHDdata

nhd_lines = '../Examples/data/NHDflowlines.shp'
newlines = '../Examples/data/added_lines2.shp'
sfrlines = '../Examples/data/SFRlines2.shp'
grid = '../Examples/data/grid2.shp'
domain = '../Examples/data/domain2.shp'
#dem = 'D:/ATLData/BR/BadRiver/grid/dem/dem_utm_ft/w001001.adf'
dem = '/Users/aleaf/Documents/BR/BadRiver/grid/dem/dem_utm_ft/w001001.adf'

#nhdpathGL = '/Users/aleaf/Documents/NHDPlus/NHDPlusGL/NHDPlus04/NHDPlusAttributes/'
nhdpathGL = '../Examples/data/'

pfvaa_files = nhdpathGL + 'PlusFlowlineVAA.dbf'
pf_files = nhdpathGL + 'PlusFlow.dbf'
elevslope_files = nhdpathGL + 'elevslope.dbf'
nhd_lines_proj4 = get_proj4(nhd_lines)

if not os.path.isdir('temp'):
    os.makedirs('temp')

nhd = NHDdata(NHDFlowline=nhd_lines, PlusFlowlineVAA=pfvaa_files, PlusFlow=pf_files,
              elevslope=elevslope_files,
              mf_grid=grid, mf_units='feet',
              model_domain=domain,
              flowlines_proj4=nhd_lines_proj4)

nhd.to_sfr()
nhd.write_linework_shapefile(basename='../Examples/data/SFRlines2')

trim_buffer=20
#nsegstart = nhd.m1.segment.astype(int).max() + 1

lns = lines(newlines, mf_grid=grid, mf_grid_node_col='node', model_domain=domain)
lns.get_end_elevs_from_dem(dem)
m1, m2 = lns.append2sfr(nhd.m1, routing_tol=200)

nhd.m1 = nhd.m1.append(m1)
nhd.m2 = nhd.m2.append(m2)

nhd.renumber_segments()
nhd.write_linework_shapefile(basename='temp/junk.shp')
nhd.write_tables(basename='temp/junk')

j=2