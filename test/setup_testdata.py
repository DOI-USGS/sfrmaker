"""Script to set up test data from larger files
"""
from GISio import shp2df, df2shp


xul, yul = 520487.3, 1194668.3
nrow, ncol = 20, 20
dxy = 5280 * .3048
buf = 1e4
bounds = xul - buf, \
         yul - dxy * nrow - buf, \
         xul + dxy * ncol + buf, \
         yul + buf

df = shp2df('/Users/aleaf/Documents/MAP/repos/sfr_output/preprocessed/flowlines_gt20km/flowlines_gt20km_edited.shp',
            filter=bounds)
df2shp(df, 'data/map_test_flowlines.shp', epsg=5070)