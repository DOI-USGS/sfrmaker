"""
Make a SFR package for a prexisting MODFLOW model,
using a shapefile to specify the grid
"""
from pathlib import Path

import flopy
import numpy as np

fm = flopy.modflow
from sfrmaker import Lines, StructuredGrid
from sfrmaker.utils import assign_layers

# NHDPlus input files (see the input requirements in the SFRmaker readme file
# (Note that multiple datasets can be supplied as lists
# when the SFR area covers multiple drainage basins)
data_dir = 'tylerforks'
pfvaa_files = ['{}/PlusFlowlineVAA.dbf'.format(data_dir)]
plusflow_files = ['{}/PlusFlow.dbf'.format(data_dir)]
elevslope_files = ['{}/elevslope.dbf'.format(data_dir)]
flowlines = ['{}/NHDFlowlines.shp'.format(data_dir)]

# DEM for sampling streambed top elevations
dem = '{}/dem_26715.tif'.format(data_dir)

# output folder
outdir = Path('temp/')
outdir.mkdir(exist_ok=True)

# make an instance of the sfrmaker.lines class from NHDPlus data
# use a shapefile of the model grid to filter what is read in
# (only flowlines within the bounding box of the grid)
lns = Lines.from_nhdplus_v2(NHDFlowlines=flowlines,
                            PlusFlowlineVAA=pfvaa_files,
                            PlusFlow=plusflow_files,
                            elevslope=elevslope_files,
                            bbox_filter='{}/grid.shp'.format(data_dir))

# load the model using flopy
m = fm.Modflow.load('tf.nam', model_ws='{}/tylerforks'.format(data_dir))

# make a sfrmaker.StructuredGrid instance from a shapefile of the model grid
# for a structured grid, attribute fields with row and column information must be passed
# active_area is a polygon that specifies the area where SFR reaches will be populated
grd = StructuredGrid.from_shapefile(shapefile='{}/grid.shp'.format(data_dir),
                                    icol='i',
                                    jcol='j',
                                    active_area='{}/active_area.shp'.format(data_dir)
                                    )

# from the lines and StructuredGrid instances, make a sfrmaker.sfrdata instance
# (lines are intersected with the model grid and converted to reaches, etc.)
sfr = lns.to_sfr(grd, model=m)
sfr.set_streambed_top_elevations_from_dem(dem, dem_z_units='meters')

# assign layers to the sfr reaches
# if the model bottom is changed, write out a new one and update model
botm = m.dis.botm.array.copy()
layers, new_botm = assign_layers(sfr.reach_data, botm_array=botm)
sfr.reach_data['k'] = layers
if new_botm is not None:
    botm[-1] = new_botm
    np.savetxt('{}/external/botm{}.dat'.format(m.model_ws,
                                               m.nlay - 1),
               new_botm, fmt='%.2f')
    sfr.modflow_sfr2.parent.dis.botm = botm

# run Flopy diagnostics (refresh the modflow_sfr2 instance first)
sfr.create_modflow_sfr2(model=m)
sfr.modflow_sfr2.check()

# write the SFR package
# turn on writing of SFR water balance ascii output by specifying unit no. for istcb2
sfr.write_package(istcb2=223)  # writes a sfr file to the model workspace
m.write_name_file()  # write new version of name file with sfr package

# wite shapefiles for visualization
sfr.export_cells(outdir + 'example_cells.shp')
sfr.export_outlets(outdir + 'example_outlets.shp')
sfr.export_lines(outdir + 'example_lines.shp')
sfr.export_routing(outdir + 'example_routing.shp')
