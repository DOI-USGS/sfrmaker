import os
import numpy as np
import pytest
import flopy.modflow as fm
import sfrmaker

#TODO: make tests more rigorous

@pytest.fixture(scope='module')
def lines_from_NHDPlus(datapath):
    pfvaa_files = ['{}/badriver/PlusFlowlineVAA.dbf'.format(datapath)]
    plusflow_files = ['{}/badriver/PlusFlow.dbf'.format(datapath)]
    elevslope_files = ['{}/badriver/elevslope.dbf'.format(datapath)]
    flowlines = ['{}/badriver/NHDflowlines.shp'.format(datapath)]

    lns = sfrmaker.lines.from_NHDPlus_v2(NHDFlowlines=flowlines,
                                PlusFlowlineVAA=pfvaa_files,
                                PlusFlow=plusflow_files,
                                elevslope=elevslope_files,
                                filter='{}/badriver/grid.shp'.format(datapath))
    return lns


@pytest.fixture(scope='module')
def active_area_shapefile(datapath):
    return '{}/badriver/active_area.shp'.format(datapath)


@pytest.fixture(scope='module')
def model(model_grid):
    m = fm.Modflow('example', model_ws=outdir)
    dis = fm.ModflowDis(m, nlay=1, nrow=nrow, ncol=nrow,
                        top=1000, botm=0)
    m.sr = model_grid
    return model


@pytest.fixture(scope='module')
def sfrmaker_grid_from_sr(model_grid, active_area_shapefile):
    grid = sfrmaker.StructuredGrid.from_sr(model_grid,
                                          active_area=active_area_shapefile)
    return grid


@pytest.fixture(scope='module')
def sfrmaker_grid_from_shapefile(grid_shapefile):
    #grid = sfrmaker.StructuredGrid.from_sr(model_grid,
    #                                       active_area=active_area_shapefile)
    #return grid
    pass

#@pytest.mark.parameterize('grid', [sfrmaker_grid_from_sr,
#                                   sfrmaker_grid_from_shapefile
#                                   ])
def test_make_sfr(outdir, sfrmaker_grid_from_sr,
                  lines_from_NHDPlus,
                  active_area_shapefile):

    sfr = lines_from_NHDPlus.to_sfr(grid=sfrmaker_grid_from_sr)

    sfr.reach_data['strtop'] = sfr.interpolate_to_reaches('elevup', 'elevdn')
    sfr.get_slopes()

    sfr.write_package(outdir + 'example.sfr')
    sfr.export_cells(outdir + 'example_cells.shp')
    sfr.export_outlets(outdir + 'example_outlets.shp')
    sfr.export_lines(outdir + 'example_lines.shp')
    sfr.export_routing(outdir + 'example_routing.shp')
