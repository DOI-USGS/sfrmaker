import os
import numpy as np
import pytest
import flopy
import flopy.modflow as fm
import sfrmaker


@pytest.fixture(scope='module')
def lines_from_NHDPlus(datapath):
    pfvaa_files = ['{}/badriver/PlusFlowlineVAA.dbf'.format(datapath)]
    plusflow_files = ['{}/badriver/PlusFlow.dbf'.format(datapath)]
    elevslope_files = ['{}/badriver/elevslope.dbf'.format(datapath)]
    flowlines = ['{}/badriver/NHDFlowlines.shp'.format(datapath)]

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
def model_grid():
    nrow, ncol = 112, 160
    dxy = 250
    sr = flopy.utils.SpatialReference(delr=np.ones(ncol)*dxy,
                                      delc=np.ones(nrow)*dxy,
                                      lenuni=1,
                                      xll=682688, yll=5139052, rotation=0,
                                      proj4_str='+init=epsg:26715')
    return sr


@pytest.fixture(scope='module')
def model(model_grid):
    m = fm.Modflow('example', model_ws=outdir)
    dis = fm.ModflowDis(m, nlay=1, nrow=nrow, ncol=nrow,
                        top=1000, botm=0)
    m.sr = model_grid


def test_make_sfr(outdir, model_grid,
                  lines_from_NHDPlus,
                  active_area_shapefile):
    grd = sfrmaker.StructuredGrid.from_sr(model_grid,
                                                      active_area=active_area_shapefile)
    sfr = lines_from_NHDPlus.to_sfr(grd)

    sfr.reach_data['strtop'] = sfr.interpolate_to_reaches('elevup', 'elevdn')
    sfr.get_slopes()

    sfr.write_package(outdir + 'example.sfr')
    sfr.export_cells(outdir + 'example_cells.shp')
    sfr.export_outlets(outdir + 'example_outlets.shp')
    sfr.export_lines(outdir + 'example_lines.shp')
    sfr.export_routing(outdir + 'example_routing.shp')
