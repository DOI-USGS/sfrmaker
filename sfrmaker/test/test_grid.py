# TODO: add unit tests for grid.py
from rasterio import Affine

import flopy
import pytest
import sfrmaker

fm = flopy.modflow
from ..gis import get_authority_crs
from ..grid import StructuredGrid
from ..units import convert_length_units


def test_structuredgrid_from_flopy_mg():
    # make a flopy modelgrid instance
    # that represents the model grid
    data_dir = 'examples/tylerforks'
    m = fm.Modflow.load('tf.nam', model_ws='{}/tylerforks'.format(data_dir), load_only=['DIS'])
    mg = flopy.discretization.StructuredGrid(delr=m.dis.delr.array * .3048,  # cell spacing along a row
                                             delc=m.dis.delc.array * .3048,  # cell spacing along a column
                                             xoff=682688, yoff=5139052,  # lower left corner of model grid
                                             angrot=0,  # grid is unrotated
                                             crs=26715
                                             # projected coordinate system of model (UTM NAD27 zone 15 North)
                                             )

    grd2 = StructuredGrid.from_modelgrid(mg,
                                         active_area='{}/active_area.shp'.format(data_dir)
                                         )
    assert grd2.xul == 682688
    assert grd2.yul == 5139052 + mg.delc.sum()
    assert grd2.dx == mg.delr[0]
    assert grd2.dy == mg.delc[0]
    assert grd2.crs == get_authority_crs(26715)


def test_grid_epsg(shellmound_sfrmaker_grid):
    assert shellmound_sfrmaker_grid.crs.srs == 'EPSG:5070'


@pytest.fixture(scope='function')
def tyler_forks_grid_from_shapefile(tyler_forks_grid_shapefile, tylerforks_active_area_shapefile):
    grid = sfrmaker.StructuredGrid.from_shapefile(tyler_forks_grid_shapefile,
                                                  node_col='node',
                                                  icol='i',
                                                  jcol='j',
                                                  active_area=tylerforks_active_area_shapefile)
    return grid


@pytest.fixture(scope='function', params=[True, False])
def tyler_forks_grid_from_shapefile(request, tyler_forks_grid_shapefile, tylerforks_active_area_shapefile):
    if request.param:
        active_area = tylerforks_active_area_shapefile
    else:
        active_area = None
    grid = sfrmaker.StructuredGrid.from_shapefile(tyler_forks_grid_shapefile,
                                                  node_col='node',
                                                  icol='i',
                                                  jcol='j',
                                                  active_area=active_area)
    return grid


def test_structuredgrid_from_shapefile(tyler_forks_grid_from_shapefile, tylerforks_sfrmaker_grid_from_flopy):
    grid = tyler_forks_grid_from_shapefile
    grid_flopy = tylerforks_sfrmaker_grid_from_flopy
    assert grid.uniform
    assert grid == grid_flopy


@pytest.mark.skip(reason='not completed')
def test_unstructuredgrid_from_shapfile(tyler_forks_grid_shapefile,
                                        tylerforks_sfrmaker_grid_from_flopy):
    # TODO: test creating unstructured grid from same shapefile
    # with no row or column information passed
    pass