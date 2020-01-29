# TODO: add unit tests for grid.py
import flopy
fm = flopy.modflow
from ..gis import crs
from ..grid import StructuredGrid
from ..units import convert_length_units


def test_structuredgrid_from_flopy_mg():
    # make a flopy.utils.reference.SpatialReference instance
    # that represents the model grid
    data_dir = 'Examples/data/badriver'
    m = fm.Modflow.load('tf.nam', model_ws='{}/tylerforks'.format(data_dir), load_only=['DIS'])
    mg = flopy.discretization.StructuredGrid(delr=m.dis.delr.array * .3048,  # cell spacing along a row
                                             delc=m.dis.delc.array * .3048,  # cell spacing along a column
                                             xoff=682688, yoff=5139052,  # lower left corner of model grid
                                             angrot=0,  # grid is unrotated
                                             proj4='epsg:26715'
                                             # projected coordinate system of model (UTM NAD27 zone 15 North)
                                             )

    grd2 = StructuredGrid.from_modelgrid(mg,
                                         active_area='{}/active_area.shp'.format(data_dir)
                                         )
    assert grd2.xul == 682688
    assert grd2.yul == 5139052 + mg.delc.sum()
    assert grd2.dx == mg.delr[0]
    assert grd2.dy == mg.delc[0]
    assert grd2.crs == crs(epsg=26715)
