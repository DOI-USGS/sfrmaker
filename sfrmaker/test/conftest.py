import os
import numpy as np
import pandas as pd
import pytest
import flopy
import flopy.modflow as fm
from flopy.utils import SpatialReference


@pytest.fixture(scope="session")
def datapath():
    return 'Examples/data'


@pytest.fixture(scope="session")
def outdir():
    # output folder
    outdir = 'sfrmaker/test/temp/'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    return outdir


@pytest.fixture(scope="function")
def sfr_test_numbering():
    rd = pd.DataFrame()
    rd['i'] = [3, 4, 5,
              7, 8, 9,
              0, 1, 2,
              4, 4, 5,
              0, 0, 0,
              3, 4, 5,
              0, 1, 2,
              4, 5, 6,
              2, 2, 2]
    rd['j'] = [0, 1, 2,
              6, 6, 6,
              6, 6, 6,
              3, 4, 5,
              9, 8, 7,
              6, 6, 6,
              0, 0, 0,
              6, 6, 6,
              9, 8, 7]
    rd['iseg'] = sorted(list(range(1, 10)) * 3)
    rd['ireach'] = [1, 2, 3] * 9

    sd = pd.DataFrame()
    sd['nseg'] = range(1, 10)
    sd['outseg'] = [4, 0, 6, 8, 3, 8, 1, 2, 8]
    return rd, sd


@pytest.fixture(scope='module')
def tyler_forks_grid_shapefile(datapath):
    return '{}/badriver/grid.shp'.format(datapath)


@pytest.fixture(scope='module')
def tylerforks_model(datapath):
    m = fm.Modflow.load('tf.nam', model_ws='{}/badriver/tylerforks'.format(datapath))
    return m


@pytest.fixture(scope='module')
def tylerforks_model_grid(tyler_forks_grid_shapefile, tylerforks_model):
    m = tylerforks_model
    sr = flopy.utils.SpatialReference(delr=m.dis.delr.array,  # cell spacing along a row
                                      delc=m.dis.delc.array,  # cell spacing along a column
                                      lenuni=1,  # model units of feet
                                      xll=682688, yll=5139052,  # lower left corner of model grid
                                      rotation=0,  # grid is unrotated
                                      proj4_str='+init=epsg:26715'
                                      # projected coordinate system of model (UTM NAD27 zone 15 North)
                                      )
    m.sr = sr
    return sr

