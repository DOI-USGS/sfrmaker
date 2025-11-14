import os
import numpy as np
import pytest
from gisutils import get_authority_crs
from sfrmaker.gis import get_bbox


def test_get_bbox(project_root_path):
    shapefile = os.path.join(project_root_path, 'examples/tylerforks/grid.shp')
    crs = get_authority_crs(4269)
    bbox = get_bbox(shapefile, dest_crs=crs)
    assert np.allclose(bbox, (-90.62442575352304, 46.37890212020774, -90.46249896050521, 46.458360301848685))


#@pytest.fixture(scope='module')
#def intersected(tylerforks_sfrmaker_grid_from_flopy, tylerforks_lines_from_NHDPlus):
#    #results = intersect()
#    pass


@pytest.mark.skip(reason='still working on faster intersection method')
def test_intersect():
    from rtree import index
    pass