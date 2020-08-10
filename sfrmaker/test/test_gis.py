import os
import numpy as np
import pytest
from sfrmaker.gis import CRS, get_bbox, intersect


# basic test that different input options don't crash CRS.__init__
@pytest.mark.parametrize('kwargs', ({},
                                    {'crs_dict': {'proj': 'utm',
                                                  'zone': 16,
                                                  'datum': 'NAD83',
                                                  'units': 'm',
                                                  'no_defs': None,
                                                  'type': 'crs'}},
                                    {'epsg': 26916},
                                    {'proj_str': 'epsg:26916'},
                                    {'prjfile': 'sfrmaker/test/data/shellmound/flowlines.prj'},
                                    {'crs_dict': {'epsg': 26916},
                                     'epsg': 26916,
                                     'proj_str': 'epsg:26916',
                                     'prjfile': 'sfrmaker/test/data/shellmound/flowlines.prj'
                                     },
                                    {'wkt': ('PROJCS["unnamed",GEOGCS["WGS 84",DATUM["unknown",\
                                             SPHEROID["WGS84",6378137,298.257223563]],PRIMEM["Greenwich",0],\
                                             UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]]],\
                                             PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],\
                                             PARAMETER["central_meridian",-90],PARAMETER["scale_factor",0.9996],\
                                             PARAMETER["false_easting",520000],PARAMETER["false_northing",-4480000],\
                                             UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],\
                                             AXIS["Northing",NORTH]]')},
                                    {'crs_dict': {'ellps': 'WGS84', 
                                                  'k': 0.9996, 
                                                  'lat_0': 0, 
                                                  'lon_0': -90, 
                                                  'no_defs': True, 
                                                  'proj': 'tmerc', 
                                                  'units': 'm', 
                                                  'x_0': 520000, 
                                                  'y_0': -4480000}}
)
                         )
def test_CRS(kwargs):
    # test init with no arguments
    if 'wkt' in kwargs:
        j=2
    crs = CRS(**kwargs)
    assert isinstance(crs, CRS)
    if len(kwargs) == 0:
        assert not crs.is_valid
    else:
        assert crs.is_valid
    assert isinstance(crs.__repr__(), str)
    j=2


def test_crs_eq():
    crs_4269_proj = CRS(proj_str='+proj=longlat +datum=NAD83 +no_defs ')
    crs_26715_epsg = CRS(epsg=26715)
    crs_26715_epsg_proj = CRS(proj_str='epsg:26715')
    crs_26715_proj = CRS(proj_str='+proj=utm +zone=15 +datum=NAD27 +units=m +no_defs ')
    crs_26715_prj = CRS(prjfile='examples/tylerforks/grid.shp')
    assert crs_4269_proj != crs_26715_epsg
    assert crs_4269_proj != crs_26715_epsg_proj
    assert crs_4269_proj != crs_26715_proj
    assert crs_4269_proj != crs_26715_prj
    assert crs_26715_epsg == crs_26715_epsg_proj
    assert crs_26715_epsg == crs_26715_proj
    assert crs_26715_epsg == crs_26715_prj


def test_crs_units():
    crs_4269_proj = CRS(proj_str='+proj=longlat +datum=NAD83 +no_defs ')
    assert crs_4269_proj.length_units == 'degree'
    crs_26715_epsg = CRS(epsg=26715)
    assert crs_26715_epsg.length_units == 'meters'
    crs_26715_epsg_proj = CRS(proj_str='epsg:26715')
    assert crs_26715_epsg_proj.length_units == 'meters'
    crs_26715_proj = CRS(proj_str='+proj=utm +zone=15 +datum=NAD27 +units=m +no_defs ')
    assert crs_26715_proj.length_units == 'meters'
    crs_26715_prj = CRS(prjfile='examples/tylerforks/grid.shp')
    assert crs_26715_prj.length_units == 'meters'


def test_is_valid():
    """
    With pyproj 2, all Proj instances are valid
    (error is raised in construction if not)
    https://github.com/pyproj4/pyproj/issues/304
    """
    crs_5070_epsg = CRS(epsg=5070)
    assert crs_5070_epsg.is_valid


@pytest.mark.xfail(reason="Invalid CRS class can't be instantiated")
def test_invalid():
    junk = CRS(proj_str='junk')
    assert not junk.is_valid


def test_crs_get_proj_str():
    crs_5070_epsg = CRS(epsg=5070)
    assert crs_5070_epsg.proj_str == 'EPSG:5070'


def test_get_bbox(project_root_path):
    shapefile = os.path.join(project_root_path, 'examples/tylerforks/grid.shp')
    crs = CRS(epsg=4269)
    bbox = get_bbox(shapefile, dest_crs=crs)
    assert np.allclose(bbox, (-90.62442575352304, 46.37890212020774, -90.46249896050521, 46.458360301848685))


@pytest.mark.skip(reason='still working on faster intersection method')
@pytest.fixture(scope='module')
def intersected(tylerforks_sfrmaker_grid_from_flopy, tylerforks_lines_from_NHDPlus):
    #results = intersect()
    pass


@pytest.mark.skip(reason='still working on faster intersection method')
def test_intersect():
    from rtree import index
    pass