import os
import numpy as np
import pandas as pd
from shapely.geometry import Point
from ..gis import get_proj_str, df2shp, shp2df, shp_properties, crs


def test_get_proj_str(tmpdir):
    proj_str = '+proj=tmerc +lat_0=0 +lon_0=-90 +k=0.9996 +x_0=520000 +y_0=-4480000 +datum=NAD83 +units=m +no_defs '
    f = os.path.join(tmpdir, 'junk.shp')
    df2shp(pd.DataFrame({'id': [0],
                         'geometry': [Point(0, 0)]
                         }),
           f, proj_str=proj_str)
    proj42 = get_proj_str(f.replace('shp', 'prj'))
    assert proj42 == proj_str


def test_shp_properties():
    df = pd.DataFrame({'reach': [1], 'value': [1.0], 'name': ['stuff']}, index=[0])
    df = df[['name', 'reach', 'value']].copy()
    assert [d.name for d in df.dtypes] == ['object', 'int64', 'float64']
    assert shp_properties(df) == {'name': 'str', 'reach': 'int', 'value': 'float'}


def test_integer_dtypes(tmpdir):

    # verify that pandas is recasting numpy ints as python ints when converting to dict
    # (numpy ints invalid for shapefiles)
    d = pd.DataFrame(np.ones((3, 3)), dtype=int).astype(object).to_dict(orient='records')
    for i in range(3):
        assert isinstance(d[i][0], int)

    df = pd.DataFrame({'r': np.arange(100), 'c': np.arange(100)})
    f = '{}/ints.dbf'.format(tmpdir)
    df2shp(df, f)
    df2 = shp2df(f)
    assert np.all(df == df2)


def test_boolean_dtypes(tmpdir):

    df = pd.DataFrame([False, True]).transpose()
    df.columns = ['true', 'false']
    f = '{}/bool.dbf'.format(tmpdir)
    df2shp(df, f)
    df2 = shp2df(f, true_values='True', false_values='False')
    assert np.all(df == df2)


def test_crs_eq():

    crs_4269_proj = crs(proj_str='+proj=longlat +datum=NAD83 +no_defs ')
    crs_26715_epsg = crs(epsg=26715)
    crs_26715_epsg_proj = crs(proj_str='+init=epsg:26715')
    crs_26715_proj = crs(proj_str='+proj=utm +zone=15 +datum=NAD27 +units=m +no_defs ')
    crs_26715_prj = crs(prjfile='Examples/data/badriver/grid.shp')
    assert crs_4269_proj != crs_26715_epsg
    assert crs_4269_proj != crs_26715_epsg_proj
    assert crs_4269_proj != crs_26715_proj
    assert crs_4269_proj != crs_26715_prj
    assert crs_26715_epsg == crs_26715_epsg_proj
    assert crs_26715_epsg == crs_26715_proj
    assert crs_26715_epsg == crs_26715_prj


def test_get_proj_str():

    crs_5070_epsg = crs(epsg=5070)
    assert crs_5070_epsg.proj_str == '+init=epsg:5070 +no_defs'
