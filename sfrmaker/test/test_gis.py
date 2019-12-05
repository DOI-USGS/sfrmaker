import os
import numpy as np
import pandas as pd
import pyproj
from shapely.geometry import Point
from gisutils import get_proj_str, df2shp, shp2df
from gisutils.shapefile import shp_properties
from ..gis import crs


def test_crs_eq():
    crs_4269_proj = crs(proj_str='+proj=longlat +datum=NAD83 +no_defs ')
    crs_26715_epsg = crs(epsg=26715)
    crs_26715_epsg_proj = crs(proj_str='epsg:26715')
    crs_26715_proj = crs(proj_str='+proj=utm +zone=15 +datum=NAD27 +units=m +no_defs ')
    crs_26715_prj = crs(prjfile='Examples/data/badriver/grid.shp')
    assert crs_4269_proj != crs_26715_epsg
    assert crs_4269_proj != crs_26715_epsg_proj
    assert crs_4269_proj != crs_26715_proj
    assert crs_4269_proj != crs_26715_prj
    assert crs_26715_epsg == crs_26715_epsg_proj
    assert crs_26715_epsg == crs_26715_proj
    assert crs_26715_epsg == crs_26715_prj


def test_isvalid():
    """
    With pyproj 2, all Proj instances are valid
    (error is raised in construction if not)
    https://github.com/pyproj4/pyproj/issues/304
    """
    crs_5070_epsg = crs(epsg=5070)
    assert crs_5070_epsg.is_valid

    junk = crs(proj_str='junk')
    assert not junk.is_valid


def test_crs_get_proj_str():
    crs_5070_epsg = crs(epsg=5070)
    assert crs_5070_epsg.proj_str.replace('+init=', '') == 'epsg:5070 +no_defs'


def test_rtree():
    from rtree import index