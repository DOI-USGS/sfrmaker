from collections.abc import Iterable
import os
from packaging import version
from pathlib import Path
import warnings
import time
import traceback
import fiona
import numpy as np
import pyproj
from shapely.geometry import shape, Polygon, box
from shapely.ops import unary_union
from shapely.validation import make_valid
import gisutils
from gisutils import df2shp, shp2df, project, get_shapefile_crs, get_authority_crs
import sfrmaker

if version.parse(gisutils.__version__) < version.parse('0.2.2'):
    warnings.warn('Automatic reprojection functionality requires gis-utils >= 0.2.2'
                  '\nPlease pip install --upgrade gis-utils')


def get_crs(prjfile=None, crs=None, **kwargs):
    if 'epsg' in kwargs:
        warnings.warn(
        "epsg argument is deprecated, "
        "use crs instead",
        PendingDeprecationWarning,
    )
        crs = kwargs['epsg']
    if 'proj_str' in kwargs:
        warnings.warn(
        "proj_str argument is deprecated, "
        "use crs instead",
        PendingDeprecationWarning,
    )
        crs = kwargs['proj_str']
    if crs is not None:
        crs = get_authority_crs(crs)
    elif prjfile is not None:
        prjfile_crs = get_shapefile_crs(prjfile)
        if (crs is not None) and (crs != prjfile_crs):
                raise ValueError('Different coordinate reference systems '
                                 f'in crs argument and supplied projection file: {prjfile}\n'
                                 f'\nuser supplied crs: {crs}  !=\ncrs from projection file: {prjfile_crs}'
                                 )
        else:
            crs = prjfile_crs
    return crs


def build_rtree_index(geom):
    """Builds an :class:`rtree.index.Index` (spatial index) object. 
    Useful for multiple intersections with same index.

    Parameters
    ==========
    geom : list
        list of shapely geometry objects
    Returns
        idx : rtree spatial index object
    """
    from rtree import index

    # build spatial index for items in geom1
    print('\nBuilding spatial index...')
    ta = time.time()
    idx = index.Index()
    for i, g in enumerate(geom):
        idx.insert(i, g.bounds)
    print("finished in {:.2f}s".format(time.time() - ta))
    return idx


def export_reach_data(reach_data, grid, filename,
                      nodes=None, geomtype='Polygon'):
    """Generic method for exporting data to a shapefile; joins
    attributes in reach_data to geometries in grid using node numbers.
    """
    assert grid is not None, "need grid attribute for export"
    if nodes is not None:
        keep = [True if n in nodes else False for n in reach_data.node]
        rd = reach_data.loc[keep].copy()
    else:
        rd = reach_data.copy()
    assert isinstance(grid, sfrmaker.grid.Grid), "grid needs to be an sfrmaker.Grid instance"
    assert np.array_equal(grid.df.node.values, np.arange(grid.size))
    assert np.array_equal(grid.df.node.values, grid.df.index.values)
    polygons = grid.df.loc[rd.node, 'geometry'].values
    if geomtype.lower() == 'polygon':
        rd['geometry'] = polygons
    elif geomtype.lower() == 'point':
        rd['geometry'] = [p.centroid for p in polygons]
    else:
        raise ValueError('Unrecognized geomtype "{}"'.format(geomtype))
    df2shp(rd, filename, crs=grid.crs)


def intersect_rtree(geom1, geom2, index=None):
    """Intersect features in geom1 with those in geom2. For each feature in geom2, return a list of
     the indices of the intersecting features in geom1.

    Parameters
    ----------
    geom1 : list
        list of shapely geometry objects
    geom2 : list
        list of shapely geometry objects to be intersected with features in geom1
    index : rtree spatial index
        Option to use an r-tree spatial index that has already been created for geom1.
        The :func:`sfrmaker.gis.build_rtree_index` function can be used to
        create a spatial index for a list of shapely geometriy objects.
        by default, None.

    Returns
    -------
    A list of the same length as geom2; containing for each feature in geom2,
    a list of indicies of intersecting geometries in geom1.
    """

    #make certain that the objects in geom1 and geom2 are all considered valid shapely geometries using "make_valid()"
    for i in range(len(geom1)): geom1[i] = make_valid(geom1[i])
    for i in range(len(geom2)): geom2[i] = make_valid(geom2[i])
    
    if index is None:
        idx = build_rtree_index(geom1)
    else:
        idx = index
    isfr = []
    print('\nIntersecting {} features...'.format(len(geom2)))
    ta = time.time()
    for pind, geom in enumerate(geom2):
        print('\r{}'.format(pind + 1), end='')
        # test for intersection with bounding box of each feature in geom2 using spatial index
        inds = [i for i in idx.intersection(geom.bounds)]
        # test each feature inside the bounding box for intersection with the geometry
        inds = [i for i in inds if geom1[i].intersects(geom)]
        isfr.append(inds)
    print("\nfinished in {:.2f}s".format(time.time() - ta))
    return isfr


def intersect(geom1, geom2):
    """Same as intersect_rtree, except without spatial indexing. Fine for smaller datasets,
    but scales by 10^4 with the side of the problem domain.

    Parameters
    ----------
    geom1 : list
        list of shapely geometry objects
    geom2 : list
        list of shapely geometry objects to be intersected with features in geom1

    Returns
    -------
    A list of the same length as geom2; containing for each feature in geom2,
    a list of indicies of intersecting geometries in geom1.
    """

    #make certain that the objects in geom1 and geom2 are all considered valid shapely geometries using "make_valid()"
    for i in range(len(geom1)): geom1[i] = make_valid(geom1[i])
    for i in range(len(geom2)): geom2[i] = make_valid(geom2[i])
    
    isfr = []
    ngeom1 = len(geom1)
    print('Intersecting {} features...'.format(len(geom2)))
    ta = time.time()
    for i, g in enumerate(geom2):
        print('\r{}'.format(i + 1), end='')
        intersects = np.array([r.intersects(g) for r in geom1])
        inds = list(np.arange(ngeom1)[intersects])
        isfr.append(inds)
    print("\nfinished in {:.2f}s".format(time.time() - ta))
    return isfr


def parse_units_from_proj_str(proj_str):
    units = None
    from pyproj import CRS
    crs = CRS.from_string(proj_str)
    try:
        # need this because preserve_units doesn't seem to be
        # working for complex proj strings.  So if an
        # epsg code was passed, we have no choice, but if a
        # proj string was passed, we can just parse it

        if "EPSG" in proj_str.upper():
            import pyproj
            from pyproj import CRS
            crs = CRS.from_epsg(4326)
            crs = pyproj.Proj(proj_str,
                              preseve_units=True,
                              errcheck=True)
            proj_str = crs.srs
        else:
            proj_str = proj_str
        # http://proj.org/parameters.html#units
        # from proj source code
        # "us-ft", "0.304800609601219", "U.S. Surveyor's Foot",
        # "ft", "0.3048", "International Foot",
        if "units=m" in proj_str:
            units = "meters"
        elif "units=ft" in proj_str or \
                "units=us-ft" in proj_str or \
                "to_meters:0.3048" in proj_str:
            units = "feet"
        return units
    except:
        pass


def read_polygon_feature(feature, dest_crs=None, feature_crs=None):
    """Read a geometric feature from a shapefile, shapely geometry object,
    or collection of shapely geometry objects. Reproject to dest_crs
    if the feature is in a different CRS.

    Parameters
    ----------
    feature : shapely Polygon, list of Polygons, or shapefile path
            Polygons must be in same CRS as linework; shapefile
            features will be reprojected if their crs is different.
    dest_crs : instance of sfrmaker.crs
        Output CRS for the feature.

    Returns
    -------
    feature : shapely geometry object
    """
    if isinstance(feature, str) or isinstance(feature, Path):
        feature_crs = get_shapefile_crs(feature)
        geoms = shp2df(feature)['geometry'].values
        feature = unary_union(geoms)
    elif isinstance(feature, Iterable):
        if isinstance(feature[0], dict):
            try:
                feature = [shape(f) for f in feature]
            except Exception as ex:
                print(ex)
                print("Supplied dictionary doesn't appear to be valid GeoJSON.")
        feature = unary_union(feature)
    elif isinstance(feature, dict):
        try:
            feature = shape(feature)
        except Exception as ex:
            print(ex)
            print("Supplied dictionary doesn't appear to be valid GeoJSON.")
    elif isinstance(feature, Polygon):
        pass
    else:
        raise TypeError("Unrecognized feature input.")
    if feature_crs is not None and dest_crs is not None and feature_crs != dest_crs:
        feature = project(feature, feature_crs, dest_crs)
    return feature.buffer(0)


def get_bbox(feature, dest_crs):
    """Get bounding box for a Polygon feature.

    Parameters
    ----------
    feature : str (shapefile path), shapely Polygon or GeoJSON
    dest_crs  : obj
        Coordinate reference system of the head observation locations.
        A Python int, dict, str, or :class:`pyproj.crs.CRS` instance
        passed to :meth:`pyproj.crs.CRS.from_user_input`

        Can be any of:
          - PROJ string
          - Dictionary of PROJ parameters
          - PROJ keyword arguments for parameters
          - JSON string with PROJ parameters
          - CRS WKT string
          - An authority string [i.e. 'epsg:4326']
          - An EPSG integer code [i.e. 4326]
          - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
          - An object with a `to_wkt` method.
          - A :class:`pyproj.crs.CRS` class

        By default, epsg:4269
    """
    if isinstance(feature, str) or isinstance(feature, Path):
        with fiona.open(feature) as src:
            l, b, r, t = src.bounds
        bbox_src_crs = box(*src.bounds)
        shpcrs = get_shapefile_crs(feature)
        if dest_crs is not None and shpcrs != dest_crs:
            bbox_dest_crs = project(bbox_src_crs, shpcrs, dest_crs)
            l, b, r, t = bbox_dest_crs.bounds
        filter = (l, b, r, t)
    elif isinstance(feature, Polygon):
        filter = feature.bounds
    elif isinstance(feature, dict):
        try:
            filter = shape(feature).bounds
        except Exception as ex:
            print(ex)
            print("Supplied dictionary doesn't appear to be valid GeoJSON.")
    return filter


