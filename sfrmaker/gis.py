import collections
import os
import time
import traceback
import fiona
import numpy as np
import pyproj
from fiona.crs import from_epsg, from_string, to_string
from shapely.geometry import shape, Polygon, box
from shapely.ops import unary_union
from gisutils import get_proj_str, df2shp, shp2df, project
import sfrmaker


class crs:

    def __init__(self, crs=None,
                 epsg=None, proj_str=None, prjfile=None):

        self._crs = crs
        self.epsg = epsg
        self._proj_str = proj_str
        self._length_units = None
        self.prjfile = prjfile
        if prjfile is not None:
            assert os.path.exists(prjfile), "{} not found.".format(prjfile)
            self._proj_str = get_proj_str(prjfile)
        if not self.is_valid:
            print("Warning, coordinate reference system {}\nis invalid. "
                  "\nPlease supply a valid epsg code, proj4 string, "
                  "or ESRI projection file".format(self.__repr__()))

    @property
    def crs(self):
        if self._crs is None:
            if self._proj_str is not None:
                self._crs = from_string(self._proj_str)
            elif self.epsg is not None:
                self._crs = from_epsg(self.epsg)
        return self._crs

    @property
    def is_valid(self):
        if isinstance(self.pyproj_Proj, pyproj.Proj):
            return True
        else:
            return False

    @property
    def length_units(self):
        if self._length_units is None:
            return parse_units_from_proj_str(self.proj_str)

    @property
    def proj_str(self):
        if self._proj_str is None and self.crs is not None:
            self._proj_str = to_string(self._crs)
        elif self._proj_str is None and self.prjfile is not None:
            self._proj_str = get_proj_str(self.prjfile)
        return self._proj_str

    @property
    def pyproj_Proj(self):
        """pyproj Proj instance. Used for comparing proj4 strings
        that represent the same CRS in different ways."""
        if self.proj_str is not None:
            try:
                return pyproj.Proj(self.proj_str)
            except Exception as ex:
                traceback.print_exception(type(ex), ex, ex.__traceback__)
                return

    def _reset(self):
        self._proj_str = None
        self._crs = None

    def __setattr__(self, key, value):
        if key == "crs":
            super(crs, self).__setattr__("_crs", value)
            self._proj_str = None
        elif key == "proj_str":
            super(crs, self).__setattr__("_proj_str", value)
            self._crs = None
        elif key == "length_units":
            super(crs, self).__setattr__("_length_units", value)
        else:
            super(crs, self).__setattr__(key, value)

    def __eq__(self, other):
        if not isinstance(other, crs):
            return False
        if other.pyproj_Proj is not None and \
                other.pyproj_Proj != self.pyproj_Proj:
            return False
        # if other.epsg is not None and other.epsg == self.epsg:
        #    return True
        # if other.proj_str is not None and other.proj_str == self.proj_str:
        #    return True
        # if other.proj_str != self.proj_str:
        #    return False
        return True

    def __repr__(self):
        return self.proj_str


def build_rtree_index(geom):
    """Builds an rtree index. Useful for multiple intersections with same index.

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
    epsg = grid.crs.epsg
    proj_str = grid.crs.proj_str
    if geomtype.lower() == 'polygon':
        rd['geometry'] = polygons
    elif geomtype.lower() == 'point':
        rd['geometry'] = [p.centroid for p in polygons]
    else:
        raise ValueError('Unrecognized geomtype "{}"'.format(geomtype))
    df2shp(rd, filename, epsg=epsg, proj_str=proj_str)


def intersect_rtree(geom1, geom2, index=None):
    """Intersect features in geom1 with those in geom2. For each feature in geom2, return a list of
     the indices of the intersecting features in geom1.

    Parameters:
    ----------
    geom1 : list
        list of shapely geometry objects
    geom2 : list
        list of shapely polygon objects to be intersected with features in geom1
    index :
        use an index that has already been created

    Returns:
    -------
    A list of the same length as geom2; containing for each feature in geom2,
    a list of indicies of intersecting geometries in geom1.
    """
    if index is None:
        idx = build_rtree_index(geom1)
    else:
        idx = index
    isfr = []
    print('\nIntersecting {} features...'.format(len(geom2)))
    ta = time.time()
    for pind, poly in enumerate(geom2):
        print('\r{}'.format(pind + 1), end='')
        # test for intersection with bounding box of each polygon feature in geom2 using spatial index
        inds = [i for i in idx.intersection(poly.bounds)]
        # test each feature inside the bounding box for intersection with the polygon geometry
        inds = [i for i in inds if geom1[i].intersects(poly)]
        isfr.append(inds)
    print("\nfinished in {:.2f}s".format(time.time() - ta))
    return isfr


def intersect(geom1, geom2):
    """Same as intersect_rtree, except without spatial indexing. Fine for smaller datasets,
    but scales by 10^4 with the side of the problem domain.

    Parameters:
    ----------
    geom1 : list
        list of shapely geometry objects
    geom2 : list
        list of shapely polygon objects to be intersected with features in geom1

    Returns:
    -------
    A list of the same length as geom2; containing for each feature in geom2,
    a list of indicies of intersecting geometries in geom1.
    """

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
    try:
        # need this because preserve_units doesn't seem to be
        # working for complex proj strings.  So if an
        # epsg code was passed, we have no choice, but if a
        # proj string was passed, we can just parse it
        if "EPSG" in proj_str.upper():
            import pyproj

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


def read_polygon_feature(feature, dest_crs, feature_crs=None):
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
    if isinstance(feature, str):
        with fiona.open(feature) as src:
            feature_crs = crs(src.crs)
        geoms = shp2df(feature)['geometry'].values
        feature = unary_union(geoms)
    elif isinstance(feature, collections.Iterable):
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
    if feature_crs is not None and feature_crs != dest_crs:
        feature = project(feature, feature_crs.proj_str, dest_crs.proj_str)
    return feature.buffer(0)


def get_bbox(feature, dest_crs):
    """Get bounding box for a Polygon feature.

    Parameters
    ----------
    feature : str (shapefile path), shapely Polygon or GeoJSON
    dest_crs : proj str
        Desired output coordinate system (shapefiles only)
    """
    if isinstance(feature, str):
        with fiona.open(feature) as src:
            l, b, r, t = src.bounds
            bbox_src_crs = box(*src.bounds)
            shpcrs = crs(proj_str=to_string(src.crs))
        if dest_crs is not None and shpcrs != dest_crs:
            bbox_dest_crs = project(bbox_src_crs, shpcrs.proj_str, dest_crs.proj_str)
            l, b, r, t = bbox_dest_crs.bounds
            # x, y = project([(l, b),
            #                (r, t)], shpcrs.proj_str, dest_crs.proj_str)
            # filter = (x[0], y[0], x[1], y[1])
            # else:
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


