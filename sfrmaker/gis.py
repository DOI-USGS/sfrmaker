import os
import collections
import time
from functools import partial
import fiona
from shapely.ops import transform, unary_union
from shapely.geometry import shape, mapping
import pyproj
import numpy as np
import pandas as pd
from fiona.crs import from_epsg, from_string, to_string

class crs:

    def __init__(self, crs=None,
                 epsg=None, proj4=None, prjfile=None):

        self._crs = crs
        self.epsg = epsg
        self._proj4 = proj4
        self._length_units = None
        self.prjfile = prjfile
        if prjfile is not None:
            assert os.path.exists(prjfile), "{} not found.".format(prjfile)
            self._proj4 = get_proj4(prjfile)

    @property
    def crs(self):
        if self._crs is None:
            if self._proj4 is not None:
                self._crs = from_string(self._proj4)
            elif self.epsg is not None:
                self._crs = from_epsg(self.epsg)
        return self._crs

    @property
    def length_units(self):
        if self._length_units is None:
            return parse_units_from_proj4(self.proj4)

    @property
    def proj4(self):
        if self._proj4 is None and self._crs is not None:
            self._proj4 = to_string(self._crs)
        elif self._proj4 is None and self.prjfile is not None:
            self._proj4 = get_proj4(self.prjfile)
        return self._proj4

    def _reset(self):
        self._proj4 = None
        self._crs = None

    def __setattr__(self, key, value):
        if key == "crs":
            super(crs, self).__setattr__("_crs", value)
            self._proj4 = None
        elif key == "proj4":
            super(crs, self).__setattr__("_proj4", value)
            self._crs = None
        elif key == "length_units":
            super(crs, self).__setattr__("_length_units", value)
        else:
            super(crs, self).__setattr__(key, value)

    def __eq__(self, other):
        if not isinstance(other, crs):
            return False
        if other.proj4 != self.proj4:
            return False
        return True

def get_proj4(prj):
    """Get proj4 string for a projection file

    Parameters
    ----------
    prj : string
        Shapefile or projection file

    Returns
    -------
    proj4 string (http://trac.osgeo.org/proj/)

    """
    '''
    Using fiona (couldn't figure out how to do this with just a prj file)
    from fiona.crs import to_string
    c = fiona.open(shp).crs
    proj4 = to_string(c)
    '''
    # using osgeo
    from osgeo import osr

    prjfile = prj[:-4] + '.prj' # allows shp or prj to be argued
    try:
        with open(prjfile) as src:
            prjtext = src.read()
        srs = osr.SpatialReference()
        srs.ImportFromESRI([prjtext])
        proj4 = srs.ExportToProj4()
        return proj4
    except:
        pass

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

def intersect_rtree(geom1, geom2):
    """Intersect features in geom1 with those in geom2. For each feature in geom2, return a list of
     the indices of the intersecting features in geom1.

    Parameters:
    ----------
    geom1 : list or rtree spatial index object
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
    if isinstance(geom1, list):
        idx = build_rtree_index(geom1)
    else:
        idx = geom1
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
    print("\nfinished in {:.2f}s\n".format(time.time() - ta))
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
        print('\r{}'.format(i+1), end='')
        intersects = np.array([r.intersects(g) for r in geom1])
        inds = list(np.arange(ngeom1)[intersects])
        isfr.append(inds)
    print("\nfinished in {:.2f}s\n".format(time.time() - ta))
    return isfr

def parse_units_from_proj4(proj4_str):
    units = None
    try:
        # need this because preserve_units doesn't seem to be
        # working for complex proj4 strings.  So if an
        # epsg code was passed, we have no choice, but if a
        # proj4 string was passed, we can just parse it
        if "EPSG" in proj4_str.upper():
            import pyproj

            crs = pyproj.Proj(proj4_str,
                              preseve_units=True,
                              errcheck=True)
            proj_str = crs.srs
        else:
            proj_str = proj4_str
        # http://proj4.org/parameters.html#units
        # from proj4 source code
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

def project(geom, projection1, projection2):
    """Reproject a shapely geometry object to new coordinate system

    Parameters
    ----------
    geom: shapely geometry object, list of shapely geometry objects,
          list of (x, y) tuples, or (x, y) tuple.
    projection1: string
        Proj4 string specifying source projection
    projection2: string
        Proj4 string specifying destination projection
    """
    # check for x, y values instead of shapely objects
    if isinstance(geom, tuple):
        return np.squeeze([projectXY(geom[0], geom[1], projection1, projection2)])

    if isinstance(geom, collections.Iterable):
        geom = list(geom) # in case it's a generator
        geom0 = geom[0]
    else:
        geom0 = geom

    if isinstance(geom0, tuple):
        a = np.array(geom)
        x = a[:, 0]
        y = a[:, 1]
        return np.squeeze([projectXY(x, y, projection1, projection2)])

    # transform shapely objects
    # enforce strings
    projection1 = str(projection1)
    projection2 = str(projection2)

    # define projections
    pr1 = pyproj.Proj(projection1, errcheck=True, preserve_units=True)
    pr2 = pyproj.Proj(projection2, errcheck=True, preserve_units=True)

    # projection function
    # (see http://toblerity.org/shapely/shapely.html#module-shapely.ops)
    project = partial(pyproj.transform, pr1, pr2)

    # do the transformation!
    if isinstance(geom, collections.Iterable):
        return [transform(project, g) for g in geom]
    return transform(project, geom)


def projectXY(x, y, projection1, projection2):
    """Project x and y coordinates to different crs

    Parameters
    ----------
    x: scalar or 1-D array
    x: scalar or 1-D array
    projection1: string
        Proj4 string specifying source projection
    projection2: string
        Proj4 string specifying destination projection
    """
    projection1 = str(projection1)
    projection2 = str(projection2)

    # define projections
    pr1 = pyproj.Proj(projection1, errcheck=True, preserve_units=True)
    pr2 = pyproj.Proj(projection2, errcheck=True, preserve_units=True)

    return pyproj.transform(pr1, pr2, x, y)

def read_feature(feature, dest_crs):
    """Read a geometric feature from a shapefile, shapely geometry object,
    or collection of shapely geometry objects. Reproject to dest_crs
    if the feature is read from a shapefile in a different CRS.

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
        geoms = shp2df(feature)['geomety'].values
        feature = unary_union(geoms)
        if feature_crs != dest_crs:
            feature = project(feature, feature_crs.proj4, dest_crs.proj4)
    elif isinstance(feature, collections.Iterable):
        feature = unary_union(feature)
    return feature

def shp2df(shplist, index=None, index_dtype=None, clipto=[], filter=None,
           true_values=None, false_values=None, layer=None,
           skip_empty_geom=True):
    """Read shapefile/DBF, list of shapefiles/DBFs, or File geodatabase (GDB)
     into pandas DataFrame.

    Parameters
    ----------
    shplist : string or list
        of shapefile/DBF name(s) or FileGDB
    index : string
        Column to use as index for dataframe
    index_dtype : dtype
        Enforces a datatype for the index column (for example, if the index field is supposed to be integer
        but pandas reads it as strings, converts to integer)
    clipto : list
        limit what is brought in to items in index of clipto (requires index)
    filter : tuple (xmin, ymin, xmax, ymax)
        bounding box to filter which records are read from the shapefile.
    true_values : list
        same as argument for pandas read_csv
    false_values : list
        same as argument for pandas read_csv
    layer : str
        Layer name to read (if opening FileGDB)
    skip_empty_geom : True/False, default True
        Drops shapefile entries with null geometries.
        DBF files (which specify null geometries in their schema) will still be read.

    Returns
    -------
    df : DataFrame
        with attribute fields as columns; feature geometries are stored as
    shapely geometry objects in the 'geometry' column.
    """
    if isinstance(shplist, str):
        shplist = [shplist]
    if not isinstance(true_values, list) and true_values is not None:
        true_values = [true_values]
    if not isinstance(false_values, list) and false_values is not None:
        false_values = [false_values]
    if len(clipto) > 0 and index:
        clip = True
    else:
        clip = False

    df = pd.DataFrame()
    for shp in shplist:
        print("\nreading {}...".format(shp))
        shp_obj = fiona.open(shp, 'r', layer=layer)

        if index is not None:
            # handle capitolization issues with index field name
            fields = list(shp_obj.schema['properties'].keys())
            index = [f for f in fields if index.lower() == f.lower()][0]

        attributes = []
        # for reading in shapefiles
        meta = shp_obj.meta
        if meta['schema']['geometry'] != 'None':

            if clip:  # limit what is brought in to items in index of clipto
                for line in shp_obj.filter(bbox=filter):
                    props = line['properties']
                    if not props[index] in clipto:
                        continue
                    props['geometry'] = line.get('geometry', None)
                    attributes.append(props)
            else:
                for line in shp_obj.filter(bbox=filter):
                    props = line['properties']
                    props['geometry'] = line.get('geometry', None)
                    attributes.append(props)
            print('--> building dataframe... (may take a while for large shapefiles)')
            shp_df = pd.DataFrame(attributes)
            # reorder fields in the DataFrame to match the input shapefile
            if len(attributes) > 0:
                shp_df = shp_df[list(attributes[0].keys())]

            # handle null geometries
            if len(shp_df) == 0:
                print('Empty dataframe! No features were read.')
                if filter is not None:
                    print('Check filter {} for consistency \
with shapefile coordinate system'.format(filter))
            geoms = shp_df.geometry.tolist()
            if geoms.count(None) == 0:
                shp_df['geometry'] = [shape(g) for g in geoms]
            elif skip_empty_geom:
                null_geoms = [i for i, g in enumerate(geoms) if g is None]
                shp_df.drop(null_geoms, axis=0, inplace=True)
                shp_df['geometry'] = [shape(g) for g in shp_df.geometry.tolist()]
            else:
                shp_df['geometry'] = [shape(g) if g is not None else None
                                      for g in geoms]

        # for reading in DBF files (just like shps, but without geometry)
        else:
            if clip:  # limit what is brought in to items in index of clipto
                for line in shp_obj:
                    props = line['properties']
                    if not props[index] in clipto:
                        continue
                    attributes.append(props)
            else:
                for line in shp_obj:
                    attributes.append(line['properties'])
            print('--> building dataframe... (may take a while for large shapefiles)')
            shp_df = pd.DataFrame(attributes)
            # reorder fields in the DataFrame to match the input shapefile
            if len(attributes) > 0:
                shp_df = shp_df[list(attributes[0].keys())]

        shp_obj.close()
        if len(shp_df) == 0:
            continue
        # set the dataframe index from the index column
        if index is not None:
            if index_dtype is not None:
                shp_df[index] = shp_df[index].astype(index_dtype)
            shp_df.index = shp_df[index].values

        df = df.append(shp_df)

        # convert any t/f columns to numpy boolean data
        if true_values is not None or false_values is not None:
            replace_boolean = {}
            for t in true_values:
                replace_boolean[t] = True
            for f in false_values:
                replace_boolean[f] = False

            # only remap columns that have values to be replaced
            cols = [c for c in df.columns if c != 'geometry']
            for c in cols:
                if len(set(replace_boolean.keys()).intersection(set(df[c]))) > 0:
                    df[c] = df[c].map(replace_boolean)

    return df