import os
import collections
import time
from functools import partial
import shutil
import fiona
from shapely.ops import transform, unary_union
from shapely.geometry import shape, mapping, Polygon, box, Point
import pyproj
import numpy as np
import pandas as pd
from fiona.crs import from_epsg, from_string, to_string

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

    @property
    def crs(self):
        if self._crs is None:
            if self._proj_str is not None:
                self._crs = from_string(self._proj_str)
            elif self.epsg is not None:
                self._crs = from_epsg(self.epsg)
        return self._crs

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
            return pyproj.Proj(self.proj_str)

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
        #if other.epsg is not None and other.epsg == self.epsg:
        #    return True
        #if other.proj_str is not None and other.proj_str == self.proj_str:
        #    return True
        #if other.proj_str != self.proj_str:
        #    return False
        return True

    def __repr__(self):
        return self.proj_str


def get_proj_str(prj):
    """Get proj_str string for a projection file

    Parameters
    ----------
    prj : string
        Shapefile or projection file

    Returns
    -------
    proj_str (https://proj.org/)

    """
    '''
    Using fiona (couldn't figure out how to do this with just a prj file)
    from fiona.crs import to_string
    c = fiona.open(shp).crs
    proj_str = to_string(c)
    '''
    # using osgeo
    from osgeo import osr

    prjfile = prj[:-4] + '.prj' # allows shp or prj to be argued
    try:
        with open(prjfile) as src:
            prjtext = src.read()
        srs = osr.SpatialReference()
        srs.ImportFromESRI([prjtext])
        proj_str = srs.ExportToProj4()
        return proj_str
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
    assert np.array_equal(grid.df.node.values, np.arange(grid.size))
    assert np.array_equal(grid.df.node.values, grid.df.index.values)
    polygons = grid.df.loc[rd.node, 'geometry'].values
    if geomtype.lower() == 'polygon':
        rd['geometry'] = polygons
    elif geomtype.lower() == 'point':
        rd['geometry'] = [p.centroid for p in polygons]
    else:
        raise ValueError('Unrecognized geomtype "{}"'.format(geomtype))
    df2shp(rd, filename, epsg=grid.crs.epsg, proj_str=grid.crs.proj_str)


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
        print('\r{}'.format(i+1), end='')
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


def project(geom, projection1, projection2):
    """Reproject a shapely geometry object to new coordinate system

    Parameters
    ----------
    geom: shapely geometry object, list of shapely geometry objects,
          list of (x, y) tuples, or (x, y) tuple.
    projection1: string
        proj string specifying source projection
    projection2: string
        proj string specifying destination projection
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
        proj string specifying source projection
    projection2: string
        proj string specifying destination projection
    """
    projection1 = str(projection1)
    projection2 = str(projection2)

    # define projections
    pr1 = pyproj.Proj(projection1, errcheck=True, preserve_units=True)
    pr2 = pyproj.Proj(projection2, errcheck=True, preserve_units=True)

    return pyproj.transform(pr1, pr2, x, y)


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
            #x, y = project([(l, b),
            #                (r, t)], shpcrs.proj_str, dest_crs.proj_str)
            #filter = (x[0], y[0], x[1], y[1])
            #else:
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
            if filter is not None:
                print('filtering on bounding box {}, {}, {}, {}...'.format(*filter))
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
                continue
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


def shp_properties(df):

    newdtypes = {'bool': 'str',
                 'object': 'str',
                 'datetime64[ns]': 'str'}

    # fiona/OGR doesn't like numpy ints
    # shapefile doesn't support 64 bit ints,
    # but apparently leaving the ints alone is more reliable
    # than intentionally downcasting them to 32 bit
    # pandas is smart enough to figure it out on .to_dict()?
    for c in df.columns:
        if c != 'geometry':
            df[c] = df[c].astype(newdtypes.get(df.dtypes[c].name,
                                               df.dtypes[c].name))
        if 'int' in df.dtypes[c].name:
            if np.max(np.abs(df[c])) > 2**31 -1:
                df[c] = df[c].astype(str)

    # strip dtypes to just 'float', 'int' or 'str'
    def stripandreplace(s):
        return ''.join([i for i in s
                        if not i.isdigit()]).replace('object', 'str')
    dtypes = [stripandreplace(df[c].dtype.name)
              if c != 'geometry'
              else df[c].dtype.name for c in df.columns]
    properties = collections.OrderedDict(list(zip(df.columns, dtypes)))
    return properties


def df2shp(dataframe, shpname, geo_column='geometry', index=False,
           retain_order=False,
           prj=None, epsg=None, proj_str=None, crs=None):
    '''
    Write a DataFrame to a shapefile
    dataframe: dataframe to write to shapefile
    geo_column: optional column containing geometry to write - default is 'geometry'
    index: If true, write out the dataframe index as a column
    retain_order : boolean
        Retain column order in dataframe, using an OrderedDict. Shapefile will
        take about twice as long to write, since OrderedDict output is not
        supported by the pandas DataFrame object.

    --->there are four ways to specify the projection....choose one
    prj: <file>.prj filename (string)
    epsg: EPSG identifier (integer)
    proj_str: pyproj style projection string definition
    crs: crs attribute (dictionary) as read by fiona
    '''
    # first check if output path exists
    if os.path.split(shpname)[0] != '' and not os.path.isdir(os.path.split(shpname)[0]):
        raise IOError("Output folder doesn't exist")

    # check for empty dataframe
    if len(dataframe) == 0:
        raise IndexError("DataFrame is empty!")

    df = dataframe.copy()  # make a copy so the supplied dataframe isn't edited

    # reassign geometry column if geo_column is special (e.g. something other than "geometry")
    if geo_column != 'geometry':
        df['geometry'] = df[geo_column]
        df.drop(geo_column, axis=1, inplace=True)

    # assign none for geometry, to write a dbf file from dataframe
    Type = None
    if 'geometry' not in df.columns:
        df['geometry'] = None
        Type = 'None'
        mapped = [None] * len(df)

    # reset the index to integer index to enforce ordering
    # retain index as attribute field if index=True
    df.reset_index(inplace=True, drop=not index)

    # enforce character limit for names! (otherwise fiona marks it zero)
    # somewhat kludgey, but should work for duplicates up to 99
    df.columns = list(map(str, df.columns))  # convert columns to strings in case some are ints
    overtheline = [(i, '{}{}'.format(c[:8], i)) for i, c in enumerate(df.columns) if len(c) > 10]

    newcolumns = list(df.columns)
    for i, c in overtheline:
        newcolumns[i] = c
    df.columns = newcolumns

    properties = shp_properties(df)
    del properties['geometry']

    # set projection (or use a prj file, which must be copied after shp is written)
    # alternatively, provide a crs in dictionary form as read using fiona
    # from a shapefile like fiona.open(inshpfile).crs

    if epsg is not None:
        from fiona.crs import from_epsg
        crs = from_epsg(int(epsg))
    elif proj_str is not None:
        from fiona.crs import from_string
        crs = from_string(proj_str)
    elif crs is not None:
        pass
    else:
        pass

    if Type != 'None':
        for g in df.geometry:
            try:
                Type = g.type
            except:
                continue
        mapped = [mapping(g) for g in df.geometry]

    schema = {'geometry': Type, 'properties': properties}
    length = len(df)

    if not retain_order:
        props = df.drop('geometry', axis=1).astype(object).to_dict(orient='records')
    else:
        props = [collections.OrderedDict(r) for i, r in df.drop('geometry', axis=1).astype(object).iterrows()]
    print('writing {}...'.format(shpname))
    with fiona.collection(shpname, "w", driver="ESRI Shapefile", crs=crs, schema=schema) as output:
        for i in range(length):
            output.write({'properties': props[i],
                          'geometry': mapped[i]})

    if prj is not None:
        """
        if 'epsg' in prj.lower():
            epsg = int(prj.split(':')[1])
            prjstr = getPRJwkt(epsg).replace('\n', '') # get rid of any EOL
            ofp = open("{}.prj".format(shpname[:-4]), 'w')
            ofp.write(prjstr)
            ofp.close()
        """
        try:
            print('copying {} --> {}...'.format(prj, "{}.prj".format(shpname[:-4])))
            shutil.copyfile(prj, "{}.prj".format(shpname[:-4]))
        except IOError:
            print('Warning: could not find specified prj file. shp will not be projected.')
