"""
Module for working with model grids. For examples of how to easily set
up a StructuredGrid instance, see :ref:`Basic usage of SFRmaker in a scripting context`.
"""

from pathlib import Path
import warnings

import fiona
import flopy
import numpy as np
import pandas as pd
from rasterio import Affine
from rasterio import features
from shapely.geometry import Polygon, shape, box
from shapely.ops import unary_union
from gisutils import shp2df, df2shp, get_shapefile_crs
from .gis import get_crs, read_polygon_feature, \
    build_rtree_index, intersect

fm = flopy.modflow


class Grid:
    """Base class for model grids. Has methods and attributes
    that are common to both Structured and Unstructured Grids. Not
    meant to be called directly.

    """
    units_dict = {0: 'undefined', 'feet': 1, 'meters': 2}

    def __init__(self, df,
                 model_units='undefined', crs_units=None,
                 bounds=None, crs=None, prjfile=None, **kwargs):

        self.df = df

        # coordinate projection stuff
        self.model_units = model_units
        self.crs = get_crs(crs=crs, prjfile=prjfile, **kwargs)
        if crs_units is not None:
            warnings.warn('The crs_units argument is deprecated. Units are now read '
                          'from the crs information via the pyproj.crs.CRS class')

        # spatial index for intersecting
        self._idx = None

        # set the active area where streams will be simulated
        self._bounds = bounds
        self._active_area = None
        self._active_area_defined_by = None

    def __setattr__(self, key, value):
        if key == "active_area":
            self._set_active_area(value)
        else:
            super(Grid, self).__setattr__(key, value)

    def __repr__(self):
        s = 'Model grid information\n'
        if isinstance(self, StructuredGrid):
            s += 'structured grid\nnnodes: {:,d}\n'.format(len(self.df))
            for dim in ['nlay', 'nrow', 'ncol']:
                s += '{}: {:d}\n'.format(dim, self.__dict__[dim])
        else:
            s += 'unstructured grid\n'
            s += 'nnodes: {:,d}\n'.format(len(self.df))
        s += 'model length units: {}\n'.format(self.model_units)
        s += 'crs: {}\n'.format(self.crs)
        s += 'bounds: {:.2f}, {:.2f}, {:.2f}, {:.2f}\n'.format(*self.bounds)
        s += 'active area defined by: {}'.format(self._active_area_defined_by)
        s += '\n'
        return s

    def __eq__(self, other):
        if not isinstance(other, Grid):
            print('not an sfrmaker.Grid instance')
            return False
        if other._structured != self._structured:
            print('different grid types')
            return False
        if other.size != self.size:
            print('grid sizes not equal!')
            return False
        if other.crs != self.crs:
            print('crs {} is not equal to {}!'.format(other.crs, self.crs))
            return False
        if not np.allclose(other.bounds, self.bounds):
            return False
        if self._structured:
            if other.nrow != self.nrow:
                return False
            if other.ncol != self.ncol:
                return False
            if other.rotation != self.rotation:
                return False
        if not np.array_equal(self.isfr, other.isfr):
            print('idomain arrays are not equal!')
            return False
        return True

    @property
    def active_area(self):
        """Shapely Polygon delinating area where SFR will be simulated."""
        return self._active_area

    @property
    def bounds(self):
        if self._bounds is None:
            allX = []  # all x coordinates
            allY = []  # all y coordinatess
            geoms = self.df.geometry.tolist()  # cell polygons
            for g in geoms:
                gx, gy = g.exterior.coords.xy
                allX += gx
                allY += gy
            allX = np.array(allX)
            allY = np.array(allY)
            self._bounds = allX.min(), allY.min(), allX.max(), allY.max()
        return self._bounds

    @property
    def size(self):
        return len(self.df)

    @property
    def spatial_index(self):
        """Rtree index for intersecting features with model grid."""
        if self._idx is None:
            self._idx = build_rtree_index(self.df.geometry.tolist())
        return self._idx

    @property
    def lenuni(self):
        return self.units_dict.get(self.model_units, 0)

    def _set_active_area(self, feature=None):
        """Establish a polygon that defines the portion of the
        grid where streams will be represented.

        feature : shapely Polygon, list of Polygons, or shapefile path
            Polygons must be in same CRS as linework; shapefile
            features will be reprojected if their crs is different.
        isfr : list or ndarray of boolean values
            Length must be equal to nrow * ncol, or the number of nodes in a layer (
            Indicates whether or not a particular cell can have an SFR reach.
            (0 or False indicates no SFR).
        """
        if feature is not None:
            self._active_area = read_polygon_feature(feature, self.crs)

            # if a polygon feature is supplied but all cells are active
            # set isfr from polygon
            if self.df.isfr.sum() == len(self.df):
                self._set_isfr_from_active_area()
                if isinstance(feature, str) or isinstance(feature, Path):
                    self._active_area_defined_by = feature
                else:
                    self._active_area_defined_by = 'supplied Polygon feature(s)'
            else:
                self._active_area_defined_by = 'isfr array'

        # otherwise if no feature is supplied but all cells are active
        # leave active area as none
        # main reason not to handle this the same as below
        # is that unary_unions of a lot of cells are slow
        # and with all cells active,
        # don't have to worry about inactive cells getting intersected
        elif self.df.isfr.sum() == len(self.df):
            self._active_area_defined_by = 'all cells'

        # no feature was supplied; set from isfr values
        else:
            self.create_active_area_polygon_from_isfr()
            self._active_area_defined_by = 'isfr array'

    def _set_isfr_from_active_area(self):
        """Intersect model grid cells with active area polygon,
        assign isfr = 1 to cells that intersect."""
        # intersections = intersect_rtree(self.df.geometry.tolist(),
        #                                [self.active_area],
        #                                index=self.spatial_index)
        print('setting isfr values...')
        intersections = intersect(self.df.geometry.tolist(),
                                  [self.active_area])
        self.df.sort_values(by='node', inplace=True)
        self.df['isfr'] = 0
        self.df.loc[np.squeeze(intersections), 'isfr'] = 1

    def create_active_area_polygon_from_isfr(self):
        """The StructuredGrid and UnstructuredGrid classes
        have their own ways of doing this."""
        return

    def get_node(self, k, i, j):
        return k * self.nrow * self.ncol + i * self.ncol + j

    def write_active_area_shapefile(self, outshp='active_area.shp'):
        if self._active_area is None:
            self.create_active_area_polygon_from_isfr()
        assert isinstance(self._active_area, Polygon), \
            "active area didn't get set correctly (not a shapely Polygon)"
        df = pd.DataFrame({'geometry': [self._active_area],
                           'description': ['Active area where SFR will be applied.']})
        df2shp(df, outshp, crs=self.crs)

    def write_grid_shapefile(self, outshp='grid.shp'):
        df2shp(self.df, outshp, crs=self.crs)


class StructuredGrid(Grid):
    """Class representing a model grid that has a row/column structure.


    Parameters
    ----------
    df : DataFrame
        Pandas DataFrame that is the primary container for information about the model grid.
        Must have the following columns:

        ========= ===  ==================================================
        k         int  model layer (zero-based)
        i         int  model row (zero-based)
        j         int  model column (zero-based)
        isfr      int  flag indicating whether the cell can have SFR reaches
                       (0=False, 1=True)
        geometry  obj  shapely :class:`Polygons <Polygon>` of model cells
        ========= ===  ==================================================

    xul : float, optional
        Upper left corner of the grid x-coordinate. Only used for creating
        the :attr:`transform` attribute, by default None
    yul : [type], optional
        Upper left corner of the grid y-coordinate.  Only used for creating
        the :attr:`transform` attribute, by default None
    dx : float, optional
        Uniform spacing in the x-direction (if the grid is uniform), 
        Only used for creating
        the :attr:`transform` attribute, by default None
    dy : float, optional
        Uniform spacing in the x-direction (if the grid is uniform), 
        Only used for creating
        the :attr:`transform` attribute, by default None
    rotation : float, optional
        Grid rotation angle in degrees, counter-clockwise
        about the origin, by default 0. Only used for creating
        the :attr:`transform` attribute, by default None
    uniform : bool, optional
        Optional flag indicating the grid is uniform, 
        by default None
    model_units : str, optional, {'meters', 'feet', ..}
        Model length units, by default 'undefined'
    crs_units : str, optional, {'meters', 'feet', ..}
        Coordinate reference system length. Usually these
        are read from the CRS information below,
        by default None
    bounds : tuple, optional
        (left, bottom, top, right) edges of the grid bounding box, if known.
        Otherwise, the :attr:`bounds` attribute is computed from the
        shapely :class:`Polygons <Polygon>` in the :attr:`Grid DataFrame <Grid.df>` attribute.
        by default None
    active_area : shapely Polygon, list of Polygons, or shapefile path, optional
        Polygon defining the active portion of the model grid.
        Polygons must be in same CRS as linework; shapefile features will be reprojected if their crs is different.
        by default None, in which case the entire grid is assumed to be active.
    crs : obj, optional
        Coordinate reference object for ``df``.
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
    prjfile: str, optional
        ESRI-style projection file with coordinate reference information for ``df``. 
    """
    _structured = True

    def __init__(self, df,
                 xul=None, yul=None, dx=None, dy=None, rotation=0.,
                 uniform=None,
                 model_units='undefined', crs_units=None,
                 bounds=None, active_area=None,
                 crs=None, prjfile=None, **kwargs):

        Grid.__init__(self, df, model_units=model_units, crs_units=crs_units,
                      bounds=bounds, active_area=active_area,
                      crs=crs, prjfile=prjfile, **kwargs)

        # structured grid parameters
        self.xul = xul
        self.yul = yul
        self.rotation = rotation
        self.nrow = self.df.i.max() + 1
        self.ncol = self.df.j.max() + 1

        # uniform structured grid parameters
        self._uniform = uniform  # whether grid is uniform or not
        self.dx = dx
        self.dy = dy

        self.nlay = df.k.max() + 1
        self.nrow = df.i.max() + 1
        self.ncol = df.j.max() + 1

        self._set_active_area(active_area)

    @property
    def isfr(self):
        return np.reshape(self.df.isfr.values,
                          (self.nrow, self.ncol)).astype(int)

    @property
    def uniform(self):
        """Check if cells are uniform by comparing their areas."""
        if self._uniform is None:
            areas = [g.area for g in self.df.geometry]
            self._uniform = np.allclose(areas, np.mean(areas), rtol=0.01)
        return self._uniform

    @property
    def transform(self):
        """Rasterio-style affine transform object.
        https://www.perrygeo.com/python-affine-transforms.html
        """
        if self.uniform:
            for param in ['dx', 'rotation', 'xul', 'dy', 'yul']:
                if self.__dict__[param] is None:
                    print('This method requires a uniform grid and '
                          'specification of xul, yul, dx, dy, and rotation.')
                    return
            return Affine(self.dx, 0., self.xul,
                          0., -self.dy, self.yul) * Affine.rotation(-self.rotation)

    def create_active_area_polygon_from_isfr(self):
        """Convert 2D numpy array representing active area where
        SFR will be simulated (isfr) to a polygon (if multiple
        polygons area created, the largest one by area is retained).
        """
        if self.transform is not None:
            # vectorize the raster
            shapes = features.shapes(self.isfr.astype(np.int32), transform=self.transform)
            # convert the shapes corresponding to raster values of 1 to shapely objects
            shapes = list(shapes)
            shapes = [shape(s[0]) for s in list(shapes) if s[1] == 1]
            # get the shape with the largest area
            areas = [s.area for s in shapes]
            self._active_area = shapes[np.argmax(areas)]
        else:
            self._active_area = unary_union(self.df.loc[self.df.isfr == 1, 'geometry'])

    @classmethod
    def from_json(cls, jsonfile, active_area=None, isfr=None,
                  crs=None, prjfile=None, **kwargs):
        from sfrmaker.fileio import load_modelgrid
        grid = load_modelgrid(jsonfile)
        return cls.from_modelgrid(grid, active_area=active_area, isfr=isfr,
                                  crs=crs, prjfile=prjfile, **kwargs)

    @classmethod
    def from_sr(cls, **kwargs):
        raise AttributeError("the from_sr method is depricated. Please use from_modelgrid,"
                             "as demonstrated in the examples.")

    @classmethod
    def from_modelgrid(cls, mg=None, active_area=None, isfr=None,
                       crs=None, prjfile=None):
        """Create StructureGrid class instance from a
        flopy.discretization.StructuredGrid instance."""
        i, j = np.indices((mg.nrow, mg.ncol))
        vertices = mg._cell_vert_list(i.ravel(), j.ravel())
        polygons = [Polygon(v) for v in vertices]
        df = pd.DataFrame({'node': range(0, len(vertices)),
                           'i': sorted(list(range(mg.nrow)) * mg.ncol),
                           'j': list(range(mg.ncol)) * mg.nrow,
                           'geometry': polygons
                           }, columns=['node', 'i', 'j', 'geometry'])

        crs = get_crs(prjfile=prjfile, crs=mg.crs)

        if isfr is not None:
            # if a 3D array is supplied for isfr, convert to 2D
            # (retain all i, j locations with at least one active layer;
            # assuming that top of each highest-active cell represents the land surface)
            if len(isfr.shape) == 3:
                isfr = np.any(isfr == 1, axis=0).astype(int)
                df['isfr'] = isfr.ravel()
            elif isfr.shape == (mg.nrow, mg.ncol):
                df['isfr'] = isfr.ravel()
            else:
                assert isfr.size == len(df), \
                    "isfr must be of shape (nlay, nrow, ncol), (nrow, ncol) or (nrow * ncol,)"
                df['isfr'] = isfr
        uniform = False
        dx, dy = None, None
        if len(set(mg.delc)) == 1 and len(set(mg.delr)) == 1:
            dx = mg.delr[0]
            dy = mg.delc[0]
            uniform = True
        bounds = mg.extent[0], mg.extent[2], mg.extent[1], mg.extent[3]

        # upper left corner
        x0 = mg.xyedges[0][0]
        y0 = mg.xyedges[1][0]
        xul, yul = mg.get_coords(x0, y0)

        return cls.from_dataframe(df, uniform=uniform,
                                  xul=xul, yul=yul, dx=dx, dy=dy,
                                  rotation=mg.angrot,
                                  bounds=bounds, active_area=active_area,
                                  crs=crs)

    @classmethod
    def from_shapefile(cls, shapefile=None,
                       node_col='node', kcol='k', icol='i', jcol='j',
                       isfr_col='isfr',
                       active_area=None,
                       crs=None, prjfile=None, **kwargs):

        if crs is None:
            crs = get_shapefile_crs(shapefile)
        crs = get_crs(prjfile=prjfile, crs=crs, **kwargs)

        with fiona.open(shapefile) as src:
            bounds = src.bounds

        df = shp2df(shapefile)
        assert 'geometry' in df.columns, "No feature geometries found in {}.".format(shapefile)

        return cls.from_dataframe(df, node_col=node_col, kcol=kcol, icol=icol, jcol=jcol,
                                  isfr_col=isfr_col,
                                  bounds=bounds, active_area=active_area,
                                  crs=crs)

    @classmethod
    def from_dataframe(cls, df=None, uniform=None,
                       kcol='k', icol='i', jcol='j',
                       isfr_col='isfr',
                       geometry_column='geometry',
                       active_area=None,
                       crs=None, prjfile=None, **kwargs):

        assert geometry_column in df.columns, \
            "No feature geometries found in dataframe column '{}'".format(geometry_column)

        assert icol in df.columns, "No icol='{}' not found".format(icol)
        assert jcol in df.columns, "No jcol='{}' not found".format(jcol)

        crs = get_crs(prjfile=prjfile, crs=crs, **kwargs)

        # set layer column
        if kcol in df.columns:
            df['k'] = df[kcol]
        else:
            df['k'] = 0

        # set i, j columns
        df['i'] = df[icol]
        df['j'] = df[jcol]
        # convert to zero-based if one-based
        for dim in ['k', 'i', 'j']:
            if df[dim].min() == 1:
                df[dim] -= 1
        nrow, ncol = df.i.max() + 1, df.j.max() + 1
        # add missing cells if the grid is incomplete
        if len(df) < (nrow * ncol):
            df.index = ncol * df.i + df.j
            df = df.reindex(np.arange(nrow * ncol))
            i, j = np.indices((nrow, ncol))
            df['k'] = 0
            df['i'] = i.ravel()
            df['j'] = j.ravel()
            missing_cells = df['geometry'].isnull()
            df.loc[missing_cells, 'geometry'] = box(0, 0, 0, 0)
            if isfr_col in df.columns:
                df.loc[missing_cells, isfr_col] = 0
            else:
                df['isfr'] = 1
                df.loc[missing_cells, 'isfr'] = 0
        df.sort_values(by=['k', 'i', 'j'], inplace=True)
        df['node'] = ncol * df.i + df.j

        # only retain first instance of each node
        # (subsequent instances would be underlying layers for same cell)
        # (this should also work for grids that are missing cells in a given layer;
        # as long as the cell is in at least one layer)
        if len(set(df.node).difference(df.groupby('k').get_group(0).node)) != 0:
            # this may take awhile for large grids
            df = df.groupby('node').first()
            df['node'] = df.index  # put node back in columns

        if isfr_col in df.columns:
            df['isfr'] = df[isfr_col].astype(int)
        else:
            df['isfr'] = 1

        return cls(df, active_area=active_area,
                   uniform=uniform,
                   crs=crs, **kwargs)


class UnstructuredGrid(Grid):
    """Class representing an unstructured model grid."""
    _structured = False

    def __init__(self, df,
                 model_units='undefined', crs_units=None,
                 bounds=None, active_area=None,
                 crs=None, prjfile=None, **kwargs):
        Grid.__init__(self, df, model_units=model_units, crs_units=crs_units,
                      bounds=bounds, active_area=active_area,
                      crs=crs, prjfile=prjfile, **kwargs)

        assert 'node' in df.columns, \
            "DataFrame df must have a 'node' column for identifying model cells."

        self._set_active_area(active_area)

    def create_active_area_polygon_from_isfr(self):
        """Create active area polygon from union of cells where isfr=1.
        """
        print('Creating active area polygon from shapely.ops.unary_union '
              'of cells with isfr=1. '
              'This will take a while for large grids. To avoid this step,'
              'supply a shapefile or shapely polygon of the SFR domain when'
              'instantiating the grid objec.')
        geoms = self.df.geometry.values[self.df.isfr == 1]
        self._active_area = unary_union(geoms)

    @classmethod
    def from_dataframe(cls, df=None,
                       node_col='node',
                       isfr_col='isfr',
                       geometry_column='geometry',
                       model_units='feet', active_area=None,
                       crs=None, prjfile=None, **kwargs):

        assert geometry_column in df.columns, \
            "No feature geometries found in dataframe column '{}'".format(geometry_column)

        crs = get_crs(prjfile=prjfile, crs=crs, **kwargs)

        if node_col in df.columns:
            df['node'] = df[node_col]
            # convert to zero-based if one-based
            if df.node.min() == 1:
                df['node'] -= 1

            df.sort_values(by=['node'], inplace=True)
            # enforce consecutive node numbering
            df['node'] = np.arange(len(df))
        else:
            df['node'] = np.arange(len(df))

        if isfr_col in df.columns:
            df['isfr'] = df[isfr_col].astype(int)
        else:
            df['isfr'] = 1
        return cls(df, active_area=active_area,
                   model_units=model_units,
                   crs=crs, **kwargs)
