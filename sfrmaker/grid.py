import os
import numpy as np
import pandas as pd
import fiona
from shapely.ops import unary_union
from shapely.geometry import Polygon
import flopy
from .gis import shp2df, get_proj4, crs, read_polygon_feature, \
    build_rtree_index, intersect_rtree, intersect
fm = flopy.modflow

class grid:
    """Class representing model grid.

    Parameters
    ----------
    model_units : str ('feet' or 'meters')
        Computation units units of model
    crs_units : str ('feet' or 'meters')
        Units of coordinate reference system for grid.
        (optional; otherwise inferred from epsg, proj4, or prjfile inputs.

    """
    units_dict = {'feet': 1, 'meters': 2}

    def __init__(self, df, structured=True,
                 model_units='feet', crs_units=None,
                 bounds=None, active_area=None,
                 epsg=None, proj4=None, prjfile=None):

        self.df = df
        self.structured = structured
        if structured:
            self.nlay = df.k.max() + 1
            self.nrow = df.i.max() + 1
            self.ncol = df.j.max() + 1
        else:
            assert 'node' in df.columns, "Unstructured grids require node numbers."
            self.nlay, self.nrow, self.ncol = None, None, None

        self.model_units = model_units
        self.crs = crs(epsg=epsg, proj4=proj4, prjfile=prjfile)
        if crs_units is not None:
            self.crs._length_units = crs_units

        # spatial index for intersecting
        self._idx = None

        # set the active area where streams will be simulated
        self._bounds = bounds
        self._active_area = None
        self._active_area_defined_by = None
        self._set_active_area(active_area)

    def __setattr__(self, key, value):
        if key == "active_area":
            self._set_active_area(value)
        else:
            super(grid, self).__setattr__(key, value)

    def __repr__(self):
        s = 'Model grid information\n'
        if self.structured:
            s += 'structured grid\nnnodes: {:,d}\n'.format(len(self.df))
            for dim in ['nlay', 'nrow', 'ncol']:
                s += '{}: {:d}\n'.format(dim, self.__dict__[dim])
        else:
            s += 'unstructured grid\n'
            s += 'nnodes: {:,d}\n'.format(len(self.df))
        s += 'model length units: {}\n'.format(self.model_units)
        s += 'crs: {}'.format(self.crs)
        s += 'bounds: {:.2f}, {:.2f}, {:.2f}, {:.2f}\n'.format(*self.bounds)
        s += 'active area defined by: {}'.format(self._active_area_defined_by)
        s += '\n'
        return s

    @property
    def active_area(self):
        """Shapely Polygon delinating area where SFR will be simulated."""
        return self._active_area

    @property
    def bounds(self):
        if self._bounds is None:
            allX = []  # all x coordinates
            allY = []  # all y coordinatess
            geoms = self.df.geometry.tolist() # cell polygons
            for g in geoms:
                gx, gy = g.exterior.coords.xy
                allX += gx
                allY += gy
            allX = np.array(allX)
            allY = np.array(allY)
            self._bounds = allX.min(), allY.min(), allX.max(), allY.max()
        return self._bounds

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
                if isinstance(feature, str):
                    self._active_area_defined_by = feature
                else:
                    self._active_area_defined_by = 'supplied Polygon feature(s)'
            else:
                self._active_area_defined_by = 'isfr array'

        # otherwise if no feature is supplied but all cells are active
        # leave active area as none
        elif self.df.isfr.sum() == len(self.df):
            self._active_area_defined_by = 'all cells'

        # no feature was supplied but some cells are inactive
        else:
            print('Creating active area polygon from unary_union of cells with isfr=1. '
                  'This will take a while for large grids. To avoid this step,'
                  'supply a shapefile or shapely polygon of the SFR domain when'
                  'instantiating the grid objec.')
            geoms = self.df.geometry.values[self.df.isfr == 1]
            self._active_area = unary_union(geoms)
            self._active_area_defined_by = 'isfr array'

    def _set_isfr_from_active_area(self):
        """Intersect model grid cells with active area polygon,
        assign isfr = 1 to cells that intersect."""
        #intersections = intersect_rtree(self.df.geometry.tolist(),
        #                                [self.active_area],
        #                                index=self.spatial_index)
        print('setting isfr values...')
        intersections = intersect(self.df.geometry.tolist(),
                                        [self.active_area])
        self.df.sort_values(by='node', inplace=True)
        self.df['isfr'] = 0
        self.df.loc[np.squeeze(intersections), 'isfr'] = 1

    def get_node(self, k, i, j):
        return k * self.nrow * self.ncol + i * self.ncol + j

    @staticmethod
    def from_sr(sr=None, active_area=None, isfr=None,
                epsg=None, proj4=None, prjfile=None):

        vertices = sr.vertices
        polygons = [Polygon(v) for v in sr.vertices]
        df = pd.DataFrame({'node': range(0, len(vertices)),
                           'i': sorted(list(range(sr.nrow)) * sr.ncol),
                           'j': list(range(sr.ncol)) * sr.nrow,
                           'geometry': polygons
                           }, columns=['node', 'i', 'j', 'geometry'])
        if epsg is None:
            epsg = sr.epsg
        if proj4 is None:
            proj4 = sr.proj4_str
        if prjfile is not None:
            proj4 = get_proj4(prjfile)
        if isfr is not None:
            assert isfr.size == len(df), \
                "isfr column must be of same length as the number of nodes in a layer."
            df['isfr'] = isfr.ravel()
        return grid.from_dataframe(df,
                                   model_units=sr.model_length_units,
                                   bounds=sr.bounds, active_area=active_area,
                                   epsg=epsg, proj4=proj4, prjfile=prjfile)

    @staticmethod
    def from_shapefile(shapefile=None,
                       node_col='node', kcol='k', icol='i', jcol='j',
                       isfr_col='isfr',
                       active_area=None,
                       epsg=None, proj4=None, prjfile=None):

        if prjfile is None:
            prjfile = shapefile.replace('.shp', '.prj')
            prjfile = prjfile if os.path.exists(prjfile) else None
        with fiona.open(shapefile) as src:
            bounds = src.bounds

        df = shp2df(shapefile)
        assert 'geometry' in df.columns, "No feature geometries found in {}.".format(shapefile)

        return grid.from_dataframe(df, node_col=node_col, kcol=kcol, icol=icol, jcol=jcol,
                                   isfr_col=isfr_col,
                                   bounds=bounds, active_area=active_area,
                                   epsg=epsg, proj4=proj4, prjfile=prjfile)

    @staticmethod
    def from_dataframe(df=None,
                       node_col='node', kcol='k', icol='i', jcol='j',
                       isfr_col='isfr',
                       geometry_column='geometry',
                       model_units='feet', active_area=None,
                       epsg=None, proj4=None, prjfile=None, **kwargs):

        assert geometry_column in df.columns, \
            "No feature geometries found in dataframe column '{}'".format(geometry_column)

        structured = False # grids will be treated as unstructured unless row/column information is found
        if icol in df.columns or jcol in df.columns:
            assert icol in df.columns, "No icol={} not found".format(icol)
            assert jcol in df.columns, "No icol={} not found".format(jcol)
            if kcol in df.columns:
                df['k'] = df[kcol]
            else:
                df['k'] = 0
            df['i'] = df[icol]
            df['j'] = df[jcol]
            # convert to zero-based if one-based
            for dim in ['k', 'i', 'j']:
                if df[dim].min() == 1:
                    df[dim] -= 1
            nrow, ncol = df.i.max()+1, df.j.max()+1
            df.sort_values(by=['k', 'i', 'j'], inplace=True)
            df['node'] = ncol * df.i + df.j
            # only retain first instance of each node
            # (subsequent instances would be underlying layers for same cell)
            # (this should also work for grids that are missing cells in a given layer;
            # as long as the cell is in at least one layer)
            if len(set(df.node).difference(df.groupby('k').get_group(0).node)) != 0:
                # this may take awhile for large grids
                df = df.groupby('node').first()
                df['node'] = df.index # put node back in columns
            structured = True

        elif node_col in df.columns:
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

        return grid(df, active_area=active_area,
                    model_units=model_units, structured=structured,
                    epsg=epsg, proj4=proj4, prjfile=prjfile, **kwargs)


