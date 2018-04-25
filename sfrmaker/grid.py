import os
import numpy as np
import pandas as pd
import fiona
from shapely.ops import unary_union
from shapely.geometry import Polygon
import flopy
from .gis import shp2df, get_proj4, crs, read_polygon_feature
fm = flopy.modflow

class grid:
    """Class representing model grid.

    Parameters
    ----------
    model_units : str ('ft' or 'm')
        Computation units units of model
    crs_units : str ('ft' or 'm')
        Units of coordinate reference system for grid.
        (optional; otherwise inferred from epsg, proj4, or prjfile inputs.

    """
    units_dict = {'ft': 1, 'm': 2}

    def __init__(self, df, structured=True,
                 model_units='ft', crs_units=None,
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

        # set the active area where streams will be simulated
        self._bounds = bounds
        self._active_area = None
        self._set_active_area(active_area)

    def __setattr__(self, key, value):
        if key == "active_area":
            self._set_active_area(value)
        else:
            super(grid, self).__setattr__(key, value)

    @property
    def active_area(self):
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
    def lenuni(self):
        return self.units_dict.get(self.model_units, 0)

    def _set_active_area(self, feature=None, ibound_array=None):
        """Establish a polygon that defines the portion of the
        grid where streams will be represented.
        """
        if feature is not None:
            self._active_area = read_polygon_feature(feature, self.crs)

        # include all i, j locations with at least one active cell
        elif ibound_array is not None:
            print('Creating active area from the ibound array... ' +
                  '(make take a while for large grids)')
            if self.structured:
                active = ((ibound_array > 0).sum(axis=0) > 0)
                isactive = active[self.df.i.values,
                                  self.df.j.values]
            else:
                active = ibound_array > 0
                isactive = active[grid.node.values]

            self.df['isactive'] = isactive
            geoms = grid.df.geometry.values[isactive]
            self._active_area = unary_union(geoms)

    def get_node(self, k, i, j):
        return k * self.nrow * self.ncol + i * self.ncol + j

    @staticmethod
    def from_sr(sr=None, active_area=None,
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
        return grid.from_dataframe(df,
                                   bounds=sr.bounds, active_area=active_area,
                                   epsg=epsg, proj4=proj4, prjfile=prjfile)

    @staticmethod
    def from_shapefile(shapefile=None,
                       node_col='node', kcol='k', icol='i', jcol='j',
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
                                   bounds=bounds, active_area=active_area,
                                   epsg=epsg, proj4=proj4, prjfile=prjfile)

    @staticmethod
    def from_dataframe(df=None,
                       node_col='node', kcol='k', icol='i', jcol='j',
                       geometry_column='geometry',
                       model_units='ft', active_area=None,
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

        return grid(df, active_area=active_area,
                    model_units=model_units, structured=structured,
                    epsg=epsg, proj4=proj4, prjfile=prjfile, **kwargs)


