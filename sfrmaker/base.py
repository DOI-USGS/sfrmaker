import os
from pathlib import Path
import time
from shapely.geometry import LineString
try:
    import flopy
    fm = flopy.modflow
    mf6 = flopy.mf6
except:
    flopy = False
from gisutils import df2shp
from .gis import export_reach_data
from .grid import StructuredGrid
from .units import get_length_units


class DataPackage:

    package_type = None  # overloaded by SFRData, RivData, etc.

    def __init__(self, grid=None, model=None, isfr=None,
                 model_length_units="undefined", model_time_units='d',
                 package_name=None,
                 **kwargs):
        """
        Base class with shared methods and attributes for model
        package input. Meant to be inherited by SFRData and RivData and
        not called directly.

        Parameters
        ----------
        grid :
        model :
        model_length_units :
        model_time_units :
        package_name :
        kwargs :

        """

        if flopy and isinstance(grid, flopy.discretization.grid.Grid):
            print('\nCreating grid class instance from flopy modelgrid...')
            ta = time.time()
            grid = StructuredGrid.from_modelgrid(grid, isfr=isfr)
            print("grid class created in {:.2f}s\n".format(time.time() - ta))

        # attributes
        self._crs = None
        self._model = None
        self._tables_path = 'tables/'  # default location for writing tables
        self._shapefiles_path = 'shps/'  # default location for writing shapefiles

        # print grid information to screen
        print(grid)
        self.grid = grid
        # units
        self.model_length_units = get_length_units(model_length_units, grid, model)
        self.model_time_units = model_time_units

        self.package_name = package_name

    def write_shapefiles(self, basename=None):
        """Write shapefiles illustrating all aspects of a boundary package.
        """
        if basename is None:
            output_path = Path(self._shapefiles_path)
            output_path.mkdir(exist_ok=True)
            basename = self.package_name
        else:
            output_path, basename = os.path.split(basename)
            basename, _ = os.path.splitext(basename)
        basename = basename.replace(self.package_type, '').strip('_') + '_{}'.format(self.package_type)
        for datatype in 'cells', 'outlets', 'lines', 'routing', 'period_data', 'observations':
            export_method_name = 'export_{}'.format(datatype)
            export_method = getattr(self, export_method_name, None)
            if export_method is None:
                print('{} not supported for '.format(self.__class__))
                continue
            if not callable(export_method):
                export_method = getattr(DataPackage, export_method_name)
            output_shapefile_name = os.path.normpath('{}/{}_{}.shp'.format(output_path, basename, datatype))
            export_method(output_shapefile_name)

        if self.package_type == 'sfr':
            inlets_shapefile = os.path.normpath('{}/{}_sfr_inlets.shp'.format(output_path, basename))
            self.export_transient_variable('flow', inlets_shapefile)

    def export_cells(self, filename=None, nodes=None, geomtype='polygon'):
        """Export shapefile of model cells with stream reaches."""
        if filename is None:
            filename = '{}_{}_cells.shp'.format(self.package_name, self.package_type)
        if self.package_type == 'sfr':
            data = self.reach_data
        else:
            data = self.stress_period_data
        export_reach_data(data, self.grid, filename,
                          nodes=nodes, geomtype=geomtype)

    def export_lines(self, filename=None):
        """Export shapefile of linework"""
        if filename is None:
            filename = '{}_{}_cells.shp'.format(self.package_name, self.package_type)
        if self.package_type == 'sfr':
            data = self.reach_data
        else:
            data = self.stress_period_data
        assert 'geometry' in data.columns and \
               isinstance(data.geometry.values[0], LineString), \
            "No LineStrings in reach_data.geometry"
        df2shp(data, filename, crs=self.grid.crs)

    def export_period_data(self, filename=None, geomtype='point'):
        """Export point shapefile showing locations of period data
        in a MODFLOW-6 SFR package (e.g. inflows, runoff, etc.)

        Parameters
        ----------
        f : str, filename
        geomtype : str ('point' or 'polygon')
            write the locations as points at the cell centers, or polygons
            of the model cells containing the period data.

        """
        if self.package_type != 'sfr':
            return self.export_cells(filename=filename, geomtype=geomtype)

        data = self.period_data.dropna(axis=1).sort_values(by=['per', 'rno'])
        if len(data) == 0:
            print('No period data to export!')
            return

        nodes = dict(zip(self.reach_data.rno, self.reach_data.node))
        for var in ['evaporation', 'inflow', 'rainfall', 'runoff', 'stage']:
            if var in data.columns:
                # pivot the segment data to segments x periods with values of varname
                aggfunc = 'mean'  # how to aggregate multiple instances of rno/per combinations
                if var in ['inflow', 'runoff']:
                    aggfunc = 'sum'
                df = data.reset_index().pivot_table(index='rno', columns='per', values=var,
                                                    aggfunc=aggfunc).reset_index()
                # rename the columns to indicate stress periods
                df.columns = ['rno'] + ['{}{}'.format(i, var) for i in range(df.shape[1] - 1)]
                df['node'] = [nodes[rno] for rno in df['rno']]
                if filename is None:
                    filename = self.package_name + '_{}_period_data_{}.shp'.format(self.package_type,
                                                                                   var)
                elif var not in filename:
                    filename = os.path.splitext(filename)[0] + '_{}.shp'.format(var)
                export_reach_data(df, self.grid, filename, geomtype=geomtype)