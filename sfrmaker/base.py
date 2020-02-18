import time
import flopy
from .grid import StructuredGrid
from .units import convert_length_units, get_length_units, itmuni_values, lenuni_values


class DataPackage:

    def __init__(self, grid=None, sr=None,
                 model=None, isfr=None,
                 model_length_units="undefined", model_time_units='d',
                 package_name=None,
                 **kwargs):
        """
        Base class with shared methods and attributes for model
        package input (SFRData and RivData).

        Parameters
        ----------
        grid :
        model :
        model_length_units :
        model_time_units :
        package_name :
        kwargs :

        """

        if grid is None and sr is not None:
            print('\nCreating grid class instance from flopy SpatialReference...')
            ta = time.time()
            grid = StructuredGrid.from_sr(sr, isfr=isfr)
            print("grid class created in {:.2f}s\n".format(time.time() - ta))
        elif isinstance(grid, flopy.discretization.grid.Grid):
            print('\nCreating grid class instance from flopy modelgrid...')
            ta = time.time()
            grid = StructuredGrid.from_modelgrid(grid, isfr=isfr)
            print("grid class created in {:.2f}s\n".format(time.time() - ta))

        # attributes
        self._crs = None
        self._model = None

        # print grid information to screen
        print(grid)
        self.grid = grid
        # units
        self.model_length_units = get_length_units(model_length_units, grid, model)
        self.model_time_units = model_time_units

        self.package_name = package_name
