from .gis import crs
from .grid import StructuredGrid, UnstructuredGrid
from .lines import lines
from .sfrdata import sfrdata
from .mf5to6 import mf6sfr
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
