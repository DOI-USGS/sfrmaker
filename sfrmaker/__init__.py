from ._version import get_versions
from .gis import CRS
from .grid import StructuredGrid, UnstructuredGrid
from .lines import Lines
from .mf5to6 import Mf6SFR
from .rivdata import RivData
from .sfrdata import SFRData

__version__ = get_versions()['version']
del get_versions

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
