import warnings
warnings.filterwarnings("ignore")

from ._version import get_versions
from sfrmaker.grid import StructuredGrid, UnstructuredGrid
from sfrmaker.lines import Lines
from sfrmaker.mf5to6 import Mf6SFR
from sfrmaker.rivdata import RivData
from sfrmaker.sfrdata import SFRData

__version__ = get_versions()['version']
del get_versions

from . import _version
__version__ = _version.get_versions()['version']
