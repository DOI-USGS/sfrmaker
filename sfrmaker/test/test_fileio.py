from flopy.discretization import StructuredGrid
from sfrmaker.fileio import load_modelgrid


def test_load_grid():
    gridfile = 'sfrmaker/test/data/shellmound/shellmound/shellmound_grid.json'
    modelgrid = load_modelgrid(gridfile)
    assert isinstance(modelgrid, StructuredGrid)