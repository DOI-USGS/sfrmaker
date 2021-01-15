import flopy

import sfrmaker


def test_shellmound_simulation(test_data_path, shellmound_simulation):
    assert isinstance(shellmound_simulation, flopy.mf6.MFSimulation)


def test_shellmound_model(shellmound_model):
    assert isinstance(shellmound_model, flopy.mf6.MFModel)


def test_shellmound_grid(shellmound_grid):
    assert isinstance(shellmound_grid, flopy.discretization.StructuredGrid)


def test_lines_from_shapefile(lines_from_shapefile):
    assert isinstance(lines_from_shapefile, sfrmaker.Lines)
