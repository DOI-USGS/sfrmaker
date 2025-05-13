import os
from pathlib import Path
import shutil
import copy
import platform
import numpy as np
import flopy
import flopy.modflow as fm
import pandas as pd
import pytest
import sfrmaker


# stuff for handling slow tests
def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture(scope="session")
def project_root_path():
    filepath = os.path.split(os.path.abspath(__file__))[0]
    return Path(os.path.normpath(os.path.join(filepath, '../../')))


@pytest.fixture(scope="session")
def datapath():
    """Example datasets for users."""
    return Path('examples/')


@pytest.fixture(scope="session")
def test_data_path():
    """Smaller datasets for faster test execution."""
    return Path('sfrmaker/test/data')


@pytest.fixture(scope="session", autouse=True)
def outdir(project_root_path):
    folder = project_root_path / 'sfrmaker/test/temp/'
    reset = True
    if reset:
        shutil.rmtree(folder, ignore_errors=True)
        folder.mkdir()
    return folder


@pytest.fixture(scope="session")
def bin_path(project_root_path):
    bin_path = os.path.join(project_root_path, "bin")
    if "linux" in platform.platform().lower():
        bin_path = os.path.join(bin_path, "linux")
    elif "darwin" in platform.platform().lower():
        bin_path = os.path.join(bin_path, "mac")
    else:
        bin_path = os.path.join(bin_path, "win")
    return Path(bin_path)


@pytest.fixture(scope="session")
def mf6_exe(bin_path):
    return os.path.join(bin_path, 'mf6')


@pytest.fixture(scope="session")
def mfnwt_exe(bin_path):
    return os.path.join(bin_path, 'mfnwt')


@pytest.fixture(scope="function")
def sfr_test_numbering():
    rd = pd.DataFrame()
    rd['i'] = [3, 4, 5,
               7, 8, 9,
               0, 1, 2,
               4, 4, 5,
               0, 0, 0,
               3, 4, 5,
               0, 1, 2,
               4, 5, 6,
               2, 2, 2]
    rd['j'] = [0, 1, 2,
               6, 6, 6,
               6, 6, 6,
               3, 4, 5,
               9, 8, 7,
               6, 6, 6,
               0, 0, 0,
               6, 6, 6,
               9, 8, 7]
    rd['iseg'] = sorted(list(range(1, 10)) * 3)
    rd['ireach'] = [1, 2, 3] * 9

    sd = pd.DataFrame()
    sd['nseg'] = range(1, 10)
    sd['outseg'] = [4, 0, 6, 8, 3, 8, 1, 2, 8]
    return rd, sd


@pytest.fixture(scope="function")
def sfr_testdata(sfr_test_numbering):
    rd, sd = sfr_test_numbering
    sd['width1'] = 1
    sd['width2'] = 1
    return sfrmaker.SFRData(reach_data=rd, segment_data=sd)


@pytest.fixture(scope='module')
def shellmound_simulation(test_data_path, outdir):
    sim_data_folder = '{}/shellmound/shellmound'.format(test_data_path)
    sim_ws = os.path.join(outdir, 'shellmound')
    #shutil.copytree(sim_data_folder, sim_ws, dirs_exist_ok=True)
    # dirs_exist_ok new with python 3.8
    if os.path.exists(sim_ws):
        shutil.rmtree(sim_ws)
    shutil.copytree(sim_data_folder, sim_ws)
    sim = flopy.mf6.MFSimulation.load('mfsim', 'mf6', 'mf6', sim_ws=sim_ws)
    #sim.set_sim_path(sim_ws)
    return sim


@pytest.fixture(scope='module')
def get_shellmound_model(shellmound_simulation, outdir):
    model = shellmound_simulation.get_model('shellmound')
    model.dis.length_units = 'meters'
    return model


@pytest.fixture(scope='function')
def shellmound_model(get_shellmound_model):
    return copy.deepcopy(get_shellmound_model)


@pytest.fixture(scope='module')
def get_shellmound_grid(get_shellmound_model):
    m = copy.deepcopy(get_shellmound_model)
    mg = flopy.discretization.StructuredGrid(delr=m.dis.delr.array,  # cell spacing along a row
                                             delc=m.dis.delc.array,  # cell spacing along a column
                                             xoff=500955.0, yoff=1176285.0,  # lower left corner of model grid
                                             angrot=0,  # grid is unrotated
                                             crs=5070
                                             )
    m._modelgrid = mg
    return mg

@pytest.fixture(scope='function')
def shellmound_grid(get_shellmound_grid):
    return copy.deepcopy(get_shellmound_grid)

@pytest.fixture(scope='function')
def shellmound_sfrmaker_grid(shellmound_grid):
    grid = sfrmaker.StructuredGrid.from_modelgrid(shellmound_grid,
                                                  active_area=None, isfr=None)
    return grid


@pytest.fixture(scope='module')
def tylerforks_active_area_shapefile(datapath):
    return '{}/tylerforks/active_area.shp'.format(datapath)


@pytest.fixture(scope='function')
def tyler_forks_grid_shapefile(datapath, tylerforks_sfrmaker_grid_from_flopy):
    grid = tylerforks_sfrmaker_grid_from_flopy
    shapename = '{}/tylerforks/grid.shp'.format(datapath)
    grid.write_grid_shapefile(shapename)
    return shapename


@pytest.fixture(scope='module')
def get_tylerforks_model(datapath, outdir):
    m = fm.Modflow.load('tf.nam', model_ws='{}/tylerforks/tylerforks'.format(datapath))
    model_ws = os.path.join(outdir, 'tylerforks')
    if not os.path.isdir(model_ws):
        os.makedirs(model_ws)
    m.model_ws = model_ws
    return m


@pytest.fixture(scope='function')
def tylerforks_model(get_tylerforks_model):
    m = copy.deepcopy(get_tylerforks_model)
    return m


@pytest.fixture(scope='module')
def tylerforks_model_grid(get_tylerforks_model):
    m = get_tylerforks_model
    mg = flopy.discretization.StructuredGrid(delr=np.round((m.dis.delr.array * .3048).astype(np.float64), 4),  # cell spacing along a row
                                             delc=np.round((m.dis.delc.array * .3048).astype(np.float64), 4),  # cell spacing along a column
                                             xoff=682688, yoff=5139052,  # lower left corner of model grid
                                             angrot=0,  # grid is unrotated
                                             crs=26715
                                             # projected coordinate system of model (UTM NAD27 zone 15 North)
                                             )
    m.modelgrid = mg
    return copy.deepcopy(mg)


@pytest.fixture(scope='module')
def tylerforks_sfrmaker_grid_from_flopy(tylerforks_model_grid, tylerforks_active_area_shapefile):
    grid = sfrmaker.StructuredGrid.from_modelgrid(mg=tylerforks_model_grid,
                                                  active_area=tylerforks_active_area_shapefile)
    return copy.deepcopy(grid)


@pytest.fixture(scope='module')
def get_tylerforks_lines_from_NHDPlus(datapath):
    pfvaa_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/PlusFlowlineVAA.dbf'.format(datapath)]
    plusflow_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/PlusFlow.dbf'.format(datapath)]
    elevslope_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/elevslope.dbf'.format(datapath)]
    flowlines = ['{}/tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp'.format(datapath)]

    lns = sfrmaker.Lines.from_nhdplus_v2(NHDFlowlines=flowlines,
                                         PlusFlowlineVAA=pfvaa_files,
                                         PlusFlow=plusflow_files,
                                         elevslope=elevslope_files,
                                         bbox_filter='{}/tylerforks/grid.shp'.format(datapath))
    return lns


@pytest.fixture(scope='function')
def tylerforks_lines_from_NHDPlus(get_tylerforks_lines_from_NHDPlus):
    return copy.deepcopy(get_tylerforks_lines_from_NHDPlus)


@pytest.fixture(scope='module')
def get_tylerforks_sfrdata(get_tylerforks_model, get_tylerforks_lines_from_NHDPlus,
                       tylerforks_sfrmaker_grid_from_flopy):
    m = get_tylerforks_model
    tylerforks_lines_from_NHDPlus = copy.deepcopy(get_tylerforks_lines_from_NHDPlus)
    # from the lines and StructuredGrid instances, make a sfrmaker.tylerforks_sfrdata instance
    # (lines are intersected with the model grid and converted to reaches, etc.)
    sfrdata = tylerforks_lines_from_NHDPlus.to_sfr(grid=tylerforks_sfrmaker_grid_from_flopy,
                                                   model=m)
    return sfrdata


@pytest.fixture(scope='function')
def tylerforks_sfrdata(get_tylerforks_sfrdata):
    return copy.deepcopy(get_tylerforks_sfrdata)


@pytest.fixture(scope='module')
def get_lines_from_shapefile(test_data_path):
    flowlines_file = '{}/shellmound/flowlines.shp'.format(test_data_path)
    lns = sfrmaker.Lines.from_shapefile(flowlines_file,
                                        id_column='COMID',
                                        routing_column='tocomid',
                                        arbolate_sum_column2='asum_calc',
                                        asum_units='km',
                                        width1_column='width1',
                                        width2_column='width2',
                                        width_units='feet',
                                        up_elevation_column='elevupsmo',
                                        dn_elevation_column='elevdnsmo',
                                        elevation_units='feet',
                                        name_column='GNIS_NAME',
                                        )
    return lns


@pytest.fixture(scope='function')
def lines_from_shapefile(get_lines_from_shapefile):
    return copy.deepcopy(get_lines_from_shapefile)


@pytest.fixture(scope='module')
def get_shellmound_sfrdata(get_shellmound_model, get_lines_from_shapefile,
                       get_shellmound_grid):
    m = copy.deepcopy(get_shellmound_model)
    lines = copy.deepcopy(get_lines_from_shapefile)
    shellmound_grid = copy.deepcopy(get_shellmound_grid)
    # from the lines and StructuredGrid instances, make a sfrmaker.sfrdata instance
    # (lines are intersected with the model grid and converted to reaches, etc.)
    sfrdata = lines.to_sfr(grid=shellmound_grid,
                                          model=m)
    return sfrdata

@pytest.fixture(scope='function')
def shellmound_sfrdata(get_shellmound_sfrdata):
    return copy.deepcopy(get_shellmound_sfrdata)

@pytest.fixture(params=['shellmound_sfrdata',
                        'tylerforks_sfrdata'])
def sfrdata(request,
            shellmound_sfrdata,
            tylerforks_sfrdata):
    return {'shellmound_sfrdata': shellmound_sfrdata,
            'tylerforks_sfrdata': tylerforks_sfrdata}[request.param]

@pytest.fixture(scope='function')
def neversink_lines_from_nhdplus_hr(datapath):
    NHDPlusHR_paths = [f'{datapath}/neversink_rondout/NHDPLUS_HR_1.gdb', 
                       f'{datapath}/neversink_rondout/NHDPLUS_HR_2.gdb']
    boundary_file = f'{datapath}/neversink_rondout/Model_Extent.shp'
    drop_fcodes = [42803]
    nhdhr_epsg = 4269

    lns = sfrmaker.Lines.from_nhdplus_hr(NHDPlusHR_paths,
                                        bbox_filter=boundary_file,
                                        drop_fcodes=drop_fcodes,
                                        crs=nhdhr_epsg)
    return lns