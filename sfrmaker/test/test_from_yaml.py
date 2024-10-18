"""Test configuration file input to Sfrmaker.
"""
import os
import yaml
import shutil
import numpy as np
import pandas as pd
from shapely.geometry import box
import pytest
import gisutils
import sfrmaker
import flopy
from sfrmaker.fileio import read_mf6_block


@pytest.fixture
def shellmound_active_area(shellmound_grid, outdir):
    """Make a shapefile of the shellmound bounding box."""
    l, r, b, t, = shellmound_grid.extent
    bbox = box(l, b, r, t)
    df = pd.DataFrame({'geometry': [bbox],
                       'id': [0]})
    out_shapefile = os.path.join(outdir, 'shellmound', 'shellmound_bbox.shp')
    gisutils.df2shp(df, out_shapefile, crs=5070)
    return out_shapefile


@pytest.fixture
def shellmound_config(shellmound_active_area, outdir):
    """Make a version of the example config file with
    shellmound active area."""
    config_file = 'sfrmaker/test/data/shellmound/shellmound_config.yml'
    with open(config_file) as src:
        cfg = yaml.load(src, Loader=yaml.Loader)
    return config_file, cfg


@pytest.fixture(autouse=True)
def keep_cwd():
    """Reset the working directory after a test.
    """
    wd = os.getcwd()
    yield wd  # provide the fixture value
    print("reverting working directory from {} to {}".format(os.getcwd(), wd))
    os.chdir(wd)


def test_from_yaml(shellmound_config):
    """Basic test of setting up an sfr parkage from a yaml file,
    without a model.
    """
    output_config_file, cfg = shellmound_config
    # test the example as-is (but with shellmound extent)
    # from a yaml file
    sfrdata = sfrmaker.SFRData.from_yaml(output_config_file)
    assert isinstance(sfrdata, sfrmaker.SFRData)
    sfr_package_file = 'sfrmaker/test/temp/shellmound/shellmound.sfr'
    assert os.path.exists(sfr_package_file)
    assert os.path.exists('sfrmaker/test/temp/shellmound/shellmound.sfr.obs')
    assert os.path.exists('sfrmaker/test/temp/shellmound/shps/shellmound_sfr_lines.shp')
    assert os.path.exists('sfrmaker/test/temp/shellmound/tables/shellmound_sfr_reach_data.csv')
    options = read_mf6_block(sfr_package_file, 'options')
    assert 'save_flows' in options
    assert options['budget'] == ['fileout', 'shellmound.sfr.cbc']
    assert options['stage'] == ['fileout', 'shellmound.sfr.stage.bin']
    assert options['obs6'] == ['filein', 'shellmound.sfr.obs']

    # test specification of multiple observation types in yaml file
    sfr_obs_file = 'sfrmaker/test/temp/shellmound/shellmound.sfr.obs'
    obsdata = read_mf6_block(sfr_obs_file, 'begin continuous fileout')
    key = list(obsdata.keys())[0]
    obsdata = obsdata[key]
    for line in obsdata:
        if line.startswith('#'):
            continue
        if line.endswith('continuous'):
            break
        obstype = line.split(' ')[0].split('-')[1]
        assert obstype in {'flow', 'stage'}
        

def get_package_version(sfr_package_file):
    with open(sfr_package_file) as src:
        header = next(src)
        if 'MODFLOW-6' in header:
            return 'mf6'
        else:
            return 'mf2005'


@pytest.mark.parametrize('package_version', ('mf6', 'mf2005'))
@pytest.mark.parametrize('model_cfg', ({},
                                       {'simulation': {'sim_name': 'shellmound',
                                                      'sim_ws': 'shellmound'},
                                       'model': {'modelname': 'shellmound'
                                                 }}
                                       )
                         )
def test_from_config(shellmound_config, package_version, model_cfg):
    """Test setting up an sfr parkage from a configuration dictionary,
        - with different package versions specified
          (verify that observations get written regardless)
        - with and without an existing model
    """
    output_config_file, cfg = shellmound_config

    # dict still has same pathing as config file
    os.chdir(os.path.split(output_config_file)[0])
    cfg['package_version'] = package_version
    cfg.update(model_cfg)
    sfrdata = sfrmaker.SFRData.from_yaml(cfg)
    package_file = os.path.join(cfg['output_path'], cfg['package_name'] + '.sfr')

    if len(model_cfg) == 0:
        assert get_package_version(package_file) == package_version
    else:
        assert sfrdata.model.version == 'mf6'
    # verify that observation input was written, even if no model
    if package_version == 'mf6' and sfrdata.model is None:
        assert os.path.exists('../../temp/shellmound/shellmound.sfr.obs')
        os.remove('../../temp/shellmound/shellmound.sfr')
        os.remove('../../temp/shellmound/shellmound.sfr.obs')
    # if there is a model, sfr package and obs input get written to model_ws
    elif sfrdata.model is not None:
        assert os.path.exists('shellmound/shellmound.sfr.obs')
        os.remove('shellmound/shellmound.sfr')
        os.remove('shellmound/shellmound.sfr.obs')
    else:
        assert os.path.exists('../../temp/shellmound/shellmound.gage')
        assert os.path.exists('../../temp/shellmound/shellmound.gage.namefile_entries')
        os.remove('../../temp/shellmound/shellmound.gage')
        os.remove('../../temp/shellmound/shellmound.gage.namefile_entries')

    # check that add_outlets works
    add_outlets = cfg['options']['add_outlets']
    if not isinstance(add_outlets, list):
        add_outlets = [add_outlets]
    for outlet_id in add_outlets:
        assert sfrdata.reach_data.loc[sfrdata.reach_data.line_id == outlet_id,
                                      'outseg'].sum() == 0


@pytest.mark.parametrize('config_file,dem', (
        ('examples/tylerforks/tf_sfrmaker_config.yml', True),
        ('examples/tylerforks/tf_sfrmaker_config2.yml', True),
        ('examples/tylerforks/tf_sfrmaker_config2.yml', False)
)
                         )
def test_tylerforks_from_config(config_file, dem):
    """Test setting up an sfr parkage from a configuration dictionary,
        - with a (georeferenced) MODFLOW-NWT model specified
        - with a model grid defined by a shapefile (tf_sfrmaker_config2.yml)
        - with and without an active area
        - with NHDPlus hydrography by files (tf_sfrmaker_config2.yml)
        - with NHDPlus hydrography with single path (tf_sfrmaker_config.yml)
        - with and without DEM
    """
    with open(config_file) as src:
        cfg = yaml.load(src, Loader=yaml.Loader)

    # dict still has same pathing as config file
    os.chdir(os.path.split(config_file)[0])
    for folder in 'tables/', 'shps/':
        if os.path.exists(folder):
            shutil.rmtree(folder)
    if not dem:
        del cfg['dem']
    sfrdata = sfrmaker.SFRData.from_yaml(cfg)
    if 'model' in cfg:
        assert os.path.exists('tylerforks/tf.sfr')
    else:
        assert os.path.exists('tf.sfr')
    assert os.path.exists('tables/tf_sfr_segment_data.csv')
    assert os.path.exists('tables/tf_sfr_reach_data.csv')
    assert os.path.exists('shps/tf_sfr_cells.shp')
    assert os.path.exists('shps/tf_sfr_lines.shp')
    assert os.path.exists('shps/tf_sfr_outlets.shp')
    assert os.path.exists('shps/tf_sfr_routing.shp')
    assert len(sfrdata.reach_data) >= 155  # number of reaches with active area defined
    rd = sfrdata.reach_data.dropna(subset=['asum'], axis=0)
    assert np.allclose(rd.strtop.values.min(), 1125, rtol=0.01)
    assert sfrdata.grid.crs.srs == 'EPSG:26715'
    
    # check that slopes were updated after sampling streambed tops from the DEM
    reach_data = sfrdata.reach_data.copy()
    reach_data.index = reach_data['rno']
    path = sfrmaker.routing.find_path(sfrdata.rno_routing, 1)[:-1]
    strtop = reach_data.loc[path, 'strtop']
    rchlen = reach_data.loc[path, 'rchlen']
    slopes = -np.diff(strtop)/rchlen[:-1]
    slopes[slopes > sfrdata._maximum_slope] = sfrdata._maximum_slope
    slopes[slopes < sfrdata._minimum_slope] = sfrdata._minimum_slope
    outlets = reach_data.loc[path, 'outreach'] == 0
    slopes[outlets] = sfrdata._default_slope
    assert np.allclose(slopes.values, reach_data.loc[path, 'slope'].values[:-1])
    
    # check that inflows were written correctly
    if 'inflows' in cfg:
        m1 = sfrdata.modflow_sfr2.parent
        m2 = flopy.modflow.Modflow(model_ws='.')
        nper = sfrdata.segment_data['per'].max() + 1
        dis = flopy.modflow.ModflowDis(m2, 
            nrow=sfrdata.grid.nrow, ncol=sfrdata.grid.ncol, 
            nlay=sfrdata.grid.nlay, nper=nper)
        sfr = flopy.modflow.ModflowSfr2.load('tf.sfr', m2)
        for i in range(nper):
            assert np.allclose(sfr.segment_data[i]['flow'], 
                               sfrdata.modflow_sfr2.segment_data[i]['flow'])
        original_inflow_data = pd.read_csv(cfg['inflows']['filename'])
        for i in range(nper):
            per_col = cfg['inflows']['period_column']
            q_col = cfg['inflows']['data_column']
            q_orig = original_inflow_data.loc[original_inflow_data[per_col] == i, q_col].values[0]
            idx = np.where(sfr.segment_data[i]['flow'] > 0)[0][0]
            q_input = sfr.segment_data[i]['flow'][idx]
            assert np.allclose(q_input, q_orig)
