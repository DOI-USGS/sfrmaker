import copy
import os

import numpy as np
import pandas as pd
import geopandas as gpd
import pytest

import sfrmaker
from sfrmaker.gis import get_authority_crs
from gisutils import project
from sfrmaker import Lines
from .test_grid import tyler_forks_grid_from_shapefile


@pytest.fixture(scope='module')
def dem(datapath):
    return os.path.join(datapath, 'tylerforks/dem_26715.tif')


def test_attribute_units(get_tylerforks_lines_from_NHDPlus, tylerforks_sfrdata):
    tylerforks_lines_from_NHDPlus = copy.deepcopy(get_tylerforks_lines_from_NHDPlus)
    assert tylerforks_lines_from_NHDPlus.elevation_units == 'meters'
    assert tylerforks_lines_from_NHDPlus.width_units == 'meters'
    assert tylerforks_sfrdata.model_length_units == 'feet'
    lines_mean_strtop = tylerforks_lines_from_NHDPlus.df[['elevup', 'elevdn']].mean().mean()
    rd = tylerforks_sfrdata.reach_data
    # some lines along the model boundary aren't included 
    # in the test subset of the elevslope database
    # since this case isn't using the DEM, they end up with elevations of 0
    sfrdata_mean_strtop = rd.loc[rd['strtop'] != 0, 'strtop'].mean()
    assert np.allclose(sfrdata_mean_strtop / lines_mean_strtop, 3.28, rtol=0.2)


def test_output_dir(tylerforks_model, outdir):
    assert tylerforks_model.model_ws == os.path.join(outdir, 'tylerforks')


def test_sample_elevations_different_proj(dem, tylerforks_sfrdata, datapath):
    sfr = tylerforks_sfrdata
    sampled_elevs1 = sfr.sample_reach_elevations(dem, method='buffers', smooth=True)
    sampled_elevs1 = np.array(list(sampled_elevs1.values()))

    reach1_geom = sfr.reach_data.geometry[0]
    crs1 = sfr.crs
    crs2 = get_authority_crs(3070)
    sfr._crs = crs2
    sfr.reach_data['geometry'] = project(sfr.reach_data['geometry'].values, crs1, crs2)
    reach1_geom_5070 = sfr.reach_data.geometry[0]

    # verify that the reaches were reprojected
    assert reach1_geom.intersection(reach1_geom_5070).area == 0
    sampled_elevs2 = sfr.sample_reach_elevations(dem, method='buffers', smooth=True)
    sampled_elevs2 = np.array(list(sampled_elevs2.values()))
    rms_error = np.sqrt(np.mean((sampled_elevs2 - sampled_elevs1) ** 2))
    assert rms_error < 0.5  # not sure why the elevations don't match better

    # verify that at least the first reach is the same
    reach1_geom_projected_back_100buffer = project(reach1_geom_5070, crs2, crs1).buffer(100)
    assert np.allclose(reach1_geom_projected_back_100buffer.area, reach1_geom.buffer(100).area)


# ugly work-around for fixtures not being supported as test parameters yet
# https://github.com/pytest-dev/pytest/issues/349
@pytest.fixture(params=['tylerforks_sfrmaker_grid_from_flopy',
                        'tyler_forks_grid_from_shapefile'])
def grid(request,
         tylerforks_sfrmaker_grid_from_flopy,
         tyler_forks_grid_from_shapefile):
    return {'tylerforks_sfrmaker_grid_from_flopy': tylerforks_sfrmaker_grid_from_flopy,
            'tyler_forks_grid_from_shapefile': tyler_forks_grid_from_shapefile}[request.param]


def test_lines_from_NHDPlus(tylerforks_lines_from_NHDPlus):
    lines = tylerforks_lines_from_NHDPlus
    
    tf = lines.df.name == 'Tyler Forks'
    lines_pr = project(lines.df.geometry, 4269, 26915)
    line_lengths = np.array([g.length for g in lines_pr])
    expected_asum1s = lines.df['asum2'] - line_lengths
    # add dropna due to some lines along boundary 
    # not being in the PFVAA subset used for the test
    assert np.allclose(lines.df['asum1'].dropna(), expected_asum1s.dropna(), atol=10)
    assert np.all(lines.df.loc[tf, 'asum1'].dropna() > 90000)
    assert isinstance(lines, Lines)


def test_make_sfr(outdir,
                  grid,
                  tylerforks_model,
                  tylerforks_lines_from_NHDPlus,
                  tylerforks_active_area_shapefile,
                  dem, mfnwt_exe):
    m = tylerforks_model
    sfr = tylerforks_lines_from_NHDPlus.to_sfr(grid=grid,
                                               model=m)
    assert sfr.crs == grid.crs  
    sfr.set_streambed_top_elevations_from_dem(dem, elevation_units='meters')

    # verify that the minimum width on the Tyler Forks
    # which has substantial drainage upstream of the test case domain
    # is at least 30 feet (that the first reach in the model was not 
    # assigned the default width of 1 but was instead computed from a non-zero asum value)
    rd = sfr.reach_data.dropna(subset=['asum'], axis=0)
    min_tf_width = rd.loc[rd.name == 'Tyler Forks', 'width'].min()
    assert min_tf_width > 30.
    
    botm = m.dis.botm.array.copy()
    layers, new_botm = sfrmaker.utils.assign_layers(sfr.reach_data, botm_array=botm)
    sfr.reach_data['k'] = layers
    if new_botm is not None:
        botm[-1] = new_botm
        np.savetxt('{}/external/botm{}.dat'.format(m.model_ws,
                                                   m.nlay - 1),
                   new_botm, fmt='%.2f')
        sfr.modflow_sfr2.parent.dis.botm = botm

    sfr.create_modflow_sfr2(model=m)
    sfr.modflow_sfr2.check()
    sfr.write_package(istcb2=223)  # writes a sfr file to the model workspace
    m.write_input()
    #m.write_name_file()  # write new version of name file with sfr package

    # wite shapefiles for visualization
    sfr.export_cells(outdir / 'example_cells.shp')
    sfr.export_outlets(outdir / 'example_outlets.shp')
    sfr.export_lines(outdir / 'example_lines.shp')
    sfr.export_routing(outdir / 'example_routing.shp')

    # run modflow
    if sfrmaker.utils.exe_exists(mfnwt_exe):
        m.exe_name = mfnwt_exe
        try:
            success, buff = m.run_model(silent=False)
        except:
            pass
        if not success:
            list_file = m.lst.fn_path
            with open(list_file) as src:
                list_output = src.read()
        assert success, 'model run did not terminate successfully:\n{}'.format(list_output)
        return m
