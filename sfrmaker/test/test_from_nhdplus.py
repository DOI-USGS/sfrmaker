import os
import shutil

import numpy as np
import pytest

import sfrmaker
from sfrmaker.checks import reach_elevations_decrease_downstream
from sfrmaker.gis import crs
from gisutils import project
from sfrmaker import Lines


@pytest.fixture(scope='module')
def dem(datapath):
    return os.path.join(datapath, 'badriver/dem_26715.tif')


@pytest.fixture(scope='function')
def sfrmaker_grid_from_shapefile(tyler_forks_grid_shapefile, tylerforks_active_area_shapefile):
    grid = sfrmaker.StructuredGrid.from_shapefile(tyler_forks_grid_shapefile,
                                                  node_col='node',
                                                  icol='i',
                                                  jcol='j',
                                                  active_area=tylerforks_active_area_shapefile)
    return grid


def test_attribute_units(tylerforks_lines_from_NHDPlus, tylerforks_sfrdata):
    assert tylerforks_lines_from_NHDPlus.attr_height_units == 'meters'
    assert tylerforks_lines_from_NHDPlus.attr_length_units == 'meters'
    assert tylerforks_sfrdata.model_length_units == 'feet'
    lines_mean_strtop = tylerforks_lines_from_NHDPlus.df[['elevup', 'elevdn']].mean().mean()
    sfrdata_mean_strtop = tylerforks_sfrdata.reach_data.strtop.mean()
    assert np.allclose(sfrdata_mean_strtop / lines_mean_strtop, 3.28, atol=0.01)


@pytest.mark.parametrize('method', ['cell polygons', 'buffers'])
def test_sample_elevations(dem, tylerforks_sfrdata, datapath, method):
    sfr = tylerforks_sfrdata
    sampled_elevs = sfr.sample_reach_elevations(dem, method=method, smooth=True)
    sfr.reach_data['strtop'] = [sampled_elevs[rno] for rno in sfr.reach_data['rno']]
    assert reach_elevations_decrease_downstream(sfr.reach_data)


def test_output_dir(tylerforks_model, outdir):
    assert tylerforks_model.model_ws == os.path.join(outdir, 'tylerforks')


def test_sample_elevations_different_proj(dem, tylerforks_sfrdata, datapath):
    sfr = tylerforks_sfrdata
    sampled_elevs1 = sfr.sample_reach_elevations(dem, method='buffers', smooth=True)
    sampled_elevs1 = np.array(list(sampled_elevs1.values()))

    reach1_geom = sfr.reach_data.geometry[0]
    crs1 = sfr.crs
    crs2 = crs(epsg=3070)
    sfr._crs = crs2
    sfr.reach_data['geometry'] = project(sfr.reach_data['geometry'].values, crs1.proj_str, crs2.proj_str)
    reach1_geom_5070 = sfr.reach_data.geometry[0]

    # verify that the reaches were reprojected
    assert reach1_geom.intersection(reach1_geom_5070).area == 0
    sampled_elevs2 = sfr.sample_reach_elevations(dem, method='buffers', smooth=True)
    sampled_elevs2 = np.array(list(sampled_elevs2.values()))
    rms_error = np.sqrt(np.mean((sampled_elevs2 - sampled_elevs1) ** 2))
    assert rms_error < 0.5  # not sure why the elevations don't match better

    # verify that at least the first reach is the same
    reach1_geom_projected_back_100buffer = project(reach1_geom_5070, crs2.proj_str, crs1.proj_str).buffer(100)
    assert np.allclose(reach1_geom_projected_back_100buffer.area, reach1_geom.buffer(100).area)


def test_structuredgrid_from_shapefile(sfrmaker_grid_from_shapefile, tylerforks_sfrmaker_grid_from_flopy):
    grid = sfrmaker_grid_from_shapefile
    grid_flopy = tylerforks_sfrmaker_grid_from_flopy
    assert grid == grid_flopy


def test_unstructuredgrid_from_shapfile(tyler_forks_grid_shapefile,
                                        tylerforks_sfrmaker_grid_from_flopy):
    # TODO: test creating unstructured grid from same shapefile
    # with no row or column information passed
    pass


# ugly work-around for fixtures not being supported as test parameters yet
# https://github.com/pytest-dev/pytest/issues/349
@pytest.fixture(params=['tylerforks_sfrmaker_grid_from_flopy',
                        'sfrmaker_grid_from_shapefile'])
def grid(request,
         tylerforks_sfrmaker_grid_from_flopy,
         sfrmaker_grid_from_shapefile):
    return {'tylerforks_sfrmaker_grid_from_flopy': tylerforks_sfrmaker_grid_from_flopy,
            'sfrmaker_grid_from_shapefile': sfrmaker_grid_from_shapefile}[request.param]


def test_lines_from_NHDPlus(tylerforks_lines_from_NHDPlus):
    lines = tylerforks_lines_from_NHDPlus
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
    sfr.set_streambed_top_elevations_from_dem(dem, dem_z_units='meters')

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
    sfr.export_cells(outdir + 'example_cells.shp')
    sfr.export_outlets(outdir + 'example_outlets.shp')
    sfr.export_lines(outdir + 'example_lines.shp')
    sfr.export_routing(outdir + 'example_routing.shp')

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
