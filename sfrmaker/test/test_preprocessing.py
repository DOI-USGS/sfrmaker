"""
Test the preprocessing module
"""
import os
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
from shapely.geometry import box
import pytest
from gisutils import df2shp, shp2df
from sfrmaker.checks import check_monotonicity
from sfrmaker.preprocessing import cull_flowlines, preprocess_nhdplus, clip_flowlines_to_polygon, edit_flowlines


@pytest.fixture(scope='module')
def test_data_path(project_root_path):
    return os.path.join(project_root_path, 'sfrmaker/test/data/shellmound/')


@pytest.fixture(scope='module')
def active_area(outfolder):
    active_area_tuple = -90.55, 33.5, -90.16, 33.86
    active_area_poly = box(*active_area_tuple)
    df = pd.DataFrame({'geometry': [active_area_poly], 'id': [0]})
    active_area = os.path.join(outfolder, 'active_area.shp')
    df2shp(df, active_area, crs=4269)
    return active_area_tuple


@pytest.fixture(scope='module')
def outfolder(outdir):
    outfolder = os.path.join(outdir, 'preprocessing')
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    return outfolder


@pytest.fixture(autouse=True)
def use_outfolder(outdir, outfolder):
    wd = os.getcwd()
    os.chdir(outfolder)
    yield outfolder
    os.chdir(wd)


# test with tuple bounding box and with shapefile
@pytest.fixture(scope='module')
def clipped_flowlines(test_data_path, outfolder, active_area):
    #nhdpaths = ['/Users/aleaf/Documents/NHDPlus/NHDPlusMS/NHDPlus08']
    nhdpaths = [os.path.join(test_data_path, 'NHDPlus08')]
    active_area = os.path.join(outfolder, 'active_area.shp')
    results = cull_flowlines(nhdpaths,
                             asum_thresh=None, intermittent_streams_asum_thresh=None,
                             cull_invalid=False, cull_isolated=False,
                             active_area=active_area,
                             outfolder=outfolder)
    return results


@pytest.fixture(scope='module')
def culled_flowlines(test_data_path, outfolder, active_area):
    nhdpaths = [os.path.join(test_data_path, 'NHDPlus08')]
    results = cull_flowlines(nhdpaths,
                             asum_thresh=20, intermittent_streams_asum_thresh=50,
                             cull_invalid=True, cull_isolated=True,
                             active_area=active_area,
                             outfolder=outfolder)
    return results


@pytest.fixture(scope='module')
def preprocessed_flowlines(test_data_path, culled_flowlines, outfolder, project_root_path):

    kwargs = culled_flowlines.copy()
    #kwargs['demfile'] = os.path.join(test_data_path, 'meras_100m_dem.tif')
    kwargs['demfile'] = os.path.join(test_data_path, 'meras_30m_dem.tif')
    #kwargs['demfile'] = os.path.join(project_root_path, 'examples/meras/dem_min_elevs_1000.tif')
    kwargs['dem_length_units'] = 'feet'
    kwargs['narwidth_shapefile'] = os.path.join(test_data_path, 'NARwidth.shp')
    kwargs['waterbody_shapefiles'] = os.path.join(test_data_path,
                                                  'NHDPlus08/NHDSnapshot/Hydrography/NHDWaterbody.shp')
    kwargs['asum_thresh'] = 20.
    kwargs['width_from_asum_a_param'] = 0.0592
    kwargs['width_from_asum_b_param'] = 0.5127
    kwargs['known_connections'] = {17955195: 17955197,
                                   17955197: 17955185,
                                   17954979: 17954993,
                                   17954993: 17955075
                                   }
    kwargs['logger'] = None
    kwargs['output_length_units'] = 'meters'
    kwargs['outfolder'] = outfolder
    kwargs['project_epsg'] = 5070
    preprocessed_flowlines = preprocess_nhdplus(**kwargs)

    # check that the known_connections were routed correctly
    for comid, tocomid in kwargs['known_connections'].items():
        assert preprocessed_flowlines.loc[comid, 'tocomid'] == tocomid

    out_shapefile = os.path.join(outfolder, 'preprocessed_flowlines.shp')
    df2shp(preprocessed_flowlines, out_shapefile, crs=5070)
    return preprocessed_flowlines


def test_cull_flowlines(clipped_flowlines, culled_flowlines, test_data_path,
                                    outfolder, active_area):
    nhdpaths = [os.path.join(test_data_path, 'NHDPlus08')]
    source_nhdfiles = [os.path.join(nhdpaths[0], 'NHDSnapshot/Hydrography/NHDFlowline.shp'),
                       os.path.join(nhdpaths[0], 'NHDPlusAttributes/PlusFlowlineVAA.dbf'),
                       os.path.join(nhdpaths[0], 'NHDPlusAttributes/PlusFlow.dbf'),
                       os.path.join(nhdpaths[0], 'NHDPlusAttributes/elevslope.dbf')
                       ]
    original_sizes = np.array([os.path.getsize(f) for f in source_nhdfiles])
    results = clipped_flowlines
    clipped_sizes = np.array([os.path.getsize(f) for f in clipped_flowlines.values()])
    culled_sizes = np.array([os.path.getsize(f) for f in culled_flowlines.values()])
    assert np.all(culled_sizes > 0)
    assert np.all(original_sizes > 0)
    assert np.all(culled_sizes <= clipped_sizes)
    assert np.all(clipped_sizes <= original_sizes)

    results2 = cull_flowlines(nhdpaths,
                              asum_thresh=None, intermittent_streams_asum_thresh=None,
                              cull_invalid=False, cull_isolated=False,
                              active_area=active_area,
                              outfolder=outfolder)
    assert results == results2
    sizes2 = np.array([os.path.getsize(f) for f in results.values()])
    assert np.all(sizes2 <= original_sizes)
    assert results != culled_flowlines


def test_preprocess_nhdplus(preprocessed_flowlines):
    fl = preprocessed_flowlines

    # check some more connections
    connections = {17955689: 17955711,
                   17956683: 17956745
                   }
    for comid, tocomid in connections.items():
        assert fl.loc[comid, 'tocomid'] == tocomid

    # these lines should have been dropped
    should_have_been_dropped = {17956691, 17955223}
    for comid in should_have_been_dropped:
        assert comid not in fl.index

    # these lines should be included
    should_be_included = {17955197, 17956745}
    for comid in should_be_included:
        assert comid in fl.index

    # check that arbolate sums at known connections are correct
    assert np.allclose(fl.loc[17955197, 'asum_calc'], 650., rtol=0.1)

    # check that arbolate sums increase monotonically downstream
    assert check_monotonicity(fl.index, fl.tocomid, fl.asum_calc, decrease=False)

    # verify that for the lines that have narwidth estimates,
    # the mean narwidth width and asum widths are within 20%
    # (should hopefully catch unit conversion mistakes)
    has_nw = ~fl.narwd_n.isna()
    np.allclose(fl.loc[has_nw, 'width1asum'].mean(),
                fl.loc[has_nw, 'narwd_mean'].mean(), rtol=0.2)


def test_preprocess_nhdplus_no_zonal_stats(culled_flowlines, preprocessed_flowlines,
                                           test_data_path, outfolder):

    kwargs = culled_flowlines.copy()
    # kwargs['demfile'] = os.path.join(test_data_path, 'meras_100m_dem.tif')
    kwargs['run_zonal_statistics'] = False
    kwargs['flowline_elevations_file'] = Path(outfolder, 'flowlines_gt20km_buffers.shp')
    kwargs['narwidth_shapefile'] = os.path.join(test_data_path, 'NARwidth.shp')
    kwargs['waterbody_shapefiles'] = os.path.join(test_data_path,
                                                  'NHDPlus08/NHDSnapshot/Hydrography/NHDWaterbody.shp')
    kwargs['asum_thresh'] = 20.
    kwargs['width_from_asum_a_param'] = 0.0592
    kwargs['width_from_asum_b_param'] = 0.5127
    kwargs['known_connections'] = {17955195: 17955197,
                                   17955197: 17955185,
                                   17954979: 17954993,
                                   17954993: 17955075
                                   }
    kwargs['logger'] = None
    kwargs['output_length_units'] = 'meters'
    kwargs['outfolder'] = outfolder
    kwargs['project_epsg'] = 5070
    preprocessed_flowlines2 = preprocess_nhdplus(**kwargs)

    # verify that the same result is produced
    # when reusing the shapefile output from zonal statistics
    pd.testing.assert_frame_equal(preprocessed_flowlines, preprocessed_flowlines2)


@pytest.mark.timeout(30)  # projection issues will cause zonal stats to hang
def test_preprocess_nhdplus_no_narwidth(test_data_path, culled_flowlines, outfolder):
    kwargs = culled_flowlines.copy()
    kwargs['demfile'] = os.path.join(test_data_path, 'meras_100m_dem.tif')
    kwargs['narwidth_shapefile'] = None
    kwargs['waterbody_shapefiles'] = None
    kwargs['asum_thresh'] = 20.
    kwargs['known_connections'] = None
    kwargs['logger'] = None
    kwargs['outfolder'] = outfolder
    kwargs['project_epsg'] = 5070
    preprocess_nhdplus(**kwargs)


def test_clip_flowlines(preprocessed_flowlines, test_data_path):
    clipped = clip_flowlines_to_polygon(preprocessed_flowlines,
                                        os.path.join(test_data_path, 'active_area.shp'),
                                        simplify_tol=100, logger=None)
    assert len(clipped) < len(preprocessed_flowlines)


@pytest.mark.parametrize('flowlines', ('preprocessed_flowlines.shp',
                                       None))
def test_edit_flowlines(flowlines, preprocessed_flowlines, test_data_path):
    if flowlines is None:
        flowlines = preprocessed_flowlines
    flowline_edits_file = os.path.join(test_data_path, 'flowline_edits.yml')
    edited_flowlines = edit_flowlines(flowlines,
                                      flowline_edits_file, logger=None)
    with open(flowline_edits_file) as src:
        cfg = yaml.load(src, Loader=yaml.Loader)
    # verify that flowlines were dropped
    assert not any(set(cfg['drop_flowlines']).intersection(edited_flowlines.COMID))
    # verify routing changes
    for comid, tocomid in cfg['reroute_flowlines'].items():
        assert edited_flowlines.loc[comid, 'tocomid'] == tocomid
    add_flowlines = shp2df(os.path.join(test_data_path, 'yazoo.shp'))
    assert not any(set(add_flowlines.comid).difference(edited_flowlines.index))
    if isinstance(flowlines, str):
        assert os.path.exists(flowlines[:-4] + '.prj')