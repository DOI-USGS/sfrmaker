"""
Test the preprocessing module
"""
import os
import yaml
import numpy as np
import pandas as pd
from shapely.geometry import box
import pytest
from gisutils import df2shp, shp2df
from sfrmaker.preprocessing import cull_flowlines, preprocess_nhdplus, clip_flowlines_to_polygon, edit_flowlines


@pytest.fixture(scope='module')
def test_data_path(project_root_path):
    return os.path.join(project_root_path, 'sfrmaker/test/data/shellmound/')


@pytest.fixture(scope='module')
def active_area():
    return -90.55, 33.5, -90.16, 33.78


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
    active_area_poly = box(*active_area)
    df = pd.DataFrame({'geometry': [active_area_poly], 'id': [0]})
    active_area = os.path.join(outfolder, 'active_area.shp')
    df2shp(df, active_area, epsg=4269)
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
def preprocessed_flowlines(test_data_path, culled_flowlines, outfolder):

    kwargs = culled_flowlines.copy()
    kwargs['demfile'] = os.path.join(test_data_path, 'meras_100m_dem.tif')
    kwargs['dem_length_units'] = 'feet'
    kwargs['narwidth_shapefile'] = os.path.join(test_data_path, 'NARwidth.shp')
    kwargs['waterbody_shapefiles'] = os.path.join(test_data_path,
                                                  'NHDPlus08/NHDSnapshot/Hydrography/NHDWaterbody.shp')
    kwargs['asum_thresh'] = 20.
    kwargs['width_from_asum_a_param'] = 0.0592
    kwargs['width_from_asum_b_param'] = 0.5127
    kwargs['known_connections'] = None
    kwargs['logger'] = None
    kwargs['output_length_units'] = 'meters'
    kwargs['outfolder'] = outfolder
    kwargs['project_epsg'] = 5070
    preprocessed_flowlines = preprocess_nhdplus(**kwargs)
    out_shapefile = os.path.join(outfolder, 'preprocessed_flowlines.shp')
    df2shp(preprocessed_flowlines, out_shapefile, epsg=5070)
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
    assert np.allclose(original_sizes, clipped_sizes)

    results2 = cull_flowlines(nhdpaths,
                              asum_thresh=None, intermittent_streams_asum_thresh=None,
                              cull_invalid=False, cull_isolated=False,
                              active_area=active_area,
                              outfolder=outfolder)
    assert results == results2
    sizes2 = np.array([os.path.getsize(f) for f in results.values()])
    assert np.allclose(original_sizes, sizes2)
    assert results != culled_flowlines


def test_preprocess_nhdplus(preprocessed_flowlines):
    fl = preprocessed_flowlines

    # verify that for the lines that have narwidth estimates,
    # the mean narwidth width and asum widths are within 20%
    # (should hopefully catch unit conversion mistakes)
    has_nw = ~fl.narwd_n.isna()
    np.allclose(fl.loc[has_nw, 'width1asum'].mean(),
                fl.loc[has_nw, 'narwd_mean'].mean(), rtol=0.2)


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