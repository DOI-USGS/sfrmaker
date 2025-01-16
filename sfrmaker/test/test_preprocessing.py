"""
Test the preprocessing module
"""
import os
from pathlib import Path
import sys
import yaml
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
import pytest
from gisutils import df2shp, shp2df
from sfrmaker.checks import check_monotonicity
from sfrmaker.preprocessing import (get_flowline_routing,
                                    cull_flowlines, 
                                    preprocess_nhdplus, 
                                    clip_flowlines_to_polygon, 
                                    edit_flowlines,
                                    swb_runoff_to_csv,
                                    preprocess_nhdplus_hr_flowlines,
                                    preprocess_nhdplus_hr_waterbodies
                                    )


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
def preprocessed_flowlines_module(test_data_path, culled_flowlines, outfolder, project_root_path):

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
    kwargs['dest_crs'] = 5070
    preprocessed_flowlines = preprocess_nhdplus(**kwargs)

    # check that the known_connections were routed correctly
    for comid, tocomid in kwargs['known_connections'].items():
        assert preprocessed_flowlines.loc[comid, 'tocomid'] == tocomid

    out_shapefile = os.path.join(outfolder, 'preprocessed_flowlines.shp')
    df2shp(preprocessed_flowlines, out_shapefile, crs=5070)
    return preprocessed_flowlines


@pytest.fixture(scope='function')
def preprocessed_flowlines(preprocessed_flowlines_module):
    return preprocessed_flowlines_module.copy()
    

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


#@pytest.mark.skipif(sys.version_info[:2] == (3, 11), reason="inexplicable negative asum values")
def test_preprocess_nhdplus_no_zonal_stats(culled_flowlines, preprocessed_flowlines,
                                           test_data_path, outfolder):

    kwargs = culled_flowlines.copy()
    preprocessed_flowlines = preprocessed_flowlines.copy()
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
    kwargs['dest_crs'] = 5070
    preprocessed_flowlines2 = preprocess_nhdplus(**kwargs)

    # verify that the same result is produced
    # when reusing the shapefile output from zonal statistics
    pd.testing.assert_frame_equal(preprocessed_flowlines, preprocessed_flowlines2)
    
    # test manual updating of COMID end elevations
    # (e.g. from measured stage data)
    # assume mean stage at the Money, MS gage
    # (near the upstream end of COMID 17991438)
    # has been measured at 37.1 m (~1.2 meters above the min DEM elevation)
    kwargs['update_up_elevations'] = {17991438: 37.1,
                                }
    # and that mean stage near Greenwood
    # (near the downstream end of COMID 18047242)
    # has been measured at 37.0 m (~1.1 meters above the min DEM elevation)
    kwargs['update_dn_elevations'] = {18047242: 37.0,
                                }
    preprocessed_flowlines3 = preprocess_nhdplus(**kwargs)
    
    assert preprocessed_flowlines3.loc[17991438, 'elevupsmo'] == 37.1
    assert preprocessed_flowlines3.loc[18047242, 'elevdnsmo'] == 37.0


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
    kwargs['dest_crs'] = 5070
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
    if isinstance(flowlines, str) or isinstance(flowlines, Path):
        assert os.path.exists(flowlines[:-4] + '.prj')
        

def test_get_flowline_routing(datapath, project_root_path):
    
    # change the directory to the root level
    # (other tests in this module use the output folder)
    wd = os.getcwd()
    os.chdir(project_root_path)
    NHDPlus_paths = [f'{datapath}/tylerforks/NHDPlus/']
    plusflow_files = [f'{datapath}/tylerforks/NHDPlus/NHDPlusAttributes/PlusFlow.dbf']

    mask = f'{datapath}/tylerforks/grid.shp'
    df = get_flowline_routing(NHDPlus_paths=NHDPlus_paths,
                              mask=mask)
    assert np.array_equal(df.columns, ['FROMCOMID', 'TOCOMID'])
    df2 = get_flowline_routing(PlusFlow=plusflow_files)
    pd.testing.assert_frame_equal(df2.loc[df2['FROMCOMID'].isin(df['FROMCOMID'])].head(),
                                  df.head())
    os.chdir(wd)
    

@pytest.mark.parametrize('limit_runoff_to_area,start_datetime,end_datetime,expected', (
    (None, None, None, (5e5, 9e5)),
    ('runoff_area.shp', '2016-01-01', '2018-01-01', (3000, 1870)),
))    
def test_swb_runoff_to_csv(test_data_path, limit_runoff_to_area, 
                           start_datetime, end_datetime,
                           expected, outdir):
    
    test_data_path = Path(test_data_path)
    if limit_runoff_to_area is not None:
        limit_runoff_to_area = test_data_path / limit_runoff_to_area
    swb_netcdf_output = test_data_path / 'runoff__1999-01-01_to_2018-12-31__989_by_661.nc'
    nhdplus_catchments_file = test_data_path / 'NHDPlus08/NHDPlusCatchment/Catchment.shp'
    outfile = Path(outdir, 'swb_runoff_by_nhdplus_comid.csv')
    aggregated = swb_runoff_to_csv(
        swb_netcdf_output, nhdplus_catchments_file,
        runoff_output_variable='runoff', 
        catchment_id_col='FEATUREID',
        limit_runoff_to_area=limit_runoff_to_area,
        output_length_units='meters',
        include_xy_in_output=True, xy_crs=5070,
        outfile=outfile)
    df = pd.read_csv(outfile)
    df['time'] = pd.to_datetime(df['time'])
    #cat = gpd.read_file(nhdplus_catchments_file)
    # model bounds
    xoffset, yoffset = 500955, 1176285
    nrow, ncol = 30, 35
    dxy = 1000
    x1 = xoffset + ncol * dxy
    y1 = yoffset + nrow * dxy
    within_model = (df.x.values > xoffset) & (df.x.values < x1) & \
                   (df.y.values > yoffset) & (df.y.values < y1)
    df = df.loc[within_model]
    mean_monthly_runoff = df.groupby(df['time'].dt.year)['runoff_m3d'].sum().mean()/12
    # no idea if this is the right number but similar to test results for modflow-setup
    # and serves as a benchmark in case the output changes
    assert np.allclose(mean_monthly_runoff, expected[0], rtol=0.2)
    
    # verify aggregrated runoff for a catchment 
    # that intersects a single SWB cell
    # (this will break if the test data change)
    if limit_runoff_to_area is None:
        loc = aggregated.index.get_level_values(1) == 17955907
        expected = aggregated.loc[loc, 'area_m2']
        # expected amount is the average inches per day
        # converted to meters, then multiplied by
        # SWB cell area of 1 km x 1 km
        assert np.allclose(aggregated.loc[loc, 'runoff_m3d'].mean(), 
                           0.008788843639195 * 1e6 * (.3048/12))
    else:
        # in this case, the expected amount is the same as
        # the aggregated runoff for the single interected SWB cell
        
        assert np.allclose(aggregated['runoff_m3d'].mean(), 
                           0.119300499558449 * 1e6 * (.3048/12))
    
    # test with "rejected net infiltration added"
    swb_rejected_net_inf_output = test_data_path / \
        'irrtest_1000mrejected_net_infiltration__1999-01-01_to_2020-12-31__989_by_661.nc'
    outfile2 = Path(outdir, 'swb_runoff_w_netinf_by_nhdplus_comid.csv')
    swb_runoff_to_csv(swb_netcdf_output, nhdplus_catchments_file,
                      runoff_output_variable='runoff', 
                      swb_rejected_net_inf_output=swb_rejected_net_inf_output,
                      catchment_id_col='FEATUREID',
                      limit_runoff_to_area=limit_runoff_to_area,
                      start_datetime=start_datetime,
                      end_datetime=end_datetime,
                      output_length_units='meters',
                      include_xy_in_output=False,
                      outfile=outfile2)
    df2 = pd.read_csv(outfile2)
    df2['time'] = pd.to_datetime(df2['time'])
    mean_monthly_runoff2 = df2.groupby(df2['time'].dt.year)['runoff_m3d'].sum().mean()/12
    # similar to above; should be a bigger number 
    # with the rejected net inf added
    # not sure if this is right but will fail if the input 
    # or the code changes substantially
    assert np.allclose(mean_monthly_runoff2, expected[1], rtol=0.2)


def test_cull_flowlines2(project_root_path):
    outfolder = Path(project_root_path, 'examples/Notebooks/output')
    outfolder.mkdir(exist_ok=True)

    results = cull_flowlines(NHDPlus_paths=[Path(project_root_path, 'examples/tylerforks/NHDPlus/')],
                            asum_thresh=5, intermittent_streams_asum_thresh=10,
                            cull_invalid=True, cull_isolated=False,
                            active_area=Path(project_root_path, 'examples/tylerforks/active_area.shp'),
                            outfolder=outfolder)
    
    
def test_preprocess_nhdplus_hr_flowlines(project_root_path, outdir):
    data_path = project_root_path / 'examples/neversink_rondout'
    nhdplus_paths = [data_path / 'NHDPlus_HR_1.gdb',  #'../source_data/NHDPLUS_H_19020302_HU8_GDB/NHDPLUS_H_19020302_HU8_GDB.gdb'
                     data_path / 'NHDPlus_HR_2.gdb'
    ]
    model_outline = data_path / 'Model_Extent.shp'    
    outfile = outdir / 'preprocessed_flowlines.shp'
    wb_outfile = outdir / 'preprocessed_waterbodies.shp'

    # drop lines representing storm sewers and aquaducts, etc.
    keep_fcodes = {46003,  # intermittent streams
                   46006,  # perennial streams
                   55800,  # artificial path (thru lakes)
                   33600,  # canal/ditch (includes many wetlands, and the MS River)
                   33400  # connector (some wetlands)
                   }
    
    # drop these IDs and all IDs above them
    drop_ids = {'10000700059483'}
    
    preprocess_nhdplus_hr_flowlines(nhdplus_paths, active_area=model_outline,
                          keep_fcodes=keep_fcodes,
                          drop_ids_upstream=drop_ids,
                          drop_isolated=False,
                          dest_crs=26918, outfile=outfile)
    
    df = gpd.read_file(outfile)
    assert df.crs == 26918
    assert not set(df['FCode']).difference(keep_fcodes)
    assert '10000700059483' not in df['NHDPlusID'].values
    # the next line up shouldn't be there either
    assert '10000700020952' not in df['NHDPlusID'].values


def test_preprocess_nhdplus_hr_waterbodies(project_root_path, outdir):
    nhdplus_path = project_root_path / 'sfrmaker/test/data/nhdplus_hr_waterbodies.shp',  #'../source_data/NHDPLUS_H_19020302_HU8_GDB/NHDPLUS_H_19020302_HU8_GDB.gdb'    
    outfile = outdir / 'preprocessed_waterbodies.shp'
  
    # drop these waterbodies, regardless of size
    drop_waterbodies = {'75004400013339'}
    
    expected_lakes = {#'75004400013854', 
                      '75004400011923', 
                      '75004400012773'}
    
    preprocess_nhdplus_hr_waterbodies(nhdplus_path, 
                                      active_area=(-151.02, 60.64855, -150.96778, 60.67559), 
                                      drop_waterbodies=drop_waterbodies,
                                      min_areasqkm=0.05,
                                      dest_crs=26905, outfile=outfile)
    df = gpd.read_file(outfile)
    df.crs == 26905
    df['NHDPlusID'] = df['NHDPlusID'].astype(int).astype(str)
    assert set(df['NHDPlusID']) == expected_lakes
    # this lake is < 0.05 km2; should have been culled
    assert '75004400012864' not in df['NHDPlusID'].values
    for nhdplusid in drop_waterbodies:
        assert int(nhdplusid) not in df['NHDPlusID'].values
