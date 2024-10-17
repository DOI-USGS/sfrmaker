"""Unit tests for test_nhdplus_utils.py"""
import pytest
from sfrmaker.nhdplus_utils import (
    load_nhdplus_v2,
    read_nhdplus,
    read_nhdplus_hr
)


@pytest.mark.parametrize('nhdplus_files,bbox_filter,expected_nrows,drop_fcodes', (
    ('neversink_rondout/NHDPlus_HR_1.gdb', 
     (-74.68061523310072, 41.74556802490179, -74.44311875132166, 41.848862965337595),
     132, None),
    ('neversink_rondout/NHDPlus_HR_1.gdb', 
     'examples/neversink_rondout/Model_Extent.shp',
     75, 55800),
    (['neversink_rondout/NHDPlus_HR_1.gdb', 'neversink_rondout/NHDPlus_HR_2.gdb'],
     'examples/neversink_rondout/Model_Extent.shp', 386, None
     )
)
)
def test_read_nhdplus_hr(nhdplus_files, bbox_filter, expected_nrows, drop_fcodes,
                         datapath):

    if not isinstance(nhdplus_files, str):
        NHDPlusHR_paths = [datapath / f for f in nhdplus_files]
    else:
        NHDPlusHR_paths = datapath / nhdplus_files
    results = read_nhdplus_hr(NHDPlusHR_paths, bbox_filter=bbox_filter, 
                              drop_fcodes=drop_fcodes)
    expected_cols = ['NHDPlusID', 'ArbolateSu','StreamOrde', 'MaxElevSmo', 
                     'MinElevSmo', 'Divergence', 'ToNHDPID']
    assert not set(expected_cols).difference(results.columns)
    assert len(results) == expected_nrows
    # NHDPlusIDs should be treated as strings
    assert all([type(s) == str for s in results['NHDPlusID'].astype(int).astype(str)])
    assert all([type(s) == str for s in results['ToNHDPID'].astype(int).astype(str)])


@pytest.mark.parametrize('nhdplus_files,bbox_filter,expected_nrows', (
    ('tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp', 
     (-90.59552642527598, 46.37324928457199, -90.45520721646749, 46.43603727904705),
     27),
    ('tylerforks/NHDPlus/NHDPlusAttributes/elevslope.dbf', None, 101),
    ('tylerforks/NHDPlus/NHDPlusAttributes/PlusFlow.dbf', None, 101),
    ('tylerforks/NHDPlus/NHDPlusAttributes/PlusFlowlineVAA.dbf', None, 101),
    (['tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp', 
      'tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp'],
     'tylerforks/active_area.shp', 27*2
     )
)
)
def test_read_nhdplus(nhdplus_files, bbox_filter, expected_nrows,
                      datapath):
    if not isinstance(nhdplus_files, str):
        nhdplus_paths = [datapath / f for f in nhdplus_files]
    else:
        nhdplus_paths = datapath / nhdplus_files
    if isinstance(bbox_filter, str):
        bbox_filter = datapath / bbox_filter
    results = read_nhdplus(nhdplus_paths, bbox_filter=bbox_filter)
    assert len(results) == expected_nrows
    # NHDPlusIDs should be treated as strings
    comid_col = [c for c in results.columns if c.lower() in {'comid', 'froncomid', 'tocomid'}]
    assert all([type(s) == str for s in results[comid_col].astype(int).astype(str)])


@pytest.mark.parametrize('kwargs,bbox_filter', (
    ({'NHDPlus_paths': 'tylerforks/NHDPlus'}, 'tylerforks/active_area.shp'),
    ({'NHDFlowlines': 'tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp',
      'elevslope': 'tylerforks/NHDPlus/NHDPlusAttributes/elevslope.dbf',
      'PlusFlowlineVAA': 'tylerforks/NHDPlus/NHDPlusAttributes/PlusFlowlineVAA.dbf',
      'PlusFlow': 'tylerforks/NHDPlus/NHDPlusAttributes/PlusFlow.dbf'
      }, 'tylerforks/active_area.shp')
))
def test_load_nhdplus_v2(kwargs, bbox_filter, datapath):
    kwargs = {k: datapath / v for k, v in kwargs.items()}
    if isinstance(bbox_filter, str):
        bbox_filter = datapath / bbox_filter
    results = load_nhdplus_v2(**kwargs, bbox_filter=bbox_filter)
    # all comids and tocomids should be strings
    assert all([type(s) == str for s in results['COMID'].astype(int).astype(str)])
    assert all([type(s) == str for comids in results['COMID'].astype(int).astype(str) 
         for s in comids ])
