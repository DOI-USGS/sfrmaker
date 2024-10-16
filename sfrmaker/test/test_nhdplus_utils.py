"""Unit tests for test_nhdplus_utils.py"""
import pytest
from sfrmaker.nhdplus_utils import (
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