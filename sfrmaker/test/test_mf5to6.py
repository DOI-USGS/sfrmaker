# TODO: add unit tests for mf5to6.py
import copy
import filecmp
import os
import platform
import numpy as np
import pandas as pd
import pytest
import flopy
import sfrmaker


@pytest.fixture(scope='function')
def shellmound_ModflowSfr2(shellmound_sfrdata):
    return shellmound_sfrdata.modflow_sfr2


@pytest.fixture(scope='function')
def mf6sfr_instance_ModflowSfr2(shellmound_ModflowSfr2):
    return sfrmaker.Mf6SFR(shellmound_ModflowSfr2)


@pytest.fixture(scope='function')
def mf6sfr_instance_SFRdata(shellmound_sfrdata):
    return sfrmaker.Mf6SFR(SFRData=shellmound_sfrdata)


@pytest.mark.parametrize('options', (None, 
                                     ['save_flows',
                                      'BUDGET FILEOUT sfr.cbc',
                                      'STAGE FILEOUT sr.stage.bin',
                                     ],
                                     ['auxiliary line_id',
                                      'unit_conversion 1.0']
                                     ),
                         )
def test_init(options, shellmound_ModflowSfr2):
    mf6sfr = sfrmaker.Mf6SFR(shellmound_ModflowSfr2, options=options)
    options_block = mf6sfr.options_block.strip('\n').split('\n')
    expected_keys = ['BEGIN']
    if options is not None:
        expected_keys = expected_keys[:1] + [' '.join(item.split()[:-1]) if len(item.split()) > 1 else item 
                        for item in options]
    if 'unit_conversion' not in expected_keys:
        expected_keys.append('unit_conversion')
    if 'auxiliary' not in expected_keys:
        expected_keys.append('auxiliary')
    expected_keys.append('END')
    keys = [' '.join(item.split()[:-1]) if len(item.split()) > 1 else item.strip() 
                    for item in options_block]
    # verify that there aren't any duplicate entries
    assert len(set(keys)) == len(keys)
    assert keys == expected_keys


def test_connectiondata(mf6sfr_instance_SFRdata, outdir):
    """Test handling of unconnected reaches"""
    mf6sfr = copy.copy(mf6sfr_instance_SFRdata)
    # remove a routing connection and then rebuild connectiondata
    remove_reach = 1
    del mf6sfr.graph[remove_reach]
    mf6sfr._connections = None
    assert remove_reach in mf6sfr.connections
    # each reach should have an entry in connectiondata
    # including unconnected reaches
    expected_connectiondata_entries = np.arange(len(mf6sfr.connections.keys())) + 1
    assert np.allclose(expected_connectiondata_entries, sorted(mf6sfr.connections.keys()))
    outfile = os.path.join(outdir, 'junk.sfr')
    mf6sfr.write_file(filename=outfile)


def test_packagedata_aux(mf6sfr_instance_SFRdata):
    mf6sfr = mf6sfr_instance_SFRdata
    packagedata = mf6sfr._get_packagedata()
    assert 'line_id' in packagedata.columns
    

def test_segment_data_to_perioddata(shellmound_sfrdata):
    # shellmound_sfrdata = copy.deepcopy(shellmound_sfrdata)
    shellmound_sfrdata.segment_data.loc[2, 'flow'] = 100
    perdata = shellmound_sfrdata.period_data
    assert isinstance(perdata, pd.DataFrame)
    assert len(perdata) == 1
    assert perdata.index.names == ['per', 'rno']
    assert set(perdata.columns) == \
           {'iseg', 'inflow', 'status',
            'ireach', 'icalc'}


def test_idomain(mf6sfr_instance_ModflowSfr2, shellmound_model):
    ibound = mf6sfr_instance_ModflowSfr2.idomain
    idomain = shellmound_model.dis.idomain.array
    assert np.array_equal(ibound, idomain)


@pytest.mark.xfail(platform.system() == 'Windows',
                   reason="for some reason filecmp.cmp fails on Windows but not OSX or Linux")
@pytest.mark.parametrize('external_files_path', ('external', None))
def test_write(shellmound_sfrdata, mf6sfr_instance_SFRdata, outdir, external_files_path):
    mf6sfr = mf6sfr_instance_SFRdata
    outfile = os.path.join(outdir, 'junk.sfr')
    outfile2 = os.path.join(outdir, 'junk2.sfr')
    outfile3 = os.path.join(outdir, 'junk3.sfr')
    
    if external_files_path is not None:
        full_external_files_path = os.path.join(outdir, external_files_path)
        if not os.path.isdir(full_external_files_path):
            os.makedirs(full_external_files_path)

    _, outfile2_name = os.path.split(outfile2)
    mf6sfr.write_file(filename=outfile,
                      options=['save_flows',
                               'BUDGET FILEOUT {}.cbc'.format(outfile2_name),
                               'STAGE FILEOUT {}.stage.bin'.format(outfile2_name),
                               ],
                      external_files_path=external_files_path
                      )

    shellmound_sfrdata.write_package(outfile2, version='mf6', 
                                     external_files_path=external_files_path)
    with open(outfile2) as src:
        text = src.read().replace('junk2_packagedata.dat', 'junk_packagedata.dat')
        with open(outfile3, 'w') as dest:
            dest.write(text)
    assert filecmp.cmp(outfile, outfile3)
    
    if external_files_path:
        for version in 'junk', 'junk2':
            assert os.path.exists(os.path.join(os.path.split(outfile2)[0],
                                            external_files_path,
                                            '{}_packagedata.dat'.format(version)))


def test_convert_mf5_to_6(datapath):
    model = flopy.modflow.Modflow.load('tf_with_sfr.nam', model_ws=datapath / 'tylerforks/tylerforks')
    mf5to6 = sfrmaker.Mf6SFR(model.sfr)
    
    # check that aux columns, including default line_id, are included in packagedata
    # otherwise MODFLOW 6 with throw an error
    if 'auxiliary' in mf5to6.options_block:
        cols = mf5to6.options_block.split('auxiliary')[1].strip().split('\n')[0].split()
        for col in cols:
            assert col in mf5to6.packagedata.columns
