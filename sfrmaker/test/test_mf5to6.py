# TODO: add unit tests for mf5to6.py
import filecmp
import os

import numpy as np
import pandas as pd
import pytest

import sfrmaker


@pytest.fixture(scope='function')
def shellmound_ModflowSfr2(shellmound_sfrdata):
    return shellmound_sfrdata.modflow_sfr2


@pytest.fixture(scope='function')
def mf6sfr_instance(shellmound_ModflowSfr2):
    return sfrmaker.mf6sfr(shellmound_ModflowSfr2)


def test_segment_data_to_perioddata(shellmound_sfrdata):
    # shellmound_sfrdata = copy.deepcopy(shellmound_sfrdata)
    shellmound_sfrdata.segment_data.loc[2, 'flow'] = 100
    perdata = shellmound_sfrdata.period_data
    assert isinstance(perdata, pd.DataFrame)
    assert len(perdata) == 1
    assert set(perdata.columns) == \
           {'iseg', 'inflow', 'rno', 'status',
            'ireach', 'icalc', 'per'}


def test_idomain(mf6sfr_instance, shellmound_model):
    ibound = mf6sfr_instance.idomain
    idomain = shellmound_model.dis.idomain.array
    assert np.array_equal(ibound, idomain)


def test_write(shellmound_sfrdata, mf6sfr_instance, outdir):
    mf6sfr = mf6sfr_instance
    outfile = os.path.join(outdir, 'junk.sfr')
    outfile2 = os.path.join(outdir, 'junk2.sfr')
    mf6sfr.write_file(filename=outfile,
                      options=['save_flows',
                               'BUDGET FILEOUT {}.cbc'.format(outfile2),
                               'STAGE FILEOUT {}.stage.bin'.format(outfile2),
                               ]
                      )

    shellmound_sfrdata.write_package(outfile2, version='mf6')
    assert filecmp.cmp(outfile, outfile2)
