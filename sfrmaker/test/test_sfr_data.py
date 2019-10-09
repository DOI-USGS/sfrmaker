import os
import pytest
import numpy as np
import pandas as pd
import flopy


def test_const(shellmound_sfrdata, sfr_testdata):
    assert shellmound_sfrdata._lenuni == 2  # meters
    assert shellmound_sfrdata._itmuni == 4  # days
    assert shellmound_sfrdata.const == 86400

    sfr_testdata._model_length_units = 'feet'
    sfr_testdata._model_time_units = 'days'
    assert sfr_testdata.const == 86400 * 1.486


@pytest.fixture(scope="function")
def mf6sfr(shellmound_sfrdata, shellmound_model):
    return shellmound_sfrdata.create_mf6sfr(model=shellmound_model)


@pytest.fixture(scope="function")
def mf6sfr_outfile(mf6sfr, shellmound_model, outdir):
    m = shellmound_model
    mf6sfr.write()
    outfile = os.path.join(m.model_ws, mf6sfr.filename)
    return outfile


def test_create_mf6sfr(mf6sfr, shellmound_sfrdata, shellmound_model):
    assert isinstance(mf6sfr, flopy.mf6.ModflowGwfsfr)
    assert hasattr(shellmound_model, 'sfr')

    # test that packagedata has same information as sfrmaker.sfrdata
    packagedata = pd.DataFrame(mf6sfr.packagedata.array)
    k, i, j = list(zip(*packagedata['cellid'].tolist()))
    packagedata['k'] = k
    packagedata['i'] = i
    packagedata['j'] = j
    packagedata.drop('cellid', axis=1, inplace=True)
    for col in packagedata.columns:
        if packagedata[col].dtype == np.object:
            packagedata[col] = pd.to_numeric(packagedata[col])
    reach_data = shellmound_sfrdata.reach_data
    assert np.array_equal(packagedata['rno'].values + 1, reach_data['rno'].values)
    packagedata['iseg'] = reach_data['iseg']
    segment_data = shellmound_sfrdata.segment_data
    packagedata_groups = packagedata.groupby('iseg').mean()
    for col in packagedata.columns:
        if col in shellmound_sfrdata.mf5names:
            mf2005col = shellmound_sfrdata.mf5names[col]
            if mf2005col in reach_data.columns:
                c1 = packagedata[col].values
                c2 = reach_data[mf2005col].values
                if col == 'rno':
                    c1 += 1
                assert np.allclose(c1, c2, rtol=0.001), \
                    print('values in packagedata.{} != sfrdata.reach_data.{}'.format(col, mf2005col))
            elif mf2005col in segment_data.columns:
                assert np.allclose(packagedata_groups[col].values, segment_data[mf2005col].values, rtol=0.001), \
                    print('values in packagedata.{} != sfrdata.segment_data.{}'.format(col, mf2005col))


def test_flopy_mf6sfr_outfile(mf6sfr, mf6sfr_outfile):
    assert os.path.exists(mf6sfr_outfile)
    # check that the written package data matches in memory package data
    pkdata1 = mf6sfr.packagedata.array.copy()
    mf6sfr.load()
    for col in pkdata1.dtype.names:
        if col == 'cellid':
            continue
        elif pkdata1[col].dtype == np.object:
            c1 = pd.to_numeric(pkdata1[col])
            c2 = pd.to_numeric(mf6sfr.packagedata.array[col])
        else:
            c1 = pkdata1[col]
            c2 = mf6sfr.packagedata.array[col]
        assert np.allclose(c1, c2)


def test_ibound_representation_of_idomain(shellmound_sfrdata, shellmound_model):
    ibound = shellmound_sfrdata.ModflowSfr2.parent.bas6.ibound.array
    idomain = shellmound_model.dis.idomain.array
    assert np.array_equal(ibound, idomain)


