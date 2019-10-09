#TODO: add unit tests for mf5to6.py
import filecmp
import os
import numpy as np
import pytest
import sfrmaker


@pytest.fixture(scope='module')
def shellmound_ModflowSfr2(shellmound_sfrdata):
    return shellmound_sfrdata.ModflowSfr2


@pytest.fixture(scope='function')
def mf6sfr_instance(shellmound_ModflowSfr2):
    return sfrmaker.mf6sfr(shellmound_ModflowSfr2)


def test_idomain(mf6sfr_instance, shellmound_model):
    ibound = mf6sfr_instance.idomain
    idomain = shellmound_model.dis.idomain.array
    assert np.array_equal(ibound, idomain)


def test_write(shellmound_sfrdata, mf6sfr_instance, outdir):
    mf6sfr = mf6sfr_instance
    outfile = os.path.join(outdir, 'junk.sfr')
    mf6sfr.write_file(filename=outfile)

    outfile2 = os.path.join(outdir, 'junk2.sfr')
    shellmound_sfrdata.write_package(outfile2,
                            version='mf6'
                            )
    assert filecmp.cmp(outfile, outfile2)