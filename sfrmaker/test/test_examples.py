import os
from subprocess import Popen, PIPE, STDOUT
import pytest
import sfrmaker
from sfrmaker.fileio import read_mf6_block


@pytest.fixture(autouse=True)
def keep_cwd():
    """Reset the working directory after a test.
    """
    wd = os.getcwd()
    yield wd  # provide the fixture value
    print("reverting working directory from {} to {}".format(os.getcwd(), wd))
    os.chdir(wd)


@pytest.mark.slow
@pytest.mark.parametrize('script', ('examples/meras/make_sfr.py',))
def test_meras_example(script):
    """Basic test of MERAS example, which takes awhile to execute.
    """
    path, scriptname = os.path.split(script)
    os.chdir(path)
    sfrdata = sfrmaker.SFRData.from_yaml('meras_sfrmaker_config.yml')
    #ival = os.system('python {}'.format(scriptname))
    #assert ival == 0, 'could not run {}'.format(scriptname)
