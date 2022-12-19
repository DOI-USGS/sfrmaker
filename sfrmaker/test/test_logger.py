"""
Tests for the SFRmaker logger module.
"""
import os
from pathlib import Path
import pytest
from sfrmaker.logger import Logger


@pytest.fixture(autouse=True)
def outfolder(outdir):
    wd = os.getcwd()
    folder = os.path.join(outdir, 'logging')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    os.chdir(folder)
    yield folder
    os.chdir(wd)


@pytest.fixture()
def existing_log_file(outfolder):
    logfile = os.path.join(outfolder, 'sfrmaker.logger')
    with open(logfile, 'w') as dest:
        dest.write('existing logger file\n')
    yield dest


@pytest.mark.parametrize('filename,mode', ((None, 'w'),
                                           ('sfrmaker.logger', 'w'),  # new logger file
                                           ('sfrmaker.logger', 'a'),  # append to existing logger file
                                           (existing_log_file, 'junk')  # work with an open file handle
                                      ))
def test_init(filename, mode, existing_log_file, outfolder):
    if isinstance(filename, str) or isinstance(filename, Path):
        filename = os.path.join(outfolder, filename)
        logger = Logger(filename, mode)
    elif filename is not None:
        existing_log_file = open(existing_log_file.name, 'a')
        logger = Logger(existing_log_file)
    else:
        logger = Logger()
    assert Path(logger.filename).resolve() == Path(outfolder) / 'sfrmaker.logger'
    assert not logger.f.closed
    logger.f.close()
    with open(logger.filename) as src:
        lines = src.readlines()
    if mode == 'w':
        assert lines[0].split()[0] == 'SFRmaker'
    else:
        assert lines[0] == 'existing logger file\n'
    j=2