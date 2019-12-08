import os
import shutil
import glob
import pytest


def included_notebooks():
    include = ['Examples/Notebooks']
    files = []
    for folder in include:
        files += glob.glob(os.path.join(folder, '*.ipynb'))
    return sorted(files)


@pytest.fixture(params=included_notebooks(), scope='module')
def notebook(request):
    return request.param


# even though test runs locally on Windows 10, and on Travis
@pytest.mark.xfail(os.environ.get('APPVEYOR') == 'True',
                   reason="jupyter kernel has timeout issue on appveyor for some reason")
def test_notebook(notebook, outdir):
    # run autotest on each notebook
    path, fname = os.path.split(notebook)
    cmd = ('jupyter ' + 'nbconvert '
           '--ExecutePreprocessor.timeout=600 '
           '--ExecutePreprocessor.kernel_name=test '
           '--to ' + 'notebook '
           '--execute ' + '{} '.format(notebook) +
           '--output-dir ' + '{} '.format(outdir) +
           '--output ' + '{}'.format(fname))
    ival = os.system(cmd)
    assert ival == 0, 'could not run {}'.format(fname)
