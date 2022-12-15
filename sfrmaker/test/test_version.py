"""
check for consistancy in package version
"""
import os
from packaging import version
import pandas as pd
import sfrmaker
import pytest


def get_readme_version(project_root_path):
    readme_file = os.path.join(project_root_path, 'Readme.md')
    with open(readme_file) as src:
        for line in src:
            if '# version' in line.lower():
                version_info = version.parse(line.strip().split()[2])
                break
    return version_info


def get_changelog_version(project_root_path):
    file = os.path.join(project_root_path, 'docs/source/release-history.rst')
    with open(file) as src:
        for line in src:
            if line.strip().lower().startswith('version'):
                _, version_info, date = line.strip().split()
                version_info = version.parse(version_info)
                date = date.strip('()')
                try:
                    pd.Timestamp(date)
                    break
                except:
                    continue
    return version_info


@pytest.mark.skipif(os.environ.get('GITHUB_ACTIONS') == 'true',
                    reason=("Clone depth of 1 doesn't allow last "
                           "release number to be detected."))
def test_version(project_root_path):
    version_info = version.parse(sfrmaker.__version__)
    readme_version = get_readme_version(project_root_path)
    changelog_version = get_changelog_version(project_root_path)
    assert version_info.base_version == changelog_version.base_version
    assert readme_version.major == version_info.major
    assert readme_version.minor == version_info.minor
