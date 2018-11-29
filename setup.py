from setuptools import setup

long_description = \
"""Python package for automating construction of stream flow routing networks from hydrography data.
"""

setup(name="sfrmaker",
      description=long_description,
      long_description=long_description,
      author="Andrew Leaf",
      author_email='aleaf@usgs.gov',
      url='https://github.com/aleaf/SFRmaker/refactor',
      download_url='https://github.com/aleaf/SFRmaker/archive/refactor.zip',
      license='New BSD',
      platforms='Windows, Mac OS-X',
      packages=["sfrmaker"],
      version="0.1")