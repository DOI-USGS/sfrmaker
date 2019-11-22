import sys

sys.path.append('..')


def test_crs(datapath):
    from sfrmaker import crs
    crsobj = crs(prjfile='{}/badriver/NHDflowlines.shp'.format(datapath))
