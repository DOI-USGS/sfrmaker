import sys

sys.path.append('..')


def test_crs(datapath):
    from sfrmaker import CRS
    crsobj = CRS(prjfile='{}/badriver/NHDflowlines.shp'.format(datapath))
