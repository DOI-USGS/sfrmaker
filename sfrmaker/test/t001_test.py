import sys

sys.path.append('..')


def test_crs(datapath):
    from sfrmaker import CRS
    crsobj = CRS(prjfile='{}/tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp'.format(datapath))
