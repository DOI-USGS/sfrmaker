import sys
sys.path.append('..')


def test_crs():

    from sfrmaker import crs
    crsobj = crs(prjfile='data/WI_flowlines.prj')


if __name__ == '__main__':
    test_crs()