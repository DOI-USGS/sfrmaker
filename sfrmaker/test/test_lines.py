import platform
import pytest
from gisutils import get_authority_crs, get_shapefile_crs
import sfrmaker
from sfrmaker.checks import is_to_one
from sfrmaker.nhdplus_utils import load_nhdplus_v2, get_prj_file


@pytest.fixture
def nhdplus_dataframe(datapath):
    pfvaa_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/PlusFlowlineVAA.dbf'.format(datapath)]
    plusflow_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/PlusFlow.dbf'.format(datapath)]
    elevslope_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/elevslope.dbf'.format(datapath)]
    flowlines = ['{}/tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp'.format(datapath)]
    df = load_nhdplus_v2(NHDFlowlines=flowlines, PlusFlowlineVAA=pfvaa_files,
                         PlusFlow=plusflow_files, elevslope=elevslope_files,
                         bbox_filter='{}/tylerforks/grid.shp'.format(datapath),
                         )
    return df


@pytest.mark.parametrize('one_to_many,crs,prjfile', [
    (True, 4269, None),
    (False, None, True),
    # test should fail if the supplied crs and prjfile are inconsistent
    pytest.param(True, 5070, True, marks=pytest.mark.xfail), 
     ])
def test_from_dataframe(nhdplus_dataframe, one_to_many, 
                        crs, prjfile, datapath):
    df = nhdplus_dataframe
    if one_to_many:
        # add some fake distributaries to a few flowlines
        id = df.COMID.values[0]
        inds = df.index.values[:3]
        modified_toids = [[toids[0], id] for toids in df.loc[inds, 'tocomid']] + df.loc[df.index.values[3]:, 'tocomid'].tolist()
        df['tocomid'] = modified_toids
        assert not is_to_one(df.tocomid)
    flowlines = ['{}/tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp'.format(datapath)]
    if prjfile is not None:
        prjfile = get_prj_file(NHDFlowlines=flowlines)

    # convert arbolate sums from km to m
    df['asum2'] = df.ArbolateSu * 1000

    # convert comid end elevations from cm to m
    if 'MAXELEVSMO' in df.columns:
        df['elevup'] = df.MAXELEVSMO / 100.
    if 'MINELEVSMO' in df.columns:
        df['elevdn'] = df.MINELEVSMO / 100.
    lines = sfrmaker.Lines.from_dataframe(df, id_column='COMID',
                                          routing_column='tocomid',
                                          name_column='GNIS_NAME',
                                          asum_units='meters',
                                          elevation_units='meters',
                                          prjfile=prjfile, crs=crs)
    # verify that values (to ids) in original routing dictionary are scalars
    assert is_to_one(lines._original_routing)
    assert lines.df.crs == lines.crs
    if crs is not None:
        assert lines.crs == get_authority_crs(crs)
    elif prjfile is not None:
        assert lines.crs == get_shapefile_crs(prjfile)


def test_write_shapefile(tylerforks_lines_from_NHDPlus, test_data_path):
    lines = tylerforks_lines_from_NHDPlus
    outshp = test_data_path / 'lines.shp'
    lines.write_shapefile(outshp)
    assert outshp.exists()


@pytest.mark.skipif(platform.system() == 'Linux', reason="inscrutable pyogrio.errors.DataSourceError")
def test_load_nhdplus_hr(neversink_lines_from_nhdplus_hr):
    
    lines = neversink_lines_from_nhdplus_hr
    assert isinstance(lines, sfrmaker.lines.Lines)
    assert is_to_one(lines._original_routing)


@pytest.mark.skipif(platform.system() == 'Linux', reason="inscrutable pyogrio.errors.DataSourceError")
@pytest.mark.parametrize('kwargs', (
    {'drop_ftypes': [428], 'drop_NHDPlusIDs': [10000200240966]},
))
def test_load_nhdplus_hr_options(datapath, kwargs):
    NHDPlusHR_paths = [f'{datapath}/neversink_rondout/NHDPLUS_HR_1.gdb', 
                       f'{datapath}/neversink_rondout/NHDPLUS_HR_2.gdb']
    boundary_file = f'{datapath}/neversink_rondout/Model_Extent.shp'

    lns = sfrmaker.Lines.from_nhdplus_hr(NHDPlusHR_paths,
                                        bbox_filter=boundary_file,
                                        **kwargs)
    assert not any(set(lns.df.id).intersection(kwargs['drop_NHDPlusIDs']))
    # these two NHDPlusIDs have FType == 428
    assert not any(set(lns.df.id).intersection([10000700047982, 10000200046339]))


@pytest.mark.parametrize('crs', 
                         [None, 
                          4269, 
                          # should fail of a different CRS is input
                          # (than the projection file)
                          pytest.param(5070, marks=pytest.mark.xfail)
                          ]
                         )
def test_from_nhdplus_v2(datapath, crs):
    pfvaa_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/PlusFlowlineVAA.dbf'.format(datapath)]
    plusflow_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/PlusFlow.dbf'.format(datapath)]
    elevslope_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/elevslope.dbf'.format(datapath)]
    flowlines = ['{}/tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp'.format(datapath)]
    df = load_nhdplus_v2(NHDFlowlines=flowlines, PlusFlowlineVAA=pfvaa_files,
                         PlusFlow=plusflow_files, elevslope=elevslope_files,
                         bbox_filter='{}/tylerforks/grid.shp'.format(datapath),
                         crs=crs
                         )
    assert len(df) > 0
    if crs is not None:
        assert df.crs == get_authority_crs(crs)


@pytest.mark.parametrize('crs,prjfile', [
    (5070, None),
    (None, 'examples/meras/flowlines.shp'),
    # test should fail if the supplied crs and prjfile are inconsistent
    # (this preprocessed shapefile is in 5070)
    pytest.param(4269, 'examples/meras/flowlines.shp', marks=pytest.mark.xfail), 
     ])
def test_from_shapefile(crs, prjfile):
    lines = sfrmaker.Lines.from_shapefile(
        shapefile='examples/meras/flowlines.shp',
        id_column='COMID',  # arguments to sfrmaker.Lines.from_shapefile
        routing_column='tocomid',
        width1_column='width1',
        width2_column='width2',
        up_elevation_column='elevupsmo',
        dn_elevation_column='elevdnsmo',
        name_column='GNIS_NAME',
        width_units='feet',  # units of source data
        elevation_units='feet',  # units of source data
        prjfile=prjfile, crs=crs
        )
    assert lines.df.crs == lines.crs
    if crs is not None:
        assert lines.crs == get_authority_crs(crs)
    elif prjfile is not None:
        assert lines.crs == get_shapefile_crs(prjfile)
        

def test_lines_to_crs(datapath):
    """This test is for the intermittent issue where reprojection 
    with pyproj fails due to an SSL certificate issue. This problem is 
    difficult to reproduce, possibly because pyproj catches data from the PROJ Network,
    so this test may only be meaningful in the situation of a new pyproj install
    on an internal network requiring an SSL certificate bundle. In the 
    situation where pyproj can't access the PROJ Network 
    (resulting in inf values on reprojection), 
    SFRmaker will pyproj.network.set_network_enabled(False) and try again. Therefore,
    reprojection as below shouldn't fail, but the results may potentially be less
    accurate than they would otherwise.
    """
    NHDPlus_paths = datapath / 'tylerforks/NHDPlus'
    lines = sfrmaker.Lines.from_nhdplus_v2(NHDPlus_paths=NHDPlus_paths)
    lines.to_crs(32616)    
        