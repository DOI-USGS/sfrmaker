import os
import time

import pandas as pd
from gisutils import shp2df
from .gis import crs, get_bbox


def get_prj_file(NHDPlus_paths=None, NHDFlowlines=None):
    if NHDPlus_paths is not None:
        if isinstance(NHDPlus_paths, str):
            NHDPlus_paths = [NHDPlus_paths]
        return os.path.join(NHDPlus_paths[0], 'NHDSnapshot/Hydrography/NHDFlowline.prj')
    elif NHDFlowlines is not None:
        if isinstance(NHDFlowlines, str):
            NHDFlowlines = [NHDFlowlines]
        return NHDFlowlines[0][:-4] + '.prj'


def get_nhdplus_v2_filepaths(NHDPlus_paths):
    print('for basins:')
    if isinstance(NHDPlus_paths, str):
        NHDPlus_paths = [NHDPlus_paths]
    for path in NHDPlus_paths:
        print(path)
    NHDFlowlines = [os.path.join(f, 'NHDSnapshot/Hydrography/NHDFlowline.shp')
                    for f in NHDPlus_paths]
    PlusFlowlineVAA = [os.path.join(f, 'NHDPlusAttributes/PlusFlowlineVAA.dbf')
                       for f in NHDPlus_paths]
    PlusFlow = [os.path.join(f, 'NHDPlusAttributes/PlusFlow.dbf')
                for f in NHDPlus_paths]
    elevslope = [os.path.join(f, 'NHDPlusAttributes/elevslope.dbf')
                 for f in NHDPlus_paths]
    for paths in NHDFlowlines, PlusFlowlineVAA, PlusFlow, elevslope:
        for f in paths:
            if not os.path.exists(f):
                raise FileNotFoundError(f)
    return NHDFlowlines, PlusFlowlineVAA, PlusFlow, elevslope


def get_nhdplus_v2_routing(PlusFlow_file,
                           from_col='FROMCOMID', to_col='TOCOMID'):
    """Read PlusFlow file and return the routing
    information as a dictionary of to:from COMID numbers.
    """
    fname, ext = os.path.splitext(PlusFlow_file)
    if ext in ['.shp', '.dbf']:
        df = shp2df(PlusFlow_file)
    elif ext == '.csv':
        df = pd.read_csv(PlusFlow_file)
    else:
        raise Exception("Unrecognized file-type for PlusFlow table: {}".format(PlusFlow_file))
    flowline_routing = dict(zip(df[from_col], df[to_col]))
    comids = set(df[from_col])
    flowline_routing = {k: v if v in comids else 0
                        for k, v in flowline_routing.items()}
    return flowline_routing


def load_nhdplus_v2(NHDPlus_paths=None,
                    NHDFlowlines=None, PlusFlowlineVAA=None, PlusFlow=None, elevslope=None,
                    filter=None,
                    epsg=None, proj_str=None, prjfile=None):
    """
    Parameters
    ==========
    NHDFlowlines : str or list of strings.
        Shapefile or list of NHDFlowline shapefiles containing
        feature geometries (line arcs) for stream network. Must contain
        the following attribute fields:
        COMID : common identifier number
    PlusFlowlineVAA : str or list of strings.
        DBF file or list of DBF files with NHDPlus attribute information.
        Must contain the following attribute fields:
        COMID : common identifier number
    PlusFlow : str or list of strings.
        DBF file or list of DBF files with NHDPlus routing information.
        Must contain the following attribute fields:
        COMID : common identifier number
    elevslope : str or list of strings.
        DBF file or list of DBF files with end elevations for each
        line arc in NHDFlowlines. Must contain the following attribute fields:
        COMID : common identifier number
    filter : tuple, str (filepath), shapely Polygon or GeoJSON polygon
        Bounding box (tuple) or polygon feature of model stream network area.
        Shapefiles will be reprojected to the CRS of the flowlines; all other
        feature types must be supplied in same CRS as flowlines.
    """
    print("\nloading NHDPlus v2 hydrography data...")
    ta = time.time()

    if NHDPlus_paths is not None:
        NHDFlowlines, PlusFlowlineVAA, PlusFlow, elevslope = \
            get_nhdplus_v2_filepaths(NHDPlus_paths)

    # get crs information from flowline projection file
    if prjfile is None:
        prjfile = get_prj_file(NHDPlus_paths, NHDFlowlines)
    nhdcrs = crs(epsg=epsg, proj_str=proj_str, prjfile=prjfile)

    # ensure that filter bbox is in same crs as flowlines
    # get filters from shapefiles, shapley Polygons or GeoJSON polygons
    if filter is not None and not isinstance(filter, tuple):
        filter = get_bbox(filter, dest_crs=nhdcrs)

    fl_cols = ['COMID',  # 'FCODE', 'FDATE', 'FLOWDIR',
               # 'FTYPE', 'GNIS_ID',
               'GNIS_NAME', 'LENGTHKM',
               # 'REACHCODE', 'RESOLUTION', 'WBAREACOMI',
               'geometry']
    pfvaa_cols = ['ArbolateSu',  # 'Hydroseq', 'DnHydroseq',
                  'StreamOrde',  # 'LevelPathI',
                  ]
    elevs_cols = ['MAXELEVSMO', 'MINELEVSMO']

    # read flowlines and attribute tables into dataframes
    fl = read_nhdplus(NHDFlowlines, bbox_filter=filter)
    pfvaa = read_nhdplus(PlusFlowlineVAA)
    pf = shp2df(PlusFlow)
    elevs = read_nhdplus(elevslope)

    # join flowline and attribute dataframes
    df = fl[fl_cols].copy()
    df = df.join(pfvaa[pfvaa_cols], how='inner')
    df = df.join(elevs[elevs_cols], how='inner')
    print("\nload finished in {:.2f}s".format(time.time() - ta))

    # add routing information from PlusFlow table;
    df['tocomid'] = get_tocomids(pf, df.index.tolist())
    return df


def get_tocomids(pf, fromcomid_list):
    print('\nGetting routing information from NHDPlus Plusflow table...')
    ta = time.time()

    # setup local variables and cull plusflow table to comids in model
    comids = fromcomid_list
    pf = pf.loc[(pf.FROMCOMID.isin(comids)) |
                (pf.TOCOMID.isin(comids))].copy()

    # subset PlusFlow entries for comids that are not in fromcomid_list
    # comids may be missing because they are outside of the model
    # or if the fromcomid_list dataset was edited (resulting in breaks in the routing)
    missing_tocomids = ~pf.TOCOMID.isin(comids) & (pf.TOCOMID != 0)
    missing = pf.loc[missing_tocomids, ['FROMCOMID', 'TOCOMID']].copy()
    # recursively crawl the PlusFlow table
    # to try to find a downstream comid that is in fromcomid_list
    missing['nextCOMID'] = [find_next_comid(tc, pf, comids)
                            for tc in missing.TOCOMID]
    pf.loc[missing_tocomids, 'TOCOMID'] = missing.nextCOMID

    # set any remaining comids not in fromcomid_list to zero
    # (outlets or inlets from outside model)
    pf.loc[~pf.FROMCOMID.isin(comids), 'FROMCOMID'] = 0
    tocomid = pf.TOCOMID.values
    fromcomid = pf.FROMCOMID.values
    tocomids = [tocomid[fromcomid == c].tolist() for c in comids]
    print("finished in {:.2f}s\n".format(time.time() - ta))
    return tocomids


def find_next_comid(comid, pftable, comids, max_levels=10):
    """Crawls the PlusFlow table to find the next downstream comid that
    is in the set comids. Looks up subsequent downstream comids to a
    maximum number of iterations, specified by max_levels (default 10).
    """
    pftable = pftable.copy()
    nextocomid = [comid]
    comids = set(comids)
    for i in range(max_levels):
        nextocomid = pftable.loc[pftable.FROMCOMID.isin(nextocomid), 'TOCOMID'].tolist()
        if len(set(nextocomid).intersection(comids)) > 0:
            # if more than one comid is found, simply take the first
            # (often these will be in different levelpaths,
            # so there is no way to determine a preferred routing path)
            return list(set(nextocomid).intersection(comids))[0]
    return 0


def read_nhdplus(shpfiles, bbox_filter=None,
                 index_col='comid'):
    # read shapefile into dataframe and find the index column
    df = shp2df(shpfiles, filter=bbox_filter)
    if len(df) > 0:
        index_col = [c for c in df.columns if c.lower() == index_col]
        if len(index_col) == 0:
            if isinstance(shpfiles, list):
                shpfiles = '\n'.join(shpfiles)
            raise IndexError('No {} column found in: \n{}'.format(index_col,
                                                                  shpfiles))
        else:
            df.index = df[index_col[0]]
        return df
