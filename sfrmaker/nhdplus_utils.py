import fiona
from fiona.crs import to_string
from .gis import shp2df, crs, project, get_bbox

def load_NHDPlus_v2(NHDFlowlines=None, PlusFlowlineVAA=None, PlusFlow=None, elevslope=None,
                    filter=None,
                    epsg=None, proj4=None, prjfile=None):
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
    # get crs information from flowline projection file
    if prjfile is None:
        prjfile = NHDFlowlines if not isinstance(NHDFlowlines, list) else NHDFlowlines[0]
    nhdcrs = crs(epsg=epsg, proj4=proj4, prjfile=prjfile)

    # ensure that filter bbox is in same crs as flowlines
    if filter is not None and not isinstance(filter, tuple):
        filter = get_bbox(filter, nhdcrs)

    fl_cols = ['COMID',  # 'FCODE', 'FDATE', 'FLOWDIR',
               # 'FTYPE', 'GNIS_ID',
               'GNIS_NAME', #'LENGTHKM',
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

    # add routing information from PlusFlow table;
    df['tocomid'] = get_tocomids(pf, df.index.tolist())
    return df

def get_tocomids(pf, fromcomid_list):
    print('\nGetting routing information from NHDPlus Plusflow table...')
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
    index_col = [c for c in df.columns if c.lower() == index_col]
    if len(index_col) == 0:
        if isinstance(shpfiles, list):
            shpfiles = '\n'.join(shpfiles)
        raise IndexError('No {} column found in: \n{}'.format(index_col,
                                                              shpfiles))
    else:
        df.index = df[index_col[0]]
    return df