import os
from pathlib import Path
import time
import warnings
import pandas as pd
import geopandas as gpd
from gisutils import shp2df, get_shapefile_crs
from .gis import get_bbox, get_crs


def get_prj_file(NHDPlus_paths=None, NHDFlowlines=None):
    if NHDPlus_paths is not None:
        if isinstance(NHDPlus_paths, str) or isinstance(NHDPlus_paths, Path):
            NHDPlus_paths = [NHDPlus_paths]
        return os.path.join(NHDPlus_paths[0], 'NHDSnapshot/Hydrography/NHDFlowline.prj')
    elif NHDFlowlines is not None:
        if isinstance(NHDFlowlines, str) or isinstance(NHDFlowlines, Path):
            NHDFlowlines = [NHDFlowlines]
        return NHDFlowlines[0][:-4] + '.prj'


def get_nhdplus_v2_filepaths(NHDPlus_paths, 
                             raise_not_exist_error=True):
    print('for basins:')
    if isinstance(NHDPlus_paths, str) or isinstance(NHDPlus_paths, Path):
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
            if raise_not_exist_error and not os.path.exists(f):
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
                    filter=None, bbox_filter=None, crs=None,
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
    bbox_filter : tuple, str (filepath), shapely Polygon or GeoJSON polygon
        Bounding box (tuple) or polygon feature of model stream network area.
        Shapefiles will be reprojected to the CRS of the flowlines; all other
        feature types must be supplied in same CRS as flowlines.
    crs : obj
        Coordinate reference system of the NHDPlus data. Only needed if
        the data do not have a valid ESRI projection (.prj) file.
        A Python int, dict, str, or :class:`pyproj.crs.CRS` instance
        passed to :meth:`pyproj.crs.CRS.from_user_input`

        Can be any of:
          - PROJ string
          - Dictionary of PROJ parameters
          - PROJ keyword arguments for parameters
          - JSON string with PROJ parameters
          - CRS WKT string
          - An authority string [i.e. 'epsg:4326']
          - An EPSG integer code [i.e. 4326]
          - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
          - An object with a `to_wkt` method.
          - A :class:`pyproj.crs.CRS` class

        By default, None
    """
    print("\nloading NHDPlus v2 hydrography data...")
    
    if filter is not None:
        raise ValueError("The 'filter' argument is deprecated; use 'bbox_filter' instead")
    
    ta = time.time()

    if NHDPlus_paths is not None:
        NHDFlowlines, PlusFlowlineVAA, PlusFlow, elevslope = \
            get_nhdplus_v2_filepaths(NHDPlus_paths)

    # get crs information from flowline projection file
    #crs_from_shapefile = get_shapefile_crs(NHDFlowlines)
    crs = get_crs(prjfile=NHDFlowlines, epsg=epsg, proj_str=proj_str, crs=crs)

    # ensure that filter bbox is in same crs as flowlines
    # get filters from shapefiles, shapley Polygons or GeoJSON polygons
    if bbox_filter is not None and not isinstance(bbox_filter, tuple):
        bbox_filter = get_bbox(bbox_filter, dest_crs=crs)

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
    fl = read_nhdplus(NHDFlowlines, bbox_filter=bbox_filter)
    pfvaa = read_nhdplus(PlusFlowlineVAA)
    pf = shp2df(PlusFlow)
    elevs = read_nhdplus(elevslope)

    # join flowline and attribute dataframes
    fl.columns = [c.upper() for c in list(fl)]  # added this, switch all to upper case
    fl = fl.rename(columns={"GEOMETRY": "geometry"})  # added this (switch GEOMETRY back to lower case)
    df = fl[fl_cols].copy()
    df = df.join(pfvaa[pfvaa_cols], how='left')
    df = df.join(elevs[elevs_cols], how='left')
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


def read_nhdplus(shapefiles, bbox_filter=None,
                 index_col='comid'):
    # read shapefile into dataframe and find the index column
    if isinstance(shapefiles, str) or isinstance(shapefiles, Path):
        shapefiles = [shapefiles]
    dfs = []
    for i, f in enumerate(shapefiles):
        df = gpd.read_file(f, bbox=bbox_filter)
        if len(dfs) > 0:
            if df.crs != dfs[-1].crs:
                # we could simply reproject, 
                # but if the CRS are different among sets of NHDPlus data, 
                # there might be other issues
                raise ValueError(f'{f} has a different CRS than {shapefiles[i-1]}')
        dfs.append(df)
    df = pd.concat(dfs)
        
    #df = shp2df(shpfiles, filter=bbox_filter)
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


def read_nhdplus_hr(NHDPlusHR_path, bbox_filter=None, drop_fcodes=None):
    ta = time.time()
    print('reading {}...'.format(NHDPlusHR_path))
    #  read NHDFLowlines from NHDPlusHR_path (NHDPlus HR OpenFileGDB)
    fl = gpd.read_file(NHDPlusHR_path, driver='OpenFileGDB', layer='NHDFlowline')
    #  get crs information from flowlines
    fl_crs = fl.crs
    
    if filter is not None:
        print('filtering flowlines...')
    
    #  ensure that filter bbox is in same crs as flowlines
    #  get filters from shapefiles, shapley Polygons or GeoJSON polygons
    if bbox_filter is not None:
        if bbox_filter is not isinstance(bbox_filter, tuple):
            bbox_filter = get_bbox(bbox_filter, dest_crs=fl_crs)
        
        #  filter to bbox using geopandas spatial indexing
        fl = fl.cx[bbox_filter[0]:bbox_filter[2], bbox_filter[1]:bbox_filter[3]]
        
    #  read NHDPlusFlowlineVAA file from NHDPlusHR_path (NHDPlus HR OpenFileGDB) and merge with flowlines
    flvaa = gpd.read_file(NHDPlusHR_path, driver='OpenFileGDB', layer='NHDPlusFlowlineVAA')
    fl = fl.merge(flvaa[['NHDPlusID', 'ArbolateSu','StreamOrde', 'MaxElevSmo', 'MinElevSmo', 'Divergence']],
                  on='NHDPlusID', how='left'
               )
    
    # read NHDPlusFlow file from NHDPlusHR_path (NHDPlus HR OpenFileGDB) 
    pf = gpd.read_file(NHDPlusHR_path, driver='OpenFileGDB', layer='NHDPlusFlow')
    
    #  Remove features classified as minor divergence pathways (Divergence == 2)
    #  from PlusFlow table
    pf_routing_dict = get_hr_routing(pf, fl)
    
    #  Add routing information from PlusFlow table.
    #  Set any remaining comids not in fromcomid_list to zero
    #  (outlets or inlets from outside model)
    fl['ToNHDPID'] = [pf_routing_dict[i] if i in pf_routing_dict else 0.0 for i in fl.NHDPlusID]
    print("finished in {:.2f}s\n".format(time.time() - ta))
    return fl


def get_hr_routing(pf, fl):
    '''
    build NHDPlus HR routing dictionary for connecting FromNHDPID with ToNHDPID
    using NHDPlusFlow and NH
    '''
    print('\nGetting routing information from NHDPlus HR Plusflow table...')
    ta = time.time()
    
    # merge divergence data info to Plusflow dataframe
    pf = pf.merge(fl[['Divergence', 'NHDPlusID']], left_on='ToNHDPID', 
                      right_on = 'NHDPlusID', how='outer')
    pf.rename(columns={'Divergence':'Divergence_ToNHDPID'}, inplace=True)

    # build routing dict excluding Divergece to == 2 (minor divergence path)
    pf_routing_dict = dict(zip(pf.loc[pf.Divergence_ToNHDPID != 2, 'FromNHDPID'], 
                               pf.loc[pf.Divergence_ToNHDPID != 2, 'ToNHDPID']))
    
    print("finished in {:.2f}s\n".format(time.time() - ta))
    return pf_routing_dict


def load_nhdplus_hr(NHDPlusHR_paths, bbox_filter=None, 
                    drop_fcodes=None, drop_ftypes=None, drop_NHDPlusIDs=None):
    """
    Parameters
    ==========
    NHDPlusHR_paths : path (or list of paths) to the NHDPlus High Resolution HU-4 
        Subregion  file geodatabase (.gdb) to include, assuming the file structure 
        is the same as when downloaded from the USGS National Map Downloader tool 
        (v2.0) website (https://apps.nationalmap.gov/downloader/#/).
    filter : tuple, str (filepath), shapely Polygon or GeoJSON polygon
        Bounding box (tuple) or polygon feature of model stream network area.
        Shapefiles will be reprojected to the CRS of the flowlines; all other
        feature types must be supplied in same CRS as flowlines.
    drop_fcodes: int or list of ints, optional
            fcode or list of NHDFlowline FCodes to drop from network. 
    drop_ftypes: int or list of ints, optional
            ftype or list of NHDFlowline Ftypes to drop from network. 
    drop_NHDPlusIDs: int or list of ints, optional
            NHDPlusID or list of NHDFlowline NHDPlusIDs to drop from network. 
    crs : obj
        Coordinate reference system of the NHDPlus data. Only needed if
        the data do not have a valid ESRI projection (.prj) file.
        A Python int, dict, str, or :class:`pyproj.crs.CRS` instance
        passed to :meth:`pyproj.crs.CRS.from_user_input`

        Can be any of:
          - PROJ string
          - Dictionary of PROJ parameters
          - PROJ keyword arguments for parameters
          - JSON string with PROJ parameters
          - CRS WKT string
          - An authority string [i.e. 'epsg:4326']
          - An EPSG integer code [i.e. 4326]
          - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
          - An object with a `to_wkt` method.
          - A :class:`pyproj.crs.CRS` class
        By default, None
        epsg: int, optional
            EPSG code identifying Coordinate Reference System (CRS)
            for features in the input shapefile.
        proj_str: str, optional
            proj_str string identifying CRS for features in the input shapefile.
        prjfile: str, optional
            File path to projection (.prj) file identifying CRS
            for features in the input shapefile. By default,
            the projection file included with the input shapefile
            will be used.

    Returns
    ==========
    df : dataframe
    """
    print("loading NHDPlus HR hydrography data...")
    ta = time.time()
    
    #  Read if using more than one HUC-4 OpenFileGDB or one passed as list
    if isinstance(NHDPlusHR_paths, list):
        dfs = [read_nhdplus_hr(p, bbox_filter=bbox_filter) for p in NHDPlusHR_paths]
        #  make sure that all FLowlines have the same CRS
        assert all(df.crs == dfs[0].crs for df in dfs), 'NHDPlusHR OpenFileGDBs have have different CRSs'
        #  concat into single df
        df = gpd.GeoDataFrame(pd.concat(dfs, ignore_index=True), crs=dfs[0].crs)
        #  get OpenFileGDB crs
        crs = df.crs
    
    #  Read if using one HUC-4 FileGDP passed as str
    if isinstance(NHDPlusHR_paths, str) or isinstance(NHDPlusHR_paths, Path):
        df = read_nhdplus_hr(NHDPlusHR_paths, filter = filter)
        #  get OpenFileGDB crs
        crs = df.crs
   
    #  Option to drop specified FCodes
    if drop_fcodes is not None:    
        df = df.loc[~df.FCode.isin(drop_fcodes)]
    if drop_ftypes is not None:
        df = df.loc[~df.FType.isin(drop_ftypes)]
    if drop_NHDPlusIDs is not None:
        df = df.loc[~df.NHDPlusID.isin(drop_NHDPlusIDs)]
        
    keep_cols = ['NHDPlusID', 'ToNHDPID', 'ArbolateSu',
               'geometry', 'StreamOrde',
               'MaxElevSmo', 'MinElevSmo', 'GNIS_Name']
    
    #  Make final dataframe with only the info needed for lines
    df = df[keep_cols].copy()

    print("\nload finished in {:.2f}s\n".format(time.time() - ta))

    return df, crs
