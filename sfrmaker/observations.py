"""
Functions for handling observations of SFR package output.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from gisutils import shp2df
try:
    import flopy
    fm = flopy.modflow
except:
    flopy = False
from .gis import get_shapefile_crs, project
from .fileio import read_tables
from .routing import get_next_id_in_subset


def add_observations(sfrdata, data, flowline_routing=None,
                     obstype=None, sfrlines_shapefile=None,
                     rno_column_in_sfrlines='rno',
                     x_location_column=None,
                     y_location_column=None,
                     line_id_column=None,
                     rno_column=None,
                     obstype_column=None,
                     obsname_column='site_no'):
    """Add SFR observations to the observations DataFrame
    attribute of an sfrdata instance. Observations can
    by located on the SFR network by specifying reach number
    directly (rno_column), by x, y location (x_column_in_data and y_column in data),
    or by specifying the source hydrography lines that they are located on
    (line_id_column).

    Parameters
    ----------
    sfrdata : sfrmaker.SFRData instance
        SFRData instance with reach_data table attribute. To add observations from x, y coordinates,
        the reach_data table must have a geometry column with LineStrings representing each reach, or
        an sfrlines_shapefile is required. Reach numbers are assumed to be in an 'rno' column.
    data : DataFrame, path to csv file, or list of DataFrames or file paths
        Table with information on the observation sites to be located. Must have
        either reach numbers (rno_column), line_ids (line_id_column),
        or x and y locations (x_column_in_data and y_column_in_data).
    obstype : str or list-like (optional)
        Type(s) of observation to record, for MODFLOW-6 (default 'downstream-flow'; see
        MODFLOW-6 IO documentation for more details). Alternatively, observation
        types can be specified by row in data, using the obstype_column_in_data argument.
    x_location_column : str (optional)
        Column in data with site x-coordinates (in same CRS as SFR network).
    y_location_column : str (optional)
        Column in data with site y-coordinates (in same CRS as SFR network).
    sfrlines_shapefile : str (optional)
        Shapefile version of SFRdata.reach_data. Only needed if SFRdata.reach_data doesn't
        have LineString geometries for the reaches.
    rno_column_in_sfrlines : str (optional)
        Column in sfrlines with reach numbers for matching lines with reaches in sfrdata, or
        reach numbers assigned to observation sites. (default 'rno')
    line_id_column : str
        Column in data matching observation sites to line_ids in the source hydrography data.
    rno_column : str
        Column in data matching observation sites to reach numbers in the SFR network.
    flowline_routing : dict
        Optional dictionary of routing for source hydrography. Only needed
        if locating by line_id, and SFR network is a subset of the full source
        hydrography (i.e. some lines were dropped in the creation of the SFR packge,
        or if the sites are inflow points corresponding to lines outside of the model perimeter).
        In this case, observation points referenced to line_ids that are missing from the SFR
        network are placed at the first reach corresponding to the next downstream line_id
        that is represented in the SFR network.
    obstype_column : str (optional)
        Column in data with MODFLOW-6 observation types. For adding observations of different types.
        If obstype and obstype_column_in_data are none, the default of 'downstream-flow' will be used.
    obsname_column : str
        Column in data with unique identifier (e.g. site number or name) for observation sites.


    Notes
    -----
    Sites located by line_id (source hydrography) will be assigned to the last reach in the
    segment corresponding to the line_id. Locating by x, y or reach number is more accurate.

    """
    sfrd = sfrdata
    reach_data = sfrdata.reach_data.copy()

    # allow input via a list of tables or single table
    data = read_tables(data, dtype={obsname_column: object})#, dtype={obsname_column: object})
    # need a special case if allowing obsname_column to also be identifier
    if obsname_column == rno_column:
        obsname_column = f'obsnamecol_{obsname_column}'
        data[obsname_column] = data[rno_column].astype(object)
    elif obsname_column == line_id_column:
        obsname_column = f'obsnamecol_{obsname_column}'
        data[obsname_column] = data[line_id_column].astype(object)
    else:
        data[obsname_column] = data[obsname_column].astype(object)
    assert data[obsname_column].dtype == object

    # read reach geometries from a shapefile
    if sfrlines_shapefile is not None:
        sfrlines = shp2df(sfrlines_shapefile)
        geoms = dict(zip(sfrlines[rno_column_in_sfrlines], sfrlines['geometry']))
        reach_data['geometry'] = [geoms[rno] for rno in reach_data['rno']]

    # if no reach number is provided
    msg = "Observation sites need reach number, (x,y) coordinates, or source hydrography IDs"
    if rno_column not in data.columns:

        rno_column = 'rno'

        # get reach numbers by x, y location of sites
        if x_location_column in data.columns and y_location_column in data.columns:
            locs = locate_sites(data, reach_data,
                                x_column_in_data=x_location_column,
                                y_column_in_data=y_location_column,
                                reach_id_col='rno',  # reach number column in reach_data
                                site_number_col=obsname_column
                                )
            data[rno_column] = locs['rno']

        # get reach number from site locations in source hydrography (line_ids)
        elif line_id_column in data.columns:
            data[line_id_column] = data[line_id_column].astype(int).astype(str)
            # map NHDPlus COMIDs to reach numbers
            if flowline_routing is None:
                line_id = dict(zip(reach_data.iseg, reach_data.line_id))
                sfr_routing = sfrdata.segment_routing.copy()

                # routing for source hydrography
                flowline_routing = {line_id.get(k, 0): line_id.get(v, 0)
                                    for k, v in sfr_routing.items()}
            # get the last reach in each segment
            r1 = reach_data.sort_values(by=['iseg', 'ireach'], axis=0).groupby('iseg').last()
            line_id_rno_mapping = dict(zip(r1['line_id'], r1['rno']))
            line_ids = get_next_id_in_subset(r1.line_id, flowline_routing,
                                             data[line_id_column])
            data[rno_column] = [line_id_rno_mapping[lid] for lid in line_ids]

        else:
            raise ValueError(msg)

    # create observations dataframe
    obsdata = pd.DataFrame(columns=sfrd.observations.columns)

    # remove duplicate locations
    data = data.groupby(rno_column).first().reset_index()
    obsdata['rno'] = data[rno_column]

    # segment and reach info
    iseg_ireach = dict(list(zip(reach_data.rno, zip(reach_data.iseg, reach_data.ireach))))
    obsdata['iseg'] = [iseg_ireach[rno][0] for rno in obsdata.rno]
    obsdata['ireach'] = [iseg_ireach[rno][1] for rno in obsdata.rno]
    for col in ['rno', 'iseg', 'ireach']:
        obsdata[col] = obsdata[col].astype(int)

    if obstype is not None:
        if isinstance(obstype, str):
            obsdata['obstype'] = obstype
            obsdata['obsname'] = data[obsname_column].astype(str)
        else:
            obstypes = obstype
            dfs = []
            for obstype in obstypes:
                df = obsdata.copy()
                df['obstype'] = obstype
                obsnme_suffix = obstype.split('-')[-1]
                df['obsname'] = [f"{obsnme}-{obsnme_suffix}" 
                                 for obsnme in data[obsname_column].astype(str)]
                dfs.append(df)
            obsdata = pd.concat(dfs)
    elif obstype_column in data.columns:
        obsdata['obstype'] = data[obstype_column]
        # check if base observation names (with just site info) are unique
        # if not, formulate name based on type as well
        site_names = data[obsname_column].astype(str)
        if len(set(site_names)) < len(data):
            obsname_suffixes = [obstype.split('-')[-1] for obstype in obsdata['obstype']]
            obsdata['obsname'] = [f"{obsnme}-{obsnme_suffix}" for obsnme, obsnme_suffix 
                                  in zip(site_names, obsname_suffixes)]
        else:
            obsdata['obsname'] = site_names
    else:
        obsdata['obstype'] = 'downstream-flow'
        obsdata['obsname'] = data[obsname_column].astype(str)

    return obsdata


def get_closest_reach(x, y, sfrlines,
                      rno_column='rno'):
    """Get the SFR reach number closest to a point feature.

    Parameters
    ----------
    x : scalar or list of scalars
        x-coordinate(s) of point feature(s)
    y : scalar or list of scalars
        y-coordinate(s) or point feature(s)
    sfrlines: dataframe
        DataFrame containing a geometry column with SFR line arcs,
        and a column rno_column with unique numbers for each reach.
    rno_column: str
        Column with unique number for each reach. default "rno"
    threshold : numeric
        Distance threshold (in CRS units). Only return reaches within
        this distance.

    Returns
    -------
    rno : int or list of ints
        Reach numbers for reaches closest to each location
        defined by x, y.
    distance: int or list of ints
        Distance from each site in x, y to the corresponding 
        closest SFR reach.
    """
    scalar = False
    if np.isscalar(x):
        scalar = True
        x = [x]
    if np.isscalar(y):
        y = [y]

    geoms = sfrlines.geometry.values
    rno = sfrlines[rno_column].values

    allX = []  # all x coordinates in sfrlines
    allY = []  # all y coordinates in sfrlines
    all_rno = []  # reach number for each x, y point in sfrlines
    for i, g in enumerate(geoms):
        if 'Multi' not in g.type:
            g = [g]
        for part in g:
            gx, gy = part.coords.xy
            allX += gx
            allY += gy
            all_rno += [rno[i]]*len(gx)
    allX = np.array(allX)
    allY = np.array(allY)
    all_rno = np.array(all_rno)
    assert len(all_rno) == len(allX)

    rno = []
    distance = []
    for i in range(len(x)):
        distances = np.sqrt((allX - x[i]) ** 2 +
                            (allY - y[i]) ** 2)
        idx = np.argmin(distances)
        rno.append(all_rno[idx])
        distance.append(np.min(distances))
    if scalar:
        return rno[0], distance[0]
    else:
        return rno, distance


def locate_sites(site_data,
                 reach_data,
                 active_area_shapefile=None,
                 x_column_in_data=None,
                 y_column_in_data=None,
                 reach_id_col='rno',
                 ireach_col='ireach',
                 iseg_col='iseg',
                 site_number_col='site_no',
                 keep_columns=None,
                 perimeter_buffer=1000,
                 distance_threshold=1000
                 ):
    """Get SFR reach locations corresponding to x, y points
    (e.g. measurement site locations).

    Parameters
    ----------
    site_data: (Geo)DataFrame or shapefile
        DataFrame or shapefile with point locations and attribute data for
        stream flow observation sites. Point locations can be specified
        in a DataFrame by either x_column_in_data and y_column_in_data, or
        a 'geometry' column of shapely points. If shapefiles or GeoDataFrames are provided
        for both site_data and reach_data, they can be in any CRS, but both must have valid
        CRS references (.prj file or .crs attribute); otherwise site_data and
        reach_data are assumed to be in the same CRS.
    reach_data: (Geo)DataFrame or shapefile
        SFRData.reach_data DataFrame, or shapefile equivalent
        with line-arcs representing all segments and/or reaches.
        If shapefiles or GeoDataFrames are provided
        for both site_data and reach_data, they can be in any CRS, but both must have valid
        CRS references (.prj file or .crs attribute); otherwise site_data and
        reach_data are assumed to be in the same CRS.
    active_area_shapefile: ESRI shapefile or shapely polygon (optional)
        Shapefile or polygon, in same CRS as sfr_lines_shapefile,
        defining areal extent (perimeter) of SFR network.
    x_column_in_data : str (optional)
        Column in data with site x-coordinates (in same CRS as SFR network).
    y_column_in_data : str (optional)
        Column in data with site y-coordinates (in same CRS as SFR network).
    reach_id_col: str
        Column with 1-based unique number for each stream line-arc. default "rno"
    ireach_col: str
        Column with 1-based reach number within segment (i.e. in the modflow-2005 SFR data model),
        by default, 'ireach'
    iseg_col: str
        Column with 1-based reach segment number (i.e. in the modflow-2005 SFR data model),
        by default, 'iseg'
    site_number_col : str
        Name of column in sites_shapefile with number identifying each
        site to be located. default "site_no"
    keep_columns: list of strings
        List of columns in sites_shapefile to retain when
        writing output_csv_file and output_shape_file.
    perimeter_buffer : scalar
        Exclude flows within this distance of perimeter defined
        by active_area. For example, a value of 1000 would
        mean that sites must be at least 1 km inside of the active area perimeter to
        be included.
    distance_threshold : scalar
        Only consider sites within this distance of a stream line-arc.


    Returns
    -------
    locs : DataFrame

    """
    sfr_crs = None
    locs_crs = None
    # read in sfr lines
    if not isinstance(reach_data, pd.DataFrame):
        sfrlines = gpd.read_file(reach_data)
        sfr_crs = sfrlines.crs
    elif isinstance(reach_data, pd.DataFrame):
        sfrlines = reach_data.copy()
        sfr_crs = getattr(sfrlines, 'crs', None)
    else:
        raise TypeError('Datatype for reach_data not understood: {}'.format(reach_data))
    if reach_id_col not in sfrlines.columns:
        sfrlines.index = np.arange(len(sfrlines)) + 1
        sfrlines[reach_id_col] = np.arange(len(sfrlines)) + 1
    else:
        sfrlines.index = sfrlines[reach_id_col]

    # sites to locate
    if not isinstance(site_data, pd.DataFrame):
        locs = gpd.read_file(site_data)
        if isinstance(site_data, list):
            locs_crs = get_shapefile_crs(site_data[0])
        else:
            locs_crs = get_shapefile_crs(site_data)
        locs['site_no'] = locs[site_number_col]  # str_ids(locs.site_no)
    elif isinstance(site_data, pd.DataFrame):
        locs = site_data.copy()
        locs_crs = getattr(locs, 'crs', None)
    else:
        raise TypeError('Datatype for site_data not understood: {}'.format(site_data))

    # to_crs if crs are available
    if locs_crs is not None and sfr_crs is not None:
        locs['geometry'] = project(locs.geometry.values, locs_crs, sfr_crs)

    # get the x and y coordinates
    if x_column_in_data is not None and y_column_in_data is not None:
        x = locs[x_column_in_data]
        y = locs[y_column_in_data]
    else:
        x = [p.x for p in locs.geometry]
        y = [p.y for p in locs.geometry]

    ids, distances = get_closest_reach(x, y, sfrlines,
                                       rno_column=reach_id_col)
    reach_id_col = reach_id_col.lower()
    locs[reach_id_col] = ids
    locs['distance'] = distances
    if iseg_col in sfrlines.columns:
        locs[iseg_col] = sfrlines.loc[ids, iseg_col].values
    if ireach_col in sfrlines.columns:
        locs[ireach_col] = sfrlines.loc[ids, ireach_col].values
    locs = locs.loc[locs['distance'] <= distance_threshold]

    # cull observations at or outside of model perimeter
    # to only those along model perimeter
    if active_area_shapefile is not None:
        active_area = active_area_shapefile
        if not isinstance(active_area_shapefile, Polygon):
            active_area = gpd.read_file(active_area_shapefile).geometry[0]
        perimeter = active_area.exterior.buffer(perimeter_buffer)
        perimeter_inside_buffer = Polygon(perimeter.interiors[0])

        keep = []
        for rn in locs[reach_id_col]:
            geom = sfrlines.loc[rn, 'geometry']
            keep.append(geom.within(perimeter_inside_buffer))
    else:
        keep = slice(None)

    if keep_columns is None:
        keep_columns = locs.columns.tolist()
    for c in [reach_id_col, iseg_col, ireach_col, 'geometry']:
        if c not in keep_columns and c in locs.columns:
            keep_columns.append(c)

    locs = locs.loc[keep, keep_columns]
    return locs


def write_gage_package(location_data,
                       gage_package_filename=None,
                       gage_namfile_entries_file=None,
                       model=None,
                       sitename_col='site_no',
                       gage_package_unit=25,
                       start_gage_unit=200):
    """

    Parameters
    ----------
    location_data : pandas.DataFrame
        Table of observation locations. Must have columns:
        'iseg': segment number
        'ireach': reach number
        sitename_col: specified by sitename_col argument (default 'site_no')
    gage_package_filename : str or pathlike
        Name for the gage package file, which will be written to the 
        model workspace (``model.model_ws``) of the supplied flopy model instance. 
        This must also be manually added to the MODFLOW .nam file, for example::
        
            GAGE          25 br_trans.gage
            
        Or, a new namefile can be written from the supplied flopy model instance.
    gage_namfile_entries_file : str or pathlike
        Namefile entries for the gage package output files will be written
        to this file. The contents of this file must then be manually
        copied and pasted into the MODFLOW .nam file.
        By default, None, in which case this file is named after gage_package_filename.
    model : flopy model instance
        Used for writing the gage package input file via flopy.
    sitename_col : str
        Unique name or number for each gage site.
        By default, 'site_no'
    gage_package_unit : int
        MODFLOW unit number to assign to gage package.
        By default, 25
    start_gage_unit : int
        Modflow unit numbers for each gage package output file 
        will be assigned sequentially starting at this number.
        By default, 200

    Returns
    -------

    """
    if gage_package_filename is None:
        if model is not None:
            gage_package_filename = '{}.gage'.format(model.modelname)
        else:
            gage_package_filename = 'observations.gage'
    if gage_namfile_entries_file is None:
        gage_namfile_entries_file = gage_package_filename + '.namefile_entries'
    if model is not None:
        gage_package_filename = os.path.join(model.model_ws,
                                             os.path.split(gage_package_filename)[1])
        gage_namfile_entries_file = os.path.join(model.model_ws,
                                                 os.path.split(gage_namfile_entries_file)[1])
    # read in stream gage observation locations from locate_flux_targets_in_SFR.py
    df = location_data.copy()
    df.sort_values(by=[sitename_col], inplace=True)
    df['gagefile'] = [f'{sitename}.ggo' for sitename in df[sitename_col]]

    if model is not None and flopy:
        # create flopy gage package object
        gage_data = fm.ModflowGage.get_empty(ncells=len(df))
        gage_data['gageloc'] = df['iseg']
        gage_data['gagerch'] = df['ireach']
        gage_data['unit'] = np.arange(start_gage_unit, start_gage_unit + len(df))
        gag = fm.ModflowGage(model, numgage=len(gage_data),
                             gage_data=gage_data, files=df.gagefile.tolist(),
                             unitnumber=gage_package_unit
                             )
        gag.fn_path = str(gage_package_filename)
        gag.write_file()
    else:
        raise NotImplementedError('writing a gage package without Flopy')

    with open(gage_namfile_entries_file, 'w') as output:
        for i, f in enumerate(df.gagefile.values):
            output.write('DATA  {:d} {}\n'.format(gage_data.unit[i], f))
    return gag


def write_mf6_sfr_obsfile(observation_locations,
                          filename,
                          sfr_output_filename, digits=5,
                          print_input=True):
    """Write MODFLOW-6 observation input for the SFR package.

    Parameters
    ----------
    observation_locations : DataFrame or str (filepath to csv file)
        line_id : unique identifier for observation locations
    filename : str
        File path for MODFLOW-6 SFR observation input file
    sfr_output_filename : str
        File path that SFR observation output file
    digits : int
        the number of significant digits with which simulated values
        are written to the output file.
    print_input : bool
        keyword to indicate that the list of observation information
        will be written to the listing file immediately after it is read.

    Returns
    -------
    writes filename
    """
    locs = observation_locations
    if isinstance(observation_locations, str) or isinstance(observation_locations, Path):
        locs = pd.read_csv(observation_locations)

    with open(filename, 'w') as output:
        output.write('BEGIN OPTIONS\n  DIGITS {:d}\n'.format(digits))
        if print_input:
            output.write('  PRINT_INPUT\n')
        output.write('END OPTIONS\n')
        output.write('BEGIN CONTINUOUS FILEOUT {}\n'.format(sfr_output_filename))
        output.write('# obsname  obstype  rno\n')
        for i, r in locs.iterrows():
            output.write('  {}  {}  {:d}\n'.format(r.obsname, r.obstype, r.rno))
        output.write('END CONTINUOUS\n')
    print('wrote {}'.format(filename))
