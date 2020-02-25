"""
Functions for handling observations of SFR package output.
"""
import os
import numpy as np
import pandas as pd
try:
    import flopy
    fm = flopy.modflow
except:
    flopy = False

from .routing import get_next_id_in_subset


def add_observations(sfrdata, data, flowline_routing=None,
                     obstype=None,
                     line_id_column_in_data=None,
                     rno_column_in_data=None,
                     obstype_column_in_data=None,
                     obsname_column_in_data='site_no'):
    """Add SFR observations to the observations DataFrame
    attribute of an sfrdata instance.

    """
    sfrd = sfrdata
    if isinstance(data, str):
        data = pd.read_csv(data)
    elif isinstance(data, pd.DataFrame):
        data = data.copy()
    else:
        raise Exception('Unrecognized input type for data:\n{}'.format(data))

    rd = sfrdata.reach_data
    if rno_column_in_data is None:
        assert line_id_column_in_data in data.columns, \
            "Data need an id column so observation locations can be mapped to reach numbers"
        # map NHDPlus COMIDs to reach numbers
        if flowline_routing is None:
            line_id = dict(zip(rd.iseg, rd.line_id))
            sfr_routing = sfrdata.segment_routing.copy()

            # routing for source hydrography
            flowline_routing = {line_id.get(k, 0): line_id.get(v, 0)
                                for k, v in sfr_routing.items()}
        rno_column_in_data = 'rno'
        r1 = sfrd.reach_data.loc[sfrd.reach_data.ireach == 1]
        line_id_rno_mapping = dict(zip(r1['line_id'], r1['rno']))
        line_ids = get_next_id_in_subset(r1.line_id, flowline_routing,
                                         data[line_id_column_in_data])
        data[rno_column_in_data] = [line_id_rno_mapping[lid] for lid in line_ids]
    else:
        assert rno_column_in_data in data.columns, \
            "Data to add need reach number, or flowline routing information is needed."

    # create observations dataframe
    obsdata = pd.DataFrame(columns=sfrd.observations.columns)

    # remove duplicate locations
    data = data.groupby(rno_column_in_data).first().reset_index()
    obsdata['rno'] = data[rno_column_in_data]

    # segment and reach info
    iseg_ireach = dict(list(zip(rd.rno, zip(rd.iseg, rd.ireach))))
    obsdata['iseg'] = [iseg_ireach[rno][0] for rno in obsdata.rno]
    obsdata['ireach'] = [iseg_ireach[rno][1] for rno in obsdata.rno]
    for col in ['rno', 'iseg', 'ireach']:
        obsdata[col] = obsdata[col].astype(int)

    if obstype is not None:
        obsdata['obstype'] = obstype
    elif obstype_column_in_data in data.columns:
        obsdata['obstype'] = data[obstype_column_in_data]
    else:
        obsdata['obstype'] = 'downstream-flow'
    obsdata['obsname'] = data[obsname_column_in_data].astype(str)
    sfrdata._observations = sfrdata.observations.append(obsdata).reset_index(drop=True)
    # enforce dtypes (pandas doesn't allow an empty dataframe to be initialized with more than one specified dtype)
    for df in [obsdata, sfrdata._observations]:
        for col in ['rno', 'iseg', 'ireach']:
            df[col] = df[col].astype(int)
    return obsdata


def write_gage_package(location_data,
                       gage_package_filename=None,
                       gage_namfile_entries_file=None,
                       model=None,
                       obsname_col='obsname',
                       gage_package_unit=25,
                       start_gage_unit=228):
    """

    Parameters
    ----------
    location_data : pandas.DataFrame
        Table of observation locations. Must have columns:
        'iseg': segment number
        'ireach': reach number
        obsname_col: specified by obsname_col argument (default 'obsname')
    gage_package_filename :
    gage_namfile_entries_file :
    model :
    obsname_col :
    gage_package_unit :
    start_gage_unit :

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
    df.sort_values(by=[obsname_col], inplace=True)
    df['gagefile'] = ['{}.ggo'.format(obsname) for obsname in df.obsname]

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
        gag.fn_path = gage_package_filename
        gag.write_file()
    else:
        j=2

    with open(gage_namfile_entries_file, 'w') as output:
        for i, f in enumerate(df.gagefile.values):
            output.write('DATA  {:d} {}\n'.format(gage_data.unit[i], f))


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
    if isinstance(observation_locations, str):
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
