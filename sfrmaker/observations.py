"""
Functions for handling observations of SFR package output.
"""
import pandas as pd

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

    if rno_column_in_data is None:
        assert line_id_column_in_data in data.columns, \
            "Data need an id column so observation locations can be mapped to reach numbers"
        # map NHDPlus COMIDs to reach numbers
        if flowline_routing is None:
            rd = sfrdata.reach_data
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

    # add inflows to period_data
    obsdata = pd.DataFrame(columns=sfrd.observations.columns)

    # remove duplicate locations
    data = data.groupby(rno_column_in_data).first().reset_index()
    obsdata['rno'] = data[rno_column_in_data]
    if obstype is not None:
        obsdata['obstype'] = obstype
    else:
        assert obstype_column_in_data in data.columns, "Need to specify observation type."
        obsdata['obstype'] = data[obstype_column_in_data]
    obsdata['obsname'] = data[obsname_column_in_data].astype(str)
    sfrdata._observations = sfrdata.observations.append(obsdata).reset_index(drop=True)
    sfrdata._observations['rno'] = sfrdata.observations.rno.astype(int)
    return obsdata


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
