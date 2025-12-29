"""Code for converting MODFLOW-2005 style SFR Package input to MODFLOW-6 style SFR package input.
"""
import os
from copy import copy
import warnings
import numpy as np
import pandas as pd
import sfrmaker
from sfrmaker.reaches import interpolate_to_reaches


class Mf6SFR:
    """Class for writing MODFLOW-6 SFR package input
    from a (MODFLOW-2005 style) flopy.modflow.ModflowSfr2 instance or
    an sfrmaker.SFRData instance.

    Parameters
    ----------
    ModflowSfr2 : flopy.modflow.ModflowSfr2 instance, optional
        Input SFR dataset
    SFRData : sfrmaker.SFRData instance, optional
        Input SFR dataset
    period_data : DataFrame, optional
        DataFrame of MODFLOW-6-style stress period data, as made by
        :meth:`SFRData.add_to_perioddata`. Only needed if SFRData isn't supplied,
        by default None
    idomain : ndarray, optional
        3D numpy array designating active cells (idomain==1).
        SFR reaches in inactive cells will be written with 'none' in the cellid field.
        by default None
    options : list, optional
        List of strings to write to the MODFLOW-6 SFR options block. For example::

                options=['save_flows',
                         'BUDGET FILEOUT model.sfr.cbc',
                         'STAGE FILEOUT model.sfr.stage.bin']

        An appropriate unit_conversion is written by default.
        See MODFLOW-6 documentation for other options.
        By default None.

    auxiliary_line_numbers : bool, optional
        If true, add 'line_id' as an auxiliary variable to the options block
        and write hydrography line IDs to the packagedata block in the auxiliary
        'line_id' column, by default True.
    """
    # convert from ModflowSfr to mf6
    mf6names = {'rno': 'rno',
                'node': 'cellid',
                'rchlen': 'rlen',
                'slope': 'rgrd',
                'strtop': 'rtp',
                'strthick': 'rbth',
                'strhc1': 'rhk',
                'roughch': 'man',
                'flow': 'inflow',
                'pptsw': 'rainfall',
                'etsw': 'evaporation',
                'runoff': 'runoff',
                'depth1': 'depth1',  # need these for distributing stage
                'depth2': 'depth2'}
    mf5names = {v: k for k, v in mf6names.items()}

    cols = ['rno', 'cellid', 'k', 'i', 'j', 'rlen', 'rwid', 'rgrd',
            'rtp', 'rbth', 'rhk', 'man', 'ncon', 'ustrf', 'ndv', 'idomain', 'line_id']
    def __init__(self, ModflowSfr2=None, SFRData=None,
                 period_data=None, idomain=None,
                 options=None, auxiliary_line_numbers=True):

        # instantiate with SFRData instance instead of ModflowSfr2 instance
        # allows auxiliary variables from SFRData.reach_data
        # (that aren't allowed in ModflowSfr2.reach_data) to be written
        if SFRData is not None:
            ModflowSfr2 = SFRData.modflow_sfr2

        # check for other packages
        if idomain is None:
            try:
                self.idomain = ModflowSfr2.parent.bas6.ibound.array.copy()
                # constant head cells (-1) also made active
                self.idomain[self.idomain != 0] = 1
            except:
                txt = 'Warning: BAS6 package not found. '
                txt += 'Cannot check for reaches in inactive cells. '
                txt += 'Converted SFR package may not run with MODFLOW 6.'
                print(txt)
                self.idomain = None
        else:
            self.idomain = idomain

        # copy modflow_sfr2 object and enforce sorting
        self.ModflowSfr2 = copy(ModflowSfr2)
        self.ModflowSfr2.segment_data[0].sort(order='nseg')
        self.ModflowSfr2.reach_data.sort(order=['iseg', 'ireach'])

        self.structured = self.ModflowSfr2.parent.structured
        self.unit_conversion = ModflowSfr2.const
        self.nreaches = len(ModflowSfr2.reach_data)
        self.nper = ModflowSfr2.nper

        # mf6 options block
        self.auxiliary_line_numbers = auxiliary_line_numbers
        self.options_block = options

        # dataframes of reach and segment data from modflow_sfr2
        if SFRData is not None:
            self.rd = SFRData.reach_data
            self.sd = SFRData.segment_data
            self._period_data = SFRData.period_data
        else:
            self.rd = pd.DataFrame(ModflowSfr2.reach_data)
            self.rd.rename(columns={'reachID': 'rno'}, inplace=True)
            self.sd = self._get_segment_dataframe()
            # period data
            self._period_data = period_data

        # package data
        self._package_data = None

        # connection info (doesn't support diversions)
        self.graph = dict(zip(self.rd.rno, self.rd.outreach))
        self._graph_r = None  # filled by properties
        self.outlets = None
        self._connections = None

        # diversions
        self.diversions = None



        self.dimensions_block = '\nBEGIN Dimensions\n  NREACHES {:d}\nEND Dimensions\n'.format(self.nreaches)

    @property
    def graph_r(self):
        if self._graph_r is None:
            outreaches = np.unique(self.rd.outreach)
            self._graph_r = {o: self.rd.rno[self.rd.outreach == o].tolist()
                             for o in outreaches}
            self.outlets = self._graph_r[0]
            del self._graph_r[0]
        return self._graph_r

    @property
    def options_block(self):
        return self._options_block
    
    @options_block.setter
    def options_block(self, options):
        # stuff to write
        options_block = '\nBEGIN Options\n'
        if options is not None:
            for opt in options:
                options_block += '  {}\n'.format(opt)
        if 'unit_conversion' not in options_block:
            options_block += '  unit_conversion  {}\n'.format(self.unit_conversion)
        if 'auxiliary' not in options_block and self.auxiliary_line_numbers:
            options_block += '  auxiliary line_id\n'
        options_block += 'END Options\n'
        self._options_block = options_block

    @property
    def packagedata(self):
        if self._package_data is None:
            self._package_data = self._get_packagedata()
        return self._package_data

    @property
    def connections(self):
        if self._connections is None:
            connections = {}
            for rno in np.arange(self.nreaches) + 1:
                outreach = self.graph.get(rno)
                #for rno, v in self.graph.items():
                cnk = []
                # downstream reach
                if outreach not in {0, None}:  # outlets aren't explicit in MF6
                    cnk += [-outreach]
                # upstream reaches
                if rno in self.graph_r.keys():
                    cnk += list(self.graph_r[rno])
                connections[rno] = cnk
            self._connections = connections
        return self._connections

    @property
    def period_data(self):
        if self._period_data is None:
            self._period_data = self._get_period_data()
        return self._period_data

    def _segment_data2reach_data(self, var):
        reach_values = []
        sd0 = self.ModflowSfr2.segment_data[0]
        if len(np.unique(sd0[var])) > 1:
            for i in range(len(sd0)):
                seg_value = sd0[i][var]
                nreaches = np.sum(self.rd.iseg == i + 1)
                reach_values += [seg_value]
        else:
            reach_values = np.ones(len(self.rd)) * sd0[var][0]
        return reach_values

    def _get_segment_dataframe(self):
        sd = pd.DataFrame()
        dfs = []
        for k, v in self.ModflowSfr2.segment_data.items():
            df = pd.DataFrame(v)
            df['per'] = k
            dfs.append(df)
        sd = pd.concat(dfs, axis=0)
        keepcols = (sd.sum(axis=0) > 0) | np.isin(sd.columns.values, ['per'])
        return sd.loc[:, keepcols]

    def _get_packagedata(self):
        # [rno, cellid, rlen, rwid, rgrd, rtp, rbth, rhk, man, ncon,
        # ustrf, ndv, aux, boundname]
        print('converting reach and segment data to package data...')
        # rwid = self.modflow_sfr2._interpolate_to_reaches('width1', 'width2')
        rwid = interpolate_to_reaches(self.rd, self.sd,
                                      'width1', 'width2',
                                      reach_data_group_col='iseg',
                                      segment_data_group_col='nseg'
                                      )
        man = self._segment_data2reach_data('roughch')

        packagedata = pd.DataFrame()
        for k, v in self.mf6names.items():
            if k in self.rd.columns:
                packagedata[v] = self.rd[k]

        packagedata['rwid'] = rwid
        packagedata['man'] = man
        packagedata['ncon'] = [len(self.connections.get(k, [])) for k in packagedata.rno]
        packagedata['ustrf'] = 1.
        packagedata['ndv'] = 0

        if self.structured:
            packagedata.drop('cellid', axis=1, inplace=True)
            for dim in ['k', 'i', 'j']:
                packagedata[dim] = self.rd[dim]

        if self.idomain is not None:
            packagedata['idomain'] = self.idomain[packagedata.k, packagedata.i, packagedata.j]
        else:
            packagedata['idomain'] = 1

        if 'auxiliary' in self.options_block:
            aux_variables = self.options_block.split('auxiliary')[1].split('\n')[0].split()
            for var in aux_variables:
                if var in self.rd.columns:
                    packagedata[var] = self.rd[var]
                elif var == 'line_id':  # included by default
                    # if no line_id column; use SFR2 segment
                    packagedata[var] = self.rd['iseg']
                else:
                    raise ValueError(f"Auxiliary variable {var} not in packagedata!")
                    
        cols = [c for c in self.cols if c in packagedata.columns]
        return packagedata[cols].sort_values(by='rno')

    def _get_period_data(self):
        print('converting segment data to period data...')
        return segment_data_to_period_data(self.sd, self.rd)

    def write_file(self, filename=None, outpath='', options=None,
                   external_files_path=None):
        """Write a MODFLOW-6 format SFR package file.

        Parameters
        ----------
        filename : str, optional
            SFR package filename. Default setting is to use the
            ModflowSfr2.file_name attribute for the ModflowSfr2 instance
            entered on init of the Mf6SFR class, by default None.
        outpath : str, optional
            Path to write sfr file (with ModflowSfr2.file_name) to. 
            Usually this is the simulation workspace. 
            Only used if filename is None.
        options : list, optional
            List of strings to write to the MODFLOW-6 SFR options block. For example::

                options=['save_flows',
                         'BUDGET FILEOUT model.sfr.cbc',
                         'STAGE FILEOUT model.sfr.stage.bin']

            An appropriate unit_conversion is written by default.
            See MODFLOW-6 documentation for other options.
            By default None.
        external_files_path : str, optional
            Path for writing an external file for packagedata, relative to the location of the SFR package file.
            If specified, an open/close statement referencing the file is written to the packagedata block.
            By default, None (packagedata table is written to the SFR package file)

        Raises
        ------
        OSError
            If an invalid external_files_path is specified.
        """        
        if filename is not None:
            outfile = filename
            outpath = os.path.split(filename)[0]
        else:
            outfile = os.path.join(outpath, self.ModflowSfr2.file_name[0])

        if options is not None:
            self.options_block = options
            
        if external_files_path is not None:
            full_external_files_path = os.path.join(outpath, external_files_path)
            if not os.path.isdir(full_external_files_path):
                raise OSError("external_files_path doesn't exist:\n{}".format(full_external_files_path))

        header = "# MODFLOW-6 SFR input; created by SFRmaker v. {}".format(sfrmaker.__version__)
        with open(outfile, 'w', newline=""
                  ) as output:
            output.write(header + '\n')
            output.write(self.options_block)
            output.write(self.dimensions_block)

            output.write('\nBEGIN Packagedata\n')
            writepakdata = self.packagedata.copy()
            for c in ['cellid', 'k', 'i', 'j']:
                if c in writepakdata.columns:
                    writepakdata[c] += 1  # convert indices to 1-based
                    writepakdata[c] = writepakdata[c].astype(str)
            # fill in NONEs for reaches in inactive cells
            inactive = writepakdata.idomain != 1
            if not 'cellid' in writepakdata.columns:
                writepakdata.loc[inactive, 'k'] = ''
                writepakdata.loc[inactive, 'i'] = 'NONE'
                writepakdata.loc[inactive, 'j'] = ''
            else:
                writepakdata.loc[inactive, 'cellid'] = 'NONE'

            columns = list(writepakdata.columns)
            columns[0] = '#{}'.format(columns[0])
            writepakdata.columns = columns
            if external_files_path:
                
                # filename for the external package data file
                packagedata_outfile = '{}_packagedata.dat'.format(os.path.splitext(os.path.split(outfile)[-1])[0])
                
                # path to package data file relative to sfr file
                packagedata_rel_path = os.path.join(external_files_path, packagedata_outfile)
                output.write('  open/close {}\n'.format(packagedata_rel_path))
                
                # path to package data file relative to cwd
                packagedata_outfile = os.path.join(full_external_files_path, packagedata_outfile)
                writepakdata.drop('idomain', axis=1).to_csv(packagedata_outfile, sep=' ', index=False)
                print('wrote {}'.format(packagedata_outfile))
                
            else:
                writepakdata.drop('idomain', axis=1).to_csv(output, sep=' ', index=False)
            output.write('END Packagedata\n')

            output.write('\nBEGIN Connectiondata\n')
            #for i in range(1, self.nreaches + 1):
            #    if i in self.connections.keys():
            #        output.write('  {} {}\n'.format(i, ' '.join(map(str, self.connections[i]))))
            for rno, connections in self.connections.items():
                connections = ' '.join(map(str, self.connections[rno]))
                output.write(f'  {rno} {connections}\n')
            output.write('END Connectiondata\n')

            # skip the diversions block for now
            if self.period_data is not None:
                periods = self.period_data.groupby('per')
                for per, group in periods:
                    output.write('\nBEGIN Period {}\n'.format(per + 1))
                    group = group.replace('ACTIVE', np.nan)
                    assert np.all(group.index.get_level_values(0) == per)
                    # drop the period from the index
                    group.reset_index(level=0, drop=True, inplace=True)
                    datacols = {'inflow', 'manning', 'rainfall', 'evaporation', 'runoff', 'stage'}
                    datacols = [c for c in group.columns if c in datacols]
                    group = group.loc[:, datacols]
                    group.stack().to_csv(output, sep=' ', index=True, header=False)
                    output.write('END Period {}\n'.format(per + 1))
        print('wrote {}'.format(outfile))


class mf6sfr(Mf6SFR):
    def __init__(self, *args, **kwargs):
        warnings.warn("The 'mf6sfr' class was renamed to Mf6SFR to better follow pep 8.",
                      DeprecationWarning)
        Mf6SFR.__init__(self, *args, **kwargs)


def segment_data_to_period_data(segment_data, reach_data):
    """Convert modflow-2005 style segment data to modflow-6 period data.
    """
    idx_cols = ['per', 'iseg', 'ireach', 'rno', 'icalc']
    variable_cols = ['status', 'evaporation', 'inflow', 'rainfall', 'runoff', 'stage']
    sd = segment_data.copy()
    prd = segment_data.copy()
    rd = reach_data.copy()
    cols = [k for k, v in Mf6SFR.mf6names.items() if k in prd.columns]
    cols += ['per', 'nseg']
    cols.remove('roughch')

    prd = prd.loc[:, cols]
    # drop all of the zero values from s.p. 0 (default)
    datacols = set(prd.columns.difference({'per', 'nseg'}))
    for c in datacols:  # explicitly convert zeros for each data column
        prd.loc[(prd.per == 0), c] = prd.loc[(prd.per == 0), c].replace(0., np.nan)
    prd.rename(columns=Mf6SFR.mf6names, inplace=True)

    # drop columns that have no values
    prd.dropna(axis=1, how='all', inplace=True)

    # drop rows that have no values
    cols = set(variable_cols).intersection(prd.columns)
    prd.dropna(axis=0, how='all', subset=cols, inplace=True)
    if len(prd) == 0:
        return pd.DataFrame(columns=idx_cols + variable_cols)

    # get rno for reach 1s
    segroups = rd.groupby('iseg')
    reach1IDs = segroups.first()['rno'].to_dict()
    prd['rno'] = [reach1IDs[seg] for seg in prd['nseg']]
    # prd['rainfall'] = np.random.randn(len(prd))

    #  distribute other variables to all reaches
    cols = {'inflow', 'manning', 'rainfall', 'evaporation', 'runoff',
            'depth1', 'depth2'}
    cols = set(prd.columns).intersection(cols)
    icalc = dict(zip(zip(sd.nseg, sd.per), sd.icalc))
    if len(cols) == 0:
        return None
    elif len(cols) > 0:
        print('distributing values to reaches:\n')
        for c in cols.difference({'depth1', 'depth2'}):
            print('{} -> {}...'.format(Mf6SFR.mf5names[c], c))
        if np.min(list(icalc.values())) < 1:
            print('strtop, depth1 & depth2 -> stage for icalc<1...')
        reaches = []
        # iterate through segments in prd (multiple periods)
        for i, r in prd.iterrows():
            seg, per = int(r.nseg), int(r.per)
            df = pd.DataFrame(segroups.get_group(seg)[['rno', 'rchlen', 'ireach']])
            df.rename(columns={'reachID': 'rno'}, inplace=True)
            df['iseg'] = seg
            df['lenfrac'] = df.rchlen / df.rchlen.sum()
            df['per'] = per
            df['icalc'] = icalc[(seg, per)]
            for c in cols:
                if c == 'runoff':
                    df[c] = r[c] * df.lenfrac  # length-weighted mean
                elif c == 'inflow':  # only assign to reach1
                    inflow = np.zeros(len(df)) * np.nan
                    inflow[0] = r.inflow
                    df['inflow'] = inflow
                elif c in {'manning', 'rainfall', 'evaporation'}:
                    df[c] = r[c]  # values already normalized to area

            # distribute depth values from stress period data
            # this will only populate stages where depth1/depth2 were >0
            if icalc[(seg, per)] < 1 and 'depth1' in cols:
                df['status'] = 'SIMPLE'

                # interpolate depth1 and depth2 to reaches
                dist = (np.cumsum(df.rchlen) - 0.5 * df.rchlen).values
                fp = [r.depth1, r.depth2]
                xp = [dist[0], dist[-1]]
                depth = np.interp(dist, xp, fp)
                strtop = segroups.get_group(r.nseg).strtop.values
                df['stage'] = strtop + depth
            else:
                df['status'] = 'ACTIVE'

            # drop out rows that didn't have any data
            df.dropna(subset=cols.difference({'depth1', 'depth2'}), inplace=True)
            reaches.append(df)
        distributed = pd.concat(reaches)
    distributed['per'] = distributed.per.astype(int)
    distributed.index = distributed.rno
    distributed.index.name = 'rno_idx'
    distributed.sort_values(by=['per', 'rno'], inplace=True)

    # for icalc=0 reaches with no depth specified in per 0
    # set stage from strtop; add to distributed period data
    strtop = rd[['rno', 'strtop', 'iseg', 'ireach']].copy()
    strtop.rename(columns={  # 'reachID': 'rno',
        'strtop': 'stage'}, inplace=True)
    strtop['icalc'] = [icalc[(s, 0)] for s in rd.iseg]
    strtop.index = strtop.rno
    # cull to only include per=0 reaches where icalc=0 and depth wasn't specified
    per0reaches = distributed.loc[distributed.per == 0, 'rno']
    strtop = strtop.loc[~strtop.rno.isin(per0reaches) & (strtop.icalc < 1)]
    strtop['per'] = 0
    strtop['status'] = 'SIMPLE'

    if len(strtop) > 0:
        distributed = pd.concat([distributed, strtop], axis=0)
    distributed.sort_values(by=['per', 'rno'], inplace=True)

    # rearrange the columns
    arrangecols = ['per', 'rno', 'status']
    if 'stage' in distributed.columns: arrangecols += ['stage']
    for col in cols.difference({'depth1', 'depth2'}):
        arrangecols.append(col)
    arrangecols += ['iseg', 'ireach', 'icalc']
    distributed = distributed[arrangecols]
    distributed['iseg'] = distributed.iseg.astype(int)
    return distributed


def cellids_to_kij(cellids, drop_inactive=True):
    """Unpack tuples of MODFLOW-6 cellids (k, i, j) to
    lists of k, i, j values; ignoring instances
    where cellid is None (unconnected cells).

    Parameters
    ----------
    cellids : sequence of (k, i, j) tuples
    drop_inactive : bool
        If True, drop cellids == 'none'. If False,
        distribute these to k, i, j.

    Returns
    -------
    k, i, j : 1D numpy arrays of integers
    """
    active = np.array(cellids) != 'none'
    if drop_inactive:
        k, i, j = map(np.array, zip(*cellids[active]))
    else:
        k = np.array([cid[0] if cid != 'none' else None for cid in cellids])
        i = np.array([cid[1] if cid != 'none' else None for cid in cellids])
        j = np.array([cid[2] if cid != 'none' else None for cid in cellids])
    return k, i, j