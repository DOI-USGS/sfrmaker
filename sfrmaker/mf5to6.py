import os
import numpy as np
import pandas as pd
from copy import copy
from .utils import interpolate_to_reaches
import sfrmaker


class mf6sfr:

    # convert from ModflowSfr to mf6
    mf6names = {'reachID': 'rno',
                'node': 'cellid',
                'rchlen': 'rlen',
                'slope': 'rgrd',
                'strtop': 'rtp',
                'strthick': 'rbth',
                'strhc1': 'rhk',
                'roughch': 'manning',
                'flow': 'inflow',
                'pptsw': 'rainfall',
                'etsw': 'evaporation',
                'runoff': 'runoff',
                'depth1': 'depth1', # need these for distributing stage
                'depth2': 'depth2'}
    mf5names = {v :k for k, v in mf6names.items()}

    cols = ['rno', 'cellid', 'k', 'i', 'j', 'rlen', 'rwid', 'rgrd',
            'rtp', 'rbth', 'rhk', 'man', 'ncon', 'ustrf', 'ndv', 'idomain']

    def __init__(self, ModflowSfr2, options=['print_input',
                                             'save_flows']):

        # check for other packages
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

        # copy ModflowSfr2 object and enforce sorting
        self.ModflowSfr2 = copy(ModflowSfr2)
        self.ModflowSfr2.segment_data[0].sort(order='nseg')
        self.ModflowSfr2.reach_data.sort(order=['iseg', 'ireach'])

        self.structured = self.ModflowSfr2.parent.structured
        self.const = ModflowSfr2.const
        self.nreaches = len(ModflowSfr2.reach_data)
        self.nper = ModflowSfr2.nper

        # dataframes of reach and segment data from ModflowSfr2
        self.rd = pd.DataFrame(ModflowSfr2.reach_data)
        self.sd = self._get_segment_dataframe()

        # package data
        self._package_data = None

        # connection info (doesn't support diversions)
        self.graph = dict(zip(self.rd.reachID, self.rd.outreach))
        sd0 = pd.DataFrame(ModflowSfr2.segment_data[0])
        self._graph_r = None # filled by properties
        self.outlets = None
        self._connections = None

        # diversions
        self.diversions = None

        # period data
        self._period_data = None

        # stuff to write
        self.options_block = '\nBEGIN Options\n'
        for opt in options:
            self.options_block += '  {}\n'.format(opt)
        self.options_block += '  unit_conversion  {}\n'.format(ModflowSfr2.const)
        self.options_block += 'END Options\n'

        self.dimensions_block = '\nBEGIN Dimensions\n  NREACHES {:d}\nEND Dimensions\n'.format(self.nreaches)

    @property
    def graph_r(self):
        if self._graph_r is None:
            outreaches = np.unique(self.rd.outreach)
            self._graph_r = {o: self.rd.reachID[self.rd.outreach == o].tolist()
                             for o in outreaches}
            self.outlets = self._graph_r[0]
            del self._graph_r[0]
        return self._graph_r

    @property
    def packagedata(self):
        if self._package_data is None:
            self._package_data = self._get_packagedata()
        return self._package_data

    @property
    def connections(self):
        if self._connections is None:
            connections = {}
            for k, v in self.graph.items():
                cnk = []
                # downstream reach
                if v != 0: # outlets aren't explicit in MF6
                    cnk += [-v]
                # upstream reaches
                if k in self.graph_r.keys():
                    cnk += list(self.graph_r[k])
                if len(cnk) > 0:
                    connections[k] = cnk
            self._connections = connections
        return self._connections

    @property
    def perioddata(self):
        if self._period_data is None:
            self._period_data = self._get_perioddata()
        return self._period_data

    def _segment_data2reach_data(self, var):
        reach_values = []
        sd0 = self.ModflowSfr2.segment_data[0]
        if len(np.unique(sd0[var])) > 1:
            for i in range(len(sd0)):
                seg_value = sd0[i][var]
                nreaches = np.sum(self.rd.iseg == i+ 1)
                reach_values += [seg_value]
        else:
            reach_values = np.ones(len(self.rd)) * sd0[var][0]
        return reach_values

    def _get_segment_dataframe(self):
        sd = pd.DataFrame()
        for k, v in self.ModflowSfr2.segment_data.items():
            df = pd.DataFrame(v)
            df['per'] = k
            sd = sd.append(df)
        keepcols = (sd.sum(axis=0) > 0) | np.in1d(sd.columns.values, ['per'])
        return sd.loc[:, keepcols]

    def _get_packagedata(self):
        # [rno, cellid, rlen, rwid, rgrd, rtp, rbth, rhk, man, ncon,
        # ustrf, ndv, aux, boundname]
        print('converting reach and segment data to package data...')
        #rwid = self.ModflowSfr2._interpolate_to_reaches('width1', 'width2')
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

        cols = [c for c in self.cols if c in packagedata.columns]
        return packagedata[cols].sort_values(by='rno')

    def _get_perioddata(self):
        print('converting segment data to period data...')

        sd = self.sd.copy()
        cols = [k for k, v in self.mf6names.items() if k in sd.columns]
        cols += ['per', 'nseg']
        cols.remove('roughch')

        sd = sd.loc[:, cols]
        # drop all of the zero values from s.p. 0 (default)
        datacols = set(sd.columns.difference({'per', 'nseg'}))
        for c in datacols:  # explicitly convert zeros for each data column
            sd.loc[(sd.per == 0), c] = sd.loc[(sd.per == 0), c].replace(0., np.nan)
        sd.head()
        sd.rename(columns=self.mf6names, inplace=True)

        # get rno for reach 1s
        segroups = self.rd.groupby('iseg')
        reach1IDs = segroups.first()['reachID'].to_dict()
        sd['rno'] = [reach1IDs[i + 1] for i in sd.index]
        # sd['rainfall'] = np.random.randn(len(sd))

        #  distribute other variables to all reaches
        cols = {'inflow', 'manning', 'rainfall', 'evaporation', 'runoff',
                'depth1', 'depth2'}
        cols = set(sd.columns).intersection(cols)
        icalc = dict(zip(zip(self.sd.nseg, self.sd.per), self.sd.icalc))
        if len(cols) == 0:
            return None
        elif len(cols) > 0:
            print('distributing values to reaches:\n')
            for c in cols.difference({'depth1', 'depth2'}):
                print('{} -> {}...'.format(self.mf5names[c], c))
            if np.min(list(icalc.values())) < 1:
                print('strtop, depth1 & depth2 -> stage for icalc<1...')
            reaches = []
            # iterate through segments in sd (multiple periods)
            for i, r in sd.iterrows():
                seg, per = int(r.nseg), int(r.per)
                df = pd.DataFrame(segroups.get_group(seg)[['reachID', 'rchlen', 'ireach']])
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
        strtop = self.rd[['reachID', 'strtop', 'iseg', 'ireach']].copy()
        strtop.rename(columns={'reachID': 'rno', 'strtop': 'stage'}, inplace=True)
        strtop['icalc'] = [icalc[(s, 0)] for s in self.rd.iseg]
        strtop.index = strtop.rno
        # cull to only include per=0 reaches where icalc=0 and depth wasn't specified
        per0reaches = distributed.loc[distributed.per == 0, 'rno']
        strtop = strtop.loc[~strtop.rno.isin(per0reaches) & (strtop.icalc < 1)]
        strtop['per'] = 0
        strtop['status'] = 'SIMPLE'

        distributed = distributed.append(strtop)
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

    def write_file(self, filename=None, outpath=''):
        if filename is not None:
            outfile = filename
        else:
            outfile = os.path.join(outpath, self.ModflowSfr2.file_name[0] + '6')
        header = "# MODFLOW-6 SFR input; created by SFRmaker v. {}".format(sfrmaker.__version__)
        with open(outfile, 'w') as output:
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
            inactive = writepakdata.idomain == 0
            if not 'cellid' in writepakdata.columns:
                writepakdata.loc[inactive, 'k'] = ''
                writepakdata.loc[inactive, 'i'] = 'NONE'
                writepakdata.loc[inactive, 'j'] = ''
            else:
                writepakdata.loc[inactive, 'cellid'] = 'NONE'

            columns = list(writepakdata.columns)
            columns[0] = '#{}'.format(columns[0])
            writepakdata.columns = columns
            writepakdata.drop('idomain', axis=1).to_csv(output, sep=' ', index=False)
            output.write('END Packagedata\n')

            output.write('\nBEGIN Connectiondata\n')
            for i in range(1, self.nreaches + 1):
                if i in self.connections.keys():
                    output.write('  {} {}\n'.format(i, ' '.join(map(str, self.connections[i]))))
            output.write('END Connectiondata\n')

            # skip the diversions block for now
            if self.perioddata is not None:
                periods = self.perioddata.groupby('per')
                for per in range(self.nper):
                    output.write('\nBEGIN Period {}\n'.format(per + 1))
                    grp = periods.get_group(per).replace('ACTIVE', np.nan)
                    assert np.array_equal(grp.index.values, grp.rno.values)
                    datacols = {'inflow', 'manning', 'rainfall', 'evaporation', 'runoff', 'stage'}
                    datacols = datacols.intersection(grp.columns)
                    grp = grp.loc[:, datacols]
                    grp.stack().to_csv(output, sep=' ', index=True, header=False)
                    output.write('END Period {}\n'.format(per + 1))
        print('wrote {}'.format(outfile))
