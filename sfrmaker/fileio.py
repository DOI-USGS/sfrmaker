import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
try:
    import flopy
except:
    flopy = False

from sfrmaker.utils import get_input_arguments


def load_mf2005_package(f, model=None):

    if model.verbose:
        sys.stdout.write('loading sfr2 package file...\n')

    tabfiles = False
    tabfiles_dict = {}
    transroute = False
    reachinput = False
    structured = model.structured
    if nper is None:
        nper = model.nper
        nper = 1 if nper == 0 else nper  # otherwise iterations from 0, nper won't run

    openfile = not hasattr(f, 'read')
    if openfile:
        filename = f
        f = open(filename, 'r')

    # Item 0 -- header
    while True:
        line = f.readline()
        if line[0] != '#':
            break

    options = None
    if model.version == "mfnwt" and "options" in line.lower():
        options = OptionBlock.load_options(f, ModflowSfr2)

    else:
        query = ("reachinput", "transroute", "tabfiles",
                 "lossfactor", "strhc1kh", "strhc1kv")
        for i in query:
            if i in line.lower():
                options = OptionBlock(line.lower().strip(),
                                      ModflowSfr2, block=False)
                break

    if options is not None:
        line = f.readline()
        # check for 1b in modflow-2005
        if "tabfile" in line.lower():
            t = line.strip().split()
            options.tabfiles = True
            options.numtab = int(t[1])
            options.maxval = int(t[2])
            line = f.readline()

        # set varibles to be passed to class args
        transroute = options.transroute
        reachinput = options.reachinput
        tabfiles = isinstance(options.tabfiles, np.ndarray)
        numtab = options.numtab if tabfiles else 0

    # item 1c
    nstrm, nss, nsfrpar, nparseg, const, dleak, ipakcb, istcb2, \
    isfropt, nstrail, isuzn, nsfrsets, \
    irtflg, numtim, weight, flwtol, option = _parse_1c(line,
                                                       reachinput=reachinput,
                                                       transroute=transroute)

    # item 2
    # set column names, dtypes
    names = _get_item2_names(nstrm, reachinput, isfropt, structured)
    dtypes = [d for d in ModflowSfr2.get_default_reach_dtype().descr
              if d[0] in names]

    lines = []
    for i in range(abs(nstrm)):
        line = f.readline()
        line = line_parse(line)
        ireach = tuple(map(float, line[:len(dtypes)]))
        lines.append(ireach)

    tmp = np.array(lines, dtype=dtypes)
    # initialize full reach_data array with all possible columns
    reach_data = ModflowSfr2.get_empty_reach_data(len(lines))
    for n in names:
        reach_data[n] = tmp[
            n]  # not sure if there's a way to assign multiple columns

    # zero-based convention
    inds = ['k', 'i', 'j'] if structured else ['node']
    _markitzero(reach_data, inds)

    # items 3 and 4 are skipped (parameters not supported)
    # item 5
    segment_data = {}
    channel_geometry_data = {}
    channel_flow_data = {}
    dataset_5 = {}
    aux_variables = {}  # not sure where the auxiliary variables are supposed to go
    for i in range(0, nper):
        # Dataset 5
        dataset_5[i] = _get_dataset(f.readline(), [-1, 0, 0, 0])
        itmp = dataset_5[i][0]
        if itmp > 0:
            # Item 6
            current = ModflowSfr2.get_empty_segment_data(nsegments=itmp,
                                                         aux_names=option)
            # container to hold any auxiliary variables
            current_aux = {}
            # these could also be implemented as structured arrays with a column for segment number
            current_6d = {}
            current_6e = {}
            # print(i,icalc,nstrm,isfropt,reachinput)
            for j in range(itmp):
                dataset_6a = _parse_6a(f.readline(), option)
                current_aux[j] = dataset_6a[-1]
                dataset_6a = dataset_6a[:-1]  # drop xyz
                icalc = dataset_6a[1]
                # link dataset 6d, 6e by nseg of dataset_6a
                temp_nseg = dataset_6a[0]
                # datasets 6b and 6c aren't read under the conditions below
                # see table under description of dataset 6c,
                # in the MODFLOW Online Guide for a description
                # of this logic
                # https://water.usgs.gov/ogw/modflow-nwt/MODFLOW-NWT-Guide/sfr.htm
                dataset_6b, dataset_6c = (0,) * 9, (0,) * 9
                if not (isfropt in [2, 3] and icalc == 1 and i > 1) and \
                        not (isfropt in [1, 2, 3] and icalc >= 2):
                    dataset_6b = _parse_6bc(f.readline(), icalc, nstrm,
                                            isfropt,
                                            reachinput, per=i)
                    dataset_6c = _parse_6bc(f.readline(), icalc, nstrm,
                                            isfropt,
                                            reachinput, per=i)
                current[j] = dataset_6a + dataset_6b + dataset_6c

                if icalc == 2:
                    # ATL: not sure exactly how isfropt logic functions for this
                    # dataset 6d description suggests that this line isn't read for isfropt > 1
                    # but description of icalc suggest that icalc=2 (8-point channel) can be used with any isfropt
                    if i == 0 or nstrm > 0 and not reachinput or isfropt <= 1:
                        dataset_6d = []
                        for _ in range(2):
                            dataset_6d.append(
                                _get_dataset(f.readline(), [0.0] * 8))
                            # dataset_6d.append(list(map(float, f.readline().strip().split())))
                        current_6d[temp_nseg] = dataset_6d
                if icalc == 4:
                    nstrpts = dataset_6a[5]
                    dataset_6e = []
                    for _ in range(3):
                        dataset_6e.append(
                            _get_dataset(f.readline(), [0.0] * nstrpts))
                    current_6e[temp_nseg] = dataset_6e

            segment_data[i] = current
            aux_variables[j + 1] = current_aux
            if len(current_6d) > 0:
                channel_geometry_data[i] = current_6d
            if len(current_6e) > 0:
                channel_flow_data[i] = current_6e

        if tabfiles and i == 0:
            for j in range(numtab):
                segnum, numval, iunit = map(int,
                                            f.readline().strip().split())
                tabfiles_dict[segnum] = {'numval': numval, 'inuit': iunit}

        else:
            continue

    if openfile:
        f.close()

    # determine specified unit number
    unitnumber = None
    filenames = [None, None, None]
    if ext_unit_dict is not None:
        for key, value in ext_unit_dict.items():
            if value.filetype == ModflowSfr2.ftype():
                unitnumber = key
                filenames[0] = os.path.basename(value.filename)

            if ipakcb > 0:
                if key == ipakcb:
                    filenames[1] = os.path.basename(value.filename)
                    model.add_pop_key_list(key)

            if abs(istcb2) > 0:
                if key == abs(istcb2):
                    filenames[2] = os.path.basename(value.filename)
                    model.add_pop_key_list(key)

    return ModflowSfr2(model, nstrm=nstrm, nss=nss, nsfrpar=nsfrpar,
                       nparseg=nparseg, const=const, dleak=dleak,
                       ipakcb=ipakcb, istcb2=istcb2,
                       isfropt=isfropt, nstrail=nstrail, isuzn=isuzn,
                       nsfrsets=nsfrsets, irtflg=irtflg,
                       numtim=numtim, weight=weight, flwtol=flwtol,
                       reach_data=reach_data,
                       segment_data=segment_data,
                       dataset_5=dataset_5,
                       channel_geometry_data=channel_geometry_data,
                       channel_flow_data=channel_flow_data,
                       reachinput=reachinput, transroute=transroute,
                       tabfiles=tabfiles, tabfiles_dict=tabfiles_dict,
                       unit_number=unitnumber, filenames=filenames,
                       options=options)


def read_tables(data, **kwargs):
    # allow input via a list of tables or single table
    input_data = data
    if not isinstance(input_data, list):
        input_data = [input_data]
    dfs = []
    for item in input_data:
        if isinstance(item, str) or isinstance(item, Path):
            dfs.append(pd.read_csv(item, **kwargs))
        elif isinstance(item, pd.DataFrame):
            item = item.copy()
            if 'dtype' in kwargs:
                for col, dtype in kwargs['dtype'].items():
                    item[col] = item[col].astype(dtype)
            dfs.append(item.copy())
        else:
            raise Exception('Unrecognized input type for data:\n{}'.format(item))
    data = pd.concat(dfs).reset_index(drop=True)
    return data


def load_json(jsonfile):
    """Convenience function to load a json file; replacing
    some escaped characters."""
    with open(jsonfile) as f:
        return json.load(f)


def load_modelgrid(filename):
    """Create a MFsetupGrid instance from model config json file."""
    cfg = load_json(filename)
    rename = {'xll': 'xoff',
              'yll': 'yoff',
              'epsg': 'crs'
              }
    for k, v in rename.items():
        if k in cfg:
            cfg[v] = cfg.pop(k)
    if np.isscalar(cfg['delr']):
        cfg['delr'] = np.ones(cfg['ncol'])* cfg['delr']
    if np.isscalar(cfg['delc']):
        cfg['delc'] = np.ones(cfg['nrow']) * cfg['delc']
    kwargs = get_input_arguments(cfg, flopy.discretization.StructuredGrid)
    return flopy.discretization.StructuredGrid(**kwargs)


def read_mf6_block(filename, blockname):
    blockname = blockname.lower()
    data = {}
    read = False
    per = None
    with open(filename) as src:
        for line in src:
            line = line.lower().strip().replace('\\', '/')
            if 'begin' in line and blockname in line:
                if blockname == 'period':
                    per = int(line.split()[-1])
                    data[per] = []
                elif blockname == 'continuous':
                    fname = line.split()[-1]
                    data[fname] = []
                elif blockname == 'packagedata':
                    data['packagedata'] = []
                else:
                    blockname = line.split()[-1]
                    data[blockname] = []
                read = blockname
                continue
            if 'end' in line and blockname in line:
                per = None
                read = False
                #break
            if read == 'options':
                line = line.split()
                data[line[0]] = line[1:]
            elif read == 'packages':
                pckg, fname, ext = line.split()
                data[pckg] = fname
            elif read == 'period':
                data[per].append(' '.join(line.split()))
            elif read == 'continuous':
                data[fname].append(' '.join(line.split()))
            elif read == 'packagedata':
                data['packagedata'].append(' '.join(line.split()))
            elif read == blockname:
                data[blockname].append(' '.join(line.split()))
    return data