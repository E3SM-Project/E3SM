#!/usr/bin/env python
"""

Utility functions to read and manipulate JIGSAW meshes.

Phillip J. Wolfram
04/06/2017
"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np


def readmsh(fname):
    """
    Reads JIGSAW msh structure and produces a dictionary with values.

    Phillip J. Wolfram
    09/22/2017
    """

    dataset = {}
    datavals = {}
    datavals['HEADER'] = ';'
    datavals['ARRAY'] = None
    with open(fname) as f:
        line = f.readline()
        while line:
            if line[0] == '#':
                datavals['HEADER'] += line[1:] + ';'
                line = f.readline()
                continue
            if '=' in line:
                datavals, dataset = _store_datavals(datavals, dataset)
                if 'COORD' in line:
                    name = 'COORD' + line.split('=')[1][0]
                    datavals[name] = line.split(';')[-1]
                else:
                    vals = line.split('=')
                    value = vals[1] if ';' in vals[1] else int(vals[1])
                    datavals[vals[0]] = value
                line = f.readline()
                continue

            # just numbers
            arrayvals = np.asarray(line.split(';'), dtype='f8')
            if datavals['ARRAY'] is None:
                datavals['ARRAY'] = [arrayvals]
            else:
                datavals['ARRAY'].append(arrayvals)
            line = f.readline()
            continue
        datavals, dataset = _store_datavals(datavals, dataset)

    return dataset


def _store_datavals(datavals, dataset):  # {{{

    if datavals['ARRAY'] is not None:
        # remove empty data
        if np.all(datavals['ARRAY'] == np.array(None, dtype='object')):
            datavals.pop('ARRAY')
        for key in [aval for aval in datavals.keys()
                    if aval in ['HEADER', 'MSHID', 'NDIMS']]:
            if key in dataset:
                dataset[key] += datavals[key]
            else:
                dataset[key] = datavals[key]
            datavals.pop(key)
        entryname = [aval for aval in datavals.keys() if aval not in [
            'ARRAY']]

        if 'TRI' in entryname[0]:
            dtype = 'i'
        else:
            dtype = 'f8'
        datavals['ARRAY'] = np.asarray(datavals['ARRAY'], dtype=dtype)

        # decided to throw away "index" from msh because it isn't truly a
        # real number
        dataset[entryname[0]] = datavals['ARRAY']
        datavals = {}
        datavals['ARRAY'] = None

    return datavals, dataset  # }}}
