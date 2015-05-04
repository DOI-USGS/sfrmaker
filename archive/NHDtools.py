__author__ = 'aleaf'
'''
Classes for working with information in NHDPlus v2

* ``Flowlines`` - class for working with NHD "Flowlines" hydrography
'''

import GISio
import numpy as np
import pandas as pd

class Flowlines(object):

    def __init__(self, filename):
        '''
        ``filename``: (string or list) shapefile(s) to load
        '''
        self.filename = filename
        self.df = GISio.shp2df(self.filename, index='COMID', geometry=True)

    def add_fcode(self, fcodefile):

        self.fcode = GISio.shp2df(fcodefile)

    def add_PlusFlowVAA(self, pfvaafile, **kwargs):

        self.pfvaa = GISio.shp2df(pfvaafile, index='COMID')
        self.df = self.df.join(self.pfvaa, **kwargs)

    def portion_perennial(self, bounds=None):
        '''
        returns new dataframe of COMIDs classified as perennial (True) or ephemeral (False)
        ``bounds (optional)``: (polygon shapefile) encompassing streams to evaluate
        '''

        streams = self.df[(self.df.FCODE == 46006) | (self.df.FCODE == 46003)].copy()
        streams['perennial'] = [True if r.FCODE == 46006 else False for i, r in streams.iterrows()]

        if bounds:

            bound = GISio.shp2df(bounds, geometry=True).iloc[0].geometry
            intersected = [g.intersects(bound) for g in streams.geometry]
            streams = streams.ix[intersected]

        return streams


