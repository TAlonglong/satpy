#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2018 Trygve Aspenes

# Author(s):

#   Trygve Aspenes <trygveas@met.no>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""SAFE SAR L2 OCN format."""

import logging
import os
import xml.etree.ElementTree as ET
from pyresample import geometry

from satpy.readers.file_handlers import BaseFileHandler
from satpy import CHUNK_SIZE

import numpy as np
import xarray as xr

logger = logging.getLogger(__name__)


class SAFENC(BaseFileHandler):
    """Measurement file reader."""

    def __init__(self, filename, filename_info, filetype_info):
        super(SAFENC, self).__init__(filename, filename_info,
                                     filetype_info)

        self._start_time = filename_info['start_time']
        self._end_time = filename_info['end_time']

        self._polarization = filename_info['polarization']

        self.lats = None
        self.lons = None
        self._shape = None
        self.area = None

        self.nc = xr.open_dataset(filename,
                                  decode_cf=True,
                                  mask_and_scale=False,
                                  chunks={'owiAzSize': CHUNK_SIZE,
                                          'owiRaSize': CHUNK_SIZE})
        self.nc = self.nc.rename({'owiAzSize': 'y'})
        self.nc = self.nc.rename({'owiRaSize': 'x'})
        print self.nc
        print "WindDirection ", self.nc['owiWindDirection']
        print "WindSpeed ", self.nc['owiWindSpeed']
        self.filename = filename

    def get_dataset(self, key, info):
        """Load a dataset."""
        logger.debug("REader %s %s", key, info)

        logger.debug('Reading keyname %s.', key.name)
        if key.name in ['owiLat', 'owiLon']:
            logger.debug('Constructing coordinate arrays ll.')

            if self.lons is None or self.lats is None:
                self.lons = self.nc['owiLon']
                self.lats = self.nc['owiLat']
            if key.name == 'owiLat':
                res = self.lats
            else:
                res = self.lons
            res.attrs = info
        else:
            logger.debug("Read data")
            res = self.nc[key.name]
            res = xr.DataArray(res, dims=['y', 'x'])
            res.attrs.update(info)
            if '_FillValue' in res.attrs:
                res = res.where(res != res.attrs['_FillValue'])
                res.attrs['_FillValue'] = np.nan

            print "DATA:", self.nc[key.name]
            print "END"

        if 'missionName' in self.nc.attrs:
            res.attrs.update({'platform_name': self.nc.attrs['missionName']})

        print "res.shape: ",res.shape
        if not self._shape:
            self._shape = res.shape

        return res

    @property
    def start_time(self):
        return self._start_time

    @property
    def end_time(self):
        return self._end_time
