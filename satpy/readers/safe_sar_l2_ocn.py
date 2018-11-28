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


def dictify(r, root=True):
    """Convert an ElementTree into a dict."""
    if root:
        return {r.tag: dictify(r, False)}
    d = {}
    if r.text and r.text.strip():
        try:
            return int(r.text)
        except ValueError:
            try:
                return float(r.text)
            except ValueError:
                return r.text
    for x in r.findall("./*"):
        print x, x.tag
        if x.tag in d and not isinstance(d[x.tag], list):
            d[x.tag] = [d[x.tag]]
            d[x.tag].append(dictify(x, False))
        else:
            d[x.tag] = dictify(x, False)
    return d


class SAFEXML(BaseFileHandler):
    """XML file reader for the SAFE format."""

    def __init__(self, filename, filename_info, filetype_info):
        print 'SAFEXML init'
        super(SAFEXML, self).__init__(filename, filename_info, filetype_info)


        self.root = ET.parse(self.filename)
        #rt = self.root.getroot()
        #for coordinates in rt.findall('gml:coordinates'):
        #    print coordinates
        #print 'After coordinates'
        #print dictify(self.root.getroot())
        #self.hdr = {}
        #if header_file is not None:
        #    self.hdr = header_file.get_metadata()
        #    print 'self.hdr', self.hdr
        print "SAFEXML END INIT"

    def get_metadata(self):
        """Convert the xml metadata to dict."""
        print "get_metadata"
        return dictify(self.root.getroot())

    def get_dataset(self, key, info):
        print "get_dataset XML"
        return

        #    @property
#    def start_time(self):
#        return self._start_time

#    @property
#    def end_time(self):
#        return self._end_time


class SAFENC(BaseFileHandler):
    """Measurement file reader."""

    def __init__(self, filename, filename_info, filetype_info):
        print "INIT SAFENC"
        super(SAFENC, self).__init__(filename, filename_info,
                                      filetype_info)

        #self.manifest = manifest_fh
        #print "manifest_fh ", manifest_fh
        #self.manifest.get_metadata()

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
        print self.nc['owiWindDirection']
        self.filename = filename
        print "END INIT"
        #self.get_gdal_filehandle()

    def get_dataset(self, key, info):
        """Load a dataset."""
        logger.debug("REader %s %s",key, info)
        #if self._polarization != key.polarization:
        #    return

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

        #print self.nc.attrs
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

#    def get_area_def(self, ds_id):
        #print 'ds_id: ',ds_id
        #print 'ds_id.name: ',ds_id.name

        #data = self[ds_id.name]

#        proj_dict = {#'init': 'epsg:4326'}
#            'proj': 'eqc',
#            'datum': 'WGS84',
#            'ellps': 'WGS84',
#            'no_defs': True
#        }

#        from pyproj import Proj, transform
#        inProj = Proj(init='epsg:4326')
#        outProj = Proj('+proj=eqc +datum=WGS84 +ellps=WGS84')
#        lon_ll, lat_ll = 2.300587, 66.669945 # 8.146225, 66.177345 # 2.300587, 66.669945
#        x_ll, y_ll = transform(inProj, outProj, lon_ll, lat_ll)
#        print x_ll, y_ll
#        lon_ur, lat_ur = 9.148409, 67.643173 # 2.944688, 68.154190 # 9.148409, 67.643173
#        x_ur, y_ur = transform(inProj, outProj, lon_ur, lat_ur)
#        print x_ur, y_ur

#        res_x = (x_ur - x_ll)/166.
#        res_y = (y_ur - y_ll)/266.
#        print "RES: ", res_x, res_y
        # A string of 4 lon, lat coordinate pairs which describe
        # the corners of the image. The string is of the form:
        # lon,lat lon,lat lon,lat lon,lat
        # The coordinates must appear in the following order:
        # last line first pixel, last line last pixel, first line last
        # pixel, first line first pixel
        # <gml:coordinates>66.177345,8.146225 66.669945,2.300587 68.154190,2.944688 67.643173,9.148409</gml:coordinates>

        # area_extent: (x_ll, y_ll, x_ur, y_ur)
        # where

        # x_ll: projection x coordinate of lower left corner of lower left pixel
        # y_ll: projection y coordinate of lower left corner of lower left pixel
        # x_ur: projection x coordinate of upper right corner of upper right pixel
        # y_ur: projection y coordinate of upper right corner of upper right pixel

#        area_extent = (x_ll, y_ll, x_ur, y_ur)

#        area = geometry.AreaDefinition(
#            'sar_ocn_area',
 #           'name_of_proj',
#            'id_of_proj',
#            proj_dict,
#            int(266),
#            int(166),
#            area_extent
#        )

#        return area
