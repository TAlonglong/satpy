#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2017.

# Author(s):

#   Trygve Aspenes <trygveas@met.no>

# This file is part of satpy.

# satpy is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# satpy is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# satpy.  If not, see <http://www.gnu.org/licenses/>.
"""MITIFF writer objects for creating MITIFF files from `Dataset` objects.

"""

import logging

import numpy as np

#import pymitiff.mitiff as mt

from satpy.utils import ensure_dir
from satpy.writers import ImageWriter

import os
from TiffImagePlugin import IMAGEDESCRIPTION
#from TiffImagePlugin import COMPRESSION
#from TiffImagePlugin import ZIPQUALITY

LOG = logging.getLogger(__name__)

KELVIN_TO_CELSIUS = -273.15

class MITIFFWriter(ImageWriter):

    def __init__(self, floating_point=False, tags=None, **kwargs):
        ImageWriter.__init__(self,
                             default_config_filename="writers/mitiff.yaml",
                             **kwargs)

        self.tags = self.info.get("tags",
                                            None) if tags is None else tags
        if self.tags is None:
            self.tags = {}
        elif not isinstance(self.tags, dict):
            # if it's coming from a config file
            self.tags = dict(tuple(x.split("=")) for x in self.tags.split(","))

        self.config = {}

    def _load_mitiff_config(self, mitiff_config_file):
        #Load a mitiff yaml config file
        import yaml
        with open(mitiff_config_file, 'r') as stream:
            try:
                self.config = yaml.load(stream)
                import pprint
                print type(self.config)
                pp = pprint.PrettyPrinter(indent=4)
                pp.pprint(self.config)
            except yaml.YAMLError as exc:
                print(exc)
                raise
                
        #return config

    def save_datasets(self, datasets, **kwargs):
        """Save all datasets to one or more files.
        """
        LOG.debug("Starting in save_datasetsssssssssssssssss ... ")
        LOG.debug("kwargs: {}".format(kwargs))
        try:
            #self._load_mitiff_config(os.path.join(self.ppp_config_dir,"mitiff-config.yaml"))

            if type(kwargs["sensor"]) not in (tuple, list, set):
                kwargs['sensor'] = kwargs['sensor'].replace("/","-")
                #if kwargs['sensor'] not in self.config:
                #    LOG.error("Sensor {} not defined in config. Go fix your config!".format(kwargs['sensor']))
                #    return False
            else:
                for i,sensor in enumerate(kwargs["sensor"]):
                    kwargs["sensor"][i] = sensor.replace("/","-")
                    #if sensor not in self.config:
                    #    LOG.error("Sensor {} not defined in config. Go fix your config!".format(sensor))
                    #    return False

            image_description = self._make_image_description(datasets, **kwargs)
            print "File pattern {}".format(self.file_pattern)
            kwargs['name']  ="shallalal"
            kwargs['start_time'] = datasets[0].info['start_time']
            print kwargs
            print self.get_filename(**kwargs)
            gen_filename = self.get_filename(**kwargs)#self._generate_filename(datasets, **kwargs)
            self._save_datasets_as_mitiff(datasets, image_description, gen_filename, **kwargs)
        except:
            raise

        return gen_filename

#    def save_dataset(self, dataset, filename=None, fill_value=None, overlay=None, **kwargs):
#        """Saves the *dataset* to a given *filename*.
#        """
#        LOG.debug("Starting in save_dataset ... ")
#        with tifffile.TiffWriter(output_fn, **tifargs) as tif:
#            tif.save(image_data, **args)


#    def save_image(self, img, filename=None, floating_point=False, **kwargs):
#        """Save the image to the given *filename* in mitiff format.
#        `floating_point` allows the saving of
#        'L' mode images in floating point format if set to True.
#        """

#        LOG.debug("About to start save_image in mitiff")
#        filename = filename or self.get_filename(**img.info)
#        mt.save(img, filename, **kwargs)
#        LOG.debug("save_image complete in mitiff")

#    def _make_calibration_string(self, datasets)
#        """Make the calibration string to insert in as a tifftag description.
#        """

    def _make_image_description(self, datasets, **kwargs):
        #generate image desdcription for mitiff.
        """
        Satellite: NOAA 18
        Date and Time: 06:58 31/05-2016
        SatDir: 0
        Channels:   6 In this file: 1-VIS0.63 2-VIS0.86 3(3B)-IR3.7 4-IR10.8 5-IR11.5 6(3A)-VIS1.6
        Xsize:  4720
        Ysize:  5544
        Map projection: Stereographic
        Proj string: +proj=stere +lon_0=0 +lat_0=90 +lat_ts=60 +ellps=WGS84 +towgs84=0,0,0 +units=km +x_0=2526000.000000 +y_0=5806000.000000
        TrueLat: 60 N
        GridRot: 0
        Xunit:1000 m Yunit: 1000 m
        NPX: 0.000000 NPY: 0.000000
        Ax: 1.000000 Ay: 1.000000 Bx: -2526.000000 By: -262.000000

        Satellite: <satellite name>
        Date and Time: <HH:MM dd/mm-yyyy>
        SatDir: 0
        Channels:   <number of chanels> In this file: <channels names in order>
        Xsize:  <number of pixels x>
        Ysize:  <number of pixels y>
        Map projection: Stereographic
        Proj string: <proj4 string with +x_0 and +y_0 which is the positive distance from proj origo to the lower left corner of the image data>
        TrueLat: 60 N
        GridRot: 0
        Xunit:1000 m Yunit: 1000 m
        NPX: 0.000000 NPY: 0.000000
        Ax: <pixels size x in km> Ay: <pixel size y in km> Bx: <left corner of upper right pixel in km> By: <upper corner of upper right pixel in km>
     
        
        if palette image write special palette
        if normal channel write table calibration:
        Table_calibration: <channel name>, <calibration type>, [<unit>], <no of bits of data>, [<calibration values space separated>]\n\n
        """
    
        translate_platform_name = {'metop01': 'Metop-B',
                                   'metop02': 'Metop-A',
                                   'metop03': 'Metop-C'}
        translate_channel_name = {'avhrr-3' : {'1':'1',
                                               '2':'2',
                                               '3a':'3',
                                               '4':'4',
                                               '5':'5',
                                               '3b':'6'}
                                  }

        if kwargs['platform_name']:
            _platform_name = translate_platform_name.get(kwargs['platform_name'],kwargs['platform_name'])

        _image_description = ''
        _image_description.encode('utf-8')

        _image_description += ' Satellite: '
        if ( _platform_name != None ):
            _image_description += _platform_name
    
        _image_description += '\n'
        
        _image_description += ' Date and Time: '
        #Select earliest start_time
        first = True
        earliest = 0
        for dataset in datasets:
            if first:
                earliest = dataset.info['start_time']
            else:
                if dataset.info['start_time'] < earliest:
                    earliest = datset.info['start_time']
            first=False
        print earliest
        _image_description += earliest.strftime("%H:%M %d/%m-%Y\n")
       
        _image_description += ' SatDir: 0\n'
    
        _image_description += ' Channels: '
        if type(kwargs["sensor"]) not in (tuple, list, set):
            kwargs["sensor"] = [kwargs["sensor"]]

        print "datasets in make_image_desc: {}".format(datasets)

        _image_description += str(len(datasets[0].info['prerequisites']))

        _image_description += ' In this file: '
        #tcn = translate_channel_name.get(kwargs['sensor'][0])

        #for dataset in datasets:
        for ch in datasets[0].info['prerequisites']:
            print ch
            try:
                _image_description += datasets[0].info['metadata_requirements'][ch]['alias']
            except KeyError:
                _image_description += ch

            _image_description += ' '
        
        #Replace last char(space) with \n
        _image_description = _image_description[:-1]
        _image_description += '\n'
       
        _image_description += ' Xsize: '
        _image_description += str(dataset.shape[0]) + '\n'
    
        _image_description += ' Ysize: '
        _image_description += str(dataset.shape[0]) + '\n'
    
        _image_description += ' Map projection: Stereographic\n'
        _image_description += ' Proj string: ' + datasets[0].info['area'].proj4_string
        if not all( datum in datasets[0].info['area'].proj4_string for datum in ['datum','towgs84']):
            _image_description += ' +towgs84=0,0,0'

        if not 'units' in datasets[0].info['area'].proj4_string:
            _image_description += ' +units=km'
        
        #Need to use center of lower left pixel. Subtract half a pixel size
        #image_description += ' +x_0=%.6f' % (-datasets[0].info['area'].area_extent[0]-datasets[0].info['area'].pixel_size_x/2.)
        #image_description += ' +y_0=%.6f' % (-datasets[0].info['area'].area_extent[1]-datasets[0].info['area'].pixel_size_y/2.)
        _image_description += ' +x_0=%.6f' % (-datasets[0].info['area'].area_extent[0]+datasets[0].info['area'].pixel_size_x)
        _image_description += ' +y_0=%.6f' % (-datasets[0].info['area'].area_extent[1]+datasets[0].info['area'].pixel_size_y)
    
        _image_description += '\n'
        _image_description += ' TrueLat: 60N\n'
        _image_description += ' GridRot: 0\n'
    
        _image_description += ' Xunit:1000 m Yunit: 1000 m\n'

        _image_description += ' NPX: %.6f' % (0)
        _image_description += ' NPY: %.6f' % (0) + '\n'

        _image_description += ' Ax: %.6f' % (datasets[0].info['area'].pixel_size_x/1000.)
        _image_description += ' Ay: %.6f' % (datasets[0].info['area'].pixel_size_y/1000.)
        #But this ads up to upper left corner of upper left pixel.
        #But need to use the center of the pixel. Therefor use the center of the upper left pixel.
        _image_description += ' Bx: %.6f' % (datasets[0].info['area'].area_extent[0]/1000. + datasets[0].info['area'].pixel_size_x/1000./2.) #LL_x
        _image_description += ' By: %.6f' % (datasets[0].info['area'].area_extent[3]/1000. - datasets[0].info['area'].pixel_size_y/1000./2.) #UR_y
        _image_description += '\n'
    
        LOG.debug("Area extent: {}".format(datasets[0].info['area'].area_extent))

        for ch in datasets[0].info['prerequisites']:
            found_channel = False
            print ch
                    
            palette=False
            #Make calibration.
            if palette:
                raise NotImplementedError("Mitiff palette saving is not implemented.")
            else:
                _image_description += 'Table_calibration: '
                try:
                    _image_description += datasets[0].info['metadata_requirements'][ch]['alias']
                except KeyError:
                    _image_description += ch

                _reverse_offset = 0.;
                _reverse_scale = 1.;

                #FIXME need to correlate the configured calibration and the calibration for the dataset.
                if datasets[0].info['metadata_requirements'][ch]['calibration'] == 'RADIANCE':
                    raise NotImplementedError("Mitiff radiance calibration not implemented.")
                    #_image_description += ', Radiance, '
                    #_image_description += '[W/m²/µm/sr]'
                    #_decimals = 8
                elif datasets[0].info['metadata_requirements'][ch]['calibration'] == 'brightness_temperature':
                    _image_description += ', BT, '
                    #_image_description += u'[\u2103]'
                    _image_description += u'[C]'

                    _reverse_offset = 255.;
                    _reverse_scale = -1.;
                    _decimals = 2
                elif datasets[0].info['metadata_requirements'][ch]['calibration'] == 'reflectance':
                    _image_description += ', Reflectance(Albedo), '
                    _image_description += '[%]'
                    _decimals = 2

                else:
                    LOG.warning("Unknown calib type. Must be Radiance, Reflectance or BT.")

                #How to format string by passing the format
                #http://stackoverflow.com/questions/1598579/rounding-decimals-with-new-python-format-function
            
                _image_description += ', 8, [ '
                for val in range(0,256):
                    #Comma separated list of values
                    #calib.append(boost::str(boost::format("%.8f ") % (prod_chan_it->min_val + (val * (prod_chan_it->max_val - prod_chan_it->min_val)) / 255.)));
                    _image_description += '{0:.{1}f} '.format((float(datasets[0].info['metadata_requirements'][ch]['min-val']) + ( (_reverse_offset + _reverse_scale*val) * ( float(datasets[0].info['metadata_requirements'][ch]['max-val']) - float(datasets[0].info['metadata_requirements'][ch]['min-val'])))/255.),_decimals)
                
                _image_description += ']\n\n'
                    
        return _image_description

    def _generate_filename(self, datasets, **kwargs):
        """Generate filename for the mitiff file
           I think this need a config input like a trollsift config or something
        """
        from trollsift.parser import compose
        filename = None
        if type(kwargs["sensor"]) not in (tuple, list, set):
            kwargs["sensor"] = [kwargs["sensor"]]

        _info = datasets[0].info
        _info.update({'area_id':datasets[0].info['area'].area_id})
        try:
            #filename = os.path.join(kwargs['output_dir'], compose(self.config[kwargs['sensor'][0]][0]['file-name'],datasets[0].info))
            filename = os.path.join(kwargs['output_dir'], compose(self.config[kwargs['sensor'][0]][0]['file-name'],_info))
        except:
            LOG.error("Failed to compose filename for sensor: {} with config: {} and output_dir: {}".format(kwargs['sensor'][0],self.config[kwargs['sensor'][0]][0]['file-name'], kwargs['output_dir']))
            raise

        return filename

    def _save_datasets_as_mitiff(self, datasets, image_description, gen_filename, **kwargs):
        """Put all togehter and save as a tiff file with the special tag making it a 
           mitiff file.
        """

        if type(kwargs["sensor"]) not in (tuple, list, set):
            kwargs["sensor"] = [kwargs["sensor"]]

        from libtiff import TIFF

        tif = TIFF.open(gen_filename, mode ='w')
        
        tif.SetField(IMAGEDESCRIPTION, str(image_description))
        
        #for ch in self.config[kwargs['sensor'][0]][0]['channels']:
        for i,ch in enumerate(datasets[0].info['prerequisites']):
            found_channel = False
            print ch
            #for dataset in datasets:
            #if ch['name'] == dataset.info['name']:
            #datasets[0].info['metadata_requirements'][ch]['alias']
            #Need to scale the data set to mitiff 0-255. 0 is no/missing/bad data.
            LOG.debug("min %f max %f value" % (float(datasets[0].info['metadata_requirements'][ch]['min-val']),float(datasets[0].info['metadata_requirements'][ch]['max-val'])))
            reverse_offset = 0.
            reverse_scale = 1.
            if datasets[0].info['metadata_requirements'][ch]['calibration'] == 'brightness_temperature':
            #if ch['calibration'] == "brightness_temperature":
                reverse_offset = 255.
                reverse_scale = -1.
                datasets[0].data[i] += KELVIN_TO_CELSIUS
                #print "after: ",datasets[0].data[i]
                
            LOG.debug("Reverse offset: %f reverse scale: %f" % ( reverse_offset,reverse_scale))

            _mask = datasets[0].mask[i]
            _data = np.clip(datasets[0].data[i], float(datasets[0].info['metadata_requirements'][ch]['min-val']),float(datasets[0].info['metadata_requirements'][ch]['max-val']))

            data=reverse_offset + reverse_scale*(((_data-float(datasets[0].info['metadata_requirements'][ch]['min-val']))/(float(datasets[0].info['metadata_requirements'][ch]['max-val']) - float(datasets[0].info['metadata_requirements'][ch]['min-val'])))*255.)

            data[_mask] = 0

            tif.write_image(data.astype(np.uint8), compression='deflate')
            found_channel = True
            #break

            if not found_channel:
                LOG.debug("Could not find configured channel in read data set. Fill with empty.")
                try:
                    fill_channel = np.zeros(datasets[0].data.shape,dtype=np.uint8)
                except IndexError as ie:
                    logger.error("Index out out bounds for channels[0].data.shape. Try area instead ... {}".format(ie))
                    fill_channel = np.zeros(datasetnotcorrect.area.shape,dtype=np.uint8)
                
                #import pdb
                #pdb.set_trace()

                tif.write_image(fill_channel, compression='deflate')
        tif.close

    def _find_config(self, dataset, **kwargs):

        for config_ch in self.config[kwargs['sensor'][0]]:
            if config_ch['type'] == 'metno-diana':
                for ch in config_ch['channels']:
                    if ch['name'] == dataset.info['name']:
                        return ch;
                else:
                    LOG.debug("product type not handeled {}".format(config_ch['type']))
        return None
