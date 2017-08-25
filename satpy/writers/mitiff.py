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

import pymitiff.mitiff as mt

from satpy.utils import ensure_dir
from satpy.writers import ImageWriter

LOG = logging.getLogger(__name__)

class MITIFFWriter(ImageWriter):

    def __init__(self, floating_point=False, tags=None, **kwargs):
        ImageWriter.__init__(self,
                             default_config_filename="writers/mitiff.cfg",
                             **kwargs)

        self.tags = self.config_options.get("tags",
                                            None) if tags is None else tags
        if self.tags is None:
            self.tags = {}
        elif not isinstance(self.tags, dict):
            # if it's coming from a config file
            self.tags = dict(tuple(x.split("=")) for x in self.tags.split(","))

    def save_image(self, img, filename=None, floating_point=False, **kwargs):
        """Save the image to the given *filename* in mitiff format.
        `floating_point` allows the saving of
        'L' mode images in floating point format if set to True.
        """

        LOG.debug("About to start save_image in mitiff")
        filename = filename or self.get_filename(**img.info)
        mt.save(img, filename, **kwargs)
        LOG.debug("save_image complete in mitiff")

