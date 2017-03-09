#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
# CoCoPy - A python toolkit for rotational spectroscopy
#
# Copyright (c) 2013 by David Schmitz (david.schmitz@chasquiwan.de).
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of 
# this software and associated documentation files (the “Software”), to deal in the 
# Software without restriction, including without limitation the rights to use, 
# copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
# Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH 
# THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
# MIT Licence (http://mit-license.org/)
#
################################################################################

import visa
import numpy as np

def get_resource_manager():
    """Return the PyVISA Resource Manager, creating an instance if necessary.
    :rtype: visa.ResourceManager
    """
    global _resource_manager
    if _resource_manager is None:
        _resource_manager = visa.ResourceManager()
    return _resource_manager


class visaInstrument():
    
    MANUFACTURER_ID = None
    MODEL_CODE = None
    __resource_manager = None
     
    def __init__(self, resource_name):

        self.__resource_manager = get_resource_manager()
        try:
            resource_info = self.__resource_manager.resource_info(resource_name)
        except visa.VisaIOError:
            raise ValueError('The resource name is invalid')

        self.resource_name = resource_name

        # The resource will be created when the driver is initialized.
        #: :type: pyvisa.resources.MessageBasedResource
        self.resource = None
        self.timeout = 5000
        self.encoding = 'utf-8'
        

    def initialize(self):
        self.resource = get_resource_manager().open_resource(self.resource_name)
        self.timeout = self.resource.timeout
        self.encoding = self.resource.encoding
        
    def finalize(self):
        self.resource.close()

    def query(self, message):
        return self.resource.query(message, delay=self.delay)

    def query_values(self, message):
        return self.resource.query(message, delay=self.delay)

    def query_waveform(self, message):
        return self.resource.query_ascii_values(message, converter=u'f', separator=u', ', container=list(), delay=None)

    def parse_query(self, message, parser):
        h = self.query(message)
        if format:
            h = parser(h)
        return h

    def write(self, message, termination=None, encoding=None):
        return self.resource.write(message, termination, encoding)

    def write_ascii_values(self, message, termination=None, encoding=None):
        return self.resource.write(message, termination, encoding)

    def write_values(self, message, termination=None, encoding=None):
        return self.resource.write(message, termination, encoding)

    def clear(self):
        self.resource.clear()
        
    def set_encoding(self, encoding):
        self.encoding = encoding
        
    def set_delay(self, delay):
        self.delay = delay

    def read(self, termination=None, encoding=None):
        return self.resource.read(termination, encoding)
        
    def read_ascii_values(self, termination=None, encoding=None):
        return self.resource.read(termination, encoding)
        
    def read_values(self, termination=None, encoding=None):
        return self.resource.read(termination, encoding)


class NiDAQmx():
    def __init__():
        