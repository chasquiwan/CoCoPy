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

import re

'''
Todo:
    - 1 yo
    - 2
    - 3
'''

def readout(fname='', calctype='opt'):
    try:
        f = open(fname, 'r')
    except IOError:
        print 'Cannot open: ', fname
    else:
        if calctype == 'opt' or calctype == 'optfreq':
            dip = 0
            freq = 0
            rot = 0
            mull = 0
            longstr = ''
            flag = False
            properties = dict()
            for line in f:
                if 'This is the Gaussian(R) 03 program' in line:
                    ver = '03'
                if 'This is part of the Gaussian(R) 09 program.' in line:
                    ver = '09'
                if 'Rotational constants (GHZ):' in line and rot == 1:
                    properties['rotA'] = float(line[30:45]) * 1000.
                    properties['rotB'] = float(line[45:60]) * 1000.
                    properties['rotC'] = float(line[60:]) * 1000.
                    rot = 0
                if '---------------------------------------------------------------------' in line:
                    rot = 1
                if dip == 1:
                    if ver == '09':
                        properties['dipA'] = float(line[13:26])
                        properties['dipB'] = float(line[39:52])
                        properties['dipC'] = float(line[65:78])
                        properties['dipTot'] = float(line[91:])
                    if ver == '03':
                        properties['dipA'] = float(line[6:17])
                        properties['dipB'] = float(line[23:34])
                        properties['dipC'] = float(line[40:51])
                        properties['dipTot'] = float(line[57:])

                    dip = 0
                if ' Dipole moment (field-independent basis, Debye):' in line:
                    dip = 1
                if 'Harmonic frequencies (cm**-1)' in line:
                    freq = 1
                    properties['freq'] = list([])
                    properties['intens'] = list([])
                if freq == 1 and 'Frequencies --' in line:
                    properties['freq'].append(float(line[15:30]))
                    if len(line) > 30:
                        properties['freq'].append(float(line[30:55]))
                    if len(line) > 55:
                        properties['freq'].append(float(line[55:]))
                if freq == 1 and 'IR Inten    --' in line:
                    properties['intens'].append(float(line[15:30]))
                    if len(line) > 30:
                        properties['intens'].append(float(line[30:55]))
                    if len(line) > 55:
                        properties['intens'].append(float(line[55:]))
                if 'Mulliken atomic charges:' in line:
                    mull = 1
                if 'SCF Done:' in line:
                    properties['energy'] = float(line.split('=')[1].split('A.U.')[0])
                if flag == True:
                    if line != '\n':
                        longstr += line
                    else:
                        flag = False
                if '(Enter /usr/product/gaussian/g09/l9999.exe)' in line:
                    flag = True
            
            if '\MP2=' in longstr:
                h = re.findall('\\\\MP2=[+-]?\d*\.\d*\\\\', longstr)
                properties['energy'] = float(re.findall('[+-]?\d*\.\d*', h[0])[0])
            
            return properties
            
        if calctype == 'optscan':
            dip = 0
            freq = 0
            rot = 0
            properties = dict(energy=[], rotA=[], rotB=[], rotC=[])
            for line in f:
                if 'This is the Gaussian(R) 03 program' in line:
                    ver = '03'
                if 'This is part of the Gaussian(R) 09 program.' in line:
                    ver = '09'
                if 'Rotational constants (GHZ):' in line and rot == 1:
                    rotA = float(line[30:45]) * 1000.
                    rotB = float(line[45:60]) * 1000.
                    rotC = float(line[60:]) * 1000.
                    rot = 0
                if '---------------------------------------------------------------------' in line:
                    rot = 1
                if 'SCF Done:' in line:
                    if ver == '09':
                        energy = float(line[23:40])
                    if ver == '03':
                        energy = float(line[26:43])
                if '-- Stationary point found.' in line:
                    properties['energy'].append(energy)
                    properties['rotA'].append(rotA)
                    properties['rotB'].append(rotB)
                    properties['rotC'].append(rotC)
            return properties

        f.close()
            
        