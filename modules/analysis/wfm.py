#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
# CoCoPy - A python toolkit for rotational spectroscopy
#
# Copyright (c) 2015 by David Schmitz (david.schmitz@chasquiwan.de).
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

'''
Todo:
    - 1 Write comments
    - 2
    - 3
'''

import numpy as np
import struct as su
import glob as gl
import os

class fidReadoutWfm:
    '''
    The fidReadoutWfm can be used to read out Tektronix waveform files (vers. 3)
    '''

    meas_params = dict(N_points=0, max_t=0.0, delta_t=0.0, fname='', name='',
                       trigger_dl=0.0)
    fid = np.array([[]])

    def __init__(self, fname = 'unnamed.wfm'):
        '''
        Readout of a Tektronix waveform file. Filename needs to be specified.

        Parameters
        ----------
        fname (str): Filename of waveform file ('unnamed.txt').

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        self.meas_params['fname'] = fname
        try:
            if fname == '':
                fname = self.meas_params['fname']
            else:
                self.meas_params['fname'] = fname
            fd = open(fname, 'rb')
        except IOError:
            print 'Cannot open: ', self.meas_params['fname']
        else:

            # Waveform static information

            params = dict();

            params['ByteOrder']                 =  su.unpack('H', fd.read(2))[0]
            params['VersionNum']                =  su.unpack('8s', fd.read(8))[0]
            params['NumDigitsInByteCount']      =  su.unpack('b', fd.read(1))[0]
            params['NumBytesToEOF']             =  su.unpack('i', fd.read(4))[0]
            params['NumBytesPerPoint']          =  su.unpack('B', fd.read(1))[0]
            params['ByteOffsetToCurveBuffer']   =  su.unpack('i', fd.read(4))[0]
            params['HorZoomScale']              =  su.unpack('i', fd.read(4))[0]
            params['HorZoomPos']                =  su.unpack('f', fd.read(4))[0]
            params['VerZoomScale']              =  su.unpack('d', fd.read(8))[0]
            params['VerZoomPos']                =  su.unpack('f', fd.read(4))[0]
            params['WaveformLabel']             =  su.unpack('32s', fd.read(32))[0]
            params['N']                         =  su.unpack('I', fd.read(4))[0]
            params['HeaderSize']                =  su.unpack('H', fd.read(2))[0]

            # WFM header

            header = dict()

            header['SetType']                   =  su.unpack('i', fd.read(4))[0]
            header['WfmCnt']                    =  su.unpack('I', fd.read(4))[0]
            jnk                                 =  fd.read(36)
            header['DataType']                  =  su.unpack('i', fd.read(4))[0]
            jnk                                 =  fd.read(8)
            header['AccuWfmCnt']                =  su.unpack('I', fd.read(4))[0]
            jnk                                 =  fd.read(30)

            # Explicit Dimension 1/2

            ed = dict()

            for n in range(0,2):

                ed['DimScale' + str(n+1)]       =  su.unpack('d', fd.read(8))[0]
                ed['DimOffset' + str(n+1)]      =  su.unpack('d', fd.read(8))[0]
                ed['DimSize' + str(n+1)]        =  su.unpack('I', fd.read(4))[0]
                ed['Units' + str(n+1)]          =  su.unpack('20s', fd.read(20))[0]
                ed['DimExtentMin' + str(n+1)]   =  su.unpack('d', fd.read(8))[0]
                ed['DimExtentMax' + str(n+1)]   =  su.unpack('d', fd.read(8))[0]
                ed['DimResolution' + str(n+1)]  =  su.unpack('d', fd.read(8))[0]
                ed['DimRefPoint' + str(n+1)]    =  su.unpack('d', fd.read(8))[0]
                ed['Format' + str(n+1)]         =  su.unpack('i', fd.read(4))[0]
                ed['StorageType' + str(n+1)]    =  su.unpack('i', fd.read(4))[0]
                jnk                             =  fd.read(20)
                ed['UserScale' + str(n+1)]      =  su.unpack('d', fd.read(8))[0]
                ed['UserUnits' + str(n+1)]      =  su.unpack('20s', fd.read(20))[0]
                ed['UserOffset' + str(n+1)]     =  su.unpack('d', fd.read(8))[0]
                ed['PointDensity' + str(n+1)]   =  su.unpack('d', fd.read(8))[0]
                ed['HRef' + str(n+1)]           =  su.unpack('d', fd.read(8))[0]
                ed['TrigDelay' + str(n+1)]      =  su.unpack('d', fd.read(8))[0]

            # Implicit Dimension 1/2

            md = dict()

            for n in range(0,2):

                md['DimScale' + str(n+1)]       =  su.unpack('d', fd.read(8))[0]
                md['DimOffset' + str(n+1)]      =  su.unpack('d', fd.read(8))[0]
                md['DimSize' + str(n+1)]        =  su.unpack('I', fd.read(4))[0]
                md['Units' + str(n+1)]          =  su.unpack('20s', fd.read(20))[0]
                jnk                             =  fd.read(16)
                md['DimResolution' + str(n+1)]  =  su.unpack('d', fd.read(8))[0]
                md['DimRefPoint' + str(n+1)]    =  su.unpack('d', fd.read(8))[0]
                md['Spacing' + str(n+1)]        =  su.unpack('I', fd.read(4))[0]
                md['UserScale' + str(n+1)]      =  su.unpack('d', fd.read(8))[0]
                md['UserUnits' + str(n+1)]      =  su.unpack('20s', fd.read(20))[0]
                md['UserOffset' + str(n+1)]     =  su.unpack('d', fd.read(8))[0]
                md['PointDesity' + str(n+1)]    =  su.unpack('d', fd.read(8))[0]
                md['HRef' + str(n+1)]           =  su.unpack('d', fd.read(8))[0]
                md['TrigDelay' + str(n+1)]      =  su.unpack('d', fd.read(8))[0]

            # Time Base 1/2 Information

            tb = dict()

            for n in range(0,2):

                tb['RPointSpacing' + str(n+1)]  =  su.unpack('I', fd.read(4))[0]
                tb['Sweep' + str(n+1)]          =  su.unpack('i', fd.read(4))[0]
                tb['TypeOfBase' + str(n+1)]     =  su.unpack('i', fd.read(4))[0]

            # Waveform Update Specification

            UpSpecs = dict()

            UpSpecs['RPOffset']                 =  su.unpack('I', fd.read(4))[0]
            UpSpecs['TTOffset']                 =  su.unpack('d', fd.read(8))[0]
            UpSpecs['FracSec']                  =  su.unpack('d', fd.read(8))[0]
            UpSpecs['GmtSec']                   =  su.unpack('i', fd.read(4))[0]

            # WFM Curve Information

            CurveInf = dict()

            jnk                                 =  fd.read(10)
            CurveInf['PrechargeStartOffset']    =  su.unpack('I', fd.read(4))[0]
            CurveInf['DataStartOffset']         =  su.unpack('I', fd.read(4))[0]
            CurveInf['PostchargeStartOffset']   =  su.unpack('I', fd.read(4))[0]
            CurveInf['PostchargeStopOffset']    =  su.unpack('I', fd.read(4))[0]
            CurveInf['EndOfCurveBufferOffset']  =  su.unpack('I', fd.read(4))[0]

            # Get the Curve Size in Byte and Points

            params['CurveSizeInBytes']          = CurveInf['PostchargeStartOffset'] - CurveInf['DataStartOffset']
            params['CurveSize']                 = params['CurveSizeInBytes'] / int(params['NumBytesPerPoint'])
            jnk                                 = fd.read(CurveInf['DataStartOffset'])

            # Extract the curve to an array
            if params['NumBytesPerPoint'] == 1:
                curve                           = np.array(su.unpack(str(params['CurveSize']) + 'b', fd.read(params['CurveSizeInBytes'])))
            elif params['NumBytesPerPoint'] == 2:
                curve                           = np.array(su.unpack(str(params['CurveSize']) + 'h', fd.read(params['CurveSizeInBytes'])))
            elif params['NumBytesPerPoint'] == 4:
                curve                           = np.array(su.unpack(str(params['CurveSize']) + 'f', fd.read(params['CurveSizeInBytes'])))
            elif params['NumBytesPerPoint'] == 8:
                curve                           = np.array(su.unpack(str(params['CurveSize']) + 'd', fd.read(params['CurveSizeInBytes'])))
            else:
                print 'Measurement file corrupted. Contact your locale Tektronix vendor and get the things fixed'
            # Create the Time Vector and merge both vectors to an array
            # Create the Time Vector and merge both vectors to an array

            y = ed['DimOffset1'] + ed['DimScale1'] * curve
            t = md['DimScale1'] * np.arange(0, params['CurveSize'])

            y = np.reshape(y, (len(y), 1))
            t = np.reshape(t, (len(t), 1))

            self.fid = np.hstack((t, y))
            self.meas_params['N_points'] = params['CurveSize']
            self.meas_params['max_t'] = t[len(t)-1, 0]
            self.meas_params['delta_t'] = t[1, 0]-t[0, 0]
            self.meas_params['trigger_dl'] = md['DimOffset1']
            self.meas_params['averages'] = header['AccuWfmCnt']

            fd.close()

class fidReadoutTxt:
    '''
    help
    '''
    meas_params = dict(N_points=0, max_t=0.0, delta_t=0.0, fname='', name='',
                       trigger_dl=0.0)
    fid = np.array([[]])

    def __init__(self, fname = 'unnamed.txt', header=6):
        '''
        Readout of a FID text waveform file (two column text file). Filename
        needs to be specified.

        Parameters
        ----------
        fname (str): Filename of waveform file ('unnamed.txt').

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        self.meas_params['fname'] = fname
        try:
            if fname == '':
                fname = self.meas_params['fname']
            else:
                self.meas_params['fname'] = fname
            f = open(fname, 'r+')
            fid_h = np.zeros((header,2))
            i = 0
            for line in f:
                fid_h[i,0] = float(line.split('\t')[-2])
                fid_h[i,1] = float(line.split('\t')[-1])
                i += 1
                if i == header:
                    break
            f.close()

        except IOError:
            print 'Cannot open: ', self.meas_params['fname']
        else:
            fid = np.loadtxt(self.meas_params['fname'], skiprows=header)
            fid = np.vstack((fid_h, fid))

            self.meas_params['N_points'] = len(fid)
            self.meas_params['max_t'] = fid[len(fid)-1, 0]
            self.meas_params['delta_t'] = fid[1, 0]-fid[0, 0]
            self.meas_params['trigger_dl'] = fid[0, 0]
            self.fid = fid

class specReadoutTxt:
    '''
    help
    '''
    meas_params = dict(N_points=0, min_f=0.0, max_f=0.0, fname='', name='',
                       delta_f=0.0)
    spec = np.array([[]])

    def __init__(self, fname = 'unnamed.txt'):
        '''
        Readout of a spectrum text file (two column text file). Filename needs
        to be specified.

        Parameters
        ----------
        fname (str): Filename of waveform file ('unnamed.txt').

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        self.meas_params['fname'] = fname
        try:
            if fname == '':
                fname = self.meas_params['fname']
            else:
                self.meas_params['fname'] = fname
            fd = open(fname, 'r')
            fd.close()
        except IOError:
            print 'Cannot open: ', self.meas_params['fname']
        else:
            spec = np.genfromtxt(self.meas_params['fname'])

            self.meas_params['N_points'] = len(spec)
            self.meas_params['max_f'] = spec[len(spec)-1, 0]
            self.meas_params['delta_f'] = spec[1, 0]-spec[0, 0]
            self.meas_params['min_f'] = spec[0, 0]
            self.spec = spec

class linesReadoutTxt:
    '''
    help
    '''
    meas_params = dict(nlines=0, min_f=0.0, max_f=0.0, fname='', name='')
    lines = np.array([[]])

    def __init__(self, fname = 'unnamed.txt'):
        '''
        Readout of a linelist text file (one or two column text file). Filename
        needs to be specified.

        Parameters
        ----------
        fname (str): Filename of waveform file ('unnamed.txt').

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        self.meas_params['fname'] = fname
        try:
            if fname == '':
                fname = self.meas_params['fname']
            else:
                self.meas_params['fname'] = fname
            fd = open(fname, 'r')
            fd.close()
        except IOError:
            print 'Cannot open: ', self.meas_params['fname']
        else:
            lines = np.loadtxt(self.meas_params['fname'])
            if len(lines.shape) == 1:
                h = np.ones(len(lines))
                h = h.reshape((len(h), 1))
                lines = np.hstack((lines.reshape((len(lines), 1)), h))

            self.meas_params['nlines'] = len(lines)
            self.meas_params['max_f'] = np.max(lines[:,0])
            self.meas_params['min_f'] = np.min(lines[:,0])
            self.lines = lines

def files_folder(ext='wfm', sort='time', key=''):
    '''
    Returns the files within a folder with a defined file extension 'ext' and
    sorted by time or name.

    Parameters
    ----------
    ext (str): File extension ('wfm')
    sort (str): Algorithm to use for sorting: 'time' or 'name' ('time')
    key (str): filename has to contain that substring

    Returns
    -------
    filelist (list, str): sorted filelist

    Notes
    -----
    none
    '''

    filelist = gl.glob('*.'+ext)
    filelist = [ x for x in filelist if key in x ]

    if sort == 'time':
        filelist.sort(key=lambda x: os.path.getmtime(x))
    elif sort == 'name':
        filelist.sort()

    return filelist
