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

import numpy as np
import matplotlib.pyplot as plt

'''
Todo:
    - 1 Implement Coadd
    - 2 123
    - 3 
'''
class chirp:
    '''
    The chirp class tries to cover everything, which is related to creating different
    types of waveforms, which can be used with the AWG
    '''
    
    def __init__(self, sampling_rate=24.0, fname = 'unnamed.txt'):
        '''
        Starts an empty chirp instance.

        Parameters
        ----------
        fname (str): Filename of waveform file ('unnamed.txt').
        sampling_rate (float): sampling rate of the produced waveform in GS/s (24.0)
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        self.length = [1.E-6]
        self.sampling_rate = sampling_rate * 1.0E9
        self.delta_t = 1./self.sampling_rate
        self.s_freq = [2000.0E6]
        self.e_freq = [8000.0E6]
        self.mono_freq = [400.0E6]
        self.total_length = 1.E-6
        
        self.t_abs = np.array([])
        self.sample_points = 1
        
        self.marker1 = np.array([])
        self.marker2 = np.array([])

        self.sig = np.array([])
        
        self.fft_sig = np.array([])
        self.fft_freq = np.array([])
        
        self.fname = 'unnamed.txt'
        
        self.fname = fname
        
                
    def save(self, fname=''):
        '''
        Saves the waveform as an textfile. Creates the markers first.

        Parameters
        ----------
        fname (str): Filename of waveform file ('').
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        if fname != '':
            self.fname = fname
            
        self.marker()
        
        output = np.column_stack((self.sig, self.marker1, self.marker2))
        np.savetxt(self.fname, output, fmt='%.18e,%d,%d')

    def chirp(self, length = [1.0], s_freq = [2000.0], e_freq = [8000.0], phase=[0]):
        '''
        Creates a chirp. If a list is used for the input parameters a stacked chirp
        will be produced. Use zero for start and end frequency if you want to
        introduce zero padding. All lists need to be the same length.
        Updates the values: self.sig, self.sample_points, self.t_abs, self.length,
        self.s_freq, self.e_freq, self.phase

        Parameters
        ----------
        length (list, float): length of the different parts of the chirp in mus ([1.0])
        s_freq (list, float): starting frequency of the chirp in MHz ([2000.0])
        e_freq (list, float): ending frequency of the chirp in MHz ([8000.0])
        phase (list, float): starting phase of the chirp in radians ([0.0])
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        if type(length) != list:
            length = list([length])
        if type(s_freq) != list:
            s_freq = list([s_freq])
        if type(e_freq) != list:
            e_freq = list([e_freq])
        if type(phase) != list:
            phase = list([phase])
            
        self.length = np.array(length) * 1.0E-6
        self.s_freq = np.array(s_freq) * 1.0E6
        self.e_freq = np.array(e_freq) * 1.0E6
        self.phase = phase
        
        if len(self.s_freq) != len(self.e_freq) or len(self.s_freq) != len(self.length):
            print('s_freq, e_freq and length need to have the same length')
        else:
            i = 0
            sig = np.array([])
            for x in self.length:
                t_abs = np.arange(0, x, self.delta_t)
                sample_points = len(t_abs)
                if self.s_freq[i] == 0.0 and self.e_freq[i] == 0:
                    sig = np.hstack((sig, np.zeros(sample_points, dtype=np.float)))
                else:
                    sig = np.hstack((sig, np.sin((self.s_freq[i] + (self.e_freq[i]-self.s_freq[i])/self.length[i] * t_abs * .5) * 2.0 * np.pi * t_abs + self.phase[i])))
                i += 1
            
            self.sig = sig
            self.sample_points = len(sig)
            self.t_abs = np.arange(0, len(sig)) * self.delta_t
            self.total_length=len(self.sig)*self.delta_t

    def single(self, length = [1.0], mono_freq = [2000.0], phase=[0.0]):
        '''
        Creates a single frequency pulse. If a list is used for the input parameters 
        a stacked single frequency pulse will be produced. Use zero for the frequency 
        if you want to introduce zero padding. All lists need to be the same length.
        Updates the values: self.sig, self.sample_points, self.t_abs, self.length,
        self.s_freq, self.e_freq, self.phase

        Parameters
        ----------
        length (list, float): length of the different parts of the chirp in mus ([1.0])
        mono_freq (list, float): frequency of the pulse in MHz ([2000.0])
        phase (list, float): starting phase of the pulse in radians ([0.0])
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        if type(length) != list:
            length = list([length])
        if type(mono_freq) != list:
            mono_freq = list([mono_freq])
        if type(phase) != list:
            phase = list([phase])
        
        self.length = np.array(length) * 1.0E-6
        self.mono_freq = np.array(mono_freq) * 1.0E6
        self.phase = phase
        
        if len(self.mono_freq) != len(self.length):
            print('mono_freq and length need to have the same length')
        else:
            i = 0
            sig = np.array([])
            for x in self.length:
                t_abs = np.arange(0, x, self.delta_t)
                sample_points = len(t_abs)
                if self.mono_freq[i] == 0.0:
                    sig = np.hstack((sig, np.zeros(sample_points, dtype=np.float)))
                else:
                    sig = np.hstack((sig, np.sin(t_abs * 2 * np.pi * self.mono_freq[i] + phase[i])))
                
                i += 1
            
            self.sig = sig
            self.sample_points = len(sig)
            self.t_abs = np.arange(0, len(sig)) * self.delta_t
            self.total_length=len(self.sig)*self.delta_t
            
    def cut(self, sample_points=0, length=0.0):
        '''
        Cuts the signal to a certain number of sample_points or to a certain length

        Parameters
        ----------
        sample_points (int): number of sample points the signal will be cut to (0)
        length (float): length of signal in mus the original signal will be cut to (0.0)
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        if sample_points != 0 or length != 0.0:
            if sample_points != 0:
                h = sample_points
            elif length != 0.0:
                h = (length*1.E-6)/self.delta_t
            
            self.sig = self.sig[0:h]
            self.t_abs = self.t_abs[0:h]
            self.sample_points=len(self.sig)
            self.total_length=len(self.sig)*self.delta_t
            
    def plot(self, img = 1):
        '''
        Plots the pulse and the FFT of it.

        Parameters
        ----------
        img (int): plot image tag (1)
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        self.fft()                
        plt.ion()
        plt.figure(img)
        plt.subplots_adjust(wspace=1, hspace=1)
        
        plt.subplot(211)
        plt.title('Chirped Pulse in the Time Domain', size=11)
        plt.ylabel('Intensity [a. U.]')
        plt.xlabel('Time [mus]')
        plt.plot(self.t_abs*10**6, self.sig, '-')
                
        plt.subplot(212)
        plt.title('Chirped Pulse in the Frequency Domain', size=11)
        plt.ylabel('Intensity [a. U.]')
        plt.xlabel('Frequency [GHz]')
        plt.axis([0, max(self.fft_freq)*10**-9, 0, max(self.fft_sig)])
        plt.plot(self.fft_freq*10**-9, self.fft_sig)
        
        plt.show()

                
    def invert(self):
        '''
        Inverts the signal/pulse.

        Parameters
        ----------
        none
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        self.sig = self.sig * -1
             
        
    def fft(self):
        '''
        Performs the FFT of the signal/pulse.

        Parameters
        ----------
        none
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        h_sig = np.absolute(np.fft.rfft(self.sig)) ** 2 / len(self.sig)
        self.fft_sig = h_sig[0:int(len(self.sig)/2)]
        h_freq = np.fft.fftfreq(self.sample_points, self.delta_t)
        self.fft_freq = h_freq[0:int(len(self.sig)/2)]
        
        
    def marker(self, marker_up=0.01, marker_down=0.002):
        '''
        Creates the marker.

        Parameters
        ----------
        marker_up (float): time in mus the marker is up == 1 (0.01)
        marker_down (float): time in mus the marker is down == 0, after the up
        time (0.002)
        
        Returns
        -------
        none
        
        Notes
        -----
        none
        
        '''
        if marker_up > marker_down:
            start_down = int((marker_down*1.E-6)/self.delta_t)
            start_up = int((marker_up*1.E-6)/self.delta_t) + start_down
            
            self.marker1 = np.zeros(self.sample_points, dtype=np.int)
            self.marker2 = np.zeros(self.sample_points, dtype=np.int)
            
            self.marker1[self.sample_points-start_up:self.sample_points-start_down] = 1
            self.marker2[self.sample_points-start_up:self.sample_points-start_down] = 1
