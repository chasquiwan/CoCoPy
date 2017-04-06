#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
# CoCoPy - A python toolkit for rotational spectroscopy
#
# Copyright (c) 2016 by David Schmitz (david.schmitz@chasquiwan.de).
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
#from analysis.util import ProgressBar
import util
import itertools as it
import os
from tempfile import TemporaryFile
import scipy.constants as con
from scipy.integrate import ode
from scipy.misc import comb


'''
Todo:
    1. Comment all
    2. module load
'''

class Bloch:

    def __init__(self, prefix='test', levels=2, parts=2, dt=1.E-11):
        '''
        This class tries to cover everything, which is related to Bloch
        equations.

        Parameters
        ----------
        prefix (str): Prefix for all temporary filenames ('test')
        levels (int): number of energy levels involved
        parts (int): number of different seperations in time??
        dt (float): timestep for the evaluation

        Returns
        -------
        Bloch instance

        Notes
        -----
        none
        '''

        self.egy = {}
        self.dip = {}
        self.pop = {}
        self.om0 = {}

        self.pulse_scheme = []

        self.pulse = np.array([])
        self.freq = {}
        self.intens = {}
        self.phase = {}

        self.DEBYE = 3.335E-30
        self.ANGULAR = 2.*np.pi*1.E6

        self.dt = dt

        self.t0 = 0.
        self.rho0 = np.array([], dtype=complex)

        self.t = np.array([])
        self.rho = np.array([], dtype=complex)

        self.LABELS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        tr = int(comb(levels, 2))

        self.params = {
            'parts': parts, 'levels': levels, 'tr': tr, 'cut': False,
            'Npoints': 1000
            }

        self.pulse = np.zeros(parts)

        for i in range(levels):
            self.egy[self.LABELS[i]] = 0.
            self.pop[self.LABELS[i]] = 0.

        for i in self.get_keys():
            self.dip[i] = 0.0 + 0.0j
            self.freq[i] = np.zeros(parts)
            self.intens[i] = np.zeros(parts)
            self.phase[i] = np.zeros(parts)

    def set_starting_values(self, t0=0., rho0=np.array([])):
        self.t0 = t0
        if len(rho0):
            self.rho0 = rho0
        else:
            self.rho0 = np.zeros((self.params['levels'],self.params['levels']), dtype=complex)
            self.rho0[0,0] = 1.+0.0j

    def set_freq(self, freq=np.array([])):
        if len(freq):
            i = 0
            for x in freq:
                self.freq[self.get_keys()[i]] = np.ones(self.params['parts']) * x * self.ANGULAR
                i += 1

    def set_intens(self, intens=np.array([])):
        if len(intens):
            i = 0
            for x in intens:
                self.intens[self.get_keys()[i]] = np.ones(self.params['parts']) * x
                i += 1


    def set_phase(self, phase=np.array([])):
        if len(phase):
            i = 0
            for x in phase:
                self.phase[self.get_keys()[i]] = np.ones(self.params['parts']) * x
                i += 1

    def set_dip(self, dip=np.array([])):
        if len(dip):
            i = 0
            for x in dip:
                self.dip[self.get_keys()[i]] = x * self.DEBYE
                i += 1

    def set_egy(self, egy=np.array([]), unit='J'):
        if len(egy):
            i = 0; h = 0.; zero = 0.
            for x in egy:
                self.egy[self.LABELS[i]] = util.convert_energy((x - zero), unit, 'J')
                if h > x:
                    print 'Energies needs to be in ascending order'
                if i == 0:
                    zero = x
                h = x
                i += 1

            for x in self.get_keys():
                h_upper = self.egy[x[1].upper()]
                h_lower = self.egy[x[0].upper()]
                self.om0[x] = util.convert_energy((h_upper - h_lower), 'J', 'Hz') / 1.E6


    def set_pop(self, pop=np.array([])):
        if len(pop):
            i = 0
            for x in pop:
                self.pop[self.LABELS[i]] = x
                i += 1

    def set_pulse(self, pulse=np.array([])):
        if len(pulse):
            self.pulse = pulse

    def get_keys(self):
        keys = []
        for x in it.combinations(self.LABELS[:self.params['levels']], 2):
            keys.append((x[0]+x[1]).lower())

        return keys

    def get_indices(self):
        indices = []
        for x in it.combinations(range(self.params['levels']), 2):
            indices.append(x)

        return indices


    def set_init(self, i):
        arg = list([self.params['levels'], self.params['tr']])
        for x in self.LABELS[:self.params['levels']]:
            arg.append(self.egy[x])
        for x in self.get_keys():
            arg.append(self.dip[x])
        for x in self.get_keys():
            arg.append(self.intens[x][i])
        for x in self.get_keys():
            arg.append(self.freq[x][i])
        for x in self.get_keys():
            arg.append(self.phase[x][i])

        return arg


    def bloch_solver(self, Npoints=-1):
        #TODO: add Npoints in params
        r = ode(self.f).set_integrator('zvode')

        for i in range(self.params['parts']):

            # Set parameters and starting values at t0
            if i == 0:
                r.set_initial_value(self.rho0.flatten(), self.t0).set_f_params(self.set_init(i))
                t = [self.t0]
                rho = [self.rho0]
                tpp = 0.

            # Set parameters for the next pulse
            else:
                r.set_initial_value(r.y, r.t).set_f_params(self.set_init(i))
                tpp += self.pulse[i-1]

            r.integrate(r.t + self.dt)

            # Integrate and solve ODEs
            while r.successful() and r.t <= tpp + self.pulse[i]:
                r.integrate(r.t + self.dt)
                if Npoints == -1 or r.t > sum(self.pulse) - Npoints*self.dt:
                    rho.append(r.y.reshape((self.params['levels'], self.params['levels'])))
                    t.append(r.t)

        rho_h = np.zeros((self.params['levels'], self.params['levels'], len(rho)), dtype=complex)
        i = 0
        for x in rho:
            rho_h[:,:,i] = x
            i += 1

        #convert list to array
        if Npoints == -1:
            self.t = np.array(t)
            self.rho = rho_h
            return np.array(t), rho_h
        else:
            self.t = np.array(t[1:])
            self.rho = rho_h[:,:,1:]
            return np.array(t[1:]), rho_h[:,:,1:]


    def save_dump(self, rho, t, tmpfile_t, tmpfile_rho):
        tmpfile_t.seek(0); tmpfile_rho.seek(0)
        t_h = np.load(tmpfile_t); rho_h = np.load(tmpfile_rho)
        tmpfile_t.seek(0); tmpfile_rho.seek(0)

        rho_hh = np.zeros((self.params['levels'], self.params['levels'], len(rho)), dtype=complex)
        i = 0
        for x in rho:
            rho_hh[:,:,i] = x
            i += 1

        rho = np.concatenate((rho_h, rho_hh), axis=2)
        t = np.append(t_h, np.array(t))
        np.save(tmpfile_t, t); np.save(tmpfile_rho, rho)


    def bloch_solver_dump(self, Npoints=10000):

        tmpfile_t = TemporaryFile()
        tmpfile_rho = TemporaryFile()

        r = ode(self.f).set_integrator('zvode')

        for i in range(self.params['parts']):

            # Set parameters and starting values at t0
            if i == 0:
                r.set_initial_value(self.rho0.flatten(), self.t0).set_f_params(self.set_init(i))
                t = [self.t0]
                rho = [self.rho0]
                tpp = 0.

                np.save(tmpfile_rho, np.array(rho).reshape((3,3,1)))
                np.save(tmpfile_t, np.array(t))

            # Set parameters for the next pulse
            else:
                r.set_initial_value(r.y, r.t).set_f_params(self.set_init(i))
                tpp += self.pulse[i-1]

            r.integrate(r.t + self.dt)

            # Integrate and solve ODEs
            j = 0
            while r.successful() and r.t <= tpp + self.pulse[i]:
                r.integrate(r.t + self.dt)

                if j == Npoints:
                    self.save_dump(rho, t, tmpfile_t, tmpfile_rho)
                    rho = []
                    t = []
                    j = 0

                rho.append(r.y.reshape((self.params['levels'], self.params['levels'])))
                t.append(r.t)
                j += 1

        self.save_dump(rho, t, tmpfile_t, tmpfile_rho)

        tmpfile_t.seek(0); tmpfile_rho.seek(0)
        rho = np.load(tmpfile_rho);
        t = np.load(tmpfile_t);
        tmpfile_rho.close(); tmpfile_t.close()

        self.t = t
        self.rho = rho

        return t, rho


    def f(self, t, y, arg):

        # Make it human readable
        n   = arg[0]
        tr  = arg[1]
        E   = arg[2:n+2]
        dip = arg[n+2+tr*0:n+2+tr*1]
        F   = arg[n+2+tr*1:n+2+tr*2]
        o   = arg[n+2+tr*2:n+2+tr*3]
        phi = arg[n+2+tr*3:n+2+tr*4]

        # Setting up the Hamiltonian
        H = np.zeros((n,n), dtype=complex)
        for i in range(n):
            H[i,i] = E[i]

        # ab, ac, ad, .... bc, bd .... ,cd...
        i = 0
        for j in it.combinations(range(n), 2):
            H[j[0],j[1]] = - dip[i] * F[i] / 2. * (np.exp(-1.j*(o[i]*t+phi[i])) + np.exp(1.j*(o[i]*t+phi[i])))
            H[j[1],j[0]] = np.conjugate(H[j[0],j[1]])
            i += 1

        # Setting up the density matrix
        rho = y.reshape((n, n))

        # Calculate commutator
        dt_rho = -1.j/con.hbar * ( np.dot(H, rho) - np.dot(rho, H) )

        return dt_rho.flatten()

    def conv_phase(self, phase):
        if phase > np.pi:
            return -np.pi + np.abs(phase - np.pi)
        elif phase < -np.pi:
            return np.pi - np.abs(phase + np.pi)
        else:
            return phase


    def get_phase(self, trans='ab', freq=0., pol=True):

        if len(self.rho):
            t = self.t

            if freq == 0.:
                freq = self.om0[trans] * self.ANGULAR
            else:
                freq *= self.ANGULAR

            if pol == True:
                signal = self.get_pol(trans)[-100:]

            else:
                signal = self.get_coh(trans)[-100:].real

            t = t[-100:]
            ref = self.generate_sin(freq, 0., t)


            ts_zero = self.zero_crossing(t, signal)
            tr_zero = self.zero_crossing(t, ref)

            return self.conv_phase((tr_zero-ts_zero)*freq)


    def zero_crossing(self, t, signal):
        h=1.; t_cross = 0.
        for i in range(100):
            if np.sign(signal[i]*h) < 0.:
                tria = np.abs(signal[i])/(np.abs(h) + np.abs(signal[i]))
                if h < signal[i]:
                    t_cross = t[i] - tria*self.dt
            h = signal[i]

        if t_cross:
            return t_cross


    def get_intens(self, trans='ab', pol=True):

        if len(self.rho):
            t = self.t

            if pol == True:
                signal = self.get_pol(trans)[-100:]
            else:
                signal = self.get_coh(trans)[-100:].real

            return np.max(np.abs(signal))


    def generate_sin(self, freq, phase, t=np.array([])):
        if len(t) == 0:
            t = self.t
            return (np.exp(-1.j*(freq*t-np.pi/2.-phase)) + np.exp(1.j*(freq*t-np.pi/2.-phase))).real


    def get_coh(self, trans='ab'):
        h = np.where(np.array(self.get_keys()) == trans)[0]
        j = self.get_indices()[h]

        if len(self.rho):
            return self.rho[j[0],j[1],:]


    def get_pol(self, trans='ab'):
        if len(self.rho):
            return (self.get_coh(trans) * np.conjugate(self.dip[trans])).real


    def get_pop(self, level='A'):
        h = self.LABELS.index(level)

        if len(self.rho):
            return self.rho[h,h,:].real

    def generate_pulse_scheme(self):

        self.pulse = np.array([])

        for x in get_keys(self):
            self.freq[x] = np.array([])
            self.intens[x] = np.array([])
            self.phase[x] = np.array([])

        self.pulse_scheme = sorted(self.pulse_scheme)
        self.pulse = []

        for i in range(len(self.pulse_scheme)):
            for j in range(len(self.pulse_scheme))[i:]:
                if self.pulse_scheme[j][0] != self.pulse_scheme[i][0]:
                    self.pulse = np.append((self.pulse, [self.pulse_scheme[j][0] - self.pulse_scheme[i][0]]))

            if i == 0:
                for x in get_keys():
                    if x == self.pulse_scheme[i][1]:
                        self.freq[x] = np.append((self.freq[x], [self.pulse_scheme[i][2]]))
                        self.intens[x] = np.append((self.intens[x], [self.pulse_scheme[i][3]]))
                        self.phase[x] = np.append((self.phase[x], [self.pulse_scheme[i][4]]))

                    else:
                        self.freq[x] = np.append((self.freq[x], [0.]))
                        self.intens[x] = np.append((self.intens[x], [0.]))
                        self.phase[x] = np.append((self.phase[x], [0.]))

            elif self.pulse_scheme[i] == self.pulse_scheme[i-1]:
                for x in get_keys():
                    if x == self.pulse_scheme[i][1]:
                        self.freq[x][-1] = self.pulse_scheme[i][2]
                        self.intens[x][-1] = self.pulse_scheme[i][3]
                        self.phase[x][-1] = self.pulse_scheme[i][4]

            else:
                for x in get_keys():
                    if x == self.pulse_scheme[i][1]:
                        self.freq[x] = np.append((self.freq[x], [self.pulse_scheme[i][2]]))
                        self.intens[x] = np.append((self.intens[x], [self.pulse_scheme[i][3]]))
                        self.phase[x] = np.append((self.phase[x], [self.pulse_scheme[i][4]]))

                    else:
                        self.freq[x] = np.append((self.freq[x], [self.freq[x][-1]]))
                        self.intens[x] = np.append((self.intens[x], [self.intens[x][-1]]))
                        self.phase[x] = np.append((self.phase[x], [self.phase[x][-1]]))






            self.pulse = np.append((self.pulse, [self.pulse_scheme[i][0] - self.pulse_scheme[i-1][0]]))
            self.freq[transition]

    def single_pulse(self, length=1., freq=2000.):
        pass

    def zero_pulse(self, length=1.):
        pass

    def plot_pulse_scheme(self):
        pass

    def f_decay(self):
        pass

    def bloch_solver_decay(self):
        pass

    def rect_pulse(self, start, duration, frequency, transition, intensity, phase = 0):
        end = start + duration
        T = 1./(frequency * 1.E6)

        tmp = np.mod(duration, T)

        start = [start, transition, frequency, intensity, phase]
        end = [end, transition, 0., 0., 0.]

        self.pulse_scheme.append(start)
        self.pulse_scheme.append(end)

    def chirped_pulse(self, start, duration, start_frequency, end_frequency, transition, intensity, phase = 0):
        t = np.arange(0, duration + self.dt)
        freq = np.linspace(start_frequency, stop_frequency, len(t))

        for i in np.arange(len(t)):
            instant_phase = phase + 2 * np.pi (stop_frequency-start_frequency) * 1E6 / (2. * duration * 1.E-6 * (t * 1E-6) ** 2)
            self.pulse_scheme.append([t, transition, freq, intensity, instant_phase])
