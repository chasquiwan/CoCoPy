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
import string as strg
from StringIO import StringIO
import re
import os
import scipy.constants as constants
import analysis.util as util
import copy
import sys

'''
Todo:
    - 1 Write comments
    - 2
    - 3
'''

class spfitPrediction:

    def __init__(self, fname = 'test', name = 'Test'):
        '''
        This class tries to cover everything, which is related to a SPFIT/SPCAT
        program suite. This means it comes with file parsers and writers. It also
        includes some specialized plotting function and autoassignment methods.

        Starts an empty spfitPrediction instance. Filename and name can be
        specified.
        '''

        self.FILE_EXT = {
            'cat': '.cat', 'par': '.par', 'lin': '.lin', 'egy': '.egy',
            'var': '.var', 'fit': '.fit', 'out': '.out', 'strf': '.str',
            'intf': '.int', 'bak': '.bak'
            }

        self.VAR_PARAMS_LABEL1 = ['title', 'npar', 'nline', 'nitr', 'nxpar',
            'thresh', 'errtst', 'frac', 'cal']

        self.VAR_PARAMS_LABEL2 = ['char', 'spind', 'nvib', 'knmin', 'knmax',
            'ixx', 'iax', 'wtpl', 'wtmn', 'vsym', 'ewt', 'diag', 'xopt']

        self.VAR_PARAMS_FORMAT2 = ['', ':d', ':d', ':d', ':d', '', ':d', ':.2f',
            ':.2f', ':d', '', ':d', '']

        self.INT_PARAMS_LABEL = ['title', 'flags', 'tag', 'qrot', 'fbgn',
            'fend', 'str0', 'str1', 'fqlim', 'temp', 'maxv']

        self.CAT_DT = np.dtype([('freq', 'f8'), ('err', 'f8'), ('lgint', 'f8'),
            ('dr', 'i2'), ('elo', 'f8'), ('gup', 'i2'), ('tag', 'i2'),
            ('qnfmt', 'i2'), ('qn', 'S24')])

        self.CAT_DEL = [13, 8, 8, 2, 10, 3, 7, 4, 24]

        self.LIN_DT = np.dtype([('qn', 'S36'), ('rest', 'S100')])

        self.LIN_DEL = [36, 100]

        self.STR_DT = np.dtype([('freq', 'f8'), ('dipole', 'f8'),
            ('qnfmt', 'i2'), ('qn', 'S24'), ('items', 'i2')])

        self.STR_DEL = [15, 15, 5, 24, 10]

        self.EGY_DT = np.dtype([('iblk', 'i2'), ('indx', 'i2'), ('egy', 'f8'),
            ('pmix', 'f8'), ('err', 'f8'), ('udn', 'S4'), ('qn', 'S18')])

        self.EGY_DEL = [6, 6, 18, 18, 12, 4, 18]


        self.var_params = {
            'title': '', 'npar': 3, 'nline': 1000, 'nitr': 99, 'nxpar': 0,
            'thresh': 1.0E-3, 'errtst': 1.0E3, 'frac': 1.0E0, 'cal': 1.0,
            'char': 'a', 'spind': 1, 'nvib': 1, 'knmin': 0, 'knmax': 15,
            'ixx': 0, 'iax': 0, 'wtpl': 0.0, 'wtmn': 0.0, 'vsym': 0, 'ewt': 0,
            'diag': 0, 'xopt': 0
            }

        self.var_fitpar = {
            'idpar': list(), 'par': np.array([]), 'erpar': np.array([]),
            'label': list(), 'nfitpar': 0
            }

        self.par_params = {
            'title': '', 'npar': 3, 'nline': 1000, 'nitr': 99, 'nxpar': 0,
            'thresh': 1.0E-3, 'errtst': 1.0E3, 'frac': 1.0E0, 'cal': 1.0,
            'char': 'a', 'spind': 1, 'nvib': 1, 'knmin': 0, 'knmax': 15,
            'ixx': '0', 'iax': 0, 'wtpl': 0.0, 'wtmn': 0.0, 'vsym': '0',
            'ewt': '0', 'diag': 0, 'xopt': '0'
            }

        self.par_fitpar = {
            'idpar': list(), 'par': np.array([]), 'erpar': np.array([]),
            'label': list(), 'nfitpar': 0
            }

        self.int_params = {
            'title': '', 'flags': '0111', 'tag': 666, 'qrot': 0.0, 'fbgn': 0,
            'fend': 99, 'str0': -10.0, 'str1': -10.0, 'fqlim': 8.5, 'temp': .3,
            'maxv': 999
            }

        self.int_dipole = {
            'idip': list(), 'dipole': np.array([]), 'comment': list(),
            'ndip': 0
            }

        self.cat_content = {
            'freq': np.array([]), 'err': np.array([]), 'lgint': np.array([]),
            'dr': np.array([]), 'elo': np.array([]), 'gup': np.array([]),
            'tag': np.array([]), 'qnfmt': np.array([]), 'qn': list()
            }

        self.lin_content = {
            'qn': list(), 'freq': np.array([]), 'err': np.array([]),
            'wt': np.array([]), 'intens': np.array([])
            }

        self.str_content = {
            'freq': np.array([]), 'dipole': np.array([]), 'qnfmt': np.array([]),
            'qn': list(), 'items': np.array([])
            }

        self.egy_content = {
            'iblk': np.array([]), 'indx': np.array([]), 'egy': np.array([]),
            'pmix': np.array([]), 'err': np.array([]), 'udn': np.array([]),
            'qn': list()
            }

        self.fit_content = {
            'err_mw_rms': 1.0, 'err_mw_avg': 1.0, 'err_old_rms': 1.0,
            'err_new_rms': 1.0, 'iteration': 1, 'reject': False,
            'N_reject_lines': 0, 'N_fitted_lines': 0, 'N_initial_lines': 0
            }

        self.mol_params = {
            'nqn': 3, 'nline': 100, 'nvib': 1, 'nvarpar': 5, 'linerr': 0.02,
            'flag_linerr': 1, 'flag_linwt': 0, 'nintpar': 9, 'var_ex': 0,
            'par_ex': False, 'int_ex': False, 'lin_ex': False, 'cat_ex': False,
            'egy_ex': False, 'fit_ex': False, 'out_ex': False, 'bak_ex': 0,
            'str_ex': False, 'name': 'test', 'fname': 'test'
            }

        self.meas_spec = {
            'spec': np.array([[]]), 'linelist': np.array([[]]),
            'qn': np.array([[]]), 'scale': 1.0
            }

        self.pred_spec = {
            'spec': np.array([[]]), 'linelist': np.array([[]]),
            'qn': np.array([[]]), 'scale': 1.0
            }

        self.energy_states = {
            'qn': np.array([[]]), 'egy': np.array([]),
            'population': np.array([])
            }

        self.exe = {'spcat': 'spcat', 'spfit': 'spfit'}

        self.mol_params['name'] = name
        self.mol_params['fname'] = fname

        if sys.platform == 'win32':
            self.exe['spcat'] = 'spcat.exe'
            self.exe['spfit'] = 'spfit.exe'


    def write_pvar(self, fname='', trg='var'):
        '''
        Creates and writes the var or a par file and uses the values stored in
        var_params, var_fitpar. Not for normal use. Use write_par(), write_var()
        instead.

        Parameters
        ----------
        fname (str): filename of the file to write without the .xxx extension
            ('')
        trg (str): defines if a par or a var file should be created ('var')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + self.FILE_EXT[trg], 'w')
        except IOError:
            self.mol_params[trg+'_ex'] = False
            print 'Cannot open: ', self.mol_params['fname'] + self.FILE_EXT[trg]
        else:
            if trg == 'var':
                params = self.var_params
                fitpar = self.var_fitpar
            else:
                params = self.par_params
                fitpar = self.par_fitpar

            self.mol_params[trg+'_ex'] = True
            filestr = self.mol_params['name'] + '\n'
            filestr += '{npar:d}\t{nline:d}\t{nitr:d}\t{nxpar:d}\t{thresh:.2e}\t{errtst:.2e}\t{frac:.2e}\t{cal:.2f}\n'.format(**params)
            i = 0
            for key in self.VAR_PARAMS_LABEL2:
                if i >= self.mol_params['nvarpar']:
                    filestr += ''
                else:
                    filestr += ('{' + key + self.VAR_PARAMS_FORMAT2[i] + '}\t').format(**params)
                i += 1
            filestr += '\n'

            i = 0
            while i < fitpar['nfitpar']:

                idpar = strg.rjust(fitpar['idpar'][i], 15)
                par = strg.rjust('{0:.8e}'.format(fitpar['par'][i]), 20)
                erpar = strg.rjust('{0:.8e}'.format(fitpar['erpar'][i]), 20)
                label = '   /' + strg.lstrip(fitpar['label'][i])
                filestr += idpar + par + erpar + label + '\n'

                i += 1
            f.write(filestr)
            f.close()


    def write_var(self, fname=''):
        '''
        Creates a var file and uses the values stored in var_params, var_fitpar.

        Parameters
        ----------
        fname (str): filename of the file to write without the .var extension
            ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        self.write_pvar(fname, trg='var')


    def write_par(self, fname=''):
        '''
        Creates a var file and uses the values stored in par_params, par_fitpar.

        Parameters
        ----------
        fname (str): filename of the file to write without the .par extension
            ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        self.write_pvar(fname, trg='par')


    def write_lin(self, fname='', sign_fig=4):
        '''
        Creates a lin file and uses the values stored in lin_content.

        Parameters
        ----------
        fname (str): filename of the file to write without the .lin extension
            ('')
        sign_fig (int): number of significant figures for the lin file

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + self.FILE_EXT['lin'], 'w')
        except IOError:
            self.mol_params['lin_ex'] = False
            print 'Cannot open: ', self.mol_params['fname'] + self.FILE_EXT['lin']
        else:
            self.mol_params['lin_ex'] = True
            i = 0
            filestr = ''
            for x in self.lin_content['qn']:
                filestr += '{0}'.format(x)
                h = ' {0:.'+str(sign_fig)+'f}'
                filestr += h.format(np.round(self.lin_content['freq'][i], sign_fig))
                if len(self.lin_content['err']) > 0:
                    filestr += strg.ljust('', 5) + '{0:.3f}'.format(self.lin_content['err'][i])
                else:
                    filestr += strg.ljust('', 5) + '{0:.3f}'.format(self.mol_params['linerr'])

                if len(self.lin_content['wt']) > 0:
                    filestr += strg.ljust('', 5) + '{0:.2f}'.format(self.lin_content['wt'][i])
                else:
                    filestr += strg.ljust('', 5) + '1.0'
                if len(self.lin_content['intens']) > 0:
                    filestr += '    /' + str(self.lin_content['intens'][i]) + '\n'
                else:
                    filestr += '    /0.0\n'
                i += 1
            f.write(filestr)
            f.close()


    def write_int(self, fname=''):
        '''
        Creates a int file and uses the values stored in int_params and
        int_dipole.

        Parameters
        ----------
        fname (str): filename of the file to write without the .int extension
            ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + self.FILE_EXT['intf'], 'w')
        except IOError:
            self.mol_params['int_ex'] = False
            print 'Cannot open: ', self.mol_params['fname'] + self.FILE_EXT['intf']
        else:
            self.mol_params['int_ex'] = True
            filestr = self.mol_params['name'] + "\n"
            filestr += "{flags}\t{tag:d}\t{qrot:.2f}\t{fbgn:d}\t{fend:d}\t{str0:.2f}\t{str1:.2f}\t{fqlim:.2f}\t{temp:.2f}\t{maxv:d}\n".format(**self.int_params)
            i = 0
            while i < self.int_dipole['ndip']:
                idip = strg.rjust(self.int_dipole['idip'][i], 10)
                dipole = strg.rjust('{0:.4}'.format(self.int_dipole['dipole'][i]), 10)
                comment = '  /' + strg.lstrip(self.int_dipole['comment'][i])
                filestr += idip + dipole + comment + '\n'
                i += 1

            f.write(filestr)
            f.close()


    def read_pvar(self, fname='', src='var'):
        '''
        Reads a var or a par file and stores the contens in the variables
        var_params, var_fitpar.
        Not for normal use. Use read_par(), read_var() instead.

        Parameters
        ----------
        fname (str): filename of the file to read from without the .xxx
            extension ('')
        src (str): defines if the par or the var file should be read ('var')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + self.FILE_EXT[src], 'r')
        except IOError:
            self.mol_params[src+'_ex'] = False
            print 'Cannot open: ', self.mol_params['fname'] + self.FILE_EXT[src]
        else:
            if src == 'var':
                params = self.var_params
                fitpar = self.var_fitpar
            else:
                params = self.par_params
                fitpar = self.par_fitpar

            flag = 0

            self.mol_params[src+'_ex'] = True
            i = 0
            for line in f:
                if re.match('^\s*\n$', line) != None:
                    i += 1
                    continue
                if i == 0:
                    params['title'] = strg.strip(line[0:50])
                if i == 1:
                    if strg.find(line, ',') != -1:
                        h = np.genfromtxt(StringIO(strg.strip(line)), dtype=str, delimiter=',')
                    else:
                        h = np.genfromtxt(StringIO(strg.strip(line)), dtype=str)

                    j = 1
                    for x in h:
                        if x == '':
                            j += 1
                            continue

                        if self.VAR_PARAMS_LABEL1[j] in ('npar', 'nline', 'nitr', 'nxpar'):
                            params[self.VAR_PARAMS_LABEL1[j]] = int(x)
                        else:
                            params[self.VAR_PARAMS_LABEL1[j]] = float(x)
                        j += 1

                if i == 2:
                    if strg.find(line, ',') != -1:
                        h = np.genfromtxt(StringIO(strg.strip(line)), dtype=str, delimiter=',')
                    else:
                        h = np.genfromtxt(StringIO(strg.strip(line)), dtype=str)

                    j = 0
                    for x in h:
                        if x == '':
                            j += 1
                            continue
                        if self.VAR_PARAMS_LABEL2[j] in ('char', 'ixx', 'vsym', 'ewt', 'xopt'):
                            params[self.VAR_PARAMS_LABEL2[j]] = x
                        elif self.VAR_PARAMS_LABEL2[j] in ('wtpl', 'wtpn'):
                            params[self.VAR_PARAMS_LABEL2[j]] = float(x)
                        else:
                            params[self.VAR_PARAMS_LABEL2[j]] = int(x)
                        j += 1

                if i >= 3:
                    if flag == 0:
                        fitpar['idpar']=list(); fitpar['par']=np.array([]); fitpar['erpar']=np.array([]); fitpar['label']=list(); fitpar['nfitpar']=0
                        flag = 1
                    if re.match('^\s\s.*$', line) == None:
                        i += 1
                        continue
                    fitpar['nfitpar'] += 1
                    h = np.genfromtxt(StringIO(strg.strip(line)), dtype=str)
                    fitpar['idpar'].append(h[0])
                    fitpar['par'] = np.append(fitpar['par'], float(h[1]))
                    if len(h) > 2:
                        fitpar['erpar'] = np.append(fitpar['erpar'], float(h[2]))
                    if len(h) > 3:
                        g = ''
                        for x in np.arange(3, len(h)):
                            g += " " + h[x]

                        fitpar['label'].append(g[2:])
                i += 1
            f.close()


    def read_var(self, fname=''):
        '''
        Reads the var file and stores the contents in the variables var_params,
        var_fitpar.

        Parameters
        ----------
        fname (str): filename of the file to read from without the .var
            extension ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        self.read_pvar(fname, 'var')

    def read_par(self, fname=''):
        '''
        Reads the par file and stores the contents in the variables par_params, par_fitpar.

        Parameters
        ----------
        fname (str): filename of the file to read from without the .par extension ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        self.read_pvar(fname, 'par')


    def read_int(self, fname=''):
        '''
        Reads the int file and stores the contents in the variables int_dipole,
        int_params.

        Parameters
        ----------
        fname (str): filename of the file to read from without the .int
            extension ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + self.FILE_EXT['intf'], 'r')
        except IOError:
            self.mol_params['int_ex'] = False
            print 'Cannot open: ', self.mol_params['fname'] + self.FILE_EXT['intf']
        else:
            self.mol_params['int_ex'] = True
            i = 0
            for line in f:
                if re.match('^\s*\n$', line) != None:
                    i += 1
                    continue
                if i == 0:
                    self.int_params['title'] = strg.strip(line)
                if i == 1:
                    if strg.find(line, ',') != -1:
                        h = np.genfromtxt(StringIO(strg.strip(line)), dtype=str, delimiter=',')
                    else:
                        h = np.genfromtxt(StringIO(strg.strip(line)), dtype=str)
                    self.mol_params['nintpar'] = len(h)+1

                    j = 1
                    for x in h:
                        if x == '':
                            j += 1
                            continue
                        if self.INT_PARAMS_LABEL[j] == 'flags':
                            self.int_params[self.INT_PARAMS_LABEL[j]] = x
                        elif self.INT_PARAMS_LABEL[j] in ('fbgn', 'fend', 'tag', 'maxv'):
                            self.int_params[self.INT_PARAMS_LABEL[j]] = int(x)
                        else:
                            self.int_params[self.INT_PARAMS_LABEL[j]] = float(x)
                        j += 1

                if i >= 2:
                    self.int_dipole['ndip'] += 1
                    h = np.genfromtxt(StringIO(strg.strip(line)), dtype=str)
                    self.int_dipole['idip'].append(h[0])
                    self.int_dipole['dipole'] = np.append(self.int_dipole['dipole'], float(h[1]))
                    if len(h) > 2:
                        g = ''
                        for x in np.arange(2, len(h)):
                            g += ' ' + h[x]
                        self.int_dipole['comment'].append(g[1:])
                i += 1
            f.close()

    def read_lin(self, fname=''):
        '''
        Reads the lin file and stores the contents in the variable lin_content.

        Parameters
        ----------
        fname (str): filename of the file to read from without the .lin
            extension ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + self.FILE_EXT['lin'], 'r')
        except IOError:
            self.mol_params['lin_ex'] = False
            print 'Cannot open: ', self.mol_params['fname'] + self.FILE_EXT['lin']
        else:
            self.mol_params['lin_ex'] = True
            data = np.genfromtxt(f, dtype=self.LIN_DT, delimiter=self.LIN_DEL)
            freq = np.array([])
            err = np.array([])
            wt = np.array([])
            intens = np.array([])

            for x in data['rest']:
                g = strg.split(x, '/')
                h = np.genfromtxt(StringIO(strg.strip(g[0])), dtype=float)
                freq = np.append(freq, h[0])


                if len(h) > 1:
                    err = np.append(err, h[1])
                else:
                    err = np.append(err, 0.0)

                if len(h) > 2:
                    wt = np.append(wt, h[2])
                else:
                    wt = np.append(wt, 1.0)

                if len(g) > 1:
                    intens = np.append(intens, np.genfromtxt(StringIO(strg.strip(g[1])), dtype=float))
                else:
                    intens = np.append(intens, 1.)


            self.lin_content['qn']=data['qn']
            self.lin_content['freq']=freq
            self.lin_content['err']=err
            self.lin_content['wt']=wt
            self.lin_content['intens']=intens
            f.close()
            i = 0
            check = '1'
            h = list()
            while check != '':
                check = strg.rstrip(data['qn'][0])[i:i+3]
                if check != '':
                    h.append(check)
                i += 3
            self.mol_params['nqn'] = int(len(h) / 2)

            self.mol_params['nline'] = len(data)

    def read_cat(self, fname=''):
        '''
        Reads the cat file and stores the contents in the variable cat_content.

        Parameters
        ----------
        fname (str): filename of the file to read from without the .cat
            extension ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + self.FILE_EXT['cat'], 'r')
        except IOError:
            self.mol_params['cat_ex'] = False
            print 'Cannot open: ', self.mol_params['fname'] + self.FILE_EXT['cat']
        else:
            self.mol_params['cat_ex'] = True
            data = np.genfromtxt(f, dtype=self.CAT_DT, delimiter=self.CAT_DEL)

            for key in self.cat_content.iterkeys():
                self.cat_content[key]=data[key]

            f.close()

            i = 0
            check = '1'
            h = list()
            while check != '':
                check = strg.rstrip(data['qn'][0])[i:i+2]
                if check != '' and check != '  ':
                    h.append(check)
                i += 2
            self.mol_params['nqn'] = int(len(h) / 2)

    def read_egy(self, fname=''):
        '''
        Reads the egy file and stores the contents in the variable egy_content.

        Parameters
        ----------
        fname (str): filename of the file to read from without the .egy
            extension ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + self.FILE_EXT['egy'], 'r')
        except IOError:
            self.mol_params['egy_ex'] = False
            print 'Cannot open: ', self.mol_params['fname'] + self.FILE_EXT['egy']
        else:
            self.mol_params['egy_ex'] = True
            data = np.genfromtxt(f, dtype=self.EGY_DT, delimiter=self.EGY_DEL)

            for key in self.egy_content.iterkeys():
                self.egy_content[key]=data[key]

            f.close()

    def read_str(self, fname=''):
        '''
        Reads the str file and stores the contents in the variable str_content.

        Parameters
        ----------
        fname (str): filename of the file to read from without the .str
            extension ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + self.FILE_EXT['strf'], 'r')
        except IOError:
            self.mol_params['str_ex'] = False
            print 'Cannot open: ', self.mol_params['fname'] + self.FILE_EXT['strf']
        else:
            self.mol_params['str_ex'] = True
            data = np.genfromtxt(f, dtype=self.STR_DT, delimiter=self.STR_DEL)

            for key in self.str_content.iterkeys():
                self.str_content[key]=data[key]

            f.close()

    def read_out(self, fname=''):
        '''
        Not yet implemented.

        Parameters
        ----------
        fname (str): filename of the file to read from without the .out extension ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        self.written = 1

    def read_fit(self, fname=''):
        '''
        Reads the fit file and stores the contents in the variable fit_content.

        Parameters
        ----------
        fname (str): filename of the file to read from without the .fit
            extension ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''

        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + self.FILE_EXT['fit'], 'r')
        except IOError:
            self.mol_params['fit_ex'] = False
            print 'Cannot open: ', self.mol_params['fname'] + self.FILE_EXT['fit']
        else:
            self.mol_params['fit_ex'] = True

            for line in f:
                if strg.find(line, 'Lines rejected from fit') != -1:
                    self.fit_content['reject'] = True
                    self.fit_content['N_reject_lines'] = int(line.split('Lines')[0])
                if re.search('^\s*\d+:', line) is not None:
                    self.fit_content['N_initial_lines'] = int(re.search('^\s*\d+:', line).group()[:-1])
                if strg.find(line, 'MICROWAVE AVG') != -1:
                    self.fit_content['err_mw_avg'] = float(strg.strip(line[20:32]))
                if strg.find(line, 'MICROWAVE RMS') != -1:
                    self.fit_content['err_mw_rms'] = float(strg.strip(line[20:32]))
                if strg.find(line, 'END OF ITERATION') != -1:
                    self.fit_content['iteration'] = int(strg.strip(line[17:20]))
                    self.fit_content['err_old_rms'] = float(strg.strip(line[42:60]))
                    self.fit_content['err_new_rms'] = float(strg.strip(line[61:75]))
            f.close()
            self.fit_content['N_fitted_lines'] = self.fit_content['N_initial_lines'] - self.fit_content['N_reject_lines']


    def cal_part_function(self, temp = 0.0):
        '''
        Calculates the partition sum Q using the energies from the egy file.

        Parameters
        ----------
        temp (float): temperature in Kelvin

        Returns
        -------
        Returns the partition sum Q.

        Notes
        -----
        none
        '''

        if temp == 0.0:
            temp = self.int_params['temp']
        if len(self.egy_content['egy']) != 0:
            Q = 0.
            i = 0
            egy_lev = util.convert_energy(self.egy_content['egy'], 'wavenumber', 'joule')

            for x in egy_lev:
                Q += (2*int(self.egy_content['qn'][i][:3]) + 1) * np.exp(-x/(constants.k*temp))
                i += 1
            return Q

    def make_linel_lin(self):
        '''
        Creates a linelist of the assigned lines from the lin file. The linelist
        is stored in the variable self.meas_spec['linelist'].

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

        if len(self.lin_content['freq']) == 0:
            print "No data available!"
        else:
            self.meas_spec['qn'] = standard_qn(self.lin_content['qn'], src='lin')
            self.meas_spec['linelist'] = np.hstack((self.lin_content['freq'].reshape((len(self.lin_content['freq'])), 1), self.lin_content['intens'].reshape((len(self.lin_content['intens']), 1))))

    def make_linel_cat(self):
        '''
        Creates a linelist of the predicted lines from the cat file. The linelist
        is stored in the variable self.pred_spec['linelist'].

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

        if len(self.cat_content['freq']) == 0:
            print "No data available!"
            print len(self.cat_content['freq'])
        else:
            self.pred_spec['qn'] = standard_qn(self.cat_content['qn'], src='cat')
            freq = self.cat_content['freq']
            intens = (10 ** (self.cat_content['lgint'])) * self.pred_spec['scale']
            self.pred_spec['linelist'] = np.hstack((freq.reshape((len(freq), 1)), intens.reshape((len(intens), 1))))

    def make_spec_cat(self):
        '''
        Creates a plottable spectrum of the predicted lines from the cat file.
        The linelist is saved in the variable self.pred_spec['spec'].

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

        self.make_linel_cat()
        if len(self.pred_spec['linelist']) != 0:
            self.pred_spec['spec'] = lines_to_spec(self.pred_spec['linelist'])

    def make_spec_lin(self):
        '''
        Creates a plottable spectrum of the assigned lines from the lin file.
        The linelist is saved in the variable self.meas_spec['spec'].

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

        self.make_linel_lin()
        if len(self.meas_spec['linelist']) != 0:
            self.meas_spec['spec'] = lines_to_spec(self.meas_spec['linelist'])

    def get_transitions(self, trans=([]), pos=([]), src='cat'):
        '''
        Extracts and returns the index of the cat linelist matching a specific
        quantum number combination.

        Parameters
        ----------
        trans (list, int): quantum numbers ([])
        pos (list, int): position of the quantum number ([])
        scr (str): only the cat source is implemented yet, so you should not
            change anythin here ('cat')

        Returns
        -------

        Notes
        -----
        none
        '''

        if len(trans) > 0:
            if src=='cat':
                self.make_linel_cat()
                std_qn = standard_qn(self.cat_content['qn'], 'cat')
                index = list()
                i = 0
                for x in trans:
                    j = 0
                    for y in std_qn:
                        if np.array_equal(y[pos[i]], x) == True:
                            index.append(j)
                        j+=1
                    i+=1
                index = set(index)
                index = sorted(index)

                return index


    def cat_run(self, fname=''):
        '''
        Runs the spcat program in the working directory. The excecutable has to
        be in this directory.

        Parameters
        ----------
        fname (str): filename for the set of files to use, without the extension
            ('')

        Returns
        -------
        check (boolean): TODO

        Notes
        -----
        none
        '''

        check = -1
        if fname == '':
            fname = self.mol_params['fname']
        else:
            self.mol_params['fname'] = fname
        check = os.system(self.exe['spcat'] + ' ' + fname + self.FILE_EXT['var'] + ' > ' + fname + '.dmp')

        return check

    def fit_run(self, fname=''):
        '''
        Runs the spcat program in the working directory. The excecutable has to
        be in this directory.

        Parameters
        ----------
        fname (str): filename for the set of files to use, without the extension
            ('')

        Returns
        -------
        check (boolean): TODO

        Notes
        -----
        none
        '''

        check = -1
        if fname == '':
            fname = self.mol_params['fname']
        else:
            self.mol_params['fname'] = fname
        check = os.system(self.exe['spfit'] + ' ' + fname + self.FILE_EXT['par'] + ' > ' + fname + '.dmp')

        return check

    def fit_min(self):
        '''
        Writes the par, int and lin file and runs the spfit program afterwards
        in the working directory. Finally it reads in the fit and the par file.
        The excecutable has to be in this directory or defined globally.

        Parameters
        ----------
        none

        Returns
        -------
        check (boolean): TODO

        Notes
        -----
        none
        '''

        self.write_par()
        self.write_int()
        self.write_lin()
        check = self.fit_run()
        self.read_fit()
        self.read_par()

        return check


    def cat_min(self):
        '''
        Writes the var and int file and runs the spcat program afterwards
        in the working directory. Finally it reads in the cat file.
        The excecutable has to be in this directory or defined globally.

        Parameters
        ----------
        none

        Returns
        -------
        check (boolean): TODO

        Notes
        -----
        none
        '''

        self.int_params['flags'] = '0000';
        self.write_var()
        self.write_int()
        check = self.cat_run()
        self.read_cat()

        return check


    def cat_max(self):
        '''
        Writes the var and int file and runs the spcat program afterwards
        in the working directory. Finally it reads in the cat, egy and str
        files. The excecutable has to be in this directory or defined globally.

        Parameters
        ----------
        none

        Returns
        -------
        check (boolean): TODO

        Notes
        -----
        none
        '''

        self.write_var()
        self.write_int()
        check = self.cat_run()
        self.read_cat()
        self.read_egy()
        self.read_str()

        return check


    def refine_lines(self, spectrum, thresh=0.050, weight=False):
        '''
        Compares the assigned lines with a spectrum and refines the lines by
        lookingfor the absolute maximum intensity or by also weighting the
        neighboring points of the maximum.

        Parameters
        ----------
        spectrum (float): 2 column array, frequencies in the first column and
            intensities in the second column

        Returns
        -------
        freq (float): list of the refined frequencies
        intens ( float): list of the intensities

        Notes
        -----
        none
        '''

        return refine_lines(self.lin_content['freq'], spectrum, thresh, weight)

    def autoassign(self, peaks, freq_thresh=10.0, int_thresh=10.0, finesse=10.0, lin=True):
        '''
        Compares the assigned lines with a spectrum and refines the lines by looking
        for the absolute maximum intensity or by also weighting the neighboring
        points of the maximum.

        Parameters
        ----------
        spectrum (np.array (2D), float): 2 column array with frequencies in the
            first column and intensities in the second column

        Returns
        -------
        freq (list, float): list of the refined frequencies
        intens (list, float): list of the intensities

        Notes
        -----
        none
        '''

        def remove_dublicates(seq1, seq2, seq3):
            seq1_new = []
            seq2_new = []
            seq3_new = []
            i = 0
            for e in seq1:
                if e not in seq1_new:
                    seq1_new.append(e)
                    seq2_new.append(seq2[i])
                    seq3_new.append(seq3[i])

                i += 1
            return seq1_new, seq2_new, seq3_new

        i = 0
        self.make_linel_cat()
        freq = np.array([])
        intens = np.array([])
        qn = list()

        if len(peaks) == 0 or len(self.pred_spec['qn']) == 0 or len(self.pred_spec['linelist']):
            for x in self.pred_spec['linelist'][:]:   # Go through the list of predictions
                succ = 0
                h_thresh = freq_thresh
                while succ == 0:
                    a = len(np.where(np.abs( peaks[:,0] - x[0]) < h_thresh)[0]) # Search for a peak nearby
                    if a == 1:  # Only one peak within the range freq_threshold_help
                        y = peaks[np.where(np.abs(peaks[:,0] - x[0]) < h_thresh)][0][0]
                        y_int = peaks[np.where(np.abs(peaks[:,0] - x[0]) < h_thresh)][0][1]
                        if (y_int/x[1] if y_int > x[1] else x[1]/y_int)  < int_thresh:
                            freq = np.append(freq, y)
                            intens = np.append(intens, y_int)
                            qn.append(self.pred_spec['qn'][i])
                        succ = 1
                    elif a == 0:    # No peak found -> go to next peak in prediction list
                        succ = 1
                    else:   # More than one peak within the range freq_threshold_help -> make freq_threshold_help smaller
                        h_thresh -= float(h_thresh)/float(finesse)
                i+=1

            freq, intens, qn = remove_dublicates(freq, intens, qn)

            if lin == True:
                self.lin_content['freq'] = freq
                self.lin_content['qn'] = convert_qn_lin(qn)
                self.lin_content['err'] = np.array([])
                self.lin_content['wt'] = np.array([])
                self.lin_content['intens'] = np.array([])
            return freq, intens, qn
        else:
            print 'No prediction data OR no peaks available'


    def plot(self, x_lim=[0.0, 0.0], img=1, scale=1.0, offset=0.0, shift=0.0,
    src='cat', tr_type='all', branch='pq', label='', inline = False,
    rasterized = False):
        '''
        Plots a spectrum of the predicted (cat) or assigned lines (lin). The
        frequency range can be specified. The intensities can be scaled by a
        scaling factor and also shifted by a specific amount. Beside that you
        can choose from different types of lines (a-type, b-type, c-type) and
        include or exclude q-type transitions.
        TODO: fig, ax

        Parameters
        ----------
        x_lim (list, float): two item list defining the start and the end
            frequency of the plot range, if both numbers are zero
        img (int): number of the image the plot should appear (1)
        scale (float): the intensities are scaled by this value (1.0)
        shift (float): the x-axis is shifted by this value (0.0)
        scr (str): source of the data: cat or lin (cat)
        tr_type (str): transition type to plot: all, a, b, c, ab, abc, etc.
        branch (str): plots only a certain branch: p or q (pq)
        label (str): label of the data in the plot ('')

        Returns
        -------
        none

        Notes
        -----
        none
        '''
        #TODO: FIX return ax and fig

        import matplotlib.pyplot as plt

        if src == 'cat':
            if len(self.pred_spec['spec']) > 0:
                plt.ion()
                fig, ax = plt.subplots()
                flag = 0

                if tr_type == 'all':
                    w = self.pred_spec['linelist']
                    v = self.pred_spec['qn']
                if 'a' in tr_type and tr_type != 'all':
                    v, w = get_abc_transitions(self.pred_spec['qn'], self.pred_spec['linelist'], 'a')
                    flag = 1
                if 'b' in tr_type:
                    if flag == 1:
                        v_h, w_h = get_abc_transitions(self.pred_spec['qn'], self.pred_spec['linelist'], 'b')
                        v = np.vstack((v, v_h))
                        w = np.vstack((w, w_h))
                    else:
                        v, w = get_abc_transitions(self.pred_spec['qn'], self.pred_spec['linelist'], 'b')
                        flag = 1
                if 'c' in tr_type:
                    if flag == 1:
                        v_h, w_h = get_abc_transitions(self.pred_spec['qn'], self.pred_spec['linelist'], 'c')
                        v = np.vstack((v, v_h))
                        w = np.vstack((w, w_h))
                    else:
                        v, w = get_abc_transitions(self.pred_spec['qn'], self.pred_spec['linelist'], 'c')
                        flag = 1

                if branch == 'p':
                    v, w = q_branch(v, w, 'p')
                elif branch == 'q':
                    v, w = q_branch(v, w, 'q')

                w = lines_to_spec(w)
                if x_lim[0] != x_lim[1]:
                    ax.xlim(x_lim)
                ax.plot(w[:,0]+shift, w[:,1]*scale+offset, label=label, rasterized=rasterized)
                ax.xlabel('Frequency [MHz]')
                ax.ylabel('Intensity [arb. U.]')
                ax.legend()
                if inline == False:
                    fig.show()

        if src == 'lin':
            if len(self.meas_spec['spec']) > 0:
                if tr_type == 'all':
                    plt.ion()
                    fig, ax = plt.subplots()
                    ax.plot(self.meas_spec['spec'][:,0]+shift, self.meas_spec['spec'][:,1]*scale, label=label)
                    if x_lim[0] != x_lim[1]:
                        ax.xlim(x_lim)
                    if inline == False:
                        fig.show()

class spfitAsym(spfitPrediction):
    '''
    The Spectrum class tries to cover everything, which is related to a SPFIT/SPCAT
    program suite. This means it comes with file parsers and writers. It also includes
    some specialized plotting function and autoassignment methods.
    '''
    def __init__(self, fname = 'test', name='Test', A=2000., B=700., C=600., dipA = 1. , dipB = 1. , dipC = 1., temp=1.5):
        '''
        Starts an empty spfitPrediction instance. Filename and Name can be specified.
        '''
        spfitPrediction.__init__(self)

        self.mol_params['name'] = name
        self.mol_params['fname'] = fname

        self.var_fitpar['idpar'] = ['10000', '20000', '30000']
        self.var_fitpar['par'] = np.array([A, B, C])
        self.var_fitpar['erpar'] = np.array([A/10., B/10., C/10.])
        self.var_fitpar['label'] = ['A', 'B', 'C']
        self.var_fitpar['nfitpar'] = 3
        self.par_fitpar = copy.copy(self.var_fitpar)

        self.lin_content['qn'] = convert_qn_lin(np.array([[1,1,1,0,0,0]]))
        self.lin_content['freq'] = np.array([1000.])

        self.int_dipole['idip'] = ['1', '2', '3']
        self.int_dipole['dipole'] = np.array([dipA, dipB, dipC])
        self.int_dipole['comment'] = ['a', 'b', 'c']

        self.int_params['temp'] = temp
        self.int_params['qrot'] = 5.3311E6*np.sqrt(self.int_params['temp']**3/(A*B*C))
        self.int_dipole['ndip'] = 3

################################################################################
def reshape_linelist(lines, intensity=1.0):
    '''
    Reshapes a linelist: 1D -> 2D

    Parameters
    ----------
    lines (list, np.array, float): 1 or 2 column array with frequencies in the
    first column and intensities in the second column or a list. If only one
    dimension is specified (only frequencies) then the intensities are set to 1.0
    intensity (float): standard intensities if only the frequencies are specified

    Returns
    -------
    linelist (2-col. array, float): the generated linelist file

    Notes
    -----
    none
    '''

    h = np.array(lines, copy=True)
    if len(h.shape) == 1:
        g = np.ones(len(h))
        g = g.reshape((len(g), 1))
        h = np.hstack((h.reshape((len(h), 1)), g))

        return h

    if len(h.shape) == 2:

        return h

def lines_to_spec(lines, eps = .0001):
    '''
    Converts a linelist to a spectrum file

    Parameters
    ----------
    lines (list, np.array, float): 1 or 2 column array with frequencies in the
    first column and intensities in the second column or a list. If only one
    dimension is specified (only frequencies) then the intensities are set to 1.0
    eps (float): linewidth in MHz (0.0001)

    Returns
    -------
    spectrum (2-col. array, float): the generated spectrum file

    Notes
    -----
    none
    '''

    h = reshape_linelist(lines)
    i = 0
    h = np.split(h, len(h))

    while i < len(h):   # Puts a zero next to each line
        h1 = [[h[i][0][0] - eps, 0.]]
        h2 = [[h[i][0][0] + eps, 0.]]
        #h[i][0][1] = intensity
        h[i] = np.vstack((h1, h[i], h2))
        i+=1

    return np.vstack(h[:])

def convert_qn_cat_lin(qn_in):
    '''
    Converts the quantum numbers from the cat to the lin format

    Parameters
    ----------
    qn_in (str): quantum numbers in the cat format

    Returns
    -------
    qn (str): quantum numbers in the lin format

    Notes
    -----
    none
    '''

    if len(qn_in) == 0:
        print "No data available!"
    else:
        return convert_qn_lin(standard_qn(qn_in))


def convert_qn_lin(qn_in):
    '''
    Converts standard quantum numbers to the lin format

    Parameters
    ----------
    qn_in: (array, int): 2D array of the quantum numbers in the standard format

    Returns
    -------
    qn_in (array, str): 2D array of the quantum numbers in the lin format

    Notes
    -----
    none
    '''

    if len(qn_in) == 0:
        print "No data available!"
    else:
        qn_lin = list()
        for x in qn_in:
            h = ''
            for y in x:
                if y >= 10 and y < 100:
                    h += ' {0:d}'.format(y)
                elif y >= 100:
                    h += '{0:d}'.format(y)
                else:
                    h += '  {0:d}'.format(y)
            h = strg.ljust(h, 36)
            qn_lin.append(h)
        return qn_lin

def standard_qn(qn_in, src='cat'):
    '''
    Converts the quantum numbers extrcted from the lin or cat file into the
    standard quantum number format

    Parameters
    ----------
    qn_in: (array, str): 2D array of the quantum numbers in the lin or cat format
    src (str): source of the data: cat or lin (cat)

    Returns
    -------
    qn_in (array, str): 2D array of the quantum numbers in the standard format

    Notes
    -----
    none
    '''

    if src == 'cat':
        inc = 2
    elif src == 'lin':
        inc = 3
    flag = 0

    if len(qn_in) == 0:
        print 'No data available!'
    else:
        qn = np.array([[]])
        for x in qn_in:
            i = 0
            check = '1'
            h = list()
            while check != '':
                check = strg.rstrip(x)[i:i+inc]
                if src == 'cat' and check != '' and check != '  ':
                    h.append(int(check))
                if src == 'lin' and check != '':
                    h.append(int(check))
                i += inc
            if flag == 0:
                qn = np.array(h, dtype=int)
                flag = 1
            else:
                qn = np.vstack((qn, np.array(h, dtype=int)))
        return qn

def get_abc_transitions(qn_in, lines, tr_type='a'):
    '''
    Extracts and returns a specific transition type from a linelist.

    Parameters
    ----------
    qn_in (array, in): quantum numbers in the standard format
    lines (list, np.array, float): 1 or 2 column array with frequencies in the
        first column and intensities in the second column or a list.
    tr_type (str): transition type to extrct: a, b, c (a)

    Returns
    -------
    qn (array, int): extracted quantum numbers for the specific transition type
    lin (array, float): linelist of the extracted transitions

    Notes
    -----
    none
    '''

    lin = reshape_linelist(lines)
    tr_ty = np.array([])

    if len(qn_in[0]) == 6:
        up = 3
    elif len(qn_in[0]) == 8:
        up = 4
    elif len(qn_in[0]) == 10:
        up = 5

    for x in qn_in:
        if np.mod(abs(x[1]-x[1+up]), 2) == 0 and np.mod(abs(x[2]-x[2+up])-1, 2) == 0:
            tr_ty = np.append(tr_ty, 'a')
        elif np.mod(abs(x[1]-x[1+up])-1, 2) == 0 and np.mod(abs(x[2]-x[2+up])-1, 2) == 0:
            tr_ty = np.append(tr_ty, 'b')
        elif np.mod(abs(x[1]-x[1+up])-1, 2) == 0 and np.mod(abs(x[2]-x[2+up]), 2) == 0:
            tr_ty = np.append(tr_ty, 'c')

    tr_ty = np.where(tr_ty == tr_type)

    return qn_in[tr_ty], lin[tr_ty]

def q_branch(qn_in, lines, branch = 'p'):
    '''
    Extracts and returns transitions of a certain branch (p or q).

    Parameters
    ----------
    qn_in (array, in): quantum numbers in the standard format
    lines (list, np.array, float): 1 or 2 column array with frequencies in the
    first column and intensities in the second column or a list.
    branch (str): specifies the branch of lines to be returned (p)

    Returns
    -------
    qn (array, int): extracted quantum numbers for the specific transition type
    lin (array, float): linelist of the extracted transitions

    Notes
    -----
    none
    '''

    lin = reshape_linelist(lines)
    tr_ty = np.array([])

    if len(qn_in[0]) == 6:
        up = 3
    elif len(qn_in[0]) == 8:
        up = 4
    elif len(qn_in[0]) == 10:
        up = 5

    for x in qn_in:
        if x[0] - x[up] == 0:
            tr_ty = np.append(tr_ty, 'q')
        else:
            tr_ty = np.append(tr_ty, 'p')

    if branch == 'q':
        tr_ty = np.where(tr_ty == 'q')
    else:
        tr_ty = np.where(tr_ty == 'p')

    return qn_in[tr_ty], lin[tr_ty]

def fit_min(spfitPred):
    '''
    Writes the par, int and lin file and runs the spfit program afterwards
    in the working directory. Finally it reads in the fit and the par file.
    The excecutable has to be in this directory or defined globally.

    Parameters
    ----------
    spfitPred (spfitPrediction): spfit prediction object

    Returns
    -------
    returns the modified spfit prediction object (spfitPred)

    Notes
    -----
    none
    '''

    spfitPred.write_par()
    spfitPred.write_int()
    spfitPred.write_lin()
    spfitPred.fit_run()
    spfitPred.read_fit()
    spfitPred.read_par()

    return spfitPred

def cat_min(spfitPred):
    '''
    Writes the var and int file and runs the spcat program afterwards
    in the working directory. Finally it reads in the cat file.
    The excecutable has to be in this directory or defined globally.

    Parameters
    ----------
    spfitPred (spfitPrediction): spfit prediction object

    Returns
    -------
    returns the modified spfit prediction object (spfitPred)

    Notes
    -----
    none
    '''

    spfitPred.int_params['flags'] = '0000';
    spfitPred.write_var()
    spfitPred.write_int()
    spfitPred.cat_run()
    spfitPred.read_cat()

    return spfitPred

def cat_max(spfitPred):
    '''
    Writes the var and int file and runs the spcat program afterwards
    in the working directory. Finally it reads in the cat, egy and str
    files. The excecutable has to be in this directory or defined globally.

    Parameters
    ----------
    spfitPred (spfitPrediction): spfit prediction object

    Returns
    -------
    returns the modified spfit prediction object (spfitPred)

    Notes
    -----
    none
    '''
    
    spfitPred.write_var()
    spfitPred.write_int()
    spfitPred.cat_run()
    spfitPred.read_cat()
    spfitPred.read_egy()
    spfitPred.read_str()

    return spfitPred


def refine_lines(linelist, spectrum, thresh=0.050, weight=False):
    '''
    Compares the assigned lines with a spectrum and refines the lines by looking
    for the absolute maximum intensity or by also weighting the neighboring
    points of the maximum.

    Parameters
    ----------
    spectrum (float): 2 column array, frequencies in the first column and
    intensities in the second column

    Returns
    -------
    freq (float): list of the refined frequencies
    intens ( float): list of the intensities

    Notes
    -----
    none
    '''

    if len(linelist) == 0:
        print "No data available!"
    else:
        delta_f = spectrum[1,0]- spectrum[0,0]
        freq = list()
        intens = list()
        for x in linelist:
            if x < spectrum[-1,0]:
                h = spectrum[np.where(np.abs(spectrum[:,0] - x) < thresh)]
                peak_index = int(np.where(h[:, 1] == np.max(h[:, 1]))[0][0])
                peak = h[peak_index, 1]
                peak_freq = h[peak_index, 0]
                if peak_index == 0:
                    peak_low = 0.
                    peak_up = h[peak_index + 1, 1]
                elif peak_index == len(h) - 1:
                    peak_up = 0.
                    peak_low = h[peak_index - 1, 1]
                else:
                    peak_up = h[peak_index + 1, 1]
                    peak_low = h[peak_index - 1, 1]

                if weight == True:
                    if peak_up > peak_low:
                        shift = (peak-peak_up)/(peak-peak_low)
                        shift = (1.-shift)*delta_f/2.
                    else:
                        shift = (peak-peak_low)/(peak-peak_up)
                        shift = -(1.-shift)*delta_f/2.

                    peak_freq += shift
            else:
                peak_freq = x
                peak = 0.0
            freq.append(peak_freq)
            intens.append(peak)

        return freq, intens
