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
import pylab as plt
import string as strg

'''
Todo:
    - 1
    - 2
    - 3
'''

class ramMolecule:
    '''
    help
    '''
        
    def __init__(self, name = 'Test', fname = 'test'):
        '''
        help
        '''    
        self.pred_params = dict(temp=0.3, str_cutoff=0.001, int_cutoff=0.001, 
                                low_freq=1000.000, up_freq=2000.000, units='MHz', 
                                J_min=1, J_max=8, vt_min=0, vt_max=2, E_lowest=1.0)
                                
        self.prediction = np.array([])
        
        self.mol_params = dict(nline=100, linerr=0.02, flag_linerr=1, nintpar=9, 
                               pred_ex=0, input_ex=0, output_ex=0, name='test', 
                               fname='test')
        
        self.mol_params['name'] = name
        self.mol_params['fname'] = fname

        def read_pred(self, fname=''):
            '''
            help
            '''
            h = dict(sym_up='E2', m_up=0, J_up=0, Ka_up=0, Kc_up=0, sym_low='E2', 
                     m_down=0, J_down=0, Ka_down=0, Kc_down=0, intens=1.0E0, 
                     pred_freq=1000.00, pred_freq_err=0.0, E_low=1.0, strength=1.0, 
                     meas_freq=1000.00, meas_freq_err=0.0, pred_meas_diff=0.0, 
                     used_fit=1)

            try:
                if fname == '':
                    fname = 'predictvt0.txt'
                else:
                    self.mol_params['fname'] = fname
                f = open(fname, 'r')
            except IOError:
                self.mol_params['pred_ex'] = 0
                print 'Cannot open: ', fname
            else:
                i = 0
                for line in f:
                    if i == 0:
                        self.pred_params['temp'] = float(line[12:23])
                        self.pred_params['str_cutoff'] = float(line[42:57])
                        self.pred_params['int_cutoff'] = float(line[74:84])
                    if i == 1:
                        self.pred_params['low_freq'] = float(line[20:31])
                        self.pred_params['up_freq'] = float(line[47:59])
                    if i == 2:
                        self.pred_params['units'] = strg.replace(line[11:15], '\n', '')
                    if i == 3:
                        self.pred_params['J_min'] = int(line[5:11])
                        self.pred_params['J_max'] = int(line[15:21])
                        self.pred_params['vt_min'] = int(line[27:33])
                        self.pred_params['vt_max'] = int(line[38:44])
                    if i == 4:
                        self.pred_params['E_lowest'] = float(line[9:36])
                    if i > 9 and line != '\n':
                        h['sym_up'] = line[0:2]
                        h['m_up'] = int(line[2:5])
                        h['J_up'] = int(line[5:9])
                        h['Ka_up'] = int(line[9:13])
                        h['Kc_up'] = int(line[13:17])
                        h['sym_low'] = line[17:23]
                        h['m_low'] = int(line[23:26])
                        h['J_low'] = int(line[26:30])
                        h['Ka_low'] = int(line[30:34])
                        h['Kc_low'] = int(line[34:38])
                        h['intens'] = float(line[38:51])
                        h['pred_freq'] = float(line[51:65])
                        h['pred_freq_err'] = float(line[66:74])
                        h['E_low'] = float(line[75:87])
                        h['strength'] = float(line[87:97])
                        if line[97:111] != '' and line[97:111] != '\n':
                            h['meas_freq'] = float(line[97:111])
                        else:
                            h['meas_freq']=''
                        if line[112:119] != '':
                            h['meas_freq_err'] = float(line[112:119])
                        else:
                            h['meas_freq_err']=''
                        if line[120:131] != '':
                            h['pred_meas_diff'] = float(line[120:131])
                        else:
                            h['pred_meas_diff']=''
                        if line[131:134] != '':
                            h['used_fit'] = int(line[131:134])
                        else:
                            h['used_fit']=''
                        
                        if i == 10:
                            self.prediction = np.array([h])
                        else:
                            self.prediction = np.concatenate((self.prediction, [h]))
                i+=1
    
            f.close()            
    # Slicing: output[23]['intens']                            
            
