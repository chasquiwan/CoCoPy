import numpy as np
import string as strg
import spfit
import os

class xiamParams:

    def __init__(self):
        
        label = list(('ncyc', 'c_print', 'aprint', 'xprint', 'ints', 'reduc', 'freq_l', 
                                     'freq_h', 'limit', 'temp', 'orger', 'c_eval', 'dfreg', 'maxm',
                                     'woods', 'ndata', 'nfold', 'spin', 'ntop', 'adjf', 'rofit', 'defer',
                                     'eps', 'weigf', 'convg', 'c_lambda', 'fitscl', 'svderr'))        
        self.gen_params = dict(ncyc=0, c_print=3, aprint=9, xprint=4, ints=3, reduc=0, freq_l=2.0, 
                               freq_h=8.5, limit=0.1, temp=5.0, orger=0, c_eval=0, dfreg=0, maxm=8,
                               woods=33, ndata=0, nfold=3, spin=1, ntop=1, adjf=0, rofit=0.0, defer=1.0E-5,
                               eps=1.0E-12, weigf=0.0, convg=0.999, c_lambda=1.0E-5, fitscl=0, svderr=0, label=label, ex_zero=0)                               

        label = list(('BJ', 'BK', 'B_m', 'DJ', 'DJK', 'DK', 'dj', 'dk', 
                                        'H_J', 'HJK', 'HKJ', 'H_K', 'h_j', 'hjk', 'h_k', 'C_p', 'C_z',
                                        'C_m'))              
        self.rr_fit_params = dict(BJ=1.0, BK=1.0, B_m=1.0, DJ=0.0, DJK=0.0, DK=0.0, dj=0.0, dk=0.0, 
                               H_J=0.0, HJK=0.0, HKJ=0.0, H_K=0.0, h_j=0.0, hjk=0.0, h_k=0.0, C_p=0.0, C_z=0.0,
                               C_m=0.0, label=label, ex_zero=1)                               

        label = list(('chi_zz', 'chi_m', 'chi_xy', 'chi_xz', 'chi_yz'))
        self.qc_fit_params = dict(chi_zz=1.0, chi_m=0.0, chi_xy=0.0, chi_xz=0.0, chi_yz=0.0, label=label, ex_zero=1)
        
        label = list(('C_p', 'C_z', 'C_m'))
        self.sr_fit_params = dict(C_p=0.1, C_z=0.0, C_m=0.0, label=label, ex_zero=1)        
        
        label = list(('V1n', 'V2n', 'F', 'rho', 'beta', 'gamma', 'Dpi2J', 'Dpi2K', 'Dpi2_m',
                                         'Dpi3J', 'F0', 'delta', 'epsil'))          
        self.ir_fit_params = dict(V1n=1.0, V2n=0.0, F=0.0, rho=0.0, beta=1.0, gamma=0.0, Dpi2J=0.0, Dpi2K=0.0, Dpi2_m=0.0,
                                  Dpi3J=0.0, F0=0.0, delta=0.0, epsil=0.0, label=label, ex_zero=1)                                  
                                  
        label = list(('mu_x', 'mu_y', 'mu_z'))
        self.dm_fit_params = dict(mu_x=0.0, mu_y=0.0, mu_z=0.0, label=label, ex_zero=0)
        
        self.params = list((self.gen_params, self.rr_fit_params, self.qc_fit_params, 
                            self.sr_fit_params, self.ir_fit_params, self.dm_fit_params))
    

class xiamPrediction:
    '''
    help
    '''
        
    def __init__(self,  fname = 'test', name = 'Test'):
             
        self.params = xiamParams()
        
        self.fit = list(('BJ', 'BK', 'B_m'))
        self.dqu = list(('BJ', 'BK', 'B_m'))
        self.dqx = list(('BJ', 'BK', 'B_m'))
        
        self.tor_symm = dict(label=list(('/A', '/E')), value=list(('S 0', 'S 1')))
        
        self.tor_state = list(('V 0'))
        
        self.mol_params = dict(nqn=3, nline=100, nvib=1, nvarfit=5, linerr=0.02, 
                               xi_ex=0, xo_ex=0, name='test', fname='test', gen=1, rr=1, qc=1, sr=1, ir=1)
            
        self.meas_spec = dict(spec=np.array([[]]), linelist=np.array([[]]), qn=np.array([[]]),
                              scale=1.0)
    
        self.pred_spec = dict(spec=np.array([[]]), linelist=np.array([[]]), qn=np.array([[]]), 
                              scale=1.0)


        self.mol_params['name'] = name
        self.mol_params['fname'] = fname

        
    def write_xi(self, fname=''):
        '''
        
        Parameters
        ----------
                
        Returns
        -------
        
        Notes
        -----
        none
        '''
        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + '.xi', 'w+')
        except IOError:
            self.mol_params['xi_ex'] = 0
            print 'Cannot open: ', self.mol_params['fname'] + '.xi'
        else:
            self.mol_params['xi_ex'] = 1
            # 1. Block
            filestr = self.mol_params['name'] + '\n\n'
            # 2. and 3. Block
            filestr += self.write_params()
            filestr += '\n'
            # 4. Block
            filestr += self.write_fit()
            filestr += '\n'
            # 5. Block
            filestr += self.write_symm()
            filestr += '\n'

            f.write(filestr)
            f.close()


    def write_params(self):

        filestr = ''
        for p in self.params.params:
            for key in p['label']:
                h_key = convert_label(key)
                if p['ex_zero'] == 1:
                    if float(p[key]) != 0.0:
                        filestr += (strg.ljust(str(h_key), 12) + str(p[key]) + '\n')
                else:
                    filestr += (strg.ljust(str(h_key), 12) + str(p[key]) + '\n')
                if key == 'svderr':
                    filestr += '\n'
        return filestr

    def write_fit(self):

        filestr = ''
        for f in self.fit:
            f = convert_label(f)
            filestr += 'fit ' + str(f) + '\n'
        for f in self.dqu:
            f = convert_label(f)
            filestr += 'dqu ' + str(f) + '\n'
        for f in self.dqx:
            f = convert_label(f)
            filestr += 'dqx ' + str(f) + '\n'
        return filestr
        
    def write_symm(self):

        filestr = ''
        i=0
        for s in self.tor_symm['label']:
            filestr += str(s) + '  ' + self.tor_symm['value'][i] + '\n'
            i += 1
        return filestr
        
    def write_state(self):

        filestr = ''
        i=0
        for s in self.tor_state:
            filestr += str(s) + '  ' + self.tor_symm['value'][i] + '\n'
            i += 1
        return filestr

    def import_from_spfit(self, fname='test', src='lin'):
        '''
        
        Parameters
        ----------
        fname (str):
        src (str): 
                
        Returns
        -------
        none
        
        Notes
        -----
        none

        '''
        f = open(fname+'.xi', 'w')
        a = spfit.spfitPrediction(fname=fname)
        if src == 'lin':
            a.read_lin()
            a.make_linel_lin()
            qn = a.meas_spec['qn']
            freq = a.meas_spec['linelist']
        if src == 'cat':
            a.read_cat()
            a.make_linel_cat()
            qn = a.pred_spec['qn']
            freq = a.pred_spec['linelist']
        filestr = ''
        i = 0
        for x in qn:
            filestr += (strg.rjust(str(x[0]), 2) + strg.rjust(str(x[1]), 3) + strg.rjust(str(x[2]), 3) + '  ')
            filestr += (strg.rjust(str(x[3]), 3) + strg.rjust(str(x[4]), 3) + strg.rjust(str(x[5]), 3))
            filestr += ('  S 1 =  ' + '{0:1.8f}'.format(freq[i][0]/1000.) + ' Err 2.0E-6 \n')
            i += 1
        
        f.write(filestr)
        f.close()


        
    def read_xo(self, fname=''):
        '''
        Reads the predicted lines of a XIAM prediction. Currently it only can read
        the predicted lines for an asymmetric top without quadrupole coupling and
        one top (3 fold)        
        
        Parameters
        ----------
        fname (str): Filename of xxx.xo file
                
        Returns
        -------
        linelist (dict): This dictionary contains the frequencies (freq), 
        intensities (intens), symmetries (sym) and quantum numbers (qn) of the
        predicted spectrum
        
        Notes
        -----
        none

        '''
        try:
            if fname == '':
                fname = self.mol_params['fname']
            else:
                self.mol_params['fname'] = fname
            f = open(fname + '.xo', 'r')
        except IOError:
            print 'Cannot open: ', fname
        else:
            mode = -1
            check = False
            flag = False 
            linelist = dict(qn=np.array([]), sym=np.array([]), freq=np.array([]), intens=np.array([]), obs=np.array([]))
            for line in f:
                if ' ints         2' in line:
                    mode = 2
                if ' ints         3' in line:
                    mode = 3
                elif 'End at Cycle' in line:
                    if flag == True:
                        mode = 0
                    else:
                        mode = 0
                elif 'Recalculation of the spectrum' in line:
                    flag = True
                if mode == 2:
                    if check == True:
                        if 'total is the product of Linestr.' in line:
                            check = False
                        if check == True:
                            if line != '\n':
                                linelist['qn'] = [line[0:19]] if len(linelist['qn']) == 0 else np.vstack((linelist['qn'], line[0:19]))
                                linelist['freq'] = [float(line[34:43])] if len(linelist['freq']) == 0 else np.vstack((linelist['freq'], float(line[34:43])))
                                linelist['intens'] = [float(line[43:52])] if len(linelist['intens']) == 0 else np.vstack((linelist['intens'], float(line[43:52])))
                                if 'S 1' in line[22:30]:
                                    linelist['sym'] = ['A'] if len(linelist['sym']) == 0 else np.vstack((linelist['sym'], 'A'))
                                else:
                                    linelist['sym'] = ['E'] if len(linelist['sym']) == 0 else np.vstack((linelist['sym'], 'E'))
                                
                    if '--' in line:
                        check = True
                        
                if mode == 3:
                    if check == True:
                        if 'total is the product of Linestr.' in line:
                            check = False
                        if check == True:
                            if line != '\n' and 'S ' in line[22:25]:
                                if line[0:19] != strg.ljust('', 19):
                                    linelist['qn'] = [line[0:19]] if len(linelist['qn']) == 0 else np.vstack((linelist['qn'], line[0:19]))
                                    qn = line[0:19]
                                else:
                                    linelist['qn'] = np.vstack((linelist['qn'], qn))
                                linelist['freq'] = [float(line[39:48])] if len(linelist['freq']) == 0 else np.vstack((linelist['freq'], float(line[39:48])))
                                linelist['intens'] = [float(line[61:67])] if len(linelist['intens']) == 0 else np.vstack((linelist['intens'], float(line[61:67])))
                                linelist['sym'] = [int(line[24:25])] if len(linelist['sym']) == 0 else np.vstack((linelist['sym'], int(line[24:25])))
                                 
                    if '--' in line:
                        check = True
                        
                if mode == 0:
                    if check == True:
                        if  'Maximum (obs-calc)/err' in line:
                            check = False
                        if check == True:
                            linelist['qn'] = [line[5:23]] if len(linelist['qn']) == 0 else np.vstack((linelist['qn'], line[5:23]))
                            linelist['freq'] = [float(line[28:40])] if len(linelist['freq']) == 0 else np.vstack((linelist['freq'], float(line[28:40])))
                            linelist['intens'] = [1.0] if len(linelist['intens']) == 0 else np.vstack((linelist['intens'], 1.0))
                            linelist['obs'] = [float(line[50:62])] if len(linelist['obs']) == 0 else np.vstack((linelist['obs'], float(line[50:62])))
                            if '/A' in line[23:28]:
                                linelist['sym'] = ['A'] if len(linelist['sym']) == 0 else np.vstack((linelist['sym'], 'A'))
                            else:
                                linelist['sym'] = ['E'] if len(linelist['sym']) == 0 else np.vstack((linelist['sym'], 'E'))
                    if 'calc/GHz' in line:
                        check = True

                if mode == 1:
                    if check == True:
                        if  'Maximum (obs-calc)/err' in line:
                            check = False
                        if check == True:
                            linelist['qn'] = [line[5:23]] if len(linelist['qn']) == 0 else np.vstack((linelist['qn'], line[5:23]))
                            linelist['freq'] = [float(line[28:40])] if len(linelist['freq']) == 0 else np.vstack((linelist['freq'], float(line[28:40])))
                            linelist['intens'] = [float(line[84:91])] if len(linelist['intens']) == 0 else np.vstack((linelist['intens'], float(line[84:91])))
                            linelist['obs'] = [float(line[50:62])] if len(linelist['obs']) == 0 else np.vstack((linelist['obs'], float(line[50:62])))
                            if '/A' in line[23:28]:
                                linelist['sym'] = ['A'] if len(linelist['sym']) == 0 else np.vstack((linelist['sym'], 'A'))
                            else:
                                linelist['sym'] = ['E'] if len(linelist['sym']) == 0 else np.vstack((linelist['sym'], 'E'))
                    if 'calc/GHz' in line:
                        check = True
                        

            linelist['freq'] *= 1000. 
            linelist['obs'] *= 1000. 
            return linelist
        
    def run(self, fname=''):
        check = -1
        if fname == '':
            fname = self.mol_params['fname']
        else:
            self.mol_params['fname'] = fname
        check = os.system('xiam < ' + fname + '.xi > ' + fname + '.xo')
        return check



def convert_label(label):
    if strg.find(label, 'c_') > -1:
        h_label = label[2:]
    elif strg.find(label, '_m') > -1:
        h_label = label[0:-2] + '-'
    elif strg.find(label, '_p') > -1:
        h_label = label[0:-2] + '+'    
    else:
        h_label = label
        
    return h_label


#def xiam_readout_diff(filename):
#
#    f = open(filename, "r")
#    check = 0
#    trans = np.array([" "])
#    final = np.array([[" ", "  "]])
#    
#    for line in f:
#        if "End at Cycle" in line:
#            check = 1
#        if "Maximum" in line:
#            check = 0
#        if check == 1:
#            trans = np.append( trans, [line], axis = 1 )
#        
#    trans = np.split( trans, [4], axis = 0 )[1]
#    
#    i = 0
#    for x in trans:
#        if i % 2 == 0:    
#            qn_low = "  " + x[16] + " " + x[19] + " " + x[22]
#            freqA = x[31:40]
#        if i % 2 > 0:
#            freqE = x[31:40]
#            diff_freq = str( float( freqE) - float( freqA ) )
#            final = np.append( final, [[qn_low, diff_freq]], axis = 0 )
#        i += 1
#        
#    final = np.split( final, [1], axis = 0 )[1]
#    
#    return final