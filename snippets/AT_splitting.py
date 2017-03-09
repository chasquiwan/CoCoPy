# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 14:40:54 2015

@author: davidschmitz
"""

# Load some modules
import scipy.constants as con
from scipy.integrate import ode
import copy as cp
import tempfile
import numpy as np
import sys

tempfile.tempdir = sys.argv[2]

# Setting up the differential equations by solving the commutator
def f_3(t, y, arg):
    
    # Make it human readable
    E_a = arg[0]; E_b = arg[1]; E_c = arg[2] 
    mu_ab = arg[3]; mu_bc = arg[4]; mu_ac = arg[5]
    F_ab = arg[6]; F_bc = arg[7]; F_ac = arg[8]
    o_ab = arg[9]; o_bc = arg[10]; o_ac = arg[11]
    phi_ab = arg[12]; phi_bc = arg[13]; phi_ac = arg[14]
    
    # Set up Hamiltonian
    H = np.zeros((3,3), dtype=np.complex)
    H[0,0] = E_a; H[1,1] = E_b; H[2,2] = E_c;
    
    H[0,1] = - mu_ab * F_ab / 2. * (np.exp(-1.j*(o_ab*t+phi_ab)) + np.exp(1.j*(o_ab*t+phi_ab))); H[1,0] = np.conjugate(H[0,1])
    H[1,2] = - mu_bc * F_bc / 2. * (np.exp(-1.j*(o_bc*t+phi_bc)) + np.exp(1.j*(o_bc*t+phi_bc))); H[2,1] = np.conjugate(H[1,2])
    H[0,2] = - mu_ac * F_ac / 2. * (np.exp(-1.j*(o_ac*t+phi_ac)) + np.exp(1.j*(o_ac*t+phi_ac))); H[2,0] = np.conjugate(H[0,2])
    
    # Set up density matrix
    rho = np.zeros((3,3), dtype=np.complex)
    rho[0,0] = y[0]; rho[0,1] = y[1]; rho[0,2] = y[2]
    rho[1,0] = y[3]; rho[1,1] = y[4]; rho[1,2] = y[5]
    rho[2,0] = y[6]; rho[2,1] = y[7]; rho[2,2] = y[8]

    # Calculate commutator
    dt_rho = -1.j/con.hbar * ( np.dot(H, rho) - np.dot(rho, H) )
    return [[dt_rho[0,0]], [dt_rho[0,1]], [dt_rho[0,2]], 
            [dt_rho[1,0]], [dt_rho[1,1]], [dt_rho[1,2]], 
            [dt_rho[2,0]], [dt_rho[2,1]], [dt_rho[2,2]]]
            
# rho is the density matrix at every timestep
# t is the time vector
# This function solves the Bloch equations and carries the full output
def bloch_solver_3(par, dump=True, cut=False, Ndpoints=1000, Ncpoints=1000):

    para = cp.deepcopy(par)
    if dump == True:
        tmpfile_t = tempfile.TemporaryFile()
        tmpfile_rho = tempfile.TemporaryFile()
        cut = False
    
    def set_init(para, i):
        return [para['E_a'], para['E_b'], para['E_c'], para['mu_ab'], para['mu_bc'], para['mu_ac'],
                para['F_ab'][i], para['F_bc'][i], para['F_ac'][i], para['o_ab'][i], para['o_bc'][i], para['o_ac'][i],
                para['phi_ab'][i], para['phi_bc'][i], para['phi_ac'][i]]
    
    def save_dump(rho, t, tmpfile_t, tmpfile_rho):
        tmpfile_t.seek(0); tmpfile_rho.seek(0)
        t_h = np.load(tmpfile_t); rho_h = np.load(tmpfile_rho)
        tmpfile_t.seek(0); tmpfile_rho.seek(0)
        rho = np.vstack((rho_h, rho[:-1])); t = np.append(t_h, t[:-1])
        np.save(tmpfile_t, t); np.save(tmpfile_rho, rho)
    
    para['o_ab'] *= (2.*np.pi*1.E6); para['o_bc'] *= (2.*np.pi*1.E6); para['o_ac'] *= (2.*np.pi*1.E6)
    para['mu_ab'] *= 3.335E-30; para['mu_bc'] *= 3.335E-30; para['mu_ac'] *= 3.335E-30


    nb = len(para['tp'])
    r = ode(f_3).set_integrator('zvode')

    for i in np.arange(nb):
        
        # Set parameters and starting values at t0
        if i == 0:
            r.set_initial_value(para['rho0'], para['t0']).set_f_params(set_init(para, 0))
            t = np.array(para['t0'])
            rho = np.array([[[para['rho0'][0], para['rho0'][1], para['rho0'][2]], 
                          [para['rho0'][3], para['rho0'][4], para['rho0'][5]], 
                          [para['rho0'][6], para['rho0'][7], para['rho0'][8]]]])
            tpp = 0.
            
            if dump == True:
                np.save(tmpfile_rho, rho); np.save(tmpfile_t, t)

            
        # Set parameters for the next pulse
        else:
            r.set_initial_value(r.y, r.t).set_f_params(set_init(para, i))
            tpp += para['tp'][i-1]
            
        # Integrate and solve ODEs
        j = 0
        while r.successful() and r.t <= tpp + para['tp'][i]:
            r.integrate(r.t + para['dt'])
            
            if dump == True and j == Ndpoints:
                save_dump(rho, t, tmpfile_t, tmpfile_rho)
                rho = np.array([rho[-1]])
                t = np.array([t[-1]])
                j = 0
                
            if cut == False or r.t > np.sum(para['tp']) - Ncpoints*para['dt']:
                rho = np.vstack((rho, [[[r.y[0], r.y[1], r.y[2]], 
                                     [r.y[3], r.y[4], r.y[5]], 
                                     [r.y[6], r.y[7], r.y[8]]]]))
                t = np.append(t, r.t)
            
            j += 1
        
    if dump == True:
        save_dump(rho, t, tmpfile_t, tmpfile_rho)
        tmpfile_t.seek(0); tmpfile_rho.seek(0)
        rho = np.load(tmpfile_rho); rho = rho[1:] 
        t = np.load(tmpfile_t); t = t[1:]
        tmpfile_rho.close(); tmpfile_t.close()
            
    if cut == False:
        return t, rho
    else:
        return t[1:], rho[1:]

par = dict({ 'dt': 4.E-11, 'tp': np.array([2.50E-8, 50E-6]), 't0': 0.,
             'mu_ab': 1.+0.j, 'mu_bc': 1.+0.j, 'mu_ac':  1.+0.j,
             'E_a': 0., 'E_b': 2.460372e-24, 'E_c': 2.674924e-24,
             'o_ab': np.array([3713.17, 3713.17]), 'o_bc': np.array([323.8, 323.8]), 'o_ac': np.array([4036.97, 4036.97]),
             'F_ab': np.array([2000., 0.]), 'F_bc': np.array([0., 25.]), 'F_ac': np.array([0., 0.]),
             'phi_ab': np.array([0.,0.]), 'phi_bc': np.array([0.,0.]), 'phi_ac': np.array([0.,0.]),
             'rho0': np.array([.5 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, 
                          0.0 + 0.0j, 0.3 + 0.0j, 0.0 + 0.0j, 
                          0.0 + 0.0j, 0.0 + 0.0j, 0.2 + 0.0j]),
             })

det = float(sys.argv[1])

par['o_bc'] = np.array([323.8, 323.8]) + det
t, rho_n = bloch_solver_3(par, dump=True, Ndpoints=10000)
fid_n = np.hstack((np.reshape(t,(len(t), 1)), np.reshape(rho_n[:,0,1].real,(len(t), 1))))
fid_n = fid_n[-1250000:]
np.save('long'+str(int(det*1000.)), fid_n)

del fid_n; del rho_n; del t
