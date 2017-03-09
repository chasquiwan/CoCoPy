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
from __future__ import print_function

import numpy as np
import io
import scipy.constants as constants
import sys
import glob


data_path = '/Volumes/DATA/'
lines = dict()
lines['acetone'] = [5253.063, 5269.071, 5270.904, 5276.057, 7343.774, 7398.953, 7399.964, 7453.183]
lines['background'] = [
2050.00, 2060.00, 2070.00, 2090.00, 2100.00, 2187.51, 2200.01, 2250.00, 2260.00, 2270.00,
2280.00, 2290.00, 2310.00, 2330.00, 2339.99, 2349.99, 2359.99, 2369.99, 2379.99, 2390.01,
2409.99, 2420.01, 2437.51, 2450.00, 2460.00, 2470.00, 2480.00, 2500.00, 2540.02, 2549.99,
2559.99, 2569.99, 2600.01, 2650.01, 2660.00, 2670.00, 2690.00, 2700.00, 2750.00, 2754.99,
2799.99, 2810.01, 2812.51, 2890.00, 2900.00, 2920.00, 2960.00, 2964.99, 2969.99, 2979.99,
2989.99, 2999.99, 3005.01, 3010.01, 3020.01, 3030.01, 3045.01, 3050.01, 3060.01, 3090.00,
3100.00, 3110.00, 3125.00, 3187.49, 3245.01, 3250.01, 3360.00, 3399.99, 3437.51, 3490.00,
3495.00, 3500.00, 3515.00, 3520.00, 3710.00, 3750.00, 3799.99, 3809.99, 4125.00, 4687.51,
5840.00, 5860.00, 5870.00, 5874.99, 5879.99, 6000.00, 6010.00, 6020.00, 2010.00, 2020.00,
6089.99, 6130.01, 6170.01, 6249.98, 6250.00, 6330.01, 6370.01, 6410.00, 6620.00, 6630.00,
6640.00, 7812.51, 2300.00, 4500.00, 7500.00, 2130.00, 2240.00, 2110.00]

'''
Todo:
    - 1 Write comments
    - 2 Write parser for out file
    - 3
'''

def stylePlot(data='spectrum'):

    import matplotlib.pyplot as plt

    ax = plt.gca()

    if data == 'spectrum':
        xlabel = r'Frequency [MHz]'
        ylabel = r'Intensity [arb. u.]'
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

    elif data == 'fid':
        xlabel = r'Time [$\mu s$]'
        ylabel = r'Intensity [arb. u.]'
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ax.ticklabel_format(axis='x', style='plain', useOffset=False)

    elif data == 'intensity':
        xlabel = r'Frequency [MHz]'
        ylabel = r'Intensity [arb. u.]'
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

    elif data == 'energy':
        xlabel = r'Angle [$\degree$]'
        ylabel = r'Energy [kcal/mol]'
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-4,4), useOffset=False)

    elif data == 'phase':
        xlabel = r'Frequency [MHz]'
        ylabel = r'Phase [radians]'
        plt.ylim(-3.5, 3.5)
        plt.yticks([-np.pi, -np.pi/2, 0., np.pi/2, np.pi], [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'])

    plt.title(r'Plot', fontsize=20, family='serif', position=(0.5, 1.05))
    plt.xlabel(xlabel, fontsize=20, horizontalalignment='center')
    plt.ylabel(ylabel, fontsize=20)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    ax.yaxis.set_tick_params(pad=8)
    ax.xaxis.set_tick_params(pad=8)
    ax.xaxis.set_label_coords(0.5, -0.1)
    ax.yaxis.set_label_coords(-0.07, .5)

def pubPlot(width = 384):

    import matplotlib.pyplot as plt

    fig_width_pt = width  # Get this from LaTeX using \the\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height =fig_width*golden_mean       # height in inches
    fig_size = [fig_width,fig_height]
    params = {'backend': 'pdf',
              'axes.labelsize': 10,
              'text.fontsize': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': fig_size,
              'font.size': 8,
              'font.family': 'serif',
              'font.serif': 'serif',
              'pdf.fonttype':42}
    plt.rcParams.update(params)
    plt.figure(1)
    plt.clf()
    plt.axes([0.125,0.2,0.95-0.125,0.95-0.2])


def execute_notebook(nbfile):
    # execute_notebook("loaddata.ipynb")

    from IPython import nbformat

    with io.open(nbfile) as f:
        nb = nbformat.read(f, 'json')

    ip = get_ipython()

    for cell in nb.worksheets[0].cells:
        if cell.cell_type != 'code':
            continue
        ip.run_cell(cell.input)

def convert_energy(energy, origin='hartree', dest='wavenumber'):
    #http://docs.scipy.org/doc/scipy/reference/constants.html#module-scipy.constants
    energy_h = 0.
    if origin == 'hartree':
        energy_h = constants.physical_constants['hartree-joule relationship'][0] * energy
    elif origin == 'wavenumber' or origin == 'cm-1':
        energy_h = constants.physical_constants['inverse meter-joule relationship'][0] * energy * 100
    elif origin == 'joule' or origin == 'J':
        energy_h = 1. * energy
    elif origin == 'hertz' or origin == 'Hz':
        energy_h = constants.physical_constants['hertz-joule relationship'][0] * energy
    elif origin == 'eV':
        energy_h = constants.physical_constants['electron volt-joule relationship'][0] * energy
    elif origin == 'wavelength' or origin == 'm':
        energy_h = constants.c * constants.h / energy
    elif origin == 'kCal/mol':
        energy_h = energy * 4184.0 / constants.Avogadro
    elif origin == 'kJ/mol':
        energy_h = energy * 1000. / constants.Avogadro

    if dest == 'hartree':
        energy_h = constants.physical_constants['joule-hartree relationship'][0] * energy_h
    elif dest == 'wavenumber' or dest == 'cm-1':
        energy_h = constants.physical_constants['joule-inverse meter relationship'][0] * energy_h / 100.
    elif dest == 'joule' or dest == 'J':
        energy_h = 1. * energy_h
    elif dest == 'hertz' or dest == 'Hz':
        energy_h = constants.physical_constants['joule-hertz relationship'][0] * energy_h
    elif dest == 'eV':
        energy_h = constants.physical_constants['joule-electron volt relationship'][0] * energy_h
    elif dest == 'wavelength' or dest == 'm':
        energy_h = constants.c * constants.h / energy_h
    elif dest == 'kCal/mol':
        energy_h = energy_h / 4184.0 * constants.Avogadro
    elif dest == 'kJ/mol':
        energy_h = energy_h / 1000. * constants.Avogadro

    return energy_h

def convert_rot(value):
    #http://docs.scipy.org/doc/scipy/reference/constants.html#module-scipy.constants
    return constants.h / (8. * np.pi**2 * value * 1.E6 * constants.atomic_mass) * 1.E20

class ProgressBar:
    def __init__(self, iterations):
        self.iterations = iterations
        self.prog_bar = '[]'
        self.fill_char = '*'
        self.width = 50
        self.__update_amount(0)

    def animate(self, iter):
        print('\r', self, end='')
        sys.stdout.flush()
        self.update_iteration(iter + 1)

    def update_iteration(self, elapsed_iter):
        self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
        self.prog_bar += '  %d of %s complete' % (elapsed_iter, self.iterations)

    def __update_amount(self, new_amount):
        percent_done = int(round((new_amount / 100.0) * 100.0))
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
        pct_string = '%d%%' % percent_done
        self.prog_bar = self.prog_bar[0:pct_place] + \
            (pct_string + self.prog_bar[pct_place + len(pct_string):])

    def __str__(self):
        return str(self.prog_bar)

def draw_cube(ax, size, pos=[0,0,0], edges = False, color = 'red', alpha=0.1):

    import matplotlib.pyplot as plt
    import matplotlib.colors as co
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from itertools import product

    h = list(np.array([-1., 1.])*size)

    tupleList=np.array([x for x in product(h,h,h)])
    tupleList[:,0] +=pos[0]; tupleList[:,1] += pos[1]; tupleList[:,2] += pos[2]

    vertices = [[0, 1, 3, 2], [0, 1, 5, 4], [0, 2, 6, 4], [2, 3, 7, 6], [4, 5, 7, 6], [1, 3, 7, 5]]

    cc = lambda arg: co.colorConverter.to_rgba(arg, alpha=alpha)

    poly3d = [[tupleList[vertices[ix][iy]] for iy in range(len(vertices[0]))] for ix in range(len(vertices))]
    if edges == True:
        ax.scatter(tupleList[:,0],tupleList[:,1],tupleList[:,2], c=color)
    collection = Poly3DCollection(poly3d, linewidths=.2, linestyles='-', edgecolors=['grey'], facecolors = [cc('red')])

    #j = co.ColorConverter()
    #face_color = co.rgb2hex(j.to_rgb(color))
    #collection.set_facecolor(face_color)
    ax.add_collection3d(collection)
