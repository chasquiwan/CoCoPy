import sys
from PyQt4 import QtGui, QtCore
from GUI.CompactGUI import CompactMainWindow, CompactMeasurementMainWindow, CompactMeasurementDialog, CompactSettingsDialog
import numpy

class Compact(CompactMainWindow):

    def f(self):
            
            self.fig.clear()
            self.ax = self.fig.add_subplot(1,1,1)
            x = numpy.arange(0., 1., .0001)
            y = numpy.sin(2. * x * self.i)
            self.ax.plot(x,y)
            self.fig.canvas.draw()
            self.i += 1.
            #print self.i

    def __init__(self):
        super(Compact, self).__init__()
        '''
        self.ui = CompactMainWindow()
        self.ui.setupUi(self)
        '''
        self.i = 0.

        self.ctimer = QtCore.QTimer()
        self.ctimer.timeout.connect(self.f)
        self.show()
        
        self.ctimer.start(10)

        
        
            