from PyQt4 import QtGui, QtCore
from Ui.CompactMainWindow import Ui_CompactMainWindow
from Ui.CompactMeasurementMainWindow import Ui_CompactMeasurementMainWindow
from Ui.CompactMeasurementDialog import Ui_CompactMeasurementDialog
from Ui.CompactSettingsDialog import Ui_CompactSettingsDialog

from matplotlib.figure import Figure
#from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt4agg import (FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
#from matplotlib.backends import qt4_compat
        
           
class CompactMainWindow(QtGui.QMainWindow, Ui_CompactMainWindow):   #or whatever Q*class it is
    def __init__(self):
        super(CompactMainWindow, self).__init__()
        self.setupUi(self)

        self.fig = Figure((3.0, 2.0), dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.centralWidget)
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.canvas.setFocus()
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.centralWidget)
        self.verticalLayout.addWidget(self.canvas)  # the matplotlib canvas
        self.verticalLayout.addWidget(self.mpl_toolbar)
        #self.measDia.clicked.connect(self.openMeasDia)
        self.measMain.clicked.connect(self.openMeasMain)
      
    def openMeasMain(self):
        # here put the code that creates the new window and shows it.
        self.MeasMain = CompactMeasurementMainWindow()
        #self.MeasMain.setupUi(self)
        self.MeasMain.show()




class CompactMeasurementMainWindow(QtGui.QMainWindow, Ui_CompactMeasurementMainWindow):   #or whatever Q*class it is
    def __init__(self, parent=None):
        super(CompactMeasurementMainWindow, self).__init__(parent)
        self.setupUi(self)
        self.fig = Figure((3.0, 2.0), dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.centralWidget)
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.canvas.setFocus()
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.centralWidget)
        self.verticalLayout.addWidget(self.canvas)  # the matplotlib canvas
        self.verticalLayout.addWidget(self.mpl_toolbar)

#    def create_child(self):   #here should go your mbutton1

class CompactMeasurementDialog(QtGui.QDialog, Ui_CompactMeasurementDialog):   #or whatever Q*class it is
    def __init__(self, parent=None):
        super(CompactMeasurementDialog, self).__init__(parent)
        self.setupUi(self)

class CompactSettingsDialog(QtGui.QMainWindow, Ui_CompactSettingsDialog):   #or whatever Q*class it is
    def __init__(self, parent=None):
        super(CompactSettingsDialog, self).__init__(parent)
        self.setupUi(self)
