# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainwindow.ui'
#
# Created: Fri Jun 26 15:10:59 2015
#      by: PyQt4 UI code generator 4.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_CompactMainWindow(object):
    def setupUi(self, CompactMainWindow):
        CompactMainWindow.setObjectName(_fromUtf8("CompactMainWindow"))
        CompactMainWindow.resize(509, 417)
        self.centralWidget = QtGui.QWidget(CompactMainWindow)
        self.centralWidget.setObjectName(_fromUtf8("centralWidget"))
        self.measDia = QtGui.QPushButton(self.centralWidget)
        self.measDia.setGeometry(QtCore.QRect(200, 60, 221, 31))
        self.measDia.setObjectName(_fromUtf8("measDia"))
        self.measMain = QtGui.QPushButton(self.centralWidget)
        self.measMain.setGeometry(QtCore.QRect(200, 90, 221, 21))
        self.measMain.setObjectName(_fromUtf8("measMain"))
        self.widget = QtGui.QWidget(self.centralWidget)
        self.widget.setGeometry(QtCore.QRect(39, 139, 431, 201))
        self.widget.setObjectName(_fromUtf8("widget"))
        self.verticalLayoutWidget = QtGui.QWidget(self.widget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(-1, -1, 431, 201))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        CompactMainWindow.setCentralWidget(self.centralWidget)
        self.menuBar = QtGui.QMenuBar(CompactMainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 509, 22))
        self.menuBar.setObjectName(_fromUtf8("menuBar"))
        self.menuDatei = QtGui.QMenu(self.menuBar)
        self.menuDatei.setObjectName(_fromUtf8("menuDatei"))
        self.menuSettings = QtGui.QMenu(self.menuBar)
        self.menuSettings.setObjectName(_fromUtf8("menuSettings"))
        CompactMainWindow.setMenuBar(self.menuBar)
        self.mainToolBar = QtGui.QToolBar(CompactMainWindow)
        self.mainToolBar.setObjectName(_fromUtf8("mainToolBar"))
        CompactMainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QtGui.QStatusBar(CompactMainWindow)
        self.statusBar.setObjectName(_fromUtf8("statusBar"))
        CompactMainWindow.setStatusBar(self.statusBar)
        self.actionSettings = QtGui.QAction(CompactMainWindow)
        self.actionSettings.setObjectName(_fromUtf8("actionSettings"))
        self.menuSettings.addAction(self.actionSettings)
        self.menuBar.addAction(self.menuDatei.menuAction())
        self.menuBar.addAction(self.menuSettings.menuAction())

        self.retranslateUi(CompactMainWindow)
        QtCore.QMetaObject.connectSlotsByName(CompactMainWindow)

    def retranslateUi(self, CompactMainWindow):
        CompactMainWindow.setWindowTitle(_translate("CompactMainWindow", "MainWindow", None))
        self.measDia.setText(_translate("CompactMainWindow", "Measurement Dialog", None))
        self.measMain.setText(_translate("CompactMainWindow", "Measurement Main Window", None))
        self.menuDatei.setTitle(_translate("CompactMainWindow", "Datei", None))
        self.menuSettings.setTitle(_translate("CompactMainWindow", "Settings", None))
        self.actionSettings.setText(_translate("CompactMainWindow", "Settings", None))

