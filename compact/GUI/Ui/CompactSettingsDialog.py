# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dialog.ui'
#
# Created: Fri Jun 26 15:11:46 2015
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

class Ui_CompactSettingsDialog(object):
    def setupUi(self, CompactSettingsDialog):
        CompactSettingsDialog.setObjectName(_fromUtf8("CompactSettingsDialog"))
        CompactSettingsDialog.resize(400, 300)
        self.buttonBox = QtGui.QDialogButtonBox(CompactSettingsDialog)
        self.buttonBox.setGeometry(QtCore.QRect(180, 250, 192, 32))
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.lineEdit = QtGui.QLineEdit(CompactSettingsDialog)
        self.lineEdit.setGeometry(QtCore.QRect(70, 150, 113, 21))
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))

        self.retranslateUi(CompactSettingsDialog)
        QtCore.QMetaObject.connectSlotsByName(CompactSettingsDialog)

    def retranslateUi(self, CompactSettingsDialog):
        CompactSettingsDialog.setWindowTitle(_translate("CompactSettingsDialog", "Dialog", None))
        self.lineEdit.setText(_translate("CompactSettingsDialog", "Settings", None))

