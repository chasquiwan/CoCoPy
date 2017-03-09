import sys
from PyQt4 import QtGui
from CompactMain import Compact

def main():
    app = QtGui.QApplication(sys.argv)
    com = Compact()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()