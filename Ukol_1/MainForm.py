# Form implementation generated from reading ui file 'First.ui'
#
# Created by: PyQt6 UI code generator 6.6.1
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt6 import QtCore, QtGui, QtWidgets
from Algorithms import Algorithms
from Load import load_polygons

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(parent=MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.Canvas = Draw(parent=self.centralwidget)
        self.Canvas.setObjectName("Canvas")
        self.horizontalLayout.addWidget(self.Canvas)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(parent=MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 26))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(parent=self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuInput = QtWidgets.QMenu(parent=self.menubar)
        self.menuInput.setObjectName("menuInput")
        self.menuAnalyze = QtWidgets.QMenu(parent=self.menubar)
        self.menuAnalyze.setObjectName("menuAnalyze")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(parent=MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(parent=MainWindow)
        self.toolBar.setObjectName("toolBar")
        MainWindow.addToolBar(QtCore.Qt.ToolBarArea.TopToolBarArea, self.toolBar)
        self.actionOpen = QtGui.QAction(parent=MainWindow)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("images/icons/open_file.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionOpen.setIcon(icon)
        self.actionOpen.setObjectName("actionOpen")
        self.actionExit = QtGui.QAction(parent=MainWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("images/icons/exit.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionExit.setIcon(icon1)
        self.actionExit.setObjectName("actionExit")
        self.actionPolygon = QtGui.QAction(parent=MainWindow)
        self.actionPolygon.setCheckable(True)
        self.actionPolygon.setObjectName("actionPolygon")
        self.actionDira = QtGui.QAction(parent=MainWindow)
        self.actionDira.setCheckable(True)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("images/icons/pointpol.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionDira.setIcon(icon2)
        self.actionDira.setObjectName("actionDira")
        self.actionClear = QtGui.QAction(parent=MainWindow)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap("images/icons/clear_all.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionClear.setIcon(icon3)
        self.actionClear.setObjectName("actionClear")
        self.actionRay_Crossing_Algorithm = QtGui.QAction(parent=MainWindow)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap("images/icons/ray.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionRay_Crossing_Algorithm.setIcon(icon4)
        self.actionRay_Crossing_Algorithm.setObjectName("actionRay_Crossing_Algorithm")
        self.actionWinding_number_Algorithm = QtGui.QAction(parent=MainWindow)
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap("images/icons/winding.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionWinding_number_Algorithm.setIcon(icon5)
        self.actionWinding_number_Algorithm.setObjectName("actionWinding_number_Algorithm")
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionExit)
        self.menuInput.addAction(self.actionDira)
        self.menuInput.addSeparator()
        self.menuInput.addAction(self.actionClear)
        self.menuAnalyze.addAction(self.actionRay_Crossing_Algorithm)
        self.menuAnalyze.addAction(self.actionWinding_number_Algorithm)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuInput.menuAction())
        self.menubar.addAction(self.menuAnalyze.menuAction())
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionOpen)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionDira)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionRay_Crossing_Algorithm)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionWinding_number_Algorithm)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionClear)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionExit)

        self.retranslateUi(MainWindow)
        self.actionOpen.triggered.connect(self.openClick) # type: ignore
        self.actionClear.triggered.connect(self.clear) # type: ignore
        self.actionWinding_number_Algorithm.triggered.connect(self.windingNumber) # type: ignore
        self.actionRay_Crossing_Algorithm.triggered.connect(self.rayCrossing) # type: ignore
        self.actionDira.triggered.connect(self.makeHole) # type: ignore
        self.actionExit.triggered.connect(MainWindow.close) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def openClick(self):
        """Open file dialog and choose *.JSON file"""
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(caption="Open File", directory="Data/.", filter="JSON file (*.json; *.geojson)")
        if filename == "":
            return None
        self.Canvas.pols = load_polygons(filename)
        
    def rayCrossing(self):
        """Ray Crossinng algorithm on the click"""
        # Get data
        q = self.Canvas.getQ()
        pol = self.Canvas.getPols()
        
        # Point was not choosen
        if q.x() <= 0 or q.y() <= 0 or not pol:
            return None
        
        # Run analysis
        messagebox = QtWidgets.QMessageBox()
        messagebox.setWindowTitle("Analyze point and polygon position")
        
        # Point inside/outside
        alg = Algorithms()
        self.Canvas.intersect = []
        # Choose polygons with point in their minmax box
        indexes = alg.preProcessPolygons(q, self.Canvas.pols)

        # Iterate throught choosen polygon
        for idx in indexes:
            if alg.rayCrossingAlgorithm(q, self.Canvas.pols[idx]):
                self.Canvas.intersect.append(idx)
        
        if self.Canvas.intersect:
            messagebox.setText("Inside")
        else:
            messagebox.setText("Outside")
        
        # Show analysis
        messagebox.exec()
        
    def windingNumber(self):
        """Winding number algorithm on the click"""
        """Ray Crossinng algorithm on the click"""
        # Get data
        q = self.Canvas.getQ()
        pol = self.Canvas.getPols()
        
        # Point was not choosen
        if q.x() <= 0 or q.y() <= 0 or not pol:
            return None
        
        # Run analysis
        messagebox = QtWidgets.QMessageBox()
        messagebox.setWindowTitle("Analyze point and polygon position")
        
        # Point inside/outside
        alg = Algorithms()
        self.Canvas.intersect = []
        # Choose polygons with point in their minmax box
        indexes = alg.preProcessPolygons(q, self.Canvas.pols)

        # Iterate throught choosen polygon
        for idx in indexes:
            if alg.windingNumber(q, self.Canvas.pols[idx]):
                self.Canvas.intersect.append(idx)
        
        if self.Canvas.intersect:
            messagebox.setText("Inside")
        else:
            messagebox.setText("Outside")
        
        # Show analysis
        messagebox.exec()
    
    def clear(self):
        """Clear Canvas"""
        self.Canvas.clearData()
        
    def makeHole(self):
        """Switch drawing between point and hole"""
        self.Canvas.switchDrawing()
        
    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Point and polygon position"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuInput.setTitle(_translate("MainWindow", "Input"))
        self.menuAnalyze.setTitle(_translate("MainWindow", "Analyze"))
        self.toolBar.setWindowTitle(_translate("MainWindow", "toolBar"))
        self.actionOpen.setText(_translate("MainWindow", "Open"))
        self.actionOpen.setToolTip(_translate("MainWindow", "Open file"))
        self.actionExit.setText(_translate("MainWindow", "Exit"))
        self.actionExit.setToolTip(_translate("MainWindow", "Close application"))
        self.actionPolygon.setText(_translate("MainWindow", "Polygon"))
        self.actionDira.setText(_translate("MainWindow", "Make Hole"))
        self.actionClear.setText(_translate("MainWindow", "Clear"))
        self.actionRay_Crossing_Algorithm.setText(_translate("MainWindow", "Ray Crossing Algorithm"))
        self.actionWinding_number_Algorithm.setText(_translate("MainWindow", "Winding number Algorithm"))
from draw import Draw


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())
