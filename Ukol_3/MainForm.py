# Form implementation generated from reading ui file 'MainForm.ui'
#
# Created by: PyQt6 UI code generator 6.6.1
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt6 import QtCore, QtGui, QtWidgets
from Algorithms import Algorithms
from Settings import *
from load import loadPoints


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
        self.menuAnalysis = QtWidgets.QMenu(parent=self.menubar)
        self.menuAnalysis.setObjectName("menuAnalysis")
        self.menuView = QtWidgets.QMenu(parent=self.menubar)
        self.menuView.setObjectName("menuView")
        self.menuClear = QtWidgets.QMenu(parent=self.menubar)
        self.menuClear.setObjectName("menuClear")
        self.menuSettings = QtWidgets.QMenu(parent=self.menubar)
        self.menuSettings.setObjectName("menuSettings")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(parent=MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(parent=MainWindow)
        self.toolBar.setObjectName("toolBar")
        MainWindow.addToolBar(QtCore.Qt.ToolBarArea.TopToolBarArea, self.toolBar)
        self.actionOpen = QtGui.QAction(parent=MainWindow)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("./images/icons/open_file.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionOpen.setIcon(icon)
        self.actionOpen.setObjectName("actionOpen")
        self.actionClose = QtGui.QAction(parent=MainWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("./images/icons/exit.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionClose.setIcon(icon1)
        self.actionClose.setObjectName("actionClose")
        self.actionCreate_DTM = QtGui.QAction(parent=MainWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("./images/icons/triangles2.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionCreate_DTM.setIcon(icon2)
        self.actionCreate_DTM.setObjectName("actionCreate_DTM")
        self.actionCreate_contour_lines = QtGui.QAction(parent=MainWindow)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap("./images/icons/contours2.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionCreate_contour_lines.setIcon(icon3)
        self.actionCreate_contour_lines.setObjectName("actionCreate_contour_lines")
        self.actionAnalyze_slope = QtGui.QAction(parent=MainWindow)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap("./images/icons/slope2.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionAnalyze_slope.setIcon(icon4)
        self.actionAnalyze_slope.setObjectName("actionAnalyze_slope")
        self.actionAnalyze_exposition = QtGui.QAction(parent=MainWindow)
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap("./images/icons/orientation2.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionAnalyze_exposition.setIcon(icon5)
        self.actionAnalyze_exposition.setObjectName("actionAnalyze_exposition")
        self.actionDTM = QtGui.QAction(parent=MainWindow)
        self.actionDTM.setCheckable(True)
        self.actionDTM.setObjectName("actionDTM")
        self.actionContour_line = QtGui.QAction(parent=MainWindow)
        self.actionContour_line.setCheckable(True)
        self.actionContour_line.setObjectName("actionContour_line")
        self.actionSlope = QtGui.QAction(parent=MainWindow)
        self.actionSlope.setCheckable(True)
        self.actionSlope.setObjectName("actionSlope")
        self.actionExposition = QtGui.QAction(parent=MainWindow)
        self.actionExposition.setCheckable(True)
        self.actionExposition.setObjectName("actionExposition")
        self.actionClear_results = QtGui.QAction(parent=MainWindow)
        icon6 = QtGui.QIcon()
        icon6.addPixmap(QtGui.QPixmap("./images/icons/clear.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionClear_results.setIcon(icon6)
        self.actionClear_results.setObjectName("actionClear_results")
        self.actionClear_all = QtGui.QAction(parent=MainWindow)
        icon7 = QtGui.QIcon()
        icon7.addPixmap(QtGui.QPixmap("./images/icons/clear_all.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionClear_all.setIcon(icon7)
        self.actionClear_all.setObjectName("actionClear_all")
        self.actionParameters = QtGui.QAction(parent=MainWindow)
        icon8 = QtGui.QIcon()
        icon8.addPixmap(QtGui.QPixmap("./images/icons/settings.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionParameters.setIcon(icon8)
        self.actionParameters.setObjectName("actionParameters")
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionClose)
        self.menuAnalysis.addAction(self.actionCreate_DTM)
        self.menuAnalysis.addSeparator()
        self.menuAnalysis.addAction(self.actionCreate_contour_lines)
        self.menuAnalysis.addAction(self.actionAnalyze_slope)
        self.menuAnalysis.addAction(self.actionAnalyze_exposition)
        self.menuView.addAction(self.actionDTM)
        self.menuView.addAction(self.actionContour_line)
        self.menuView.addAction(self.actionSlope)
        self.menuView.addAction(self.actionExposition)
        self.menuClear.addAction(self.actionClear_results)
        self.menuClear.addAction(self.actionClear_all)
        self.menuSettings.addAction(self.actionParameters)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuAnalysis.menuAction())
        self.menubar.addAction(self.menuView.menuAction())
        self.menubar.addAction(self.menuClear.menuAction())
        self.menubar.addAction(self.menuSettings.menuAction())
        self.toolBar.addAction(self.actionOpen)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionCreate_DTM)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionCreate_contour_lines)
        self.toolBar.addAction(self.actionAnalyze_slope)
        self.toolBar.addAction(self.actionAnalyze_exposition)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionClear_results)
        self.toolBar.addAction(self.actionClear_all)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionParameters)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionClose)

        # Settings
        self.settings = QtWidgets.QDialog()
        self.ui = Ui_Settings()
        self.ui.setupUi(self.settings)

        self.retranslateUi(MainWindow)
        self.actionOpen.triggered.connect(self.openClick) # type: ignore
        self.actionCreate_DTM.triggered.connect(self.createDTMClick) # type: ignore
        self.actionAnalyze_slope.triggered.connect(self.analyzeSlopeClick) # type: ignore
        self.actionAnalyze_exposition.triggered.connect(self.analyzeAspectClick) # type: ignore
        self.actionClear_results.triggered.connect(self.clearClick) # type: ignore
        self.actionClear_all.triggered.connect(self.clearAllClick) # type: ignore
        self.actionDTM.triggered.connect(self.viewDTMClick) # type: ignore
        self.actionExposition.triggered.connect(self.viewExpositionClick) # type: ignore
        self.actionContour_line.triggered.connect(self.viewContourLineClick)
        self.actionSlope.triggered.connect(self.viewSlopeClick)
        self.actionCreate_contour_lines.triggered.connect(self.createContourLinesClick)


        self.actionClose.triggered.connect(MainWindow.close) # type: ignore

        self.actionParameters.triggered.connect(self.setParameters) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def openClick(self):
        """Open file dialog and choose *.JSON file"""
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(caption="Open File", directory="Data/.", filter="JSON file (*.json; *.geojson)")
        if filename == "":
            return None

        self.Canvas.clearAll()
        self.Canvas.points, self.Canvas.border = loadPoints(filename)

    def createDTMClick(self):
        """Create Delaunay triangulation on click"""

        # Get points
        points = self.Canvas.getPoints()

        # Points has not been loaded yet
        if len(points) == 0:
            return None

        # Create triangulation
        alg = Algorithms()
        dt = alg.createDT(points)

        # Clip dt by border
        dt = alg.clipDt(dt, self.Canvas.getBorder())

        # Upadate DT
        self.Canvas.setDT(dt)

        # Check menu item
        self.actionDTM.setChecked(True)
        
        # Repaint screen
        self.Canvas.repaint()


    def createContourLinesClick(self):
        """Create Contour Lines on click"""

        # Get Delaunay triangulation
        alg = Algorithms()

        # Do we have a triangulation
        dt = self.Canvas.getDT()

        # No triangulation constructed
        if not dt:
            # Get points
            points = self.Canvas.getPoints()

            # Points has not been loaded yet
            if len(points) == 0:
                return None

            # Create DT
            dt = alg.createDT(points)

            # Clip dt by border
            dt = alg.clipDt(dt, self.Canvas.getBorder())

            # Set results
            self.Canvas.setDT(dt)

        # Get contour lines parameters from dialog window
        zmin = float(self.ui.minVal.text())
        zmax = float(self.ui.maxVal.text())
        dz = float(self.ui.intVal.text())

        # Create contour lines
        contours = alg.createContourLines(dz, zmin, zmax, self.Canvas.getDT())

        # Set results
        self.Canvas.setContours(contours)
        
        # Check menu item
        self.actionContour_line.setChecked(True)
        
        # Repaint screen
        self.Canvas.repaint()

    def analyzeSlopeClick(self):
        """Analyze slope on click"""

        # Get Delaunay triangulation
        alg = Algorithms()

        # Do we have a triangulation
        dt = self.Canvas.getDT()

        # No triangulation constructed
        if not dt:
            # Get points
            points = self.Canvas.getPoints()

            # Points has not been loaded yet
            if len(points) == 0:
                return None

            # Create DT
            dt = alg.createDT(points)

            # Clip dt by border
            dt = alg.clipDt(dt, self.Canvas.getBorder())

            # Set results
            self.Canvas.setDT(dt)

        # Analyze dtm slope
        slope = alg.analyzeDTMSlope(self.Canvas.getDT())

        # Set results
        self.Canvas.setDTMSlope(slope)

        # Check menu item
        self.actionSlope.setChecked(True)
        
        # Repaint screen
        self.Canvas.repaint()
        
    def analyzeAspectClick(self):
        """Create Aspect on click"""

        # Get Delaunay triangulation
        alg = Algorithms()

        # Do we have a triangulation
        dt = self.Canvas.getDT()

        # No triangulation constructed
        if not dt:
            # Get points
            points = self.Canvas.getPoints()

            # Points has not been loaded yet
            if len(points) == 0:
                return None

            # Create DT
            dt = alg.createDT(points)

            # Clip dt by border
            dt = alg.clipDt(dt, self.Canvas.getBorder())

            # Set results
            self.Canvas.setDT(dt)

        # Analyze dtm aspect
        aspect = alg.analyzeDTMAspect(self.Canvas.getDT())

        # Set results
        self.Canvas.setDTMAspect(aspect)

        # Check menu item
        self.actionExposition.setChecked(True)

        # Repaint screen
        self.Canvas.repaint()
        
    def clearClick(self):
        """Clear analysis"""

        # Clear results
        self.Canvas.clearResults()

        # Repaint screen
        self.Canvas.repaint()

    def clearAllClick(self):
        """Clear all"""

        # Clear all data
        self.Canvas.clearAll()

        # Repaint screeen
        self.Canvas.repaint()

    def viewDTMClick(self):
        """Change visibility of Delaunay triangulation"""

        # Enable/disable DTM
        self.Canvas.setViewDT(self.actionDTM.isChecked())

        # Upadate
        self.Canvas.update()

    def viewContourLineClick(self):
        """Change visibility of Contour Lines"""

        # Enable/disable Contour lines
        self.Canvas.setViewContourLine(self.actionContour_line.isChecked())

        # Upadate
        self.Canvas.update()

    def viewSlopeClick(self):
        """Change visibility of Slope"""

        # Enable/disable Slope
        self.Canvas.setViewSlope(self.actionSlope.isChecked())

        # Upadate
        self.Canvas.update()

    def viewExpositionClick(self):
        """Change visibility of Exposition"""

        # Enable/disable Exposition
        self.Canvas.setViewAspect(self.actionExposition.isChecked())

        # Upadate
        self.Canvas.update()

    def setParameters(self):
        """Set parameters for creating Contour Lines using Dialog window"""

        self.settings.show()

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "DTM Analysis"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuAnalysis.setTitle(_translate("MainWindow", "Analysis"))
        self.menuView.setTitle(_translate("MainWindow", "View"))
        self.menuClear.setTitle(_translate("MainWindow", "Clear"))
        self.toolBar.setWindowTitle(_translate("MainWindow", "toolBar"))
        self.actionOpen.setText(_translate("MainWindow", "Open"))
        self.actionClose.setText(_translate("MainWindow", "Close"))
        self.actionCreate_DTM.setText(_translate("MainWindow", "Create DTM"))
        self.actionCreate_contour_lines.setText(_translate("MainWindow", "Create contour lines"))
        self.actionAnalyze_slope.setText(_translate("MainWindow", "Analyze slope"))
        self.actionAnalyze_exposition.setText(_translate("MainWindow", "Analyze exposition"))
        self.actionDTM.setText(_translate("MainWindow", "DTM"))
        self.actionContour_line.setText(_translate("MainWindow", "Contour line"))
        self.actionSlope.setText(_translate("MainWindow", "Slope"))
        self.actionExposition.setText(_translate("MainWindow", "Exposition"))
        self.actionClear_results.setText(_translate("MainWindow", "Clear results"))
        self.actionClear_all.setText(_translate("MainWindow", "Clear all"))
from Draw import Draw 


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())