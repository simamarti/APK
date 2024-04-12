# Form implementation generated from reading ui file 'form.ui'
#
# Created by: PyQt6 UI code generator 6.6.1
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt6 import QtCore, QtGui, QtWidgets
from Algorithms import Algorithms
from Load import load_buildings
from math import inf
# C:\Users\simam\Onedrive\Dokumenty\GitHub\APK\Ukol_2
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(parent=MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.Canvas = Draw(parent=self.centralwidget)
        self.Canvas.setObjectName("Canvas")
        self.verticalLayout.addWidget(self.Canvas)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(parent=MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 26))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(parent=self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuSimplify = QtWidgets.QMenu(parent=self.menubar)
        self.menuSimplify.setObjectName("menuSimplify")
        self.menuView = QtWidgets.QMenu(parent=self.menubar)
        self.menuView.setObjectName("menuView")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(parent=MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(parent=MainWindow)
        self.toolBar.setObjectName("toolBar")
        MainWindow.addToolBar(QtCore.Qt.ToolBarArea.TopToolBarArea, self.toolBar)
        self.actionOpen = QtGui.QAction(parent=MainWindow)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("Images/open_file.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionOpen.setIcon(icon)
        self.actionOpen.setObjectName("actionOpen")
        self.actionClose_2 = QtGui.QAction(parent=MainWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("Images/exit.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionClose_2.setIcon(icon1)
        self.actionClose_2.setObjectName("actionClose_2")
        self.actionMinimum_Area_Enclosing_Rectangle = QtGui.QAction(parent=MainWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("Images/maer.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionMinimum_Area_Enclosing_Rectangle.setIcon(icon2)
        self.actionMinimum_Area_Enclosing_Rectangle.setObjectName("actionMinimum_Area_Enclosing_Rectangle")
        self.actionPCA = QtGui.QAction(parent=MainWindow)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap("Images/pca.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionPCA.setIcon(icon3)
        self.actionPCA.setObjectName("actionPCA")
        self.actionClear_results = QtGui.QAction(parent=MainWindow)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap("Images/clear_ch.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionClear_results.setIcon(icon4)
        self.actionClear_results.setObjectName("actionClear_results")
        self.actionClear_all = QtGui.QAction(parent=MainWindow)
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap("Images/clear_er.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionClear_all.setIcon(icon5)
        self.actionClear_all.setObjectName("actionClear_all")
        self.actionLongest_Edge = QtGui.QAction(parent=MainWindow)
        icon6 = QtGui.QIcon()
        icon6.addPixmap(QtGui.QPixmap("Images/longestedge.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionLongest_Edge.setIcon(icon6)
        self.actionLongest_Edge.setObjectName("actionLongest_Edge")
        self.actionWall_Average = QtGui.QAction(parent=MainWindow)
        icon7 = QtGui.QIcon()
        icon7.addPixmap(QtGui.QPixmap("Images/wa.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionWall_Average.setIcon(icon7)
        self.actionWall_Average.setObjectName("actionWall_Average")
        self.actionWeighted_Bisector = QtGui.QAction(parent=MainWindow)
        icon8 = QtGui.QIcon()
        icon8.addPixmap(QtGui.QPixmap("Images/weightedbisector.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionWeighted_Bisector.setIcon(icon8)
        self.actionWeighted_Bisector.setObjectName("actionWeighted_Bisector")
        self.actionJarvis_Scan = QtGui.QAction(parent=MainWindow)
        self.actionJarvis_Scan.setCheckable(True)
        self.actionJarvis_Scan.setChecked(True)
        icon9 = QtGui.QIcon()
        icon9.addPixmap(QtGui.QPixmap("Images/Jarvis_scan.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionJarvis_Scan.setIcon(icon9)
        self.actionJarvis_Scan.setObjectName("actionJarvis_Scan")
        self.actionGrahamScan = QtGui.QAction(parent=MainWindow)
        self.actionGrahamScan.setCheckable(True)
        icon10 = QtGui.QIcon()
        icon10.addPixmap(QtGui.QPixmap("Images/Lída_scan.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionGrahamScan.setIcon(icon10)
        self.actionGrahamScan.setObjectName("actionL_da_ob_lka")
        self.actionValidation = QtGui.QAction(parent=MainWindow)
        icon11 = QtGui.QIcon()
        icon11.addPixmap(QtGui.QPixmap("Images/Validate.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.actionValidation.setIcon(icon11)
        self.actionValidation.setObjectName("actionValidation")
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionClose_2)
        self.menuSimplify.addAction(self.actionMinimum_Area_Enclosing_Rectangle)
        self.menuSimplify.addAction(self.actionPCA)
        self.menuSimplify.addAction(self.actionLongest_Edge)
        self.menuSimplify.addAction(self.actionWall_Average)
        self.menuSimplify.addAction(self.actionWeighted_Bisector)
        self.menuSimplify.addSeparator()
        self.menuSimplify.addAction(self.actionJarvis_Scan)
        self.menuSimplify.addAction(self.actionGrahamScan)
        self.menuView.addAction(self.actionValidation)
        self.menuView.addSeparator()
        self.menuView.addAction(self.actionClear_results)
        self.menuView.addSeparator()
        self.menuView.addAction(self.actionClear_all)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuSimplify.menuAction())
        self.menubar.addAction(self.menuView.menuAction())
        self.toolBar.addAction(self.actionOpen)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionLongest_Edge)
        self.toolBar.addAction(self.actionMinimum_Area_Enclosing_Rectangle)
        self.toolBar.addAction(self.actionPCA)
        self.toolBar.addAction(self.actionWall_Average)
        self.toolBar.addAction(self.actionWeighted_Bisector)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionJarvis_Scan)
        self.toolBar.addAction(self.actionGrahamScan)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionValidation)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionClear_results)
        self.toolBar.addAction(self.actionClear_all)
        self.toolBar.addAction(self.actionClose_2)

        self.retranslateUi(MainWindow)
        self.actionOpen.triggered.connect(self.openClick) # type: ignore
        self.actionMinimum_Area_Enclosing_Rectangle.triggered.connect(self.mbrClick) # type: ignore
        self.actionClear_all.triggered.connect(self.clearAllClick) # type: ignore
        self.actionPCA.triggered.connect(self.pcaClick) # type: ignore
        self.actionClear_results.triggered.connect(self.clearClick) # type: ignore
        self.actionClose_2.triggered.connect(MainWindow.close) # type: ignore
        self.actionJarvis_Scan.triggered.connect(self.jarvisScan) # type: ignore
        self.actionGrahamScan.triggered.connect(self.grahamScan) # type: ignore
        self.actionLongest_Edge.triggered.connect(self.longestEdge) # type: ignore
        self.actionWall_Average.triggered.connect(self.wallAverage) # type: ignore
        self.actionWeighted_Bisector.triggered.connect(self.weightedBisector) # type: ignore
        self.actionValidation.triggered.connect(self.validation) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
    
    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Simplify buildings"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuSimplify.setTitle(_translate("MainWindow", "Simplify"))
        self.menuView.setTitle(_translate("MainWindow", "View"))
        self.toolBar.setWindowTitle(_translate("MainWindow", "toolBar"))
        self.actionOpen.setText(_translate("MainWindow", "Open file"))
        self.actionClose_2.setText(_translate("MainWindow", "Exit"))
        self.actionMinimum_Area_Enclosing_Rectangle.setText(_translate("MainWindow", "Minimum Area Enclosing Rectangle"))
        self.actionPCA.setText(_translate("MainWindow", "PCA"))
        self.actionClear_results.setText(_translate("MainWindow", "Clear results"))
        self.actionClear_all.setText(_translate("MainWindow", "Clear all"))
        self.actionLongest_Edge.setText(_translate("MainWindow", "Longest Edge"))
        self.actionWall_Average.setText(_translate("MainWindow", "Wall Average"))
        self.actionWeighted_Bisector.setText(_translate("MainWindow", "Weighted Bisector"))
        self.actionJarvis_Scan.setText(_translate("MainWindow", "Jarvis Scan"))
        self.actionGrahamScan.setText(_translate("MainWindow", "Graham Scan"))
        self.actionValidation.setText(_translate("MainWindow", "Validation"))
        
    def openClick(self):
        """Open file dialog and choose *.JSON file"""
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(caption="Open File", directory="Data/.", filter="JSON file (*.json; *.geojson)")
        if filename == "":
            return None
        
        self.Canvas.clearData()
        print(f"Canvas height: {self.Canvas.height()}, canvas width: {self.Canvas.width()}")
        self.Canvas.buildings = load_buildings(filename)
        
    def mbrClick(self):
        
        a = Algorithms()
        #Get building
        buildings = self.Canvas.getBuildings()
        
        for building in buildings:
            #Simplify building
            simplify = None
            if self.Canvas.jarvis:
                simplify = a.createMBR(building.building, a.jarvisScan)
            else:
                simplify = a.createMBR(building.building, a.grahamScan)
            
            building.setBuildingGeneralize(simplify)     
        
        #Repaint screen
        self.Canvas.repaint()
        
    def pcaClick(self):
        
        a = Algorithms()
        
        #Get building
        buildings = self.Canvas.getBuildings()
        
        #Simplify building
        
        for building in buildings:
            #Simplify building
            pca = a.createERPCA(building)
            
            #Update Generalize Building
            building.setBuildingGeneralize(pca)    
            
        #Repaint screen
        self.Canvas.repaint()
    
    def clearClick(self):
        self.Canvas.clearmbr()
    
        #Repaint screen
        self.Canvas.repaint()
        
    def clearAllClick(self):
        self.Canvas.clearData()
        
        #Repaint screen
        self.Canvas.repaint()
    
    def jarvisScan(self):
        self.Canvas.jarvis = True
        self.actionGrahamScan.setChecked(False)
    
    def grahamScan(self):
        self.Canvas.jarvis = False
        self.actionJarvis_Scan.setChecked(False)
    
    def longestEdge(self):
        
        a = Algorithms()

        #Get building
        buildings = self.Canvas.getBuildings()
        
        for building in buildings:
            #Simplify building
            simplify = a.longestEdge(building.building)      
                  
            #Update Generalize Building
            building.building_generalize = simplify
        
        #Repaint screen
        self.Canvas.repaint()
    
    def wallAverage(self):
        
        a = Algorithms()
        
        #Get building
        buildings = self.Canvas.getBuildings()
        
        for building in buildings:
            #Simplify building
            simplify = a.wallAverage(building.building)      
            
            #Update Generalize Building
            building.building_generalize = simplify
        
        #Repaint screen
        self.Canvas.repaint()
    
    def weightedBisector(self):
        
        a = Algorithms()
        
        #Get buildings
        buildings = self.Canvas.getBuildings()
        
        for building in buildings:
            #Simplify building
            maer = a.weightedBisector(building.building)
            
            #Update Generalize Building
            building.setBuildingGeneralize(maer)    
        
        #Repaint screen
        self.Canvas.repaint()

    def validation(self):
        alg = Algorithms()
        
        text = alg.validation(self.Canvas.buildings)
        if text == "":
            return
        messagebox = QtWidgets.QMessageBox()
        messagebox.setWindowTitle("Validation")
        messagebox.setText(text)
            
        # Show analysis
        messagebox.exec()
        
    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Simplify buildings"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuSimplify.setTitle(_translate("MainWindow", "Simplify"))
        self.menuView.setTitle(_translate("MainWindow", "View"))
        self.toolBar.setWindowTitle(_translate("MainWindow", "toolBar"))
        self.actionOpen.setText(_translate("MainWindow", "Open file"))
        self.actionClose_2.setText(_translate("MainWindow", "Exit"))
        self.actionMinimum_Area_Enclosing_Rectangle.setText(_translate("MainWindow", "Minimum Area Enclosing Rectangle"))
        self.actionPCA.setText(_translate("MainWindow", "PCA"))
        self.actionClear_results.setText(_translate("MainWindow", "Clear results"))
        self.actionClear_all.setText(_translate("MainWindow", "Clear all"))
        self.actionLongest_Edge.setText(_translate("MainWindow", "Longest Edge"))
        self.actionWall_Average.setText(_translate("MainWindow", "Wall Average"))
        self.actionWeighted_Bisector.setText(_translate("MainWindow", "Weighted Bisector"))
        self.actionJarvis_Scan.setText(_translate("MainWindow", "Jarvis Scan"))
        self.actionGrahamScan.setText(_translate("MainWindow", "Graham Scan"))
from Draw import Draw


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())