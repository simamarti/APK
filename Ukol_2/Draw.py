from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from Building import Building

class Draw(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # self.building = QPolygonF([QPointF(100, 100), QPointF(100, 200), QPointF(200, 150)])
        self.buildings : list[Building] = []
        
        self.jarvis = True
       
    def getBuildings(self) -> list:
        # Return analyzed polygons
        return self.buildings
       
    def setJarvis(self):
        self.jarvis = not self.jarvis
         
    def clearData(self):
        # Clear building
        self.building = []
        
        # Clear MBR
        self.mbr = []
        
        # Repaint screen
        self.repaint()
        
    def clearmbr(self):
        self.mbr = []
        
    def mousePressEvent(self, e: QMouseEvent):
        
        # Get Coordinates
        x = e.position().x()
        y = e.position().y()
        
        # Add vertex
        p = QPointF(x, y)
                
        # Add point to building
        # self.building.append(p)
        
        # Repaint screen
        self.repaint()
        
    def paintEvent(self, e: QPaintEvent):
        # Draw situation
        
        # Create new draphic object
        qp = QPainter(self)
        
        # Start drawing
        qp.begin(self)
        
        # Set atributes
        
        
        # Draw buildings
        for building in self.buildings:
            # Set Atributes for bulding
            qp.setPen(Qt.GlobalColor.black)
            qp.setBrush(Qt.GlobalColor.yellow)
            
            qp.drawPolygon(building.building)
    
            # Set atributes for simplified building
            qp.setPen(Qt.GlobalColor.red)
            qp.setBrush(Qt.GlobalColor.transparent)
            
            # Draw mbr
            qp.drawPolygon(building.building_generalize)
        
        # Set graphical atributes for ch
        # Draw convex hull
        
        # Set graphical atributes for mbr
        # Draw Minimum bounding rectangle
        
        # Draw point
        qp.setPen(Qt.GlobalColor.blue)
        
        # End drawing
        qp.end()