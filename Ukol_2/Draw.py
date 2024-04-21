from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from Building import Building

class Draw(QWidget):
    """ A class for drawing
    
    Attributes
    ----------
    self.buildings : list[Building]
        list of buildings (object of class Building)
    self.jarvis : bool
        True - Jarvis Scan is used for construct Convex Hull, otherwise Graham Scan is used
    
    Methods
    -------
    getBuldings():
        Getter for access to buildings
        
    setJarvis():
        Choose method for creating Convex Hull
    
    cleardata()
        Clear all data
        
    clearmbr()
        Clear generalized buildings
    
    paintEvent(e)
        Method for redraw Canvas
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.buildings : list[Building] = []
        self.jarvis = True
       
    def getBuildings(self) -> list:
        """Getter for access to buildings"""
        return self.buildings
       
    def setJarvis(self):
        """Choose method for creating Convex Hull"""
        self.jarvis = not self.jarvis
         
    def clearData(self):
        """Clear all data"""
        
        self.buildings = []
        
    def clearmbr(self):
        """Clear generalized buildings"""
        
        for building in self.buildings:
            building.setBuildingGeneralize(QPolygonF())
        
    def paintEvent(self, e: QPaintEvent):
        """Method for redraw Canvas"""
        
        # Create new draphic object
        qp = QPainter(self)
        
        # Start drawing
        qp.begin(self)        
        
        # Draw buildings
        for building in self.buildings:
            # Set Atributes for bulding
            qp.setPen(Qt.GlobalColor.black)
            qp.setBrush(Qt.GlobalColor.yellow)
            
            qp.drawPolygon(building.getBuilding())
    
            # Set atributes for simplified building
            qp.setPen(Qt.GlobalColor.red)
            qp.setBrush(Qt.GlobalColor.transparent)
            
            # Draw generalized building
            qp.drawPolygon(building.getBuildingGeneralize())
        
        # End drawing
        qp.end()