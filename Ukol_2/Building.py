from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from math import inf

class Building(QWidget):
    """ A class to represent a polygon
    
    Attributes
    ----------
    self.buildings : QPolygonF
    self.bulding_generalize : QPolygonF
    
    Methods
    -------
    __len__():
        Overloaded build-in function len()

    adVertex(q):
        Add vertex to polygon
        
    getBuilding():
        Getter for building
        
    getBuildingGeneralize():
        Getter for generalized building
        
    setBuildingGeneralize(pol):
        Setter for update generalized building
    """
    def __init__(self) -> None:
        """Constructor for creating polygon"""
        self.building : QPolygonF =  QPolygonF()
        self.buildingGeneralize : QPolygonF = QPolygonF()

    def __len__(self) -> int:
        """Overloaded build-in function len()"""
        return len(self.building)
    
    def addVertex(self, q : QPointF) -> None:
        """Add vertex to polygon
        
        Parameters
        ----------
        q : QPointF
            Point which will be added to polygon
        """
        self.building.append(q)
    
    def getBuilding(self) -> QPolygonF:
        """Getter for building"""
        
        return self.building
    
    def getBuildingGeneralize(self) -> QPolygonF:
        """Getter for generalized building"""
        
        return self.buildingGeneralize
        
    def setBuildingGeneralize(self, pol : QPolygonF) -> None:
        """Setter for update generalized building
        
        Parameters
        ----------
        pol : QPolygonF
            Polygon of the generalized building
        """
        self.buildingGeneralize = pol