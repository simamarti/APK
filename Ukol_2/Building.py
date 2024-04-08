from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from math import inf

class Building(QWidget):
    
    def __init__(self) -> None:
        """Constructor for creating polygon"""
        self.building : QPolygonF =  QPolygonF()
        self.building_generalize : QPolygonF = QPolygonF()

    def __len__(self) -> int:
        """Overloaded build-in function len()"""
        return len(self.building)
    
    def addVertex(self, q : QPointF) -> None:
        """Add vertex to polygon
        
        Parameters
        ----------
        q : QPainterpath
            Point which will be added to polygon
        """
        self.building.append(q)
    
    def setBuildingGeneralize(self, pol : QPolygonF) -> None:
        self.building_generalize = pol