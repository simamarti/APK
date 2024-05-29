from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QPolygonF
from QPoint3DF import QPoint3DF

class Triangle:
    """ A class to represent Triangle
    
    Attributes
    ----------
    self.vertices : QPolygonF
    self.slope : float
    self.exposition : float
    
    Methods
    -------
    getVertices():
        Get list of verticies of Triangle
    getAspect():
        Get Aspect of Triangle
        
    getSlope():
        Get Slope of Triangle
    """

    def __init__(self, p1 : QPoint3DF, p2 : QPoint3DF, p3 : QPoint3DF, slope : float, exposition : float) -> None:

        self.vertices = QPolygonF()
        self.vertices.append(p1)
        self.vertices.append(p2)
        self.vertices.append(p3)

        self.slope = slope
        self.exposition = exposition

    def getVertices(self):
        """Get list of verticies of Triangle"""

        return self.vertices

    def getAspect(self):
        """Get Aspect of Triangle"""

        return self.exposition

    def getSlope(self):
        """Get Slope of Triangle"""

        return self.slope