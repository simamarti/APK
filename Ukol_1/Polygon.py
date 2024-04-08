from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from math import inf

class Polygon(QWidget):
    """ A class to represent a polygon
    
    Attributes
    ----------
    None
    
    Methods
    -------
    __len__():
        Overloaded build-in function len()

    __getitem__(idx):
        Overloaded square brackets
        
    adVertex(q):
        Add vertex to polygon
        
    isInMMB(q):
        is point in minmax box of polygon?
    
    updateMMB(q)
        update Minmax box
        
    convertPolToPath()
        Convert QPolyfonF to QPainterPath
    """
    def __init__(self) -> None:
        """Constructor for creating polygon"""
        self.verticies : QPolygonF =  QPolygonF()
        self.path = QPainterPath()
        self.min_x : float = inf
        self.max_x : float = -inf
        self.min_y : float = inf
        self.max_y : float = -inf

    def __len__(self) -> int:
        """Overloaded build-in function len()"""
        return len(self.verticies)
    
    def __getitem__(self, idx : int) -> QPointF:
        """Overloaded square brackets"""
        return self.verticies[idx]
    
    def addVertex(self, q : QPointF) -> None:
        """Add vertex to polygon
        
        Parameters
        ----------
        q : QPainterpath
            Point which will be added to polygon
        """
        self.verticies.append(q)
        
    def isInMMB(self, q : QPointF) -> bool:
        """Decide if point is in minmax box of the polygon
        
        Parameters
        ----------
        q : QPointF
            Point which will be analyzed
        
        Returns
        -------
        True: Point is in the minmax box
        False: Poin is not in the minmax box
        """
        return (    self.min_x <= q.x() and self.max_x >= q.x() and 
                    self.min_y <= q.y() and self.max_y >= q.y())
    
    def updateMMB(self, q : QPointF) -> None:
        """Update minmax box
        
        Parameters
        ----------
        q : QPointF
            Point which will be analyzed
        """
        self.min_x = min(self.min_x, q.x())
        self.min_y = min(self.min_y, q.y())
        self.max_x = max(self.max_x, q.x())
        self.max_y = max(self.max_y, q.y())
        
    def convertPolToPath(self) -> None:
        """Convert QPolyfonF to QPainterPath
        
        Parameters
        ----------
        None
        """
        self.verticies.append(self.verticies[0])        # Connect last point with first point     
        self.path.addPolygon(self.verticies)
