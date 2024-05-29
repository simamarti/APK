from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QPolygonF
from QPoint3DF import QPoint3DF

class Edge:
    """ A class to represent Edge
    
    Attributes
    ----------
    self.start : QPoint3DF
    self.end : QPoint3DF
    
    Methods
    -------
    getStart():
        Get starting point of the Edge
    getEnd():
        Get ending point of the Edge
        
    switchOriantation():
        Switch orientation of Edge
        
    __eq__(other)
        Overloaded equal operator
    """

    def __init__(self, start : QPoint3DF, end : QPoint3DF) -> None:
        self.start = start
        self.end = end

    def getStart(self):
        """Get starting point of the Edge"""

        return self.start

    def getEnd(self):
        """Get ending point of the Edge"""

        return self.end

    def switchOriantation(self):
        """Switch orientation of Edge"""

        return Edge(self.end, self.start)

    def __eq__(self, other) -> bool:
        """Overloaded equal operator"""

        return (self.start == other.start and self.end == other.end)