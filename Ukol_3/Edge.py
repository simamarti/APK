from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QPolygonF
from QPoint3DF import QPoint3DF

class Edge:
    def __init__(self, start : QPoint3DF, end : QPoint3DF) -> None:
        self.start = start
        self.end = end
    
    def getStart(self):
        return self.start
    
    def getEnd(self):
        return self.end
    
    def switchOriantation(self):
        return Edge(self.end, self.start)
    
    def __eq__(self, other) -> bool:
        
        return (self.start == other.start and self.end == other.end)
    
    