from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QPolygonF

class QPoint3DF(QPointF):
    def __init__(self, x : float, y : float, z : float) -> None:
        super().__init__(x, y)
        self.z = z
        
    def __eq__(self, other) -> bool:
        return (self.x() == other.x() and self.y() == other.y() and self.getZ() == other.getZ())
    
    def getZ(self):
        return self.z