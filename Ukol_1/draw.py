from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *

class Draw(QWidget):
    """ A class for drawing on the Canvas
    
    Attributes
    ----------
    None
    
    Methods
    -------
    switchDrawing(q, polygons):
        Decide what will be drawn (point/polygon) 

    getPols():
        get polygons
    
    getQ():
        get point
        
    clearData()
        clear Canvas
        
    mousePressEvent(e)
        acions when mouse is pressed
        
    paintEvent(e)
        draw on the Canvas
    """
    def __init__(self, *args, **kwargs):
        """Constructor for class Draw"""
        super().__init__(*args, **kwargs)
        self.q = QPointF(-100, -100)
        self.pols = []
        self.scale = 1
        self.dx = 0
        self.dy = 0
        self.addVertex = False
        self.intersect = []
    
    def switchDrawing(self):
        """Decide what will be drawn (point/polygon)"""
        # Change what will be drawn
        self.addVertex = not self.addVertex
    
    def getPols(self) -> list:
        """Get polygons"""        
        # Return polygons
        return self.pols
    
    def getQ(self) -> QPointF:
        """Get point"""        
        # Return analyzed point
        return self.q
    
    def clearData(self):
        """Clear Canvas"""        
        self.pols = []
        self.intersect = []
        self.q.setX(-100)
        self.q.setY(-100)
        
        # Repaint screen
        self.repaint()
            
    def mousePressEvent(self, e: QMouseEvent):
        """Actions when mouse is presssed 
            
        Parameters
        ----------
        e : QMouseEvent
            class decribes mouse event
        """
        # Get Coordinates
        x = e.position().x()
        y = e.position().y()
        
        if self.addVertex: 
            # ZRUŠIT??
            # Create temporary point
            p = QPointF(x, y)
                
            # Add point to polygon
            self.pols.append(p)
            print(len(self.pols))
        else:
            self.q.setX(x)
            self.q.setY(y)
            
            self.intersect = []
        
        # Repaint screen
        self.repaint()
        
    def paintEvent(self, e: QPaintEvent):
        """Draw on the Canvas
            
        Parameters
        ----------
        e : QMPaintEvent
            class decribes paint event
        """
        # Draw situation
        
        # Create new draphic object
        qp = QPainter(self)
        
        # Start drawing
        qp.begin(self)
        
        # Set atributes
        qp.setPen(Qt.GlobalColor.red)
        qp.setBrush(Qt.GlobalColor.green)
        
        # Draw polygons
        for polygon in self.pols:
            qp.drawPolygon(polygon.verticies)
        
        if self.intersect:
            qp.setPen(Qt.GlobalColor.yellow)
            qp.setBrush(Qt.GlobalColor.red)
            for pol in self.intersect:
                qp.drawPolygon(self.pols[pol].verticies)
        elif self.pols:
            qp.setPen(Qt.GlobalColor.red)
            qp.setBrush(Qt.GlobalColor.green)
            for pol in self.intersect:
                qp.drawPolygon(self.pols[pol].verticies)
            
        # Draw point
        qp.setPen(Qt.GlobalColor.yellow)
        qp.setBrush(Qt.GlobalColor.yellow)
        r = 5
        qp.drawEllipse(int(self.q.x() - r), int(self.q.y() - r), 2*r, 2*r)
        # End drawing
        qp.end()