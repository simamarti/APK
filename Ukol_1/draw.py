from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from Polygon import Polygon

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
    def __init__(self, *args, **kwargs) -> None:
        """Constructor for class Draw"""
        super().__init__(*args, **kwargs)
        self.q : QPointF = QPointF(-100, -100)
        self.pols : list[Polygon] = []
        self.hole : QPolygonF = QPolygonF()
        self.addHole : bool = False
        self.intersect : list[int] = []
    
    def switchDrawing(self) -> None:
        """Decide what will be drawn (point/hole)"""
        # Change what will be drawn
        self.addHole = not self.addHole
        
        # Color all to green (point is not in any polygons)
        self.intersect = []
        self.repaint()
    
    def getPols(self) -> list:
        """Get polygons"""        
        # Return polygons
        return self.pols
    
    def getQ(self) -> QPointF:
        """Get point"""        
        # Return analyzed point
        return self.q
    
    def clearData(self) -> None:
        """Clear Canvas"""        
        self.pols = []
        self.hole = QPolygonF()
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
        
        if self.addHole: 
            # „Delete“ point from Canvas
            self.q.setX(-100)
            self.q.setY(-100)
            
            # Create temporary point
            p = QPointF(x, y)
                
            # Add point to hole
            self.hole.append(p)
            
            # Create path from hole and subtracted from polygons
            if len(self.hole) > 1:
                path_hole = QPainterPath()
                path_hole.addPolygon(self.hole)
                for pol in self.pols:
                    pol.path = pol.path.subtracted(path_hole)
        else:               # or point will be inserted
            self.q.setX(x)
            self.q.setY(y)
            self.intersect = []         # Clear array with indexes of highlighted polygons

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
            qp.drawPath(polygon.path)

        if self.intersect:                          # Draw polygons with point
            qp.setPen(Qt.GlobalColor.yellow)
            qp.setBrush(Qt.GlobalColor.red)
            for pol in self.intersect:
                qp.drawPath(self.pols[pol].path)
        elif self.pols:                             # or redraw all polygon to green
            qp.setPen(Qt.GlobalColor.red)
            qp.setBrush(Qt.GlobalColor.green)
            for pol in self.intersect:
                qp.drawPath(self.pols[pol].path)
            
        # Draw point
        qp.setPen(Qt.GlobalColor.yellow)
        qp.setBrush(Qt.GlobalColor.yellow)
        r = 5
        qp.drawEllipse(int(self.q.x() - r), int(self.q.y() - r), 2*r, 2*r)
        # End drawing
        qp.end()