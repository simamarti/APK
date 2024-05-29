from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from Edge import Edge
from QPoint3DF import QPoint3DF
from Triangle import Triangle
from random import random
from math import pi

class Draw(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.points : list[QPoint3DF] = []
        self.border : list[QPoint3DF] = []
        self.dt : list[Edge] = []
        self.contours : list[Edge] = []
        self.dtmSlope : list[Triangle] = []
        self.dtmAspect : list[Triangle] = []
        self.viewDT = True
        self.viewContourLines = True
        self.viewSlope = True
        self.viewAspect = True
        
    def clearAll(self):
        # Clear points
        self.points.clear()
        
        # Clear results
        self.clearResults()
        
        # Repaint screen
        self.repaint()
        
    def clearResults(self):
        # Clear DT
        self.dt.clear()
        
        # Clear contours
        self.contours.clear()
        
        # Clear slope
        self.dtmSlope.clear()
        
        # Clear aspect
        self.dtmAspect.clear()
        
        # Repaint screen
        self.repaint()
        
    def getPoints(self):
        return self.points
    
    def getBorder(self):
        return self.border
    
    def setBorder(self, list : list[QPoint3DF]):
        self.points = list
    
    def getDT(self):
        # Return DT
        return self.dt
           
    def setDT(self, dt : list[Edge]) -> None:
        self.dt = dt
        
    def getDTMSlope(self):
        return self.dtmSlope
    
    def setDTMSlope(self, dtm_slope : list[Triangle]):
        self.dtmSlope = dtm_slope

    def getDTMAspect(self):
        return self.dtmAspect
    
    def setDTMAspect(self, dtm_aspect : list[Triangle]):
        self.dtmAspect = dtm_aspect

    def setContours(self, contours : list[Edge]) -> None:
        self.contours = contours
        
    def setViewDT(self, viewDT):
        self.viewDT = viewDT
        
    def setViewContourLine(self, viewContourLines):
        self.viewContourLines = viewContourLines
        
    def setViewSlope(self, viewSlope):
        self.viewSlope = viewSlope
        
    def setViewAspect(self, viewAspect):
        self.viewAspect = viewAspect
        
    def mousePressEvent(self, e: QMouseEvent):
        
        # Get Coordinates
        x = e.position().x()
        y = e.position().y()

        # Generate random height
        zmin = 150
        zmax = 1600
        z = random() * (zmax - zmin) + zmin

        # Add new point
        p = QPoint3DF(x, y, z)
                
        # Add point to building
        self.points.append(p)
        
        # Repaint screen
        self.repaint()
        
    def paintEvent(self, e: QPaintEvent):
        # Draw situation

        # Create new draphic object
        qp = QPainter(self)
        
        # Start drawing
        qp.begin(self)
    
        if self.viewAspect:
            # Draw aspect
            qp.setPen(Qt.GlobalColor.gray)
            for t in self.dtmAspect:
                
                # Map Aspect from <-pi; pi> to <0; 2pi>
                aspect = t.getAspect() + pi
                
                # Convert radians to degrees
                aspect = aspect*180/pi
                
                # Corvert Aspect to color using HSV color model
                color = QColor.fromHsv((int(aspect) - int(aspect)%22 + 11)%359, 255, 255, 255)
                qp.setBrush(color)

                # Draw triangle
                qp.drawPolygon(t.getVertices())
                
        if self.viewSlope:
            # Draw slope
            qp.setPen(Qt.GlobalColor.gray)
            for t in self.dtmSlope:
                slope = t.getSlope()
                
                # Convert slope to color
                mju = 510/pi
                col = int(255 - mju*slope)
                color = QColor(col, col, col)
                qp.setBrush(color)

                # Draw triangle
                qp.drawPolygon(t.getVertices())
        
        # Set atributes
        qp.setBrush(Qt.GlobalColor.transparent)
        
        if self.viewContourLines:
            qp.setBrush(Qt.GlobalColor.green)
            # Draw contour lines
            for e in self.contours:
                qp.drawLine(int(e.getStart().x()), int(e.getStart().y()), int(e.getEnd().x()), int(e.getEnd().y()))
                
        if self.viewDT:
            # Set atributes
            qp.setPen(QColor(255, 165, 0))
            
            # Draw triangulation
            for edge in self.dt:
                qp.drawLine(int(edge.getStart().x()), int(edge.getStart().y()), int(edge.getEnd().x()), int(edge.getEnd().y()))
        
        
        # Draw points
        # Set atributes
        qp.setPen(Qt.GlobalColor.black)
        qp.setBrush(Qt.GlobalColor.yellow)
        
        r = 2
        for p in self.points:
            qp.drawEllipse(int(p.x()) - r, int(p.y()) - r, 2*r, 2*r)
        
        # End drawing
        qp.end()