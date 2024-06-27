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
    """ A class to represent drawing
    
    Attributes
    ----------
    self.points : list[QPoint3DF]
    self.border : list[QPoint3DF]
    self.hole : list[QPoint3DF]
    self.dt : list[Edge]
    self.contours : list[Edge]
    self.dtmSlope : list[Triangle]
    self.dtmAspect : list[Triangle]
    self.createBorderSwitch : bool
    self.createHoleSwitch : bool
    self.viewDT : bool
    self.viewContourLines : bool
    self.viewSlope : bool
    self.viewAspect : bool
    
    Methods
    -------
    clearAll():
        Clear all things on the screen
    clearResults():
        Clear analysis on the screen
        
    getPoints():
        Get points
        
    getBorder()
        Get border points of region
    
    getHole()
        Get Hole
        
    setBorder(list : list[QPoin3DF])
        Set border points
    
    switchHole():
        switch between creating hole in points
    
    getSwtichHole():
        Get hole status
        
    switchBorder():
        switch between creating border of points
    
    getSwtichBorder():
        Get border status
        
    getDT()
        Get Delaunay Triangulation
        
    setDT(dt : list[Edge])
        Set Delaunay triangulation
    
    getDTMSlope()
        Get Slope of Delaunay Triangulation
        
    setDTMSlope(dtm_slope : list[Triangle])
        Set Slope of Delaunay triangulation
    
    getDTMAspect()
        Get Aspect of Delaunay Triangulation
        
    setDTMAspect(dtm_aspect : list[Triangle])
        Set Aspect of Delaunay triangulation
        
    setContours(contours : list[Edge])
        Set Contour Lines of Delaunay triangulation
        
    setViewDT(viewDT : bool)
        Set visibility of Delaunay Triangulation
        
    setViewContourLine(viewContourLine : bool)
        Set visibility of Contour Lines
        
    setViewSlope(slope : bool)
        Set visibility of Slope
        
    setViewAspect(aspect : bool)
        Set visibility of Aspect
        
    mousePressEvent(e: QMouseEvent)
        Method for mouse click
        
    paintEvent(e: QPaintEvent)
        Method for paint event
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.points : list[QPoint3DF] = []
        self.border : list[QPoint3DF] = []
        self.hole : list[QPoint3DF] = []
        self.dt : list[Edge] = []
        self.contours : list[Edge] = []
        self.dtmSlope : list[Triangle] = []
        self.dtmAspect : list[Triangle] = []
        self.createHoleSwitch : bool = False
        self.createBorderSwitch : bool = False
        self.viewDT = True
        self.viewContourLines = True
        self.viewSlope = True
        self.viewAspect = True

    def clearAll(self):
        """Clear all things on the screen"""

        # Clear points
        self.points.clear()
        
        # Clear border
        self.border.clear()
        
        # Clear hole
        self.hole.clear()

        # Clear results
        self.clearResults()

        # Repaint screen
        self.repaint()

    def clearResults(self):
        """Clear analysis on the screen"""

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
        """Get points"""

        return self.points

    def getBorder(self):
        """Get border points of region"""
        return self.border

    def getHole(self):
        """Get Hole"""
        
        return self.hole
        
    def setBorder(self, list : list[QPoint3DF]):
        """Set border points"""

        self.points = list

    def switchHole(self):
        """switch between creating hole in points"""

        self.createHoleSwitch = not self.createHoleSwitch
    
    def getSwitchHole(self):
        """Get Hole status"""
        return self.createHoleSwitch
    
    def switchBorder(self):
        """switch between creating border of points"""

        self.createBorderSwitch = not self.createBorderSwitch
    
    def getSwtichBorder(self):
        """Get Border status"""
        return self.createBorderSwitch
       
    def getDT(self):
        """Get Delaunay Triangulation"""

        # Return DT
        return self.dt

    def setDT(self, dt : list[Edge]) -> None:
        """Set Delaunay triangulation"""

        self.dt = dt

    def getDTMSlope(self):
        """Get Slope of Delaunay Triangulation"""

        return self.dtmSlope

    def setDTMSlope(self, dtm_slope : list[Triangle]):
        """Set Slope of Delaunay triangulation"""

        self.dtmSlope = dtm_slope

    def getDTMAspect(self):
        """Get Aspect of Delaunay Triangulation"""

        return self.dtmAspect

    def setDTMAspect(self, dtm_aspect : list[Triangle]):
        """Set Aspect of Delaunay triangulation"""

        self.dtmAspect = dtm_aspect

    def setContours(self, contours : list[Edge]) -> None:
        """Set Contour Lines of Delaunay triangulation"""

        self.contours = contours

    def setViewDT(self, viewDT):
        """Set visibility of Delaunay Triangulation"""

        self.viewDT = viewDT

    def setViewContourLine(self, viewContourLines):
        """Set visibility of Contour Lines"""

        self.viewContourLines = viewContourLines

    def setViewSlope(self, viewSlope):
        """Set visibility of Slope"""

        self.viewSlope = viewSlope

    def setViewAspect(self, viewAspect):
        """Set visibility of Aspect"""

        self.viewAspect = viewAspect

    def mousePressEvent(self, e: QMouseEvent):
        """Method for mouse click"""

        # Get Coordinates
        x = e.position().x()
        y = e.position().y()

        # Generate random height
        zmin = 150
        zmax = 1600
        z = random() * (zmax - zmin) + zmin
        
        # Add new point
        p = QPoint3DF(x, y, z)

        if self.getSwitchHole():
            
            # Add point to hole
            self.hole.append(p)
            
        elif self.getSwtichBorder():
            
            # Add point to border
            self.border.append(p)
        else:
            
            # Add point to DTM
            self.points.append(p)

        # Repaint screen
        self.repaint()

    def paintEvent(self, e: QPaintEvent):
        """Method for paint event"""

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
                hue = 0
                if aspect > 337.5 or aspect <= 22.5:    hue = 0
                elif aspect <= 67.5:                    hue = 45
                elif aspect <= 112.5:                   hue = 90
                elif aspect <= 157.5:                   hue = 135
                elif aspect <= 202.5:                   hue = 180
                elif aspect <= 247.5:                   hue = 225
                elif aspect <= 292.5:                   hue = 270
                elif aspect <= 337.5:                   hue = 315

                color = QColor.fromHsv(hue, 255, 255, 255)
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
            qp.setPen(QColor(139,69,19))
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


        # Draw hole
        # Set atributes
        qp.setPen(Qt.GlobalColor.red)
        qp.setBrush(QColor(240, 240, 240))

        n = len(self.border)
        for idx in range(n):
            p1 = self.border[idx]
            p2 = self.border[(idx + 1)%n]

            qp.drawLine(int(p1.x()), int(p1.y()), int(p2.x()), int(p2.y()))
            
        qp.setPen(Qt.GlobalColor.blue)

        n = len(self.hole)
        pol =QPolygonF()
        for idx in range(n):
            p = self.hole[idx]

            pol.append(QPointF(p.x(), p.y()))
        
        qp.drawPolygon(pol)
        # End drawing
        qp.end()