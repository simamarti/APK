from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from numpy.linalg import norm
from numpy import array, finfo, float64
from Polygon import Polygon
from math import acos, pi

class Algorithms:
    """ A class with algorithms
    
    Attributes
    ----------
    None
    
    Methods
    -------
    preProcessPolygons(q, polygons):
        Choose polygons in whitch q is in their minmax box 

    rayCrossingAlgorithm(q, pol):
        Ray Crossing algorithm
    
    windingNumber(q, pol):
        Winding Number algorithm
        
    angleTwoLine(q, startPoint, endPoint):
        Compute angle between two line
        
    pointEdgePosition(q, startPoint, endPoint):
        is line on the left of the line
    """
    def __init__(self):
        pass
    
    def preProcessPolygons(self, q : QPointF, polygons : list[Polygon]) -> list:
        """Choose polygons in whitch q is in their minmax box 
            
        Parameters
        ----------
        q : QPointF
            Analyzed point
        polygons : list [Polygon]
            list of all loaded polygons
            
        Returns
        -------
        inBB : list
            list of indexes of polygons whose minmax box contain point
        """
        inBB = []
        for idx, polygon in enumerate(polygons):
            if polygon.isInMMB(q):
                inBB.append(idx)
        return inBB
    
    def rayCrossingAlgorithm(self, q : QPointF, polygon : Polygon) -> int:
        """Ray Crossing algorithm 
            
        Parameters
        ----------
        q : QPointF
            Analyzed point
        pol : QPolygonF
            polygon
            
        Returns
        -------
        0: Point is in the polygon
        1: Point is outside the polygon
        """
        # Split QPainterPath to QPolygonF
        parts = polygon.path.toSubpathPolygons()
        
        # Inicialize amount of intersection
        k = 0
        eps = 1.0e-5        # epsilon for comparing floats
        
        # Iterate through all paths
        for path in parts:

            # Amount of vertices
            n = len(path)
            
            # process all segments
            for i in range(n):
                
                # Reduce coordinates
                xir = path[i].x() - q.x()
                xiir = path[(i+1)%n].x() - q.x()
                
                yir = path[i].y() - q.y()
                yiir = path[(i+1)%n].y() - q.y()
                
                # Point is vertex
                if abs(xir) < eps and abs(yir) < eps:
                    return 1

                # Ray intersects vertex
                if abs(yir) < eps:
                    iyr = path[(i-1)%n].y() - q.y()
                    if iyr*yiir < 0:
                        k = k + 1
                        
                # Point is on the Edge
                p_1 = array([xir, yir])
                p_2 = array([xiir, yiir])
                if abs(norm(-p_1) + norm(-p_2) - norm(p_2 - p_1)) < eps:
                    return 1
                
                # Suitable segment?
                if (((yiir > eps) and (yir <= eps)) 
                    or ((yir > eps) and (yiir <= eps))):
                    # Compute inersecion
                    xm = (xir*yiir-xiir*yir)/(yir - yiir)
                    if xm > eps:
                        k += 1
                    
        # Point inside polygon       
        return k%2
    
    def windingNumber(self, q : QPointF, polygon : Polygon) -> int:     
        """Winding number algorithm 
            
        Parameters
        ----------
        q : QPointF
            Analyzed point
        pol : Polygon
            polygon
            
        Returns
        -------
        1: Point is in the polygon
        0: Point is outside the polygon
        """
        # Split QPainterPath to QPolygonF
        parts = polygon.path.toSubpathPolygons()
        
        # Inicialize Winding Number
        
        eps = 1.0e-5                    # Epsilon for comparing floats
        wn = 0                         # In how many paths is the point
        # Iterate through all paths
        for path in parts:  
                      
            omega = 0
            n = len(path)
            # process all segments
            for i in range(n):
                
                startPoint = path[i]        # Starting point of edge
                endPoint = path[(i+1)%n]    # Ending point of edge
                
                # Point is vertex
                if startPoint == q or endPoint == q:
                    return 1
                
                # Compute determinant for analyze position point relative to the line
                det = self.pointEdgePosition(q, startPoint, endPoint)
                
                # Compute Winding Number
                angle = self.angleTwoLine(q, startPoint, endPoint)
                if det < -eps:    
                    omega += angle
                elif det > eps:                 
                    omega -= angle
                elif abs(angle - pi) < eps:           # Point is on the Edge
                    return 1
            if abs(abs(omega) - 2*pi) < eps:             # Point is inside the path
                wn += 1

        if wn == 1:                                     # if point is only inside one path
            return 1
        return 0
    
    def angleTwoLine(self, q : QPointF, startPoint : QPointF, endPoint : QPointF) -> float:
        """Compute angle between two line
            
        Parameters
        ----------
        q : QPointF
            Analyzed point
        startPoint : QPointF
            coordinates of starting point of the Edge
        endPoint : QPointF
            coordinates of ending point of the Edge
        Returns
        -------
        agle between lines
        """
        ux = startPoint.x() - q.x(); uy = startPoint.y() - q.y()
        vx = endPoint.x() - q.x(); vy = endPoint.y() - q.y()
        nu = (ux**2 + uy**2)**0.5; nv = (vx**2 + vy**2)**0.5
        nom = ux*vx + uy*vy; denom = nu*nv
        
        afi = min(max(nom/denom, -1), 1)
        return abs(acos(afi))
    
    def pointEdgePosition(self, q : QPointF, startPoint : QPointF, endPoint : QPointF) -> float:
        """is line on the left of the line
            
        Parameters
        ----------
        q : QPointF
            Analyzed point
        startPoint : QPointF
            coordinates of starting point of the Edge
        endPoint : QPointF
            coordinates of ending point of the Edge
        Returns
        -------
        Position of point
        """
        eps = finfo(float64).eps        # Epsilon for comparing float640
        ux = endPoint.x() - startPoint.x(); uy = endPoint.y() - startPoint.y()
        vx = q.x() - startPoint.x(); vy = q.y() - startPoint.y()
        return (ux*vy - uy*vx)