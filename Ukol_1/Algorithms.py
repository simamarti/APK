from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from numpy.linalg import norm
from numpy import array, finfo, float64
from Polygon import Polygon

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
        parts = polygon.path.toSubpathPolygons()
        # Inicialize amount of intersection
        k = 0
        
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
                if not xir and not yir:
                    return 1

                # Ray intersects vertex
                if not yir:
                    iyr = path[(i-1)%n].y() - q.y()
                    if iyr*yiir < 0:
                        k = k + 1
                        
                # Point is on the Edge
                p_1 = array([xir, yir])
                p_2 = array([xiir, yiir])
                if (norm(-p_1) + norm(-p_2) - norm(p_2 - p_1)) < finfo(float64).eps:
                    return 1
                
                # Suitable segment?
                if ((yiir > 0) and (yir <= 0)) or ((yir > 0) and (yiir <= 0)):
                    # Compute inersecion
                    xm = (xir*yiir-xiir*yir)/(yir - yiir)
                    if xm > 0:
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
        parts = polygon.path.toSubpathPolygons()
        # Inicialize amount of intersection
        k = 0
        
        for path in parts:
            
            n = len(path)
            # process all segments
            for i in range(n):
                
                # Reduce coordinates
                xir = path[i].x() - q.x()
                xiir = path[(i+1)%n].x() - q.x()
                
                yir = path[i].y() - q.y()
                yiir = path[(i+1)%n].y() - q.y()
                
                # Point is vertex
                if not xir and not yir:
                    return 1

                # Ray intersects vertex
                if not yir:
                    iyr = path[(i-1)%n].y() - q.y()
                    if iyr*yiir < 0:
                        k = k + 1
                        
                # Point is on the Edge
                p_1 = array([xir, yir])
                p_2 = array([xiir, yiir])
                if (norm(-p_1) + norm(-p_2) - norm(p_2 - p_1)) < finfo(float64).eps:
                    return 1
                
                if ((yiir > 0) and (yir <= 0)):
                    # Compute inersecion
                    xm = (xir*yiir-xiir*yir)/(yir - yiir)
                    if xm > 0:
                        k += 1
                elif ((yir > 0) and (yiir <= 0)):
                    # Compute inersecion
                    xm = (xir*yiir-xiir*yir)/(yir - yiir)
                    if xm > 0:
                        k -= 1
        return k