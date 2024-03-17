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
    
    def rayCrossingAlgorithm(self, q : QPointF, pol : QPolygonF) -> bool:
        """Ray Crossing algorithm 
            
        Parameters
        ----------
        q : QPointF
            Analyzed point
        pol : QPolygonF
            polygon
            
        Returns
        -------
        True: Point is in the polygon
        False: Point is outside the polygon
        """
        # Inicialize amount of intersection
        k = 0
        
        # Amount of vertices
        n = len(pol)
        
        # process all segments
        for i in range(n):
            
            # Reduce coordinates
            xir = pol[i].x() - q.x()
            xiir = pol[(i+1)%n].x() - q.x()
            
            yir = pol[i].y() - q.y()
            yiir = pol[(i+1)%n].y() - q.y()
            
            # Point is vertex
            if not xir and not yir:
                print("Point is vertex")
                return 1

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
        return bool(k%2)