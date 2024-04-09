from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QPolygonF
from PyQt6.QtWidgets import *
from math import acos, cos, sin, sqrt, inf, atan2
from numpy import array, cov
from numpy.linalg import svd
from Building import Building
from copy import deepcopy

class Algorithms:
    
    def __init__(self):
        pass
    
    def errorWindow(self, text : str) -> None:
        messagebox = QtWidgets.QMessageBox()
        messagebox.setWindowTitle(text)
        messagebox.setText("Polygon contains less than 3 points. " + text + " cannot be created")
          
        # Show analysis
        messagebox.exec()
        
    def getTwoLineAngle(self, p1 : QPointF, p2 : QPointF, p3 : QPointF, p4 : QPointF) -> float:
        # Get two line angle
        ux = p2.x() - p1.x(); uy = p2.y() - p1.y()
        vx = p4.x() - p3.x(); vy = p4.y() - p3.y()
        
        dot = ux*vx + uy*vy
        nu = (ux**2 + uy**2)**0.5
        nv = (vx**2 + vy**2)**0.5
        
        afi = min(max(dot/(nu*nv), -1), 1)
        
        return acos(afi)
    
    def jarvisScan(self, pol : Building) -> QPolygonF:
        # Convex Hull constructed using Jarvis scan
        ch = QPolygonF()

        # Find pivot 1
        polygon = pol.building
        q = min(polygon,  key = lambda k: k.y())
        
        # Find pivot 2
        s = min(polygon,  key = lambda k: k.x())
        
        # Initialize last 2 points CH        
        qj = q
        qj1 = QPointF(s.x() - 10, q.y())

        # Add pivot to ch
        ch.append(q)
        
        # Process all points
        while True:
            # Maximum and its index
            omegaMax = 0
            indexMax = -1
            
            for i in range(len(polygon)):
                
                if qj != polygon[i]:
                    omega = self.getTwoLineAngle(qj, qj1, qj, polygon[i])
                    
                    # Update maximum
                    if omega > omegaMax:
                        omegaMax = omega
                        indexMax = i
                    
            # Add new point po convex hull
            ch.append(polygon[indexMax])
            
            
            # We found pivot again
            if polygon[indexMax] == q: 
                break
            
            # Update last segment
            qj1 = qj
            qj = polygon[indexMax]
    
        return ch
        
    def lidaScan(self, pol : QPolygonF) -> QPolygonF:
        raise NotImplementedError("Method Lída Scan has not been implemented yet")
    
    def mmb(self, pol : QPolygonF) -> QPolygonF:
        # Compute min_max box
        
        # Compute points with extreme coordinates
        px_min = min(pol, key = lambda k : k.x()); px_max = max(pol, key = lambda k : k.x())
        py_min = min(pol, key = lambda k : k.y()); py_max = max(pol, key = lambda k : k.y())
        
        # Compute min-max box points
        v1 = QPointF(px_min.x(), py_min.y())
        v2 = QPointF(px_max.x(), py_min.y())
        v3 = QPointF(px_max.x(), py_max.y())
        v4 = QPointF(px_min.x(), py_max.y())
        
        # Create min max box
        box = QPolygonF([v1, v2, v3, v4])
        
        return box
        
    def rotate(self, pol : QPolygonF, sig : float) -> QPolygonF:
        # Rotate polygon by given angle
        
        polr = QPolygonF()
        
        # Process all points
        for p in pol:
            # Rotae point
            xr = cos(sig)*p.x() - sin(sig)*p.y()
            yr = sin(sig)*p.x() + cos(sig)*p.y()
            
            # Create roated point
            p_r = QPointF(xr, yr)
            
            # Add point to polygon
            polr.append(p_r)
            
        return polr
    
    def getArea(sell, pol : QPolygonF) -> float:
        # Return polygon area
        # polygon length
        n = len(pol)
        suma = 0
        
        # Processing of Vertices
        for idx in range(n):
            suma += pol[idx].x()*(pol[(idx + 1)%n].y() - pol[(idx - 1 + n)%n].y())
            
        return abs(suma)/2
    
    def l2(self, p1 : QPointF, p2 : QPointF) -> float:
        # Euclidian distance
        
        return sqrt((p1.x() - p2.x())**2 + (p1.y() - p2.y())**2)
        
    def searchDiagonals(self, pol: QPolygonF):
        
        first = [pol[0], pol[1]]; d1 = self.l2(pol[0], pol[1])
        second = [pol[0], pol[2]]; d2 = self.l2(pol[0], pol[2])
        n = len(pol)
        
        for i in range(n):
            for j in range(i+1, n):
                d = self.l2(pol[i], pol[j])
                if d > d1:
                    d1 = d
                    second = [first[0], first[1]]
                    first = [pol[i], pol[j]]
                elif d > d2:
                    second = [pol[i], pol[j]]

        return first, d1, second, d2
    
    def weightedMean(self, sigma1 : float, d1 : float, sigma2 : float, d2 : float) -> float:
        nom = sigma1*d1 + sigma2*d2
        denom = d1 + d2
        
        return nom/denom
    
    def slope(self, diagonal : list[QPointF, QPointF]) -> float:
        
        p1 = diagonal[0]; p2 = diagonal[1]
        p3 = QPointF(0, 0); p4 = QPointF(100, 0)              #x axis
        if p2.y() < p1.y():
            tmp = deepcopy(p1)
            p1 = deepcopy(p2)
            p2 = deepcopy(tmp)
        
        sigma = self.getTwoLineAngle(p1, p2, p3, p4)
            
        return sigma
        
    def resizeRectangle(self, rect : QPolygonF, build : Building) -> QPolygonF:
        # Resize rectangle to fit area of the building
        
        # Compute areas
        Ab = self.getArea(build.building)
        A = self.getArea(rect)
        
        # Compute ratio
        k = Ab/A
        
        # Center of mass
        tx = (rect[0].x() + rect[1].x() + rect[2].x() + rect[3].x())/4
        ty = (rect[0].y() + rect[1].y() + rect[2].y() + rect[3].y())/4
        
        # Vectors
        u1x = rect[0].x() - tx; u1y = rect[0].y() - ty
        u2x = rect[1].x() - tx; u2y = rect[1].y() - ty
        u3x = rect[2].x() - tx; u3y = rect[2].y() - ty
        u4x = rect[3].x() - tx; u4y = rect[3].y() - ty
        
        # New vertices
        p1 = QPointF(tx + sqrt(k)*u1x, ty + sqrt(k)*u1y)
        p2 = QPointF(tx + sqrt(k)*u2x, ty + sqrt(k)*u2y) 
        p3 = QPointF(tx + sqrt(k)*u3x, ty + sqrt(k)*u3y) 
        p4 = QPointF(tx + sqrt(k)*u4x, ty + sqrt(k)*u4y)
        
        # Append QPointF to QPolygonF
        rectR = QPolygonF([p1, p2, p3, p4])
        
        return rectR
    
    def createMBR(self, pol : QPolygonF, convexHull) -> QPolygonF:
        # Create minimum bounding rectangle
        
        if len(pol) < 3:
            return QPolygonF()
        
        # Compute convex hull
        ch = convexHull(pol)
        
        # Initialization of mmb_min, area_min
        mmb_min = self.mmb(ch)
        area_min = self.getArea(mmb_min)
        sigma_min = 0
        
        n = len(ch)
        # Process all segments in ch
        for idx in range(n):
            dx = ch[(idx+1)%n].x() - ch[idx].x()
            dy = ch[(idx+1)%n].y() - ch[idx].y()
            
            # Direction
            sigma = atan2(dy, dx)
            
            # Rotate convex hull by -sigma
            ch_rot = self.rotate(ch, -sigma)    # vstup QPolygonF, výstup QPolygonF
            
            # Find mmb and its area
            mmb_rot = self.mmb(ch_rot)      # vstup QPolygonF, výstup - QPolygonF
            area_rot = self.getArea(mmb_rot)
            
            # Is it a better approximation?
            if area_rot < area_min:
                mmb_min = mmb_rot
                area_min = area_rot
                sigma_min = sigma
        
        # Back rotation
        mmb_unrot = self.rotate(mmb_min, sigma_min) # vstup - QPolygonF
        
        # Resize rectangle
        mmb_res = self.resizeRectangle(mmb_unrot, pol)
        
        return mmb_res
    
    def createERPCA(self, pol : Building) -> QPolygonF:
        # Create enclosing rectangle using PCA

        if len(pol) < 3:
            return QPolygonF()
        # Create lis of coordinates
        x = []; y = []
        
        # Add coordinates to lists
        for point in pol.building:
            x.append(point.x())
            y.append(point.y())
        
        # Create numpy array    
        P = array([x, y])
        
        # Compute covariance matrix
        C = cov(P)
        
        # Compute Singular Value Decomposition
        [U, S, V] = svd(C)
        
        # Compute sigma
        sigma = atan2(U[0][1], U[0][0])
        
        # Rotate polygon by minus sigma
        pol_unrot = self.rotate(pol.building, -sigma)
        
        # Find min-max box
        mmb = self.mmb(pol_unrot)
        
        # rotate min-max box (create emclosing rectangle)
        er = self.rotate(mmb, sigma)
        
        # Resize enclosing rectangle
        er_r = self.resizeRectangle(er, pol)
        
        return er_r
        
    def longestEdge(self, pol : QPolygonF) -> QPolygonF:
        raise NotImplementedError("Method longest Edge has not been implemented yet")
    
    def wallAverage(self, pol : QPolygonF) -> QPolygonF:
        raise NotImplementedError("Method Wall Average has not been implemented yet")

    def weightedBisector(self, pol : QPolygonF) -> QPolygonF:
    
        if len(pol) < 3:
            return QPolygonF()
        
        # Search two longest diagonals
        first, d1, second, d2 = self.searchDiagonals(pol)
        
        # Compute slope of diagonals
        sigma1 = self.slope(first)
        sigma2 = self.slope(second)

        # Compute weighted slope     
        sigma = self.weightedMean(sigma1, d1, sigma2, d2)

        # Rotate polygon by minus sigma
        pol_unrot = self.rotate(pol, -sigma)
        
        # Find min-max box
        mmb = self.mmb(pol_unrot)
        
        # rotate min-max box (create emclosing rectangle)
        er = self.rotate(mmb, sigma)
        
        return er