from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QPolygonF
from Edge import Edge
from QPoint3DF import QPoint3DF
from Triangle import Triangle
from math import inf, sqrt, acos, atan2
from numpy import arange
from numpy.linalg import norm
from numpy import array


class Algorithms:
    
    def __init__(self):
        pass
    
    def getPointLinePosition(self, p : QPoint3DF, p1 : QPoint3DF, p2 : QPoint3DF):
        
        ux = p2.x() - p1.x(); uy = p2.y() - p1.y()
        vx = p.x() - p1.x(); vy = p.y() - p1.y()
        
        eps = 10e-10
        # Compute test
        t = ux*vy - uy*vx
        # Point in the left halfplane
        if t > eps:
            return 1
        
        # Point in the right halfplane
        if t < -eps:
            return 0
        
        # Point on the line
        return -1
        
    def getNearestPoint(self, q : QPoint3DF, points : list[QPoint3DF]) -> QPoint3DF:
        # Find point nearest to q
        
        pNearest = None
        dNearest = inf
        
        # Process all points
        for point in points:
            if point != q:
                
                #Compute distance
                dx = q.x() - point.x()
                dy = q.y() - point.y()
                dist = sqrt(dx**2 + dy**2)
                if dist < dNearest:
                    dNearest = dist
                    pNearest = point
        return pNearest
    
    def getTwoLineAngle(self, p1 : QPointF, p2 : QPointF, p3 : QPointF, p4 : QPointF) -> float:
        # Get two line angle
        ux = p2.x() - p1.x(); uy = p2.y() - p1.y()
        vx = p4.x() - p3.x(); vy = p4.y() - p3.y()
        
        dot = ux*vx + uy*vy
        nu = (ux**2 + uy**2)**0.5
        nv = (vx**2 + vy**2)**0.5
        
        afi = min(max(dot/(nu*nv), -1), 1)
        
        return acos(afi)
    
    def getDelaunayPoint(self, start : QPoint3DF, end : QPoint3DF, points : list[QPoint3DF]):
        # Delaunay point to edge
        pDt = None
        angleMax = 0

        # Process all points
        for point in points:
            
            # Point p is different from q
            if start != point and end != point:
                
                # Point in the left halfplane
                if self.getPointLinePosition(point, start, end) == 1:

                    angle = self.getTwoLineAngle(point, start, point, end)
                    
                    # Update maximum
                    if angle > angleMax:
                        angleMax = angle
                        pDt = point
        return pDt
    
    def updateAEL(self, e : Edge, ael : list[Edge]) -> list[Edge]:
        e_op = e.switchOriantation()
        if e_op in ael:
            ael.remove(e_op)
        else:
            ael.append(e)
    
    def createDT(self, points : list[QPoint3DF]) -> list[Edge]:
        # Create Delaunay triangulation with incremental method
        
        ael : list[Edge] = []
        dt = []
        
        # Sort points by x
        
        p1 = min(points, key = lambda p: p.x())
        
        # Find nearest point
        p2 = self.getNearestPoint(p1, points)
        
        # Create edges
        e = Edge(p1, p2)
        e_op = Edge(p2, p1)
        
        # Add both edges to ael
        ael.append(e)
        ael.append(e_op)
        
        # Repeat until ael is empty
        while ael:
            # Take first edge
            e1 = ael.pop()
            
            # Change orientation
            e1_op = e1.switchOriantation()
            
            # Find optimal Delaunay point
            p_dt = self.getDelaunayPoint(e1_op.getStart(), e1_op.getEnd(), points)
            
            # Did we find a suitable point?
            if p_dt != None:
                # create remaing edges
                e2 = Edge(e1_op.getEnd(), p_dt)
                e3 = Edge(p_dt, e1_op.getStart())
                
                # Create Delaunay triangle
                dt.append(e1_op)
                dt.append(e2)
                dt.append(e3)
                
                # Update AEL
                self.updateAEL(e2, ael)
                self.updateAEL(e3, ael)
        return dt
    
    def getContourPoint(self, p1 : QPoint3DF, p2 : QPoint3DF, z : float) -> QPoint3DF:
        # Intersection of triangle and intersection plane
        
        xb = (p2.x() - p1.x())/(p2.getZ() - p1.getZ()) * (z-p1.getZ()) + p1.x()
        yb = (p2.y() - p1.y())/(p2.getZ() - p1.getZ()) * (z-p1.getZ()) + p1.y()
        
        return QPoint3DF(xb, yb, z)
    
    def createContourLines(self, dz : int, zmin : float, zmax : float, dt : list[Edge]) -> list[Edge]:
        # Create contour lines defined by interval and step
        
        contours : list[Edge] = []
        
        for i in arange(0, len(dt), 3):          # Iterate through all triangles
            # Get vertices of triangle
            p1 = dt[i].getStart(); p2 = dt[i].getEnd(); p3 = dt[i + 1].getEnd()
            
            # Get z coordinates
            z1 = p1.getZ(); z2 = p2.getZ(); z3 = p3.getZ()
            # Create all conture lines
            for z in arange(zmin, zmax, dz):
                # Compute edges' height differences
                dz1 = z - z1; dz2 = z - z2; dz3 = z - z3
                
                # Skip coplanar triangle
                if dz1 == 0 and dz2 == 0 and dz3 ==0:
                    continue
                
                #Edge p1 and p2 is colinear
                elif dz1 == 0 and dz2 == 0:
                    contours.append(dt[i])
                    
                #Edge p2 and p3 is colinear
                elif dz2 == 0 and dz3 == 0:
                    contours.append(dt[i + 1])
                    
                #Edge p3 and p1 is colinear
                elif dz3 == 0 and dz1 == 0:
                    contours.append(dt[i + 2])
                
                # Edges p1, p2 and p2, p3 intersected by plane
                elif dz1*dz2 <= 0 and dz2*dz3 <=0:
                    # Compute intersection
                    a = self.getContourPoint(p1, p2, z)
                    b = self.getContourPoint(p2, p3, z)
                    
                    # Create contour point
                    e = Edge(a, b)
                    
                    # Add edge to contour line
                    contours.append(e)
                    
                # Edges p2, p3 and p3, p1 intersected by plane
                elif dz2*dz3 <= 0 and dz3*dz1 <=0:
                    # Compute intersection
                    a = self.getContourPoint(p2, p3, z)
                    b = self.getContourPoint(p3, p1, z)
                    
                    # Create contour point
                    e = Edge(a, b)
                    
                    # Add edge to contour line
                    contours.append(e)
                    
                # Edges p3, p1 and p1, p2 intersected by plane
                elif dz3*dz1 <= 0 and dz1*dz2 <=0:
                    # Compute intersection
                    a = self.getContourPoint(p3, p1, z)
                    b = self.getContourPoint(p1, p2, z)
                    
                    # Create contour point
                    e = Edge(a, b)
                    
                    # Add edge to contour line
                    contours.append(e)

        return contours

    def computeSlope(self, p1 : QPoint3DF, p2 : QPoint3DF, p3 : QPoint3DF):
        # Compute triangle sope

        # Direction
        ux = p1.x() - p2.x(); uy = p1.y() - p2.y(); uz = p1.getZ() - p2.getZ()
        vx = p3.x() - p2.x(); vy = p3.y() - p2.y(); vz = p3.getZ() - p2.getZ()

        # Normal vector
        nx = uy * vz - vy * uz
        ny = -(ux * vz - vx * uz)
        nz = ux * vy - vx * uy

        # Norm
        norm = (nx**2 + ny**2 + nz**2)**(1/2)

        return acos(abs(nz)/norm)
    
    def computeAspect(self, p1 : QPoint3DF, p2 : QPoint3DF, p3 : QPoint3DF):
        # Direction
        ux = p1.x() - p2.x(); uy = p1.y() - p2.y(); uz = p1.getZ() - p2.getZ()
        vx = p3.x() - p2.x(); vy = p3.y() - p2.y(); vz = p3.getZ() - p2.getZ()

        # Normal vector
        nx = uy * vz - vy * uz
        ny = -(ux * vz - vx * uz)
        
        return atan2(nx, ny)
    
    def analyzeDTMSlope(self, dt : list[Edge]):
    
        dtm_slope : list[Triangle] = []

        for i in range(0, len(dt), 3):          # Iterate through all triangles
            # Get vertices of triangle
            p1 = dt[i].getStart(); p2 = dt[i].getEnd(); p3 = dt[i + 1].getEnd()

            # Get slope
            slope = self.computeSlope(p1, p2, p3)

            # Create triangle
            triangle = Triangle(p1, p2, p3, slope, 0)

            dtm_slope.append(triangle)

        return dtm_slope
    
    def analyzeDTMAspect(self, dt : list[Edge]):
    
        dtm_aspect : list[Triangle] = []

        for i in range(0, len(dt), 3):          # Iterate through all triangles
            # Get vertices of triangle
            p1 = dt[i].getStart(); p2 = dt[i].getEnd(); p3 = dt[i + 1].getEnd()

            # Get slope
            aspect = self.computeAspect(p1, p2, p3)

            # Create triangle
            triangle = Triangle(p1, p2, p3, 0, aspect)

            dtm_aspect.append(triangle)

        return dtm_aspect
    
    def rayCrossingAlgorithm(self, q : QPoint3DF, polygon : list[QPoint3DF]) -> int:
        """Ray Crossing algorithm 
            
        Parameters
        ----------
        q : QPoint3DF
            Analyzed point
        pol : list[QPoint3DF]
            polygon
            
        Returns
        -------
        0: Point is in the polygon
        1: Point is outside the polygon
        """
        
        # Inicialize amount of intersection
        k = 0
        eps = 1.0e-5        # epsilon for comparing floats
        
        # Amount of vertices
        n = len(polygon)
            
        # process all segments
        for i in range(n):
                
            # Reduce coordinates
            xir = polygon[i].x() - q.x()
            xiir = polygon[(i+1)%n].x() - q.x()
                
            yir = polygon[i].y() - q.y()
            yiir = polygon[(i+1)%n].y() - q.y()
                
            # Point is vertex
            if abs(xir) < eps and abs(yir) < eps:
                return 1

            # Ray intersects vertex
            if abs(yir) < eps:
                iyr = polygon[(i-1)%n].y() - q.y()
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
    
    def clipDt(self, dt : list[QPoint3DF], border : list[QPoint3DF]):
        """Clip Delaunay triangulation by convex/non-convex polygon"""
        print("Border: ")
        for pt in border:
            print(f"[{pt.x()}, {pt.y()}]")
            
        dt_clipped = []
        
        for i in range(0, len(dt), 3):          # Iterate through all triangles
            # Get vertices of triangle
            p1 = dt[i].getStart(); p2 = dt[i].getEnd(); p3 = dt[i + 1].getEnd()

            # Centroid of triangle
            cX = (p1.x() + p2.x() + p3.x())/3; cY = (p1.y() + p2.y() + p3.y())/3
            p = QPoint3DF(cX, cY, 0)
            
            if not self.rayCrossingAlgorithm(p, border):
                print(f"> [{p.x()}, {p.y()}]")
                dt_clipped.append(dt[i])
                dt_clipped.append(dt[i + 1])
                dt_clipped.append(dt[i + 2])
                
        return dt_clipped
    
    def getTwoLineAngle(self, p1 : QPoint3DF, p2 : QPoint3DF, p3 : QPoint3DF, p4 : QPoint3DF) -> float:
        # Get two line angle
        ux = p2.x() - p1.x(); uy = p2.y() - p1.y()
        vx = p4.x() - p3.x(); vy = p4.y() - p3.y()
        
        dot = ux*vx + uy*vy
        nu = (ux**2 + uy**2)**0.5
        nv = (vx**2 + vy**2)**0.5
        
        afi = min(max(dot/(nu*nv), -1), 1)
        
        return acos(afi)
    
    def jarvisScan(self, pol : list[QPoint3DF]) -> list[QPoint3DF]:
        # Convex Hull constructed using Jarvis scan
        ch = []

        # Find pivot 1
        q = min(pol,  key = lambda k: k.y())
        
        # Find pivot 2
        s = min(pol,  key = lambda k: k.x())
        
        # Initialize last 2 points CH        
        qj = q
        qj1 = QPoint3DF(s.x() - 10, q.y(), 0)

        # Add pivot to ch
        ch.append(q)
        
        # Process all points
        while True:
            # Maximum and its index
            omegaMax = 0
            indexMax = -1
            
            for i in range(len(pol)):
                
                if qj != pol[i]:
                    omega = self.getTwoLineAngle(qj, qj1, qj, pol[i])
                    
                    # Update maximum
                    if omega > omegaMax:
                        omegaMax = omega
                        indexMax = i
                    
            # Add new point po convex hull
            ch.append(pol[indexMax])
            
            
            # We found pivot again
            if pol[indexMax] == q: 
                break
            
            # Update last segment
            qj1 = qj
            qj = pol[indexMax]
    
        return ch