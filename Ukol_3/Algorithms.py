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
    """ A class for algorithms
    
    Attributes
    ----------
    None
    
    Methods
    -------
    getPointLinePosition(p, p1, p2):
        Analyze position of point and line
        
    getNearestPoint(q, points):
        Get nearest point to q
         
    getTwoLineAngle(p1, p2, p3, p4):
        Compute angle between two lines
        
    getDelaunayPoint(start, end, points):
        Get suitable Delaunay point    
    
    updateAEL(e, ael):
        Update AEL (Active Edge List)
    
    createDT(self, points):
        Create Delaunay triangulation
    
    getContourPoint(p1, p2, z):
        Get contour point
    
    createContourLines(dz, zmin, zmax, dt):
        Create contour lines for all region
    
    computeSlope(p1, p2, p3):
        Compute slope of Triangle
    
    computeAspect(p1, p2, p3):
        Compute Aspect of Triangle
    
    analyzeDTMSlope(dt):
        Analyze slope for all region
    
    analyzeDTMAspect(dt):
        Analyze Aspect for all region
        
    rayCrossingAlgorithm(q, polygon):
        Ray Crossing Algorithm for decide if point is inside of region
    
    intesectionTwoLines(self, p1, p2, pol):
        Check intersection of lines
        
    clipDt(self, dt, border, hole):
        Clip Delaunay triangulation by convex/non-convex polygon
    jarvisScan(pol):
        Jarvis Scan
    """

    def __init__(self):
        pass

    def getPointLinePosition(self, p : QPoint3DF, p1 : QPoint3DF, p2 : QPoint3DF):
        """Analyze position of point and line"""

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
        """Get nearest point to q"""
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
        """Compute angle between two lines"""

        # Get two line angle
        ux = p2.x() - p1.x(); uy = p2.y() - p1.y()
        vx = p4.x() - p3.x(); vy = p4.y() - p3.y()

        dot = ux*vx + uy*vy
        nu = (ux**2 + uy**2)**0.5
        nv = (vx**2 + vy**2)**0.5

        afi = min(max(dot/(nu*nv), -1), 1)

        return acos(afi)

    def getDelaunayPoint(self, start : QPoint3DF, end : QPoint3DF, points : list[QPoint3DF]):
        """Get suitable Delaunay point"""

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
        """Update AEL (Active Edge List)"""

        e_op = e.switchOriantation()
        if e_op in ael:
            ael.remove(e_op)
        else:
            ael.append(e)

    def createDT(self, points : list[QPoint3DF]) -> list[Edge]:
        """Create Delaunay triangulation"""

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
        """Get contour point"""

        # Intersection of triangle and intersection plane

        xb = (p2.x() - p1.x())/(p2.getZ() - p1.getZ()) * (z-p1.getZ()) + p1.x()
        yb = (p2.y() - p1.y())/(p2.getZ() - p1.getZ()) * (z-p1.getZ()) + p1.y()

        return QPoint3DF(xb, yb, z)

    def createContourLines(self, dz : int, zmin : float, zmax : float, dt : list[Edge]) -> list[Edge]:
        """Create contour lines for all region"""

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
        """Compute slope of Triangle"""

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
        """Compute Aspect of Triangle"""

        # Direction
        ux = p1.x() - p2.x(); uy = p1.y() - p2.y(); uz = p1.getZ() - p2.getZ()
        vx = p3.x() - p2.x(); vy = p3.y() - p2.y(); vz = p3.getZ() - p2.getZ()

        # Normal vector
        nx = uy * vz - vy * uz
        ny = -(ux * vz - vx * uz)

        return atan2(nx, ny)

    def analyzeDTMSlope(self, dt : list[Edge]):
        """Analyze slope for all region"""

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
        """Analyze Aspect for all region"""

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

    def rayCrossingAlgorithm(self, q : QPoint3DF, polygon : list[QPoint3DF]) -> bool:
        """Ray Crossing Algorithm for decide if point is inside of region
            
        Parameters
        ----------
        q : QPoint3DF
            Analyzed point
        pol : list[QPoint3DF]
            polygon
            
        Returns
        -------
        True: Point is in the polygon
        False: Point is outside the polygon
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
                return True

            # Ray intersects vertex
            if abs(yir) < eps:
                iyr = polygon[(i-1)%n].y() - q.y()
                if iyr*yiir < 0:
                    k = k + 1

            # Point is on the Edge
            p_1 = array([xir, yir])
            p_2 = array([xiir, yiir])
            if abs(norm(-p_1) + norm(-p_2) - norm(p_2 - p_1)) < eps:
                return True

            # Suitable segment?
            if (((yiir > eps) and (yir <= eps)) 
                or ((yir > eps) and (yiir <= eps))):
                # Compute inersecion
                xm = (xir*yiir-xiir*yir)/(yir - yiir)
                if xm > eps:
                    k += 1
                        
        # Point inside polygon 
        if k%2:
            return True 
        # Point outside polygon     
        return False
    
    def intesectionTwoLines(self, p1 : QPoint3DF, p2 : QPoint3DF, pol : list[QPoint3DF]) -> bool:
        """Check intersection of lines"""
        
        eps = 1.0e-5        # epsilon for comparing floats
        
        x1 = p1.x(); x2 = p2.x()
        y1 = p1.y(); y2 = p2.y()
        
        n = len(pol)
        for i in range(n):
            x3 = pol[i].x(); x4 = pol[(i + 1)%n].x()
            y3 = pol[i].y(); y4 = pol[(i + 1)%n].y()

            t1 = (x2 - x1)*(y4 - y1) - (x4 - x1)*(y2 - y1)
            t2 = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
            t3 = (x4 - x3)*(y1 - y3) - (x1 - x3)*(y4 - y3)
            t4 = (x4 - x3)*(y2 - y3) - (x2 - x3)*(y4 - y3)
        
            # Check if vertex of polygon is on the line
            u_x = x2 - x1; u_y = y2 - y1
            v_x = x3 - x1; v_y = y3 - y1
            
            t = u_x*v_y - v_x*u_y
            if abs(t) < eps:
                continue
            
            # One of the point of DT is on the border
            if (abs(x1 - x3) < eps and abs(y1 - y3) < eps) or (abs(x1 - x4) < eps and abs(y1 - y4) < eps) or\
                (abs(x2 - x3) < eps and abs(y2 - y3) < eps) or (abs(x2 - x4) < eps and abs(y2 - y4) < eps):
                    continue
                
            # Intersection exists
            if t1*t2 < eps and t3*t4 < eps:
                return True

        # Intersection doesn't exist
        return False
    
    def clipDt(self, dt : list[QPoint3DF], border : list[QPoint3DF], hole : list[QPoint3DF]):
        """Clip Delaunay triangulation by convex/non-convex polygon"""
        
        dt_clipped = []
        
        idx = 0
        for i in range(0, len(dt), 3):          # Iterate through all triangles
            # Get vertices of triangle
            p1 = dt[i].getStart(); p2 = dt[i].getEnd(); p3 = dt[i + 1].getEnd()
            
            # If all verticies of Delaunay Triangle is inside of the border (points are inside the border, edges of DT doesn't intersect border)
            if (self.rayCrossingAlgorithm(p1, border) and self.rayCrossingAlgorithm(p2, border) and self.rayCrossingAlgorithm(p3, border)) and\
                (not self.intesectionTwoLines(p1, p2, border) and not self.intesectionTwoLines(p2, p3, border) and not self.intesectionTwoLines(p3, p1, border)):
                
                # Add Triangle to final triangulation if hole doesn't exist or is outside of triangle (points are outside the hole, edges of DT doesn't intersect hole)
                if hole == [] or ((not self.rayCrossingAlgorithm(p1, hole) and not self.rayCrossingAlgorithm(p2, hole) and not self.rayCrossingAlgorithm(p3, hole)) and\
                    (not self.intesectionTwoLines(p1, p2, hole) and not self.intesectionTwoLines(p2, p3, hole) and not self.intesectionTwoLines(p3, p1, hole))):

                    dt_clipped.append(dt[i])
                    dt_clipped.append(dt[i + 1])
                    dt_clipped.append(dt[i + 2])
        
        return dt_clipped

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