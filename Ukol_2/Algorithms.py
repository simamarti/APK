from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QPolygonF
from PyQt6.QtWidgets import *
from math import acos, cos, sin, sqrt, inf, atan2, pi, floor
from numpy import array, cov
from numpy.linalg import svd
from Building import Building

class Algorithms:
    """ A class to algorithms
    
    Attributes
    ----------
    None
    
    Methods
    -------
    errorWindow(text):
        Error Window will be arised if error occurs

    analyzePointLinePosition(p, p1, p2):
        Analyze position of point and line
        
    getTwoLineAngle(p1, p2, p3, p4):
        Compute angle between two lines
        
    jarvisScan(pol):
        Jarvis Scan
        
    grahamScan(pol):
        Graham Scan
        
    mmb(pol):
        Compute min-max box for building

    rotate(pol, sig):
        Rotate building
        
    getArea(pol):
        Compute area of the building
        
    l2(p1, p2):
        Euclidian distance
        
    searchDiagonals(pol):
        Searching two longest diagonals
        
    weightedMean(sigma1, d1, sigma2, d2):
        Weighted mean
        
    slope(diagonal):
        Compute slope of the line
        
    resizeRectangle(rect, build):
        Resize rectangle to fit area of the building

    createMBR(pol, convexHull):
        Create minimum bounding rectangle
        
    createERPCA(pol):
        Create enclosing rectangle using PCA
        
    longestEdge(pol):
        Create enclosing rectangle using Longest Edge
        
    wallAverage(pol):
        Create enclosing rectangle using Wall Average
        
    weightedBisector(pol):
        Create enclosing rectangle using Weighted Bisector
        
    validation(buildings):
        Create enclosing rectangle using Weighted Bisector
    """
    def __init__(self):
        pass
    
    def errorWindow(self, text : str) -> None:
        """Error Window will be arised if error occurs
        
        Parameters
        ----------
        text : str
            Text for showing
        
        Returns
        -------
        None
        """
        messagebox = QtWidgets.QMessageBox()
        messagebox.setWindowTitle(text)
        messagebox.setText("Polygon contains less than 3 points. " + text + " cannot be created")
          
        # Show analysis
        messagebox.exec()
    
    def analyzePointLinePosition(self, p: QPointF, p1: QPointF, p2: QPointF) -> bool:
        """Analyze position of point and line
        
        Parameters
        ----------
        p : QPointF
            Analyze point
            
        p1 : QPointF
            Starting point of line
            
        p2 : QPointF
            Ending point of line
        
        Returns
        -------
        True: Point is in the left-halfplane
        False: Poin is not in the left-halfplane
        """
        # Compute vectors
        ux = p2.x() - p1.x(); uy = p2.y() - p1.y()
        vx = p.x() - p1.x(); vy = p.y() - p1.y()
        
        eps = 10e-10
        # Compute test
        t = ux*vy - uy*vx
        
        # Point in the left halfplane
        if t > eps:
            return True

        return False
      
    def getTwoLineAngle(self, p1 : QPointF, p2 : QPointF, p3 : QPointF, p4 : QPointF) -> float:
        """Compute angle between two lines
        
        Parameters
        ----------
        p1 : QPointF
            Starting point of first line
            
        p2 : QPointF
            Ending point of first line
        
        p3 : QPointF
            Starting point of second line
            
        p4 : QPointF
            Ending point of second line
            
        Returns
        -------
        angle between lines
        """
        # Get two line angle
        ux = p2.x() - p1.x(); uy = p2.y() - p1.y()
        vx = p4.x() - p3.x(); vy = p4.y() - p3.y()
        
        dot = ux*vx + uy*vy
        nu = (ux**2 + uy**2)**0.5
        nv = (vx**2 + vy**2)**0.5
        
        afi = min(max(dot/(nu*nv), -1), 1)
        
        return acos(afi)
    
    def jarvisScan(self, pol : QPolygonF) -> QPolygonF:
        """Jarvis Scan
        
        Parameters
        ----------
        pol : QPolygonF
            Building
            
        Returns
        -------
        ch : QPolygonF
            Convex Hull
        """
        # Convex Hull constructed using Jarvis scan
        ch = QPolygonF()

        # Find pivot 1
        q = min(pol,  key = lambda k: k.y())
        
        # Find pivot 2
        s = min(pol,  key = lambda k: k.x())
        
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
        
    def grahamScan(self, pol : QPolygonF) -> QPolygonF:
        """Graham Scan
        
        Parameters
        ----------
        pol : QPolygonF
            Building
            
        Returns
        -------
        ch : QPolygonF
            Convex Hull
        """
        # CH construction using Graham scan algorithm
        ch = QPolygonF()
        
        # Find pivot
        q = min(pol, key=lambda k : k.y())
        
        # Remove point from QPolygonF
        tmp = QPolygonF()
        for point in pol:
            if q != point:
                tmp.append(point)
        pol = tmp
        
        # Sort points by angle
        sorted_points = sorted(pol, key=lambda p: self.getTwoLineAngle(q, QPointF(q.x() - 10, q.y()), q, p), reverse = True)    
        
        # Initialize stack, add pivot and first point
        S = []
        S.append(q)
        S.append(sorted_points[0])
        n = len(pol)
        
        # Processing sorted points
        j= 2
        while j < len(sorted_points):
            pj = sorted_points[j]
            
            # Add point to the stack            
            if self.analyzePointLinePosition(pj, S[-2], S[-1]):
                S.append(pj)
                j += 1
            else:
                S.pop()
                
        for point in S:
            ch.append(point)
           
        return ch
    
    def mmb(self, pol : QPolygonF) -> QPolygonF:
        """Compute min-max box for building
        
        Parameters
        ----------
        pol : QPolygonF
            Building
            
        Returns
        -------
        box : QPolygonF
            min-max box
        """
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
        """Rotate building
        
        Parameters
        ----------
        pol : QPolygonF
            Building
        
        sig : float
            angle for rotating
            
        Returns
        -------
        polr : QPolygonF
            Rotated polygon
        """
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
    
    def getArea(self, pol : QPolygonF) -> float:
        """Compute area of the building
        
        Parameters
        ----------
        pol : QPolygonF
            Building
          
        Returns
        -------
        Area of the building
        """
        # polygon length
        n = len(pol)
        suma = 0
        
        # Processing of Vertices
        for idx in range(n):
            suma += pol[idx].x()*(pol[(idx + 1)%n].y() - pol[(idx - 1 + n)%n].y())
            
        return abs(suma)/2
    
    def l2(self, p1 : QPointF, p2 : QPointF) -> float:
        """Euclidian distance
        
        Parameters
        ----------
        p1 : QPointF
            First point
        
        p2 : QPointF
            Second point
            
        Returns
        -------
        Euclidian distance between points
        """
        
        return sqrt((p1.x() - p2.x())**2 + (p1.y() - p2.y())**2)
        
    def searchDiagonals(self, pol: QPolygonF) -> tuple:
        """Searching two longest diagonals
        
        Parameters
        ----------
        pol : QPolygonF
            Building
        
        Returns
        -------
        First - Longest diagonal
        d1 - length of longest diagonal
        Second - second longest diagonal
        d2 - length of second longest diagonal
        """
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
        """Weighted mean
        
        Parameters
        ----------
        sigma1 : float
            slope of first diagonal
        
        d1 : float
            lenght of first diagonal
            
        sigma2 : float
            slope of second diagonal
        
        d2 : float
            lenght of second diagonal
            
        Returns
        -------
        Weighted mean of angles
        """
        nom = sigma1*d1 + sigma2*d2
        denom = d1 + d2
        
        return nom/denom
    
    def slope(self, diagonal : list[QPointF]) -> float:
        """Compute slope of the line
        
        Parameters
        ----------
        diagonal : list[QPointF]
            diagonal
            
        Returns
        -------
        angle - slope of the diagonal
        """
        p1 = diagonal[0]; p2 = diagonal[1]
        d_x = p2.x() - p1.x()
        d_y = p2.y() - p1.y()
        # Find slope angle of longest edge
        angle = atan2(d_y, d_x)

        return angle
        
    def resizeRectangle(self, rect : QPolygonF, build : QPolygonF) -> QPolygonF:
        """Resize rectangle to fit area of the building
        
        Parameters
        ----------
        rect : QPolygonF
            generalized building (rectangle)
        
        build : QPolygonF
            original buildong
            
        Returns
        -------
        recR : QPolygonF
            resized generalized building
        """
        
        # Compute areas
        Ab = self.getArea(build)
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
        """Create minimum bounding rectangle
        
        Parameters
        ----------
        pol : QPolygonF
            building
        
        convexHull : function
            method for creating Convex Hull (Jarvis Scan, Graham Scan)
            
        Returns
        -------
        mmb_res : QPolygonF
            generalized building
        """
        
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
    
    def createERPCA(self, pol : QPolygonF) -> QPolygonF:
        """Create enclosing rectangle using PCA
        
        Parameters
        ----------
        pol : QPolygonF
            building
            
        Returns
        -------
        er_r : QPolygonF
            generalized building
        """

        if len(pol) < 3:
            return QPolygonF()
        # Create lis of coordinates
        x = []; y = []
        
        # Add coordinates to lists
        for point in pol:
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
        pol_unrot = self.rotate(pol, -sigma)
        
        # Find min-max box
        mmb = self.mmb(pol_unrot)
        
        # rotate min-max box (create emclosing rectangle)
        er = self.rotate(mmb, sigma)
        
        # Resize enclosing rectangle
        er_r = self.resizeRectangle(er, pol)
        
        return er_r
        
    def longestEdge(self, pol : QPolygonF) -> QPolygonF:
        """Create enclosing rectangle using Longest Edge
        
        Parameters
        ----------
        pol : QPolygonF
            building
            
        Returns
        -------
        er_r : QPolygonF
            generalized building
        """
        
        if len(pol) < 3:
            return QPolygonF()
        
        n = len(pol)
        # Initialize longest edge
        longest_edge = 0
        angle = 0
        
        # Process all edges
        for i in range(n):
            # Get length of the current edge
            edge_length = self.l2(pol[i], pol[(i + 1)%n])
            # Find the longest edge of the building
            if edge_length > longest_edge:
                longest_edge = edge_length
                angle = self.slope([pol[i], pol[(i + 1)%n]])
                
        # Rotate polygon by -sigma
        pol_unrot = self.rotate(pol, -angle)
        
        # Find min-max box
        mmb = self.mmb(pol_unrot)
        
        # Rotate min-max box
        er = self.rotate(mmb, angle)
        
        # Resize enclosing rectangle 
        er_r = self.resizeRectangle(er, pol)
        
        return er_r
    
    def wallAverage(self, pol : QPolygonF) -> QPolygonF:
        """Create enclosing rectangle using Wall Average
        
        Parameters
        ----------
        pol : QPolygonF
            building
            
        Returns
        -------
        er_r : QPolygonF
            generalized building
        """
        
        n = len(pol)
        if n < 3:
            return QPolygonF()
        
        # Calculate sigma
        d_x = pol[1].x() - pol[0].x()
        d_y = pol[1].y() - pol[0].y()
        sigma = atan2(d_y, d_x)
        
        # Average remainder
        av_r = 0
        
        # Process all edges
        for i in range(1, n):
            # Compute sigma for current segment
            d_x_i = pol[(i + 1)%n].x() - pol[i].x()
            d_y_i = pol[(i + 1)%n].y() - pol[i].y()
            sigma_i = atan2(d_y_i, d_x_i)
            # Compute direction difference
            dir_dif = sigma_i - sigma
            # Adjustment when direction difference is less than 0
            if dir_dif < 0:
                dir_dif += 2 * pi
            # Fraction by pi/2
            k_i = round((2 * dir_dif)/pi)
            # Calculate remainder for current segment
            rem_i = dir_dif - (k_i * (pi/2))
            # Update average remainder
            av_r += rem_i  
        
        # Compute the average remainder
        av_r = av_r/n
        
        # Average rotation
        sigma_av = sigma + av_r
        
        # Rotate polygon by -sigma
        pol_unrot = self.rotate(pol, -sigma_av)
        
        # Find min-max box
        mmb = self.mmb(pol_unrot)
        
        # Rotate min-max box
        er = self.rotate(mmb, sigma_av)
        
        # Resize enclosing rectangle 
        er_r = self.resizeRectangle(er, pol)
        
        return er_r

    def weightedBisector(self, pol : QPolygonF) -> QPolygonF:
        """Create enclosing rectangle using Weighted Bisector
        
        Parameters
        ----------
        pol : QPolygonF
            building
            
        Returns
        -------
        er_res : QPolygonF
            generalized building
        """
        
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
        
        er_res = self.resizeRectangle(er, pol)
        return er_res

    def validation(self, buildings : list[Building]) -> str:
        """Create enclosing rectangle using Weighted Bisector
        
        Parameters
        ----------
        buildings : list[Building]
            list of buildings
            
        Returns
        -------
        text : str
            text with results
        """
        text = f""
        sigma1sum = 0; sigma2sum = 0
        sigma1Perc = 0; sigma2Perc = 0
        
        # No data has been loaded or generalization has not been launched
        if not buildings or not buildings[0].getBuildingGeneralize():
            return text
        
        # Iterate through all buildings
        for building in buildings:          
            sigma1 = 0; sigma2 = 0
            
            # Compute r and k for generalized building
            rect = building.getBuildingGeneralize()
            slope = 0
            if self.l2(rect[0], rect[1]) > self.l2(rect[1], rect[2]):
                slope = self.slope([rect[0], rect[1]])
            else:
                slope = self.slope([rect[1], rect[2]])

            k = (2*slope)/pi
            r = (k - floor(k))*(pi/2)
            
            pol = building.getBuilding()
            
            # Iterate through each building
            n = len(pol)
            for idx in range(n):
                if pol[idx] == pol[(idx + 1)%n]:
                    continue
                
                slopeEdge = self.slope([pol[idx], pol[(idx + 1)%n]])
                
                ki = (2*slopeEdge)/pi
                ri = (ki - floor(ki))*(pi/2)
                
                sigma1 += (ri - r)
                sigma2 += pow(ri - r, 2)
                
            sigma2 = sqrt(sigma2)
            
            sigma1 *= (pi/(2*n)); sigma2 *= (pi/(2*n))
            sigma1 *= 180/pi; sigma2 *= 180/pi
            
            # Condition for acceptance
            if abs(sigma1) < 10: sigma1Perc += 1
            if abs(sigma2) < 10: sigma2Perc += 1
            
            sigma1sum += sigma1
            sigma2sum += sigma2
            
        sigma1sum /= len(buildings); sigma2sum /= len(buildings)
        sigma1Perc /= len(buildings); sigma2Perc /= len(buildings)

        # Generating text result
        text = f"Mean angular deviation:\t\t{sigma1sum:.2f}°\nSquare angular deviation:\t{sigma2sum:.2f}°\n"
        text += f"Percentage of buildings which met the condition (|\u03C3|<10°):\n"
        text += f"\tAngular deviation:\t\t{sigma1Perc*100:.2f} %\n\tSquare angular deviation:\t{sigma2Perc*100:.2f} %"
        return text