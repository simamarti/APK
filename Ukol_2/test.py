from PyQt6.QtWidgets import QWidget
from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QPolygonF
from copy import deepcopy
from math import acos, atan2, pi

def getTwoLineAngle(p1 : QPointF, p2 : QPointF, p3 : QPointF, p4 : QPointF) -> float:
        # Get two line angle
        ux = p2.x() - p1.x(); uy = p2.y() - p1.y()
        vx = p4.x() - p3.x(); vy = p4.y() - p3.y()
        
        dot = ux*vx + uy*vy
        nu = (ux**2 + uy**2)**0.5
        nv = (vx**2 + vy**2)**0.5
        
        afi = min(max(dot/(nu*nv), -1), 1)
        
        return acos(afi)
    
def slope(diagonal : list[QPointF, QPointF]) -> None:
        
        p1 = diagonal[0]; p2 = diagonal[1]
        p3 = QPointF(0, 0); p4 = QPointF(100, 0)              #x axis
        if p2.y() < p1.y():
            tmp = deepcopy(p1)
            p1 = deepcopy(p2)
            p2 = deepcopy(tmp)
        
        sigma = getTwoLineAngle(p1, p2, p3, p4)
        print(f"\tMetoda 1: {sigma*180/pi}°")
        d_x = p2.x() - p1.x()
        d_y = p2.y() - p1.y()
        # Find slope angle of longest edge
        angle = atan2(d_y, d_x)
        print(f"\tMetoda 2: {angle*180/pi}°")

p1 = QPointF(0, 0)
p2 = QPointF(100, 100)
print("[0, 0], [100, 100]")
slope([p1, p2])
print("[100, 100], [0, 0]")
slope([p2, p1])

p1 = QPointF(0, 0)
p2 = QPointF(-100, 100)
print("[0, 0], [-100, 100]")
slope([p1, p2])
print("[-100, 100], [0, 0]")
slope([p2, p1])