from PyQt6.QtCore import QPointF
from PyQt6.QtGui import QPolygonF
import json
from Polygon import Polygon
from math import inf

def load_polygons(path : str) -> list[QPolygonF]:
    """Load polygons
    
        Polygons from json file will be loaded. Coordinates will be scaled and translated 
        
        Parameters
        ----------
        path : str
            Path to the file
    """
    with open(path, 'r', encoding = 'utf-8') as f:
        file = json.load(f)
    polygons = file['POLYGONS']
    
    dx = 0; dy = 0; scale = 0
    min_x = inf; max_x = -inf; min_y = inf; max_y = -inf 
    
    for region in polygons:
        points = region['POINTS']
        for point in points:
            if point[0] == "N" or point[1] == "N":
                continue
            x = float(point[0])
            y = float(point[1])
            min_x = min(min_x, x)
            max_x = max(max_x, x)
            min_y = min(min_y, y)
            max_y = max(max_y, y) 
    
    width = 800
    height = 600
    offset = 50
    dx = (max_x + min_x)/2; dy = (max_y + min_y)/2
    
    if abs(max_x - min_x)/(width - offset) > abs(max_y - min_y)/(height - offset):
        scale = abs(max_x - min_x)/(width - offset)
    else:
        scale = abs(max_y - min_y)/(height - offset)
    
    painted_polygons = []    
    
    for region in polygons:
        points = region['POINTS']
        pol = Polygon()
        for point in points:
            if point[0] == "N" or point[1] == "N":
                continue
            new_x = (float(point[0]) - dx)/scale + width/2
            new_y = (-1)*((float(point[1]) - dy)/scale + height/2) + height - offset
            point = QPointF(new_x, new_y)
            pol.addVertex(point)
            pol.updateMMB(point)
        painted_polygons.append(pol)
    
    return painted_polygons, scale, dx, dy