import json
from PyQt6.QtGui import QPolygonF
from PyQt6.QtCore import QPointF
from math import inf
from Building import Building
from PyQt6 import QtWidgets

def load_buildings(path : str) -> list[QPolygonF]:
    errorText = ""
    try:
        with open(path, 'r', encoding = 'utf-8') as f:
            file = json.load(f)
    
            
        polygons = file['POLYGONS']
        
        dx = 0; dy = 0; scale = 0
        min_x = inf; max_x = -inf; min_y = inf; max_y = -inf 
        
        # Compute min_max box
        for region in polygons:
            points = region['POINTS']
            
            for p in points:
                if p[0] == 'N' or p[1] == 'N' or p[0] == 'P' or p[1] == 'P':
                    continue
                min_x = min(min_x, float(p[0]))
                max_x = max(max_x, float(p[0]))
                min_y = min(min_y, float(p[1]))
                max_y = max(max_y, float(p[1])) 
    
    except FileNotFoundError:
        errorText = "File has not bee"
        return []
    except KeyError:
        errorText = "JSON file has not proper structure."
        return []
    except PermissionError:
        errorText = "You have not permissions to open file."
        return []
    except IndexError:
        errorText = "Wrong format of data."
        return []
    except IOError:
        errorText = "Error during input."
        return []
        
    finally:   
        if errorText != "":
            messagebox = QtWidgets.QMessageBox()
            messagebox.setWindowTitle(errorText)
            messagebox.setText(errorText)
            
            # Show analysis
            messagebox.exec()
            
    width = 800
    height = 600
    offset = 50
    # Compute shifts and scalling factor
    dx = (max_x + min_x)/2; dy = (max_y + min_y)/2
    print(abs(max_x - min_x)/(width - offset))
    print(abs(max_y - min_y)/(height - offset))
    scale = max(1, abs(max_x - min_x)/(width - offset), abs(max_y - min_y)/(height - offset))

    # transform coordinates
    painted_polygons = transformPolygons(polygons, scale, dx, dy, width, height, offset)
   
    return painted_polygons

def transformPolygons(polygons : list[dict], scale : float, dx : float, dy : float, width : int, height : int, offset : int) -> list[Building]:
    """Transform polygons to new coordinates which are compatible with canvas coordinates
        
        Parameters
        ----------
        polygons : list[dict]
            list of dictionaries with points
        scale : float
            scale factor for transformation
        dx : float
            shift along the x axis
        dy : float
            shift along the y axis
        width : int
            width of canvas
        height : int
            height of canvas
        offset : int
            offset for drawing on the canvas
            
        Returns
        -------
        painted polygons : list[Polygon]
            list of polygons with transformed coordinates
    """
    painted_polygons = []    
    
    for region in polygons:
        pol = Building()
        points = region['POINTS']
        if len(points) < 3:
            continue
        
        for point in points:
            if point[0] == 'N' or point[1] == 'N' or point[0] == 'P' or point[1] == 'P':
                continue
            
            new_x = (float(point[0]) - dx)/scale + width/2
            new_y = (-1)*((float(point[1]) - dy)/scale + height/2) + height - offset
            point = QPointF(new_x, new_y)
            pol.addVertex(point)
        painted_polygons.append(pol)
        
    return painted_polygons