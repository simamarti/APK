import json
from PyQt6.QtGui import QPolygonF
from PyQt6.QtCore import QPointF
from math import inf
from Building import Building
from PyQt6 import QtWidgets

def load_buildings(path : str) -> list[Building]:
    """Load from file and transform polygons to new coordinates which are compatible with canvas coordinates
        
        Parameters
        ----------
        path : str
            relative/absolute path to the file with coordinates of buildings
            
        Returns
        -------
        painted polygons : list[Building]
            list of buildings with transformed coordinates
    """
    
    errorText = ""
    try:                
        # Avoid falling of the program due to the problems with working with IO stream
        with open(path, 'r', encoding = 'utf-8') as f:
            file = json.load(f)
    
            
        polygons = file['POLYGONS']
        
        # Initilate scaling factor, shift and min/max of coordinates
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
        errorText = "File has not been found."
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
        # After try block in the case of error arise message box with error
        if errorText != "":
            messagebox = QtWidgets.QMessageBox()
            messagebox.setWindowTitle(errorText)
            messagebox.setText(errorText)
            
            # Show analysis
            messagebox.exec()
            
    width = 800
    height = 600
    offset = 50
    scale = 0
    # Compute shifts and scalling factor
    dx = (max_x + min_x)/2; dy = (max_y + min_y)/2
    scale = max(1, abs(max_y - min_y)/(height - offset), abs(max_x - min_x)/(width - offset))  

    # transform coordinates
    painted_polygons = transformBuildings(polygons, scale, dx, dy, width, height, offset)
   
    return painted_polygons

def transformBuildings(buildings : list[dict], scale : float, dx : float, dy : float, width : int, height : int, offset : int) -> list[Building]:
    """Transform polygons to new coordinates which are compatible with canvas coordinates
        
        Parameters
        ----------
        buildings : list[dict]
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
        painted polygons : list[Building]
            list of buildings with transformed coordinates
    """
    painted_buildings = []    
    
    # Iterate through all buildings
    for region in buildings:
        pol = Building()
        points = region['POINTS']
        if len(points) < 3:
            continue
        
        #Iterate through all points in building
        for point in points:
            if point[0] == 'N' or point[1] == 'N' or point[0] == 'P' or point[1] == 'P':
                continue
            
            # Scale and move original coordinates
            new_x = (float(point[0]) - dx)/scale + width/2
            new_y = (-1)*((float(point[1]) - dy)/scale + height/2) + height - offset # In Qt origin of the coordianate system is on the upper left corner
            
            point = QPointF(new_x, new_y)   
            pol.addVertex(point)            # Save new point in building
        painted_buildings.append(pol)

    return painted_buildings