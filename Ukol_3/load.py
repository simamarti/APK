from QPoint3DF import QPoint3DF
from math import inf
import json
from PyQt6 import QtWidgets
from Algorithms import Algorithms

def loadPoints(path : str) -> list[list[QPoint3DF]]:
    """Load points from JSON file
        
        Parameters
        ----------
        path : str
            Relative address of input file
            
        Returns
        -------
        painted polygons : list[QPoint3DF]
            list of points
        
        borderpoints : list[QPoint3DF]
            list of border points
    """

    errorText = ""
    try:
        with open(path, 'r', encoding = 'utf-8') as f:
            file = json.load(f)


        points = file['POINTS']

        dx = 0; dy = 0; scale = 0
        min_x = inf; max_x = -inf; min_y = inf; max_y = -inf 

        # Compute min_max box
        for point in points:
            x = point[0]; y = point[1]; z = point[2]

            # Check if string is float (e. g. -0.123)
            if not isNumber(x) or not isNumber(y) or not isNumber(z):
                continue

            min_x = min(min_x, float(x))
            max_x = max(max_x, float(x))
            min_y = min(min_y, float(y))
            max_y = max(max_y, float(y)) 

    except FileNotFoundError:
        errorText = "File has not been found"
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

    scale = max(1, abs(max_x - min_x)/(width - offset), abs(max_y - min_y)/(height - offset))

    # transform coordinates
    painted_points = transformPoints(points, scale, dx, dy, width, height, offset)
    borderPoints = []

    # Points of border (Delaunay Triangulation of non-Convex region)
    if 'BORDER' in file:
        border = file['BORDER']
        borderPoints = transformPoints(border, scale, dx, dy, width, height, offset)
    else:
        # Set brder points as Convex Hull of points
        alg = Algorithms()
        borderPoints = alg.jarvisScan(painted_points)
        borderPoints = transformPoints(borderPoints, scale, dx, dy, width, height, offset, qpoint3df = True)

    return painted_points, borderPoints

def transformPoints(points : list[list], scale : float, dx : float, dy : float, width : int, height : int, offset : int, qpoint3df = False) -> list[QPoint3DF]:
    """Transform polygons to new coordinates which are compatible with canvas coordinates
        
        Parameters
        ----------
        points : list[list]
            list of list coordinates
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
        qpoint3df : bool (optional)
            input is list of lists or list of QPoint3DF
            
        Returns
        -------
        painted polygons : list[QPoints3DF]
            list of points with transformed coordinates
    """
    painted_points = []    
    for point in points:

        if qpoint3df:
            x = point.x(); y = point.y(); z = point.getZ()
        else:
            x = point[0]; y = point[1]; z = point[2]

        # Check if string is float (e. g. -0.123)
        if not isNumber(x) or not isNumber(y) or not isNumber(z):
            continue

        # Transform coordinates
        new_x = (float(x) - dx)/scale + width/2
        new_y = (-1)*((float(y) - dy)/scale + height/2) + height - offset
        pt = QPoint3DF(new_x, new_y, z)
        painted_points.append(pt)

    return painted_points

def isNumber(string : str | int | float) -> bool:
    """Check if string is number (float e. g. -0.123)
        
        Parameters
        ----------
        string : str | int | float
            String whitch should represent number
            
        Returns
        -------
        True: String is float
        Fase: String is not float
    """
    # Is it float or int?
    if isinstance(string, float) or isinstance(string, int):
        return True

    # Does string consist of only digits (0–9)?
    if string.isdigit():
        return True

    # Check each character if it is number or ./-
    for char in string:
        charAscii = ord(char)
        if (charAscii < 48 or charAscii > 57) and charAscii != 45 and charAscii != 46:
            return False
    return True