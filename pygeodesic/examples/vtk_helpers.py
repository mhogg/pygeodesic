import vtk
from vtkmodules.util import numpy_support as nps
from vtkmodules.numpy_interface import dataset_adapter as dsa
import numpy as np


def getPointsAndCellsFromPolydata(polydata):
    """
    Extract points and cells from polydata object

    Args:
        polydata (vtk.vtkPolyData): polydata object representing mesh
    
    Returns:
        ndarray: Array of points
        ndarrdy: Array of cells
    """    
    
    # Use numpy wrapper to easily extract points and cells
    polydata = dsa.WrapDataObject(polydata) if isinstance(polydata, vtk.vtkPolyData) else polydata

    # Get points
    points = np.array(polydata.GetPoints(), dtype=np.float64)

    # Get cells
    polygons = np.array(polydata.GetPolygons(), dtype=np.int32)
    n = polygons[0]+1
    polygons = np.resize(polygons, (polygons.size//n, n))
    polygons = polygons[:, 1:n]
    
    return points, polygons


def polydataFromPointsAndCells(points, cells):
    """
    Create new vtkPolyData object from provided points and cells

    Args:
        points (list or ndarray): 3D array of point coordinates (float type)
        cells (list or numpy array): 3D or 4D array of point connectivity for 
            each cell (int type)
    
    Returns:
        vtk.vtkPolyData: polydata representation of mesh points and cells
    """
    # Get vtk point array
    points = np.array(points)
    vtkpoints = vtk.vtkPoints()
    vtkpoints.SetData(nps.numpy_to_vtk(points, deep=True))
    # Get vtk cell array
    cells  = np.array(cells)
    number_of_cells, points_per_cell = cells.shape
    points_per_cell_array = np.full((number_of_cells,1), points_per_cell, dtype=np.int32)
    cellArray = np.hstack((points_per_cell_array, cells)).flatten()
    vtkcells = vtk.vtkCellArray()
    vtkcells.SetCells(number_of_cells, nps.numpy_to_vtk(cellArray, deep=True,
                        array_type=vtk.vtkIdTypeArray().GetDataType()))
    # Put point and cell data into a vtkPolyData object
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtkpoints)
    polydata.SetPolys(vtkcells)
    
    return polydata


def createPolyDataActor(polydata, color=(1,1,1), opacity=1.0):
    """
    Create and return a vtkActor from vtkPolyData object
    
    Args:
        polydata (vtk.vtkPolyData): polydata object representing the mesh
        color (list[int]): RBG mesh color 
        opacity (float): Value representing opacity of actor [0-1]
    
    Returns:
        vtk.vtkActor
    """
    mapper = vtk.vtkPolyDataMapper() 
    mapper.SetInputData(polydata) 
    actor = vtk.vtkActor() 
    actor.SetMapper(mapper) 
    actor.GetProperty().SetColor(color)
    actor.GetProperty().SetOpacity(opacity)
    return actor


def createPolyLineActor(pointCoords, linewidth=3, color=(1,1,1)):
    """
    Create and return a vtkActor from a list of point coordinates
    
    Args:
        pointCoords (ndarray): array of point coordinates along line
        linewidth (float): Value representing width of line actor (default=3)
        color (list[int]): RBG mesh color  
    
    Returns:
        vtk.vtkActor
    """
    # Points
    points = vtk.vtkPoints()
    for p in pointCoords:
        points.InsertNextPoint(p)
    # Lines
    lines = vtk.vtkCellArray()
    numLines = len(pointCoords)-1
    for i in range(numLines):
        lines.InsertNextCell(2)
        lines.InsertCellPoint(i)
        lines.InsertCellPoint(i+1)
    # Polydata from points and lines
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)    
    # Create and return actor
    mapper = vtk.vtkPolyDataMapper() 
    mapper.SetInputData(polydata) 
    actor = vtk.vtkActor() 
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)
    actor.GetProperty().SetLineWidth(linewidth)
    return actor


def createSphereActor(center, radius=1.5, color=(1,1,1), resolution=20):
    """
    Create and return a vtkActor based on vtkSphereSource
    
    Args:
        center (ndarray): array of center coordinates
        radius (float): value representing radius of the sphere
        color (list[int]): RBG mesh color      
        resolution (int): phi and theta mesh resolution of the sphere
    
    Returns:
        vtk.vtkActor    
    """
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(center)
    sphere.SetRadius(radius)
    sphere.SetPhiResolution(resolution)
    sphere.SetThetaResolution(resolution)
    sphereMapper = vtk.vtkPolyDataMapper()
    sphereMapper.SetInputConnection(sphere.GetOutputPort())
    sphereActor = vtk.vtkActor()
    sphereActor.GetProperty().SetColor(color)  
    sphereActor.SetMapper(sphereMapper)
    return sphereActor


class Viewer:
    """
    Convience class for creating a VTK render window and adding actors

    Attributes:
        backgroundColor (List[float]): RGBF color of renderwindow background
        windowSize (List[int]): Size of VTK renderwindow

    Methods:
        addActor(actor)
        addActors(actors)
        removeActor(actor)
        removeAllActors()
        show()

    Examples:
        >>> viewer = Viewer()
        >>> viewer.addActors([actor1, actor2, ])
        >>> viewer.show()
    """

    def __init__(self, windowSize=(800,600), backgroundColor=(82, 87, 110)):
        self.actors = []
        self._windowSize = windowSize
        self._backgroundColor = [v/255. for v in backgroundColor]

    def setup(self):
        # Renderers
        # Renderer 1 - For normal actors
        self.renderer1 = vtk.vtkRenderer()
        self.renderer1.SetBackground(*self._backgroundColor)
        self.renderer1.SetLayer(0)
        self.renderer1.SetUseDepthPeeling(True)
        # Render Window
        self.render_window = vtk.vtkRenderWindow()
        self.render_window.AddRenderer(self.renderer1)
        self.render_window.SetSize(self.windowSize)
        self.setWindowName('VTK viewer')
        # Interactor
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetRenderWindow(self.render_window)
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.iren.SetInteractorStyle(style)
        # Add actors
        self.addActorsToRenderWindow()
        self.renderer1.ResetCamera()
        self.render_window.Render()

    def addActorsToRenderWindow(self):
        """
        Add actors to render window
        For vtkFollower objects (i.e. landmark labels), also need to set the camera
        """
        for actor in self.actors:
            if isinstance(actor, vtk.vtkFollower):
                actor.SetCamera(self.renderer2.GetActiveCamera())
                self.renderer2.AddActor(actor)
            else:
                self.renderer1.AddActor(actor)

    def setWindowName(self, windowName):
        """
        Set window name
        
        Args:
            windowName (str)
        """
        # Must call render before changing window name
        # Ref: https://lorensen.github.io/VTKExamples/site/Python/Visualization/WindowTitle/
        self.render_window.Render()
        self.render_window.SetWindowName(windowName)

    def show(self):
        """
        Call to display the VTK render window
        """
        self.setup()
        self.iren.Initialize()
        self.iren.Start()

    @property
    def windowSize(self):
        """
        windowSize (List[int]): Size (W, H) of the VTK renderwindow
        """
        return self._windowSize

    @windowSize.setter
    def windowSize(self, windowSize):
        self._windowSize = windowSize

    def addActor(self, actor):
        """
        Add single VTK actor

        Args:
            actor (vtk.vtkActor)
        """
        self.actors.append(actor)

    def addActors(self, actors):
        """
        Add multiple VTK actors

        Args:
            actors (List[vtk.vtkActor])
        """
        for actor in actors:
            self.actors.append(actor)

    def removeActor(self, actor):
        """
        Remove actor
        
        Args:
            actor (vtk.vtkActor)
        """
        try:
            self.actors.remove(actor)
        except ValueError as e:
            print('Error in removeActor:', e)

    def removeAllActors(self):
        """
        Remove all actors
        """
        self.actors = []