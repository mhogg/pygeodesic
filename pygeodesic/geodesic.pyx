# -*- coding: utf-8 -*-
#!python
#cython: language_level=3

import numpy
cimport numpy
from libcpp.vector cimport vector


cdef extern from "geodesic_mesh_elements.h" namespace "geodesic":
    cdef cppclass Face:
        Face()

cdef extern from "geodesic_mesh_elements.h" namespace "geodesic":
    cdef cppclass Vertex:
        Vertex()

cdef extern from "geodesic_mesh_elements.h" namespace "geodesic":
    cdef cppclass SurfacePoint:
        SurfacePoint()
        SurfacePoint(Vertex*)
        double& x()
        double& y()
        double& z()

cdef extern from "geodesic_mesh.h" namespace "geodesic":
    cdef cppclass Mesh:
        Mesh()
        void initialize_mesh_data(vector[double]&, vector[unsigned]&)
        vector[Vertex]& vertices()
        vector[Face]& faces()

cdef extern from "geodesic_algorithm_exact.h" namespace "geodesic":
    cdef cppclass GeodesicAlgorithmExact:
        GeodesicAlgorithmExact(Mesh*)
        void propagate(vector[SurfacePoint]&, double, vector[SurfacePoint]*)
        unsigned best_source(SurfacePoint&, double&)
        void trace_back(SurfacePoint&, vector[SurfacePoint]&)
        void geodesic(SurfacePoint&, SurfacePoint&, vector[SurfacePoint]&)

cdef extern from "geodesic_algorithm_base.h" namespace "geodesic":        
    double length(vector[SurfacePoint]&)

cdef extern from "geodesic_constants_and_simple_functions.h" namespace "geodesic":
    double GEODESIC_INF


cdef class PyGeodesicAlgorithmExact:

    """
    Wrapper around C++ class GeodesicAlgoritmExact

    Usage:
    >>> import pygeodesic.geodesic as geodesic
    >>> geoalg = GeodesicAlgorithmExact(points, faces)
    >>> sourceIndex = 0
    >>> targetIndex = 100
    >>> distance, path = geoalg.geodesicDistance(sourceIndex, targetIndex)
    """

    cdef GeodesicAlgorithmExact *algorithm
    cdef Mesh mesh
    cdef vector[double] points
    cdef vector[unsigned] faces

    def __cinit__(self, numpy.ndarray[numpy.float64_t, ndim=2] _points,
                        numpy.ndarray[numpy.int32_t, ndim=2] _faces):
        """
        Variables:
          points (ndarray, type float64): Array containing the mesh point coordinates
          faces (ndarray, type int32): Array containing the nodal connectivity for each tri face in mesh
        """                        

        cdef numpy.float64_t coord
        for coord in _points.flatten():
            self.points.push_back(coord)

        cdef numpy.int32_t indx
        for indx in _faces.flatten():
            self.faces.push_back(indx)

        self.mesh.initialize_mesh_data(self.points, self.faces)
        self.algorithm = new GeodesicAlgorithmExact(&self.mesh)

    def geodesicDistance(self, int sourceIndex, int targetIndex):
        """
        Calculates the geodesic distance from the mesh vertex with index 'sourceIndex'
        to the mesh vertex with index 'targetIndex'.

        Variables:
          sourceIndex (int): index of source vertex in mesh points array
          targetIndex (int): index of target vertex in mesh points array
          
        Returns:
          path_length (double): the geodesic distance from the source vertex to the target vertex
          path_points (ndarray): the coordinates of the points that make up the path
        """
        
        cdef Py_ssize_t i       
        cdef vector[SurfacePoint] path

        # Reset source, target indices to within limits
        def resetIndicesToWithinLimits(int indx):
            cdef minIndex = 0
            cdef maxIndex = self.mesh.vertices().size()-1        
            if indx < minIndex: indx = minIndex
            if indx > maxIndex: indx = maxIndex
            return indx

        sourceIndex = resetIndicesToWithinLimits(sourceIndex)
        targetIndex = resetIndicesToWithinLimits(targetIndex)

        cdef SurfacePoint source = SurfacePoint(&self.mesh.vertices()[sourceIndex])
        cdef SurfacePoint target = SurfacePoint(&self.mesh.vertices()[targetIndex])
        self.algorithm.geodesic(source, target, path)
        
        cdef list path_points = []
        for i in range(path.size()):
            path_points.append([path[i].x(), path[i].y(), path[i].z()])

        cdef double path_length = length(path)

        return path_length, numpy.array(path_points)

    def geodesicDistances(self, numpy.ndarray[numpy.int32_t, ndim=1] source_indices=None,
                          numpy.ndarray[numpy.int32_t, ndim=1] target_indices=None,
                          double max_distance = GEODESIC_INF):
        
        """
        Calculates the distance of each target vertex from the best (closest) source vertex

        Variables:
          - source_indices (ndarray, type int32): array containing the indices of all the source vertices in the mesh
          - target_indices (ndarray, type int32): array containing the indices of all the target vertices in the mesh
        
        Returns:
          distances (ndarray): the geodesic distance of the target vertex to the closest source
          best_sources (ndarray): the index of the closest source with respect to source_indices       
        """

        cdef Py_ssize_t i
        cdef vector[SurfacePoint] all_sources
        cdef vector[SurfacePoint] stop_points
        cdef numpy.ndarray[numpy.float64_t, ndim=1] distances

        # Setup sources. Defaults to single vertex with index = 0 if None provided
        if source_indices is None:
            source_indices = numpy.arange(0, dtype=numpy.int32)
        for i in source_indices:
            all_sources.push_back(SurfacePoint(&self.mesh.vertices()[i]))
        
        # Setup targets. Defaults to all vertices if None provided
        if target_indices is None:
            for i in range(self.mesh.vertices().size()):
                stop_points.push_back(SurfacePoint(&self.mesh.vertices()[i]))
            self.algorithm.propagate(all_sources, max_distance, NULL)
        else:
            for i in target_indices:
                stop_points.push_back(SurfacePoint(&self.mesh.vertices()[i]))
            self.algorithm.propagate(all_sources, max_distance, &stop_points)
        
        # Calculate distance of each target from the best (closest) source
        distances   = numpy.zeros((stop_points.size(), ), dtype=numpy.float64)
        best_source = numpy.zeros((stop_points.size(), ), dtype=numpy.int32)
        for i in range(stop_points.size()):
            best_source[i] = self.algorithm.best_source(stop_points[i], distances[i])
        distances[distances==GEODESIC_INF] = numpy.inf

        return distances, best_source
    
    def __dealloc__(self):
        del self.algorithm


def read_mesh_from_file(filename):
    """
    Read mesh from example files

    Variables:
        filename (str): File name of mesh file
    """
    
    points = []; faces = []
    with open(filename,'r') as f:
        
        # Read header line for number of points and faces
        try:
            vals = f.readline().strip().split(' ')
            vals = [int(v) for v in vals]
            numpoints, numfaces = vals
            assert numpoints >= 3, "Number of points not >= 3"
        except AssertionError as e:
            print(f'Error reading header: {e}')
            return
        except Exception as e:
            print(f'Error reading header: {e}')
            return
        
        # Read points
        try:
            for i in range(numpoints):
                vals = f.readline().strip().replace('\t',' ').split(' ')
                vals = [float(v) for v in vals]
                points.append(vals)
            points = numpy.array(points)
        except Exception as e:
            print(f'Error reading points: {e}')
            return
        
        # Read faces
        try:
            for i in range(numfaces):
                vals = f.readline().strip().replace('\t',' ').split(' ')
                vals = [int(v) for v in vals]
                faces.append(vals)
            faces = numpy.array(faces)
        except Exception as e:
            print(f'Error reading faces: {e}')
            return            

    return points, faces
