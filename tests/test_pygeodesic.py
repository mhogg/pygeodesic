"""Basic unit tests for the pygeodesic library."""
import os
import pathlib
import pygeodesic
import pygeodesic.geodesic as geodesic
import numpy as np
import pytest


THIS_PATH = pathlib.Path(__file__).parent.resolve()
DATA_PATH = os.path.join(THIS_PATH, "data")


def test_version():
    assert isinstance(pygeodesic.__version__, str)


@pytest.fixture()
def geoalg():
    # Read the mesh to get the points and faces of the mesh
    filename = os.path.join(DATA_PATH, "flat_triangular_mesh.txt")
    points, faces = geodesic.read_mesh_from_file(filename)
    return geodesic.PyGeodesicAlgorithmExact(points, faces)


def test_distance(geoalg):
    source_index = 25
    target_index = 97
    distance, path = geoalg.geodesicDistance(source_index, target_index)
    assert np.isclose(distance, 1.697056274847714)

    known_path = np.array(
        [
            [3.1, 0.8, 0.0],
            [3.0, 0.7, 0.0],
            [2.9, 0.6, 0.0],
            [2.8, 0.5, 0.0],
            [2.7, 0.4, 0.0],
            [2.6, 0.3, 0.0],
            [2.5, 0.2, 0.0],
            [2.4, 0.1, 0.0],
            [2.3, 0.0, 0.0],
            [2.2, -0.1, 0.0],
            [2.1, -0.2, 0.0],
            [2.0, -0.3, 0.0],
            [1.9, -0.4, 0.0],
        ]
    )

    assert np.allclose(known_path, path)


def test_distances(geoalg):
    source_index = 25
    source_indices = np.array([source_index])
    distances, best_source = geoalg.geodesicDistances(source_indices)

    known_distances = np.array(
        [
            0.72111026,
            0.56568542,
            0.4472136,
            0.4,
            0.4472136,
            0.56568542,
            0.72111026,
            0.89442719,
            1.07703296,
            1.26491106,
            1.45602198,
            0.63245553,
            0.4472136,
            0.28284271,
            0.2,
            0.28284271,
            0.4472136,
            0.63245553,
            0.82462113,
            1.0198039,
            1.21655251,
            1.41421356,
            0.6,
            0.4,
            0.2,
            0.0,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            0.63245553,
            0.4472136,
            0.28284271,
            0.2,
            0.28284271,
            0.4472136,
            0.63245553,
            0.82462113,
            1.0198039,
            1.21655251,
            1.41421356,
            0.72111026,
            0.56568542,
            0.4472136,
            0.4,
            0.4472136,
            0.56568542,
            0.72111026,
            0.89442719,
            1.07703296,
            1.26491106,
            1.45602198,
            0.84852814,
            0.72111026,
            0.63245553,
            0.6,
            0.63245553,
            0.72111026,
            0.84852814,
            1.0,
            1.16619038,
            1.34164079,
            1.52315462,
            1.0,
            0.89442719,
            0.82462113,
            0.8,
            0.82462113,
            0.89442719,
            1.0,
            1.13137085,
            1.28062485,
            1.44222051,
            1.61245155,
            1.16619038,
            1.07703296,
            1.0198039,
            1.0,
            1.0198039,
            1.07703296,
            1.16619038,
            1.28062485,
            1.41421356,
            1.56204994,
            1.72046505,
            1.34164079,
            1.26491106,
            1.21655251,
            1.2,
            1.21655251,
            1.26491106,
            1.34164079,
            1.44222051,
            1.56204994,
            1.69705627,
            1.84390889,
            1.52315462,
            1.45602198,
            1.41421356,
            1.4,
            1.41421356,
            1.45602198,
            1.52315462,
            1.61245155,
            1.72046505,
            1.84390889,
            1.97989899,
            1.70880075,
            1.64924225,
            1.61245155,
            1.6,
            1.61245155,
            1.64924225,
            1.70880075,
            1.78885438,
            1.88679623,
            2.0,
            2.12602916,
        ]
    )

    assert np.allclose(known_distances, distances)
