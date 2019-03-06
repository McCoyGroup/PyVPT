from unittest import *
from .DataGenerator import DataGenerator
from ..Coordinerds.CoordinateSystems import *

class ConverterTest(TestCase):
    def test_CoordinateSet(self):
        import numpy as np
        coord_set = CoordinateSet(DataGenerator.coords(500))
        self.assertIsInstance(coord_set.coords, np.ndarray)

    def test_Loader(self):
        loaded = CoordinateSystemConverters.get_converter(CartesianCoordinates3D, ZMatrixCoordinates)
        self.assertIsInstance(loaded, CoordinateSystemConverter)

    def test_CartesianToZMatrix(self):
        coord_set = CoordinateSet(DataGenerator.coords(10))
        coord_set = coord_set.convert(ZMatrixCoordinates, use_rad = False)
        self.assertEqual(coord_set.coords.shape, (9, 6))

    def test_CartesianToZMatrixMulti(self):
        coord_set = CoordinateSet(DataGenerator.multicoords(10, 10))
        coord_set = coord_set.convert(ZMatrixCoordinates, use_rad = False)
        self.assertEqual(coord_set.coords.shape, (10, 9, 6))