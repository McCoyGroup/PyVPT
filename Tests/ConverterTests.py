from unittest import *
from .DataGenerator import DataGenerator
from ..Coordinerds.CoordinateSystems import *

class ConverterTest(TestCase):
    def test_CoordinateSet(self):
        import numpy as np
        coord_set = CoordinateSet(DataGenerator.coords(500))
        self.assertIsInstance(coord_set.coords, np.ndarray)

    def test_Loader(self):
        loaded = CoordinateSystemConverters.load_converter("ZMatrixToCartesian")
        self.assertIsInstance(loaded, CoordinateSystemConverter)

    def test_CartesianToZMatrix(self):
        coord_set = CoordinateSet(DataGenerator.coords(10))
        CoordinateSystemConverters.load_converter("ZMatrixToCartesian")
        coord_set.convert(ZMatrixCoordinates)
        self.assertEqual(coord_set.coords.shape, (10, 6))