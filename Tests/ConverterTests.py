
from .TestUtils import *
from ..Coordinerds.CoordinateSystems import *

class ConverterTest(TestCase):

    def setUp(self):
        self.n = 1000
        self.test_zmats = CoordinateSet(DataGenerator.zmats(self.n, 15), system=ZMatrixCoordinates)
        self.test_carts = CoordinateSet(DataGenerator.multicoords(self.n, 10))

    @validationTest
    def test_CoordinateSet(self):
        import numpy as np
        coord_set = CoordinateSet(DataGenerator.coords(500))
        self.assertIsInstance(coord_set.coords, np.ndarray)
    @validationTest
    def test_Loader(self):
        loaded = CoordinateSystemConverters.get_converter(CartesianCoordinates3D, ZMatrixCoordinates)
        self.assertIsInstance(loaded, CoordinateSystemConverter)
    @validationTest
    def test_CartesianToZMatrix(self):
        coord_set = CoordinateSet(DataGenerator.coords(10))
        coord_set = coord_set.convert(ZMatrixCoordinates, use_rad = False)
        self.assertEqual(coord_set.coords.shape, (9, 6))
    @validationTest
    def test_CartesianToZMatrixMulti(self):
        coord_set = CoordinateSet(DataGenerator.multicoords(10, 10))
        coord_set = coord_set.convert(ZMatrixCoordinates, use_rad = False)
        self.assertEqual(coord_set.coords.shape, (10, 9, 6))
    @validationTest
    def test_ZMatrixToCartesian(self):
        coords = self.test_zmats.convert(CartesianCoordinates3D, use_rad = False)
        self.assertEqual(coords.coords.shape, (self.n, 16, 3))
    @timeitTest(number=100)
    def test_CartToZTiming(self):
        coord_set = self.test_carts
        coord_set = coord_set.convert(ZMatrixCoordinates, use_rad = False)
        self.assertEqual(coord_set.coords.shape, (self.n, 9, 6))
    @timeitTest(number=100)
    def test_ZMatrixToCartesianTiming(self):
        coords = self.test_zmats.convert(CartesianCoordinates3D, use_rad = False)
        self.assertEqual(coords.coords.shape, (self.n, 16, 3))