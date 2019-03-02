from unittest import *
import numpy as np
from CoordinateTransformations import *

class DataGenerator:

    @staticmethod
    def get_test_array(n=50):
        return np.random.rand(n, 3)

    @staticmethod
    def get_test_transf():
        return np.random.rand(3, 3)

    @staticmethod
    def get_test_shift():
        return np.random.rand(3)

class TestAffine(TestCase):

    def assetAffineValid(self):
        transf = DataGenerator.get_test_transf()
        shift = DataGenerator.get_test_shift()
        assert isinstance(AffineTransform(transf, shift), TransformationFunction)
