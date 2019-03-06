from unittest import *
import numpy as np
from ..Coordinerds.CoordinateTransformations.TransformationUtilities import *

class DataGenerator:

    @staticmethod
    def coords(n=50):
        return np.random.rand(n, 3)

    @staticmethod
    def mats(n=1):
        return np.random.rand(n, 3, 3)

    @staticmethod
    def vecs(n=1):
        return np.random.rand(n, 3)

class TestAffineMatrix(TestCase):

    def load(self, n=10):
        self.cases = n
        self.transforms = DataGenerator.mats(n)
        self.shifts = DataGenerator.vecs(n)
        self.mats = affine_matrix(self.transforms, self.shifts)

    def test_mat_dim(self):
        self.load()
        self.assertEqual(self.mats.shape, (self.cases, 4, 4))

    def test_mat_mul(self):
        self.load()
        vecs = DataGenerator.vecs(self.cases)
        fuck_this = np.ones((self.cases, 1))
        vecs = np.concatenate((vecs, fuck_this), axis=1)
        vec_prod = mat_vec_muls(self.mats, vecs)
        vec_prod2 = np.array([ a @ b for a, b in zip(self.mats, vecs)])
        self.assertAlmostEqual(np.sum(vec_prod-vec_prod2), 0.)
