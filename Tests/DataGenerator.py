import numpy as np

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