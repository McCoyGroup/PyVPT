"""The coordinate transform class defines an architecture to transform coordinates

"""
from ..Common import *

######################################################################################################
##
##                                   CoordinateTransform Class
##
######################################################################################################

class CoordinateTransform:
    """The CoordinateTransform class provides a simple, general way to represent a
    compound coordinate transformation

    :transforms:
    """

    def __init__(self, transforms):
        self.transform_list = [self.parse_transform(tf) for tf in transforms]

    @staticmethod
    def parse_transform(tf):
        pass

    @staticmethod
    def matrix_transform(mat):
        return np.array(mat)

