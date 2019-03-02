import numpy as np
from ..CoordinateSystems.CoordinateSystemConverter import CoordinateSystemConverters as converters
from ..CoordinateSystems.CommonCoordinateSystems import CartesianCoordinates3D

######################################################################################################
##
##                                   CoordinateSet Class
##
######################################################################################################

class CoordinateSet:

    def __init__(self, coords, system = CartesianCoordinates3D):
        self.coords = np.asarray(coords)
        self.system = system


    def transform(self, tf):
        """Applies a transformation to the stored coordinates

        :param system:
        :type system:
        :return:
        :rtype:
        """
        coords = self.coords
        new_coords = tf(coords)
        return type(self)(new_coords, self.system)


    def convert(self, system, **kw):
        """Converts across coordinate systems

        :param system:
        :type system:
        :return:
        :rtype:
        """
        cosys = self.system
        converter = converters.get_converter(cosys, system)
        coords = self.coords
        new_coords = converter['converter'](coords, **kw)
        return type(self)(new_coords, system)