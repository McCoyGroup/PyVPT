from .CoordinateSystemConverter import CoordinateSystemConverter
from .CommonCoordinateSystems import CartesianCoordinates3D, ZMatrixCoordinates
# this import gets bound at load time, so unfortunately PyCharm can't know just yet
# what properties its class will have

class CartesianToZMatrixConverter(CoordinateSystemConverter):
    """A converter class for going from Cartesian coordinates to ZMatrix coordinates

    """

    @property
    def types(self):
        return (CartesianCoordinates3D, ZMatrixCoordinates)

    def canonicalize_order_list(self, ncoords, order_list):
        if order_list is None:
            order_list = zip(range(ncoords), range(-1, ncoords-1), range(n-2, ncoords-2))
        elif all((isinstance(elm, int) for elm in order_list)):
            order_list = ...
        raise NotImplemented

    def convert(self, coords, order_list, **kw):
        raise NotImplemented
