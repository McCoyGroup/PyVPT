from .CoordinateSystemConverter import CoordinateSystemConverter
from .CommonCoordinateSystems import CartesianCoordinates3D, ZMatrixCoordinates
# this import gets bound at load time, so unfortunately PyCharm can't know just yet
# what properties its class will have and will try to claim that the files don't exist

class CartesianToZMatrixConverter(CoordinateSystemConverter):
    """A converter class for going from Cartesian coordinates to ZMatrix coordinates

    """

    @property
    def types(self):
        return (CartesianCoordinates3D, ZMatrixCoordinates)

    def canonicalize_order_list(self, ncoords, order_list):
        """Normalizes the way the ZMatrix coordinates are built out

        :param ncoords:
        :type ncoords:
        :param order_list: the basic ordering to apply for the
        :type order_list: iterable or None
        :return:
        :rtype: iterator of int triples
        """
        if order_list is None:
            normalized_list = zip(
                range(ncoords),
                range(-1, ncoords-1),
                range(ncoords-2, ncoords-2)
            )
        else:
            normalized_list = [None] * len(order_list)
            for i, el in enumerate(order_list):
                if isinstance(el, int):
                    triple = (
                        el,
                        normalized_list[i-1] if i > 0 else -1,
                        normalized_list[i-2] if i > 1 else -1
                    )
                else:
                    triple = tuple(el)
                    if len(triple) < 3:
                        raise ValueError("z-matrix conversion triple {} not long enough".format(el))

                normalized_list[i] = triple

        raise normalized_list

    def convert(self, coords, ordering=None, use_rad=True, **kw):
        """

        :param coords:    array of cartesian coordinates
        :type coords:     np.ndarray
        :param use_rad:   whether to user radians or not
        :type use_rad:    bool
        :param ordering:  optional ordering parameter for the z-matrix
        :type ordering:   None or tuple of ints or tuple of tuple of ints
        :param kw:        ignore key-word arguments
        :type kw:
        """
        pass

