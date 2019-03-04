from .CoordinateSystemConverter import CoordinateSystemConverter
from .CommonCoordinateSystems import CartesianCoordinates3D, ZMatrixCoordinates
from .Utilities.VectorMath import vec_norms, vec_angles, pts_dihedrals
import numpy as np
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
                range(ncoords-2, ncoords-2),
                range(ncoords-3, ncoords-2)
            )
        else:
            normalized_list = [None] * len(order_list)
            for i, el in enumerate(order_list):
                if isinstance(el, int):
                    spec = (
                        el,
                        normalized_list[i-1][0] if i > 0 else -1,
                        normalized_list[i-2][0] if i > 1 else -1,
                        normalized_list[i-3][0] if i > 2 else -1
                    )
                else:
                    spec = tuple(el)
                    if len(spec) < 4:
                        raise ValueError("z-matrix conversion spec {} not long enough".format(el))

                normalized_list[i] = spec

        return np.array(normalized_list)

    @staticmethod
    def get_dists(points, centers):
        return vec_norms(centers-points)
    @staticmethod
    def get_angles(lefts, centers, rights):
        # need to look up again what the convention is for which atom is the central one...
        v1s = centers-lefts
        v2s = centers-rights
        return vec_angles(v1s, v2s)[0]
    @staticmethod
    def get_diheds(points, centers, seconds, thirds):
        return pts_dihedrals(centers, seconds, thirds, points)

    def convert(self, coords, ordering=None, use_rad=True, **kw):
        """The ordering should be specified like:

        [
            [n1],
            [n2, n1]
            [n3, n1/n2, n1/n2]
            [n4, n1/n2/n3, n1/n2/n3, n1/n2/n3]
            [n5, ...]
            ...
        ]

        :param coords:    array of cartesian coordinates
        :type coords:     np.ndarray
        :param use_rad:   whether to user radians or not
        :type use_rad:    bool
        :param ordering:  optional ordering parameter for the z-matrix
        :type ordering:   None or tuple of ints or tuple of tuple of ints
        :param kw:        ignored key-word arguments
        :type kw:
        """
        ncoords = len(coords)
        ol = self.canonicalize_order_list(ncoords, ordering)
        om = {
                old:new for old, new in zip(
                    (a[0] for a in ol),
                    range(ncoords)
                )
        }
        targ = [None]*(ncoords-1)
        # need to check against the cases of like 1, 2, 3 atom molecules
        # annoying but not hard
        dists = self.get_dists(
            coords[ol[1:, 0]],
            coords[ol[1:, 1]]
        )
        angles = self.get_angles(
            coords[ol[2:, 0]],
            coords[ol[2:, 1]],
            coords[ol[2:, 2]]
        )
        diheds = self.get_diheds(
            coords[ol[3:, 0]],
            coords[ol[3:, 1]],
            coords[ol[3:, 2]],
            coords[ol[3:, 3]]
        )
        final_coords = np.array(
            [
                [om[ol[1, 1]], dists[0], 0,            0,         0, 0],
                [om[ol[2, 1]], dists[0], om[ol[2, 2]], angles[0], 0, 0]
            ] + [
                [om[o[1]], d, om[o[2]], a, om[o[3]], h] for
                    o, d, a, h in zip( ol[3:], dists[2:], angles[1:], diheds)
            ]
        )
        #### should find some way to return the order, right?
        return final_coords


