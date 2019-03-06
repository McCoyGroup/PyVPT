from .CoordinateSystemConverter import CoordinateSystemConverter
from .CommonCoordinateSystems import CartesianCoordinates3D, ZMatrixCoordinates
from ..CoordinateTransformations.TransformationUtilities import vec_norms, vec_angles, pts_dihedrals, vec_normalize, vec_crosses
from ..CoordinateTransformations.TransformationUtilities import affine_matrix, merge_transformation_matrix, rotation_matrix
import numpy as np
# this import gets bound at load time, so unfortunately PyCharm can't know just yet
# what properties its class will have and will try to claim that the files don't exist



class ZMatrixToCartesianConverter(CoordinateSystemConverter):
    """A converter class for going from ZMatrix coordinates to Cartesian coordinates

    """

    @property
    def types(self):
        return (ZMatrixCoordinates, CartesianCoordinates3D)

    @staticmethod
    def zmatrix_affine_transforms(centers, vecs1, vecs2, angles, dihedrals):
        """Builds a single set of affine transformation matrices to apply to the vecs1 to get the next set of points

        :param refs1:
        :type refs1:
        :param refs2:
        :type refs2:
        :param refs3:
        :type refs3:
        :param angles:
        :type angles:
        :param dihedrals:
        :type dihedrals:
        :return:
        :rtype:
        """
        crosses = vec_crosses(vecs1, vecs2)
        rot_mats_1 = rotation_matrix(crosses, angles)
        rot_mats_2 = rotation_matrix(vecs2, dihedrals)
        transfs = affine_matrix(np.dot(rot_mats_1, rot_mats_2), centers)
        return transfs

    def build_next_points(self, refs1, dists, refs2, angles, refs3, dihedrals):
        vecs1 = dists*vec_normalize(refs2 - refs1)
        vecs2 = refs3 - refs1
        transfs = self.zmatrix_affine_transforms(refs1, vecs1, vecs2, angles, dihedrals)
        newstuff = np.dot(transfs, vecs1)
        return newstuff[:, :3]


    def convert_many(self, coordlist, origins=None, axes=None, use_rad=True, **kw):
        """Expects to get a list of configurations
        These will look like:
            [
                [anchor, dist, ref, angle, plane, dihedral ]
                ...
            ]
        **For efficiency it is assumed that all configurations have the same length**

        :param coordlist:
        :type coordlist:
        :param origins:
        :type origins:
        :param axes:
        :type axes:
        :param use_rad:
        :type use_rad:
        :param kw:
        :type kw:
        """
        coordnum = len(coordlist)

        if origins is None:
            origins = [0, 0, 0]
        if not isinstance(origins, np.ndarray):
            origins = np.array(origins)
        if len(origins) < coordnum:
            origins = np.repeat(origins, coordnum)

        if axes is None:
            axes = [1, 0, 0]
        if not isinstance(axes, np.ndarray):
            axes = np.array(axes)
        if len(axes) < coordnum:
            axes = np.repeat(axes, coordnum)
        axes = vec_normalize(axes)

        crosses = vec_crosses(axes, np.repeat([0, 1, 0], coordnum)) + np.array([0, 0, 1])

        raise NotImplementedError("Need to finish implementing this...")


    def convert(self, coords, **kw):
        """dipatches to convert_many but only pulls the first"""
        return self.convert_many(np.reshape(coords, (1,)+coords.shape), **kw)[0]


__converters__ = [ ZMatrixToCartesianConverter() ]