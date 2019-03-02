from .AffineTransform import AffineTransform
import numpy as np

######################################################################################################
##
##                                   RotationTransform Class
##
######################################################################################################

class RotationTransform(AffineTransform):
    """A simple AffineTransform implementation of the TransformationFunction abstract base class

    """

    def __init__(self, theta, axis="z", center = None):
        """

        :param theta: angle through which to rotate
        :type theta: float
        :param axis: axis about which to rotate
        :type axis: axis about which to rotate
        :param center: center point for the rotation
        :type center: None or np.array
        """

        from .Utilities.RotationMatrices import rotation_matrix
        from .Utilities.TranslationMatrices import translation_matrix
        from .Utilities.TransformationTransformations import merge_transformation_mats

        rot_mat = rotation_matrix(theta)
        if center is not None:
            translate = translation_matrix(-center)
            neg_translate = translation_matrix(center)
            rot_mat = merge_transformation_mats(merge_transformation_mats(translate, rot_mat), neg_translate)

        super().__init__(rot_mat, shift=None)

    def reverse(self):
        base_mat = self.transform
        shift = self.shift
        inverse = base_mat.transpose()
        new_shift = -np.dot(inverse, shift)
        type(self)(inverse, shift=new_shift)
