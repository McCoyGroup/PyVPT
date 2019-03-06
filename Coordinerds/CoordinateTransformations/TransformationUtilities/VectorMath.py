"""A module of useful math for handling coordinate transformations and things

"""

import numpy as np

################################################
#
#       vec_dots
#

def vec_dot_vanilla(vecs1, vecs2):
    return np.sum(vecs1*vecs2, axis=1)

def vec_dot_matty(vecs1, vecs2):
    return np.diag(np.dot(vecs1, vecs2.T))

def vec_dots(vecs1, vecs2, mode = None, arbitrary_cutoff = 1000000):
    """Computes the pair-wise dot product of two lists of vecs

    # some question as to the fastest way to compute pairwise dot products of a pile of vectors
    # the internet reccomends einsum but this seems to not perform well over huge arrays
    # a definite contender instead is to use an element wise product and a sum
    # the shitty thing is that this two operations, but the latter is O(n) so it shouldn't matter
    # probably fastest is to do a dot and take the diagonal, so that'll be the default method up to a
    # set bytecount where we'll switch to the sum and vec mode
    # alternately could introduce some chunking, but that's for the future
    # the mode argument will someday be used to allow you to control the algorithm

    :param vecs1:
    :type vecs1:
    :param vecs2:
    :type vecs2:
    """
    if vecs1.shape != vecs2.shape:
        raise ValueError("pairwise dot has to be done on things of the same size")
    if vecs1.nbytes > arbitrary_cutoff:
        return vec_dot_vanilla(vecs1, vecs2)
    else:
        return vec_dot_matty(vecs1, vecs2)

################################################
#
#       vec_norms
#

def vec_norms(vecs):
    return np.linalg.norm(vecs, axis=1)

################################################
#
#       vec_normalize
#

def vec_normalize(vecs):
    return vecs/vec_norms(vecs)

################################################
#
#       vec_crosses
#

def vec_crosses(vecs1, vecs2, normalize=False):
    crosses = np.cross(vecs1, vecs2)
    if normalize:
        crosses = crosses/vec_norms(crosses)
    return crosses

################################################
#
#       vec_cos
#

def vec_cos(vectors1, vectors2):
    """Gets the cos of the angle between two vectors

    :param vectors1:
    :type vectors1: np.ndarray
    :param vectors2:
    :type vectors2: np.ndarray
    """
    dots   = vec_dots(vectors1, vectors2)
    norms1 = vec_norms(vectors1)
    norms2 = vec_norms(vectors2)

    return dots/(norms1*norms2)

################################################
#
#       vec_sins
#

def vec_sins(vectors1, vectors2):
    """Gets the sin of the angle between two vectors

    :param vectors1:
    :type vectors1: np.ndarray
    :param vectors2:
    :type vectors2: np.ndarray
    """
    crosses= vec_crosses(vectors1, vectors2)
    norms1 = vec_norms(vectors1)
    norms2 = vec_norms(vectors2)

    return crosses/(norms1*norms2)


################################################
#
#       vec_angles
#

def vec_angles(vectors1, vectors2):
    """Gets the angles and normals between two vectors

    :param vectors1:
    :type vectors1: np.ndarray
    :param vectors2:
    :type vectors2: np.ndarray
    :return: angles and normals between two vectors
    :rtype: (np.ndarray, np.ndarray)
    """
    dots    = vec_dots(vectors1, vectors2)
    crosses = vec_crosses(vectors1, vectors2)
    norms1  = vec_norms(vectors1)
    norms2  = vec_norms(vectors2)
    norm_prod = norms1*norms2
    cos_comps = dots/norm_prod
    sin_comps = crosses/norm_prod

    return (np.arctan2(sin_comps, cos_comps), crosses)


################################################
#
#       pts_normals
#

def pts_normals(pts1, pts2, pts3, normalize=True):
    """Provides the vector normal to the plane of the three points

    :param pts1:
    :type pts1: np.ndarray
    :param pts2:
    :type pts2: np.ndarray
    :param pts3:
    :type pts3: np.ndarray
    :param normalize:
    :type normalize:
    :return:
    :rtype: np.ndarray
    """
    # should I normalize these...?
    return vec_crosses(pts2-pts1, pts3-pts1, normalize=normalize)

################################################
#
#       pts_dihedrals
#

def pts_dihedrals(pts1, pts2, pts3, pts4):
    """Provides the dihedral angle between pts4 and the plane of the other three vectors

    :param pts1:
    :type pts1: np.ndarray
    :param pts2:
    :type pts2: np.ndarray
    :param pts3:
    :type pts3: np.ndarray
    :return:
    :rtype:
    """
    # should I normalize these...?
    normals = pts_normals(pts2, pts1, pts3, normalize=False)
    off_plane_vecs = pts4-pts3
    return vec_angles(off_plane_vecs, normals)[0]


################################################
#
#       mat_vec_muls

def mat_vec_muls(mats, vecs):
    """Pairwise multiplies mats and vecs

    :param mats:
    :type mats:
    :param vecs:
    :type vecs:
    :return:
    :rtype:
    """

    vecs_2 = np.reshape(vecs, vecs.shape + (1,))
    vecs_2 = np.matmul(mats, vecs_2)
    return np.reshape(vecs_2, vecs.shape)






