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
#       vec_crosses
#

def vec_crosses(vecs1, vecs2):
    return np.cross(vecs1, vecs2)

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






