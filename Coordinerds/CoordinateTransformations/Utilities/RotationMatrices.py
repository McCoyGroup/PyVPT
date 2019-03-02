
import math, numpy as np

def rotation_matrix_basic(xyz, theta):
    """rotation matrix about x, y, or z axis

    :param xyz: x, y, or z axis
    :type xyz: str
    :param theta: counter clockwise angle in radians
    :type theta: float
    """

    axis = xyz.lower()
    if axis == "z": # most common case so it comes first
        mat = [
            [ math.cos(theta), math.sin(theta), 0.],
            [-math.sin(theta), math.cos(theta), 0.],
            [0.,               0.,              1.]
        ]
    elif axis == "y":
        mat = [
            [ math.cos(theta), 0., math.sin(theta)],
            [0.,               1.,              0.],
            [-math.sin(theta), 0., math.cos(theta)]
        ]
    elif axis == "x":
        mat = [
            [1.,               0.,              0.],
            [0.,  math.cos(theta), math.sin(theta)],
            [0., -math.sin(theta), math.cos(theta)]
        ]
    else:
        raise Exception("{}: axis '{}' invalid".format('rotation_matrix_basic', xyz))
    return np.array(mat)

#thank you SE for the nice Euler-Rodrigues imp: https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
def rotation_matrix_ER(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([
        [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]
    ])

def rotation_matrix(axis, theta):
    if type(axis) == str:
        mat_fun = rotation_matrix_basic
    else:
        mat_fun = rotation_matrix_ER

    return mat_fun(axis, theta)