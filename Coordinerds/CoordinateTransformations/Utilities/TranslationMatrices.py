import numpy as np

def translation_matrix(shift):
    return np.array(
        [
            [1., 0., 0., shift[0]],
            [0., 1., 0., shift[1]],
            [0., 0., 1., shift[2]],
            [0., 0., 0., 1.       ]
        ]
    )