
import numpy as np, scipy.sparse as sparse

class FchkForceConstants:
    """Holder class for force constants coming out of an fchk file
    Allows us to construct the force constant matrix in lazy fashion if we want

    """
    def __init__(self, fcs):
        self.fcs = fcs
        self._n = None
    def __len__(self):
        return len(self.fcs)
    def _get_n(self):
        """

        :return:
        :rtype: int
        """
        if self._n is None:
            self._n = int((-1 + np.sqrt(1 + 8*len(self)))/6) # solving 3n*3n == 2*l - 3n
        return self._n
    @property
    def n(self):
        return self._get_n()
    def _get_array(self):
        """Uses the computed n to make and symmetrize an appropriately formatted array from the lower-triangle data

        :return:
        :rtype: np.ndarray
        """
        n = self.n
        full_array = np.zeros((3*n, 3*n))
        full_array[np.tril_indices_from(full_array)] = self.fcs
        full_array = full_array + np.tril(full_array, 1).T
        return full_array
    @property
    def array(self):
        return self._get_array()

class FchkForceDerivatives:
    """Holder class for force constant derivatives coming out of an fchk file"""
    def __init__(self, derivs):
        self.derivs = derivs
        self._n = None
    def __len__(self):
        return len(self.derivs)
    def _get_n(self):
        if self._n is None:
            l = len(self)
            # had to use Mathematica to get this from the cubic poly
            #  2*(3n-6)*(3n)^2 == 2*l - 2*(3n-6)*(3n)
            l_quad = 81*l**2 + 3120*l - 5292
            l_body = (3*np.sqrt(l_quad) - 27*l - 520)
            if l_body > 0:
                l1 = l_body**(1/3)
            else:
                l1 = -(-l_body)**(1/3)
            n = (1/18)*( 10 + (2**(1/3))*( l1 - 86/l1) )
            self._n = np.ceil(n) # precision issues screw this up in python, but not in Mathematica (I think)
        return self._n
    @property
    def n(self):
        return self._get_n()
    def _get_third_derivs(self):
        # fourth and third derivs are same len
        d = self.derivs
        return d[:int(len(d)/2)]
    def _get_fourth_derivs(self):
        # fourth and third derivs are same len
        d = self.derivs
        return d[int(len(d)/2):]
    @property
    def third_derivs(self):
        return self._get_third_derivs()
    @property
    def fourth_derivs(self):
        return self._get_fourth_derivs()

    def _get_third_deriv_array(self):
        """we make a (3n-6)(3n)^2 X (3n-6)(3n)^2 array to fill from the third_derivs
        and then symmetrize it

        :return:
        :rtype: np.ndarray
        """
        n = self.n
        derivs = self.third_derivs
        dim = (3*n-6)*(3*n)**2
        full_array_1 = np.zeros((dim, dim))
        inds_1 = np.tril_indices_from(full_array_1)
        full_array_1[inds_1] = derivs
        full_array_1 = full_array_1 + np.tril(full_array_1, 1).T
        return full_array_1
    @property
    def third_deriv_array(self):
        return self._get_third_deriv_array()

    def _get_fourth_deriv_array(self):
        """We'll make our array of fourth derivs a sparse matrix with the 3n

        :return:
        :rtype: sparse.csr_matrix
        """
        # do the same for the diagonal 4th derivs
        n = self.n
        fds = self.fourth_derivs
        dim1 =(3*n-6)*(3*n)**2
        dim = (3*n-6)*dim1

        inds = np.tril_indices(dim1) # generates a tuple by default I think...
        inds_mat_0 = np.broadcast(inds[0], (3*n-6, np.shape(inds[0])))
        inds_mat_1 = np.broadcast(inds[1], (3*n-6, np.shape(inds[1])))
        inds = (
            np.reshape(inds_mat_0, (dim,)),
            np.reshape(inds_mat_1, (dim,))
        )
        fdarray = sparse.csr_matrix((), shape=((dim, dim)))
        fdarray[inds] = fds # hopefully this works...
        fdarray = fdarray + np.tril(fdarray, 1).T
        return fdarray

    @property
    def fourth_deriv_array(self):
        return self._get_fourth_deriv_array()