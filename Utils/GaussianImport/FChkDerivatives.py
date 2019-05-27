
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
            self._n = int(np.ceil(n)) # precision issues screw this up in python, but not in Mathematica (I think)
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

    @staticmethod
    def _fill_3d_tensor(n, derivs):
        """Makes and fills a 3D tensor for our derivatives

        :param n:
        :type n:
        :param derivs:
        :type derivs:
        :return:
        :rtype: np.ndarray
        """
        dim_1 = (3*n)
        mode_n = 3*n-6
        full_array_1 = np.zeros((mode_n, dim_1, dim_1))
        # set the lower triangle
        inds_1 = np.tril_indices_from(full_array_1[0])
        l_per = len(inds_1[0])
        main_ind = np.broadcast_to(np.arange(mode_n), (l_per, mode_n)).flatten(order="F")
        sub_ind_1 = np.broadcast_to(inds_1[0], (mode_n, l_per)).flatten()
        sub_ind_2 = np.broadcast_to(inds_1[1], (mode_n, l_per)).flatten()
        inds = ( main_ind, sub_ind_1, sub_ind_2 )
        full_array_1[inds] = derivs
        # set the upper triangle
        inds2 = ( main_ind, sub_ind_2, sub_ind_1 ) # basically just taking a transpose
        full_array_1[inds2] = derivs
        return full_array_1

    def _get_third_deriv_array(self):
        """we make the appropriate 3D tensor from a bunch of 2D tensors

        :return:
        :rtype: np.ndarray
        """
        n = self.n
        derivs = self.third_derivs
        return self._fill_3d_tensor(n, derivs)

    @property
    def third_deriv_array(self):
        return self._get_third_deriv_array()

    def _get_fourth_deriv_array(self):
        """We'll make our array of fourth derivs exactly the same as the third
        admittedly this should be a 4D tensor, but we only have the diagonal elements so it's just 3D

        :return:
        :rtype: sparse.csr_matrix
        """
        n = self.n
        derivs = self.fourth_derivs
        return self._fill_3d_tensor(n, derivs)

    @property
    def fourth_deriv_array(self):
        return self._get_fourth_deriv_array()