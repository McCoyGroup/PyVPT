import numpy as np, scipy.sparse as sp
from Psience import Wavefunction
from McUtils.Zachary import SparseTensor

__all__ = [
    'PerturbationTheoryWavefunction',
    'PerturbationTheoryException',
    'PerturbationTheoryHamiltonian'
]

class PerturbationTheoryException(Exception):
    pass

class PerturbationTheoryWavefunction(Wavefunction):

    def __init__(self, hamiltonian):

        super().__init__(None, hamiltonian)

    @classmethod
    def from_params(cls, *params):
        """Creates a PerturbationTheoryHamiltonian and uses that to get the terms we want"""

        raise NotImplemented

    def _perturbtation_wavefunction_coeffs(self, states, energies, expect_1):
        """Computes raw perturbation expansion of the wavefunctions

        :param states: the energies to target
        :type states: iterable[int]
        :param energies: the energies for the basis states
        :type energies: iterable[float]
        :param expect_1: the function to compute the expectation values of the first order correction to the Hamiltonian
        :type expect_1: function
        :return: coeffs
        :rtype: np.ndarray
        """

        return np.array(
                 expect_1(n, np.delete(np.arange(len(energies)), n)) /
                    energies[n] - np.delete(energies, n) for n in states
        )


    def _pertubation_energies(self,
                            states,
                            omega,
                            energies,
                            expect_1,
                            expect_2,
                            hbar = 1
                            ):
        """Computes raw perturbation theory energies

        :param states: the energies to target
        :type states: iterable[int]
        :param omega: the base harmonic frequency
        :type omega: float
        :param energies: the energies for the basis states
        :type energies: iterable[float]
        :param expect_1: the function to compute the expectation values of the first order correction to the Hamiltonian
        :type expect_1: function
        :param expect_2: the function to compute the expectation values of the first order correction to the Hamiltonian
        :type expect_2:
        :param hbar: the hbar to use -- generally probably want this to be in atomic units ?
        :type hbar:
        :return: energies
        :rtype: np.ndarray
        """

        en = [
            (
                hbar * omega * ( n + 1/2 ) +
                np.sum(
                    np.power(
                        expect_1(n, np.delete(np.arange(len(energies)), n)),
                        2
                    ) / (energies[n] - np.delete(energies, n))
                ) +
                expect_2(n, n)
            ) for n in states
        ]

        return np.array(en)

class PerturbationTheoryHamiltonian:
    """Represents the main handler used in the perturbation theory calculation

    :param coords: The current coordinates of the system (used for computing the Jacobians)
    :type coords:
    :param pot_derivs: The first, second, third, and fourth derivatives of the potential
    :type pot_derivs:
    :param modes: The cartesian coordinate normal modes
    :type modes: CoordinateSystem
    :param internals: The internal coordinate normal modes
    :type internals:
    """

    def __init__(self,
                 coords,
                 pot_derivs,
                 modes,
                 internals
                 ):
        self.coords = coords
        self.coords_n = (len(modes) + 6)
        self.mode_n = len(modes)
        self.grad, self.fcs, self.thirds, self.fourths = self._canonicalize_derivs(pot_derivs)
        self.modes = modes
        self.internals = internals

    def _canonicalize_derivs(self, derivs):

        try:
            grad, fcs, thirds, fourths = derivs
        except:
            grad, fcs, derivs = derivs
            thirds = derivs.third_deriv_array
            fourths = derivs.fourth_deriv_array

        coord_n = self.coords_n
        modes_n = self.mode_n
        if grad.shape != (coord_n,):
            raise PerturbationTheoryException("{}.{}: length of gradient array ({2[0]}) is not {3[0]}}",
                                              type(self).__name__,
                                          "_canonicalize_force_constants",
                                              grad.shape,
                                              (coord_n,)
                                              )
        if fcs.shape != (coord_n, coord_n):
            raise PerturbationTheoryException("{}.{}: dimension of force constant array ({2[0]}x{2[1]}) is not {3[0]}x{3[1]}",
                                              type(self).__name__,
                                          "_canonicalize_force_constants",
                                              fcs.shape,
                                              (coord_n, coord_n)
                                              )
        if thirds.shape != (modes_n, coord_n, coord_n):
            raise PerturbationTheoryException("{}.{}: dimension of third derivative array ({2[0]}x{2[1]}x{2[2]}) is not {3[0]}x{3[1]}x{3[2]}",
                                              type(self).__name__,
                                          "_canonicalize_derivs",
                                              thirds.shape,
                                              (modes_n, coord_n, coord_n)
                                              )
        # this might need to change in the future
        if fourths.shape != (modes_n, coord_n, coord_n):
            raise PerturbationTheoryException("{}.{}: dimension of fourth derivative array ({2[0]}x{2[1]}x{2[2]}) is not {3[0]}x{3[1]}x{3[2]}",
                                              type(self).__name__,
                                          "_canonicalize_derivs",
                                              fourths.shape,
                                              (modes_n, coord_n, coord_n)
                                              )

        return grad, fcs, thirds, fourths

    def pmatrix_ho(self, n): # we'll depend on the coefficient to be handled later
        """

        :param n:
        :type n:
        :return:
        :rtype: sp.csr_matrix
        """

        ar = np.sqrt(np.arange(1, n))
        bands = [
            [ - ar,  1 ],
            [   ar, -1 ]
        ]
        return sp.diags( [ b[0] for b in bands], [ b[1] for b in bands ] )
    def qmatrix_ho(self, n): # again, we'll depend on the coefficient to be handled later
        """

        :param n:
        :type n:
        :return:
        :rtype: sp.csr_matrix
        """

        ar = np.sqrt(np.arange(1, n))
        bands = [
            [ ar,  1 ],
            [ ar, -1 ]
        ]
        return sp.diags( [ b[0] for b in bands], [ b[1] for b in bands ] )

    def _operator_submatrix(self, funcs, dims, totinds, inds):
        """Returns the operator submatrix for a product operator like piQjpk or whatever

        :param funcs: the functions that take a dimension size and return a matrix for the suboperator
        :type funcs:
        :param dims: dimensions of each coordinate (e.g. (5, 8, 2, 9))
        :type dims: np.ndarray
        :param totinds: the total list of indices of which inds is a subset
        :type totinds: np.ndarray
        :param inds: the list of indices
        :type inds: np.ndarray
        :return:
        :rtype:
        """
        import functools as ft

        pieces = [ None ] * len(totinds)
        for f, i in zip(funcs, inds):
            if pieces[i] is None:
                pieces[i] = f(dims[i])
            else:
                pieces[i] = pieces[i].dot(f(dims[i]))

        for j in np.setdiff1d(totinds, inds):
            pieces[j] = sp.identity(dims[j])

        return ft.reduce(sp.kron, pieces)

    def product_operator_tensor(self, funcs, dims):
        """Generates the tensor created from the product of funcs over the dimensions dims

        :param funcs:
        :type funcs:
        :param dims:
        :type dims:
        :return:
        :rtype:
        """

        ndim = (len(dims),)*len(funcs)
        base_tensor = np.transpose(np.indices(ndim, dtype=object), np.roll(np.arange(len(funcs))))
        tot_inds = np.arange(len(dims))
        news_boy = lambda inds, ti=tot_inds, f=funcs, d=dims: self._operator_submatrix(f, d, ti, inds)
        news_boys = np.apply_along_axis(news_boy, -1, base_tensor)

        return SparseTensor(news_boys)

    def QQQ(self, coefficients, dimensions, qmatrix = None):
        #TODO: get the coefficients in there...
        if qmatrix is None:
            qmatrix = self.qmatrix_ho

        return self.product_operator_tensor((qmatrix, qmatrix, qmatrix), dimensions)

    # TODO: gotta actually include the coefficients somehow...?
    def pQp(self,
            coefficients,
            dimensions,
            pmatrix = None,
            qmatrix = None
            ):

        if pmatrix is None:
            pmatrix = self.pmatrix_ho
        if qmatrix is None:
            qmatrix = self.qmatrix_ho

        return self.product_operator_tensor((pmatrix, qmatrix, pmatrix), dimensions)

    def _compute_h1(self, gmatrix_derivs, V_derivs, pQp, QQQ, n, m, dimensions):
        """

        :param gmatrix_derivs:
        :type gmatrix_derivs: np.ndarray
        :param V_derivs:
        :type V_derivs: np.ndarray
        :param pQp:
        :type pQp:
        :param QQQ:
        :type QQQ:
        :param n:
        :type n:
        :param m:
        :type m:
        :param dimensions:
        :type dimensions:
        :return:
        :rtype:
        """
        # we'll vectorize the computation of these components in a bit...


        dimensions = np.asarray(dimensions)
        bits = np.cumprod(dimensions) # should ignore the last one I think...
        n_i = np.dot(bits, np.asarray(n))
        m_i = np.dot(bits, np.asarray(m))

        pQp_shape = pQp.shape
        piQkpj = np.array([
            pQp[i][j][k][n_i][m_i] for k in pQp_shape[2] for j in pQp_shape[1] for i in pQp_shape[0]
        ]) # this could be quite slow but the only alternative is to use tons of memory...

        ke = np.dot(gmatrix_derivs.flatten(), piQkpj)

        QQQ_shape = QQQ.shape
        QiQjQk = np.array([
            QQQ[i][j][k][n_i][m_i] for k in QQQ_shape[2] for j in QQQ_shape[1] for i in QQQ_shape[0]
        ])

        pe = np.dot(V_derivs.flatten(), QiQjQk)

        return 1/2 * ke + 1/6 * pe



