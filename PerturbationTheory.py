import numpy as np, scipy.sparse as sp
from McUtils.Zachary.LazyTensors import SparseTensor

class PerturbationTheoryError(Exception):
    pass

class PerturbationTheoryHandler:
    """Holds all the data and methods necessary for doing PerturbationTheory

    """
    def __init__(self,
                 coords,
                 pot_derivs,
                 modes,
                 internals
                 ):
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
            raise PerturbationTheoryError("{}.{}: length of gradient array ({2[0]}) is not {3[0]}}",
                                          type(self).__name__,
                                          "_canonicalize_force_constants",
                                          grad.shape,
                                          (coord_n,)
                                          )
        if fcs.shape != (coord_n, coord_n):
            raise PerturbationTheoryError("{}.{}: dimension of force constant array ({2[0]}x{2[1]}) is not {3[0]}x{3[1]}",
                                          type(self).__name__,
                                          "_canonicalize_force_constants",
                                          fcs.shape,
                                          (coord_n, coord_n)
                                          )
        if thirds.shape != (modes_n, coord_n, coord_n):
            raise PerturbationTheoryError("{}.{}: dimension of third derivative array ({2[0]}x{2[1]}x{2[2]}) is not {3[0]}x{3[1]}x{3[2]}",
                                          type(self).__name__,
                                          "_canonicalize_derivs",
                                          thirds.shape,
                                          (modes_n, coord_n, coord_n)
                                          )
        # this might need to change in the future
        if fourths.shape != (modes_n, coord_n, coord_n):
            raise PerturbationTheoryError("{}.{}: dimension of fourth derivative array ({2[0]}x{2[1]}x{2[2]}) is not {3[0]}x{3[1]}x{3[2]}",
                                          type(self).__name__,
                                          "_canonicalize_derivs",
                                          fourths.shape,
                                          (modes_n, coord_n, coord_n)
                                          )

        return grad, fcs, thirds, fourths

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

    def QQQ(self, coefficients, dimensions, qmatrix = None):
        import functools as ft

        if qmatrix is None:
            qmatrix = self.qmatrix_ho

        dimensions = np.asarray(dimensions)
        M = len(dimensions)
        tensor = [ [ [ None ] * M ] * M ] * M # initialize the 3D thing
        for i in range(M):
            for j in range(M):
                for k in range(M):
                    pieces = [ None ] * M
                    test = {i, j, k}
                    for n, m in enumerate(dimensions): # we gotta do a 4th loop to build our KP terms
                        if n in test:
                            if i == j == k:
                                pieces[n] = qmatrix(m).dot(qmatrix(m)).dot(qmatrix(m))
                            elif i == j or i == k or j == k:
                                pieces[n] = qmatrix(m).dot(qmatrix(m))
                            elif (i == n == j) or (i == n == k) or (j == n == k):
                                pieces[n] = qmatrix(m).dot(qmatrix(m))
                            else:
                                pieces[n] = qmatrix(m)
                        else:
                            pieces[n] = sp.identity(m)

                    tensor[i][j][k] = ft.reduce(sp.kron, pieces)

        return SparseTensor(tensor, (M, M, M))

    # TODO: gotta actually include the coefficients somehow...?
    def pQp(self,
                        coefficients,
                        dimensions,
                        pmatrix = None,
                        qmatrix = None
                        ):

        import functools as ft
        if pmatrix is None:
            pmatrix = self.pmatrix_ho
        if qmatrix is None:
            qmatrix = self.qmatrix_ho

        dimensions = np.asarray(dimensions)
        M = len(dimensions)
        tensor = [ [ [ None ] * M ] * M ] * M # initialize the 3D thing
        for i in range(M):
            for j in range(M):
                for k in range(M):
                    pieces = [ None ] * M
                    test = {i, j, k}
                    for n, m in enumerate(dimensions): # we gotta do a 4th loop to build our KP terms
                        if n in test:
                            if i == j == k:
                                pieces[n] = pmatrix(m).dot(qmatrix(m)).dot(pmatrix(m))
                            elif i == k == n:
                                pieces[n] = pmatrix(m).dot(qmatrix(m))
                            elif i == j == n:
                                pieces[n] = pmatrix(m).dot(pmatrix(m))
                            elif j == k == n:
                                pieces[n] = qmatrix(m).dot(pmatrix(m))
                            elif n == k:
                                pieces[n] = qmatrix(m)
                            else:
                                pieces[n] = pmatrix(m)
                        else:
                            pieces[n] = sp.identity(m)

                    tensor[i][j][k] = ft.reduce(sp.kron, pieces)

        return SparseTensor(tensor, (M, M, M))

    def _compute_h1(self, gmatrix_derivs, V_derivs, pQp, QQQ, n, m, dimensions):
        """

        :param gmatrix_derivs:
        :type gmatrix_derivs: np.ndarray
        :param V_derivs:
        :type V_derivs:
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

        pQp_shape = pQp.shape
        dimensions = np.asarray(dimensions)
        bits = np.cumprod(dimensions) # should ignore the last one I think...
        n_i = np.dot(bits, np.asarray(n))
        m_i = np.dot(bits, np.asarray(m))
        piQkpj = np.array([
            pQp[i][j][k][n_i][m_i] for k in pQp_shape[2] for j in pQp_shape[1] for i in pQp_shape[0]
        ])

        np.dot(gmatrix_derivs.flatten(), )



