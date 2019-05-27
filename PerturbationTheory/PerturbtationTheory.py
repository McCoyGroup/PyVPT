import numpy as np

class PerturbationTheoryError(Exception):
    pass

class PerturbationExpansion:
    """Holds all the data and methods necessary for doing PerturbationTheory

    """
    def __init__(self,
                 masses,
                 modes, # a Modes object as implemented in Coordinerds
                 force_constants, # Should I find some way to embded the requisite CartesianCoordinateSystem...?
                 derivs, # the third and fourth derivatives as spit out by Gaussian
                 internals
                 ):
        self.masses = masses
        self.n = len(masses)
        self.mode_n = len(masses) - 6 # number of modes
        self.fcs = self._canonicalize_force_constants(force_constants)
        self.thirds, self.fourths = self._canonicalize_derivs(derivs)
        self.internals = internals

    def _canonicalize_force_constants(self, fcs):
        if not isinstance(fcs, np.ndarray):
            fcs = fcs.array # not type checking so that it can be extensible
        n = self.n
        if fcs.shape != (3*n, 3*n):
            raise PerturbationTheoryError("{}.{}: dimension of force constant array ({2[0]}x{2[1]}) is not {3[0]}x{3[1]}",
                                          type(self).__name__,
                                          "_canonicalize_force_constants",
                                          fcs.shape,
                                          (3*n, 3*n)
                                          )
        return fcs

    def _canonicalize_derivs(self, derivs):
        try:
            thirds, fourths = derivs
        except:
            thirds = derivs.third_deriv_array
            fourths = derivs.fourth_deriv_array
        n = self.n
        if thirds.shape != (3*n-6, 3*n, 3*n):
            raise PerturbationTheoryError("{}.{}: dimension of third derivative array ({2[0]}x{2[1]}x{2[2]}) is not {3[0]}x{3[1]}x{3[2]}",
                                          type(self).__name__,
                                          "_canonicalize_derivs",
                                          thirds.shape,
                                          (3*n-6, 3*n, 3*n)
                                          )
        if fourths.shape != (3*n-6, 3*n, 3*n):
            raise PerturbationTheoryError("{}.{}: dimension of third derivative array ({2[0]}x{2[1]}x{2[2]}) is not {3[0]}x{3[1]}x{3[2]}",
                                          type(self).__name__,
                                          "_canonicalize_derivs",
                                          thirds.shape,
                                          (3*n-6, 3*n, 3*n)
                                          )
        return (thirds, fourths)

    def internalize_force_constants(self):
        """Converts the force_constant data to an internal coordinate representation

        :return:
        :rtype:
        """
        pass

    def internalize_derivs(self):
        """Converts the derivs data to an internal coordinate representation

        :return:
        :rtype:
        """
        pass

    def build_perturbation_wavefunction(self):
        """Constructs the perturbation wavefunction from the second, third, and fourth deriv data

        :return:
        :rtype:
        """
        pass