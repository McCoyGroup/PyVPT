import numpy as np

class PotentialExpansion:
    """A class for handling expansions of an internal coordinate potential up to 4th order
    Uses Cartesian derivative matrices and the Cartesian <-> Internal normal mode Jacobian"""
    def __init__(self, derivs, transf):
        self.derivs = derivs
        self.transf = transf
        self._d1 = None
        self._d2 = None
        self._d3 = None
        self._d4 = None
        # we'll want a bunch of intermediate expressions too I think...

    @property
    def first_derivs(self):
        if self._d1 is None:
            self._d1 = self._compute_d1()
        return self._d1
    def _compute_d1(self):
        # simple mat-vec
        pass

    @property
    def second_derivs(self):
        if self._d2 is None:
            self._d2 = self._compute_d2()
        return self._d2
    def _compute_d2(self):
        # we'll use d1 and do a few more mat-vecs
        pass

    @property
    def third_derivs(self):
        if self._d3 is None:
            self._d3 = self._compute_d3()
        return self._d3
    def _compute_d3(self):
        # not the worst expression; needs a few tensordots though
        pass

    @property
    def fourth_derivs(self):
        if self._d4 is None:
            self._d4 = self._compute_d4()
        return self._d4
    def _compute_d4(self):
        # now things get rougher since we only want to iterate at most in linear time on the python side
        # it's a bunch of tensor-dots basically
        pass