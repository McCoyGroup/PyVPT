from McUtils.Zachary.LazyTensors import Tensor

########################################################################################################################
#
#                                           PotentialExpansion
#
class PotentialExpansionException(Exception):
    pass
class PotentialExpansion:
    """A class for handling expansions of an internal coordinate potential up to 4th order
    Uses Cartesian derivative matrices and the Cartesian <-> Internal normal mode Jacobian"""
    class CoordinateTransforms:
        def __init__(self, transforms):
            self._transf = [Tensor.from_array(t) for t in transforms]
        def __getitem__(self, i):
            if len(self._transf) < i:
                raise PotentialExpansionException("{}: transformations requested up to order {} but only provide to order {}".format(
                    type(self).__name__,
                    i,
                    len(self._transf)
                ))
            return self._transf[i]
    class PotentialDerivatives:
        def __init__(self, derivs):
            self.len = [Tensor.from_array(t) for t in derivs]
        def __getitem__(self, i):
            if len(self.len) < i:
                raise PotentialExpansionException("{}: derivatives requested up to order {} but only provide to order {}".format(
                    type(self).__name__,
                    i,
                    len(self.len)
                ))
            return self.len[i]
    def __init__(self, derivatives, transforms = None):
        self._derivs = self.PotentialDerivatives(derivatives) # derivatives of the potential with respect to the cartesians
        if transforms is None:
            self._transf = None
        else:
             # transformation matrices from cartesians to internals
            self._transf = self.CoordinateTransforms(transforms)
            self._d1 = None
            self._d2 = None
            self._d3 = None
            self._d4 = None
            # we may want a bunch of intermediate expressions too I think...

    def expand(self, *coords):
        coord_tensors = [Tensor.from_array(t) for t in coords]
        for i, t in coord_tensors:
            term

    @property
    def first_derivs(self):
        d = self._d1
        if d is None:
            d = self._compute_derivs(0)
        return d
    @property
    def second_derivs(self):
        d = self._d2
        if d is None:
            d = self._compute_derivs(1)
        return d
    @property
    def third_derivs(self):
        d = self._d3
        if d is None:
            d = self._compute_derivs(2)
        return d
    @property
    def fourth_derivs(self):
        d = self._d4
        if d is None:
            d = self._compute_derivs(3)
        return d

    def _compute_derivs(self, i):
        # simple mat-vec
        if self._transf is None:
            term = self._derivs[i]
        elif i == 0:
            term = self._transf[0]*self._derivs[0]
        elif i == 1:
            term = (
                    self._transf[1]*self._derivs[0] +
                    (self._transf[0] ** 2) * self._derivs[1]
            )
        elif i == 2:
            term = (
                    self._transf[2] * self._derivs[0] +
                    2 * ( self._transf[1] * self._transf[0] * self._derivs[1] )+
                    self._transf[0] * self._transf[1] * self._derivs[1] +
                    (self._transf[0] ** 3) * self._derivs[2]
            )
        elif i == 3:
            term = (
                    self._transf[3]*self._derivs[0]+
                    3*self._transf[2]*self._transf[0]*self._derivs[1]+
                    3*(self._transf[1]**2)*self._derivs[1]+
                    self._transf[0]*self._transf[2]*self._derivs[1]+
                    2*self._transf[0]*self._transf[1]*self._transf[0]*self._derivs[2]+
                    3*self._transf[1]*(self._transf[0]**2)*self._derivs[2]+
                    (self._transf[0]**2)*self._transf[1]*self._derivs[2]+
                    (self._transf[0]**4)*self._derivs[3]
            )

        else:
            raise PotentialExpansionException("{}: expansion only provided up to order {}".format(type(self).__name__, 4))

        return term