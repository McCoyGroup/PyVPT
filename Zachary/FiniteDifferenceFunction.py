"""Provides a general, convenient FiniteDifferenceFunction class to handle all of our difference FD imps

"""
import numpy as np

class FiniteDifferenceFunction:
    """The FiniteDifferenceFunction encapsulates a bunch of functionality extracted from Fornberger's
    Calculation of Wieghts in Finite Difference Formulas (https://epubs.siam.org/doi/pdf/10.1137/S0036144596322507)

    """
    def __init__(self, coefficients, widths, order, gridtype):
        self._coefficients = coefficients # used to apply the weights to function values
        # NOTE: the coefficients should be supplied without 'h' if an evenly spaced grid
        # These will be forced back in later
        self._widths = widths # used to determine how to apply the coefficients
        self._order = order # the order of the derivative
        self._gridtype = gridtype # if 'even' will dispatch to 'call_evenly_spaced' otherwise we'll use more general code

    @classmethod
    def EquispacedGridFunction(cls, n, accuracy = 4, dimension = 1):
        """Constructs a finite-difference function that computes the nth derivative to a given accuracy order

        :param deriv:
        :type deriv:
        :param accuracy:
        :type accuracy:
        :return:
        :rtype:
        """
        stencil = n + accuracy
        lefthand_coeffs = cls._even_grid_coeffs(m, 0, stencil)
        centered_coeffs = cls._even_grid_coeffs(m, stencil / 2, stencil)
        righthand_coeffs = -1*np.flip(lefthand_coeffs)
        widths = [ [0, stencil], [stencil/2, stencil/2], [stencil, 0] ]
        coeffs = [lefthand_coeffs, centered_coeffs, righthand_coeffs]
        coeffs = [coeffs for i in range(dimension)]
        widths = [widths for i in range(dimension)]
        return cls(coeffs, widths, n, "even")

    @property
    def widths(self):
        ws = self._widths
        if isinstance(ws, int):
            w = ((ws, ws))
        else:
            ws = [ (w, w) if isinstance(w, int) else w for w in ws]
        return ws

    def _evenly_spaced_FDF(self):
        """Generates a closure that applies the calculated coefficients to the case of an evenly spaced grid
         in a hopefully efficient manner

        :param h:
        :type h:
        :param function_values:
        :type function_values:
        :return:
        :rtype:
        """

        coeffs = self._coefficients
        widths = self._widths
        orders = self._order
        return self._even_grid_FDF(self._pull_f_vals, coeffs, widths, orders)

    def get_FDF(self):
        if self._gridtype == "regular":
            return self._evenly_spaced_FDF()
        else:
            return self._unevenly_spaced_FDF()

    # def __call__(self, grid_values, function_values):
    #     """This is basically just aspirational. One day this thing will detect the kind of FD it is and cache stuff
    #     and call the built FDF on the data...but today is not that day.
    #
    #     :param grid_values:
    #     :type grid_values:
    #     :param function_values:
    #     :type function_values:
    #     :return:
    #     :rtype:
    #     """
    #     ...

    @staticmethod
    def _pull_f_vals(f_slice, width_left, width_right):
        """Pulls the set of overlapping sublists out of the array slice

        :param f_slice:
        :type f_slice:
        :param width_left:
        :type width_left:
        :param width_right:
        :type width_right:
        :return:
        :rtype:
        """
        # we're gonna use the simple implementation from https://stackoverflow.com/a/43413801/5720002

        f = np.asarray(f_slice)
        # Store the shape and strides info
        shp = f.shape
        s  = f.strides
        L = width_left + width_right

        # Compute length of output array along the first axis
        nd0 = shp[0]-L+1

        # Setup shape and strides for use with np.lib.stride_tricks.as_strided
        # and get (n+1) dim output array
        shp_in = (nd0,L)+shp[1:]
        strd_in = (s[0],) + s
        return np.lib.stride_tricks.as_strided(f, shape=shp_in, strides=strd_in)

    @staticmethod
    def _even_grid_FDF(take_lists, coeffs, widths, orders):
        """For simplicity this will handle mixed derivatives simply by doing each dimension separately

        :param take_lists:
        :type take_lists:
        :param coeffs:
        :type coeffs:
        :param widths:
        :type widths:
        :return:
        :rtype:
        """

        def FDF_1D(f_vals, h, w, c, o):
            """Calculates a 1D finite difference"""
            cdata = c / (h**o)
            vals = np.asarray(f_vals)
            arrays = take_lists(vals, *w)
            fdf_vals = np.matmul(arrays, cdata)
            return fdf_vals

        if len(widths) == 1:
            def FDF(f_vals, h):
                "Calculates a 1D finite difference"
                w = widths[0]
                c = coeffs[0]
                o = orders[0]
                return FDF_1D(f_vals, h, w, c, o)
        else:
            # we do our finite differencing along multiple different axes as needed
            def FDF(f_vals, mesh_spacings):
                """Calculates the multi-dimensional FD by FD-ing each dimension in turn"""
                red_vals = f_vals
                for h, c, w, o in zip(mesh_spacings, coeffs, widths, orders):
                    red_vals = FDF_1D(red_vals, h, c, w, o)
                return red_vals

        return FDF

    @staticmethod
    def _StirlingS1_mat(n):
        # simple recursive definition of the StirlingS1 function in Mathematica
        stirlings = np.eye(n, dtype=np.int64)
        for i in range(1, n):
            for j in range(1, i+1):
                sgn = 1 if ((i-j) % 2) == 0 else -1
                stirlings[i, j] =  sgn*( (i-1)*abs(stirlings[i-1, j]) + abs(stirlings[i-1, j-1]) )
        return stirlings

    @staticmethod
    def _Binomial_mat(n):
        # simple recursive Binomial coefficients up to r, computed all at once to vectorize later ops
        # wastes space, justified by assuming a small-ish value for n
        binoms = np.eye(n, n, dtype=np.int64)
        binoms[:, 0] = 1
        rows = np.arange(2, n, dtype=np.int8)
        mids = np.ceil((rows+1)/2).astype(np.int8)
        for i, k in zip(rows, mids):
            for j in range(1, k):
                binoms[i, j] = binoms[i-1, j-1] + binoms[i-1, j]
                binoms[i, i-j] = binoms[i, j]
        return binoms

    @staticmethod
    def _GammaBinomial_list(s, n):
        # Generalized binomial gamma function
        g = np.math.gamma
        g1 = g(s+1)
        g2 = np.array([g(m+1)*g(s-m+1) for m in range(n)])
        g3 = g1/g2
        return g3

    @staticmethod
    def _Factorial_list(n):
        # I was hoping to do this in some built in way with numpy...but I guess it's not possible?
        # looks like by default things don't vectorize and just call math.factorial
        base = np.arange(n, dtype=np.int64)
        base[0] = 1
        for i in range(1, n):
            base[i] = base[i]*base[i-1]
        return base

    @classmethod
    def _even_grid_coeffs(cls, m, s, n):
        """Finds the series coefficients for x^s*ln(x)^m centered at x=1. Uses the method:

             Table[
               Sum[
                ((-1)^(r - k))*Binomial[r, k]*
                    Binomial[s, r - j] StirlingS1[j, m] (m!/j!),
                {r, k, n},
                {j, 0, r}
                ],
               {k, 0, n}
               ]
             ]

        which is shown by J.M. here: https://chat.stackexchange.com/transcript/message/49528234#49528234
        """

        n = n+1 # in J.M.'s algorithm we go from 0 to n in Mathematica -- which means we have n+1 elements
        stirlings = cls._StirlingS1_mat(n)[:, m]
        bins = cls._Binomial_mat(n)
        if isinstance(s, int):
            bges = bins[ s ]
        else:
            bges = cls._GammaBinomial_list(s, n)
        bges = np.flip(bges)
        facs = cls._Factorial_list(n)
        fcos = facs[m]/facs # factorial coefficient (m!/j!)

        import sys
        coeffs = np.zeros(n)
        for k in range(n):
            # each of these bits here should go from
            # Binomial[s, r - j] * StirlingS1[j, m] *
            bs = bges
            ss = stirlings
            fs = fcos
            bits = np.zeros(n-k)
            for r in range(k+1, n+1):
              bits[r-k-1] = np.dot(bs[-r:], ss[:r]*fs[:r])

            # (-1)^(r - k))*Binomial[r, k]
            cs = (-1)**(np.arange(n-k)) * bins[k:n, k]
            # print(bits, file=sys.stderr)
            coeffs[k] = np.dot(cs, bits)

        return coeffs

    @staticmethod
    def _uneven_spaced_weights(m, z, x):
        """Extracts the grid weights for an unevenly spaced grid based off of the algorithm outlined by
        Fronberger in https://pdfs.semanticscholar.org/8bf5/912bde884f6bd4cfb4991ba3d077cace94c0.pdf

        :param m: highest derivative order
        :type m:
        :param z: center of the derivatives
        :type z:
        :param X: grid of points
        :type X:
        """
        from .ZachLib import UnevenFiniteDifferenceWeights # loads from C extension

        x = np.asarray(x)
        return UnevenFiniteDifferenceWeights(m, z, x).T






