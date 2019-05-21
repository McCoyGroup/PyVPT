"""Provides a general, convenient FiniteDifferenceFunction class to handle all of our difference FD imps

"""
import numpy as np

##########################################################################################
#
#                                    FiniteDifferenceFunction
#
##########################################################################################

class FiniteDifferenceFunction:
    """The FiniteDifferenceFunction encapsulates a bunch of functionality extracted from Fornberger's
    Calculation of Wieghts in Finite Difference Formulas (https://epubs.siam.org/doi/pdf/10.1137/S0036144596322507)

    """
    REGULAR_GRID = "regular"
    IRREGULAR_GRID = "irregular"
    def __init__(self, coefficients, widths = None, order = None, gridtype = None, regularize_results = False):

        self._coefficients = coefficients # used to apply the weights to function values
        # NOTE: the coefficients should be supplied without 'h' if an evenly spaced grid
        # These will be forced back in later
        self._widths = self._clean_widths(widths) # used to determine how to apply the coefficients
        self._order = self._clean_order(order) # the order of the derivative
        if gridtype is None:
            gridtype = self.REGULAR_GRID
        self._gridtype = gridtype # if 'even' will dispatch to 'call_evenly_spaced' otherwise we'll use more general code
        self._function = None
        self._regr = regularize_results

    @classmethod
    def RegularGridFunction(cls, n, accuracy = 4, end_point_precision = None, dimension = None):
        """Constructs a finite-difference function that computes the nth derivative to a given accuracy order

        :param deriv:
        :type deriv:
        :param accuracy:
        :type accuracy:
        :return:
        :rtype:
        """

        # Gotta make sure we have the right dimensions for our derivatives
        if isinstance(n, int):
            if not isinstance(dimension, int):
                dimension =1
            n = [ n ] * dimension
        else:
            dimension = len(n)
        if isinstance(accuracy, int):
            accuracy = [ accuracy ] * dimension

        coeff_list = [ None ] * dimension
        width_list = [ None ] * dimension
        for i, m in enumerate(n):
            stencil = m + accuracy[i]
            if end_point_precision is None:
                end_point_precision = np.math.ceil(stencil / 2)
            outer_stencil = stencil + end_point_precision
            lefthand_coeffs = cls._even_grid_coeffs(m, 0, outer_stencil)
            centered_coeffs = cls._even_grid_coeffs(m, np.math.ceil(stencil / 2), stencil)
            righthand_coeffs = cls._even_grid_coeffs(m, outer_stencil, outer_stencil)
            stencil = len(centered_coeffs)
            widths = [
                [0, len(lefthand_coeffs)],
                [np.math.ceil(stencil/2), np.math.floor(stencil/2)],
                [len(righthand_coeffs), 0]
            ]
            coeffs = [lefthand_coeffs, centered_coeffs, righthand_coeffs]

            coeff_list[i] = coeffs
            width_list[i] = widths

        return cls(coeff_list, widths = width_list, order = n, gridtype = cls.REGULAR_GRID, regularize_results= True )

    @classmethod
    def IrregularGridFunction(cls, n, grid, stencil = 10  ):
        """Constructs a finite-difference function that computes the nth derivative with a given width

        :param deriv:
        :type deriv:
        :param accuracy:
        :type accuracy:
        :return:
        :rtype:
        """
        if stencil > len(grid):
            stencil = len(grid)

        # we're gonna do just the 1D case for now
        slices = cls._pull_f_vals(grid, 0, stencil)
        slices = np.concatenate( (
            np.tile(slices[0], np.math.floor(stencil/2)),
            slices,
            np.tile(slices[0], np.math.ceil(stencil/2)),
        ) )
        coeffs = [ cls._uneven_spaced_weights(n, z, x) for z, x in zip(grid, slices)]
        widths = [ stencil ]

        coeff_list = [ coeffs ]
        width_list = [ widths ]
        return cls(coeff_list, widths = width_list, order = n, gridtype = cls.IRREGULAR_GRID )

    def _clean_coeffs(self, cfs):
        cf1 = cfs[0]
        if isinstance(cf1, (int, float)):
            cfs = [ [ cfs ] ]
        return cfs
    def _clean_widths(self, ws):

        if isinstance(ws, (int, float)):
            ws = ( (int(ws), int(ws)) )
        else:
            w2s = [ None ] * len(ws)
            for i,w in enumerate(ws):
                if isinstance(w, (int, float)):
                    w2s[i] = (int(w), int(w))
                elif isinstance(w[0], (int, float)):
                    w2s[i] = (int(w[0]), int(w[1]))
                else:
                    w2s[i] = self._clean_widths(w)
            ws = tuple(w2s)
        return ws
    def _clean_order(self, order):
        if not isinstance(order, int):
            try:
                isinst = all((isinstance(o, int) for o in order))
            except:
                isinst = False
            if not isinst:
                raise TypeError("{}: order {} is expected to be an int or list of ints",
                                type(self).__name__, order)
        else:
            order = [ order ]
        return order

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
        return self._even_grid_FDF(self._pull_f_vals, coeffs, widths, orders, regularize=self._regr)

    def _unevenly_spaced_FDF(self):
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
        return self._uneven_grid_FDF(self._pull_f_vals, coeffs, widths, orders)

    @property
    def function(self):
        if self._function is None:
            self._function = self.get_FDF()
        return self._function
    def get_FDF(self):
        if self._gridtype != self.IRREGULAR_GRID:
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
    def _even_grid_FDF(take_lists, coeffs, widths, orders, regularize = False):
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

        def FDF_single(f_vals, h, w, c, o):
            """Calculates a 1D finite difference"""
            cdata = c / (h**o)
            vals = np.asarray(f_vals)
            arrays = take_lists(vals, *w)
            fdf_vals = np.matmul(arrays, cdata)
            return fdf_vals

        def FDF_1D(f_vals, h, widths, coeffs, o):
            """Calculates all the relevant finite differences for the difference widths and coeffs"""
            results = [FDF_single(f_vals, h, w, c, o) for w, c in zip(widths, coeffs) ]
            if len(results) == 1:
                return results[0]
            elif regularize and len(results) == 3:
                missing_points = len(f_vals) - len(results[1])
                left_pad = np.math.ceil( missing_points/2 )
                right_pad = np.math.floor( missing_points/2 )
                results = np.concatenate( ( results[0][:left_pad], results[1], results[2][-right_pad:] ) )

            return results

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
    def _uneven_grid_FDF(take_lists, coeffs, widths, orders, regularize = False):
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

        def FDF_1D(f_vals, w, c):
            """Calculates a 1D finite difference"""
            vals = np.asarray(f_vals)
            slices = take_lists(vals, 0, w)
            np.concatenate( (
                np.tile(slices[0], np.math.floor(w/2)),
                slices,
                np.tile(slices[0], np.math.ceil(w/2)),
            ) )
            fdf_vals = np.sum(slices*c, axis=1)
            return fdf_vals

        if len(widths) == 1:
            def FDF(f_vals):
                "Calculates a 1D finite difference"
                w = widths[0]
                c = coeffs[0]
                return FDF_1D(f_vals, w, c)
        else:
            # we do our finite differencing along multiple different axes as needed
            def FDF(f_vals):
                """Calculates the multi-dimensional FD by FD-ing each dimension in turn"""
                red_vals = f_vals
                for c, w in zip(coeffs, widths):
                    red_vals = FDF_1D(red_vals, c, w)
                return red_vals

        return FDF

    @staticmethod
    def _StirlingS1_mat(n):
        # simple recursive definition of the StirlingS1 function in Mathematica
        # implemented at the C level mostly just for fun
        from .ZachLib import StirlingS1 # loads from C extension

        return StirlingS1(n)

    @staticmethod
    def _Binomial_mat(n):
        # simple recursive Binomial coefficients up to r, computed all at once to vectorize later ops
        # wastes space, justified by assuming a small-ish value for n
        from .ZachLib import Binomial # loads from C extension

        return Binomial(n)

    @staticmethod
    def _old_Bin(n):
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
        sTest = s - int(s)
        if sTest == 0:
            bges = bins[ int(s) ]
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


##########################################################################################
#
#                                    finite_difference
#
##########################################################################################

def finite_difference(grid, values, order, npts = None, mode = FiniteDifferenceFunction.REGULAR_GRID):
    """

    :param grid:
    :type grid: np.ndarray
    :param values:
    :type values: np.ndarray
    :param order: order of the derivative to compute
    :type order: int or list of ints
    :param npts: number of points to use in the stencil
    :type npts: int or list of ints
    :param mode: 'regular' or 'irregular' specifying grid type
    :type mode: str
    :return:
    :rtype:
    """
    gp = np.asarray(grid)
    gv = np.asarray(values)

    if npts is None:
        npts = order + 2

    mode = mode.lower()
    if mode != FiniteDifferenceFunction.IRREGULAR_GRID:
        func = FiniteDifferenceFunction.RegularGridFunction(order, npts - order)
        mesh_spacings = np.array(
            [(lambda g:abs(g[1]-g[0]))(np.take(grid, [0, 1], axis=i)) for i in range(gp.ndim)]
        )
        fdf = func.function
        return fdf(gv, mesh_spacings)
    else:
        func = FiniteDifferenceFunction.IrregularGridFunction(order, grid, npts)
        fdf = func.function
        return fdf(gv)