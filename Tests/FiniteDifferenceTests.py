
from .TestUtils import *
from ..Utils.Zachary import *
import sys

class FiniteDifferenceTests(TestCase):

    # @debugTest
    # def test_finite_difference_2d(self):
    #     grid_1 = 4*np.math.pi / 100 * np.arange(100) + 2*np.math.pi
    #     grid = np.meshgrid(grid_1, grid_1)
    #     self.assertAlmostEqual(1, 0)#round(np.linalg.norm(dtest - sin_d3_vals), 4), 0)

    @validationTest
    def test_finite_difference(self):
        sin_grid = [
            0., 0.2513274123, 0.5026548246, 0.7539822369, 1.005309649, 1.256637061, 1.507964474,
            1.759291886, 2.010619298, 2.261946711, 2.513274123, 2.764601535, 3.015928947, 3.26725636,
            3.518583772, 3.769911184, 4.021238597, 4.272566009, 4.523893421, 4.775220833, 5.026548246,
            5.277875658, 5.52920307, 5.780530483, 6.031857895, 6.283185307
        ]
        sin_grid = np.array(sin_grid)
        sin_vals = [
            0., 0.2486898872, 0.4817536741, 0.6845471059, 0.8443279255, 0.9510565163,
            0.9980267284, 0.9822872507, 0.9048270525, 0.7705132428, 0.5877852523, 0.3681245527, 0.1253332336,
            -0.1253332336, -0.3681245527, -0.5877852523, -0.7705132428, -0.9048270525, -0.9822872507,
            -0.9980267284, -0.9510565163, -0.8443279255, -0.6845471059, -0.4817536741, -0.2486898872, 0.
        ]
        sin_vals = np.array(sin_vals)
        sin_d3_vals = np.array([
            -1., -0.9685831611, -0.87630668, -0.7289686274, -0.535826795, -0.3090169944, -0.06279051953,
            0.1873813146, 0.4257792916, 0.6374239897, 0.8090169944, 0.9297764859, 0.9921147013, 0.9921147013,
            0.9297764859, 0.8090169944, 0.6374239897, 0.4257792916, 0.1873813146, -0.06279051953, -0.3090169944,
            -0.535826795, -0.7289686274, -0.87630668, -0.9685831611, -1.
        ])
        derivs = [finite_difference(sin_grid, sin_vals, n, o) for n, o in ((1, 3), (1, 4), (2, 5), (7, 12), (7, 13)) ]
        # print(derivs, file=sys.stderr)
        dtest = finite_difference(sin_grid, sin_vals, 3, 8)
        # print(dtest, sin_d3_vals, file=sys.stderr)
        self.assertAlmostEqual(round(np.linalg.norm(dtest - sin_d3_vals), 4), 0)

    @validationTest
    def test_stirs(self):
        stir = FiniteDifferenceFunction._StirlingS1_mat(8)
        ans = np.array([
            [1,  0,    0,     0,     0,   0,    0,  0],
            [0,  1,    0,     0,     0,   0,    0,  0],
            [0, -1,    1,     0,     0,   0,    0,  0],
            [0,  2,   -3,     1,     0,   0,    0,  0],
            [0, -6,    11,   -6,     1,   0,    0,  0],
            [0,  24,  -50,    35,   -10,  1,    0,  0],
            [0, -120,  274,  -225,   85, -15,   1,  0],
            [0,  720, -1764,  1624, -735, 175, -21, 1]
        ])
        # print(stir, file=sys.stderr)
        self.assertAlmostEqual(np.round(np.linalg.norm(stir-ans), 6), 0.)
    @validationTest
    def test_bin_gs(self):
        gbin = FiniteDifferenceFunction._GammaBinomial_list(7/2, 8)
        ans = np.array([1., 3.5, 4.375, 2.1875, 0.273438, -0.0273438, 0.00683594, -0.00244141])
        # print(gbin, file=sys.stderr)
        self.assertAlmostEqual(np.round(np.linalg.norm(gbin-ans), 5), 0.)
    # @timeitTest(number=500)
    # def test_bin_speed_old(self):
    #    FiniteDifferenceFunction._old_Bin(55)
    @timeitTest(number=500)
    def test_bin_speed_new(self):
        FiniteDifferenceFunction._Binomial_mat(55)
    @validationTest
    def test_bins(self):
        gbin = FiniteDifferenceFunction._Binomial_mat(8)
        ans = np.array([
            [1, 0, 0,  0,  0,  0,  0, 0],
            [1, 1, 0,  0,  0,  0,  0, 0],
            [1, 2, 1,  0,  0,  0,  0, 0],
            [1, 3, 3,  1,  0,  0,  0, 0],
            [1, 4, 6,  4,  1,  0,  0, 0],
            [1, 5, 10, 10, 5,  1,  0, 0],
            [1, 6, 15, 20, 15, 6,  1, 0],
            [1, 7, 21, 35, 35, 21, 7, 1]
        ])
        # print(gbin, file=sys.stderr)
        self.assertAlmostEqual(np.sum(gbin-ans), 0.)
    @validationTest
    def test_facs(self):
        facs = FiniteDifferenceFunction._Factorial_list(8)
        ans = np.array([1, 1, 2, 6, 24, 120, 720, 5040])
        # print(facs, file=sys.stderr)
        self.assertAlmostEqual(np.linalg.norm(facs-ans), 0.)
    @validationTest
    def test_fd_weights(self):
        coeffs = FiniteDifferenceFunction._even_grid_coeffs(3, 7/2, 7)
        ans = [ -0.0192708, 0.259896, -2.02969, 4.92448, -4.92448, 2.02969, -0.259896, 0.0192708 ]
        coeffs2 = FiniteDifferenceFunction._even_grid_coeffs(3, 0, 7)
        ans2 = [
            -8.058333333333334, 42.53333333333333, -98.225,
            129.66666666666666, -106.04166666666667, 53.6,
            -15.408333333333333, 1.9333333333333333
        ]
        # print(coeffs2, file=sys.stderr)
        r1 = np.round(np.linalg.norm(coeffs-ans), 5)
        r2 = np.round(np.linalg.norm(coeffs2-ans2), 10)
        self.assertAlmostEqual(r1 + r2, 0.)
    @validationTest
    def test_uneven_weights(self):
        import numpy as np
        weights = FiniteDifferenceFunction._uneven_spaced_weights
        uweights =  [
            weights(2, 0, np.array([-2, -1, 0, 1, 2])),
            weights(1, 0, np.array([-3/2, -1/2, 1/2, 3/2])),
            weights(1, 1, np.array([-3, -2, -1, 0, 1]))
            #weights(3, 1/2, np.array([0, 1/3, 1, 2, 7/2, 6]))
            ]
        eweights = FiniteDifferenceFunction._even_grid_coeffs
        targ_weights = [
            eweights(2, 2,   4),
            eweights(1, 3/2, 3),
            eweights(1, 4,   4)
            #[13/10, -5832/1105, 42/5, -177/20, 288/65, -1/1340] # turns out these are wrong...I tested my weights
            # in Mathematica and they actually work great
        ]
        passed = True
        for e, w in zip(targ_weights, uweights):
            norm = np.linalg.norm(e-w[-1])
            if norm > .000001:
                passed = False
                for x in (norm, e, w[-1]):
                    print(x, file=sys.stderr)
        self.assertIs(passed, True)
