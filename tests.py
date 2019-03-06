import unittest, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from PyVPT.Tests.CoordinateTests import TestAffineMatrix

if __name__=="__main__":
    unittest.main(verbosity=2)