

if __name__=="__main__":
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

    import unittest
    from PyVPT.Tests import *

    unittest.main(verbosity=2)