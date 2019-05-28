"""These are the set of tests to the PyVPT package
Every type of tests should be its own module and should be tagged as either a debugTest, validationTest, or timingTest


"""

# we'll put the parent dir on path
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import PyVPT.Peeves.run_tests as rt # runs the tests
