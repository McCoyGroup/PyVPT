import os, sys

lib_dir = os.path.dirname(os.path.abspath(__file__))
try:
    sys.path.insert(0, lib_dir)
    from ZachLib import *
except ImportError:
    from .setup import failed
    if failed:
        raise
    from ZachLib import *
finally:
    sys.path.pop(0)