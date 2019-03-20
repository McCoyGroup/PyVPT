from distutils.core import setup, Extension
import shutil, os, sys

curdir = os.getcwd()
lib_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(lib_dir)
sysargv1 = sys.argv

sys.argv = ['build', 'build_ext', '--inplace']
module = Extension(
    'ZachLib',
    sources = [ 'ZachLib.cpp', "PyExtLib.cpp" ]
)

setup (name = 'ZachLib',
       version = '1.0',
       description = 'Performance critical parts of Zachary package',
       ext_modules = [module],
       language="c++"
       )

os.chdir(curdir)

ext = ""
target = os.path.join(lib_dir, "ZachLib")
src = None

for f in os.listdir(lib_dir):
    if f.endswith(".so"):
        ext = ".so"
        src = os.path.join(lib_dir, f)
        target += ext
        break
    elif f.endswith(".pyd"):
        ext = ".pyd"
        src = os.path.join(lib_dir, f)
        target += ext
        break

if src is not None:
    failed = False
    try:
        os.remove(target)
    except:
        pass
    os.rename(src, target)
    shutil.rmtree(os.path.join(lib_dir, "build"))
else:
    failed = True