# PyVPT

The Python Vibrational Perturbation Theory package (PyVPT) is intended to make it easy to do vibrational perturbation theory in internal coordinates in python.
 
## Included Components
 
 - **Coordinerds**
 
    The Coordinerds package is essentially a standalone package for handling many types of coordinate transformations in python. It includes a way to represent coordinates in a given basis, an extensible conversion system, and some support for providing fast coordinate transformations
     
 - **PerturbationTheory**
 
    The PerturbationTheory package *will* be the main implementation of perturbtation theory in the package. It'll work with the coordinate work under development to provide a way to do all the expansions in internal coordinates.
 
 - **Tests**
 
    Some integration with the python unittest framework exists. More will come as the package gets more fleshed out. Connects to Travis.
 
 - **TBD**
 
    Some faculty for doing nice finite differencing, Taylor series expansions, etc. will be needed. Where exactly this will go is still TBD.

## Road Map

### Coordinate Basics
 
 - [x] Set up basic project structure
 - [x] Create coordinate system representation
 - [x] Create coordinate transformation system
 - [x] Create coordinate conversion system
 - [x] Implement Cartesian <-> ZMatrix conversions
 - [ ] Generalize conversion code?
 - [ ] Implement Cartesian <-> Spherical conversions
 - [ ] Implement ZMatrix <-> Spherical conversions?

### Molecules
 - [x] Write efficient Gaussian Log File importer
 - [x] Write efficient Gaussian fchk file importer
 - [ ] Write convenience type representing molecular modes
 - [ ] Create efficient importer of all relevant molecule data from Gaussian job data
 - [ ] Create efficient constructor of all molecule data from scan (for use with psi4)
 
### Expansions
 - [x] Implement arbitrary order finite differencing (higher-dimensional extension needs testing still)
 - [ ] Provide utility type representing function over a coordinate set (useful for potentials)
 - [ ] Implement actual expression for a fourth-order potential expansion (make C-compileable...?)

### Perturbation Theory
 - [ ] Get potential expansions with Gaussian potential / force consts
 - [ ] Write code to handle directly from potential scan
 - [ ] Implement proper exprs for perturbation expansions
 
### Wavefunctions / Spectra
 - [ ] Write flexible Wavefunction metaclass that can be subclassed for DVR, DMC, and VPT2
 - [ ] Write Spectrum class that can plot and analyze itself in intuitive ways like I already have in Mathematica

### Future
 - [x] Break these packages into git submodules for others to use... (well on its way)

## Suggestions / Issues / Warnings

This is really not intended for outside use, but if there are features that should be included they can be suggested on the repo [issues page](https://github.com/McCoyGroup/PyVPT/issues)

