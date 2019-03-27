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

### Expansions
 
 - [x] Implement arbitrary order finite differencing (still needs higher-dimensional extension)
 - [ ] Provide utility type representing function over a coordinate set
 - [ ] Implement actual expression for a fourth-order potential expansion (make C-compileable...?)

### Perturbation Theory
 
 - [ ] Write efficient Gaussian Log File importer (in process)
 - [ ] Get potential expansions with Gaussian potential / force consts
 - [ ] Implement proper exprs for perturbation expansions


## Suggestions / Issues / Warnings

This is really not intended for outside use, but if there are features that should be included they can be suggested on the repo [issues page](https://github.com/McCoyGroup/PyVPT/issues)

