# ScalResNEP
![alt text](Rockets.png)
Repository containg code to compute Laplace eigenvalues and solve nonlinear eigenvalue problems
using the rational approximation of scalarized resolvents following the paper.
- Bruno, Oscar P., Manuel A. Santana, and Lloyd N. Trefethen. "Evaluation of resonances: adaptivity and AAA rational approximation of randomly scalarized boundary integral resolvents." arXiv preprint arXiv:2405.19582 (2024).

The [examples](/examples) folder contains code for reproducing some of the numerical examples from the paper, and hopefully is clear enough to extend to all of the examples in the paper. Feel free to 
raise an issue if any questions arise.

### Overview of Directories
- [Curves](Curves) Contains the code for curve parameterizations for open and closed curves.
- [Dependencies](Dependencies) Contains stand alone versions of the AAA algorithm and the butterfly and CD player NEP from nlevp. The licenses for Chebfun and and the nlevp are also included.
- [Examples](Examples) Contains the code to reproduce the examples in the paper.
- [LayerPotentials](LayerPotentials) Contains the code for the discretization of the boundary integral equations
- [ScalResAlgs](ScalResAlgs) Contains implementations of the algorithms introduced in the paper and some utility functions.


### Dependencies 
If for some reason the dependencies do not work they can be downloaded here.
- [Chebfun](https://www.chebfun.org/download/) Only the stand alone code for the aaa algorithm is required.
- [NLEVP](https://github.com/ftisseur/nlevp) Required to run the butterfly and CD player example.


