# ScalResNEP
![alt text](Rockets.png)
Repository containg code to compute Laplace eigenvalues and solve nonlinear eigenvalue problems
using the rational approximation of scalarized resolvents following the paper.
- Bruno, Oscar P., Manuel A. Santana, and Lloyd N. Trefethen. "Evaluation of resonances: adaptivity and AAA rational approximation of randomly scalarized boundary integral resolvents." arXiv preprint arXiv:2405.19582 (2024).

The [examples](/examples) folder contains code for reproducing some of the numerical examples from the paper, and hopefully is clear enough to extend to all of the examples in the paper. Feel free to 
raise an issue if any questions arise.

### Dependencies 
- [Chebfun](https://www.chebfun.org/download/) Only the stand alone code for the aaa algorithm is required.
- [NLEVP](https://github.com/ftisseur/nlevp) Required to run the butterfly and CD player example.


