### Layer Potentials
For closed curves the book by Colton and Kress refers to one of the following

- Colton, David, and Rainer Kress. Integral equation methods in scattering theory. Society for Industrial and Applied Mathematics, 2013
- Colton, David L., Rainer Kress, and Rainer Kress. Inverse acoustic and electromagnetic scattering theory. Vol. 93. Berlin: Springer, 1998.

Sometimes the method used is called the MK split as well.

The code for open curves was adapted from code by Stephane Linter. The relevant paper is

- Bruno, Oscar P., and Stéphane K. Lintner. "Second-kind integral solvers for TE and TM problems of diffraction by open arcs." Radio Science 47.06 (2012): 1-13.

Finally it is worth mentioning again that this implementations are minimal. There is no acceleration
method for the evaluation of the layer potentials implemented, (such as FMM or IFGF), no consideration
of domains with corners, and no considerations for evaluating the integral representation
formula at points very near to the boundary of the curve. All these things exist however and 
could be implemented if desired.
