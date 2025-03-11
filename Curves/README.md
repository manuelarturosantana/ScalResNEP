### Curves
Because the curve for the open curves was not written by me the API for the open and closed
curves is different. An example of how to use each is given in the examples. However, it 
should be more straight forward to extend the closed curves following the object oriented
method.

Two closed curve gotchas!

For method bie_mat gives the matrix corresponding to the exterior Dirichelt problem.
In order to solve interior problems use the method lp_mat along with appropriate
addition of the identity. See Colton and Kress chapter 3.

You tell it to discretize with n points, but it will actually use N = 2 * n points.
This is a laziness of programming to match the notation in Colton and Kress.