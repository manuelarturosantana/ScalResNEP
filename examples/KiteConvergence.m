%% Kite Convergence
% Script which gives an convergence analysis for discretization of the
% kite. Additionally this script shows how to solve the Helmholtz equation
% for an incident plane wave.


% First we create a solution of the helmholtz equation which is valid everywhere
% in the exterior of the curve by putting a point source inside the curve (here at 
% we place the point source at the point (0,0).
% Slternatively this could be changed to a plane wave of the form
% exp(1i* k * p * curve.xs) where p is the direction of the plane wave.
% In this case rather than use the sol_func as a true solution just a 
% convergence analysis should be performed.
k = 6;
sol_func = @(x) besselh(0,k * x);



% Now we compare our solution to the actual solution outside the curve.
% test_points = [-1.5, -1, -0.5, 0 , 0.5, 1, 1.5];
test_points = [-1.5, -1.2, -0.54, 0 , 0.52, 1.3, 1.3];
test_points = [test_points; 2 * ones(1, length(test_points))];
% Change this to a converged solution if using a plane wave.
true_sol = sol_func(vecnorm(test_points)).';


ns = [10:5:100];

err = [];
for N = ns

    kite = Kite(N,[],true);
    
    % Now we create the RHS for the integral equation, letting x_0 = 0
    f_vals = sol_func(vecnorm(kite.xs)).';
    
    % Now we create the discretized layer potential operator
    dl = DoubleLayer(kite);
    
    % Here we make the matrix corresponding to the boundary integral
    % equation
    mat = dl.bie_mat(k);
    phi = mat \ (f_vals);
    
    % Create the discretized representation formula for the integral
    % equation
    sol_mat = dl.sol_rep_mat(k, test_points);
    % Multiply by the density to compute the numerical solution.
    num_sol = sol_mat * phi;


    err = [err,max(abs(num_sol - true_sol))];
end

figure(1)
semilogy(2*ns,err,'*-')
xlabel("Number of Discretization Points")
ylabel("err")
