function [polsr, errs, nfevals] = sec_ref(f,pols, num_its, tol, num_workers)
    % Wrapper function which performs the secant method to refine the
    % approximation of the poles
    % Inputs:
    %      f       : Function which has poles approximately at pols
    %      pols    : The approximate poles of the function.
    %      num_its : Optional: Number of iterations in the secant method.
    %                Default 4
    %      tol     : Optional: Secant Method tolerance. Default 1e-13
    %      num_workers : Optional: Number of workers for the for loop to
    %                    which performs the secant method iteration.
    %                    Default 0
    % Warning: By default the second guess for the secant method is taken to be pol
    % + 1e-5.
    %
    % Outputs:
    %   polsr       : The refined poles of the function
    %   max_err     : The errors returned by the secant method.
    %   nfevals     : The total number of iterations + the number of poles

    if nargin == 2
        num_its = 4;
        tol = 1e-13;
        num_workers = 0;
    end

    if nargin == 3
        tol = 1e-12;
        num_workers = 0;
    end

    if nargin == 4
        num_workers = 0;
    end

    polsr = zeros(size(pols));
    errs  = zeros(size(pols));

    % parfor (pp = 1:length(pols), num_workers)
     nfevals = length(pols);
     for pp = 1:length(pols)
        [x,~, iters, x_diff] = secant_method(@(z)1/f(z), pols(pp), pols(pp) + 1e-5, num_its, tol, 1);
        polsr(pp) = x;
        errs(pp)  = x_diff;
        nfevals = nfevals + iters;
    end

end