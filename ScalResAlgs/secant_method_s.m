function [x,f_val, iters, x_diff, fevals] = secant_method_s(f, x_old, x_curr, max_iters, tol, xmin,xmax,ymin,ymax)
    % The secant method for when we don't know the derivative
    % the _s signifies that this implementation exits if the iterations
    % leaves the square in the complex plane given by xmin, xmax, ymin, ymax
    % where x and y refer to the real and complex part of the value respectively
    % Inputs:
    %   f         : The function to find the zero of.
    %   x_old     : The first initial starting point
    %   x_curr    : The second initial starting point
    %   max_iters : Maximum number of iterations before stopping 
    %   tol       : Stopping tolerance when change in root is small < tol
    % 
    %  Outputs:
    %   x      : The zero
    %   f_val  : The function value at the zero.
    %   iters  : The number of itertions ran
    %   x_diff : The change in root approximation for the last two iterations.

    f_old = f(x_old);
    f_curr = f(x_curr);

    fevals = 2;
    
    iters = 0;

    % We use the change in root approximation rather than function
    % evaluation because the change in function evaluation could converge
    % slowly
    while (iters < max_iters) && abs(x_curr - x_old) > tol
        quot = (x_curr - x_old) / (f_curr - f_old);

        x_new = x_curr - f_curr * quot;

        % To pick our new point in the secant method we
        % use which ever one has a smaller function evaluation
        if abs(f_old) > abs(f_curr)
            f_old  = f_curr;
            x_old  = x_curr;
        end

        f_curr = f(x_new);
        fevals = fevals + 1;
        f_val  = f_curr;
        x_curr = x_new;

        % If statement to stop early if the method is diverging.
        if isempty(inrectangle(x_curr,xmin,xmax,ymin,ymax))
            x_diff = NaN;
            break
        end

        x_diff = abs(x_curr - x_old);

        iters = iters +1;
    end
    
    f_val = f_curr;
    x = x_new;
end

