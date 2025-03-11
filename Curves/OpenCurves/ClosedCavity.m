% Closed Cavity
% A smooth curve which is strongly trapping.
% Taken from Nyström discretizations of boundary integral equations for the solution of 2D elastic scattering problems
% Dominguez and Turc. Journal of Comp and Applied Math 2024
classdef ClosedCavity < Curve

    methods
        function cc = ClosedCavity(n)
            % Constructor
            % The number of points is 2n -1.

            % Helper parameterizations to following the paper
            A   = @(t) sin(t) + sin(2 * t) + (1/2) * sin(3 * t);
            Ap  = @(t) cos(t) + 2 * cos(2 *t) + (3/2) * cos(3 * t);
            App = @(t) -sin(t) - 4 * sin(2 * t) - (9/2) * sin(3 * t);

            As   = @(t) -4 * sin(t) + 7 * sin(2 * t) - 6 * sin(3  *t) + 2 * sin(4 * t);
            Asp  = @(t) -4 * cos(t) + 14 * cos(2 * t) - 18 * cos(3 * t) + 8 * cos(4 * t);
            Aspp = @(t) 4 * sin(t)  - 28 * sin(2 * t) + 54 * sin(3 * t) - 32 * sin(4 * t);

            x_t   = @(t)  (1/4) * vertcat(cos(t) + 2 * cos(2 * t), A(t) /2 - As(t) / 48);
            xp_t  = @(t)  (1/4) * vertcat(-sin(t) - 4 * sin(2 * t), Ap(t) / 2 - Asp(t) / 48);
            xpp_t = @(t) (1/4) * vertcat(-cos(t) - 8 *cos(2 * t), App(t) /2 - Aspp(t) / 48);
            
            cc@Curve(n, x_t, xp_t, xpp_t);
        end 
    end % Methods
end % Classdef