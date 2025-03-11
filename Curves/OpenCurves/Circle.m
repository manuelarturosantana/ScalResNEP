% Circle Class
% This provides a good test case as the interior and exterior eigenvalues are known 
% from separation of variables.
classdef Circle < Curve

    properties
        radius; % The radius of the circle Default 1
        center; % The center of the circle Default [0;0].
    end

    methods
        function circ = Circle(n, rad, center)
            % Constructor
            % rad is the radius of the circle default 1
            % center is the center of the circle default [0;0]
            % The number of points is 2n -1.

            if nargin < 2
                rad = 1;
            end

            if nargin < 3
                center = [0;0];
            else
                % Make sure it is a column vector
                center = center(:);
            end

            x_t    = @(t) center + rad * vertcat(cos(t),sin(t));
            xp_t   = @(t) rad * vertcat(-sin(t), cos(t));
            xpp_t  = @(t) rad * vertcat(-cos(t), -sin(t));
            xppp_t = @(t) rad * vertcat(sin(t), -cos(t));
            
            circ@Curve(n, x_t, xp_t, xpp_t,xppp_t);

            circ.radius = rad;
            circ.center = center;
        end

        function bool = test_distance(circ, x, y)
            % Test if we are outside the circle.
            bool = ((x-circ.center(1))^2 + (y-circ.center(2))^2) > circ.radius^2;
        end
        
    end % Methods
end % Classdef