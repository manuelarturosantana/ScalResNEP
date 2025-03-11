% Abstract Class to hold curves.
% To create a new curve, make a new class which inherits from the Curve class. In order to 
% work with the layer potentials the curve must be positively oriented, and 
% be parameterized with t \in [0, 2pi]
classdef Curve

    properties 
        % The reason for 2n- 1 points is so we can follow 
        % Colton and Kriess Fourth edition page 93.
        n;    % Number of points in the curve is 2n
        N;    % N = 2n
        %  Note this won't include the starting point again.
        ts;    % Row array with t_j is pi * j / n for j = 0,.., 2n- 1S
        xs;    % 2 x N array of the boundar coordinates of the curve.
        xps;   % 2 x N array of derivative of the curve at each point.
        xpps;  % 2 x N array of the second derivative of the curve at each point.
        xppps; % 2 x N array of the third derivatives of the curve at each point.
        n_xps; % The norm of the derivative of the curve.
        x_t;   % Function handle taking in t_vals, returns 2 x N array of points for the curve. Assumed from 0:2 * pi
        xp_t;  % Function handle taking in t_vals returns a 2 x N array of points of the derivative.
        xpp_t; % Function handle taking in t_vals returns  2 x N array of points of the second derivative 
        xppp_t;% Function handle taking in t_vals retunrs 2 x N array of points of the third derivative.

    end 

    methods
        %
        % Superclass constructor. Each subclass should hold x_t, xp_t, xpp_t, xppp_t
        %
        function curve = Curve(n, x_t, xp_t, xpp_t, xppp_t)

            curve.n = n;
            curve.N = 2 * n;
            curve.ts = pi * (0:1:(2* curve.n - 1)) / curve.n;

            curve.x_t    = x_t; 
            curve.xp_t   = xp_t;
            
            
            curve.xs    = curve.x_t(curve.ts);
            curve.xps   = curve.xp_t(curve.ts);
            curve.n_xps = vecnorm(curve.xps);
            if nargin > 3
                curve.xpp_t  = xpp_t;
                curve.xpps  = curve.xpp_t(curve.ts);
            end

            if nargin > 4
                curve.xppp_t = xppp_t;
                curve.xppps = curve.xppp_t(curve.ts);
            end
        end

        
        function draw(curve,varargin)
            % Simple Method to plot the curve to make sure it is correct
            % varargin hold parameters for plot function.
            
            % Re link the curve with the starting point.
            xxs = horzcat(curve.xs(1,:), curve.xs(1,1));
            yys = horzcat(curve.xs(2,:), curve.xs(2,1));
            plot(xxs, yys,varargin{:})
        end

        function len = curve_len(curve)
            % Compute the curve length by trapezoidal rule
            len = integral(@(t) int_fun(curve,t),0,2*pi);
        end

        function bool = test_distance(curve,x,y)
            % Return true if the value is outside the curve
            bool = ~inpolygon(x,y,curve.xs(1,:),curve.xs(2,:));
        end
    end
end


% integrand for curve length
function vals = int_fun(curve, t)
    points = curve.xp_t(t);
    vals = sqrt(points(1,:).^2 + points(2,:).^2);
end