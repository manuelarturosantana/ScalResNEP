% Class which takes a curve and rotates it counterclockwise
classdef RotatedCurve < Curve 
    methods
        function kt = RotatedCurve(curve,theta,xshift, yshift)
            % constructor:
            %     curve  : The previous curve
            %     theta  : The angle, (in radians) to rotate it by
            %     xshift, yshift: optional arguments to also shift the curve.

            if nargin < 3
                xshift = 0;
                yshift = 0;
            end

            
            % Elements of the rotation matrix. Since these are constants with respect to 
            % t the derivatives are just rotated by theta.
            R = @(t) [cos(t) -sin(t); sin(t) cos(t)];
            % Only the term without the derivative cares about the shift.
            x_t    = @(t) R(theta) * curve.x_t(t) + vertcat(xshift * ones(1,length(t)), yshift * ones(1,length(t)));
            xp_t   = @(t) R(theta) * curve.xp_t(t);
            xpp_t  = @(t) R(theta) * curve.xpp_t(t);
            xppp_t = @(t) R(theta) * curve.xppp_t(t);

            kt@Curve(curve.n, x_t, xp_t, xpp_t,xppp_t);
        end

    end % Method
end % Class

