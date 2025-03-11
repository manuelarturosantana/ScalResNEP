% Kite Curve
classdef Kite < Curve 

    properties
        % Kite is defined as (x(t), y(t)) with
        % x(t) = Acos(t) + Bcos(2t) + C
        % y(t) = Dsin(t)
        % params is a vector [A,B,C,D]
        params
    end

    methods
        function kt = Kite(n, ps, use_CK)
            % Constructor
            % ps is a vector with [A,B,C,D]
            % use_KC: changes ps to be the coefficients in Colton and Kriess, which is 
            %         the default kite parameters used in the paper.
            if nargin > 2 && use_CK
                ps = [1, 0.65, -0.65, 1.5];
            end

            x_t    = @(t) vertcat(ps(1) * cos(t) + ps(2) * cos(2 * t) + ps(3), ps(4) * sin(t));
            xp_t   = @(t) vertcat(-ps(1) * sin(t) - 2 * ps(2) * sin(2 * t), ps(4) * cos(t));
            xpp_t  = @(t) vertcat(-ps(1) * cos(t) - 4 * ps(2) * cos(2 * t), -ps(4) * sin(t));
            xppp_t = @(t) vertcat(ps(1) * sin(t) + 8 * ps(2) *  sin(2 * t), -ps(4) * cos(t));

            kt@Curve(n, x_t, xp_t, xpp_t,xppp_t);

            kt.params = ps;
        end


    
    end % Method
end % Class