% Missile Class
% Class to implement the missile from DIFFARCS as one of my curves.
% Here we use the correct jacobian, (not the corrected one for their numerical method).
classdef Missile < Curve

    methods
        function missile = Missile(n)
            % Constructor
            % This code is Modified from Curve Parametrization in Diffarcs

            C  =  @(t) cos(t);
            S  =  @(t) sin(t);
            C2 = @(t) cos(2*t);
            S2 = @(t) sin(2*t);
            C3 = @(t) cos(3*t);
            S3 = @(t) sin(3*t);
            C4 = @(t) cos(4*t);
            S4 = @(t) sin(4*t);
            C6 = @(t) cos(6*t);
            S6 = @(t) sin(6*t);
            C8 = @(t) cos(8*t);
            S8 = @(t) sin(8*t);

            R = @(t) 0.35+0.1*C(t)+0.12*C2(t)+0.15*C3(t)+0.1*C4(t)+0.1*C6(t)+0.05*C8(t);

            Rp = @(t) -0.1*S(t)-2*0.12*S2(t)-3*0.15*S3(t)-0.1*4*S4(t)-0.1*6*S6(t)-0.05*8*S8(t);

            Rpp= @(t) -0.1*C(t)-4*0.12*C2(t)-9*0.15*C3(t)-0.1*16*C4(t)-0.1*36*C6(t)-0.05*64*C8(t);

            x  = @(t) R(t).*C(t);
            y  = @(t) R(t).*S(t);
            xp = @(t) (-R(t).*S(t)+Rp(t).*C(t));
            yp = @(t) (R(t).*C(t)+Rp(t).*S(t));

            
            x_t   = @(t) vertcat(x(t), y(t));
            xp_t  = @(t) vertcat(xp(t), yp(t));
            xpp_t = @(t) [(-R(t).*C(t)-2*Rp(t).*S(t)+Rpp(t).*C(t)); (-R(t).*S(t)+2*Rp(t).*C(t)+Rpp(t).*C(t))];

            missile@Curve(n, x_t, xp_t, xpp_t);  
        end

        function x_vals = reg_points(missile, n, mu)
            % Function to create the regularization points lambda/2 away from the edge
            % Both inside and outside.
            % Inputs:
            %   n : How many points to use in both the inside and outside regularization curve.
            %   mu: The frequency. The points will be placed a distance of the wavelength/2 away in the normal direction
            % Outputs: 
            %       x_vals: A 2 x 2n array of the new points, with each column a point. The first n are the outside curve.
            
            % Get the wavelength from the frequency
            lambda = (2 * pi) / mu;

            % Generate the t_vals
            t_vals = pi * (0:1:(2* n - 1)) / n;
            curve_vals = missile.x_t(t_vals);

            % Generate the  unit normal vectors
            normal = missile.xp_t(t_vals);
            normal = [-normal(2,:); normal(1,:)];
            normal = normal ./ vecnorm(normal);

            % Place points lambda/2 away in the direction of the normal.
            outside = curve_vals + (lambda / 2) * normal;
            inside  = curve_vals - (lambda / 2) * normal;

            x_vals = [outside, inside];
        end

        function val = aperature_size(missile, ap)
            % Given the ap "size" as requested input by DIFFARCS
            % return the length of the curve removed from the closed curve
            % Inputs:
            %     ap: The "size" of the aperature as input to DIFFARCS
            % Outputs:
            %     val: the length of the curve

            % Compute the limits in t from the DIFFARCS aperature code
            % (See Curve_Parameterization.m)
            theta0=-0.9;
            tmin=theta0+ap/2;
            tmax=theta0+2*pi-ap/2;

            % Integrate the curve length and subtract from the full curve
            % length.
            crack_len = integral(@(t) vecnorm(missile.xp_t(t)),tmin,tmax);
            val = missile.curve_len() -crack_len;
        end % aperature_size

        function ap = calc_ap_size(missile, len)
            % Function which given a length returns the DIFFARCS ap size 
            % needed to remove a portion of size len from the curve
            to_zero = @(aper) missile.aperature_size(aper) - len;
            ap = fzero(to_zero,[0, 2 * pi]);
        end

        function bool = test_distance(missile, x,y)
            % Function to test distance of a point given by diffarcs.
            [t,r]=cart2pol(x,y);
    
    
            C=cos(t);
            C2=cos(2*t);
            C3=cos(3*t);
            C4=cos(4*t);
            C6=cos(6*t);
            C8=cos(8*t);
            
            R=0.35+0.1*C+0.12*C2+0.15*C3+0.1*C4+0.1*C6+0.05*C8;
           
         if(r>abs(R)*(1.01))
             bool=1;
         else
             bool=0;
         end;
        end
        
    end % Methods
end % Classdef




















