function [pol, Z, vals] = aaanep(F,nc,center, radius,  in_circ, nworkers)
    % Function to use AAA to solve an NEP on the interior of a circle.
    % Inputs:
    %   F : A scalr function F whos poles will be approximated.
    %   nc : Number of points along the contour.
    %   center : center of the circle to evaluate f at
    %   radius : radius of the circle to evaluate the f at
    %   in_circ: If true only return the poles in the circle
    %   nworkers: Number of workers for evaluating points along the circle.
    % Outputs:
    %   pols : The pols of the rational approximation
    %   Z    : The points on the contour
    %   vals : Function values along the contour
    
    Z = center + radius * exp(2i * pi * (1:nc)/nc);
    parfor (j = 1:length(Z),nworkers)
    % for j = 1:length(Z)
        vals(j) = F(Z(j)); 
    end
    [~,pol] = aaa(vals,Z);
    if in_circ
        pol = pol(inpolygon(real(pol),imag(pol),real(Z), imag(Z)));
    end
end
