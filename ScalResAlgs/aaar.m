function [pol, deg, Z, vals] = aaar(f, p,r,n,deg)
    % Function which takes in a pole and does a local AAA approximation in
    % order to increase the approximation of the poles.
    % Experience says this works best with small input deg, a good number
    % if 4.
    % Inputs:
    %   f   : the function
    %   p   : The initial estimate for the pole
    %   r   : The radius of the circle around the pole 1e-5 is a good rule
    %         of thumb for this value.
    %   n   : The number of points to use on the circle.
    %   deg : The maximum degree of the rational approximation. 4 is a good
    %         rule of thumb
    
    if isempty(deg)
        deg = 99; % aaa default
    end
    % We evaluate at the pole and the small circle around it.
    Z = p + r * exp(2i*pi*(1:n)/n);
    for k = 1:length(Z), vals(k) = f(Z(k)); end
    [~, pol, ~, ~, ~, ~, ~, errvec] = aaa(vals,Z, 'degree',deg, 'lawson',0);
    deg = length(errvec) - 1;
end