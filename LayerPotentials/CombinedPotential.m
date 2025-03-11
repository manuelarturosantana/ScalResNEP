
classdef CombinedPotential
    % Class for the combined field layer potential. See 3.28, 3.29 in Colton and Kress

        properties
            eta; % The real coupling parameter. Pass in eta ~= 0
            sl;  % The single layer potential.
            dl;  % The double layer potential.
            curve; % The curve
        end

        methods
            function cp = CombinedPotential(curve, eta)
                % Constructor
                %   Curve:  curve object already initialized
                %   eta  : the real coupling parameter. Pass in eta ~= 0.
                
                cp.eta = eta;
                cp.sl = SingleLayer(curve);
                cp.dl = DoubleLayer(curve);
                cp.curve = cp.dl.curve;
            end

            function cp_mat = lp_mat(cp, mu)
                % Inputs:
                % mu    : mu^2 = k^2 wavenumber
                %
                % Outputs:
                %   3.28 in the new Colton and Kriess, however everything is divided by 2
                cp_mat = cp.dl.lp_mat(mu) - 1i * cp.eta * cp.sl.lp_mat(mu);

            end
          
            function com_mat = bie_mat(cp, mu)
                % Inputs:
                % mu    : mu^2 = k^2 = lambda.
                %
                % Outputs:
                %   3.29 in the new colton and kriess,
                com_mat = cp.dl.bie_mat(mu) - (1i) * cp.eta * cp.sl.bie_mat(mu);
            end
    
            function mat = sol_rep_mat(cp, mu, x_0)
                mat = cp.dl.sol_rep_mat(mu, x_0) - 1i * cp.eta * cp.sl.sol_rep_mat(mu, x_0);
            end

        end % Methods    
end % Classdef
    