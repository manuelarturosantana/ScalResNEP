% Template class for building different layer potentials
classdef LayerPotential

    properties
        curve;   % The curve to evaluate the layer on
        density; % The computed density. Not initialized on construction.
        R;       % The matrix of weights based on the curve.
        % Properties for MK split kernel only. Defaults to [] for non MK Curves
        A_log         = [];     % Matrix of logarithmic values for each t, tau
        norm_diff     = [];     % Matrix of norm difference of boundary points.
    end

    properties (Constant)
    % The Euler Mascheroni constant from vpa(eulergamma) to many more digits than necessary
       gamma = 0.5772156649015328606065120900824;
    end

    methods 

        function lp = LayerPotential(varargin)
            % Superclass constructor. If the layer potential uses an MK splitting then we 
            % need to precompute the same distance values for each one, so we do that once
            % here. 
            % Inputs:
            %   is_MK: If true, precompute the distances and A_log matrix for MK Splitting
            %   curve: The Curve object for the closed curve.
            %   reg_points: (optional) Points for regularizing the curve.
            % Note they must be passed in in this order

            % Classes which inherit from but don't explicitly call this constructor will
            % call it with no arguments. 
            if nargin > 0 && varargin{1}

                curve    = varargin{2};
                lp.curve = curve;
                
                % Build matrices which will be used for creating matrices for many values 
                % of mu. This from the MK splitting
                A_log = zeros(curve.N, curve.N);
                norm_diff  = zeros(curve.N, curve.N);
        
                % TODO: This is non-optimal runtime implementation.
                for jj = 1:curve.N
                    A_log(jj,:) = 4 * sin((curve.ts(jj) - curve.ts) / 2).^2;
                    norm_diff(jj,:) = vecnorm(curve.xs(:,jj) - curve.xs);
                end

                lp.A_log = log(A_log);
                lp.norm_diff = norm_diff;
                
                lp.R = lp.comp_R(); 
            end
        end

        function lp = set_density(lp, density)
            lp.density = density;
        end

        function vals = eval_point(lp, mu, x_0)
            % Evaluates the layer potential at given points.
            % Inputs:
            %     lp : The single layer potential
            %     mu : The wave number
            %     x_0: 2 x m array of points where we compute the solution values
            %
            % Outputs
            %   vals: m x 1 array of computed function values.
    
            mat = lp.sol_rep_mat(mu, x_0);
    
            vals = mat * lp.density
    
        end

        function R = comp_R(lp)
            % Function to compute the matrix corresponding to the R_{i,j} terms (integration weights in Colton and Kriess)
            % Output
            % R a symmetric matrix with R_{i,j} = R_|i-j|
            
            % Initialize m values from m = 1: n - 1
                ms = 1:lp.curve.n-1;
                ms_i = ms.^-1;
                R_j = [];
            
                % TODO: This implementation can be sped up with DCT.

                % Now we build the first row of R
                % See equation 3.119 in Colton and Kriess
                for jj = 0: lp.curve.N - 1
                    % Here we are subtract t_0 - t_j = -t_j
                    val = sum(ms_i .* cos(ms * jj * pi / lp.curve.n));
                    val = (-2 * pi / lp.curve.n) * val - ((-1)^jj * pi) / lp.curve.n^2;
                    R_j = [R_j, val];
                end
            
                % Build the toeplitz matrix
                R = toeplitz(R_j);
            end
            
            function sol_rep_mat(lp)
               % Build the matrix corresponding to the solution representation at the point x_0
               error("Not Implemented for Double Layer Adjoint and Hypersingular operator") 
            end
            
            function bie_mat(lp)
                error("Not Implemented for Double Layer Adjoint and Hypersingular operator") 
            end
    end

    methods (Abstract)
        lp_mat();  % Build the matrix corresponding to the discretization of the layer potential
    end

end