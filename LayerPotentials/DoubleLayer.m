
classdef DoubleLayer < LayerPotential
    % Class for creating the double layer potential and its adjoint.
    % Notation tries to follow section 3.6 of Colton and Kriess.

    properties
        normd_term     % Term corresponding to normal derivative times the gradient.
        L2_diag        % The diagonal of L2
    end

    methods
        function dl = DoubleLayer(curve)
            % Constructor
            %   Curve     : A curve object already initialized
            %   reg_points: (optional) a 2 x n vector of points which will be used for the regularization in conical newtons.
            
            % Initialize as a MK split curve
            dl = dl@LayerPotential(true, curve);   

            % Normal derivative term contained in L1, as a matrix. (x(t) - t(tau))
            % -1 normd_term is the term in L
            t1 = (curve.xs(1,:).' - curve.xs(1,:)) .* (curve.xps(2,:));
            t2 = (curve.xs(2,:).' - curve.xs(2,:)) .* (curve.xps(1,:));
            dl.normd_term = t1 - t2;
            
            % Diagonal term for L2 as a vector
            L2_diag = curve.xps(1,:) .* curve.xpps(2,:) - curve.xps(2,:) .* curve.xpps(1,:);
            dl.L2_diag = (1/(2 * pi)) * L2_diag ./ curve.n_xps.^2;
        end
        
        function DL_mat = lp_mat(dl, mu)
            % Inputs:
            % mu    : mu^2 = k^2 = lambda.
            %
            % Note that i represents varying in values t
            % and j represents varying t_j, the second aurgument of the kernl.
            
            % Output:
            % DL_mat 2n x 2n matrix corresponding the double layer potential
               
            % Use the functions below to compute the matrix A
            L1 = dl.L_1(mu);
            L2 = dl.L_2(mu, L1);
            
            % We add the negative sign since they took it outside the integral in Colton
            % and Kriess. Also divide by 2 to not use their convention.
            DL_mat = (-1/2) * (dl.R .* L1 + (pi/ dl.curve.n) * L2);
        end

        function DLA_mat = adj_lp_mat(dl,mu)
            % Computes the matrix for the adjoint of the double layer potential, which takes
            % the normal derivative in x.
            % Inputs:
            % mu    : mu^2 = k^2 = lambda.
            % Output:
            % DLA_mat 2n x 2n matrix corresponding the double layer potential adjoint.

            % The normal derivative is with respect to the fixed variable
            % down the columns of the lp_mat so we begin by taking the transpose.
            DLA_mat = dl.lp_mat(mu).';

            % Now for the gotcha! In the first equation of Colton and Kress page 92 (4th edition)
            % they don't have the Jacobian because it cancels out with the normalization factor
            % on the normal derivative. However because the normalization factor is in x,
            % and the jacobian is in y we need to multiply by both of these
            % Multiply jacobian on each element of the rows
            DLA_mat = DLA_mat .* dl.curve.n_xps;
            % Divide by normal derivative multiplication factor on each column.
            DLA_mat = DLA_mat ./ (dl.curve.n_xps.');
        end

        function DL_mat = bie_mat(dl, mu)
            %Per Colton and Kriess (the old one) Theorem
            % 3.15 By adding the identity we can make the double layer a solution of the 
            % exterior Neuman Problem
            DL_mat = lp_mat(dl, mu);
            s     = size(DL_mat);
            index = 1:s(1)+1:s(1)*s(2);
            DL_mat(index) =  DL_mat(index) + 1/2;
            %(1/2) * eye(size(DL_mat)) + DL_mat;
        end

        function A = L_1(dl, mu)
            % Function to compute L_1
            % Input: 
            %    mu^2 = k^2 = lambda.
            %
            % Output:
            % A: The matrix with A_{i,j} = L_1(i,j)
            
            % Build the matrix corresponding to L1
            A = (mu /  (2 * pi)) * dl.normd_term .* ...
            besselj(1, mu * dl.norm_diff) ./ dl.norm_diff;
            
            % While smooth the diagonal numerically blows up so we set it to its true value
            % based on bessel series representation
            s     = size(A);
            index = 1:s(1)+1:s(1)*s(2);

            A(index) = zeros(size(index));

        end
            
        function L2 = L_2(dl, mu, L1)
            % Function to return L_2
            % Inputs:
            % mu : The wave number
            % M1: The matrix M1
            %
            % Output:
            % L2: The matrix with L_2 from the splitting
            
            % Build the matrix for A
            % The normd_term is negated since the signs are switched
            L = (1i * mu / 2) * (-dl.normd_term) .* ...
                besselh(1, mu * dl.norm_diff) ./ dl.norm_diff;

            L2 = L - L1 .* dl.A_log;

            % Grab the indices of the diagonal of M2
            s     = size(L2);
            index = 1:s(1)+1:s(1)*s(2); 

            L2(index) = dl.L2_diag;
        end

        function mat = sol_rep_mat(dl, mu, x_0)
            % Create the matrix corresponding to the solution representation at the points x_0
            % Inputs: 
            %   mu  : mu^2 = k^2 Wavenumber
            %   x_0 : A 2 x n vector consisting of the points to be evaluated at
            % Outputs:
            %     mat: A n x (length(density)) matrix such that mat * density = SL.x_0

            % Compute i/2 partial H_0^(1)(mu|x(t) - x(\tau)|)/ partial eta |x'(\tau)| as a 
            % matrix with rows corresponding to each x point to be valuated at, and 
            % columes correspoding the tau_i values.
            m = size(x_0, 2);
            mat = zeros(m, dl.curve.N);
            
            for ii = 1:m
                % Compute the bessel function term
                temp = (1i/ 4) * mu * besselh(1, mu * vecnorm(x_0(:,ii) - dl.curve.xs));
                temp = temp ./ vecnorm(x_0(:,ii) - dl.curve.xs);

                % Compute the chain rule part of the gradient multiplied with normal
                temp2 = (x_0(1,ii) - dl.curve.xs(1,:)) .* dl.curve.xps(2,:);

                % Minus sign from computing the normal
                temp2 = temp2 - (x_0(2,ii) - dl.curve.xs(2,:)) .* dl.curve.xps(1,:);

                % Note the division by the norm of the tangent to make it the normal 
                % derivative cancels out with the jacobian.
                % % Divide by the length of the normal to make it the unit normal
                % temp2 = temp2 ./ vecnorm(dl.curve.xps);.* dl.curve.n_xps;

                % Multiply by the jacobian
                mat(ii,:) = temp .* temp2;
            end

            % Constants for trapezoidal rule integration.
            mat = (pi / dl.curve.n) * mat;
        end
    end % Methods
end % Classdef
    
