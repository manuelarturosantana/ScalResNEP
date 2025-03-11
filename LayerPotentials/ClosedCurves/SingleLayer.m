
classdef SingleLayer < LayerPotential
% Class for first kind single layer equatoin

    methods
        function sl = SingleLayer(curve)
            % Constructor
            %   Curve      : a curve object already initialized
            %   reg_points : (optional) 2 x n vector of points which will be used for the regularization for conical newtons. 
            
            % Initialize as a MK split curve
            sl = sl@LayerPotential(true, curve);         
        end
      
        function SL_mat = lp_mat(sl, mu)
            % Inputs:
            % mu    : mu^2 = k^2 the wave number
            %
            % Note that i represents varying in values t
            % and j represents varying t_j, the second aurgument of the kernel.
            % Works for interior and exterior Dirichelt Problems
            
            % Output:
            % SL_mat 2n x 2n matrix corresponding to the LHS of 3.44 in Colton and Kriess
            % The older book. (Here 2n is the number of points in the discretization)
            
            
            % Use the functions below to compute the matrix A
            M1 = sl.M_1(mu);
            M2 = sl.M_2(mu, M1);
            
            % Multiply by 0.5 so we don't use the KC convention, which is annoying for then
            % computing the solution representation.
            SL_mat = 0.5 * (sl.R .* M1 + (pi/ sl.curve.n) * M2);
        end

        function SL_mat = bie_mat(sl, mu)
            % The lp_mat is the same as the bie_mat
            SL_mat = sl.lp_mat(mu);
        end

        % Function for M_1
        function A = M_1(sl, mu)
            % Function to compute M_1
            % Input: 
            %   mu
            % Same as build_single_layer function
            %
            % Output:
            % A: The matrix with A_{i,j} = M_1(i,j)
            
            % Build the matrix corresponding to J_0(k|x(t) - x(tau)|)
            A = besselj(0, mu * sl.norm_diff) .* sl.curve.n_xps;
           
            A = (- 1 / (2 * pi)) * A ;
            
        end
            
        % Function for M_2
        function M2 = M_2(sl, mu, M1)
            % Function to return M_2
            % Inputs:
            % mu : The wave number
            % M1: The matrix M1
            %
            % Output:
            % A: The matrix with A_{i,j} = M_2(i,j)

            M = (1i / 2) * besselh(0, mu * sl.norm_diff) .* sl.curve.n_xps;

            M2 = M - M1 .* sl.A_log;

            % Now fix the NaN diagaonal by the real diagonal
            M2_diag = ((1i/2) - (sl.gamma / pi) - (1/pi) * log((mu/2) * sl.curve.n_xps)) .* sl.curve.n_xps;

            % Grab the indices of the diagonal of M2
            s     = size(M2);
            index = 1:s(1)+1:s(1)*s(2); 

            M2(index) = M2_diag;

        end

        function mat = sol_rep_mat(sl, mu, x_0)
            % Create the matrix corresponding to the solution representation at the points x_0
            % Inputs: 
            %   mu  : mu^2 = k^2.
            %   x_0 : A 2 x n vector consisting of the points to be evaluated at
            % Outputs:
            %     mat: A n x (length(density)) matrix such that mat * density = SL.x_0

            % Compute i/4 H_0^(1)(mu|x(t) - x(\tau)|)|x'(\tau)| as a matrix with rows corresponding 
            % to each x point to be valuated at, and columes correspoding the tau_i values.
            m = size(x_0, 2);
            mat = zeros(m, sl.curve.N);
            
            for ii = 1:m
                % Divide by 4 because we are not following the CK convention.
                temp = (1i/ 4) * besselh(0, mu * vecnorm(x_0(:,ii) - sl.curve.xs));

                mat(ii,:) = temp .* sl.curve.n_xps;
            end

            mat = (pi / sl.curve.n) * mat;
        end

        function mat = dkbie_mat(sl, mu)
            % Function which computes the first derivative of the bie mat with
            % respect to mu
            % Inputs:
            %   mu  : The wave number mu.
            
            % Compute the derivative of the Hankel Function
            mat = besselh(1, mu * sl.norm_diff) .* sl.norm_diff;

            % Now change the diagonal to be the correct one
            s     = size(mat);
            index = 1:s(1)+1:s(1)*s(2); 
            mat(index) = -(2 / pi) * 1i;
            
            % Finally multiply by the green's function constant,
            % the jacobian of the curve, and the trapezoidal rule constant.
            mat = ((1i/4) * (pi/ sl.curve.n)) * mat .* sl.curve.n_xps;
        end
       
        function [sols, x_vals, y_vals] = inner_eval(sl, mu, range, n)
            % DEPRECIATED: Was needed for an old idea, but not maintained
            % or used any more.
            % Function which evaluates the integral equation at many points inside the curve. Assuming the curve
            % is centered at zero so it can be shrunk easily.
            % Inputs:
            %   mu      : The wave number. 
            %   range   : array of values contained in (0,1). The curve will be shrunk by this much
            %             and plotted around. WARNING: This requires the curve to be centered at 0
            %   n       : How many points to calculate around each subcurve
            %
            % Outputs:
            %   sols   : A solution array of size length(range) x n. Rows correspond to each curve, cols to points on a curve.
            %   x_vals : Array same size as sols. x values for each point 
            %   y_vals : Array same size as sols. y values for each point
        
            sols   = [];
            x_vals = [];
            y_vals = [];
        
            t_vals = 0:(2 * pi /n):(2 * pi - (1/n));
        
            for ii = range
                % Get new points by shrinking the curve
                points = ii * sl.curve.x_t(t_vals);
                
                x_vals = [x_vals; points(1,:)];
                y_vals = [y_vals; points(2,:)];
        
        
                sols = [sols; sl.eval_point(mu, points)];
            end     
        end % inner_evals

    end % Methods

end % Classdef
