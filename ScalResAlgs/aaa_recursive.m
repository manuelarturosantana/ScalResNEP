
% Function which computes the adaptive scalarized resolvent algorithm
%  in a recursive fashion applied along rectangles
% by dividing diadically in each dimension. 
% For optimal speed reuses the computations from previous steps. For a full list of inputs
% see rec_params.m
%
% Notes
%  1.) one difference between this implementation and the paper. In this implementation
% the secant method is applied after every recursive step, rather than just at the end. This 
% change makes the method converge faster if AAA is returning many false
% polse. One easy way to run this is to run this as in the paper is to call
% aaa_recursive with secant_method = false, and then call sec_ref.m
%
% 2.) The functions aaa_rec and eval_f in this script contain parallel for loops. If you don't
% have access to the parallel computing tool box search fot the parfor loops and replace with
% the for loops. Note for small problems, like the butterfly this will
% actually be slower. It can be very helpful however for the high frequency
% Laplace eigenvalue problems considered.

function pol = aaa_recursive(f,rp)
    % Function which computes the poles of f using a recursive aaa application within the 
    % bounds of the box specified.
    % Inputs: 
    %   f: The scalar function to find the poles of.
    %   rp: The problem structure from rec_params.m type doc rec_params for more details on
    %       hyperparameters.

    % Initial Approximation is that there is not poles inside.
    prev_pol = [];
    prs.has_prevs = zeros(1,4);

    pol = aaa_rec(f, rp, rp.xmin, rp.xmax, rp.ymin, rp.ymax, prev_pol, prs,0);
    
end

function pol = aaa_rec(f,rp, xmin, xmax, ymin, ymax, prev_pol, prs, rec_num)
    % Recursive function which performs subdivision on AAA boxes.
    % Inputs:
    %   Same as aaa_recursive with the addition of 
    %   prev_pol: The poles from the previous iteration of the algorithm
    %   prs     : The struct containing previous function value information.
    % The box below is the key for side numbers.
    %     1
    %     _
    %   2| |4
    %     -
    %     3

    % Reduce the old poles to be only the ones inside the current square.
    prev_pol = inrectangle(prev_pol,xmin,xmax,ymin,ymax);

    if rec_num == rp.m_rec
        warning("Reached Maximum number of recursive steps")
        pol = prev_pol
        return
    end

    % Compute the points for the four sides.
    [top,    prs]  = create_points(rp.ntb, 1,  xmin, xmax, ymin, ymax, prs);
    [left,   prs]  = create_points(rp.nlr, 2,  xmin, xmax, ymin, ymax, prs);
    [bottom, prs]  = create_points(rp.ntb, 3,  xmin, xmax, ymin, ymax, prs);
    [right,  prs]  = create_points(rp.nlr, 4,  xmin, xmax, ymin, ymax, prs);

    [fvals,Z, top,left,bottom,right,fsl,fsb,fsr,fst]= eval_f(f,rp, prs,top,left,right,bottom);
    

    % link  left, bottom, top, right together in Z
    
    % We use a higher tolerance to be extra careful to avoid spurious poles.
    
    [~,pol] = aaa(fvals,Z,'tol',rp.aaa_tol,'cleanuptol',rp.aaa_ctol,'cleanupflag',rp.cleanupflag);

    % Reduce the number of poles to just the ones inside current box
    pol      = inrectangle(pol,xmin,xmax,ymin,ymax);
    if rp.use_secant
        %parfor pind = 1:length(pol)
        for pind = 1:length(pol)
            cp = pol(pind);
            [x, ~, ~, x_diff] = secant_method_s(@(z) 1 / f(z), cp,cp+rp.x1_diff,...
                rp.msec_its, rp.sec_tol, xmin,xmax,ymin,ymax);
            if x_diff > rp.fp_tol || isnan(x_diff)
                pol(pind) = nan;
            else
                pol(pind) = x;
            end
        end
        % Remove the nan values.
        pol = rmmissing(pol);
    end
    
    % check to see if we find new poles or not.
    % and force it to divide at least one time if no poles are found
    if length(prev_pol) == length(pol) && rec_num ~= 0
        return 
    % Otherwise we divide and search again.
    else
        xmid = (xmax + xmin) / 2; ymid = (ymax + ymin) / 2;
        % Create prs for the subsquare with the previous function values and points.
 
        % Bottom left square
        prs = create_prs(2,left,fsl,3,bottom,fsb);
        p1 = aaa_rec(f,rp, xmin, xmid, ymin, ymid, pol, prs, rec_num + 1);

        % Bottom right square
        prs = create_prs(3,bottom,fsb,4,right,fsr);
        p2 = aaa_rec(f,rp, xmid, xmax, ymin, ymid, pol, prs, rec_num + 1);

        % Top left square
        prs = create_prs(1,top,fst,2,left,fsl);
        p3 = aaa_rec(f,rp, xmin, xmid, ymid, ymax, pol, prs, rec_num + 1);
        
        % Top right square
        prs = create_prs(1,top,fst,4,right,fsr);
        p4 = aaa_rec(f,rp, xmid, xmax, ymid, ymax, pol, prs, rec_num + 1);

        pol = [p1(:);p2(:);p3(:); p4(:)];
    end
end



function [xs,prs] = create_points(n,sind, xmin, xmax, ymin, ymax, prs)
    % Creates the grid for a square to evaluate aaa in
    % n number of points on each side.
    % sind: Which side are we doing according to the square below
    % f   : The function to evaluate.
    % xmin, xmax, ymin, ymax: The corners of the square
    % prs: A structure which has the previous data and function values.
    % new points. The following key gives the number corresponding to side indices.
    %     1
    %     _
    %   2| |4
    %     -
    %     3
    % Output:
    %   xs  The new x values for this side
    %   fs  The function values in the correct order for this side.

    % Check if sind is in previous
    if prs.has_prevs(sind) == 1
        % Split on real part for top and bottom imaginary part for left and right
        if sind == 1 || sind == 3
            % Gotcha! Don't forget to make sure the fs the ones inside the
            % box.
            [xs, inds] = sort(prs.xvals{sind},'ComparisonMethod','real');
            fs_old = prs.fvals{sind}; fs_old = fs_old(inds);
            % Ensure that the xs are only the ones within the box.
            inds = real(xs) > xmin - 1e-10 & real(xs) < xmax + 1e-10;
            xs = xs(inds); fs_old = fs_old(inds);
            % Compute the midpoint of each term
            xsmid = real((xs(2:end) + xs(1:end-1)) / 2);
            
            if sind == 1
                xsmid = xsmid + 1i * ymax;
            else
                xsmid = xsmid + 1i * ymin;
            end
            % Save reused terms.
            [xs_new,inds] = sort([xs, xsmid],'ComparisonMethod','real');
        else
            [xs, inds] = imag_sort(prs.xvals{sind});
            fs_old = prs.fvals{sind}; fs_old = fs_old(inds);
            % Ensure that the xs are the only ones within the box.
            inds = imag(xs) > ymin - 1e-10 & imag(xs) < ymax + 1e-10;
            xs = xs(inds); fs_old = fs_old(inds);

            xsmid = imag((xs(2:end) + xs(1:end-1)) / 2);
            
            if sind == 2
                xsmid = xmin + 1i * xsmid;
            else
                xsmid = xmax + 1i * xsmid;
            end
            [xs_new,inds] = imag_sort([xs, xsmid]);
        end
         % Return only the points to be evaluated.
        xs = xsmid;
        % Save terms for later use.
        prs.fs_old{sind} = fs_old;
        prs.xs_new{sind} = xs_new;
        prs.inds{sind}   = inds;
    else
        switch sind
            case 1
                z = linspace(xmin,xmax,n); xs = z + 1i * ymax;
            case 2
                z = linspace(ymin, ymax,n); xs = xmin + 1i * z;
            case 3
                z = linspace(xmin,xmax,n); xs = z + 1i * ymin;
            case 4
                z = linspace(ymin, ymax,n); xs = xmax + 1i * z;
            otherwise
                error("create_points: side index is incorrect.")
        end
    end
end


function prs = create_prs(s1,xs1,fs1,s2,xs2,fs2)
    % Function which creates a structure containing previous information
    % Inputs:
    %   si,xsi,fsi for two sides. si corresponds to which side of the square these values 
    %              are for. (see create_points)
    % Outputs:
    %   prs: A struct with fields
    %       has_prevs a 1 x 4 array. 0 if contains previous indicies, 1 otherwise.
    %       fvals    a 1 x 4 cell array. Empty if it doesn't have previous values, otherwise contains the previous function values.
    %       xvals    a 1 x 4 cell array. Empty if it doesn't previous values, otherwise contains the previous x values.

    prs.has_prevs = zeros(1,4);
    prs.has_prevs(s1) = 1; prs.has_prevs(s2) = 1;

    prs.xvals = cell(1,4);
    prs.xvals{s1} = xs1; prs.xvals{s2} = xs2;

    prs.fvals = cell(1,4);
    prs.fvals{s1} = fs1; prs.fvals{s2} = fs2;

    % These will be intialized later to contain just the old fvalues on one side, and the
    % all the new xvalues, and also the inds which will sort the old and new f values.
    prs.xs_new = cell(1,4);
    prs.fs_old = cell(1,4);
    prs.inds   = cell(1,4);
end

function [fvals,Z, top,left,bottom,right,fsl,fsb,fsr,fst] = eval_f(f,rp,prs,top,left,right,bottom)
    % Very ugly function (hideous in fact) which allows for evaluation of all the values together.
    % The goal is to increase load balancing when for parallel processing.

    Z = [left(:);bottom(:);right(:);top(:)];


    % Returns a column vector of f evaluated at the points Z
    for ii = 1:length(Z)
    %parfor ii = 1:length(Z)
        vals(ii) = f(Z(ii));
    end

    % Now we go through each side and link up function values and terms
    % we order L B R T as is right for plotting
    if prs.has_prevs(2)
        fsmid   = vals(1:length(left)); vals = vals(length(left) + 1:end);
        fs = [prs.fs_old{2}(:); fsmid(:)];
        fsl = fs(prs.inds{2}); fvals = fsl;
        left = prs.xs_new{2};
    else
        fsl = vals(1:length(left)); fvals = fsl;
        % This is lazy programming to just cut off the part we don't need anymore.
        vals = vals(length(left) + 1:end);
    end

    if prs.has_prevs(3)
        fsmid   = vals(1:length(bottom)); vals = vals(length(bottom) + 1:end);
        fs = [prs.fs_old{3}(:); fsmid(:)];
        fsb = fs(prs.inds{3}); fvals = [fvals(:); fsb(:)];
        bottom = prs.xs_new{3};
    else
        fsb = vals(1:length(bottom)); fvals = [fvals(:); fsb(:)];
        vals = vals(length(bottom) + 1:end);
    end

    if prs.has_prevs(4)
        fsmid   = vals(1:length(right)); vals = vals(length(right) + 1:end);
        fs = [prs.fs_old{4}(:); fsmid(:)];
        fsr =fs(prs.inds{4}); fvals = [fvals(:); fsr(:)];
        right = prs.xs_new{4};
    else
        fsr = vals(1:length(right)); fvals = [fvals(:); fsr(:)];
        vals  = vals(length(right)+1:end);
    end

    if prs.has_prevs(1)
        fsmid   = vals(1:length(top));
        fs = [prs.fs_old{1}(:); fsmid(:)];
        fst = fs(prs.inds{1}); fvals = [fvals(:);fst(:)];
        top = prs.xs_new{1};
    else
        fst = vals; fvals = [fvals(:); fst(:)];
    end

    Z = [left(:);bottom(:);right(:);top(:)];

    if rp.plotls
        figure(1)
        hold on
        plot(Z,'k')
        hold off
    end

end

function [vals, inds] = imag_sort(vals)
    % Function to do comparison method imaginary with sort.
    vs = imag(vals);
    [~,inds] = sort(vs);
    vals = vals(inds);
end