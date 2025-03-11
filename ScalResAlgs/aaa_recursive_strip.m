function [pol, nfevals] = aaa_recursive_strip(f, xmin, xmax, ymin, ymax, n, m, nworkers)
    % Function which performs aaa recursive algorithm on a thin strip in the complex
    % plane. It only divides along the x direction.
    % Inputs:
    %   f :the function to be evaluated at
    %   xmin/xmax : Lower and upper real bounds to search for eigenvalues in.
    %               these do change during the iteration.
    %   ymin/ymax : Lower and upper imaginary bounds to search for eigenvalues in. These 
    %               do not change during the iteration.
    %   n : The number of points to use in x direction, 
    %   m : The number of points to use in the y direction
    %       For an even split make n odd
    %   nworkers: Optional argument to evaluate f in parallel with n
    %             workers. Pass in as 0 to not use parpool
    % Outputs:
    %     pol: The found poles in the area;
    %     nfevals: The total number of function evaluations.

        if nargin < 8
            nworkers = 0;
        end
        
        xs = linspace(xmin,xmax, n); top = xs + 1i * ymax; bottom = xs + 1i * ymin;
        ys = linspace(ymin, ymax, m); left = xmin + 1i * ys; right = xmax + 1i * ys;
        Z = [left(:);bottom(:);right(:);top(:);];
        
        % Evaluate the function at all points
        Ft = eval_f(f,top,nworkers); Fb = eval_f(f,bottom,nworkers);
        Fl = eval_f(f,left,nworkers); Fr = eval_f(f,right,nworkers);
        % Ordering is important here, make sure it goes around nicely for
        % plotting
        vals = [Fl(:);Fb(:);Fr(:);Ft(:);];

        % Plot the box
        figure(1)
        hold on
        plot(Z, 'k')
        hold off
    
        [~,pol] = aaa(vals,Z);

        % Sample the mid point y line so we aren't computing it twice per
        % iteration
        mind = ceil(n / 2); xmid = xs(mind);
        mid = xmid + 1i * linspace(ymin,ymax,n);
        Fm = eval_f(f,mid, nworkers);

        odl.xl = left; odl.xr = mid; 
        odl.xt = top(1:mind); odl.xb = bottom(1:mind);
        odl.Fl = Fl; odl.Fr = Fm; 
        odl.Ft = Ft(1:mind); odl.Fb = Fb(1:mind);
        
        odr.xl = mid; odr.xr = right; 
        odr.xt = top(mind:end); odr.xb = bottom(mind:end);
        odr.Fl = Fm; odr.Fr = Fr; 
        odr.Ft = Ft(mind:end); odr.Fb = Fb(mind:end);
     
        % Sample and pass into the left and right halves.
        [pol1, nfe1] = aaa_rss(f,xmin,xmid, ymin, ymax, n, m, pol, odl, nworkers);
        [pol2, nfe2] = aaa_rss(f,xmid,xmax, ymin, ymax, n, m, pol, odr, nworkers);
        pol = [pol1(:); pol2(:)];
        nfevals = nfe1 + nfe2 + n;
    
    end
    
    function [pol, nfe] = aaa_rss(f, xmin, xmax, ymin, ymax, n, m, prev_pol, od, nworkers)
    % Function to call recursively for the strip todo the AAA search
    % Inputs:
    %   f :the function to be evaluated at
    %   xmin/xmax : Lower and upper bounds to search for real eigenvalues in.
    %   n : The number of points to use in each part. 
    %       Warning: Here we assume that n is odd so we can just 
    %       take the midpoint when recursing
    %   prev_pol: The pols from the previous step
    %   od: Structures with zl,zr,zt,zb, and Fl,Fr, Ft,Fb representing zs
    %       and function values around the rectangle
    %   nworkers: The number of workers for the function evaluation.
    % Outputs:
    %     pol: The found poles in the area;
    %     nfevals: The total number of function evaluations.
        
        % Grab the non-changing y points.
        left = od.xl; right = od.xr; Fl = od.Fl; Fr = od.Fr;
        % Do the 1-d refinement on the top and on the bottom
        % Take the real part
        topr    = (od.xt(2:end) + od.xt(1:end-1)) / 2;
        bottomr = (od.xb(2:end) + od.xb(1:end-1)) / 2;

        Ftr = eval_f(f,topr, nworkers);
        Fbr = eval_f(f,bottomr, nworkers);

        % Now sort into the correct order
        [top, inds]   = sort([od.xt(:); topr(:)], ComparisonMethod="real"); 
        bottom = [od.xb(:); bottomr(:)]; bottom = bottom(inds);
        Ft = [od.Ft(:); Ftr(:)]; Ft = Ft(inds);
        Fb = [od.Fb(:); Fbr(:)]; Fb = Fb(inds);
  
        Z = [left(:);bottom(:);right(:);top(:);];
        

        figure(1)
        hold on
        plot(Z, 'k')
        hold off

        vals = [Fl(:);Fb(:);Fr(:);Ft(:);];
 
        [~,pol] = aaa(vals,Z);
        
        % Reduce the number of previous poles to just the ones inside the box
        prev_pol = insquare(prev_pol,xmin,xmax,ymin, ymax);
        

        % Reduce the poles to only those inside the box
        pol = insquare(pol,xmin,xmax,ymin, ymax);
        
        % check to see if we find new poles or not.
        if length(prev_pol) == length(pol) || isempty(pol)
            nfe = length(Z);
            return 
        % Otherwise we divide and search again.
        else
            % Divide across again and use the real part of the top
            mind = ceil(n / 2); xmid = real(top(mind));
            mid = xmid + 1i * linspace(ymin,ymax,n);
            Fm = eval_f(f,mid, nworkers);
    
            odl.xl = left; odl.xr = mid; 
            odl.xt = top(1:mind); odl.xb = bottom(1:mind);
            odl.Fl = Fl; odl.Fr = Fm; 
            odl.Ft = Ft(1:mind); odl.Fb = Fb(1:mind);
            
            odr.xl = mid; odr.xr = right; 
            odr.xt = top(mind:end); odr.xb = bottom(mind:end);
            odr.Fl = Fm; odr.Fr = Fr; 
            odr.Ft = Ft(mind:end); odr.Fb = Fb(mind:end);

            % Create new prevZ, prevF
            xmid = (xmax + xmin) / 2;
            [pol1, nfe1] = aaa_rss(f,xmin,xmid, ymin, ymax, n, m, pol, odl,nworkers);
            [pol2, nfe2] = aaa_rss(f,xmid,xmax, ymin, ymax, n, m, pol, odr,nworkers);
            pol = [pol1(:); pol2(:)];
            % m for the middle evaluation, 2*n for the top and bottom
            % refinement.
            nfe = nfe1 + nfe2 + m + 2 * n;
        end
    end
    
    function F = eval_f(f, Z, nworkers)
        % Helper function to evaluate F
        parfor (ii = 1:length(Z), nworkers)
            F(ii) = f(Z(ii));
        end
    end

    function val = insquare(val, xmin,xmax,ymin,ymax)
        val = val(real(val) > xmin & real(val) < xmax);
        val = val(imag(val) > ymin & imag(val) < ymax);
    end