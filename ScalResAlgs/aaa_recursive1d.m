function [pol, nfevals] = aaa_recursive1d(f, xmin, xmax, n,imag_tol)
% Function which uses the aaa search recursively on the real line to find the eigenvalues.
% Inputs:
%   f : The scalarized resolvant.
%   xmin/xmax : Lower and upper bounds to search for real eigenvalues in.
%   n : The number of points on the real line to be evaluated in each interval. 
%       Warning: In the current implementation n should be odd
%       
%   imag_tol: Tolerance on the absolute value of the imaginary part of the eigenvalues used to
%             determine if an eigenvalue is an eigenvalue or not.
% Outputs:
%     pol: The found poles in the area;
%     nfevals: The total number of function evaluations.

    Z = linspace(xmin, xmax, n);
    figure(1)
    hold on
    plot([xmin,xmin],[-0.1,0.1],'k')
    plot([xmax,xmax],[-0.1,0.1],'k')
    hold off

    % Evaluate the function at all points
    vals = [];
    for ii = 1:length(Z)
        vals(ii) = f(Z(ii));
    end
    [~,pol] = aaa(vals,Z, 'tol', 1e-11);
    xmid = (xmin + xmax) / 2;
    

    [pol1, nfe1] = aaa_r1ds(f, xmin, xmid, pol, Z, vals,imag_tol);
    [pol2, nfe2] = aaa_r1ds(f, xmid, xmax, pol, Z, vals,imag_tol);
    pol = [pol1(:); pol2(:)];
    nfevals = nfe1 + nfe2 + n;

end

function [pol, nfe] = aaa_r1ds(f, xmin, xmax, prev_pol, prev_Z, prev_F,imag_tol)
% Function to call recursively to do the aaa search
% Inputs:
%   f :the function to be evaluated at
%   xmin/xmax : Lower and upper bounds to search for real eigenvalues in.
%   prev_pol prev_Z prev_F: The pols, Zs, F values from the previous step.
% Outputs:
%   pol: The found poles in the area;
%   nfevals: The total number of function evaluations.


    figure(1)
    hold on
    plot([xmin,xmin],[-0.1,0.1],'k')
    plot([xmax,xmax],[-0.1,0.1],'k')
    % xline(xmin)
    % xline(xmax)
    hold off

    % Grab the appropriate values of Z for this cut. We use a small value
    % to make sure the end points are included. 
    loc = (prev_Z >= (xmin - 1e-14) & prev_Z <= (xmax + 1e-14));
    prev_Z = prev_Z(loc); prev_F = prev_F(loc);

    % Now we compute the midpoints to refine Z
    % Note that this computation relies on Z being in a sorted order
    Zmid = (prev_Z(2:end) + prev_Z(1:end-1)) / 2;
    
    % Evaluate the function at all points
    vals = [];
    for ii = 1:length(Zmid)
        vals(ii) = f(Zmid(ii));
    end
    % Combine previous and new values.
    vals = [vals(:); prev_F(:)];
    Z    = [Zmid(:); prev_Z(:)];
    % Sort, so we can do the the midpoint computation as above.
    [Z, inds] = sort(Z);
    vals = vals(inds);
    [~,pol] = aaa(vals,Z, 'tol', 1e-11);


    % Reduce the previous poles to only the ones near the real line.
    prev_pol = insquare(prev_pol,xmin,xmax,-imag_tol, imag_tol);

    % Reduce the poles to only those inside the box
    pol = insquare(pol,xmin,xmax,-imag_tol, imag_tol);
    
    % Plotting code to make the slides demonstrating the adaptive algo
    % figure(1)
    % hold on
    % plot(real(pol),zeros(size(pol(:))),'.b','markersize',10)
    % hold off


    % check to see if we find new poles or not.
    if length(prev_pol) == length(pol) || isempty(pol)
        nfe = length(Zmid);
        return 
    % Otherwise we divide and search again.
    else
        xmid = (xmax + xmin) / 2;
        [p1, nfe1] = aaa_r1ds(f, xmin, xmid, pol,Z,vals,imag_tol);
        [p2, nfe2] = aaa_r1ds(f, xmid, xmax, pol,Z,vals,imag_tol);
        pol = [p1(:);p2(:);];
        nfe = nfe1 + nfe2 + length(Zmid);
    end
end

function val = insquare(val, xmin,xmax,ymin,ymax)
    val = val(real(val) > xmin & real(val) < xmax);
    val = val(imag(val) > ymin & imag(val) < ymax);
end