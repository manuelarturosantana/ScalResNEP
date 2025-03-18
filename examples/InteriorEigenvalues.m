% Script which demonstrates computing real eigenvalues of the kite using
% the real line adaptive algorithm

% Limits of the interval.
xmin = 1; xmax = 10;

% Initialize the kite and integral operator
% Note for the Rocket use
% curve =  Missile(n)
% n = 150;
n = 16; %Use this for the low tolerance values.
% The syntax (n,[],true) uses the kite parameters used in the paper,
% although the implementation allows for more general parameters.
curve = Kite(n,[],true);
% dl = DoubleLayer(curve);
sl = SingleLayer(curve);





% Scalarized resolvent creation
ut = rand(1,2*n,'like',1i); v = rand(2*n,1,'like',1i);
f = @(z) ut*(sl.lp_mat(z)\v);

% If one desires to use the double layer rather than the single layer
% it can be done as follows.
% For interior problems the sign is switched for the bie_mat
% I = eye(2*n, 2*n);
% F = @(z) (1/2) * I - dl.lp_mat(z);
% f = @(z) ut*(F(z)\v);




%% Compute the poles from the Scalarized Resolvent.

figure(1)
clf

tic
% [pol, nfe] = aaa_recursive1d(f,xmin,xmax,200,1e-2,1e-11);
[pol, nfe] = aaa_recursive1d(f,xmin,xmax,500,1e-1, 1e-1); % Use this for the low tolerance example.
toc

pol = sort(pol,'ComparisonMethod','real');

%% Plot the poles
figure(1)
hold on
plot(pol,'^');
legend('off')
hold off

%% Perform the secant method on each pol to refine it.
tic
warning('off')
pols_sec = [];
max_err2 = [];
for ii = 1:length(pol)
    [poln,~,~,xdiff] = secant_method(@(z) 1/f(z),pol(ii),pol(ii) + 1e-5,4,1e-13)
    pols_sec(ii) = poln;
    max_err2(ii) = xdiff;
end
max(max_err2)
toc

figure(1)
hold on
plot(pols_sec,'ks')
hold off



%% Now compute the eigenvectors from the nullspace computation.
vs = [];
for ii = 1:length(pols_sec)
    % We use bie_mat or lp mat for exterior problems.
    evec = comp_ns((1/2) * I - dl.lp_mat(pols_sec(ii)));
    vs = [vs, evec];
    % See how close to the zero vector the computed nullspace vector gives.
    norm(((1/2) * I - dl.lp_mat(pols_sec(ii))) * evec)
end


%% Compute and visualize the eigenfunctions.
% Two notes:
%       -This is very slow, due to a lack of an acceleration method like
%        IFGF or FMM
%       - To get smooth plots like in the paper you will need to increase
%          numx and numy. To not be singular near the boundary you will
%          need to increase the number of input points in the double layer.
figure(2)
clf
tiledlayout(2,5)

for pind = 1:10
    nexttile
    % Here we get the eigenfunctions
    k = pols_sec(pind);
    xmin =-2; xmax= 2; numx = 100;
    ymin = -2; ymax = 2; numy = 100;
    sol = vs(:,pind); 
    tic
    [x,y,vals] = gen_sol(sol,k,lp,[xmin, xmax], [ymin,ymax],numx,numy, "interior");
    toc
    
    hold on
    pcolor(x,y,abs(real(vals)))
    plot([curve.xs(1,:),curve.xs(1,1)],[curve.xs(2,:),curve.xs(2,1)],'k', 'Linewidth',2);
    shading interp;
    colormap(jet)
    axis square
    title("k = " + num2str(real(k),11))
    hold off
    
end


%% DELETE ME
function [pol, nfevals] = aaa_recursive1d(f, xmin, xmax, n,imag_tol,aaa_tol)
% Function which uses the aaa search recursively on the real line to find the eigenvalues.
% Inputs:
%   f : The scalarized resolvant.
%   xmin/xmax : Lower and upper bounds to search for real eigenvalues in.
%   n : The number of points on the real line to be evaluated in each interval. 
%       Warning: In the current implementation n should be odd
%       
%   imag_tol: Tolerance on the absolute value of the imaginary part of the eigenvalues used to
%             determine if an eigenvalue is an eigenvalue or not.
%   aaa_tol: Tolerance on the aaa algorithm
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
    [~,pol] = aaa(vals,Z, 'tol', aaa_tol);
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



