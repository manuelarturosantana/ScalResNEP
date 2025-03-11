
% Script Demonstrating open curve computation of eigenvalues and
% eigenfunctions.

% x and y limits define the rectange to compute the eigenvalues inside of
% with x representing the real part and y representing the imaginary part.
xmin = 1; xmax = 10;
ymin = -0.2; ymax = 0;


n = 200;
ap = pi/4; % Controls the size of the opening.
curve = Build_Curve("circular cavity",n, ap);
lp = OpenSingleLayer(curve);
ut = randn(1,n); v = randn(n,1);
f = @(z) ut*(lp.bie_mat(z)\v);


%% Compute the eigenvalues.

figure(1)
clf

tic
[pol, nfe] = aaa_recursive_strip(f,xmin,xmax,ymin,ymax,50, 50);
run_time = toc

pol = sort(pol,"ComparisonMethod","real");

%% Plot the eigenvalues
figure(1)
clf
hold on
plot(pol,'^');
legend('off')
hold off


%% Secant method refinement of eigenvalues
tic
warning('off')
pols_sec = [];
xdiffs = [];
for ii = 1:length(pol)
    [poln,~,~,xdiff] = secant_method(@(z) 1/f(z),pol(ii),pol(ii) + 1e-5,4,1e-13);
    pols_sec(ii) = poln;
    xdiffs(ii) = xdiff;
end
toc
max(abs(xdiffs))

figure(1)
hold on
plot(pols_sec,'ks')
hold off

%% Compute the eigenvectors
num2plot = length(pol);
vs = [];
for ii = 1:num2plot 
    % We use bie_mat or lp mat for exterior problems.
    e = eig(lp.bie_mat(pols_sec(ii)));
    evec = comp_ns(lp.bie_mat(pols_sec(ii)));
    vs = [vs, evec];
    norm(lp.bie_mat(pols_sec(ii)) * evec)
end

%% Plot the eigenfunctions
% Note that this may be slow since an acceleration method such as FMM or 
% IFGF is not being used.

figure(2)
clf

tiledlayout(2,5)
for pind = 1:10
    % Here we get the eigenfunctions
    k = pols_sec(pind);
    xmin =-3; xmax= 3; numx = 150;
    ymin = -3; ymax = 3; numy = 150;


    sol = vs(:,pind); 
    tic
    [x,y,vals] = gen_sol(sol,k,lp,[xmin, xmax], [ymin,ymax],numx,numy, "open");
    toc
    
    % Now we plot the eigenfunction
    nexttile
    hold on

    pcolor(x,y,abs(vals))
    plot(curve.X,curve.Y,'k', 'Linewidth',2);
    shading interp;
    colormap(jet)
    axis square
    title("k = " + num2str(k,11))
    hold off

    axis off
    hold off
end


