% Script which demonstrates computing real eigenvalues of the kite using
% the real line adaptive algorithm

% Limits of the interval.
xmin = 1; xmax = 10;

% Initialize the kite and integral operator
% Note for the Rocket use
% curve =  Missile(n)
n = 150;
% n = 15; Use this for the low tolerance values.
% The syntax (n,[],true) uses the kite parameters used in the paper,
% although the implementation allows for more general parameters.
curve = Kite(n,[],true);
dl = DoubleLayer(curve);


% Scalarized resolvent.
ut = rand(1,2*n,'like',1i); v = rand(2*n,1,'like',1i);
I = eye(2*n, 2*n);

% For interior problems the sign is switched for the bie_mat
F = @(z) (1/2) * I - dl.lp_mat(z);
f = @(z) ut*(F(z)\v);


%% Compute the poles from the Scalarized Resolvent.

figure(1)
clf

tic
[pol, nfe] = aaa_recursive1d(f,xmin,xmax,200,1e-2,1e-11);
% [pol, nfe] = aaa_recursive1d(f,xmin,xmax,500,1e-2, 1e-2); % Use this for the low tolerance example.
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
parfor ii = 1:length(pol)
    [poln,~,~,xdiff] = secant_method(@(z) 1/f(z),pol(ii),pol(ii) + 1e-5,4,1e-13)
    pols_sec(ii) = poln
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



