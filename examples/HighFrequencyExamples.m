% x limits correspond to the real part, a y limits the imaginary part
xmin = 499.9; xmax = 500.1;
ymin = -0.3; ymax = 0;
% Number of points to use in each direction of a rectangle.
numx = 200; numy = 200; 
% Size of the discretization.
n = 4000;

% Make circle with very small opening like in the paper.
numx = 200; numy = 200;
curve = Build_Curve("circular cavity",n,pi / 100);

lp = OpenSingleLayer(curve);
ut = rand(1,n,'like',1i); v = rand(n,1,'like',1i);
f = @(z) ut*(lp.bie_mat(z)\v);

% Note for the rocket shaped cavity use
% curve = Build_Curve("missile cavity",n,0.89522806584082)
% xmin = 399; xmax = 400

% Also for the interior curves follow the script Interior_Eigenvalues.m

%% Perform the non-adaptive scalarized resolvent AAA algorithm on a square.

figure(1) 
clf
tic
xs = linspace(xmin,xmax, numx); top = xs + 1i * ymax; bottom = xs + 1i * ymin;
ys = linspace(ymin, ymax, numy); left = xmin + 1i * ys; right = xmax + 1i * ys;
Z = [left(:);bottom(:);right(:);top(:);];

% Evaluate the function at all points
Ft = eval_f(f,top); Fb = eval_f(f,bottom);
Fl = eval_f(f,left); Fr = eval_f(f,right);

vals = [Fl(:);Fb(:);Fr(:);Ft(:);];
vals = [Fl(:);Fb(:);Fr(:);Ft(:);];   
[~,pol] = aaa(vals,Z, 'tol',1e-12);

eval_time = toc
pol = insquare(pol,xmin,xmax,ymin,ymax);

%% Plot the Poles
figure(1)
clf
hold on
plot(pol,'^');
legend('off')
hold off

%% Secant Method
tic
warning('off')
pols_ref = [];
xdiffs = [];
for ii = 1:length(pol)
    ii
    xdiffs = [];
    [pol_s,~,~,xdiff] = secant_method(@(z) 1/f(z),pol(ii),pol(ii) + 1e-5,4,1e-13);
    xdiffs(ii) = xdiff;
    pols_ref = [pols_ref,pol_s];
end
secant_eval = toc

%% Grab the pole used in the paper, although any of the poles can be used.
% For the rocket use
% [~,ind] = min(pols_ref - 399.9694808817 -0.00434495360);

[~,ind] = min(pols_ref - 499.9073989141- 0.000779974959);
paper_pol = pols_ref(ind)

%% Now sort and grab the eigenvalues

% We use bie_mat or lp mat for exterior problems.
evec = comp_ns(lp.bie_mat(paper_pol));
% Validate the accuracy of the nullspace computation
norm(lp.bie_mat(paper_pol) * evec)


%% Plot and save functions
% Note this is a slow implementation which would be greatly helped by
% use of a frequency domain acceleration method such as FMM or IFGF


% Here we get the eigenfunctions
k = paper_pol;

% Rocket plotting parameters
% xmin = -0.7; xmax= 1.1; 
% ymin = -0.7; ymax = 0.7; 
% numy = 3000; numx = ceil(numy * (1.8/1.2));

% Circle
xmin = -1.1; xmax= 1.1; 
ymin = -1.1; ymax = 1.1; 
numy = 2000; numx = 2000; 
tic
[x,y,vals] = gen_sol(evec,k,lp,[xmin, xmax], [ymin,ymax],numx,numy, "open");
sol_gen_time = toc
%%

figure(1)
clf
hold on
pcolor(x,y,abs(vals))
plot(curve.X,curve.Y,'y', 'Linewidth',2);
shading interp;
colormap(parula)
axis square

title("k = " + num2str(k,11))
hold off


%% Helper functions
function F = eval_f(f, Z)
    % Helper function to evaluate F
    parfor ii = 1:length(Z)
        F(ii) = f(Z(ii));
    end
end

function val = insquare(val, xmin,xmax,ymin,ymax)
    val = val(real(val) > xmin & real(val) < xmax);
    val = val(imag(val) > ymin & imag(val) < ymax);
end

