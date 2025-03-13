% Example of computing all the interior eigenvalues of the circle between
% (1,100). Non optimal implementation in terms of speed.

% Load the bessel function roots.
Bessel1_100
xmin = 1; xmax = 100;
e = bjz(bjz < xmax & bjz > xmin);

%%
n = 400;
curve = Circle(n,1);
dl = DoubleLayer(curve);
I = eye(2*n,2*n);
 ut = randn(1,2*n); v = randn(2*n,1);
F = @(z) (1/2) * I - dl.lp_mat(z);
f = @(z) ut*(F(z)\v);

figure(1)
clf

tic
xmid = (xmin + xmax) / 2;
[pol1, nfe1] = aaa_recursive1d(f,xmin,xmid,501);
[pol2, nfe2] = aaa_recursive1d(f,xmid,xmax,501);
pol = [pol1(:); pol2(:)];
nfe = nfe1 + nfe2;
toc

%%
figure(1)
hold on
plot(real(pol),zeros(size(pol)),'^');
plot(e, zeros(size(e)),'bs')
legend('off')
hold off

%%

max_err = 0;
for ii=1:length(pol)
    % Display is a eigenvalue was missed which one it is.
    if min(abs(pol(ii) - e)) > 1e-7
       pol(ii)
    end
    max_err = max(min(abs(pol(ii) - e)), max_err);
    
end
max_err

%%
tic
max_err2 = [];
parfor ii = 1:length(pol)
    poln = secant_method(@(z) 1/f(z),pol(ii),pol(ii) + 1e-5,4,1e-13)
    max_err2(ii) = min(abs(poln - e));
end
max(max_err2)
toc