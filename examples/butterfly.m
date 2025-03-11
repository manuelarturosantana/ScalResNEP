%% Butterfly NEP example from the paper.
% Note Matlab may issue many warnings due to the secant method
% requiring inversion of near singular matrices.


[coeffs,fun,F] = nlevp('butterfly');

rng(1)
% Compute the "true" eigenvalues by linearization
e = polyeig(coeffs{:});
dim = size(coeffs{1},1); ut = randn(1,dim); v = randn(dim,1);
f = @(z) ut*(F(z)\v);

xmin = -2; xmax = 2; ymin = -2; ymax = 2;
figure(1)
clf

rp = rec_params('xlims',[xmin,xmax],'ylims',[ymin,ymax],'n',50,'aaa_tol',1e-13);
pol = aaa_recursive(f,rp);


% Check the max error
max_err = 0;
for ii=1:length(pol)
    max_err = max(min(abs(pol(ii) - e)), max_err);
end



% Plot the figure
figure(1)
MS = 'markersize';
hold on
plot(e,'.r',MS,7)
plot(pol,'ok',MS,7,'linewidth',1.5)
ax = gca();
ax.FontSize = 16;
axis square
title("Max Error " + num2str(max_err))
hold off



