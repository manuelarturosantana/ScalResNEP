% Example on subset of cd_player problem from nlevp. It has eigenvalues as
% small as 2e-4 and as large as -41.1399 in the range [-50,5]; 

rng(1)

% If the nlevp is not downloaded, just call the cd_player function from the dependencies folder
% [coeffs,fun,F] = cd_player
[coeffs,fun,F] = nlevp('cd_player');

e = polyeig(coeffs{:});



xmin = -50; xmax = 5;
% Use only the eigenvalues from the linearization in the specified range.
e = e(e < xmax & e > xmin);

dim = size(coeffs{1},1); ut = randn(1,dim); v = randn(dim,1);
f = @(z) ut*(F(z)\v);

% Compute the eigenvalues.
[pol, nfe] = aaa_recursive1d(f,xmin,xmax,250);


figure(1)
hold on
MS = 'markersize';
plot(real(pol),zeros(size(pol)), '.r',MS,7)
plot(e,zeros(size(e)),'ok',MS,7,'linewidth',1.5)
ax = gca();
ax.FontSize = 16;
axis square
hold off

% Check the Error
max_err = 0;
pol_max = 0;
for ii=1:length(e)
    if min(abs(pol - e(ii))) > max_err
        pol_max = e(ii);
        max_err = min(abs(pol - e(ii)));
    end
end
max_err
