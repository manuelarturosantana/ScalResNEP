% Example on subset of cd_player problem from nlevp. It has eigenvalues as
% small as 2e-4 and as large as 0.78 in the range [-1,1]; 

[coeffs,fun,F] = nlevp('cd_player');


rng(1)
e = polyeig(coeffs{:});
xmin = -50; xmax = 5;
e = e(e < xmax & e > xmin);

dim = size(coeffs{1},1); ut = randn(1,dim); v = randn(dim,1);
f = @(z) ut*(F(z)\v);

% figure(1)
% clf

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

max_err = 0;
pol_max = 0;
for ii=1:length(e)
    if min(abs(pol - e(ii))) > max_err
        pol_max = e(ii);
        max_err = min(abs(pol - e(ii)));
    end
    % max_err = max(min(abs(pol - e(ii))), max_err);
end
max_err