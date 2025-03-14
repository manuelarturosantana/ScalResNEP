%% Example of robustness of choice of random vector
% Here we run the butterfly test a number of times switching which random
% vectors are used each time and see if the poles are still found.
[coeffs,fun,F] = butterfly();

e = polyeig(coeffs{:});
found_all = [];
max_errs = [];

xmin = -2; xmax = 2; ymin = -2; ymax = 2;
rp = rec_params('xlims',[xmin,xmax],'ylims',[ymin,ymax],'n',100,'aaa_tol',1e-11,...
    "use_secant",false,"plotls",false);
%%
for ii = 1:1000
    ii
    % Compute the "true" eigenvalues by linearization
  
    dim = size(coeffs{1},1); ut = randn(1,dim); v = randn(dim,1);
    f = @(z) ut*(F(z)\v);
    

    tic
    pol = aaa_recursive(f,rp);
    toc
    
    % Secant Method final stage
    % Note this is here because aaa_recursive is coding slightly different
    % than in the paper. See the file for aaa_recursive for more details.
    parfor pind = 1:length(pol)
        cp = pol(pind);
        [x, ~, ~, x_diff] = secant_method_s(@(z) 1 / f(z), cp,cp+1e-5,...
            4, 1e-13, xmin,xmax,ymin,ymax);
        if x_diff > 1e-8 || isnan(x_diff)
            pol(pind) = nan;
        else
            pol(pind) = x;
        end
    end
    % Remove the nan values.
    pol = rmmissing(pol);


    % Check the max error
    max_err = 0;
    for jj=1:length(e)
        max_err = max(min(abs(e(jj) - pol)), max_err);
    end
    
    max_errs(ii) = max_err;
    npfound(ii) = length(pol);
    found_all(ii) = length(pol) == 256 && abs(max_err) < 1e-8;
end




