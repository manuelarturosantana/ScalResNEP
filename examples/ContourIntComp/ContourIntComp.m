% Script to do the comparision between this method and the contour
% integration methods. The way this loop is writen isn't the fastest
% version, but it works.
clear
rng(2)

n = 720;
curve = Kite(n/2,[],true);

lp = DoubleLayer(curve);
ncs = [40:2:140]; % The number of points along the countour
errsc = []; errsc2 = []; errsa = [];
errsc_s = []; errsa_s = []; errsa_r = [];
F = @(z) lp.bie_mat(z);
ell = 1; r = 1;
L = rand(n,ell); R = rand(n,r);

gam = 3 - 1.5i; rad = 1;
% The eigenvalue of the kite. See compute kite evals.m
true_roots = 2.299005732126516 - 1.597683594805217i;
for nc = ncs
    nc
    w = exp(2i*pi*(1:nc)/nc); z = gam+rad*w; % unit roots and quad pts 
    f = @(Z) L'*(F(Z)\R);
    vals = [];
    for k = 1:length(z), vals(k) = f(z(k)); end
    [evs,z]    = blockSS(nc,vals,  gam, rad,10);
    
    % 
    [r,pol] = aaa(vals,z,'tol',1e-12);

    pol = pol(inpolygon(real(pol),imag(pol),real(z), imag(z)));
    
    [evs2,z] = beyn1(nc,F,gam,rad,5,5);

    c_err = 0; c_err2 = 0; a_err = 0;
    c_err  = max(min(abs(evs - true_roots)), c_err);
    c_err2 = max(min(abs(evs2 - true_roots)), c_err2);
    a_err  = max(min(abs(pol - true_roots)), a_err);
    errsc = [errsc, c_err];
    errsc2 = [errsc2, c_err2];
    errsa = [errsa,a_err];
    
    % Now we do the secant method step
    polsn = [];
    for pp = 1:length(pol)
        polsn(pp) = secant_method(@(z) 1/f(z),pol(pp),pol(pp)+1e-3,3,1e-15,1);
    end

    % Now we do the aaar step
    polsr = [];
    for pp = 1:length(pol)
        polsr(pp) = aaar(f,pol(pp),1e-5, 4,1);
    end

    evsn = [];
    for jj = 1:length(evs)
        evsn(jj) = secant_method(@(z) 1/ f(z),evs(jj), evs(jj) + 1e-3,3,1e-15,1);
    end

    cs_err = 0; as_err = 0; ar_err = 0;
        cs_err = max(min(abs(evsn - true_roots)), cs_err);
        as_err  = max(min(abs(polsn - true_roots)), as_err);
        ar_err  = max(min(abs(polsr - true_roots)), ar_err);
    errsc_s = [errsc_s, cs_err];
    errsa_s = [errsa_s, as_err];
    errsa_r = [errsa_r, ar_err];
end
%%
figure(2)
clf
hold on
plot(z)
plot(pol, '^')
plot(true_roots, '*')
plot(evs,'o')
plot(evs2,'.r')
hold off


%%
lw = 3;
figure(1)
clf
hold on

plot(ncs, log10(errsc),'--','LineWidth',lw)
plot(ncs,log10(errsc2),'LineWidth',lw)
plot(ncs, log10(errsa),'LineWidth',lw)
plot(ncs,log10(errsa_s),'-.','LineWidth',lw)
plot(ncs,log10(errsa_r),'LineWidth',lw)
xlabel("Number of frequencies $k$ used along $C_2$",'FontSize',20, 'Interpreter','latex')
ylabel("Log10 Error",'FontSize',20, 'Interpreter','latex')
legend(["block SS", "Beyn 1", "Algorithm 1", "Algorithm 1 with Secant", "Algorithm 1 with Local AAA"],'FontSize',20, 'Interpreter','latex')
hold off
savefig("root" + num2str(rootind));
