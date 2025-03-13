n = 160;
curve = Kite(n/2,[],true);
lp = CombinedPotential(curve,-4);
lp = DoubleLayer(curve);
gam = 3 - 1.5i; rad = 1;
nc = 200; 
w = exp(2i*pi*(1:nc)/nc); z = gam+rad*w; % unit roots and quad pts 

F = @(z) lp.bie_mat(z);
ell = 1; r = 1;
L = rand(n,ell); R = rand(n,r);
f = @(Z) L'*(F(Z)\R);
vals = [];
for k = 1:length(z), vals(k) = f(z(k)); end
%%
[r,pol] = aaa(vals,z,'tol',1e-10);
%%
[evs,z]    = beyn2adaptive(nc,vals,  gam, rad,6);
[evs2,z,num_r] = beyn1(nc,F,gam,rad,5,5);
%%
pol = pol(inpolygon(real(pol),imag(pol),real(z),imag(z)));
pol = sort(pol,'ComparisonMethod','real');

%%
polsn = [];
for pp = 1:length(pol)
    polsn(pp) = secant_method(@(z) 1/f(z),pol(pp),pol(pp)+1e-3,3,1e-15,1);
end

%%
figure(1)
clf
hold on
plot(pol,'*')
plot(z)
hold off

