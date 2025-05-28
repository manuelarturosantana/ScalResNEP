k = 3.831705970207515;

tic, N = 100; MS = 'markersize'; LW = 'linewidth'; FS = 'fontsize'; rng(1)
IN = 'interpreter'; LT = 'latex'; LO = 'location'; SE = 'southeast';
CO = 'color'; HA = 'horizontalalignment'; RT = 'right'; blue = [0 0 1]; red = [1 0 0];
green = [0 .7 0];
v1 = rand(N,1)+1i*rand(N,1); v1 = v1/norm(v1);
v2 = rand(N,1)+1i*rand(N,1); v2 = v2/norm(v2);
ssx = (1:5)'/5; ssy = 1i*(1:3)'/3;
thetavec = 2.^(0:-.25:-6.5);
dx = .5; dy = .5;
for theta = thetavec
xmin = k-1.5*dx; xmax = k+.5*dx; ymin = -1.5*dy; ymax = .5*dy;
Z = [xmin + 1i*ymax + ssx*(xmax-xmin); xmax + 1i*ymax + ssy*(ymin-ymax);
xmax + 1i*ymin + ssx*(xmin-xmax); xmin + 1i*ymin + ssy*(ymax-ymin)];
curve = Build_Curve('circular cavity',N,2*theta); osl = OpenSingleLayer(curve);
f = @(z) v1'*(osl.bie_mat(z)\v2);
F = 0*Z; for j = 1:length(Z), F(j) = f(Z(j)); end
[r,pol] = aaa(F,Z,'tol',1e-11); [~,ii] = sort(abs(pol-k));
p = pol(ii); p = p(1:2);
erreven = p(2)-k; errodd = p(1)-k;
subplot(121)
loglog(theta,-real(erreven),'.b',MS,8), grid on, hold on
loglog(theta,-imag(erreven),'.r',MS,8)
subplot(122)
loglog(theta,-real(errodd),'.b',MS,8), grid on, hold on
loglog(theta,-imag(errodd),'.r',MS,8)
dx = -real(erreven); dy = -imag(erreven);
end
subplot(121)
xlabel('gap angle $\theta$',IN,LT), ylabel('deviation from limit')
title('Second Mode'), loglog(thetavec,.2*thetavec.^2,'--k')
axis([.02 1 1e-15 1]), loglog(thetavec,.2*thetavec.^4,'--k'), hold off
text(.12,1e-5,'$\theta^4$',IN,LT), text(.12,1.5e-2,'$\theta^2$',IN,LT)
text(.3,6e-4,'imag part',FS,8,CO,red), text(.5,2e-1,'real part',FS,8,CO,blue,HA,RT)
subplot(122)
xlabel('gap angle $\theta$',IN,LT), ylabel('deviation from limit')
title('Third Mode'), loglog(thetavec,.01*thetavec.^4,'--k')
axis([.02 1 1e-15 1]), loglog(thetavec,.001*thetavec.^8,'--k'), hold off
text(.12,1e-11,'$\theta^8$',IN,LT), text(.12,1e-5,'$\theta^4$',IN,LT)
text(.3,1e-8,'imag part',FS,8,CO,red), text(.5,3e-3,'real part',FS,8,CO,blue,HA,RT)
