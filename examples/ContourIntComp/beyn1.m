function [evs,z] = beyn1(nc, F, gam, rad, ell, r)
% Implementation of Beyn's method from the survey paper
% The Non-Linear-Eigenvalue Problem Guttel and Tissuer 2017
% It is adaptive in how it picks the parameters ell and r following
% the method suggested in beyns origonal paper.
% Inputs:
%   nc: Number of Contour points
%   F : F Function evaluated along the contour
%   gam,rad: Center and radius of the circle being used
%   n : size of the NEP
% Outputs:
%   evs: The eigenvalues
%   z  : The points along the contour.
    n = size(F(1),1);
    while true
        L = rand(n,ell); R = rand(n,r);  % probing matrices 
        w = exp(2i*pi*(1:nc)/nc); z = gam+rad*w; % unit roots and quad pts
        pbar=1;
        A = zeros(ell,r,2*pbar); % matrices of moments 
        for k = 1:nc 
            Fz = L'*(F(z(k))\R); 
            for j = 0:2*pbar-1 
                A(:,:,j+1) = A(:,:,j+1) + (w(k)^j*rad*w(k)/nc)*Fz; 
            end 
        end 
        A = A(:,:); 
        B0 = zeros(pbar*ell,pbar*r); 
        B1 = B0; 
        for j = 0:pbar-1 
            B0(1+j*ell:(j+1)*ell,:) = A(:,1+j*r:pbar*r+j*r); 
            B1(1+j*ell:(j+1)*ell,:) = A(:,1+(j+1)*r:pbar*r+(j+1)*r); 
        end 
        [V,Sig,W] = svd(B0); 
        mbar = find(diag(Sig)/Sig(1)>1e-10,1,'last'); 
        if mbar == ell
            ell = ell + 1;
            r   = r + 1;
            continue
        end
        r
        V0 = V(:,1:mbar); Sig0 = Sig(1:mbar,1:mbar); 
        W0 = W(:,1:mbar); M = (V0'*B1*W0)/Sig0; evs = eig(M); 
        evs = gam+rad*evs(abs(evs)<1);
        break
    end
end