function [evs,z] = blockSS(nc, Fvals, gam, rad, pbarstart)
% Implementation of the blockSS method fro the survey paper. Modified
% to pick the correct number of moments adaptively as suggested in 
% Beyns paper with his similar algorithm.
% The Non-Linear-Eigenvalue Problem Guttel and Tissuer 2017
% Inputs:
%   nc: Number of Contour points
%   F : F Function evaluated along the contour
% Outputs:
%   evs: The eigenvalues
%   z  : The points along the contour.
    pbar = pbarstart;
    while true
        ell = 1; r = 1;  % probing matrices 
        w = exp(2i*pi*(1:nc)/nc); z = gam+rad*w; % unit roots and quad pts 
         A = zeros(ell,r,2*pbar); % matrices of moments 
        for k = 1:nc 
            Fz = Fvals(k);
            % Fz = L'*(F(z(k))\R); 
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
        % If you look at the diagrams in Beyn of the singular values
        % splitting you see for a low number of contour points the singular
        % values don't split as much. So in practice one could develope a
        % rule of thumb to pick the magic number 1e-12, which may have to
        % be as low as 1e-3 for the smaller eigenvalues. Since that isn't
        % the point of the paper we instead just limit the number of
        % moments to half the amount of the discretization.
        mbar = find(diag(Sig)/Sig(1)>1e-12,1,'last');
        if mbar == pbar && mbar < nc/2
            pbar = pbar + 1;
            continue
        end
        mbar
        V0 = V(:,1:mbar); Sig0 = Sig(1:mbar,1:mbar); 
        W0 = W(:,1:mbar); M = (V0'*B1*W0)/Sig0; evs = eig(M); 
        evs = gam+rad*evs(abs(evs)<1);
        break
    end
end