function Res=Single_Extern(phi,x0,y0,Curve,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Single_Extern.m (function)                           %
%           Author: Stephane Lintner, ACM Caltech 2006                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Evaluate the Single-Layer Operator for a point outside the curve
%Works for both, closed and open curves
%NOTE: if  (x0,y0) is too close to the curve, phi might need to be
%zeropadded using the corresponding routine, so as to capture the 
%quasi singular behaviour of the kernel.
%For points ON the curve, use Single_Layer_Point.m instead.


r=sqrt((Curve.X-x0).^2+(Curve.Y-y0).^2);

Jacobian=sqrt(Curve.Xp.^2+Curve.Yp.^2);

Res=i/4*Curve.dt*sum(besselh(0,1,k*r).*phi.*Jacobian);
