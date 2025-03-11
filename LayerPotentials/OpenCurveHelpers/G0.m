function Res=G0(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          G0.m (function)                                %
%           Author: Stephane Lintner, ACM Caltech 2006                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is to describe the regular part in H^1_0(x)
%Used by Single_Layer.m

%Euler's constant
C=-psi(1);

%Coefficients for Taylor expansion when x is very small
a0=1+2*C*i/pi;
a2=1-2*i*(C-1)/pi;
a4=1-i*(3*C-2)/pi;


[N,M]=size(x);
Res=zeros(N,M);

%Regular definition for points not too small
Index=(abs(x)>0.0000001);
Res(Index)=besselh(0,1,x(Index))-i*2/pi*log(x(Index)/2).*besselj(0,x(Index));

%Taylor expansion for the other points
Index=logical(1-Index);
Res(Index)=a0+a2*x(Index).^2/4+x(Index).^4/64*a4;
    
