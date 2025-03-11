function Curve=Build_Curve(Curve_Flag,N,ap);
%This builds the curve according to the Flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Build_Curve.m (function)                         %
%           Author: Stephane Lintner, ACM Caltech 2006                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Intializes the curve variables. 
%Curve is an object with several fields, defined at the bottom of this
%file. 


ON=1;
OFF=0;


%if no aperture size specified, put it to zero
%The aperture is only useful for cracks.
if nargin==2

    if(strcmp(Curve_Flag,'circular cavity') ||strcmp(Curve_Flag,'elliptic cavity') ||strcmp(Curve_Flag,'missile cavity') )
    answer=inputdlg('Size of Aperture?');
    ap=str2num(answer{1,1});

    else
       ap=0;
    end;
 
 end;


if(strcmp(Curve_Flag,'circle')||strcmp(Curve_Flag,'ellipse') || strcmp(Curve_Flag,'kite') || strcmp(Curve_Flag,'missile'))
    Edge_Flag=OFF;    
else
    Edge_Flag=ON;
end;




%CLOSED CURVES FIRST
%-------------------
if(Edge_Flag==OFF)

  dt=2*pi/N;
  T=0:dt:2*pi-dt;
  T=T.';

else
    
    dt=pi/N;
    T=dt/2:dt:pi;
    T=T.';
    
end;

[X,Y,Xp,Yp,Xpp,Ypp]=Curve_Parametrization(T,Curve_Flag,Edge_Flag,ap);
Curve.N=N;
Curve.dt=dt;
Curve.T=T;
Curve.X=X;
Curve.Y=Y;
Curve.Xp=Xp;
Curve.Yp=Yp;
Curve.Xpp=Xpp;
Curve.Ypp=Ypp;
Curve.Flag=Curve_Flag;
Curve.Edge_Flag=Edge_Flag;
Curve.ap=ap;

%These last four are only used in the matrix computations
Curve.s=cos(T);
Curve.s2=1-(Curve.s).^2;
Curve.Jacobian=sqrt(Xp.^2+Yp.^2);
Curve.OverJacobian=1./(sqrt(Xp.^2+Yp.^2));
