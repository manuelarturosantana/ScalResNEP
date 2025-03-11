function [x,y,xp,yp,xpp,ypp]=Curve_Parametrization(t,Curve_Flag,Edge_Flag,ap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Curve_Parametrization.m (function)                %
%           Author: Stephane Lintner, ACM Caltech 2006                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This defines all the curves available in this code.
%To add a new curve, just copy and paste and one already available, and change the name. 
%Unfortunately, this design is not optimal, so to add (or even modify a
%parametrization) it is necessary to also update the following files: 
%1) DIFFARCS.fig, adding an extra field to the Curve Button with the new
%   string. 
%2) Test_Distance.m  if the curve is closed.
%3) Build_Curve.m, if the curve is closed, add it to the list of closed
%ones. 




ON=1;
OFF=0;

if(Edge_Flag==ON)
s=cos(t);
end;


if(nargin==5)
    ap=0.5;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HERE ARE THE CLOSED CURVES
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(Curve_Flag,'circle'))
    
  x=cos(t);
  y=sin(t);
  xp=-y;
  yp=x;
  xpp=-x;
  ypp=-y;
      
end;%End circle

 
  
if(strcmp(Curve_Flag,'ellipse') )
      
      a=2;
      b=1;
      
      x=a*cos(t);
      y=b*sin(t);
      xp=-a*sin(t);
      yp=b*cos(t);
      xpp=-a*cos(t);
      ypp=-b*sin(t);
      
end;%End ellipse
  

  if(strcmp(Curve_Flag,'kite'))
     
    %This is the Colton and Kress kite  
    x=cos(t)+0.65*cos(2*t)-0.65;
    y=1.5*sin(t);
    xp=-sin(t)-2*0.65*sin(2*t);
    yp=1.5*cos(t);
    xpp=-cos(t)-4*0.65*cos(2*t);
    ypp=-1.5*sin(t);
    
   %Replace by 'my' kite if you want a family that tends to the parabola
   %See thesis
   
   % A=0.2;
   % C=cos(t);
   % C2=cos(2*t);
   % S2=sin(2*t);
   % Ss=sin(t);
   % x=(1-A)*C+A*C2;
   % y=Ss;
   % xp=-(1-A)*Ss-2*A*S2;
   % yp=C;
   % xpp=-(1-A)*C-4*A*C2;
   % ypp=-Ss;
    
end;

%Missile (custom built...)
if(strcmp(Curve_Flag,'missile'))
   
   %theta0=-0.9;
   %tmin=theta0+ap/2;
   %tmax=theta0+2*pi-ap/2;
   Dt=2*pi;
   %t=tmin+(s+1)/2*Dt;
   C=cos(t);
   S=sin(t);
   C2=cos(2*t);
   S2=sin(2*t);
   C3=cos(3*t);
   S3=sin(3*t);
   C4=cos(4*t);
   S4=sin(4*t);
   C6=cos(6*t);
   S6=sin(6*t);
   C8=cos(8*t);
   S8=sin(8*t);
   
   R=0.35+0.1*C+0.12*C2+0.15*C3+0.1*C4+0.1*C6+0.05*C8;
  
   Rp=-0.1*S-2*0.12*S2-3*0.15*S3-0.1*4*S4-0.1*6*S6-0.05*8*S8;
   
   Rpp=-0.1*C-4*0.12*C2-9*0.15*C3-0.1*16*C4-0.1*36*C6-0.05*64*C8;
   
   x=R.*C;
   y=R.*S;
   xp=Dt/2*(-R.*S+Rp.*C);
   yp=Dt/2*(R.*C+Rp.*S);
   xpp=Dt^2/4*(-R.*C-2*Rp.*S+Rpp.*C);
   ypp=Dt^2/4*(-R.*S+2*Rp.*C+Rpp.*C);
   
   
   
end;

%%%%%%%%%%%%%%%%%%%%%%%%%
%Here are the open curves
%%%%%%%%%%%%%%%%%%%%%%%%%


if(strcmp(Curve_Flag,'strip'))
    x=s;
    y=s*0; 
    xp=s*0+1;
    yp=0*s;
    xpp=0*s;
    ypp=0*s;
    end;


       
%Spiral
if(strcmp(Curve_Flag,'spirale'))
 
    E=exp(s);
    Cs=cos(5*s);
    Ss=sin(5*s);
    
       x=E.*Cs;
       y=E.*Ss;
       xp=E.*(-5*Ss+Cs);
       yp=E.*(5*Cs+Ss);
       xpp=E.*(-25*Cs-10*Ss+Cs);
       ypp=E.*(-25*Ss+10*Cs+Ss);
    
    end;

    
    
%Half parabola
if(strcmp(Curve_Flag,'parabola'))
    x=1-2*s.^2;
    y=s;
    xp=-4*s;
    yp=1+0*s;
    xpp=-4+0.*s;
    ypp=0*s;

end;



%Cavity
if(strcmp(Curve_Flag,'circular cavity')||strcmp(Curve_Flag,'elliptic cavity') )
 
    a=1;
    b=1;
    
    if(strcmp(Curve_Flag,'elliptic cavity'))
        a=2;
    end;
    
  
    aperture=ap;
   
   tmin=-pi+aperture/4;
   tmax=pi-aperture/4;
   Dt=tmax-tmin;
   
   theta=tmin+(s+1)/2*Dt;
    C=cos(theta);
    S=sin(theta);
    x=a*S;
    y=b*C;
    xp=a*Dt/2*C;
    yp=-b*Dt/2*S;
    xpp=-a*Dt^2/4*S;
    ypp=-b*Dt^2/4*C;
    
end;

%Missile with a small crack
if(strcmp(Curve_Flag,'missile cavity'))
   
    theta0=-0.9;
   tmin=theta0+ap/2;
   tmax=theta0+2*pi-ap/2;
   Dt=tmax-tmin;
   t=tmin+(s+1)/2*Dt;
   C=cos(t);
   S=sin(t);
   C2=cos(2*t);
   S2=sin(2*t);
   C3=cos(3*t);
   S3=sin(3*t);
   C4=cos(4*t);
   S4=sin(4*t);
   C6=cos(6*t);
   S6=sin(6*t);
   C8=cos(8*t);
   S8=sin(8*t);
   
   R=0.35+0.1*C+0.12*C2+0.15*C3+0.1*C4+0.1*C6+0.05*C8;
  
   Rp=-0.1*S-2*0.12*S2-3*0.15*S3-0.1*4*S4-0.1*6*S6-0.05*8*S8;
   
   Rpp=-0.1*C-4*0.12*C2-9*0.15*C3-0.1*16*C4-0.1*36*C6-0.05*64*C8;
   
   x=R.*C;
   y=R.*S;
   xp=Dt/2*(-R.*S+Rp.*C);
   yp=Dt/2*(R.*C+Rp.*S);
   xpp=Dt^2/4*(-R.*C-2*Rp.*S+Rpp.*C);
   ypp=Dt^2/4*(-R.*S+2*Rp.*C+Rpp.*C);
   
end;


    