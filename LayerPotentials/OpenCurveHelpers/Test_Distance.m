function Res=Test_Distance(Curve,x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Test_Distance.m (function)                     %
%           Author: Stephane Lintner, ACM Caltech 2006                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For the closed curves, test whether (x,y) lies inside or outside
%NOte: this is sligthly sloppy because if the file Curve_Parametriation.m
%is changed, this needs to be updated accordingly....


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(Curve.Flag,'circle'))
    
r=sqrt(x^2+y^2);    
    
if(r>1)
    Res=1;
else
    Res=0;

end;

end;%end circle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(Curve.Flag,'ellipse'))
    
    a=2;
    b=1;

    r=sqrt((x/a)^2+(y/b)^2);
    
if(r>1)
    Res=1;
else
    Res=0;

end;

end;%end ellipse



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(Curve.Flag,'kite'))

    %First look at y
    if(abs(y)>1.5)
        Res=1;
    %if y is in the possible area, determine for x now    
    else
    t1=asin(y/1.5);  
    t2=pi-t1;
    
    x1=cos(t1)+0.65*cos(2*t1)-0.65;
    x2=cos(t2)+0.65*cos(2*t2)-0.65;
    
    a=min(x1,x2);
    b=max(x1,x2);
    
    if( x<a || x> b)
        Res=1;
    else 
        Res=0;
    end;
        
    
    end;
    
end;%end kite


if(strcmp(Curve.Flag,'missile'))
 
    [t,r]=cart2pol(x,y);
    
    
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
  
if(r>abs(R)*(1.01))
    Res=1;
else
    Res=0;
end;


end;%end missile


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(Curve.Flag,'strip') || strcmp(Curve.Flag,'parabola') || strcmp(Curve.Flag,'spirale') || strcmp(Curve.Flag,'circular_cavity'))
    Res=1;
end;

if(strcmp(Curve.Flag,'circular cavity') || strcmp(Curve.Flag,'elliptic cavity')||strcmp(Curve.Flag,'missile cavity') )
    Res=1;
end;
