function Res=Scattered_Field(Sol,Ssol,x0,y0,Curve,k,Problem_Flag,Prec_Flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Scattered_Field.m (script)                          %
%           Author: Stephane Lintner, ACM Caltech 2006                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given the solution (ie the density), compute the scattered field at
%(x0,y0)


if(strcmp(Problem_Flag,'Dirichlet'))
    
     if(strcmp(Prec_Flag,'None') || strcmp(Prec_Flag,'Calderon'))
        Res=Single_Extern(Sol,x0,y0,Curve,k);
     end;
            
    if(strcmp(Prec_Flag,'Colton-Kress Combined'))
       Res=Combined_Dirichlet_Extern(Sol,x0,y0,Curve,k);
    end;
    
else
    
    if(strcmp(Prec_Flag,'None'))
        Res=Double_Extern(Sol,x0,y0,Curve,k);
    end;
    
    if(strcmp(Prec_Flag,'Calderon'))
        Res=Double_Extern(Ssol,x0,y0,Curve,k);
    end;
    
    if(strcmp(Prec_Flag,'Colton-Kress Combined'))
        Res=Combined_Neumann_extern(Sol,Ssol,x0,y0,Curve,k);
    end;
    
    
    
end;
