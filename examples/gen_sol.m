% This script allows us to plot the interior of any open arc
% Stolen and modified from our old friend DIFFARCS
function [x,y,Val] = gen_sol(sol, k, lp,xbds,ybds,numx,numy,problem_flag)
    % Inputs:
    %   sol : The density solution
    %   k   : The wavenumber for the solution
    %   lp  : The layerpotential object. Must have solrepmat defined.
    %   xbds/ybds : Bounds on x to plot on.
    %   numx/numy : Number of x lp.sol_rep_mat(k, [xx,yy]) and y values to plot (grid resolution)
    %   problem_flag: "exterior","interior","open"
    
    warning('off','all');
    xmin = xbds(1); xmax = xbds(2); 
    ymin = ybds(1); ymax = ybds(2);
    Val=zeros(numy,numx);

    %Build the grid
    dx=(xmax-xmin)/numx;
    dy=(ymax-ymin)/numy;
    x=xmin+dx/2:dx:xmax;
    y=ymin+dy/2:dy:ymax;
    
    if strcmpi(problem_flag,"interior")
    % We can scale down assuming our curve is centered at the origin so we 
    % stay away from the boundary. Note there are methods to evaluate the
    % layer potentials close to the boundary, but they are not implemented
    % here.
        cxs = 0.98 * lp.curve.xs(1,:);
        cys = 0.98 * lp.curve.xs(2,:);
    elseif strcmpi(problem_flag,"exterior")
        cxs = lp.curve.xs(1,:);
        cys = lp.curve.xs(2,:);
    elseif strcmpi(problem_flag,"open")
        cxs = 0;
        cys = 0;
    else
        error("Invalid problem flag")
    end
    
    % Track how many function values are inside
    num_int = 0;
    % disp("parfor off")
    parfor n=1:numy
        for m=1:numx
               xx=x(m);
               yy=y(n);
               %test whether inside or not  
               if strcmpi(problem_flag,"interior")
                tt = inpolygon(xx,yy,cxs,cys);
               elseif strcmpi(problem_flag,"exterior")
                   tt = ~inpolygon(xx,yy,cxs,cys);
               else
                tt = Test_Distance(lp.curve,xx,yy);
               end
               
               if(tt==1)
                   num_int = num_int + 1;
                   %Scattered Field at that point
                   if strcmpi(problem_flag,"interior") || strcmpi(problem_flag,"exterior")
                    Val(n,m)= lp.sol_rep_mat(k, [xx;yy]) * sol;
                   else
                     Val(n,m) = Scattered_Field(sol,[],xx,yy,lp.curve,k,"Dirichlet","None");
                   end
               else
                   Val(n,m)=nan;
               end;
        end;
    end;
    non_nan = Val(~isnan(Val));
    Val = Val / norm(non_nan(:));

end
