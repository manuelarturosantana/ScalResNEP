% Class to convert the Diffarcs code for building and evaluating a single
% layer potential into my LP code, so I don't have to change any other
% code. Admittedly this will be a little messy, and not follow the
% LayerPotential format exactly. Instead it will implement the necessary
% methods for lp_sv in an efficient way. Note that the
classdef OpenSingleLayer < LayerPotential
    properties
        jacobian;      % jacobian of the curve which will be precomputed for speed! Stored as a row vector
    end

    methods
        
        function osl = OpenSingleLayer(curve)
            % Constructor
            %    curve: Curve object from OpenCurves, rather than my curve class.

            osl.curve = curve;
            % Precompute the Jacobian
            osl.jacobian = sqrt(curve.Xp.^2+curve.Yp.^2).';
        end
        
        
        function mat = bie_mat(osl, mu)
            % Use the DIFFARCS function to build the osl mat.
            mat = Build_Single_Matrix(osl.curve, mu);
        end

        function mat = lp_mat(osl, mu)
            mat = bie_mat(osl, mu);
        end

        function val = is_reg(osl)
            val = ~isempty(osl.reg_norm_diff2);
        end

    end
    
end

function Mat=Build_Single_Matrix(Curve,k)
    %This builds the matrix for the single layer operator.
    %Only available for open curves.
    
    warning off;
    
    N=Curve.N;
    Mat=zeros(N,N);
    
    
    %The smooth Kernel K2 first
    %--------------------------
    
    %Distance
    Xv=repmat(Curve.X,1,N);
    Xh=repmat(Curve.X.',N,1);
    X=Xv-Xh;
    Yv=repmat(Curve.Y,1,N);
    Yh=repmat(Curve.Y.',N,1);
    Y=Yv-Yh;
    
    R=sqrt(X.^2+Y.^2);
        
    %Regularizing the log
    s=cos(Curve.T);
    Sv=repmat(s,1,N);
    Sh=repmat(s.',N,1);
    S=Sv-Sh;
    reg=abs(R./S);
    reg(1:N+1:end)=sqrt((Curve.Xp).^2+(Curve.Yp).^2);
    
    %Bessel function for splitting
    Bess=besselj(0,k*R);
    
    %Smooth Part of the Kernel
    Mat1= G0(k*R)+ 2*sqrt(-1)/pi*log(k*reg/2).*Bess;
    
    %Line Element
    Jacobian=repmat(sqrt(Curve.Xp.^2+Curve.Yp.^2).',N,1);
    
    %Trapezoidal Rule
    Mat1=Curve.dt.*Mat1.*Jacobian;
     
    
    %The Singular Part now
    %----------------------
    
    lambda=-pi*[log(2);2./(1:N-1).']/N;
    
    Rs=zeros(2*N,1);
    m=0:N-1;
    l=0:2*N-1;
    l=l.';
    Rs=cos(l*m*pi/N)*lambda;
    
    
    Rn1=toeplitz(Rs(1:N));
    Rn2=zeros(N,N);
    for i=1:N
    Rn2(i,:)=Rs(i+1:i+N).';
    end;
    
    Rn=0.5*(Rn1+Rn2);
    
    %Rn=zeros(N,N);
    %for i=1:N
     %   for j=1:N
    
      %      Rn(i,j)=1/2*(Rs(abs(i-j)+1)+Rs(i+j));
            
       % end;
    %end;
    
    
    Mat2=2/pi*sqrt(-1)*Rn.*Jacobian.*Bess;
    
    Mat=sqrt(-1)/4*(Mat1+Mat2);
end    