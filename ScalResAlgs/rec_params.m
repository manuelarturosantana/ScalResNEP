function rp = rec_params(varargin)
    % Function which generates the parameters for the time domain scattering problem
    % Pass in as key, value pairp.
    % xlims, ylims : Limits on the real (resp imag) part to search on. Will initialize
    %                default [-1,1] for both
    %                Parameters xmin, xmax, ymin, ymax will also be in the
    %                struct
    % n            : Number of points to use on the shorter side. Will force to be odd.
    %                so that the current implementation of the subdivision works. Default 151;
    % ntb, nlr     : Number of points to use in the top bottom and left right subdivision.
    %                Defaults to zero unless passed in. Forced to be odd like n.
    %                 side (resp. and left/right side)
    % use_secant   : If true use the secant method refinement stage after each iteration. Default true
    % x1_diff      : Value added to x0 to create x1 for secant method.
    %                Default is 1e-5;
    % plotls       : If true plot the rectangles during the recursive sterp.
    % fp_tol       : Tolerance on change in x from secant method to be consider a true pole. If
    %                the secant method did not change less than this on its last step, the pole 
    %                is considered false. Default 1e-8
    % msec_its     : Maximum number of secant method iterations. Default 4
    % sec_tol      : Tolerance on step size change of secant method. Default 1e-10;
    % m_rec        : Maximum number of recursive sterp to take. Default inf;
    % aaa_tol      : Tolerance on AAA algorithm. Default 1e-10
    % aaa_ctol     : Tolerance on AAA cleanup.  Default 1e-10
    % cleanupflag  : Method of doing cleanup for aaa. 1 is the standard
    %                residue based way, and 2 is the undocumented pole-zero
    %                distance way. Default 1
    
    
    rp.n               = 151;
    rp.ntb             = false;
    rp.nlr             = false;
    rp.xlims           = [-3,3];
    rp.ylims           = [-3,3];
    rp.use_secant      = true;
    rp.x1_diff         = 1e-5;
    rp.plotls          = true;
    rp.fp_tol          = 1e-8;
    rp.msec_its        = 4;
    rp.sec_tol         = 1e-10;
    rp.m_rec           = inf;
    rp.aaa_tol         = 1e-10;
    rp.aaa_ctol        = 1e-10;
    rp.cleanupflag      = 1;
   
    
    if (mod(nargin, 2) ~= 0)
        error("rec_params, odd number of arguments, pass in key value pairp")
    end
    
    ii = 1;
    while (ii < nargin)
        key = varargin{ii};
        val = varargin{ii + 1};
        if (strcmp(key, 'n'))
            rp.n = val;
        elseif (strcmp(key, 'ntb'))
            rp.ntb = val;
        elseif (strcmp(key, 'nlr'))
            rp.nlr = val;
        elseif (strcmp(key, 'xlims'))
            rp.xlims = val;
        elseif (strcmp(key, 'ylims'))
            rp.ylims = val;
        elseif (strcmp(key, 'use_secant'))
            rp.use_secant = val;
        elseif (strcmp(key, 'x1_diff'))
            rp.x1_diff = 1e-5;
        elseif (strcmp(key, 'plotls'))
            rp.plotls = val;
        elseif (strcmp(key, 'fp_tol'))
            rp.fp_tol = val;
        elseif (strcmp(key, 'msec_its'))
            rp.msec_its = val;
        elseif (strcmp(key, 'm_rec'))
            rp.m_rec = val;
        elseif (strcmp(key,'aaa_tol'))
            rp.aaa_tol = val;
        elseif (strcmp(key,'aaa_ctol'))
            rp.aaa_ctol = val;
        elseif (strcmp(key,"cleanupflag"))
            rp.cleanupflag  = val;
        else
            err_string = strcat("Error: The following is not a valid field ", key);
            error(err_string)
        end
        % Jump by two to skip past the val.
        ii  = ii + 2;
    end

    % Store bounds by an alias.
    rp.xmin = rp.xlims(1); rp.xmax = rp.xlims(2);
    rp.ymin = rp.ylims(1); rp.ymax = rp.ylims(2);

    % Initialize left and right points if not passed in.
    if ~rp.ntb && ~rp.nlr
        rp.ntb = rp.n;
        rp.nlr = rp.n;
    end

    % Use odd number of points so subdivision is easy.
    if mod(rp.ntb,2) == 0
        rp.ntb = rp.ntb + 1;
    end

    if mod(rp.nlr,2) == 0
        rp.nlr = rp.nlr + 1;
    end


end
    