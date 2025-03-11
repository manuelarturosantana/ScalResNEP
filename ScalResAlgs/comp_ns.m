function vec = comp_ns(mat)
    % Function to compute the nullspace given a matrix which is nearly
    % singular (Singular up to numerical error). We assume that there is
    % only one zero eigenvalue and the nullspace is one dimensional.
    
    % First we compute the LU decomposition.
    [~,U] = lu(mat);
    
    % Now we compute the nullspace as in an introductory linear algebra
    % class to solve Ax = 0. We have freedom, so we pick the last column to
    % be free with value s = -1; 
    b = U(1:end-1,end);
    A = U(1:end-1,1:end-1);
    % Now we solve the smaller linear system. Since we used the negative
    % value we append -1 to the end to get our eigenvector.
    vec = A\b; vec = [vec;-1];
    vec = vec / norm(vec);
end
