% Copyright (c) 2020, Timo Betcke, Nicholas J. Higham, Volker Mehrmann,
% Gian Maria Negri Porzio, Christian Schroeder and Francoise Tisseur.
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.

% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [coeffs,fun,F] = cd_player
%CD_PLAYER      QEP from model of CD player.
%  [COEFFS,FUN,F] = nlevp('cd_player') constructs a 60-by-60 quadratic matrix
%  polynomial lambda^2*M + lambda*D + K arising in the
%  study of a CD player control task.
%  The matrices are returned in a cell array: COEFFS = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle K + lambda*D + lambda^2*M.
%  This problem has the properties pep, qep, real.

%  Reference:
%  Y. Chahlaoui and P. M. Van Dooren. A collection of benchmark examples
%  for model reduction of linear time invariant dynamical systems.
%  MIMS EPrint 2008.22, Manchester Institute for Mathematical Sciences,
%  The University of Manchester, UK, 2002.

load cd_player.mat

n = size(D,1);
M = eye(n);
K2 = K;
coeffs = {K2, D, M};
fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);

end

function varargout = nlevp_monomials(lambda,k)
%NLEVP_MONOMIALS    Evaluate monomials and their derivatives.
%  [F,FP,FPP,...] = nlevp_monomials(LAM,K), for a column vector LAM,
%   generates
%    F = [1, LAM, LAM.^2, ..., LAM.^(K-1), LAM.^K],
%    FP = [0, 1, 2*LAM, 3*LAM.^2, ..., K*LAM.^(K-1)]'  (first derivative)
%  and higher derivatives FPP, ...
%  Thus for scalar LAM: F, FP, FPP, ... are row vectors of length K+1,
%  while for vectors LAM: F, FP, FPP, ... contain a row per element of LAM.
%
%  Define a polynomial P(x) = sum_{i=0:k} c(i)*x^i.
%  Then P(LAM) can be evaluated as P(LAM) = F*c.
%  The second and third derivatives can be evaluated as
%    P'(LAM) = FP*c, P''(LAM) = FPP*c,
%  and so on.

lambda = lambda(:);
n = length(lambda);

f = [ones(n,1),cumprod(repmat(lambda,1,k),2)];

varargout = cell(1,nargout);
varargout{1} = f;
for i=2:nargout
    f = [zeros(n,1),f(:,1:k).*repmat(1:k,n,1)];
    varargout{i} = f;
end
end

function F = nlevp_handleQEP(coeffs)
%NLEVP_HANDLE_QEP  builds a function_handle for QEPs
% F= nlevp_handle_QEP(K,D,M) builds the function_handle
% F = @(z) K + lam.*D + lam.^2.*M, where K, D, M are square
% matrices of the quadratic eigenvalue problem (QEP)
% K + lam*D + lam^2*M.
%
% This is a helper function. It is called by all the QEP problems. It
% should not be called directly by the user.

  F =  @(lam) coeffs{1} + lam.*(coeffs{2} + lam.*coeffs{3});

end