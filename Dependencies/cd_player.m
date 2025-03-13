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
