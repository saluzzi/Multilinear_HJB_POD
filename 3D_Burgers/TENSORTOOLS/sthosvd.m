function [A, p, Sset] = stmlsvd(T,tau,p)
% STMLSVD  Sequentially truncated multilinear singular value decomposition.
%
%       A = stmlsvd(T,mlr,p) computes the ttensor A, as rank-mlr ST-MLSVD 
%       approximation corresponding to the order p of the tensor T.
%
%       A = stmlsvd(T,mlr) computes the ttensor A, as rank-mlr ST-MLSVD
%       approximation of T, corresponding to the order of the strongest
%       reduction of the size of the partially truncated core tensor.
%
%       NOTE: only supports dense unstructured tensors efficiently. Structured 
%       tensors may loose their structure after the first projection step.

% Sanity check.
%mlr = min( size(T), mlr );
mlr = size(T);    
% Determine, approximately, the most efficient processing order.
if nargin < 3
    ub = size(T);
    p = [];
    td = 1 : length(mlr);
    for i = 1 : length(mlr)
       [~,nextp] = min(arrayfun(@(j) ub(j)*prod(ub(td))*prod(mlr(p)), td));
       p = [p td(nextp)];
       td = [td(1:nextp-1) td(nextp+1:end)];
    end
end

% Factor matrices and core tensor.
d = ndims(T);
U = cell(1,d);
S = T;
Sset = zeros(d,ub(1));
% ST-MLSVD algorithm.
for k = 1 : d
    m = p(k);
    
    % Compute dominant subspace.
    Sm = tenmat(S,m);
    Sm = Sm.data;
    [Q,ss,~] = svd(Sm,'econ');
    l = sum(diag(ss)/ss(1,1) > tau);
    Sset(k,1:l) = diag(ss(1:l,1:l))';
    U{m} = Q(:,1:min(size(Q,2), l));
    % Project onto dominant subspace.
    S = ttm(S, U{m}, m, 't');
end

A = ttensor(S,U);
end