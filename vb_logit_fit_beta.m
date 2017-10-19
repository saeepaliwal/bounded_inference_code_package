function [w, V, invV, logdetV, E_a, L] = vb_logit_fit_beta(X, y,a0, b0, BE)
% Adapted from Jan Drugowitsch's Variational Logistic Regression Code
% Copyright (c) 2013, Jan Drugowitsch
% All rights reserved.
%% hyperprior parameters
if nargin < 3,  a0 = 1e-2;  end
if nargin < 4,  b0 = 1e-4;  end

%% pre-compute some constants
[N, D] = size(X);
max_iter = 500;
an = BE*(a0 + 0.5 * D);    gammaln_an_an = gammaln(an) + an;
t_w = 0.5 * sum(bsxfun(@times, X, y), 1)';


%% start first iteration kind of here, with xi = 0 -> lam_xi = 1/8
lam_xi = ones(N, 1) / 8;
E_a = a0 / b0;
invV = BE*E_a * eye(D) + 2 * X' * bsxfun(@times, X, lam_xi);
V = inv(invV);
w = V * t_w;
bn = BE*(b0 + 0.5 * (w' * w + trace(V)));
L_last = - N * log(2) ...
         + 0.5 * (w' * invV * w - logdet(invV)) ...
         - an / bn * b0 - an * log(bn) + gammaln_an_an;


%% update xi, bn, (V, w) iteratively
for i = 1:max_iter;
    % update xi by EM-algorithm
    xi = sqrt(sum(X .* (X * (V + w * w')), 2));
    lam_xi = lam(xi);

    % update posterior parameters of a based on xi
    bn = BE*(b0 + 0.5 * (w' * w + trace(V)));
    E_a = an / bn;

    % recompute posterior parameters of w
    invV = BE*E_a * eye(D) + 2 * X' * bsxfun(@times, X, lam_xi);
    V = inv(invV);
    logdetV = - logdet(invV);
    w = V * t_w;

    % variational bound, ingnoring constant terms for now
    L = - sum(log(1 + exp(- xi))) + sum(lam_xi .* xi .^ 2) ...
        + 0.5 * (w' * invV * w + BE*logdetV - sum(xi)) ...
        - BE*E_a * b0 - BE*an * log(bn) + BE*gammaln_an_an;

    if (L_last > L) || (abs(L_last - L) < abs(0.00001 * L))
        break
    end
    L_last = L;  
end;
if i == max_iter
    warning('Bayes:maxIter', ...
        'Bayesian logistic regression reached maximum number of iterations.');
end

%% add constant terms to variational bound
L = L - BE*gammaln(a0) + BE*a0 * log(b0);

function out = lam(xi)
% returns 1 / (4 * xi) * tanh(xi / 2)
divby0_w = warning('query', 'MATLAB:divideByZero');
warning('off', 'MATLAB:divideByZero');
out = tanh(xi ./ 2) ./ (4 .* xi);
warning(divby0_w.state, 'MATLAB:divideByZero');
% fix values where xi = 0
out(isnan(out)) = 1/8;

function out = logdet(A)
out = 2 * sum(log(diag(chol(A))), 1);
