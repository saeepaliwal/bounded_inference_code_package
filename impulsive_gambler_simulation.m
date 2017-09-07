function impulsive_gambler_simulation
addpath('./tools');
%% Bounded variational logistic regression
% Sadder but wiser
% Optimsm bias
% Overconfidence
clear
nice_colors
run_sim = 1;
if run_sim
    
    N = 100;
    
    % Random attribution
    X = rand(N,1);
    y = binornd(1,0.2,100,1);
    
    % Priors
    a0 = 1e-2;
    b0 = 1e-4;
    
    
    hold on
    % Low is unbounded, high is bounded
    BE = [0.1 1 100];
    for b = 1:length(BE)
        bb = BE(b);
        [w(b), V, invV, logdetV, E_a(b), L] = vb_logit_fit_beta(X, y,a0,b0,bb);
        p(:,b) = 1./(1+exp(-X*w(b)));
    end
    
    Xnew = randn(N,1);
    Xnew = repmat(Xnew,1,3);
    w_mat = repmat(w,N,1);
    pnew =1./(1+exp(-Xnew.*w_mat));
    
    % Specify bet size
    large = 100;
    small = 1;
    winsizes = [0 50 100 150 200];
    
    for cc = 1:100
        clear bet_size account
        account(1,1:3)= N;
        for d = 1:3
            for j = 1:N
                if pnew(j,d) > 0.55
                    bet_size(j,d) = large;
                else
                    bet_size(j,d) = small;
                end
                
                do_win = binornd(1,pnew(j,2));
                if do_win == 1
                    account(j,d) = account(j,d)-bet_size(j,d)+datasample(winsizes,1);
                else
                    account(j,d) = account(j,d)-bet_size(j,d);
                end
                account(j+1,d) = account(j,d);
            end
        end
        all_account(:,:,cc) = account;
    end
else
    load impulsive_gambler_simulation
end

%% Plot simulation
figure(101)
clf
subplot(1,2,1);
hold on
plot(X,y,'o');
plot(X,p(:,1),'Color',red);
plot(X,p(:,2),'Color',green);
plot(X,p(:,3),'k');
subplot(1,2,2);
hold on;
plot(squeeze(all_account(:,1,:)),'Color',red);
plot(squeeze(all_account(:,2,:)),'Color',green);
plot(squeeze(all_account(:,3,:)),'k');
xlim([0 100]);
purty_plot(101,'../figures/BIpaper_Figure6', 'pdf');


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
