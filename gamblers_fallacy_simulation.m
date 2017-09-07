function gamblers_fallacy_simulation
addpath('./tools');
clear
clf

red = [0.6350    0.0780    0.1840];
blue = [0    0.4470    0.7410];

for i = 1:3001
    r(i) = normrnd(0-i/3000,5);
end
r = r/100;
rc = r +1;
rc = cumprod(rc);
rc= rc-1;
plot(rc)

returns = r';

a_0 = 1;
b_0 = 1;
mu_0 = 1;
la_0 = 10;

start_trade = 1000;
final_trade = 2000;

trade_amt = 10;
% Set beta:
betas = [0.01 1];
trades = [1000 1000];
for b = 1:length(betas)
    x = returns(1:start_trade);
    be = betas(b);
    t = 1;
    for i = start_trade+1:final_trade
        x = [x; returns(i)];
        [F,H,mu_N(t,b),la_N,a_N(t,b),b_N(t,b)]  = fe_vb_1D(x,be,mu_0,la_0,a_0,b_0);
        
        sig_hat = sqrt(1/(a_N(t,b)/b_N(t,b)));
        expret = normrnd(mu_N(t,b),sig_hat);
        if t>1
            % Trade decision algorithm
            if expret > 0.05 % Buy
                trades(t,b) = trades(t-1,b)+trade_amt;
            elseif -0.5 < expret < 0.5 % Hold
                trades(t,b) = trades(t-1,b);
            else
                trades(t,b) = trades(t-1,b)-trade_amt;
            end
        end
        exprets(t,b) = expret;
        t = t+1;
    end
    performance(:,b) = trades(:,b).*returns(start_trade:final_trade-1)/100;
end


%% Financial Simulation: Figure
figure(101)
D = 10000;
M = 20;
N = length(x);
prior = normrnd(mu_0,1/la_0,D,1);
cumexp=exprets;
for j = 1:2
    cumexp(:,j) = cumexp(:,j)+1;
    cumexp(:,j) = cumprod(cumexp(:,j));
    cumexp(:,j) = cumexp(:,j)-1;
end

%%% Plot 1
subplot(2,6,1:2)
hold on;
[n1,x1] = hist(prior,M);
[n2,x2] = hist(x,M);
bar(x1,n1/sum(n1))
bar(x2,n2/sum(n2))
aa = get(gca,'child');
set(aa(1),'FaceColor',red,'EdgeColor','none');
set(aa(2),'FaceColor','k','EdgeColor','none');
xlim([-1 1.5]);
ylim([0 0.2]);

%%% Plot 2
subplot(2,6,3:4);
hold on
sig_hat = sqrt(1/(a_N(end,1)/b_N(end,1)));
[n1,x1] = hist(normrnd(mu_N(end,1),sig_hat,D,1),M);
[n2,x2] = hist(x,M);
[n3,x3] = hist(prior,M);


bar(x1,n1/sum(n1))
bar(x3,n3/sum(n3))
bar(x2,n2/sum(n2))

aa = get(gca,'child');
set(aa(3),'FaceColor',blue,'EdgeColor','none');
set(aa(2),'FaceColor','k','EdgeColor','none');
set(aa(1),'FaceColor',red,'EdgeColor','none');
ylim([0 0.2]);
xlim([-1 1.5]);

%%% Plot 3
subplot(2,6,5:6);
hold on
sig_hat = sqrt(1/(a_N(end,2)/b_N(end,2)));
[n1,x1] = hist(normrnd(mu_N(end),sig_hat,D,1),M);
[n2,x2] = hist(x,M);
[n3,x3] = hist(prior,M);

bar(x3,n3/sum(n3));
bar(x1,n1/sum(n1));
bar(x2,n2/sum(n2));

aa = get(gca,'child');
set(aa(2),'FaceColor',blue,'EdgeColor','none');
set(aa(3),'FaceColor','k','EdgeColor','none');
set(aa(1),'FaceColor',red,'EdgeColor','none');
xlim([-1 1.5]);
ylim([0 0.2]);

%%% Plot 4
subplot(2,6,7:8);

a=1;
%b=length(trades);
b=200;
plot(rc(1001+a:b+1000),'k','LineWidth',2);
xlim([0 b-a]);

%%% Plot 5
subplot(2,6,9:10);
plot(trades(a:b,:),'LineWidth',2)

%%% Plot 6
subplot(2,6,11:12);
plot(cumsum(performance(a:b,:)),'LineWidth',2);
xlim([0 b-a]);


purty_plot(101,'../figures/BIpaper_Figure1', 'png');

function [F,H,mu_N,la_N,a_N,b_N] = fe_vb(x,be,mu_0,la_0,a_0,b_0);

max_iter = 100;
N = length(x);

% Initial values for posterior
a_N = a_0; 
b_N = b_0;
mu_N = 0;

% Initialize FE
F = -inf;

for i = 1:max_iter
    
    la_N = (la_0 + be*N)*a_N/b_N;
    mu_N = ((la_0*mu_0 + be*N*mean(x))/(la_0 + be*N));
    
    a_N = (a_0 + be*N/2);
    b_N = b_0 + (1/2)*(be*sum((x-mu_N).^2) + la_0*((mu_N-mu_0)^2));

    F_old = F;
    F = fe_calc(x,be, a_N, b_N, la_N, mu_N,a_0,b_0,mu_0,la_0);
    
    % Convergence criterion
    if (F - F_old < 10e-4), break; end
    if (i == max_iter), warning('Reached %d iterations',max_iter); end

end
[F,H] = fe_calc(x,be, a_N, b_N, la_N, mu_N,a_0,b_0,mu_0,la_0);

function [F,H] = fe_calc(x,be,a_N, b_N, la_N, mu_N,a_0,b_0,mu_0,la_0)

N = size(x, 1);
d = size(x, 2);
mult=10^(round(log10(1/a_N)));
a_N = a_N*mult;
b_N = b_N*mult;

% < ln p(D|mu,tau)>_q
T1 = (N/2)*(digamma(a_N) - log(b_N) - log(2*pi)) - 0.5*(a_N/b_N)*sum((x.^2) - 2*(mu_N*x) + (mu_N^2) + (1/la_N)); % likelihood

% <ln q(mu)>_q_mu
T2 = 0.5 + 0.5*log(2*pi/la_N);

% <ln q(tau) >_q_tau
T3 = a_N - log(b_N) + gammaln(a_N) + (1-a_N)*digamma(a_N);

% <ln p(mu|tau)>_q
T4 = (N/2)*(log(la_0) + digamma(a_N) - log(b_N) - log(2*pi)) - 0.5*la_0*(a_N/b_N)*((mu_N^2) + (1/la_N) - 2*mu_N*mu_0 + (mu_0^2));

% <ln p(tau)>_q
T5 = a_0*log(b_0)-gammaln(a_0) + (a_0 - 1)*(digamma(a_N) - log(b_N)) - b_0*(a_N/b_N);

% F = < ln p(D|mu,tau)  - ln q(mu) - ln q(tau) + ln p(mu|tau) + ln(p(tau)) >_q
% Likelihood
J = T1;

% KL (constraint)
H =  T2 + T3 - T4 - T5;

F = J - (1/be)*H;



