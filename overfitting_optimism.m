function overfitting_optimism
addpath('./tools');
nice_colors;

%% Generate machine win trace
n = 10;
start_amount = 10*n;

win_prob = 0.1;
for i = 1:500
    if rand > win_prob
        win_trace(i) = -randi(n,1,1);
    else
        win_trace(i) = 10*randi(n,1,1);
    end
end
account = cumsum([start_amount win_trace]);

% Model inversions

N = length(win_trace);

a_0 = 1;
b_0 = 1;
mu_0 = -1;
la_0 = 1;

betas = [0.0001 10000];
x = account;
for i = 1:length(betas)
    be = betas(i);
    [F,H,mu_N(i),la_N,a_N(i),b_N(i)]  = fe_vb(x,be,mu_0,la_0,a_0,b_0);
    sig_hat(i) = sqrt(1/(a_N(i)/b_N(i)));
end

posterior_low_beta = normrnd(mu_N(1),sig_hat(1),1000,1);
posterior_high_beta = normrnd(mu_N(2),sig_hat(2),1000,1);

% "Within sample" performance

for j = 1:2
    for i = 1:N
        xp = normrnd(mu_N(j), sig_hat(j));
        
        if xp > 200
            bet(i,j) = randi([6 10]);
        else
            bet(i,j) = randi([1 5]);
        end
        if win_trace(i) < 0
            win_trace_within(i,j) = -bet(i,j);
        else
            win_trace_within(i,j) = 10*bet(i,j);
        end 
    end
    account_within(:,j) = cumsum([start_amount; win_trace_within(:,j)]);
end

% Out of sample performance

p_win = 0.05;

for j = 1:2
    for i = 1:length(win_trace)
        xp = normrnd(mu_N(j), sig_hat(j));
        
        if xp > 200
            bet(i,j) = randi([6 10]);
        else
            bet(i,j) = randi([1 5]);
        end
        win_prob = rand;
        if rand>p_win    
            win_trace_oos(i,j) = -bet(i,j);
        else
            win_trace_oos(i,j) = 10*bet(i,j);
        end 
    end
    account_oos(:,j) = cumsum([start_amount; win_trace_oos(:,j)]);
end

%% Generate figure
x = win_trace;
figure(101)
clf
D = 100;
M = 10;
al = 0.8;
prior = normrnd(mu_0,1/la_0,D,1);

subplot(1,3,1)

hold on;
[n2,x2] = hist(x,M);
h2 = bar(x2,n2/sum(n2));
set(h2,'FaceColor','k','EdgeColor','none','FaceAlpha',al);
set(gca,'YTick',[0 0.01 0.02 0.03],'YTickLabels',{'0.1';'0.2';'0.3';'1'})
plot([-20 -15], [0.025 0.027],'k','LineWidth',2);
plot([-20 -15], [0.0262 0.0282],'k','LineWidth',2);
xlim([-20 110]);
ylim([0 0.03]);
xlabel('Win amount')

subplot(1,3,2);
hold on
plot(account,'k','LineWidth',2)
%title('B                                                          ');
xlabel('Games')
ylabel('Account value')

subplot(1,3,3);
hold on
[n2,x2] = hist(account,M);
[n3,x3] = hist(posterior_low_beta,M);
[n4,x4] = hist(posterior_high_beta,M);

h2 = bar(x2,n2/sum(n2));
h4 = bar(x4,n4/sum(n4));
h3 = bar(x3,n3/sum(n3));

set(h2,'FaceColor','k','EdgeColor','none');
set(h3,'FaceColor',dark_grey,'EdgeColor','None');
set(h4,'FaceColor',grey,'EdgeColor','None','FaceAlpha',al);

set(gca,'YTick',[0.1 0.2 0.3 0.4],'YTickLabels',{'0.1';'0.2';'0.3';'1'})
plot([-200 -150], [0.32 0.35],'k','LineWidth',2);
plot([-200 -150], [0.335 0.365],'k','LineWidth',2);

xlim([-200 800]);
ylim([0 0.4]);
xlabel('Win amount')
legend 'True' 'Unbounded' 'Bounded'

% subplot(1,3,3);
% hold on
% plot(account_within(:,1),'Color',green,'LineWidth',2);
% plot(account_within(:,2),'Color',red,'LineWidth',2);
% xlabel('Games')
% ylabel('Account value')
% 
% subplot(2,6,10:12);
% hold on
% plot(account_oos(:,1),'Color',green,'LineWidth',2);
% plot(account_oos(:,2),'Color',red,'LineWidth',2);
% xlabel('Games')

%% Print figure 
%purty_plot(101,'../figures/Figure1','eps');

%% Persistence in gambling
betas = [0.0001 10000];
x_all = account_within;

for j = 1:500
    x_all = [x_all; account_oos(1:j,:)];
    for i = 1:length(betas)
        be = betas(i);
        x = x_all(:,i);
        [F,H,mu_NP(j,i),la_NP,a_NP(i),b_NP(i)]  = fe_vb(x,be,mu_0,la_0,a_0,b_0);
        sig_hat_P(j,i) = sqrt(1/(a_N(i)/b_N(i)));
        position(j,i) = normrnd(mu_NP(j,i),sig_hat_P(j,i));
    end
end


