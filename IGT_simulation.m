function IGT_simulation
addpath('./tools');
nice_colors;

%% Generate machine win trace

start_amount = 2000;
trials = 109;
win = [100 100 50 50];
loss = [-200 -200 -50 -50];

for i = 1:4
    trace = binornd(1,0.5,trials,1)*win(i);
    trace(trace==0) = loss(i);
    machine{i} = trace;
end

% Run simulation
a_0 = 1;
b_0 = 1;
mu_0 = 0;
la_0 = 10;

betas = [1 1000];

quantiles = [20 40 60 80 100];

for m = 1:500
    for k = 1:length(betas)
        choice = randi(4);
        trials = [];
        for i = 1:100
            game_trace(i,k) = machine{choice}(i);
            deck(i,k) = choice;
            trials = [trials i];

            if length(trials<5)
                x = machine{choice}(trials);
            else
                x = machine{choice}(trials(end-5:end));
            end
            be = betas(k);
            
            [F,H,mu(i,k),la_N,a_N(i,k),b_N(i,k)]  = ...
                fe_vb(x,be,mu_0,la_0,a_0,b_0);
  
            exp_win(i,k) = mu(i,k);
            
            if exp_win(i,k) < 50
                choice = randi(4);
                trials = [];
            end
            
            if ismember(i,quantiles)
                idx = find(quantiles==i);
                adv_cards(idx,m,k) = sum(ismember(deck(i-19:i,k),[3 4]));
                
                disadv_cards(idx,m,k) = sum(ismember(deck(i-19:i,k),[1 2]));
            end
            
        end
        
        account = [start_amount ;game_trace(1:end-1,k)];
        account = cumsum(account,1);
        
        quant_account(:,m,k) = mean(reshape(account,20,5),1)';
        quant_account_std(:,m,k) = std(reshape(account,20,5))';
        
        all_adv_cards(m,k) = sum(adv_cards(:,m,k));
        all_dis_cards(m,k) = sum(disadv_cards(:,m,k));
    end
    

end

%% Generate Plot
figure(1)
clf
subplot(1,3,1)
hold on
errorbar([20 40 60 80 100],mean(quant_account(:,:,1),2), mean(quant_account_std(:,:,1),2),'k','LineWidth',1)
errorbar([20 40 60 80 100],mean(quant_account(:,:,2),2), mean(quant_account_std(:,:,2),2),'k--','LineWidth',1);
xlim([10 110]);
xticklabels({'0-20';'20-40';'40-60';'60-80';'80-100'});
legend 'Healthy' 'Unbounded'
ylabel('Performance')

subplot(1,3,2)
hold on
h = bar([20 40 60 80 100],[ mean(adv_cards(:,:,1),2) mean(adv_cards(:,:,2),2)],'k','LineWidth',1);
set(h(1),'FaceColor','none','EdgeColor','k');
set(h(2),'FaceColor','k','EdgeColor','k');
set(gca,'XTick',[20 40 60 80 100]);
xticklabels({'0-20';'20-40';'40-60';'60-80';'80-100'});
xlim([10 110]);
legend 'Healthy' 'Unbounded'
ylabel('No. advantageous cards');

subplot(1,3,3);
h = bar([mean(all_adv_cards,1); mean(all_dis_cards,1)]');
set(h(1),'FaceColor',grey,'EdgeColor','k');
set(h(2),'FaceColor','k','EdgeColor','k');
legend 'Advantageous' 'Disadvantageous'
xticklabels({'Healthy';'Unbounded'});
purty_plot(1,'../figures/IGT_Simulation','eps');
ylabel('No. Cards')


