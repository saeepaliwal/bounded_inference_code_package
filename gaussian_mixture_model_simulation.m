function gaussian_mixture_model_simulation

addpath('./tools');
addpath(genpath('~/Dropbox/Doctorate/tools/spm12'));
nice_colors
%%
aa = mvnrnd([10 10],[1 1],100);
bb = mvnrnd([1 1],[1 1],100); 
cc = mvnrnd([1 1],[5 5],50);
dd = mvnrnd([10 10],[10 10],50);
ee = mvnrnd([6 5],[5 5],50);
ff = mvnrnd([3 12],[1 1],100);
sample = [aa; bb; cc; dd; ee; ff];
sample(:,1) = sample(:,1)+5;
sample(:,2) = (sample(:,2)-5)/100;

y = sample;
all_betas = [0.0001 1 10000];
win_mix = {};
for b = 1:3
    BE = all_betas(b);
    all_mix = {};
    for n = 1:15
        all_mix{n} = gmm_beta(y,n,0,BE);
        all_f(n) = all_mix{n}.fm;
    end
    
    idx = find(all_f==max(all_f));
    all_idx(b) = idx;
    mix = all_mix{idx};
    
    for i=1:idx
        plot(mix.state(i).m(1),mix.state(i).m(2),'rx');
    end

    win_mix{b} = mix;
end
all_idx

% Plot

figure(2);
clf
hold on;

for j = 1:3
    mix = win_mix{j};
    idx = all_idx(j);
    
    subplot(1,3,j);
     for i=1:idx
        plot(mix.state(i).m(1),mix.state(i).m(2),'kx');
     end
     
    hold on;
    plot(y(:,1),y(:,2),'.','Color',grey,'MarkerSize',10);
     
    hold on;
    plot_mix(mix,[0 25 -0.5 0.2],1,'k',0.4,0.5);
    ylim([-0.1 0.15]);

end
keyboard
%purty_plot(2,'../figures/GMM_Simulation','eps')

    
    