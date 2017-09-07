addpath('./tools');

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
scatter(sample(:,1),sample(:,2));
xlim([0 20]);

%%

%clear
%load faith
%y = faith;
y = sample;
%all_betas = [1 5 10 15];
all_betas = [0.0001 0.5 1 100];
figure(110)
clf;
hold on;
clf;
win_mix = {};
for b = 1:4
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
    hold on
    subplot(2,2,b);
    hold on;
    plot_mix(mix,[0 1000 0 1000],1,'r');
    plot(y(:,1),y(:,2),'.','MarkerSize',15);
    win_mix{b} = mix;
end
all_idx

%% Plot
%load workspace2_june6
figure(111);
hold on;

for j = 1:4
    mix = win_mix{j};
    idx = all_idx(j);
    
    subplot(2,2,j);
     for i=1:idx
        plot(mix.state(i).m(1),mix.state(i).m(2),'rx');
     end
     
    hold on;
    plot(y(:,1),y(:,2),'.','MarkerSize',15);
     
    hold on;
    plot_mix(mix,[0 25 -0.5 0.2],1,'r',0.4,0.5);
    ylim([-0.1 0.15]);

end
purty_plot(111,'BIpaper_Figure3','pdf')

    
    