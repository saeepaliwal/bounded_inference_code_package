

%%
% Loop through pictures
clear
all_betas = [0.0001 0.5 1 100];
win_mix = {};

for m = 1:10
    pic = imread(['./snowy_pictures_VPT/snowy_pictures_' sprintf('%d',m) '.jpg']);
    m
    n = 1;
    for i = 1:size(pic,1)
        for j = 1:size(pic,2)
            if mean(pic(i,j,:)) < 250
                pic_data(n,1) = j;
                pic_data(n,2) = size(pic,1)-i;
                n = n+1;
            end
            
        end
    end
    y = pic_data;
%     % Visualize data against picture
%     figure(1)
%     clf
%     subplot(1,2,1);
%     image(pic);
%     subplot(1,2,2);
%     scatter(y(:,1),y(:,2),'.');

    for b = 1:4
        BE = all_betas(b);
        all_mix = {};
        for n = 1:10
            all_mix{n} = gmm_beta(y,n,0,BE);
            all_f(n) = all_mix{n}.fm;
        end
        
        idx = find(all_f==max(all_f));
        mix = all_mix{idx};
        win_mix{m,b} = all_mix{idx};
        for i=1:idx
            plot(mix.state(i).m(1),mix.state(i).m(2),'rx');
        end
        all_idx(m,b) = idx;
    
%         subplot(2,2,b);
%         hold on;
%         plot_mix(mix,[0 200 0 200],1,'r');
%         plot(y(:,1),y(:,2),'.');

        % Response model
        if all_idx(m,b) > 2
            is_image(m,b) = 1;
        else
            is_image(m,b) = 0;
        end
    end
end

save workspace_gmm_nov22


%% Pick out FE and response model

addpath(genpath('./tools'));
addpath(genpath('../../../tools'));

red = [0.6350    0.0780    0.1840];
blue = [0    0.4470    0.7410];
truth = [1 3 4 5 6 10 11 17 19 21 22 24];
no_pic = [2 7 8 9 12 13 14 15 16 18 20 23];
for t = 1:24
    if ismember(t,truth)
    pic_exists(t) = 1;
    else  
       pic_exists(t) = 0;
    end
end

for m = 1:24
    for b = 1:4
        win_f(m,b) = win_mix{m,b}.fm;
        if b>1
        mean_prec(m,b) = mean(win_mix{m,b}.lambda);
        end
    end
end

resp_model = win_f./all_idx;
for b = 1:4
    for m = 1:24
        
        if all_idx(m,b) >= 10
            is_image(m,b) = 1;
        else
            if rand >= 0.5
                is_image(m,b) = 1;
            else
                is_image(m,b) = 0;
            end
        end
    end
    decision = is_image(:,b)-pic_exists';
    
    false_positive(b)= sum(decision>0);
    false_negative(b) = sum(decision<0);
    
end
accuracy = [false_positive' false_negative'];
accuracy(1,1) =0.1;
accuracy(4,2) = 0.1;


%%
figure(201);
subplot(1,2,1);
bar(mean(all_idx,1), 'FaceColor','k');
set(gca,'XTickLabel',{'0.01', '0.5', '1', '100'})
xlabel('\beta value');
ylabel('No. features');

subplot(1,2,2);
b = bar(accuracy);
b(1).FaceColor = blue;
b(1).EdgeColor = 'none';
b(2).FaceColor = red;
b(2).EdgeColor = 'none';
set(gca,'XTickLabel',{'0.01', '0.5', '1', '100'})
xlabel('\beta value');
ylabel('Accuracy');
purty_plot(201,'../figures/BIpaper_Figure4', 'png');

