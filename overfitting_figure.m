addpath('./tools')
nice_colors


% Training data
x = [0.1:0.1:3]';
y = x.^2 + normrnd(0,0.5,length(x),1);

% Test data
x_test = [1:0.03:1.88]';
y_test = x_test.^2 + normrnd(0,0.5,length(x_test),1);

% Linear fit
regs = regstats(y,x);
trend_line = regs.beta(1)+ regs.beta(2)*x;
trend_line_test = regs.beta(1) + regs.beta(2)*x_test;

% Quadratif fit
qf = fit(x,y,'poly2');
quad_fit = qf.p1*x.^2 + qf.p2*x + qf.p3;
quad_fit_test = qf.p1*x_test.^2 + qf.p2*x_test + qf.p3;

% Spline interpolant
[sf,gof,out] = fit(x,y,'smoothingspline');
sf_test = feval(sf,x_test);

% Within-sample error
sse_within(1) = sum((trend_line-y).^2);
sse_within(2) = sum((quad_fit-y).^2);
sse_within(3) = gof.sse;

% Out-of-sample error
sse_oos(1) = sum((trend_line_test-y_test).^2);
sse_oos(2) = sum((quad_fit_test-y_test).^2);
sse_oos(3) = sum((sf_test-y_test).^2);

%%
marker_size = 10;
line_width = 2;
font_size = 15;

figure(1);
clf
subplot(2,6,1:2)
plot(x,y,'.','Color',dark_grey,'MarkerSize',marker_size);
hold on;
plot(x,trend_line,'k-','LineWidth',line_width);
xticks([]);
yticks([]);
set(gca,'FontSize',font_size);

subplot(2,6,3:4);
plot(x,y,'.','Color',dark_grey,'MarkerSize',marker_size);
hold on
plot(x,quad_fit,'k-','LineWidth',line_width);

xticks([]);
yticks([]);
set(gca,'FontSize',font_size);

subplot(2,6,5:6);
h = plot(sf,x,y);
h(1).Color = dark_grey;
h(1).MarkerSize = marker_size; 
h(2).Color = 'k';
h(2).LineWidth = line_width; 
legend off
xticks([]);
yticks([]);
set(gca,'FontSize',font_size);

subplot(2,6,7:9)
bar(sse_within,'FaceColor',dark_grey,'EdgeColor','None');
title('Within-sample SSE')
set(gca,'FontSize',font_size);
subplot(2,6,10:12)
bar(sse_oos,'FaceColor',dark_grey,'EdgeColor','None');
title('Out-of-sample SSE')
set(gca,'FontSize',font_size);

saveas(gcf,'../figures/Overfit','epsc')
