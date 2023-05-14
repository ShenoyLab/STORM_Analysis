% plot reaction rates
% draw_hq_boxplot(g_me_stiff,g_me_soft,{'Stiff','Soft'},'\Gamma_{me}/\Gamma_{ac}')
function draw_hq_pdf(data_1,data_2,xname,yname)
figure('Position',  [500, 500, 300, 250]);
[p,x_bins] = hist(data_1,50);
plot(x_bins,smoothdata(p/sum(p)/(x_bins(2)-x_bins(1)),'gaussian',5),"Color",[234,64,37]/255,'linewidth',3);hold on
[p,x_bins] = hist(data_2,50);
plot(x_bins,smoothdata(p/sum(p)/(x_bins(2)-x_bins(1)),'gaussian',5),"Color",[241,168,59]/255,'linewidth',3);
grid on
hold on
% set(gca,'ylim',[0,1.2*max(prctile(data_1,91),prctile(data_2,91))])
xlabel(xname)
ylabel(yname)
% set(findobj(gca,'type','line'),'linew',2,'color','k')
set(gca,'LineWidth',2)
set(gca,'FontName','Arial','FontSize',18)

% set(gca,'XTick',[1,1.5])
