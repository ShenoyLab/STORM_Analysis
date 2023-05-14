% draw_hq_boxplot3(data_mixed{3,1}.spacing,data_mixed{1,1}.spacing,data_mixed{2,1}.spacing,{'Young','Old','Diseased'},'Spacing')
% draw_hq_boxplot3(data_mixed{3,1}.affinity,data_mixed{1,1}.affinity,data_mixed{2,1}.affinity,{'Young','Old','Diseased'},'Affinity(Normalized)')
% draw_hq_boxplot3(data_mixed{3,1}.radius,data_mixed{1,1}.radius,data_mixed{2,1}.radius,{'Young','Old','Diseased'},'Radius nm')
% draw_hq_boxplot3(data_mixed{3,1}.lads_thick,data_mixed{1,1}.lads_thick,data_mixed{2,1}.lads_thick,{'Young','Old','Diseased'},'LADs Thickness /nm')
% draw_hq_boxplot3(data_mixed{3,1}.g_ratio,data_mixed{1,1}.g_ratio,data_mixed{2,1}.g_ratio,{'Young','Old','Diseased'},'\Gamma_{me}/\Gamma_{ac}')
% draw_hq_boxplot3(data_mixed{3,1}.g_ac,data_mixed{1,1}.g_ac,data_mixed{2,1}.g_ac,{'Young','Old','Diseased'},'\Gamma_{me}/\Gamma_{ac}')
% draw_hq_boxplot3(data_mixed{3,1}.g_me,data_mixed{1,1}.g_me,data_mixed{2,1}.g_me,{'Young','Old','Diseased'},'\Gamma_{me}')

function draw_hq_boxplot3(data_1,data_2,data_3,xname,yname)

figure('Position',  [500, 500, 300, 250],'name',yname);

boxplot([data_1;data_2;data_3], [0*ones(length(data_1),1);ones(length(data_2),1);2*ones(length(data_3),1)], ...
       'symbol','', ...
       'Labels',xname, ...
       'BoxStyle','outline')
hold on
plot([mean(data_1),mean(data_2),mean(data_3)], 'd')
hold off
set(gca,'ylim',[0,1.1*max(max(prctile(data_1,88),prctile(data_2,88)),prctile(data_3,88))])
h = findobj(gca,'Tag','Box');
colors = [241,168,59;234,64,37;70,143,244]/255;
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',0.4);
end
ylabel(yname)
set(findobj(gca,'type','line'),'linew',2,'color','k')
set(gca,'LineWidth',2)
set(gca,'FontName','Arial','FontSize',18)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'FontName','Arial','FontSize',18,'XTickLabelRotation',0)
% set(gca,'XTick',[1,1.5])

% add p value
% [~,p,~,~]=ttest2(data_1,data_2);
% text_height = max(prctile(data_1,92),prctile(data_2,92))*1.3;
% if p <= 0.001
%     p_symbol = '***';
% elseif p <= 0.01
%     p_symbol = '**';
% elseif p <=0.05
%     p_symbol = '*';
% else
%     p_symbol = 'ns';
%     text_height = max(max(prctile(data_1,91),prctile(data_2,91)),prctile(data_3,91))*1.35;
% end
% 
% yt = get(gca, 'YTick');
% xt = get(gca, 'XTick');hold on
% plot(xt([1 2]), [1 1]*max(max(prctile(data_1,91),prctile(data_2,91)),prctile(data_3,91))*1.25, '-k','LineWidth',1);hold on
% font_size = 25;
% text(mean(xt)-length(p_symbol)/2*font_size/400,text_height,p_symbol,'FontSize',font_size)
end