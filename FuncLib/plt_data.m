% clc,clear,close all
figure()
% load dt.mat;
height_glass = [];
height_stiff = [];
height_soft = [];
affinity_glass = [];
affinity_stiff = [];
affinity_soft = [];

filter = 20;
for i = 1:14
    x = data{i,1}.lads_flat(:,1);
    y = data{i,1}.lads_flat(:,2);
    lads_thick = data{i,1}.lads_thickness;
    affinity = data{i,1}.affinity;
    yy = smoothdata(y,'gaussian',20);
    subplot(14,1,i)
    bar(x,yy,1,'FaceColor',[0.2 0.2 0.2]);hold on
    
%     xlabel('Circumference / nm')
%     ylabel('Height / nm')
    xlim([0,max(x)])
    ylim([0,1500])
    legend('off')
    if i <= 5
        
        height_glass = [height_glass;lads_thick];
        affinity_glass = [affinity_glass;affinity];
    elseif i <=9
        
        height_stiff = [height_stiff;lads_thick];
        affinity_stiff = [affinity_stiff;affinity];
        
    else
        height_soft = [height_soft;lads_thick];
        affinity_soft = [affinity_soft;affinity];
    end
end

height_mean = [mean(height_glass),mean(height_stiff),mean(height_soft)];
height_std = [std(height_glass),std(height_stiff),std(height_soft)];
neg = height_std;
pos = height_std;
figure;
bar(1:3,height_mean); hold on
errorbar(1:3,height_mean,neg,pos,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 2)

aff_mean = [mean(affinity_glass),mean(affinity_stiff),mean(affinity_soft)];
aff_std = [std(affinity_glass),std(affinity_stiff),std(affinity_soft)];
neg = aff_std;
pos = aff_std;
figure;
bar(1:3,aff_mean); hold on
errorbar(1:3,aff_mean,neg,pos,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 2)