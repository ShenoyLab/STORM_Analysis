%% summary plot of reaction rate
% ----------------------------- SUMMARY ------------------------------
function results = get_summary_plot(data)

num_of_data = length(data);
names = string([]);
for cell_idx = 1:num_of_data
%     names = [names;extractBefore(data{cell_idx,1}.name,'_')];
    names = [names;data{cell_idx,1}.name(1:find(data{cell_idx,1}.name == '_', 1, 'last')-1)];
end
num_of_treatments = length(unique(names));
sample_info = zeros(num_of_treatments,1);
temp = unique(names);
for i = 1:num_of_treatments
    sample_info(i) = length(names(names==temp(i)));
end

% assign color
colorlib = string(['r';'b';'g';'c';'m';'k']);
cmaps = strings([num_of_data,1]);
treatment_index = zeros([num_of_data,1]);
treatments = unique(names);
for i = 1:num_of_treatments
    cmaps(names == treatments(i)) = colorlib(i);
    treatment_index(names == treatments(i)) = i;
end

figure('name','summary')
set(gcf, 'Position',  [50, 50, 500, 700]) % [leftbottom_x, left_bottom_y, L, Height]
subplot(421)
for cell_idx = 1:length(data)
    cmap = cmaps(cell_idx);
    spacing = data{cell_idx,1}.spacing; d = spacing/1000;
    [p,x_bins] = hist(d,50);
    plot(x_bins,smoothdata(p/sum(p)/(x_bins(2)-x_bins(1)),'gaussian',5),cmap,'linewidth',0.75); hold on
    xlim([0,1.1]);
    title('PDF for Spacing');xlabel('Spacing /\mum');
end

subplot(422)
R_store = [];
for cell_idx = 1:length(data)
    cmap = cmaps(cell_idx);
    radius = data{cell_idx,1}.hetero_radius; R = radius/1000;
    R_store = [R_store; mean(R)];
    [p,x_bins] = hist(R,length(R)/20);
    plot(x_bins,smoothdata(p/sum(p)/(x_bins(2)-x_bins(1)),'gaussian',1),cmap,'linewidth',0.75); hold on
    title('PDF for Radius');xlabel('Radius /\mum');
    xlim([0,0.15]);
end

subplot(423) % barchart for radius
radius_mean = [];radius_std = [];
for treatment_type = 1:num_of_treatments
    radius_mean = [radius_mean; mean(R_store(treatment_index == treatment_type))];
    radius_std = [radius_std; std(R_store(treatment_index == treatment_type))];
end
bar(radius_mean); hold on
errorbar(1:num_of_treatments,radius_mean,radius_std,radius_std,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 2)
title('Bar chart for radius');ylabel('Radius /\mum');
set(gca,'xticklabel',treatments)
% radius_collect = [];
% radius_index = [];
% for cell_idx = 1:length(data)
%     radius_collect = [radius_collect;data{cell_idx,1}.hetero_radius./1000];
%     radius_index = [radius_index;(treatment_index(cell_idx)-1)*ones(length(data{cell_idx,1}.hetero_radius),1)];
% end
% boxplot(radius_collect, radius_index,'Labels',treatments)
% set(gca,'ylim',[0,0.15])

subplot(424) % barchart for lads thickness
T_store = [];
for cell_idx = 1:length(data)
%     thickness = data{cell_idx,1}.lads_flat(:,2); T = thickness/1000;
    thickness = sum(data{cell_idx,1}.lads_area)/data{cell_idx,1}.lads_flat(end,1);T = thickness/1000;
    T_store = [T_store;mean(T)];
end
thickness_mean = [];thickness_std = [];
for treatment_type = 1:num_of_treatments
    thickness_mean = [thickness_mean; mean(T_store(treatment_index == treatment_type))];
    thickness_std = [thickness_std; std(T_store(treatment_index == treatment_type))];
end
% thickness_mean = [mean(T_store(1:sample_info(1))),mean(T_store(sample_info(1)+1:end))];
% thickness_std = [std(T_store(1:sample_info(1))),std(T_store(sample_info(1)+1:end))];
bar(thickness_mean); hold on
errorbar(1:num_of_treatments,thickness_mean,thickness_std,thickness_std,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 2)
title('Bar chart for thickness');ylabel('Thickness /\mum');
set(gca,'xticklabel',treatments)


subplot(425)
g_ac_store = [];
for cell_idx = 1:length(data)
    cmap = cmaps(cell_idx);
    mu = data{cell_idx,1}.g_ac_coeff(1);
    sigma = data{cell_idx,1}.g_ac_coeff(2);
    g_ac_store = [g_ac_store;exp(mu+sigma^2/2)];
    x_hat = linspace(0,80,1000);
    p_hat = exp(-(log(x_hat)-mu).^2/(2.*sigma.^2))./(x_hat.*sigma.*sqrt(2.*pi));
    plot(x_hat, p_hat,cmap,'linewidth',0.75);hold on
    xlim([0,40])
    title('PDF for \Gamma_{ac}');xlabel('\Gamma_{ac} /\mum');
end
subplot(426)
g_me_store = [];
for cell_idx = 1:length(data)
    cmap = cmaps(cell_idx);
    mu = data{cell_idx,1}.g_me_coeff(1);
    sigma = data{cell_idx,1}.g_me_coeff(2);
    g_me_store = [g_me_store;exp(mu+sigma^2/2)];
    x_hat = linspace(0,3,10000);
    p_hat = exp(-(log(x_hat)-mu).^2/(2.*sigma.^2))./(x_hat.*sigma.*sqrt(2.*pi));
    plot(x_hat, p_hat,cmap,'linewidth',0.75);hold on
    xlim([0,1]);
    ylim([0,6]);
    title('PDF for \Gamma_{me}');xlabel('\Gamma_{me} /\mum');
end
subplot(427)
ac_mean = [];ac_std = [];
for treatment_type = 1:num_of_treatments
    ac_mean = [ac_mean; mean(g_ac_store(treatment_index == treatment_type))];
    ac_std = [ac_std; std(g_ac_store(treatment_index == treatment_type))];
end
bar(ac_mean); hold on
errorbar(1:num_of_treatments,ac_mean,ac_std,ac_std,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 2)
title('Bar chart for \Gamma_{ac}');
set(gca,'xticklabel',treatments)

subplot(428)
me_mean = [];me_std = [];
for treatment_type = 1:num_of_treatments
    me_mean = [me_mean; mean(g_me_store(treatment_index == treatment_type))];
    me_std = [me_std; std(g_me_store(treatment_index == treatment_type))];
end
bar(me_mean); hold on
errorbar(1:num_of_treatments,me_mean,me_std,me_std,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 2)
title('Bar chart for \Gamma_{me}');
set(gca,'xticklabel',treatments)

results = cell(1);
results{1}.R_store = R_store;
results{1}.T_store = T_store;
results{1}.g_ac_store = g_ac_store;
results{1}.g_me_store = g_me_store;