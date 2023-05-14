function data_mixed = get_mixed_data(data)

num_of_data = length(data);
names = string([]);
for cell_idx = 1:num_of_data
%     names = [names;extractBefore(data{cell_idx,1}.name,'_')];
    names = [names;data{cell_idx,1}.name(1:find(data{cell_idx,1}.name == '_', 1, 'last')-1)];
end
num_of_treatments = length(unique(names));

% get number of each treatment
sample_info = zeros(num_of_treatments,1);
temp = unique(names);
for i = 1:num_of_treatments
    sample_info(i) = length(names(names==temp(i)));
end

treatment_index = zeros([num_of_data,1]);
treatments = unique(names);
for i = 1:num_of_treatments
    
    treatment_index(names == treatments(i)) = i;
end

% pre-define color map
cmaps = [239,27,35;246,130,32;255,203,9;137,154,157];

data_mixed = cell(num_of_treatments,1);
cell_idx = 1;
% now mix all the data
for i = 1:num_of_treatments
    data_mixed{i,1}.name = treatments(i);
    data_mixed{i,1}.radius = [];
    data_mixed{i,1}.spacing = [];
    data_mixed{i,1}.lads_thick = [];
    data_mixed{i,1}.lads_ave = [];
    data_mixed{i,1}.affinity = [];
    for j = 1:length(find(treatment_index==i))
        data_mixed{i,1}.radius = [data_mixed{i,1}.radius;data{cell_idx,1}.hetero_radius];
        data_mixed{i,1}.spacing = [data_mixed{i,1}.spacing;data{cell_idx,1}.spacing];
%         temp1 = prctile(data{cell_idx,1}.lads_thickness,5);
%         temp2 = prctile(data{cell_idx,1}.lads_thickness,95);
%         data_mixed{i,1}.lads_thick = [data_mixed{i,1}.lads_thick;data{cell_idx,1}.lads_thickness(data{cell_idx,1}.lads_thickness>temp1 & data{cell_idx,1}.lads_thickness<temp2)];
        temp = smoothdata(data{cell_idx,1}.lads_flat(:,2),'gaussian',10);
        data_mixed{i,1}.lads_thick = [data_mixed{i,1}.lads_thick;temp(temp>1)];
        data_mixed{i,1}.affinity = [data_mixed{i,1}.affinity;data{cell_idx,1}.affinity];
%         data_mixed{i,1}.lads_ave = [data_mixed{i,1}.lads_ave; data{cell_idx,1}.ave_thickness];
        
        cell_idx = cell_idx + 1;
    end

    data_mixed{i,1}.radius_mean = mean(data_mixed{i,1}.radius);
    data_mixed{i,1}.spacing_mean = mean(data_mixed{i,1}.spacing);
    data_mixed{i,1}.lads_thick_mean = mean(data_mixed{i,1}.lads_thick);
    data_mixed{i,1}.affinity_mean = mean(data_mixed{i,1}.affinity);

end

% fitting mixture data
D = 20;
for i = 1:num_of_treatments
    g_ac = 10;
    R = data_mixed{i,1}.radius/1000;
    R = R(R<0.25);
    data_mixed{i,1}.g_me = g_ac.^2.*R.^2./(D-g_ac.*R.^2);
    data_mixed{i,1}.g_me_coeff = logNfit(data_mixed{i,1}.g_me(data_mixed{i,1}.g_me<0.08*(max(data_mixed{i,1}.g_me))), 50, [-5,3]);
    data_mixed{i,1}.g_me_mean = exp(data_mixed{i,1}.g_me_coeff(1)+(data_mixed{i,1}.g_me_coeff(2))^2/2);

    
%     % g_ac
%     data_mixed{i,1}.g_ac = D./((data_mixed{i,1}.spacing/1000).^2); 
%     data_mixed{i,1}.g_ac_coeff = logNfit(data_mixed{i,1}.g_ac(data_mixed{i,1}.g_ac<0.4*max(data_mixed{i,1}.g_ac)), 20, [3,0.2]);
%     data_mixed{i,1}.g_ac_mean = exp(data_mixed{i,1}.g_ac_coeff(1)+data_mixed{i,1}.g_ac_coeff(2)^2/2); 
%     % g_me
%     data_mixed{i,1}.g_me = data_mixed{i,1}.g_ac./(((data_mixed{i,1}.spacing./1000)./(data_mixed{i,1}.radius/1000)).^2-1);
% %     g_me = g_ac./(((spacing/1000)./hetero_radius).^2-1);
%     data_mixed{i,1}.g_ratio = data_mixed{i,1}.g_me./data_mixed{i,1}.g_ac;
%     data_mixed{i,1}.g_ratio = data_mixed{i,1}.g_ratio(data_mixed{i,1}.g_ratio>0);
%     data_mixed{i,1}.g_ratio_mean = mean(data_mixed{i,1}.g_ratio);
%     data_mixed{i,1}.g_me = data_mixed{i,1}.g_me(data_mixed{i,1}.g_me>0);
%     data_mixed{i,1}.g_me = data_mixed{i,1}.g_me(data_mixed{i,1}.g_me<=10);
%     data_mixed{i,1}.g_me_coeff = logNfit(data_mixed{i,1}.g_me(data_mixed{i,1}.g_me<0.08*(max(data_mixed{i,1}.g_me))), 15, [-5,3]); 
%     data_mixed{i,1}.g_me_mean = exp(data_mixed{i,1}.g_me_coeff(1)+data_mixed{i,1}.g_me_coeff(2)^2/2); 
end

% plotting mixture data
figure('name','Data Summary')
set(gcf, 'Position',  [50, 50, 900, 900]) % [leftbottom_x, left_bottom_y, L, Height]

% 1. PDF for radius
% subplot(331)
% figure()
% for i = 1:num_of_treatments
%     cmap = cmaps(i,:);
%     [p,x_bins] = hist(data_mixed{i,1}.radius,50);
%     plot(x_bins,smoothdata(p/sum(p)/(x_bins(2)-x_bins(1)),'gaussian',5),"Color",cmap/255,'linewidth',3);
%     grid on
%     hold on
%     xlabel('Radius /nm')
%     ylabel('PDF')
%     set(gca,'LineWidth',2)
%     set(gca,'FontName','Arial','FontSize',18)
% end
% % 2. PDF for spacing
% % subplot(332)
% figure()
% for i = 1:num_of_treatments
%     cmap = cmaps(i,:);
%     [p,x_bins] = hist(data_mixed{i,1}.spacing,50);
%     plot(x_bins,smoothdata(p/sum(p)/(x_bins(2)-x_bins(1)),'gaussian',5),"Color",cmap/255,'linewidth',3);
%     grid on
%     hold on
%     xlabel('Radius /nm')
%     ylabel('PDF')
%     set(gca,'LineWidth',2)
%     set(gca,'FontName','Arial','FontSize',18)
% end

% 3. boxplot for lads thick



