function create_reaction_map(data)

for cell_idx = 1:length(data)
    D = 1;
    radius = data{cell_idx,1}.hetero_radius;
    R = radius/1000;
    spacing = data{cell_idx,1}.spacing; 
    d = spacing/1000;
    g_ac = D./(d.^2);
    g_me = g_ac./((d./R).^2-1);
    
    % load location info
    non_lads = data{cell_idx,1}.non_lads(:,1:2);
    labels = data{cell_idx,1}.non_lads(:,3);
    num_of_groups = length(unique(labels));
    
    figure()
    subplot(121) % g_ac
    
    ub = 30;    
    g_ac_round = round(g_ac);
    g_ac_round(g_ac_round >=ub ) = ub;
    clr = jet(ub);
    clr_bar = zeros(length(g_ac_round),3);
    for i = 1:num_of_groups
        clr_bar(i,:) = clr(g_ac_round(i),:);
    end
    
    for i = 1:num_of_groups
        grp = non_lads(labels==i,:);
        bd_i = grp(boundary(grp,0.5),:);
%         plot(bd_i(:,1),bd_i(:,2),'k','linewidth',0.5); hold on
        scatter(grp(:,1),grp(:,2),2,clr_bar(i,:),'filled'); hold on
    end
    
    set(gca,'Color','k')
    title('Distribution of \Gamma_{ac}')
    axis('equal')
    
    subplot(122) % g_me
    
    ub = 30;  
    scale_factor = 10;
    g_me_round = round(scale_factor*g_me);
    g_me_round(g_me_round >=ub ) = ub;
    g_me_round(g_me_round <=0 ) = 1;
    clr = jet(ub);
    clr_bar = zeros(length(g_me_round),3);
    for i = 1:num_of_groups
        clr_bar(i,:) = clr(g_me_round(i),:);
    end
    
    for i = 1:num_of_groups
        grp = non_lads(labels==i,:);
        bd_i = grp(boundary(grp,0.5),:);
%         plot(bd_i(:,1),bd_i(:,2),'k','linewidth',0.5); hold on
        scatter(grp(:,1),grp(:,2),2,clr_bar(i,:),'filled'); hold on
    end
    
    set(gca,'Color','k')
    title('Distribution of \Gamma_{me}')
    axis('equal')
    drawnow()
end