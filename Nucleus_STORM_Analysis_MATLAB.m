% This code works as a standard code for the params extraction of 
% STORM images (@ShenoyLab)
% Notes: This file is used to extract key params of nucleus as follows:
% 1. Inner Heterochromatin Size 
% 2. LADs thickness


clc, clear, close all


% file directory
% all the H2B localization file should be stored in this directory
myDir = 'Input_LocsLib'; 
addpath(genpath(myDir))
addpath(genpath('FuncLib'))

eps = 30;
min_num = 3;

myFiles = dir(fullfile(myDir, '*.txt')); % gets all txt files in struct
data = cell(length(myFiles), 1);
for cell_idx = 1:length(myFiles)
    data{cell_idx}.name = extractBefore(myFiles(cell_idx).name,'.txt');
    locs = importdata(myFiles(cell_idx).name);
    locs = locs(:,1:2);
    locs = unique(locs,'rows');
    fprintf('Now processing %s -- Data Size: %d \n', ...
        data{cell_idx}.name, length(locs(:,1)))

    % downsampling data -- default downsample size (dss) = 4
    down_sample_scale = 4;
    down_sample_size = fix(length(locs(:,1))/down_sample_scale);
    vec = 1:length(locs(:,1));
    rng('default');
    vec = vec(randperm(length(vec)));
    I = vec(1:down_sample_size);
    locs = locs(I,:);
    clear down_sample_size vec I
    
    % get threshold (for the demo we choose 50 percentile for illustration)
    [density_threshold, storm_data] = get_percentile_threshold(locs,50);

    % apply threshold (obtain heterochromatic localization)
    Img = storm_data(:,1:2);
    density = storm_data(:,3);
    Hetero = Img(density>=density_threshold,:);

    % find nuclear boundary
    bd_old = Img(boundary(Img, 0.5),:);
    bd = smoothdata(bd_old(1:length(bd_old)-1,:),'gaussian',10); clear bd_old
    bd(end+1,:) = bd(1,:);

    % approximate size of nucleus
    nucleus_area = polyarea(bd(:,1),bd(:,2));
    approximate_radius = sqrt(nucleus_area/pi);
    fprintf('Nucleus Diameter is %.0f nm \n',approximate_radius*2)

    % dbscan
    fprintf('Now Running dbscan ... (might take a few minutes)')
    labels = dbscan(Hetero,eps,min_num);
    fprintf('  -- Done! \n')

    Hetero_flt = removerows(Hetero,'ind',find(labels == -1)); % filter background noise
    labels_flt = removerows(labels,'ind',find(labels == -1)); 
    labels_flt_shuffle = shufflelabel(labels_flt); % shuffle label

    data{cell_idx}.Hetero_flt_with_label = [Hetero_flt,labels_flt];
    numGroups = length(unique(labels_flt_shuffle));
    fprintf('Number of connected clusters %d .\n',numGroups);

    % define LADs and non-LADs domain
    lads_is = [];
    lads_index = [];
    non_lads_is = [];
    non_lads_index = [];
    f = waitbar(0, 'Starting');
    cnt = 1;
    nl_cnt = 1;
    for i=1:numGroups
        grp = Hetero_flt(labels_flt == i,:);
        % sampling from grp (speed up calculation)
        mid_point = [mean(grp(:,1)), mean(grp(:,2))];
        rid = randi([1 length(grp(:,1))],1,ceil(length(grp)/15));
        rnd_point = grp(rid,:);
        sam_points = [mid_point; rnd_point];
        dist_mat = pdist2(sam_points, bd);
        min_dist = min(min(dist_mat));
        if min_dist <= 0.05*approximate_radius && length(grp(:,1))>=35
            lads_is = [lads_is; find(labels_flt == i)]; % old index
            lads_index = [lads_index; cnt*ones(length(grp(:,1)),1)]; % new index
            cnt = cnt + 1;
         
        % filter away noise point    
        elseif length(grp(:,1))>=35
            temp_bd = grp(boundary(grp,0.3),:);
            temp_area = polyarea(temp_bd(:,1),temp_bd(:,2));
            temp_r = sqrt(temp_area/pi);
            non_lads_is = [non_lads_is; find(labels_flt == i)];
            non_lads_index = [non_lads_index; nl_cnt*ones(length(grp(:,1)),1)];
            nl_cnt = nl_cnt + 1;
            clear temp_r temp_area temp_bd
        end
    
        waitbar(i/numGroups, f, sprintf('Finding LADs regions: %d %% ...', floor(i/numGroups*100)));
        pause(0.005);
    end
    close(f)
    non_lads = Hetero_flt(non_lads_is,:);
    lads = Hetero_flt(lads_is,:);
    
    data{cell_idx,1}.lads = [lads,lads_index];
    data{cell_idx,1}.non_lads = [non_lads,non_lads_index];
    data{cell_idx,1}.lads2total = length(lads)/(length(lads)+length(non_lads));
    
    % analysis of inner heterochromatin domain size
    num_of_hetero = length(unique(non_lads_index));
    hetero_radius = zeros(num_of_hetero,1);
    hetero_center = zeros(num_of_hetero,2);
    for i=1:num_of_hetero
        grp = non_lads(non_lads_index==i,:);
        hetero_center(i,:) = [mean(grp(:,1)),mean(grp(:,2))];
        bd_i = grp(boundary(grp,0.3),:);
        area = polyarea(bd_i(:,1),bd_i(:,2));
        hetero_radius(i) = sqrt(area/pi);
        clear area grp bd_i
    end
    ave_radius = mean(hetero_radius);
    data{cell_idx,1}.hetero_center = hetero_center;
    data{cell_idx,1}.hetero_radius = hetero_radius;
    fprintf('Average Diameter is: %.2f nm; \n', ave_radius*2)

    % find spacing
    N_neighbour = 5;
    dist_mat = mink(pdist2(hetero_center,hetero_center),N_neighbour+1,2);
    spacing = sum(dist_mat,2)/N_neighbour;
    data{cell_idx,1}.spacing = spacing;
    
    %% compute lads thickness

    % create accumulate center location
    bd_acc = zeros(length(bd(:,1)),1);
    for i=2:length(bd(:,1))
        bd_acc(i) = bd_acc(i-1) + norm(bd(i,:)-bd(i-1,:));
    end

    % create unit vector
    tangents=diff(bd);
    u_tangents=tangents./sqrt(sum(tangents.^2,2));
    normals = zeros(length(bd(:,1))-1,2);
    for i=1:length(bd(:,1))-1
        normals(i,:)=[0 -1;1 0]*u_tangents(i,:)';
    end
    normals(end+1,:) = normals(1,:);

    % define LADs segments
    seg_num = 50;
    seg_len = floor(length(bd(:,1))/seg_num);

    % calculate lads thickness within each segment
    bd_idx = 1;
    indent_len = 2000;
    seg_len_store = [];
    seg_area = [];
    seg_thickness = [];
    seg_starts = [];
    seg_normals = [];
    if mod(length(bd(:,1)),seg_num)==0
        num_iters = seg_num;
    else
        num_iters = seg_num+1;
    end
    for seg_idx = 1:num_iters
        point1 = bd(bd_idx,:);
        point4 = point1 + indent_len*normals(bd_idx,:);
        seg_normals = [seg_normals;normals(bd_idx,:)];
        if seg_idx<num_iters
            bd_idx = bd_idx +seg_len;
        else
            bd_idx = bd_idx + length(bd(:,1))-seg_num*seg_len-1;
        end
        point2 = bd(bd_idx,:);
        point3 = point2 + indent_len*normals(bd_idx,:);
        
        in = inpolygon(lads(:,1),lads(:,2),...
            [point1(1);point2(1);point3(1);point4(1);point1(1)],...
            [point1(2);point2(2);point3(2);point4(2);point1(2)]);
        grp = lads(in,:);
        bd_grp = grp(boundary(grp,0.5),:);
        area = polyarea(bd_grp(:,1),bd_grp(:,2));
        seg_area = [seg_area;area];
        seg_len_store = [seg_len_store;norm(point1-point2)];
        seg_thickness = [seg_thickness;area/norm(point1-point2)];
        seg_starts = [seg_starts;point1];
        
    end
    seg_starts(end+1,:) = seg_starts(1,:);
    seg_normals(end+1,:) = seg_normals(1,:);
    data{cell_idx,1}.seg_normals = seg_normals;
    data{cell_idx,1}.seg_locs = seg_starts;
    data{cell_idx,1}.lads_seg_len = seg_len_store;
    seg_thickness(seg_thickness <= 0) = 1;
    data{cell_idx,1}.lads_seg_thickness = seg_thickness;
    fprintf('Average LADs thickness is: %.2f nm; \n\n', sum(seg_area)/sum(seg_len_store))


    figure();
    plot(bd(:,1),bd(:,2),'r-','LineWidth',2); hold on
    scatter(lads(:,1),lads(:,2),0.5,'g'); hold on
    scatter(non_lads(:,1),non_lads(:,2),0.5,'k'); hold on
    legend('off')
    for i = 1:length(seg_starts(:,1))-1
        point1 = seg_starts(i,:);
        point2 = seg_starts(i+1,:);
        point4 = point1 + seg_thickness(i)*seg_normals(i,:);
        point3 = point2 + seg_thickness(i)*seg_normals(i+1,:);
        pgon = polyshape([point1(1);point2(1);point3(1);point4(1);point1(1)],...
            [point1(2);point2(2);point3(2);point4(2);point1(2)]);
        plot(pgon,'FaceColor','k');hold on
        axis equal
    end
    drawnow()

    % visualize Voronoi density plot
    %% calculate Voronoi density
    % f = waitbar(0,'Calculating Voronoi Areas');
    x = locs(:,1); y = locs(:,2);
    dt = delaunayTriangulation(locs(:,1),locs(:,2));
    [vertices,connections] = voronoiDiagram(dt);
    voronoi_cells = cellfun(@(x) vertices(x,:),connections,'UniformOutput',false);
    voronoi_areas = cellfun(@(x) polyarea(x(:,1),x(:,2)),voronoi_cells,'UniformOutput',false);
    voronoi_areas = vertcat(voronoi_areas{:});
    Density = 1./voronoi_areas;
    
    
    %% visualize the nuclei in voronoi polygon
    logDensity = log10(Density);
    
%     min_val = prctile(logDensity,10);

    max_val = -1.9;
    logDensity(logDensity>=max_val)=max_val;
    min_val = -2.8;
%     max_val = prctile(logDensity,65);
    idx = logDensity>=min_val & logDensity<=max_val;    
    connections = connections(idx);    
    logDensity = logDensity(idx);
    
    max_connections = max(cellfun(@(x) length(x), connections));
    faces = NaN(length(connections),max_connections);
    for i = 1:length(connections)
        faces(i,1:length(connections{i})) = connections{i};
    end
    
    monitor = figure();
    set(gcf,'name','Voronoi Area Density Plot','NumberTitle','off','color','k','units','normalized','position',[0.25 0.15 0.5 0.7]);
    
    subplot = @(m,n,p) subtightplot (m, n, p, [0 0], [0 0], [0 0]);    
    subplot(1,1,1)
    hold on
    patch('Vertices',vertices,'Faces',faces,'FaceVertexCData',logDensity,'FaceColor','flat','edgecolor','none')    
    h = colorbar();    
    pbaspect([1 1 1])
    h.Position = [0.9 0.1 0.02 0.8];
    h.TickLabelInterpreter = 'latex';
    h.FontSize = 14;
    h.Color = 'w';
    set(gca,'colorscale','log')
    ylabel(h, '$Log10(Voronoi Area-nm^{2})$','FontSize',14,'interpreter','latex');
    caxis([min_val max_val])
    pixels_um_x = 3*1000;
    pixels_um_y = 0.3*1000;
    rectangle('Position',[min(x)-900 min(y)-900 pixels_um_x pixels_um_y],'facecolor','w')
    map = colormap(jet);
    
    axis equal
    L=30000;
    % xlim([min(x)-1000 max(x)+1000])
    % ylim([min(y)-1000 max(y)+1000])
    xlim([min(x)-1000 min(x)+L+1000])
    ylim([min(y)-1000 min(y)+L+1000])
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'colormap',map,'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
    
%     name = extractBefore(myFiles(cell_idx).name,'.');
%     print(gcf,['InputLocsLib/Voronoi_density_',name,'.png'],'-dpng','-r800');
    % close(monitor)
end
    