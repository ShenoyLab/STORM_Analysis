function [threshold, storm_data] = get_percentile_threshold(locs,pct)
    x = locs(:,1);
    y = locs(:,2);
    
    dt = delaunayTriangulation(x,y);
    [vertices,connections] = voronoiDiagram(dt);
    voronoi_cells = cellfun(@(x) vertices(x,:),connections,'UniformOutput',false);
    voronoi_areas = cellfun(@(x) polyarea(x(:,1),x(:,2)),voronoi_cells,'UniformOutput',false);
    voronoi_areas = vertcat(voronoi_areas{:});
    
    % cal density
    Density = 1./voronoi_areas;
    storm_data = [locs(:,1:2),Density];
    threshold = prctile(storm_data(:,3),pct);
%     fprintf('Now processing %s -- Data Size: %d \n', data{cell_idx}.name, length(locs(:,1)))
    fprintf('Density threshold is %d\n', threshold)
end