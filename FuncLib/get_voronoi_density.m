% get voronoi density
function storm_data = get_voronoi_density(locs)
x = locs(:,1);
y = locs(:,2);
f = waitbar(0,'Creating Voronoi Polygons');
dt = delaunayTriangulation(x,y);
[vertices,connections] = voronoiDiagram(dt);
voronoi_cells = cellfun(@(x) vertices(x,:),connections,'UniformOutput',false);
voronoi_areas = cellfun(@(x) polyarea(x(:,1),x(:,2)),voronoi_cells,'UniformOutput',false);
voronoi_areas = vertcat(voronoi_areas{:});

waitbar(0.5,f,'Calculating Voronoi Density');

% cal density
Density = 1./voronoi_areas;
storm_data = [locs,Density];
waitbar(1,f,'Done');
close(f)
end