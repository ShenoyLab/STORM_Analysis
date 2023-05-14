% This is the code generating high-resolution images for STORM-Voronoi
% Polygon Density

% store different value for uplim
% Stiffness 5.5e-2
% TSA 3e-2
% tenocytes 
% Fb
% H3K27me3 

clc,clear,close all
addpath(genpath('LocsLib'))
myDir = 'LocsLib/Dynamics'; % gets directory
myFiles = dir(fullfile(myDir,'*.txt')); %gets all wav files in struct
data = cell(length(myFiles),1);
for num_of_files = 1:1%length(myFiles)
    name = extractBefore(myFiles(num_of_files).name,'.'); % bug (might extract the file with the same name)
    Img = importdata(myFiles(num_of_files).name);
    get_density_plot(Img,2.51e-2,name);
end

% main function 
function get_density_plot(Img,uplim,name)

x = Img(:,1);
y = Img(:,2);

%% create Voronoi density plot
% canculate vol
f = waitbar(0,'Calculating Voronoi Areas');
dt = delaunayTriangulation(x,y);
[vertices,connections] = voronoiDiagram(dt);
voronoi_cells = cellfun(@(x) vertices(x,:),connections,'UniformOutput',false);
voronoi_areas = cellfun(@(x) polyarea(x(:,1),x(:,2)),voronoi_cells,'UniformOutput',false);
voronoi_areas = vertcat(voronoi_areas{:});

% define limit of colorbar
waitbar(0.5,f,'Calculating Voronoi Density');
% uplim = 5.5e-2;
L=30000; % 30000 for stiffness; 20000 for tsa; 30000 for Teno;

% apply the limit
Density = 1./voronoi_areas;
% Img = Img(Density>1.5e-2,:);
% writematrix([x,y,Density],['ImgLib/',name,'.txt'],'Delimiter',' ');
% max(Density)
% lolim = 10^(-3);
Density(Density>=uplim) = uplim;
% Density(Density<=lolim) = lolim;

logDensity = log10(Density);

min_val = prctile(logDensity,5);
max_val = prctile(logDensity,99);
min_val = -2.8; % 3e-2 -2.8
idx = logDensity>=min_val & logDensity<=max_val;    
connections = connections(idx);    
logDensity = logDensity(idx);

max_connections = max(cellfun(@(x) length(x), connections));
faces = NaN(length(connections),max_connections);
for i = 1:length(connections)
    faces(i,1:length(connections{i})) = connections{i};
end

waitbar(0.8,f,'Plotting Voronoi Polygons');
monitor = figure();
set(gcf,'name','Voronoi Area Density Plot','NumberTitle','off','color','k','units','normalized','position',[0.25 0.15 0.5 0.7]);

subplot = @(m,n,p) subtightplot (m, n, p, [0 0], [0 0], [0 0]);    
subplot(1,1,1)

% uimenu('Text','Save Image High Quality','ForegroundColor','b','CallBack',@save_image);
% uimenu('Text','Change Colormap Limits','ForegroundColor','b','CallBack',@change_colormap_limits);    

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
% xlim([min(x)-1000 max(x)+1000])
% ylim([min(y)-1000 max(y)+1000])
xlim([min(x)-1000 min(x)+L+1000])
ylim([min(y)-1000 min(y)+L+1000])
set(gcf, 'InvertHardCopy', 'off');
set(gca,'colormap',map,'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
waitbar(1,f,'Plotting Voronoi Polygons');
close(f)

% saveas(gcf,['ImgLib/',name,'.png'])
print(gcf,['ImgLib/ActD/Voronoi_density_',name,'.png'],'-dpng','-r800');
close(monitor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = waitbar(0,'Plotting STORM Images');
monitor2 = figure();
set(gcf,'name','Storm Plot','NumberTitle','off','color','k','units','normalized','position',[0.25 0.15 0.5 0.7]);
subplot = @(m,n,p) subtightplot (m, n, p, [0 0], [0 0], [0 0]);    
subplot(1,1,1)
% dss
down_sample_size = fix(length(x(:,1))/3);
vec = 1:length(x(:,1));
rng('default');
vec = vec(randperm(length(vec)));
I = vec(1:down_sample_size);
xx = x(I,:);
yy = y(I,:);
scatter(xx,yy,0.01,'MarkerEdgeColor','g','MarkerFaceColor','g'); hold on

pixels_um_x = 3*1000;
pixels_um_y = 0.3*1000;
rectangle('Position',[min(x)-900 min(y)-900 pixels_um_x pixels_um_y],'facecolor','w'); hold on
axis equal

% xlim([min(x)-1000 max(x)+1000])
% ylim([min(y)-1000 max(y)+1000])
xlim([min(x)-1000 min(x)+L+1000])
ylim([min(y)-1000 min(y)+L+1000])
set(gcf, 'InvertHardCopy', 'off');
set(gca,'color','k','box','on','BoxStyle','full','XColor','k','YColor','k');
waitbar(1,f,'Plotting Voronoi Polygons');
close(f)
print(gcf,['ImgLib/ActD/',['STORM_',name],'.png'],'-dpng','-r800');
pause(2)
close(monitor2)
% print(gcf,[name,'.tif'],'-dpng','-r300');
end