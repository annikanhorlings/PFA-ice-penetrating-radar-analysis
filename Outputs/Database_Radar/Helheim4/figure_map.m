%This code plots a map showing aquifer location
%Annika Horlings
%University of Washington
%15 February 2022

%% Get the data
%DEM
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Data/DEM/gimpDEM');
g = geotiffinfo('gimpdem_90m_v01.1.tif');
[A, R] = geotiffread('gimpdem_90m_v01.1.tif');
A = double(A);
dem = A;
Xmin = g.CornerCoords.X(1);
Xmax = g.CornerCoords.X(2);
Ymax = g.CornerCoords.Y(1);
Ymin = g.CornerCoords.Y(3);
xgrid = (Xmin:g.Width:Xmax);
ygrid = (Ymin:g.Height:Ymax);

%Clem's AR data
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Data/from_Clem/MATfiles');
date = struct2cell(load('date_by_year.mat'));
depth = struct2cell(load('depth_by_year.mat'));
elevation = struct2cell(load('elev_by_year.mat'));
latitude = struct2cell(load('lat_by_year.mat'));
longitude = struct2cell(load('long_by_year.mat'));

%convert lat/long to x and y coords polar stereographic, using polartereo_fwd, a function created Andy Bliss
for i = 2:length(date)
    [xclemall{i}, yclemall{i}] = polarstereo_fwd(latitude{i}, longitude{i}, 6378137.0,0.08181919,70,-45);
end

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/Helheim1');
d1r = load('database_Helheim1_RDS.mat'); extent1r = d1r.extentPFA_fromline; dates1r = d1r.dates;
xref1 = d1r.xcoordref; yref1 = d1r.ycoordref;
Cmnx1 = d1r.Cmnx; Cmxx1 = d1r.Cmxx; Cmny1 = d1r.Cmny; Cmxy1 = d1r.Cmxy;
xpfa1 = d1r.xPFA_inrange; ypfa1 = d1r.yPFA_inrange;

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/IkertivaqN1');
d2r = load('database_IkertivaqN1_RDS.mat'); extent2r = d2r.extentPFA_fromline; dates2r = d2r.dates;
xref2 = d2r.xcoordref; yref2 = d2r.ycoordref;
Cmnx2 = d2r.Cmnx; Cmxx2 = d2r.Cmxx; Cmny2 = d2r.Cmny; Cmxy2 = d2r.Cmxy;
xpfa2 = d2r.xPFA_inrange; ypfa2 = d2r.yPFA_inrange;

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/KogeBugtS1');
d3r = load('database_KogeBugtS1_RDS.mat'); extent3r = d3r.extentPFA_fromline; dates3r = d3r.dates;
xref3 = d3r.xcoordref; yref3 = d3r.ycoordref;
Cmnx3 = d3r.Cmnx; Cmxx3 = d3r.Cmxx; Cmny3 = d3r.Cmny; Cmxy3 = d3r.Cmxy;
xpfa3 = d3r.xPFA_inrange; ypfa3 = d3r.yPFA_inrange;

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/Helheim4');
d4r = load('database_Helheim4_RDS.mat'); extent4r = d4r.extentPFA_fromline; dates4r = d4r.dates;
xref4 = d4r.xcoordref; yref4 = d4r.ycoordref;
Cmnx4 = d4r.Cmnx; Cmxx4 = d4r.Cmxx; Cmny4 = d4r.Cmny; Cmxy4 = d4r.Cmxy;
xpfa4 = d4r.xPFA_inrange; ypfa4 = d4r.yPFA_inrange;

%% Plot
figure;
cmap = [pink(max(10)); lines(max(7))];
imagesc(xgrid/1000, fliplr(ygrid/1000), dem);
hold on;
for i = 1:length(xclemall)
    scatter(xclemall{i}/1000, yclemall{i}/1000, 18, cmap(i+1, :), 'filled');
    hold on;
end
%run the first part of climate_MAR_longterm.m here
axis equal;
colormap(gray);
set(gca, 'Ydir', 'Normal');
hold on;
xlim([100 400]);
ylim([-2800 -2500]);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 14);
xlabel('X Coordinate');
ylabel('Y Coordinate');
legend('2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', 'Location', 'Northwest', 'AutoUpdate','off');

hold on;
line([Cmnx1 Cmxx1],[Cmxy1 Cmny1], 'Color', [0 204/255 204/255], 'LineWidth', 3, 'LineStyle', ':');
hold on;
line([Cmnx2 Cmxx2],[Cmxy2 Cmny2], 'Color', [0 204/255 204/255], 'LineWidth', 3, 'LineStyle', ':');
hold on;
line([Cmnx3 Cmxx3],[Cmxy3 Cmny3], 'Color', [0 204/255 204/255], 'LineWidth', 3, 'LineStyle', ':');
hold on;
line([Cmnx4 Cmxx4],[Cmxy4 Cmny4], 'Color', [0 204/255 204/255], 'LineWidth', 3, 'LineStyle', ':');

% plot(xpfa1{end}, ypfa1{end}, 'r', 'LineWidth', 4, 'Color', [0 204/255 204/255]);
% hold on;
% plot(xpfa1{end}(1), ypfa1{end}(1), 'd', 'MarkerSize', 10, 'MarkerFaceColor', [0 204/255 204/255], 'MarkerEdgeColor', 'k');
% hold on;
% plot(xpfa2{end}, ypfa2{end}, 'r', 'LineWidth', 4, 'Color', [0 204/255 204/255]);
% hold on;
% plot(xpfa2{end}(1), ypfa2{end}(1), 'd', 'MarkerSize', 10, 'MarkerFaceColor', [0 204/255 204/255], 'MarkerEdgeColor', 'k');
% hold on;
% plot(xpfa3{end}, ypfa3{end}, 'r', 'LineWidth', 4, 'Color', [0 204/255 204/255]);
% hold on;
% plot(xpfa3{end}(1), ypfa3{end}(1), 'd', 'MarkerSize', 10, 'MarkerFaceColor', [0 204/255 204/255], 'MarkerEdgeColor', 'k');
% hold on;
% plot(xpfa4{end}, ypfa4{end}, 'r', 'LineWidth', 4, 'Color', [0 204/255 204/255]);
% hold on;
% plot(xpfa4{end}(end), ypfa4{end}(end), 'd', 'MarkerSize', 10, 'MarkerFaceColor', [0 204/255 204/255], 'MarkerEdgeColor', 'k');

handaxes2 = axes('position', [0.53 0.14 0.35 0.35]);
greenland('patch','facecolor', [224/255 224/255 224/255]);
set(gca,'visible','off');
hold on;
rectangle('Position',[0,-2800,300,300],'EdgeColor',[213/255 73/255 73/255],...
    'LineWidth',3);

hold on;
line([Cmnx1 Cmxx1],[Cmxy1 Cmny1], 'Color', [0 204/255 204/255], 'LineWidth', 3, 'LineStyle', ':');
hold on;
line([Cmnx2 Cmxx2],[Cmxy2 Cmny2], 'Color', [0 204/255 204/255], 'LineWidth', 3, 'LineStyle', ':');
hold on;
line([Cmnx3 Cmxx3],[Cmxy3 Cmny3], 'Color', [0 204/255 204/255], 'LineWidth', 3, 'LineStyle', ':');
hold on;
line([Cmnx4 Cmxx4],[Cmxy4 Cmny4], 'Color', [0 204/255 204/255], 'LineWidth', 3, 'LineStyle', ':');


