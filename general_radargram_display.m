%This code tries to intake RDS/AR radar data and ouput the radargrams/other info
%Annika Horlings
%28 June 2021

%% Assign input
%specify:
dir1 = '/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Data/';
dir2 = '/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Computed_vals/';

PFA = 'IkertivaqN1';
dataloc = 'all';
dataname = 'RDS_2015_edit.mat';

%% get data & computed variables
cd(dir1); cd(PFA); cd(dataloc);
p = load(dataname);

if isfield(p, 'picks') == 1
    picks = p.picks.samp2;
    lat = p.lat;
    long = p.long;
    data = p.data;
    
else
    picks = NaN;
    lat = p.lat;
    long = p.long;
    data = p.data;
end
[xcoord, ycoord] = polarstereo_fwd(lat, long, 6378137.0,0.08181919,70,-45);
    
%% map
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Data/DEM/gimpDEM');
g = geotiffinfo('gimpdem_90m_v01.1.tif'); [A, R] = geotiffread('gimpdem_90m_v01.1.tif');
A = double(A); dem = A;
Xmin = g.CornerCoords.X(1);
Xmax = g.CornerCoords.X(2);
Ymax = g.CornerCoords.Y(1);
Ymin = g.CornerCoords.Y(3);
xgrid = (Xmin:g.Width:Xmax);
ygrid = (Ymin:g.Height:Ymax);

%% Plots
fig = figure;
subplot(3, 2, 1);
imagesc(xgrid/1000, fliplr(ygrid/1000), dem); hold on;
colormap(gray);
set(gca, 'Ydir', 'Normal');
axis equal;
xlim([0 450]);
ylim([-2700 -2550]);
xlabel('X Coordinate');
ylabel('Y Coordinate');
%set(gca, 'Fontsize', 12);
set(gca, 'LineWidth', 1.5);
hold on;
plot(xcoord/1000, ycoord/1000, 'Color', 'r', 'LineWidth', 2);
hold on;
plot(xcoord(1)/1000, ycoord(1)/1000, 'rx');

subplot(3, 2, [3 4 5 6]);
imagesc(data);
hold on;
plot(picks, 'Color', 'b', 'LineWidth', 2)
colormap(gray);
caxis([min(min(p.data)) max(max(p.data))]);
%     ylim([0 600]);
title(dataname, 'Interpreter', 'none');
set(gca, 'LineWidth', 1.5);

    
    
    