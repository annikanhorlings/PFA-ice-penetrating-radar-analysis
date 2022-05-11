%This figure plots the extra climate vars from MAR that are not included in
%the paper Horlings et al. (2022), like temperature.

%Author: Annika Horlings
%University of Washington
%Last updated: 10 May 2022

%% Load in parameters
%These were already averaged over the area 
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Database_climate');
clim = load('climupstream.mat'); daysonset = clim.daysonset; 
datesmar = clim.datesmar; tempw = clim.meantempH1annualwinter; 
tempall = clim.meantempH1annual; 

%% Plot
% -winter temperature, mean annual
% -temperature, mean annual

fclim = figure;
cmap = intense(max(10));

%1 -  %winter temperature, mean annual
subplot(2, 1, 1);
plot(datesmar, tempall, '-', 'LineWidth', 2, 'Color', 'k', ...
    'MarkerFaceColor', 'k'); 
hold on;
[sddaysonset, mndecdaysonset, datesmean] = pltstats(tempall, datesmar);
cdaysonset = findchangepts(tempall, 'Statistic', 'mean');
for i = 1:length(mndecdaysonset)
    p = patch('vertices', [datesmean{i}(1), ...
        mndecdaysonset(i) - sddaysonset(i); ...
        datesmean{i}(1), mndecdaysonset(i) + sddaysonset(i); ...
        datesmean{i}(end)+1, mndecdaysonset(i) + sddaysonset(i); ...
        datesmean{i}(end)+1, mndecdaysonset(i) - sddaysonset(i)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'k', ...
          'FaceAlpha', 0.3, ...
          'EdgeColor', 'none');
end
ylabel('Annual Temperature (C)'); 
xlim([1947.5 2016.5]);
set(gca,'XTickLabel',[]);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 14);
ylim([min(tempall) - 1   max(tempall) + 1 ]);
grid on;
plot(datesmar, tempall, '-', 'LineWidth', 2, 'Color', 'k', ...
    'MarkerFaceColor', 'k'); 
hold on;
plot(ones(1, 100).*datesmar(cdaysonset), linspace(min(tempall) - ...
    0.05, max(tempall) + 0.05, 100), 'g');
hold on;
plot(linspace(datesmar(1), datesmar(cdaysonset), cdaysonset), ...
    ones(1, cdaysonset).*mean(tempall(1:cdaysonset)), 'r');
hold on;
plot(linspace(datesmar(cdaysonset), datesmar(end), ...
    length(tempall) - cdaysonset), ones(1, length(tempall)-...
    cdaysonset).*mean(tempall(cdaysonset:end)), 'r');

%2- winter T
subplot(2, 1, 2); 
plot(datesmar, tempw, '-', 'LineWidth', 2, 'Color', 'k', ...
    'MarkerFaceColor', 'k'); 
hold on;
[sddaysduration, mndecdaysduration, datesmean] = ...
    pltstats(tempw, datesmar); 
cdaysduration = findchangepts(tempw, 'Statistic', 'mean');
for i = 1:length(mndecdaysduration)
    p = patch('vertices', [datesmean{i}(1), ...
        mndecdaysduration(i) - sddaysduration(i); ...
        datesmean{i}(1), mndecdaysduration(i) + sddaysduration(i); ...
        datesmean{i}(end)+1, mndecdaysduration(i) + sddaysduration(i); ...
        datesmean{i}(end)+1, mndecdaysduration(i) - sddaysduration(i)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'k', ...
          'FaceAlpha', 0.3, ...
          'EdgeColor', 'none');
end
xlabel('Year');
ylabel('Winter Temperature (C)'); 
xlim([1947.5 2016.5]);
ylim([min(tempw) - 1  max(tempw) + 1]);
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 14);
grid on;
plot(datesmar, tempw, '-', 'LineWidth', 2, 'Color', 'k', ...
    'MarkerFaceColor', 'k'); 
plot(ones(1, 100).*datesmar(cdaysduration), linspace(min(tempw) - ...
    0.05, max(tempw) + 0.05, 100), 'g');
hold on;
plot(linspace(datesmar(1), datesmar(cdaysduration), cdaysduration), ...
    ones(1, cdaysduration).*...
    mean(tempw(1:cdaysduration)), 'r');
hold on;
plot(linspace(datesmar(cdaysduration), datesmar(end), ...
length(tempw) - cdaysduration), ones(1, length(tempw)- ...
cdaysduration).*mean(tempw(cdaysduration:end)), 'r');

%set figure size and position
set(fclim, 'Position', [560   536   705   412]);

