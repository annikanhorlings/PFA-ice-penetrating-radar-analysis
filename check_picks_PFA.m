%This code produces figures for each cresis data file
%Annika Horlings
%Jan 28 2022


%% Profiles:
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Data/Helheim1/RDS');
profiles_struct = dir('*.mat');
profiles = struct2cell(profiles_struct);
[col, lengthp] = size(profiles);
for i = 1:lengthp %get profiles
    p{i} = load(profiles{1, i});
end
for j = 1:lengthp %get desired vars from profiles
    lat{j, 1} = p{j}.lat;
    long{j, 1} = p{j}.long;
    name{j, 1} = profiles{1, j};
    TWT{j, 1} = p{j}.travel_time;
    data{j} = p{j}.data;
    picks{j, 1} = p{j}.picks.samp2;
end

%convert
for k = 1:lengthp
    [xcoord{k}, ycoord{k}] = polarstereo_fwd(lat{k}, long{k}, 6378137.0,0.08181919,70,-45);
    xcoordkm{k} = xcoord{k}/1000;
    ycoordkm{k} = ycoord{k}/1000;
end


for f = 1:lengthp
    for g = 2:length(xcoord{1, f})
        distance{f}(1) = 0;
        distance{f}(g) = sqrt((xcoordkm{f}(g -1) - xcoordkm{f}(g)).^2 + ...
            (ycoordkm{f}(g -1) - ycoordkm{f}(g)).^2) + distance{f}(g-1);
    end
end

%%
for j = 1:lengthp
    h = figure;
    subplot(2, 1, 1);
    imagesc(distance{1, j}, TWT{j}(:, 1), data{j});
    colormap(gray);
    
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 12);
    xlabel('Distance (m)');
    ylabel('Two-way travel time (us)');
    ylim([min(min(picks{j})) - 500 max(max(picks{j})) + 800]);
    title(name{j}, 'Interpreter', 'none');
    
    subplot(2, 1, 2);
    imagesc(data{j});
    colormap(gray);
    hold on;
    plot(picks{j}(1, :), 'r');
    
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 12);
    xlabel('Trace number');
    ylabel('Two-way travel time (us)');
    ylim([min(min(picks{j})) - 500 max(max(picks{j})) + 800]);
    title('with aquifer picks', 'Interpreter', 'none');
    
    filename = [name{j} '.fig'];
    saveas(h, filename);
    
    filename = [name{j} '.png'];
    saveas(h, filename);
end


%%
% figure;
% subplot(2, 1, 1);
% imagesc(data);
% colormap(gray);
% title('AR2018H1');
% set(gca, 'LineWidth', 2);
% ylim([900 1300]);
% caxis([-140 -85]);
% 
% subplot(2, 1, 2);
% imagesc(data);
% colormap(gray);
% hold on;
% plot(picks, 'r');
% title('AR2018H1  with picks');
% set(gca, 'LineWidth', 2);
% ylim([900 1300]);
% caxis([-140 -85]);

%% calculate basic statistics







