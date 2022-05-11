%This code produces a "database" for the RDS data
%Annika Horlings
%University of Washington
%8 January 2021

%% Load in the data
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Data/Radar/allpicks');

profiles_struct = dir('*.mat');
profiles = struct2cell(profiles_struct);

[m, n] = size(profiles);
for j = 1:n
    p{j} = load(profiles{1, j});
    name{j} = profiles{1, j};
end

%% Retrieve variables

for i = 1:n
    if isfield(p{i}, 'picks') == 0
    picks{i} = p{i}.samp2;
    
    else
    lat{i} = p{i}.lat;
    long{i} = p{i}.long;
    end
end

cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/all_picks/all_variables');
save('picks_RDS_all.mat', 'picks');
save('lat_RDS_all.mat', 'lat');
save('long_RDS_all.mat', 'long');
% save('trace_num_RDS_all', 'trace_num');
% save('travel_time_RDS_all', 'travel_time');


