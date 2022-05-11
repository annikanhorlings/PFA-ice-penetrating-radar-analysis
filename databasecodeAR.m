%This code produces a "database" for the AR data
%Annika Horlings
%University of Washington
%8 January 2021

%% Load in the data
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Data/AR/picks/all');

profiles_struct = dir('*.mat');
profiles = struct2cell(profiles_struct);

[m, n] = size(profiles);
for j = 1:n
    p{j} = load(profiles{1, j});
    name{j} = profiles{1, j};
end

%% Retrieve variables

for i = 1:n
    picks{i} = p{i}.picks.samp2;
    lat{i} = p{i}.lat;
    long{i} = p{i}.long;
    trace_num{i} = p{i}.trace_num;
    travel_time{i} = p{i}.travel_time;
end

cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/all_picks/all_variables');
save('picks_AR_all.mat', 'picks');
save('lat_AR_all.mat', 'lat');
save('long_AR_all.mat', 'long');
save('trace_num_AR_all', 'trace_num');
save('travel_time_AR_all', 'travel_time');



    


