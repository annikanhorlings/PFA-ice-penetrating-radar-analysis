%This script plots CRESIS data location in comparison with known regions to
%make sure its the correct data

%Annika Horlings
%University of Washington
%10 September

%% Load in the data
loc = '/Users/annikahorlings/Desktop';
cd(loc);
profiles_struct = dir('*.mat');

for i = 1:length(profiles_struct)
    profiles{i} = load(profiles_struct(i).name);
    
end

%Load the target profile
loc2 = '/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Computed_vals/IkertivaqN1/all';
cd(loc2);
lat = load('lat_RDS'); lat = lat.lat;
long = load('long_RDS'); long = long.long;


%plot location of all data
figure;
for i = 1:length(long)
    plot(long{i}, lat{i}, 'r:', 'LineWidth', 2);
    hold on;
end

for i = 1:length(profiles_struct)
    plot(profiles{i}.Longitude, profiles{i}.Latitude, 'LineWidth', 2);
    hold on;
    
    m=input('Do you want to save the name of this profile, Y/N [Y]:','s');
    
    if m=='Y'
        names{i} = profiles_struct(i).name;
    end
    
     
    n=input('Do you want to end the loop, Y/N [Y]:','s');
    if n=='Y'
        break
    end

end
xlabel('long');
ylabel('lat');

idx = find(~cellfun(@isempty,names));
for i = 1:length(idx)
    data_names_final{i} = names{idx(i)};
end

cd('/Users/annikahorlings/Desktop/');
save('data_names', 'data_names_final');




