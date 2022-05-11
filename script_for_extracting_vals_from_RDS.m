%This code produces a "database" for the RDS data at Helheim 1
%Annika Horlings
%University of Washington
%8 January 2021

%The following values are either extracted or computed:
%     data
%     dates
%     distance
%     depth
%     surface elevation
%     extent
%     elevation
%     latitude
%     longitude
%     xcoordinate
%     ycoordinate
%     picks (in TWT)
%     distance W
%     distance E
%     index W
%     index E

%% Load in the data
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Data/Helheim1/RDS');

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
    data{i} = p{i}.data;
    travel_time{i} = p{i}.travel_time;
end

%% Calculate variables:

%(1) LAT/LONG
% convert lat/long to x and y coords polar stereographic
%using polartereo_fwd, a function created Andy Bliss
for i = 1:n
    [xcoordinate{i}, ycoordinate{i}] = polarstereo_fwd(lat{i}, long{i}, 6378137.0,0.08181919,70,-45);
end

%(2)REPEAT FLIGHT INDEX
% find the indices that are within the repeated flight lines rectangle
% (determined by visual inspection)
for i = 1:n
    repeatflightindex{i} = find(ycoordinate{i}/1000 < -2582 & ycoordinate{i}/1000 > -2584 ...
        & xcoordinate{i}/1000 > 235 & xcoordinate{i}/1000 < 280);
end
for i = 1:n
    latitude_repeat{i} = lat{i}(repeatflightindex{i});
    longitude_repeat{i} = long{i}(repeatflightindex{i});
    xcoordinate_repeat{i} = xcoordinate{i}(repeatflightindex{i});
    ycoordinate_repeat{i} = ycoordinate{i}(repeatflightindex{i});
    picks_repeat{i} = picks{i}(repeatflightindex{i});
end

%(3) INDEX W/INDEX E
% calculate where we should cut off each of the flight lines
% find distance associated with each of these
for k = 1:n
    [kxwlim_repeat(k), indexwestlim_repeat(k)] = (min(abs((xcoordinate_repeat{k}/1000 - 235))));
end

for l = 1:n
    [kxelim_repeat(l), indexeastlim_repeat(l)] = (min(abs((xcoordinate_repeat{l}/1000 - 280))));
end

%(4) DISTANCE
v1 = [235*1000,-2582*1000,0]; %a line going perpendicular to the main flight lines here
v2 = [235*1000,-2585*1000,0];
for i = 1:n % calculate the distance along the flight paths
    
    for j = 1:length(xcoordinate_repeat{1, i})
        distance_repeat{i}(j) = point_to_line([xcoordinate_repeat{i}(j), ycoordinate_repeat{i}(j), 0], v1, v2);
    
    end
end

% for f = 1:n
%     distancewestlimit_repeat{f} = distance_repeat{f}(indexwestlim_repeat(f));
%     distanceeastlimit_repeat{f} = distance_repeat{f}(indexeastlim_repeat(f));
% end
% 
% for i = 1:n % start the distance at the beginning of the area of interest
%     if distancewestlimit_repeat{i} < distanceeastlimit_repeat{i}
%         distancecorrected_repeat{i} = distance_repeat{i} - distancewestlimit_repeat{i};
%     else
%         distancecorrected_repeat{i} = distance_repeat{i} - distanceeastlimit_repeat{i};
%     end
% end

%(5) AQUIFER EXTENT
%Calculate areas where there are NaN values in the picks (this is where the
%aquifer is)
for i = 1:n
    knan{i} = find(isnan(picks_repeat{i})); %Helheim 1
    kdiff{i} = find(diff(knan{i}) > 1);
    kloc{i} = zeros(1, length(picks_repeat{i}));
    kloc{i}(knan{i}) = 1;
    for j = 1:length(kloc{i}) - 5
        if (kloc{i}(1, j) + kloc{i}(1, j+1) + kloc{i}(1, j+2) + kloc{i}(1, j+3) + kloc{i}(1, j+4) < 5)
            kloc{i}(j) = 0;
        end
    end
    kloc{i}(end-5:end) = 0;
    kloc{i}(kloc{i} == 0) = NaN;
end
for i = 1:n
    extent{1, i} = distance_repeat{1, i}.*(kloc{i});
end

%(6) DATES
dates = {1993 1998 2001 2003 2005 2006 2008 2011 2012 2013 2014 2015 2017 2018};

%(7) ELEVATION (not in RDS data)


%% (8) DEPTH
%This may take awhile:
for j = 1:n
    for i = 1:length(travel_time{j})
        depth{j}(i) = TZ_firn_bdot_T(travel_time{j}(i)./(10^6), 2, -8);
        end
end

%% save 
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Computed_vals/Helheim1/all');
save('RDS_database', 'data', 'dates', 'depth', 'distance_repeat', 'distancecorrected_repeat', 'extent', 'lat', 'latitude_repeat', ...
    'long', 'longitude_repeat', 'xcoordinate', 'xcoordinate_repeat', 'ycoordinate', 'ycoordinate_repeat', 'picks', ...
    'indexwestlim_repeat', 'indexeastlim_repeat', 'repeatflightindex', 'travel_time');
