% This code loads in and manipulates Clem's pick data
% %Author: Annika Horlings
% University of Washington
% Last updated: 2021

% Note: only years 2010, 2011, 2012, 2013, 2014 and 2017 are using AR. 
%Years 2015 and 2016 were done using RDS (no AR data for these two years).

%this is for the region where the MAR climate output is averaged over
%for the paper
%% Load in
%For each files 5 columns: lat,lon,elevation,depth,flight date
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Data/from_Clem');
file10 = csvread('FirnAquiferDetections2010.csv', 1, 0);
file11 = csvread('FirnAquiferDetections2011.csv', 1, 0);
file12 = csvread('FirnAquiferDetections2012.csv', 1, 0);
file13 = csvread('FirnAquiferDetections2013.csv', 1, 0);
file14 = csvread('FirnAquiferDetections2014.csv', 1, 0);
file15 = csvread('FirnAquiferDetections2015.csv', 1, 0);
file16 = csvread('FirnAquiferDetections2016.csv', 1, 0);
file17 = csvread('FirnAquiferDetections2017.csv', 1, 0);

%segregate variables
lat10 = file10(:, 1); %latitude
lat11 = file11(:, 1);
lat12 = file12(:, 1);
lat13 = file13(:, 1);
lat14 = file14(:, 1);
lat15 = file15(:, 1);
lat16 = file16(:, 1);
lat17 = file17(:, 1);

long10 = file10(:, 2);
long11 = file11(:, 2);
long12 = file12(:, 2);
long13 = file13(:, 2);
long14 = file14(:, 2);
long15 = file15(:, 2);
long16 = file16(:, 2);
long17 = file17(:, 2);

elev10 = file10(:, 3);
elev11 = file11(:, 3);
elev12 = file12(:, 3);
elev13 = file13(:, 3);
elev14 = file14(:, 3);
elev15 = file15(:, 3);
elev16 = file16(:, 3);
elev17 = file17(:, 3);

depth10 = file10(:, 4);
depth11 = file11(:, 4);
depth12 = file12(:, 4);
depth13 = file13(:, 4);
depth14 = file14(:, 4);
depth15 = file15(:, 4);
depth16 = file16(:, 4);
depth17 = file17(:, 4);

date10 = file10(:, 5);
date11 = file11(:, 5);
date12 = file12(:, 5);
date13 = file13(:, 5);
date14 = file14(:, 5);
date15 = file15(:, 5);
date16 = file16(:, 5);
date17 = file17(:, 5);

%% segregate into basin locations

bounds{1} = [-42.6 -38.3 65.00 66.49];

for i = 1:length(bounds)
    idx10{i} = find (lat10 > bounds{i}(3) & lat10 < bounds{i}(4) & long10 > bounds{i}(1) & long10 < bounds{i}(2));
    idx11{i} = find (lat11 > bounds{i}(3) & lat11 < bounds{i}(4) & long11 > bounds{i}(1) & long11 < bounds{i}(2));
    idx12{i} = find (lat12 > bounds{i}(3) & lat12 < bounds{i}(4) & long12 > bounds{i}(1) & long12 < bounds{i}(2));
    idx13{i} = find (lat13 > bounds{i}(3) & lat13 < bounds{i}(4) & long13 > bounds{i}(1) & long13 < bounds{i}(2));
    idx14{i} = find (lat14 > bounds{i}(3) & lat14 < bounds{i}(4) & long14 > bounds{i}(1) & long14 < bounds{i}(2));
    idx15{i} = find (lat15 > bounds{i}(3) & lat15 < bounds{i}(4) & long15 > bounds{i}(1) & long15 < bounds{i}(2));
    idx16{i} = find (lat16 > bounds{i}(3) & lat16 < bounds{i}(4) & long16 > bounds{i}(1) & long16 < bounds{i}(2));
    idx17{i} = find (lat17 > bounds{i}(3) & lat17 < bounds{i}(4) & long17 > bounds{i}(1) & long17 < bounds{i}(2));
end

for i = 1:length(bounds)
    lat10basins{i} = lat10(idx10{i});
    long10basins{i} = long10(idx10{i});
    elev10basins{i} = elev10(idx10{i});
    depth10basins{i} = depth10(idx10{i});
    date10basins{i} = date10(idx10{i});
end

for i = 1:length(bounds)
    lat11basins{i} = lat11(idx11{i});
    long11basins{i} = long11(idx11{i});
    elev11basins{i} = elev11(idx11{i});
    depth11basins{i} = depth11(idx11{i});
    date11basins{i} = date11(idx11{i});
end

for i = 1:length(bounds)
    lat12basins{i} = lat12(idx12{i});
    long12basins{i} = long12(idx12{i});
    elev12basins{i} = elev12(idx12{i});
    depth12basins{i} = depth12(idx12{i});
    date12basins{i} = date12(idx12{i});
end

for i = 1:length(bounds)
    lat13basins{i} = lat13(idx13{i});
    long13basins{i} = long13(idx13{i});
    elev13basins{i} = elev13(idx13{i});
    depth13basins{i} = depth13(idx13{i});
    date13basins{i} = date13(idx13{i});
end

for i = 1:length(bounds)
    lat14basins{i} = lat14(idx14{i});
    long14basins{i} = long14(idx14{i});
    elev14basins{i} = elev14(idx14{i});
    depth14basins{i} = depth14(idx14{i});
    date14basins{i} = date14(idx14{i});
end

for i = 1:length(bounds)
    lat15basins{i} = lat15(idx15{i});
    long15basins{i} = long15(idx15{i});
    elev15basins{i} = elev15(idx15{i});
    depth15basins{i} = depth15(idx15{i});
    date15basins{i} = date15(idx15{i});
end

for i = 1:length(bounds)
    lat16basins{i} = lat16(idx16{i});
    long16basins{i} = long16(idx16{i});
    elev16basins{i} = elev16(idx16{i});
    depth16basins{i} = depth16(idx16{i});
    date16basins{i} = date16(idx16{i});
end

for i = 1:length(bounds)
    lat17basins{i} = lat17(idx17{i});
    long17basins{i} = long17(idx17{i});
    elev17basins{i} = elev17(idx17{i});
    depth17basins{i} = depth17(idx17{i});
    date17basins{i} = date17(idx17{i});
end


%% save things
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Data/from_Clem/MATfiles/by_region');

%by basin
save('lat_by_basins', 'lat10basins', 'lat11basins', 'lat12basins', 'lat13basins', 'lat14basins', 'lat15basins', ...
    'lat16basins', 'lat17basins');

save('long_by_basins', 'long10basins', 'long11basins', 'long12basins', 'long13basins', 'long14basins', 'long15basins', ...
    'long16basins', 'long17basins');

save('elev_by_basins', 'elev10basins', 'elev11basins', 'elev12basins', 'elev13basins', 'elev14basins', 'elev15basins', ...
    'elev16basins', 'elev17basins');

save('depth_by_basins', 'depth10basins', 'depth11basins', 'depth12basins', 'depth13basins', 'depth14basins', 'depth15basins', ...
    'depth16basins', 'depth17basins');

save('date_by_basins', 'date10basins', 'date11basins', 'date12basins', 'date13basins', 'date14basins', 'date15basins', ...
    'date16basins', 'date17basins');

%by year
save('lat_by_year', 'lat10', 'lat11', 'lat12', 'lat13', 'lat14', 'lat15', 'lat16', 'lat17');

save('long_by_year', 'long10', 'long11', 'long12', 'long13', 'long14', 'long15', 'long16', 'long17');

save('elev_by_year', 'elev10', 'elev11', 'elev12', 'elev13', 'elev14', 'elev15', 'elev16', 'elev17');

save('depth_by_year', 'depth10', 'depth11', 'depth12', 'depth13', 'depth14', 'depth15','depth16', 'depth17');

save('date_by_year', 'date10', 'date11', 'date12', 'date13', 'date14', 'date15','date16', 'date17');




