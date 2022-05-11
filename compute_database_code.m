% This code computes the database for the four firn aquifers for the paper
% Horlings et al. (2022).

%Author: Annika Horlings
%University of Washington
%%Last updated: 14 April 2022

%% Helheim 1 RDS
%1993 1998 2001 2003 2005 2006 2012 2013 2014 2017 2018
clear;
cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/Helheim1/RDS')
profiles_struct = dir('*.mat');
profiles = struct2cell(profiles_struct);
[m, n] = size(profiles);
for i = 1:n
    data{i} = load(profiles{1, i});
end

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/Helheim1');
Cmnx = load('Cmnx.mat'); Cmnx = Cmnx.Cmnx;
Cmxx = load('Cmxx.mat'); Cmxx = Cmxx.Cmxx;
Cmny = load('Cmny.mat'); Cmny = Cmny.Cmny;
Cmxy = load('Cmxy.mat'); Cmxy = Cmxy.Cmxy;
rdrtype = 0; %either 0 for RDS or 1 or AR
number = 11;
[xcoordinate, ycoordinate, xPFA_inrange, yPFA_inrange, extentPFA_fromline, ...
    distance_fromline, idxPFA_inrange, whereidxofPFAinrange, picks_inrange, ...
    xcoordref,  ycoordref] = project_PFA_extent(data, ...
    rdrtype, number, Cmnx, Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/Helheim1/RDS')
[dates, depth_picks, depthpicks_inrange, depth_vec, extent, latitude, ...
    longitude, xPFAir, yPFAir] = calc_database(data, rdrtype, Cmnx, ...
    Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/Helheim1');
save('database_Helheim1_RDS', 'Cmnx', 'Cmny', 'Cmxx', 'Cmxy', 'dates', 'depth_picks', 'depth_vec', 'depthpicks_inrange', ...
    'distance_fromline', 'extent', 'extentPFA_fromline', 'idxPFA_inrange', 'latitude', 'longitude', 'picks_inrange', ...
    'whereidxofPFAinrange', 'xPFA_inrange', 'xPFAir', 'xcoordinate', 'xcoordref', 'yPFA_inrange', 'yPFAir', 'ycoordinate', 'ycoordref');

%% Helheim 1 AR
%[2010 2011 2012 2013 2018]
clear;
cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/Helheim1/AR')
profiles_struct = dir('*.mat');
profiles = struct2cell(profiles_struct);
[m, n] = size(profiles);
for i = 1:n
    data{i} = load(profiles{1, i});
end

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/Helheim1');
Cmnx = load('Cmnx.mat'); Cmnx = Cmnx.Cmnx;
Cmxx = load('Cmxx.mat'); Cmxx = Cmxx.Cmxx;
Cmny = load('Cmny.mat'); Cmny = Cmny.Cmny;
Cmxy = load('Cmxy.mat'); Cmxy = Cmxy.Cmxy;
rdrtype = 1; %either 0 for RDS or 1 or AR
number = 5;
[xcoordinate, ycoordinate, xPFA_inrange, yPFA_inrange, extentPFA_fromline, ...
    distance_fromline, idxPFA_inrange, whereidxofPFAinrange, picks_inrange, ...
    xcoordref,  ycoordref] = project_PFA_extent(data, ...
    rdrtype, number, Cmnx, Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/Helheim1/AR')
[dates, depth_picks, depthpicks_inrange, depth_vec, extent, latitude, ...
    longitude, xPFAir, yPFAir] = calc_database(data, rdrtype, Cmnx, ...
    Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/Helheim1');
save('database_Helheim1_AR', 'Cmnx', 'Cmny', 'Cmxx', 'Cmxy', 'dates', 'depth_picks', 'depth_vec', 'depthpicks_inrange', ...
    'distance_fromline', 'extent', 'extentPFA_fromline', 'idxPFA_inrange', 'latitude', 'longitude', 'picks_inrange', ...
    'whereidxofPFAinrange', 'xPFA_inrange', 'xPFAir', 'xcoordinate', 'xcoordref', 'yPFA_inrange', 'yPFAir', 'ycoordinate', 'ycoordref');


%% IN1 RDS
%[1998 2001 2002 2006 2012 2014 2017]
clear;
cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/IkertivaqN1/RDS')
profiles_struct = dir('*.mat');
profiles = struct2cell(profiles_struct);
[m, n] = size(profiles);
for i = 1:n
    data{i} = load(profiles{1, i});
end

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/IkertivaqN1/');
Cmnx = load('Cmnx.mat'); Cmnx = Cmnx.Cmnx;
Cmxx = load('Cmxx.mat'); Cmxx = Cmxx.Cmxx;
Cmny = load('Cmny.mat'); Cmny = Cmny.Cmny;
Cmxy = load('Cmxy.mat'); Cmxy = Cmxy.Cmxy;
rdrtype = 0; %either 0 for RDS or 1 or AR
number = 7;
[xcoordinate, ycoordinate, xPFA_inrange, yPFA_inrange, extentPFA_fromline, ...
    distance_fromline, idxPFA_inrange, whereidxofPFAinrange, picks_inrange, ...
    xcoordref,  ycoordref] = project_PFA_extent(data, ...
    rdrtype, number, Cmnx, Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/IkertivaqN1/RDS')
[dates, depth_picks, depthpicks_inrange, depth_vec, extent, latitude, ...
    longitude, xPFAir, yPFAir] = calc_database(data, rdrtype, Cmnx, ...
    Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/IkertivaqN1');
save('database_IkertivaqN1_RDS', 'Cmnx', 'Cmny', 'Cmxx', 'Cmxy', 'dates', 'depth_picks', 'depth_vec', 'depthpicks_inrange', ...
    'distance_fromline', 'extent', 'extentPFA_fromline', 'idxPFA_inrange', 'latitude', 'longitude', 'picks_inrange', ...
    'whereidxofPFAinrange', 'xPFA_inrange', 'xPFAir', 'xcoordinate', 'xcoordref', 'yPFA_inrange', 'yPFAir', 'ycoordinate', 'ycoordref');

%% IN1 AR
%[2011 2012 2014]
clear;
cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/IkertivaqN1/AR')
profiles_struct = dir('*.mat');
profiles = struct2cell(profiles_struct);
[m, n] = size(profiles);
for i = 1:n
    data{i} = load(profiles{1, i});
end

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/IkertivaqN1/');
Cmnx = load('Cmnx.mat'); Cmnx = Cmnx.Cmnx;
Cmxx = load('Cmxx.mat'); Cmxx = Cmxx.Cmxx;
Cmny = load('Cmny.mat'); Cmny = Cmny.Cmny;
Cmxy = load('Cmxy.mat'); Cmxy = Cmxy.Cmxy;
rdrtype = 1; %either 0 for RDS or 1 or AR
number = 3;
[xcoordinate, ycoordinate, xPFA_inrange, yPFA_inrange, extentPFA_fromline, ...
    distance_fromline, idxPFA_inrange, whereidxofPFAinrange, picks_inrange, ...
    xcoordref,  ycoordref] = project_PFA_extent(data, ...
    rdrtype, number, Cmnx, Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/IkertivaqN1/AR')
[dates, depth_picks, depthpicks_inrange, depth_vec, extent, latitude, ...
    longitude, xPFAir, yPFAir] = calc_database(data, rdrtype, Cmnx, ...
    Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/IkertivaqN1');
save('database_IkertivaqN1_AR', 'Cmnx', 'Cmny', 'Cmxx', 'Cmxy', 'dates', 'depth_picks', 'depth_vec', 'depthpicks_inrange', ...
    'distance_fromline', 'extent', 'extentPFA_fromline', 'idxPFA_inrange', 'latitude', 'longitude', 'picks_inrange', ...
    'whereidxofPFAinrange', 'xPFA_inrange', 'xPFAir', 'xcoordinate', 'xcoordref', 'yPFA_inrange', 'yPFAir', 'ycoordinate', 'ycoordref');


%% Helheim4 RDS
clear;
%[2001 2002 2006 2010 2014 2017];
cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/Helheim4/RDS')
profiles_struct = dir('*.mat');
profiles = struct2cell(profiles_struct);
[m, n] = size(profiles);
for i = 1:n
    data{i} = load(profiles{1, i});
end

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/Helheim4/');
Cmnx = load('Cmnx.mat'); Cmnx = Cmnx.Cmnx;
Cmxx = load('Cmxx.mat'); Cmxx = Cmxx.Cmxx;
Cmny = load('Cmny.mat'); Cmny = Cmny.Cmny;
Cmxy = load('Cmxy.mat'); Cmxy = Cmxy.Cmxy;
rdrtype = 0; %either 0 for RDS or 1 or AR
number = 6;
[xcoordinate, ycoordinate, xPFA_inrange, yPFA_inrange, extentPFA_fromline, ...
    distance_fromline, idxPFA_inrange, whereidxofPFAinrange, picks_inrange, ...
    xcoordref,  ycoordref] = project_PFA_extent(data, ...
    rdrtype, number, Cmnx, Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/Helheim4/RDS')
[dates, depth_picks, depthpicks_inrange, depth_vec, extent, latitude, ...
    longitude, xPFAir, yPFAir] = calc_database(data, rdrtype, Cmnx, ...
    Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/Helheim4');
save('database_Helheim4_RDS', 'Cmnx', 'Cmny', 'Cmxx', 'Cmxy', 'dates', 'depth_picks', 'depth_vec', 'depthpicks_inrange', ...
    'distance_fromline', 'extent', 'extentPFA_fromline', 'idxPFA_inrange', 'latitude', 'longitude', 'picks_inrange', ...
    'whereidxofPFAinrange', 'xPFA_inrange', 'xPFAir', 'xcoordinate', 'xcoordref', 'yPFA_inrange', 'yPFAir', 'ycoordinate', 'ycoordref');

%% Koge Bugt RDS
%[1998 2001 2002 2005 2013]
clear;
cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/KogeBugtS1/RDS')
profiles_struct = dir('*.mat');
profiles = struct2cell(profiles_struct);
[m, n] = size(profiles);
for i = 1:n
    data{i} = load(profiles{1, i});
end

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/KogeBugtS1/');
Cmnx = load('Cmnx.mat'); Cmnx = Cmnx.Cmnx;
Cmxx = load('Cmxx.mat'); Cmxx = Cmxx.Cmxx;
Cmny = load('Cmny.mat'); Cmny = Cmny.Cmny;
Cmxy = load('Cmxy.mat'); Cmxy = Cmxy.Cmxy;
rdrtype = 0; %either 0 for RDS or 1 or AR
number = 5;
[xcoordinate, ycoordinate, xPFA_inrange, yPFA_inrange, extentPFA_fromline, ...
    distance_fromline, idxPFA_inrange, whereidxofPFAinrange, picks_inrange, ...
    xcoordref,  ycoordref] = project_PFA_extent(data, ...
    rdrtype, number, Cmnx, Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Data/KogeBugtS1/RDS')
[dates, depth_picks, depthpicks_inrange, depth_vec, extent, latitude, ...
    longitude, xPFAir, yPFAir] = calc_database(data, rdrtype, Cmnx, ...
    Cmxx, Cmny, Cmxy);

cd('/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database/KogeBugtS1');
save('database_KogeButS1_RDS', 'Cmnx', 'Cmny', 'Cmxx', 'Cmxy', 'dates', 'depth_picks', 'depth_vec', 'depthpicks_inrange', ...
    'distance_fromline', 'extent', 'extentPFA_fromline', 'idxPFA_inrange', 'latitude', 'longitude', 'picks_inrange', ...
    'whereidxofPFAinrange', 'xPFA_inrange', 'xPFAir', 'xcoordinate', 'xcoordref', 'yPFA_inrange', 'yPFAir', 'ycoordinate', 'ycoordref');
















