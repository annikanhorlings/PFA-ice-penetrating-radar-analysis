%This code calculates the climate and cold content from the MAR regional 
%climate model at the aquifer sites and surroundin sites for the paper, 
%given the geographic locations.
%
%Author: Annika Horlings
%Last updated: 28 February 2022
%University of Washington
%
%Outputs:
%   annual melt
%   annual winter temperatures
%   annual melt/smb ratio
%   annual melt between radargram detections
%   onset of melt - above a threshold of 0.5 mm per day
%   melt duration - from the first and last day of melt above a threshold of
%       0.5 mm per day
%   number of days - above the melt threshold of 0.5 mm per day
%   winter balance - of previous winter
%   mean annual summer temperature
%   annual rainfall
%   annual rainfall between radargram detections
%
%Data dependencies:
%MAR regional climate model output - https://arcticdata.io/catalog/view/doi%3A10.18739%2FA2H12V80V
%
%Function dependencies:
%polarstereo_fwd.m
%coldcontent.m

%% 1 Specify user parameters

% specify directories 
%%% use absolute paths
fnmar = '/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Data/MAR'; 
fnradar = '/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Computed_vals/Helheim1/all';
fnsave = '/Users/annikahorlings/Desktop/finalized_material_for_paper_final/Database_climate';

% geographic information for target sites
%%% this set of coordinates is divided as follows: (1:3) seward limit of
%%% aquifer locations; (4:6) upstream aquifer locatoin; (7:9) northern
%%% site, inland extent of aquifers there; (10:12) west of the aquifer
%%% sites in the paper; (13:15) even more west of the aquifer sites in the
%%% paper; (16:19) southern site, inland extent of aquifers there.
lonextall = [-38.928 -39.867 -41.503 -39.3954 -40.4463 -41.8572 -34.374 ...
    -34.552 -34.762 -39.907 -40.946 -42.375 -40.407 -41.446 -42.875 ...
    -43.882 -43.773 -43.670 -43.691];
latextall = [66.354 65.765 65.106 66.3610 65.8776 65.2031 68.031 67.620 ...
    67.507 66.387 65.8878 65.194 66.3887 65.8878 65.194 62.767 62.584 ...
    62.301 62.086];
lonext = lonextall(16:19); %specify which location
latext = latextall(16:19);

% dates
datesmar = 1948:1:2016;

% get date of when the radar data was taken
cd(fnradar); doyin = load('Helheim1_dayofyear.mat'); doy = doyin.doy; 
date = doyin.dv; idxr = doyin.idxofrdgms; doy = 200;

% specify file output name
fnout = 'climsouth';

%% 2 Climate analysis

% get climate RCM output
cd(fnmar); mar_struct = dir('*.nc'); N = length(mar_struct);
for b = 1:N %load MAR variables
    thisfilename{b} = string(mar_struct(b).name); 
    latitude{b} = ncread(thisfilename{b},'LAT');
    longitude{b} = ncread(thisfilename{b},'LON');
    time{b} = ncread(thisfilename{b},'TIME'); % hours since 1901-01-15 00:00:00
    smb{b} = ncread(thisfilename{b}, 'SMB');
    temp{b} = ncread(thisfilename{b}, 'TTZ');
    melt{b} = ncread(thisfilename{b}, 'ME');
    rain{b} = ncread(thisfilename{b}, 'RF');
    tempice{b} = ncread(thisfilename{b}, 'TI1'); %these use the depth "OUTLAY"
    rhos{b} = ncread(thisfilename{b}, 'RO1'); %these use the depth "OUTLAY"
    outlay{b} = ncread(thisfilename{b}, 'OUTLAY');
end
    

for ik = 1:length(datesmar)
    % lat/long to x/y coordinates of climate output
    [xcoord{ik}, ycoord{ik}] = polarstereo_fwd(latitude{ik}, ...
        longitude{ik}, 6378137.0,0.08181919,70,-45);
    xcoordk{ik} = xcoord{ik}/1000;
    ycoordk{ik} = ycoord{ik}/1000;
    [m, n] = size(ycoordk{ik});
    
    % resizing & taking the mean for each year
    meltmean{ik} = mean(melt{ik}(:, :, 1, :), 4); ...
        meltmean{ik}(find(meltmean{ik} > 9)) = NaN;
    smbmean{ik} = mean(smb{ik}(:, :, 1, :), 4);
    tempmean{ik} = mean(temp{ik}(:, :, 1, :), 4);
end


for iq = 1:length(lonext)
    % lat/long to x/y coordinates for target locations from user input
    [xrefall(iq),yrefall(iq)] = polarstereo_fwd(latext(iq), lonext(iq), 6378137.0,0.08181919,70,-45);
    xrefallk(iq) = xrefall(iq)/1000; %from m to km
    yrefallk(iq) = yrefall(iq)/1000;
    
    
    % find indices within climate output closest to target locations
    [rowidxpre{iq}, colidxpre{iq}] = find(abs(xrefallk(iq)*ones(m,n) - xcoordk{1}) < 15 & ...
        abs(yrefallk(iq)*ones(m, n) - ycoordk{1}) < 15);
    if length(rowidxpre{iq}) > 1
        rowidx(iq) = rowidxpre{iq}(1);
        colidx(iq) = colidxpre{iq}(1);
    elseif isempty(rowidxpre{iq}) == 1
        rowidx(iq) = NaN;
        colidx(iq) = NaN;
    else
        rowidx(iq) = rowidxpre{iq};
        colidx(iq) = colidxpre{iq};
    end
end
rowidxnon= rowidx(~isnan(rowidx));
colidxnon= colidx(~isnan(colidx));
rowsubscript = ind2sub(size(ycoordk), rowidxnon);
colsubscript = ind2sub(size(ycoordk), colidxnon);

% calculate climate metrics for all years

for it = 1:length(datesmar) - 1
    %%% melt:smb ratio (last winter's smb)
    smbmeanwinterbef{it + 1} = (mean(smb{it}(:, :, 1, 300:365), 4) + mean(smb{it + 1}(:, :, 1, 1:90), 4))/2;   %Oct-March
    meltsmbratio{it+1} = meltmean{it+1}./smbmeanwinterbef{it+1};
    
    for im = 1:length(rowsubscript)
        wintersmbH1{it + 1} = smbmeanwinterbef{it+1}(rowsubscript(im), colsubscript(im));
        meltsmbratioH1(it + 1) = meltsmbratio{it+1}(rowsubscript(im), colsubscript(im));
    end
    cumwintersmb{it + 1} = cumsum(wintersmbH1{it + 1});
    sumwintersmb(it + 1) = sum(wintersmbH1{it + 1});  %Oct-March
end

for ix = 1:length(datesmar)
    %%% melt metrics
    for ir = 1:length(rowsubscript)
        meltH1daily{ix}(ir, 1, 1, :) = melt{ix}(rowsubscript(ir), colsubscript(ir), 1, :); %: x 1 x 1 x 365
    end
    meanmeltH1daily{ix} = squeeze(mean(meltH1daily{ix}, 1))/1000;
    onsetidx{ix} = find(meanmeltH1daily{ix}(1:doy) > 0.0005); %0.5 is just a test
    countabovehigh(ix) = sum(meanmeltH1daily{ix}(1:doy) > 0.0005);
    if isempty(onsetidx{ix}) ==0
        daysonset(ix) = onsetidx{ix}(1); %onset above 0.5 melt
        daysduration(ix) = max(onsetidx{ix}) - min(onsetidx{ix}); %duration above 0.5 melt
    else
        daysonset(ix) = NaN;
        daysduration(ix) = NaN;
    end
     sdmeltyr{ix} = std(meanmeltH1daily{ix});
     meanofmeanmelt{ix} = mean(meanmeltH1daily{ix});
     meltcum{ix} = cumsum(meanmeltH1daily{ix}); %total
     meltsum(ix) = sum(meanmeltH1daily{ix});
     
     meltcumbef{ix} = cumsum(meanmeltH1daily{ix}(1:doy));
     meltsumbef(ix) = sum(meanmeltH1daily{ix}(1:doy));
     meltvar{ix} = std(meanmeltH1daily{ix}(1:doy));
     [meltvarmax{ix}, meltvarmaxidx{ix}] = max(meltvar{ix});
     meanmeltH1annual{ix} = nanmean(meanmeltH1daily{ix});
     summeltH1annual{ix} = sum(meanmeltH1daily{ix});

     %%% rainfall
     for ib = 1:length(rowsubscript)
         rainH1daily{ix}(ib, 1, :) = rain{ix}(rowsubscript(ib), colsubscript(ib), :); %: x 1 x 1 x 365
     end
     meanrainH1daily{ix} = squeeze(mean(rainH1daily{ix}, 1));
     onsetidxrain{ix} = find(meanrainH1daily{ix}(1:doy) > 0.5); %0.5 is just a test
     countabovehighrain{ix} = sum(meanrainH1daily{ix}(1:doy) > 0.5);
     if isempty(onsetidxrain{ix}) ==0
         daysonsetrain{ix} = onsetidxrain{ix}(1); %onset above mean melt
         daysdurationrain{ix} = max(onsetidxrain{ix}) - min(onsetidxrain{ix}); %duration above mean melt
     else
         daysonsetrain{ix} = NaN;
         daysdurationrain{ix} = NaN;
     end
    raincum{ix} = cumsum(meanrainH1daily{ix}); 
    rainsum(ix) = sum(meanrainH1daily{ix});
    raincumbef{ix} = cumsum(meanrainH1daily{ix}(1:doy));
    rainsumbef(ix) = sum(meanrainH1daily{ix}(1:doy));
    
    %%% temperature
    for iz = 1:length(rowsubscript)
        tempH1daily{ix}(iz, 1, 1, :) = temp{ix}(rowsubscript(iz), colsubscript(iz), 1, :); %: x 1 x 1 x 365
    end
    meantempH1daily{ix} = squeeze(mean(tempH1daily{ix}, 1));
    meantempH1dailywinter{ix} = squeeze(mean(tempH1daily{ix}(1:90), 1)); %winter defined as Jan-March
    meantempH1annual(ix) = nanmean(meantempH1daily{ix});
    meantempH1annualwinter(ix) = nanmean(meantempH1dailywinter{ix});
    
    %%% for calculating the cold content
    for ic = 2:length(outlay{ix})
        volumes{ix}(1) = 0;
        volumes{ix}(ic) = outlay{ix}(ic) - outlay{ix}(ic-1); 
    end
    for i = 1:length(rowsubscript)
        tempiceH1daily{ix}(i, 1, :, :) = tempice{ix}(rowsubscript(i), colsubscript(i), :, :); %: x 1 x 1 x 365
        rhosH1daily{ix}(i, 1, :, :) = rhos{ix}(rowsubscript(i), colsubscript(i), :, :); %: x 1 x 1 x 365
    end
    meantempiceH1daily{ix} = squeeze(mean(tempiceH1daily{ix}));
    meanrhoH1daily{ix} = squeeze(mean(rhosH1daily{ix}));
    for il = 1:18
        mass{ix}(il, :) = volumes{ix}(il)*meanrhoH1daily{ix}(il, :); % calculate mass in each cell
    end
    weightsmass{ix} = mass{ix}./(sum(mass{ix})); %weights for later
    weightsvolume{ix} = volumes{ix}./(sum(volumes{ix}));
    meantempiceH1annual{ix} = nanmean(meantempiceH1daily{ix}, 2); %annual unweighted tempf and rhof
    meanrhoH1annual{ix} = nanmean(meanrhoH1daily{ix}, 2);
    
    meantempiceH1dailyweighted{ix} = meantempiceH1daily{ix}.*weightsmass{ix}; %daily weighted tempf and rhof - sum these
    meanrhoH1dailyweighted{ix} = meanrhoH1daily{ix}.*(weightsvolume{ix}'*ones(1, length(meanrhoH1daily{ix})));
    
    meantempiceH1annualweighted{ix} = sum(nanmean(meantempiceH1dailyweighted{ix}, 2)); %mean across days and mean across column (already weighted)
    meanrhoH1annualweighted{ix} = sum(nanmean(meanrhoH1dailyweighted{ix}, 2));
    
    coldcontentannual(ix) = cold_content(meanrhoH1annualweighted{ix}, meantempiceH1annualweighted{ix}, 20); %in J; annual
    latentheatcontentannual(ix) = latent_heat(summeltH1annual{ix}, meanrhoH1annualweighted{ix}); %in J; annual
end


%% 3 Save results
cd(fnsave);
save(fnout, 'meltsum', 'meltsumbef', 'daysonset', 'daysduration', ...
    'meltsmbratioH1', 'countabovehigh', 'sumwintersmb', 'meantempH1annual', ...
    'meantempH1annualwinter', 'rainsumbef', 'rainsum', 'coldcontentannual', ...
    'latentheatcontentannual', 'lonext', 'latext', 'xrefallk', 'yrefallk', ...
    'datesmar');

