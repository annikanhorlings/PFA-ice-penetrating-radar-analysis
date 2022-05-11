%This script calculates the cold content & latent heat of total melt 
%production through a simple energy balance approximation, at 6 different 
%locations in southeast Greenland averaged over 3-4 nearest MAR cells.

%Author: Annika Horlings
%University of Washington
%Last updated: 1 February 2022

%Outputs:



%Dependencies:
% MAR climate reanlysis
% - MARv3.5.2 climate model simulation outputs and metadata are available 
%   from the NSF Arctic Data Center at:
%   https://arcticdata.io/catalog/view/doi:10.18739/A21Z41T0T.
% - reference coordinates 
    %which are produced here:


%% Load in
%Need total melt, temperature of firn, and density of firn
%comment out if running other sections, this is slow
cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Data/Climate/MAR');
mar_struct = dir('*.nc');
N = length(mar_struct);
for b = 1:N
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
datesmar = 1948:1:2016;    

%coordinates
for ik = 1:length(datesmar)
    [xcoord{ik}, ycoord{ik}] = polarstereo_fwd(latitude{ik}, longitude{ik}, 6378137.0,0.08181919,70,-45);
    xcoordk{ik} = xcoord{ik}/1000;
    ycoordk{ik} = ycoord{ik}/1000;
end

%resizing & taking the mean for each year
for il = 1:length(datesmar)
    meltmean{il} = mean(melt{il}(:, :, 1, :), 4); meltmean{il}(find(meltmean{il} > 9)) = NaN; %meltmean = meltmean.*(10^35);%calculate mean climate variable mean
    smbmean{il} = mean(smb{il}(:, :, 1, :), 4); %calculate mean climate variable mean
    tempmean{il} = mean(temp{il}(:, :, 1, :), 4); %calculate mean climate variable mean
end

% -----------------------------------------------
% position along reference coordinates - select which to average over
% % %1 seaward limit of aquifer location
% lonext = [-38.928 -39.867 -41.503];
% latext = [66.354 65.765 65.106];

%2 upstream aquifer location
lonext = [-39.3954 -40.4463 -41.8572];
latext = [66.3610 65.8776 65.2031];

% %3 more north site, inland extent of those aquifers
% lonext = [-34.374 -34.552 -34.762];
% latext = [68.031 67.620 67.507];

% %5 more west
% lonext = [-39.907 -40.946 -42.375];
% latext = [66.387 65.8878 65.194];

% %6 even more west
% lonext = [-40.407 -41.446 -42.875];
% latext = [66.3887 65.8878 65.194];

% %4 more south site, inland extent of those aquifers
% lonext = [-43.882 -43.773 -43.670 -43.691];
% latext = [62.767 62.584 62.301 62.086];
% --------------------------------------------------
for i = 1:length(lonext)
    [xrefall(i),yrefall(i)] = polarstereo_fwd(latext(i), lonext(i), 6378137.0,0.08181919,70,-45);
end
xrefallk = xrefall/1000;
yrefallk = yrefall/1000;

[m, n] = size(ycoordk{1});
for i = 1:length(xrefallk)
    [rowidxpre{i}, colidxpre{i}] = find(abs(xrefallk(i)*ones(m,n) - xcoordk{1}) < 15 & ...
        abs(yrefallk(i)*ones(m, n) - ycoordk{1}) < 15);
    
    if length(rowidxpre{i}) > 1
        rowidx(i) = rowidxpre{i}(1);
        colidx(i) = colidxpre{i}(1);
    elseif isempty(rowidxpre{i}) == 1
        rowidx(i) = NaN;
        colidx(i) = NaN;
    else
        rowidx(i) = rowidxpre{i};
        colidx(i) = colidxpre{i};
    end
end
rowidxnon= rowidx(~isnan(rowidx));
colidxnon= colidx(~isnan(colidx));

rowsubscript = ind2sub(size(ycoordk), rowidxnon);
colsubscript = ind2sub(size(ycoordk), colidxnon);

% Get date of when the radar data was taken
% cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Computed_vals/Helheim1/all');
% doyin = load('Helheim1_dayofyear.mat'); 
% doy = doyin.doy; date = doyin.dv; idxr = doyin.idxofrdgms;
doy = 200; %if no date or want to make consistent, just select one near July 19
%% Get limate metrics

%1 - melt to last winter's smb ratio
    for it = 1:length(datesmar) - 1
        smbmeanwinterbef{it + 1} = (mean(smb{it}(:, :, 1, 300:365), 4) ...
            + mean(smb{it + 1}(:, :, 1, 1:90), 4))/2;   %Oct-March
        meltsmbratio{it+1} = meltmean{it+1}./smbmeanwinterbef{it+1};
    end
    
for rf = 1:length(datesmar) - 1
    for i = 1:length(rowsubscript)
        wintersmbH1{rf + 1} = smbmeanwinterbef{rf+1}(rowsubscript(i), ...
            colsubscript(i));
    end
end

for ri = 1:length(datesmar) - 1
    for i = 1:length(rowsubscript)
        meltsmbratioH1{ri + 1} = meltsmbratio{ri+1}(rowsubscript(i), ...
            colsubscript(i));
    end
end

%2 & 3 - melt timing (e.g., early/late season melt) & duration
for ih = 1:length(datesmar)
     meanmeltH1daily{ih} = squeeze(mean(melt{ih}, 1));
     sdmeltyr{ih} = std(meanmeltH1daily{ih});
     meanofmeanmelt{ih} = mean(meanmeltH1daily{ih});
end
meantot = mean(cell2mat(meanofmeanmelt));
meansdtot = mean(cell2mat(sdmeltyr));
%rainfall 
% for ih = 1:length(datesmar)
%      meanrainH1daily{ih} = squeeze(mean(rain{ih}, 1));
%      sdrainyr{ih} = std(meanrainH1daily{ih});
%      meanofmeanrain{ih} = mean(meanrainH1daily{ih});
% end

for ix = 1:length(datesmar)
    for i = 1:length(rowsubscript)
        meltH1daily{ix}(i, 1, 1, :) = melt{ix}(rowsubscript(i), ...
            colsubscript(i), 1, :); %: x 1 x 1 x 365
    end
    meanmeltH1daily{ix} = squeeze(mean(meltH1daily{ix}, 1))/1000;
    onsetidx{ix} = find(meanmeltH1daily{ix}(1:doy) > 0.5); %rate of 0.5 mm per yr
    countabovehigh{ix} = sum(meanmeltH1daily{ix}(1:doy) > 0.5);
    
    if isempty(onsetidx{ix}) ==0
        daysonset{ix} = onsetidx{ix}(1); %onset above mean melt
        daysduration{ix} = max(onsetidx{ix}) - min(onsetidx{ix}); 
    else
        daysonset{ix} = NaN;
        daysduration{ix} = NaN;
    end
    
    meanmeltH1annual{ix} = nanmean(meanmeltH1daily{ix});
    summeltH1annual{ix} = sum(meanmeltH1daily{ix});
end

for ix = 1:length(datesmar)
    for i = 1:length(rowsubscript)
        rainH1daily{ix}(i, 1, :) = rain{ix}(rowsubscript(i), colsubscript(i), :); %: x 1 x 1 x 365
    end
%     meanmeltH1daily{ix} = squeeze(mean(meltH1daily{ix}, 1));
%     onsetidx{ix} = find(meanmeltH1daily{ix} > mean(meanmeltH1daily{ix})); %where the melt exceeds the mean melt for the entire year
%     daysonset{ix} = onsetidx{ix}(1); %onset above mean melt
%     daysduration{ix} = max(onsetidx{ix}) - min(onsetidx{ix}); %duration above mean melt

    meanrainH1daily{ix} = squeeze(mean(rainH1daily{ix}, 1))/1000;
    onsetidx{ix} = find(meanrainH1daily{ix}(1:doy) > 0.5); %0.5 is just a test
    countabovehigh{ix} = sum(meanrainH1daily{ix}(1:doy) > 0.5);
    
    if isempty(onsetidx{ix}) ==0
        daysonsetrain{ix} = onsetidx{ix}(1); %onset above mean melt
        daysdurationrain{ix} = max(onsetidx{ix}) - min(onsetidx{ix}); %duration above mean melt
    else
        daysonsetrain{ix} = NaN;
        daysdurationrain{ix} = NaN;
    end
end

for ix = 1:length(datesmar) %volume
    for ib = 2:length(outlay{ix})
        volumes{ix}(1) = 0;
        volumes{ix}(ib) = outlay{ix}(ib) - outlay{ix}(ib-1); %volume
    end
end

for ix = 1:length(datesmar) %ice temp & snow density % weighting
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
    
    %meantempiceH1dailyweighted{ix} = meantempiceH1daily{ix}.*weightsmass{ix}.*(weightsvolume{ix}'*ones(1, length(meanrhoH1daily{ix}))); %daily weighted tempf and rhof - sum these
    meantempiceH1dailyweighted{ix} = meantempiceH1daily{ix}.*weightsmass{ix}; %daily weighted tempf and rhof - sum these
    meanrhoH1dailyweighted{ix} = meanrhoH1daily{ix}.*(weightsvolume{ix}'*ones(1, length(meanrhoH1daily{ix}))); 
    
    meantempiceH1annualweighted{ix} = sum(nanmean(meantempiceH1dailyweighted{ix}, 2)); %mean across days and mean across column (already weighted)
    meanrhoH1annualweighted{ix} = sum(nanmean(meanrhoH1dailyweighted{ix}, 2));
    
    sumtempiceH1annualweighted{ix} = sum(sum(meantempiceH1dailyweighted{ix}, 2)); %sum temperatures for sum of cold content (????)
    
    %calculate cold content for upper 10 m (not 20 m, like above) just as
    %an experiment
    
%     weightsmass10{ix} = mass{ix}(1:16, :)./(sum(mass{ix}(1:16))); %weights for later
%     weightsvolume10{ix} = volumes{ix}(1:16)./(sum(volumes{ix}(1:16)));
%     
%     meantempiceH1dailyweighted10{ix} = meantempiceH1daily{ix}(1:16, :).*weightsmass10{ix}.*(weightsvolume10{ix}'*ones(1, length(meanrhoH1daily{ix}))); %daily weighted tempf and rhof - sum these
%     meanrhoH1dailyweighted10{ix} = meanrhoH1daily{ix}(1:16, :).*(weightsvolume10{ix}'*ones(1, length(meanrhoH1daily{ix}))); 
%     
%     meantempiceH1annualweighted10{ix} = sum(nanmean(meantempiceH1dailyweighted10{ix}, 2)); %mean across days and mean across column (already weighted)
%     meanrhoH1annualweighted10{ix} = sum(nanmean(meanrhoH1dailyweighted10{ix}, 2));
%     
%     sumtempiceH1annualweighted10{ix} = sum(sum(meantempiceH1dailyweighted10{ix}, 2)); %sum temperatures for sum of cold content (????)
end

%4 - sum & cumulative melt
for ip = 1:length(datesmar)
    meltcum{ip} = cumsum(meanmeltH1daily{ip}); %total
    meltsum{ip} = sum(meanmeltH1daily{ip});

    meltcumbef{ip} = cumsum(meanmeltH1daily{ip}(1:doy));
    meltsumbef{ip} = sum(meanmeltH1daily{ip}(1:doy));

end

for ip = 1:length(datesmar)
    raincum{ip} = cumsum(meanrainH1daily{ip}); %total
    rainsum{ip} = sum(meanrainH1daily{ip});

    raincumbef{ip} = cumsum(meanrainH1daily{ip}(1:doy));
    rainsumbef{ip} = sum(meanrainH1daily{ip}(1:doy));

end

% sum of rainfall
for ip = 1:length(datesmar)
    raincum{ip} = cumsum(meanrainH1daily{ip}); %total
    rainsum{ip} = sum(meanrainH1daily{ip});

    raincumbef{ip} = cumsum(meanrainH1daily{ip}(1:doy));
    rainsumbef{ip} = sum(meanrainH1daily{ip}(1:doy));

end

meltsummat = cell2mat(meltsum);
for ip = 4:length(idxr) - 1 %starting at the one after 2001 (ie 2003) - check dates!!!!
 %melt between last radargram and time radargram was taken
    meltsumbeftot{idxr(ip) - 8} = meltsumbef{idxr(ip) - 10} + ...
        sum(meltsummat(idxr(ip) - 10: idxr(ip) - 10));
end


%5 - sum & cumulative wintertime smb
    for ig = 1:length(datesmar) - 1
        cumwintersmb{ig + 1} = cumsum(wintersmbH1{ig + 1});
        sumwintersmb{ig + 1} = sum(wintersmbH1{ig + 1});  %Oct-March
    end
    
%6 - melt variability (std of melt) of daily & day of max value
    for iy = 1:length(datesmar)
        meltvar{iy} = std(meanmeltH1daily{iy}(1:doy));
        [meltvarmax{iy}, meltvarmaxidx{iy}] = max(meltvar{iy});
    end

%7 - melt anomaly
for jf = 1:length(datesmar) - 7
    meltsumearly(jf, :) = meltsum{jf};
end
meltsummean01to09 = mean(meltsumearly, 1);
for jf = 1:length(datesmar) - 7
    meltsuml{jf} = meltsum{jf} - meltsummean01to09;
end

% mean temperature
for ix = 1:length(datesmar)
    for i = 1:length(rowsubscript)
        tempH1daily{ix}(i, 1, 1, :) = temp{ix}(rowsubscript(i), colsubscript(i), 1, :); %: x 1 x 1 x 365
    end
    meantempH1daily{ix} = squeeze(mean(tempH1daily{ix}(150:240), 1)); %summer defined as June-Aug
    meantempH1dailywinter{ix} = squeeze(mean(tempH1daily{ix}(1:61), 1)); %winter defined as Jan-Feb
    
    %annual mean summer temp
    meantempH1annual{ix} = nanmean(meantempH1daily{ix});
    meantempH1annualwinter{ix} = nanmean(meantempH1dailywinter{ix});
    
end
%% Compute cold content and latent heat
%these are the vars of interest:
m = meanmeltH1daily; %daily melt, in m of WEQ
L = 333000; %latent heat of fusion, in J
rhow = 1000; %density of water
c = 2097; %heat capacity of ice, in J/kgK 
Tf = meantempiceH1daily; %absolute value of the mass-weighted mean firn/ice temperature, in C (need to see what depth this averages over!!!!)
rhof = 650; %the volume-weighted mean firn density, in kg (see if we can figure this as a transient out more specifically)
rhoi = 920; %ice density

%cold content daily
for ih = 1:length(Tf)
    coldcontent{ih} = 20*c*meanrhoH1dailyweighted{ih}.*meantempiceH1dailyweighted{ih}; %in J
end

%cold content annual
for ih = 1:length(Tf)
    coldcontentannual{ih} = 20*c*meanrhoH1annualweighted{ih}.*meantempiceH1annualweighted{ih}; %in J
end

%latent heat content daily
for im = 1:length(m)
    latentheatcontent{im} = L*m{im}*rhow;%in J
end

%mean latent heat content annual
for im = 1:length(m)
    latentheatcontentannual{im} = L* summeltH1annual{im}*rhow*((rhoi - meanrhoH1annualweighted{im})/rhoi);%in J; includes saturating water %
end

%For now, do this for the days:
datesdaymar = linspace(1976, 2017, 366*41);

%% Plot
cmap = intense(max(10));

figure;
subplot(3, 1, 1);
for io = 1:length(coldcontentannual)
    
     coldcontentplot(:, io) = coldcontentannual{io};
%     for ip = 1:length(coldcontent{io})
%         plot(datesdaymar(ip*(datesmar(io)- 1975)), coldcontent{io}(ip), 'k.');
%         hold on
%     end
end
    
plot(datesmar,  abs(coldcontentplot)/(10^6), '-o', 'Color', cmap(3, :), 'MarkerFaceColor',cmap(3, :));
xlabel('date');
ylabel('cold content (MJ)'); %*above mean for year 
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 14);
xlim([1947.5 2016.5]);
grid on;
set(gca, 'Ydir', 'Normal');

subplot(3, 1, 2);
for io = 1:length(latentheatcontentannual)
    
    latentheatcontentplot(io) = latentheatcontentannual{io};

%     for ip = 1:length(coldcontent{io})
%         plot(datesdaymar(ip*(datesmar(io)- 1975)), coldcontent{io}(ip), 'k.');
%         hold on
%     end
end
plot(datesmar, latentheatcontentplot/(10^6), '-o', 'Color', cmap(4, :), 'MarkerFaceColor',cmap(4, :));
xlabel('date');
ylabel('latent heat content (MJ)'); %*above mean for year 
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 14);
xlim([1947.5 2016.5]);
grid on;

subplot(3, 1, 3);
for io = 1:length(latentheatcontentannual)
    
    LCCCratio(io) = latentheatcontentannual{io}./abs(coldcontentannual{io});
end
plot(datesmar, LCCCratio, '-o', 'Color', cmap(4, :), 'MarkerFaceColor',cmap(5, :));
xlabel('date');
ylabel('LH/CC ratio'); %*above mean for year 
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 14);
xlim([1947.5 2016.5]);
grid on;

