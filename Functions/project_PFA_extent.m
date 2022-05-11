function [xcoordinate, ycoordinate, xPFA_inrange, yPFA_inrange, extentPFA_fromline, ...
     distance_fromline, idxPFA_inrange, whereidxofPFAinrange, picks_inrange,  ...
    xcoordref,  ycoordref] = project_PFA_extent(data, ...
    rdrtype, number, Cmnx, Cmxx, Cmny, Cmxy)
    %This function projects a vector of the PFA extent onto desired
    %coordinates. 

    %INPUT VARIABLES:
    %data: structure of all CRESIS profiles processed and picked by ImpDAR, 
        %exceeding 4 or more profiles total
    %rdrtype: either 0 for RDS or 1 or AR
    %number: number of profiles
    %Cmnx: min x value for perpendular line to reference extent in km
    %Cmnx: max x value for perpendular line to reference extent in km
    %Cmny: min y value for perpendular line to reference extent in km
    %Cmxy: max y value for perpendular line to reference extent in km
    
    %OUTPUT VARIABLES:
    %xPFA_inrange: x coordinates within the common range
    %yPFA_inrange: y coordinates within the common range
    %xb: x coordinates of the "base vector" onto which to project
    %yb: y coordinates of the "base vector" onto which to project
    %project_xPFAbase: x coordinates of the PFA that are projected on f(xb)
    %project_yPFAbase: y coordinates of the PFA that are projected on f(xb)
    %distancec: distance along flight lines of the identified PFA
    %extentPFA_fromline: distance from a line perpendicular to all flight
        %lines of each of the identified PFAs (i.e., the extent of PFA)
    %whereidxofPFAinrange: index of the PFA, with regard to the coordinate of 
        %PFA_inrange
        
%     EXAMPLE
%     cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Data/Helheim1/AR');
%     profiles_struct = dir('*.mat');
%     profiles = struct2cell(profiles_struct);
%     [m, n] = size(profiles);
%     for i = 1:n
%        data{i} = load(profiles{1, i});
%     end
%     cd('/Users/annikahorlings/Documents/Projects/Greenland_Firn_Aquifer/Finalized_material_for_paper/Computed_vals/Helheim1/all');
%       plv = load('perpendicular_line_vals');
%       Cmnx = plv.Cmnx/1000;
%       Cmxx = plv.Cmxx/1000;
%       Cmny = plv.Cmny/1000;
%       Cmxy = plv.Cmxy/1000;
%     rdrtype = 1;
%     number = 5;
% [xcoordinate, ycoordinate, xPFA_inrange, yPFA_inrange, extentPFA_fromline, ...
%      distance_fromline, idxPFA_inrange, whereidxofPFAinrange, picks_inrange, ...
%      xcoordref,  ycoordref] = project_PFA_extent(data, ...
%     rdrtype, number, Cmnx, Cmxx, Cmny, Cmxy);
% ---------------------------------------------------------------------------------------------------
    %FUNCTION:
    %load variables from the data
    for jl = 1:number
        latPFA{jl} = data{jl}.lat; %latitude
        lonPFA{jl} = data{jl}.long; %longitude
        picks{jl} = data{jl}.picks.samp2; %picks
    end

    %Convert the lat and lon to x and y coordinates in polar stereographic
    %projection
    for hj = 1:length(latPFA)
        [xPFA{hj}, yPFA{hj}] = polarstereo_fwd(latPFA{hj}, lonPFA{hj}, 6378137.0,0.08181919,70,-45);
        xcoordinate{hj} = xPFA{hj}/1000;
        ycoordinate{hj} = yPFA{hj}/1000;
    end
    
    %Find common coordinates of all the x and y coordinates, where:
        %x_on: the x coordinates on which to interpolate
        %y_on: the y coordinates on which to interpolate
    
    lx = length(xPFA);   %take the length of the cell
    for i = 1:lx %reorganize cells to mat
        lxi = length(xPFA{i});
        lyi = length(yPFA{i});
        xn(i, 1:lxi) = xPFA{i}/1000;
        yn(i, 1:lyi) = yPFA{i}/1000;
    end
    xn(xn==0)=NaN; %replace 0s that were filler with NaN
    yn(yn==0)=NaN;
    %find the common values for the first few
    
%     if number == 1
%         disp('Error: cannot compare common range of values');
%     end
%     
%     if number < 3 && number > 1
%         Cx = mintersect(round(xn(1, :), 1), round(xn(2, :), 1));
%         Cy = mintersect(round(yn(1, :), 1), round(yn(2, :), 1));
%         Cmxx = max(Cx); Cmnx = min(Cx); %find the min/max of the common range of values
%         Cmxy = max(Cy); Cmny = min(Cy);
%     end
%     
%     if number > 3
%         Cx = mintersect(round(xn(1, :), 1), round(xn(2, :), 1), round(xn(3, :), 1), round(xn(4, :), 1));
%         Cy = mintersect(round(yn(1, :), 1), round(yn(2, :), 1), round(yn(3, :), 1), round(yn(4, :), 1));
%         Cmxx = max(Cx); Cmnx = min(Cx); %find the min/max of the common range of values
%         Cmxy = max(Cy); Cmny = min(Cy);
%     end
    
    
    for k = 1:length(xPFA) %calculate indices that are shared & corresponding x/y coords
        idxPFA_inrange{k} = find(xn(k, :) < Cmxx & xn(k, :) > Cmnx & yn(k, :) < Cmxy & yn(k, :) > Cmny); %find the vals in the range
        lidxi = length(idxPFA_inrange{k});
        xPFA_inrange{k} = xn(k, idxPFA_inrange{k});
        yPFA_inrange{k} = yn(k, idxPFA_inrange{k});
        
        picks_inrange{k} = picks{k}(idxPFA_inrange{k}); %take the picks within that range
        
     end
%     
%     %Interpolate a line or function through the x and y coords of the PFA
%     for ik = 1:length(xPFA)
%         pPFA{ik} = polyfit(xPFA_inrange{ik}, yPFA_inrange{ik}, 2);
%         fPFA{ik} = polyval(pPFA{ik},xPFA_inrange{ik});
%         coordallPFA{ik} = vertcat(xPFA_inrange{ik}, fPFA{ik});
%     end

    %Determine a line orthogonal to the interior-seaward dimension to
    %calculate inland extent
%     
%      %Interpolate a line or function through the x and y coords that will
%      %be the basis of projecting everything onto
%      xb = linspace(Cmnx, Cmxx, 500);
%      yb = linspace(Cmny, Cmxy, 500);
%      pb = polyfit(xb, yb, 2);
%      fb = polyval(pb,xb);
%      coordallb = vertcat(xb, fb);
%      
%     %identify where the aquifer is, coordinate-wise
%     for ig = 1:length(xPFA)
%         if rdrtype == 1
%             coordallPFAextent{ig} = NaN(2, length(xPFA_inrange{ig}));
%             coordallPFAextent{ig}(:, find(~isnan(picks_inrange{ig}))) ...
%                 = coordallPFA{ig}(:, find(~isnan(picks_inrange{ig}))); %extent for AR; not the NaN
%         end
%         
%         if rdrtype == 0
%             coordallPFAextent{ig} = coordallPFA{ig};
%             coordallPFAextent{ig} = ...
%                 coordallPFA{ig}(:, find(isnan(picks_inrange{ig}))); %extent for RDS; NaN vals
%         end
%     end
%      
%     %project onto the basis coords for each PFA
%     for ij = 1:length(xPFA)
%         for ih = 1:length(coordallPFAextent{ij})
%             proj_PFA_pts{ij}(:, ih) = proj(coordallb, coordallPFAextent{ij}(:, ih)');
%         end
%         
%         proj_xPFAbase{ij} = proj_PFA_pts{ij}(1, :);
%         proj_yPFAbase{ij} = proj_PFA_pts{ij}(2, :);
%     end
%     
%     %calculate the distance along each flight line
%     for ix = 1:length(xPFA)
%             distance{ix}(1) = 0;
%         for iu = 2:length(coordallPFA{ix}(1,:))-1
%             distance{ix}(iu) = sqrt((coordallPFA{ix}(1,iu) - coordallPFA{ix}(1,iu-1)).^2 + (coordallPFA{ix}(2,iu) - coordallPFA{ix}(2,iu)).^2) + distance{ix}(iu-1);
%         end
%     end
    
   %calculate extent, from a line perpendicular to the main flight lines
    v1 = [Cmnx,Cmxy,0]; %a line going perpendicular to the main flight lines - fix!!!!!!!!!!!!
    v2 = [Cmnx,Cmny,0];
    
    for ic = 1:length(xPFA) 
        for jc = 1:length(xPFA_inrange{ic})
            distance_fromline{ic}(jc) = point_to_line([xPFA_inrange{ic}(jc), ...
                yPFA_inrange{ic}(jc), 0], v1, v2);
        end
        
        if rdrtype == 1 %AR
            isno{ic} = ~isnan(picks_inrange{ic}); %where there are not NaNs
            extentPFA_fromline{1, ic} = distance_fromline{ic}.*(isno{ic}); %extent of PFA
            extentPFA_fromline{1, ic}(extentPFA_fromline{ic} == 0) = NaN;
            whereidxofPFAinrange{ic} = isno{ic};
        end
        
        if rdrtype == 0 %RDS
            knan{ic} = find(isnan(picks_inrange{ic})); %where there are NaNs
            kdiff{ic} = find(diff(knan{ic}) > 1);
            kloc{ic} = zeros(1, length(picks_inrange{ic}));
            kloc{ic}(knan{ic}) = 1;
            for j = 1:length(kloc{ic}) - 5
                    if (kloc{ic}(1, j) + kloc{ic}(1, j+1) + kloc{ic}(1, j+2) + kloc{ic}(1, j+3) + kloc{ic}(1, j+4) < 5)
                        kloc{ic}(j) = 0;
                    end
            end
            kloc{ic}(end-5:end) = 0;
            kloc{ic}(kloc{ic} == 0) = NaN;
            extentPFA_fromline{1, ic} = distance_fromline{ic}.*(kloc{ic}); %extent of PFA   
            extentPFA_fromline{1, ic}(extentPFA_fromline{ic} == 0) = NaN;
            whereidxofPFAinrange{ic} = kloc{ic};
        
        end
    end
    
    %define a line perpendcular to reference line onto which to
    %interpolate things onto in the future
    xcoordref = linspace(Cmnx, Cmxx, 1500); %you can choose the length
    ycoordref = linspace(Cmxy, Cmny, 1500);
    
%     double check with plot
%     figure;
%     for i = 1:5
%         plot(xPFA_inrange{i}, yPFA_inrange{i});
%         hold on
%     end
%     hold on
%     plot(v2(1), v2(2), 'k*')
%     hold on
%     plot(v1(1), v1(2), 'k*')
%     hold on;
%     plot(xcoordref, ycoordref, 'k--');

   
end




