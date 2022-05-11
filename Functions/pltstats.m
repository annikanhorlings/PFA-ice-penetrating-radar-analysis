function [sddec, mndec, datesdec] = pltstats(var, datesmar)
%This function calculates the standard deviation and decadal mean of the
%MAR climate outputs for the paper for plotting.

%Author: Annika Horlings
%Date Created: 28 Feb 2022
%University of Washington
sddec(1) = std(var(1:2)); 
sddec(2) = std(var(3:12));
sddec(3) = std(var(13:22));
sddec(4) = std(var(23:32));
sddec(5) = std(var(33:42));
sddec(6) = std(var(43:52));
sddec(7) = std(var(53:62));
sddec(8) = std(var(63:69));

mndec(1) = nanmean(var(1:2)); 
mndec(2) = nanmean(var(3:12));
mndec(3) = nanmean(var(13:22));
mndec(4) = nanmean(var(23:32));
mndec(5) = nanmean(var(33:42));
mndec(6) = nanmean(var(43:52));
mndec(7) = nanmean(var(53:62));
mndec(8) = nanmean(var(63:69));

datesdec{1} = datesmar(1:2); 
datesdec{2} = datesmar(3:12);
datesdec{3} = datesmar(13:22);
datesdec{4} = datesmar(23:32);
datesdec{5} = datesmar(33:42);
datesdec{6} = datesmar(43:52);
datesdec{7} = datesmar(53:62);
datesdec{8} = datesmar(63:69);
end