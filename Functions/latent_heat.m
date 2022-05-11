function [lh] = latent_heat(melt, rhofirn)
%This function calculates the latent heat from a given meltwater input 
%balance within the firn of mean density rhofirn.
%Author: Annika Horlings
%Date created: 28 Feb 2022
%University of Washington

L = 333000; %latent heat of fusion, in J
rhow = 1000; %density of water
rhoi = 920; %ice density

lh = L*melt*rhow*((rhoi - rhofirn)/rhoi); %in J

end