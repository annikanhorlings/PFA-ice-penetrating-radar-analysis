function [cc] = cold_content(rho, tempfirn, thickness)
%This function calculates the cold content given the inputs of firn density
%and firn temperature for a given firn thickness.
%Author: Annika Horlings
%Date created: 28 Feb 2022
%University of Washington

c = 2097; %heat capacity of ice, in J/kgK 
cc = thickness*c*rho*tempfirn; %in J

end