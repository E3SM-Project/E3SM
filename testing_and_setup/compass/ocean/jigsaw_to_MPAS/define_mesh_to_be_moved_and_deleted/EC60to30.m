function [cellWidthGlobal,lon,lat] = EC60to30
% Create cell width for this mesh on a regular latitude-longitude grid.
% Outputs:
%    cellWidthGlobal - m x n array, entries are desired cell width in km
%    lon - longitude, vector of length m, with entries between -180 and 180, degrees
%    lat - latitude, vector of length n, with entries between -90 and 90, degrees

   ddeg = 1;
   lat = [ -90:ddeg: 90]';
   lon = [-180:ddeg:180]';
 
   EC60to30 = EC_CellWidthVsLat(lat);
   cellWidthGlobal = EC60to30*ones([1, length(lon)]);
