function [cellWidthGlobal,lon,lat] = Arctic60to9_Atl10deg
% Create cell width for this mesh on a regular latitude-longitude grid.
% Outputs:
%    cellWidthGlobal - m x n array, entries are desired cell width in km
%    lon - longitude, vector of length m, with entries between -180 and 180, degrees
%    lat - latitude, vector of length n, with entries between -90 and 90, degrees

   ddeg = 1;
   lat = [ -90:ddeg: 90]';
   lon = [-180:ddeg:180]';

   EC60to30 = EC_CellWidthVsLat(lat);
   QU1 = ones(size(lat));
   
   AtlNH = mergeCellWidthVsLat(lat, 30*QU1, 9*QU1, 10, 10);
   AtlGrid = mergeCellWidthVsLat(lat, EC60to30, AtlNH, 0, 1);

   PacNH = mergeCellWidthVsLat(lat, 30*QU1, 9*QU1, 45, 10);
   PacGrid = mergeCellWidthVsLat(lat, EC60to30, PacNH, 0, 1);

	 cellWidthGlobal = AtlanticPacificGrid(lon, lat, AtlGrid, PacGrid);
