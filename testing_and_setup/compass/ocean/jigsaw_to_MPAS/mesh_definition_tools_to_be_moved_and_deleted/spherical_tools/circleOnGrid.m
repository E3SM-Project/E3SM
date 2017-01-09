function cellWidthOut = circleOnGrid(lon, lat, centerLon, centerLat, radius, tanhWidth)
% circleOnGrid: combine two cell width distributions using a tanh function.
% This is intended as part of the workflow to make an MPAS global mesh.
%
% Syntax: cellWidthOut = circleOnGrid(lat, lon, cellWidthInAtlantic, cellWidthInPacific)
%
% Inputs:
%    lon - vector of length m, with entries between -180, 180, degrees
%    lat - vector of length n, with entries between -90, 90, degrees
%
% Optional inputs:
%
% Outputs:
%    cellWidthOut - m by n array, grid cell width on globe, km
%
% Example: 
%    RRS18to6 = circleOnGrid(lat,18,6)
%
% See also: 

% Author: Mark Petersen
% Los Alamos National Laboratory
% March 2018; Last revision: 3/27/2018

cellWidthOut = zeros(length(lon),length(lat));
for i=1:length(lon)
  for j=1:length(lat)
    [dist d2km]=lldistkm([centerLat, centerLon], [lat(j), lon(i)]);
    cellWidthOut(i,j) = 0.5*(-tanh((dist - radius)/tanhWidth) + 1.0);
  end
end
