function cellWidthOut = AtlanticPacificGrid(lon, lat, cellWidthInAtlantic, cellWidthInPacific)
% AtlanticPacificGrid: combine two cell width distributions using a tanh function.
% This is intended as part of the workflow to make an MPAS global mesh.
%
% Syntax: cellWidthOut = AtlanticPacificGrid(lat, lon, cellWidthInAtlantic, cellWidthInPacific)
%
% Inputs:
%    lon - vector of length m, with entries between -180, 180, degrees
%    lat - vector of length n, with entries between -90, 90, degrees
%    cellWidthInAtlantic - vector of length n, cell width in Atlantic as a function of longitude, km
%    cellWidthInPacific - vector of length n, cell width in Pacific as a function of longitude, km
%
% Optional inputs:
%
% Outputs:
%    cellWidthOut - m by n array, grid cell width on globe, km
%
% Example: 
%    RRS18to6 = RRS_CellWidthVsLat(lat,18,6)
%
% See also: 

% Author: Mark Petersen
% Los Alamos National Laboratory
% March 2018; Last revision: 3/27/2018

cellWidthOut = zeros(length(lon),length(lat));
for i=1:length(lon)
  for j=1:length(lat)
    % set to Pacific mask as default
    cellWidthOut(i,j) = cellWidthInPacific(j);
    % test if in Atlantic Basin:
    if lat(j)>65.0
      if and(lon(i)>-150.0, lon(i)<170.0)
        cellWidthOut(i,j) = cellWidthInAtlantic(j);
      end
    elseif lat(j)>20.0
      if and(lon(i)>-100.0, lon(i)<35.0)
        cellWidthOut(i,j) = cellWidthInAtlantic(j);
      end
    elseif lat(j)>0.0
      if and(lon(i)>-2*lat(j)-60.0, lon(i)<35.0)
        cellWidthOut(i,j) = cellWidthInAtlantic(j);
      end
    else
      if and(lon(i)>-60.0, lon(i)<20.0)
       cellWidthOut(i,j) = cellWidthInAtlantic(j);
      end
    end
  end
end
cellWidthOut = cellWidthOut';
