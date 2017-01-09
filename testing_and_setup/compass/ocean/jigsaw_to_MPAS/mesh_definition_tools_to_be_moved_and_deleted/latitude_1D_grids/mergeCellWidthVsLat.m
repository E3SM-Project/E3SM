function cellWidthOut = mergeCellWidthVsLat(latitude, cellWidthInSouth, cellWidthInNorth, varargin)
% mergeCellWidthVsLat: combine two cell width distributions using a tanh function.
% This is intended as part of the workflow to make an MPAS global mesh.
%
% Syntax: cellWidthOut = mergeCellWidthVsLat(latitude, cellWidthInSouth, cellWidthInNorth, latTransition, latWidthTransition)
%
% Inputs:
%    latitude - vector of length n, with entries between -90 and 90, degrees
%    cellWidthInSouth - vector of length n, first distribution
%    cellWidthInNorth - vector of length n, second distribution
%
% Optional inputs:
%    latTransition = 0; % latitude to change from cellWidthInSouth to cellWidthInNorth, degrees
%    latWidthTransition = 0; % width of latitude transition, degrees
%
% Outputs:
%    cellWidthOut - vector of length n, entries are cell width as a function of latitude

% Author: Mark Petersen
% Los Alamos National Laboratory
% March 2018; Last revision: 4/20/2018

% Assign defaults
latTransition = 0; % latitude to change from cellWidthInSouth to cellWidthInNorth, degrees
latWidthTransition = 0; % width of latitude transition, degrees

try
   latTransition = varargin{1};
   latWidthTransition = varargin{2};
end

if latWidthTransition == 0
  for j=1:length(latitude)
    if latitude(j)<latTransition
      cellWidthOut(j) = cellWidthInSouth(j);
    else
      cellWidthOut(j) = cellWidthInNorth(j);
    end
  end
else
  for j=1:length(latitude)
    weightNorth = 0.5*(tanh((latitude(j) - latTransition)/latWidthTransition) + 1.0);
    weightSouth = 1.0 - weightNorth;
    cellWidthOut(j) = weightSouth*cellWidthInSouth(j) + weightNorth*cellWidthInNorth(j);
  end
end
