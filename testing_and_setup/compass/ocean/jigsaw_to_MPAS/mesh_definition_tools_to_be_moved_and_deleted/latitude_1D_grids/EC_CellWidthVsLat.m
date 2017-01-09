function cellWidthOut = EC_CellWidthVsLat(latitude, varargin)
% EC_CellWidthVsLat - Create Eddy Closure spacing as a function of latitude.
% This is intended as part of the workflow to make an MPAS global mesh.
%
% Syntax: cellWidthOut = EC_CellWidthVsLat(latitude, cellWidthEq, cellWidthMidLat, cellWidthPole,
%                                          latPosEq, latPosPole, latTransition, 
%                                          latWidthEq, latWidthPole)
% Inputs:
%    latitude - vector of length n, with entries between -90 and 90, degrees
%
% Optional inputs:
%    % Default values for Cell width, km
%    cellWidthEq = 30.0; % Eq is equatorial latitude
%    cellWidthMidLat = 60.0; % MidLat is mid latitude
%    cellWidthPole = 35.0; % Pole is polar latitude
%    
%    % Default values for latitude positions in degrees
%    latPosEq = 15.0; % position of center of transition region
%    latPosPole = 73.0; % position of center of transition region
%    latTransition = 40; % latitude to change from Eq to Pole function
%    latWidthEq = 6.0; % width of transition region
%    latWidthPole = 9.0; % width of transition region
%    
% Outputs:
%    cellWidthOut - vector of length n, entrie are cell width as a function of latitude
%
% Example: 
%    EC60to30 = EC_CellWidthVsLat(latitude)
%    EC120to60 = EC_CellWidthVsLat(latitude,60,120,70)

% Author: Mark Petersen
% Los Alamos National Laboratory
% March 2018; Last revision: 4/20/2018

% Default values for Cell width, km
cellWidthEq = 30.0; % Eq is equatorial latitude
cellWidthMidLat = 60.0; % MidLat is mid latitude
cellWidthPole = 35.0; % Pole is polar latitude

% Default values for latitude positions in degrees
latPosEq = 15.0; % position of center of transition region
latPosPole = 73.0; % position of center of transition region
latTransition = 40; % latitude to change from Eq to Pole function
latWidthEq = 6.0; % width of transition region
latWidthPole = 9.0; % width of transition region

try
  cellWidthEq =     varargin{1}; 
  cellWidthMidLat = varargin{2}; 
  cellWidthPole =   varargin{3}; 
  latPosEq =        varargin{4}; 
  latPosPole =      varargin{5}; 
  latTransition =   varargin{6}; 
  latWidthEq =      varargin{7}; 
  latWidthPole =    varargin{8}; 
end

degToRad = pi/180.0; % convert degrees to radians
minCellWidth = min(cellWidthEq, min(cellWidthMidLat, cellWidthPole));
densityEq = (minCellWidth/cellWidthEq)^4;
densityMidLat = (minCellWidth/cellWidthMidLat)^4;
densityPole = (minCellWidth/cellWidthPole)^4;
densityEC = zeros(size(latitude));
cellWidthOut = zeros(size(latitude));
for j=1:length(latitude)
  if abs(latitude(j))<latTransition
    densityEC(j) = ((densityEq-densityMidLat) * (1.0 + tanh( (latPosEq - abs(latitude(j)))/ latWidthEq)) / 2.0) + densityMidLat;
  else
    densityEC(j) = ((densityMidLat-densityPole) * (1.0 + tanh( (latPosPole - abs(latitude(j)))/ latWidthPole)) / 2.0) + densityPole;
  end
  cellWidthOut(j) = minCellWidth/densityEC(j)^0.25;
end
