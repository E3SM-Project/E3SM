function p = params(E0,Y0);
% PARAMS  Return parameters in a structure for subglacial hydrology model.
% Example:  View default values
%   >> params()
% CODE WRITTEN BY ED BUELER: https://github.com/bueler/hydrolakes/tree/master/codes

%p.spera = 31556926.0;
p.spera = 3600.0*24.0*365.0;
p.rhoi  = 910.0;         % kg m-3
p.rhow  = 1000.0;        % kg m-3
%p.g     = 9.81;          % m s-2
p.g     = 9.80616; 


% major model parameters:
p.A  = 3.1689e-24;       % ice softness (Pa-3 s-1)
p.K  = 1.0e-2;           % m s-1   FIXME: want Kmax or Kmin according to W > Wr
p.Wr = 1.0;              % m
p.c1 = 0.500;            % m-1
p.c2 = 0.040;            % [pure]
p.c0 = p.K / (p.rhow * p.g);   % constant in velocity formula

% regularization of "W=Y" pressure closure
if nargin < 1
  p.E0 = 1.0;            % m; what is optimal?
else
  p.E0 = E0;
end

% regularizations of closing term
if nargin < 2
  p.Y0 = 0.001;            % m
else
  p.Y0 = Y0;
end

