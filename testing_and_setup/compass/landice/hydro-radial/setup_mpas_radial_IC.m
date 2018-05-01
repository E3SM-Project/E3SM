
% calculate the stuff from Ed's code
[r,W,P,h,vb] = radialsteady(true)


fname = '/Users/mhoffman/documents/mpas-git/TESTS/landice/hydro-radial/1000m/steady_state_drift_test/run_model/landice_grid.nc';
xCell = ncread(fname, 'xCell');
yCell = ncread(fname, 'yCell');

rCell = (xCell.^2 + yCell.^2).^0.5;

h0 = interp1(r,W,rCell);
h0(isnan(h0))=0.0;

P0 = interp1(r,P,rCell);
P0(isnan(P0))=0.0;

figure(66); clf; hold all
subplot(2,1,1); hold all
plot(r, W)
plot(rCell, h0, 'r*')

subplot(2,1,2); hold all
plot(r, P)
plot(rCell, P0, 'r*')


% save the IC
ncwrite(fname, 'waterThickness', h0);
ncwrite(fname, 'waterPressure', P0);

