% script test

delta=1;
lat = [-90:delta:90]';
lon = [-180:delta:180]';
RRS30to10= RRS_CellWidthVsLat(lat,30,10);
RRS18to6 = RRS_CellWidthVsLat(lat,18,6);
EC60to30 = EC_CellWidthVsLat(lat);
cell_spac = EC60to30*ones(size(lon))';
size(cell_spac);
QU1 = ones(size(lat));

AtlNH = mergeCellWidthVsLat(lat, 30*QU1, 6*QU1, 30, 10);
AtlGrid = mergeCellWidthVsLat(lat, EC60to30, AtlNH, 0, 10);
PacNH = mergeCellWidthVsLat(lat, 30*QU1, 6*QU1, 60, 10);
PacGrid = mergeCellWidthVsLat(lat, EC60to30, PacNH, 0, 10);

%globalGrid = circleOnGrid(lon, lat, -50.0, -40.0, 5000, 1000);
%imagesc(lon', lat', globalGrid')
%set(gca,'Ydir','Normal')
%colorbar
%return

subplot(3,1,1)
plot(lat, EC60to30, lat, RRS18to6)
grid on
axis([-90 90 0 62])
xlabel('latitude, degrees')
ylabel('cell size, km')
legend('EC60to30','RRS18to6')

subplot(3,1,2)
plot(lat, AtlGrid, lat, PacGrid, lat, RRS18to6)
axis([-90 90 0 62])
xlabel('latitude, degrees')
ylabel('cell size, km')
legend('Atlantic','Pacific')
grid on

globalGrid = AtlanticPacificGrid(lon, lat, AtlGrid, PacGrid);
subplot(3,1,3)
imagesc(lon', lat', globalGrid')
title('Grid cell size, km')
set(gca,'Ydir','Normal')
xlabel('longitude, degrees')
ylabel('latitude, degrees')
colorbar
return
% outer product:
lon = -180:5:175;
lonOnes = ones(size(lon));
%RRS18to6LonLat = RRS18to6*lonOnes;
%imagesc(lon,y,RRS18to6LonLat)

