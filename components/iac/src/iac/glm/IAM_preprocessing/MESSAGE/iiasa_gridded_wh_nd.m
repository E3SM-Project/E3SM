clear all

wh2000be=importdata('Wood_Nov_25/FOHarvWoodEnergy2000A2.asc',' ',6);
wh2000be=wh2000be.data;
wh2010be=importdata('Wood_Nov_25/FOHarvWoodEnergy2010A2.asc',' ',6);
wh2010be=wh2010be.data;
wh2005be=(wh2000be+wh2010be)/2;
wh2100be=importdata('Wood_Nov_25/FOHarvWoodEnergy2100A2.asc',' ',6);
wh2100be=wh2100be.data;

wh2000sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2000A2.asc',' ',6);
wh2000sw=wh2000sw.data;
wh2010sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2010A2.asc',' ',6);
wh2010sw=wh2010sw.data;
wh2005sw=(wh2000sw+wh2010sw)/2;
wh2100sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2100A2.asc',' ',6);
wh2100sw=wh2100sw.data;

years=[2005,2100];

%interpolate to annual grids

eval(['wh_range(1,:,:)=(wh2005be+wh2005sw)*1.3;']);
eval(['wh_range(2,:,:)=(wh2100be+wh2100sw)*1.3;']);

annual=years(1):years(2);
annual_wh = interp1([years(1),years(2)],wh_range,annual);

clear wh_range

for k=1:(length(annual)-1)
    eval(['wh_grid' num2str(annual(k)) '= reshape(annual_wh(k,:,:),[360,720]);']);
    dlmwrite(['processed/nodata/gfwhd.', num2str(annual(k)), '.txt'], eval(['wh_grid',num2str(annual(k))]) ,'precision','%.0f','delimiter',' ');
end;

dlmwrite('processed/nodata/gfwhd.2100.txt',(wh2100sw+wh2100be)*1.3 ,'precision','%.0f','delimiter',' ');

