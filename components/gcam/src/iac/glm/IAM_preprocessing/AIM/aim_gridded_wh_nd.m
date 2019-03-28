clear all

wh2005=importdata('6w\wood_harvested/wh_2005.txt',' ',6);
wh2005=wh2005.data;
wh2100=importdata('6w\wood_harvested/wh_2100.txt',' ',6);
wh2100=wh2100.data;

years=[2005,2100];

%interpolate to annual grids

eval(['wh_range(1,:,:)=wh2005*700*0.8*0.48*1.3/1000;']);
eval(['wh_range(2,:,:)=wh2100*700*0.8*0.48*1.3/1000;']);

annual=years(1):years(2);
annual_wh = interp1([years(1),years(2)],wh_range,annual);

clear wh_range

for k=1:(length(annual)-1)
    eval(['wh_grid' num2str(annual(k)) '= reshape(annual_wh(k,:,:),[360,720]);']);
    dlmwrite(['processed/nodata/gfwhd.', num2str(annual(k)), '.txt'], eval(['wh_grid',num2str(annual(k))]) ,'precision','%.0f','delimiter',' ');
end;

dlmwrite('processed/nodata/gfwhd.2100.txt',wh2100*700*0.8*0.48*1.3/1000 ,'precision','%.0f','delimiter',' ');

