clear all

wh2000be=importdata('Wood_Nov_25/FOHarvWoodEnergy2000A2.asc',' ',6);
wh2000be=wh2000be.data;
wh2010be=importdata('Wood_Nov_25/FOHarvWoodEnergy2010A2.asc',' ',6);
wh2010be=wh2010be.data;
wh2005be=(wh2000be+wh2010be)/2;
wh2020be=importdata('Wood_Nov_25/FOHarvWoodEnergy2020A2.asc',' ',6);
wh2020be=wh2020be.data;
wh2030be=importdata('Wood_Nov_25/FOHarvWoodEnergy2030A2.asc',' ',6);
wh2030be=wh2030be.data;
wh2040be=importdata('Wood_Nov_25/FOHarvWoodEnergy2040A2.asc',' ',6);
wh2040be=wh2040be.data;
wh2050be=importdata('Wood_Nov_25/FOHarvWoodEnergy2050A2.asc',' ',6);
wh2050be=wh2050be.data;
wh2060be=importdata('Wood_Nov_25/FOHarvWoodEnergy2060A2.asc',' ',6);
wh2060be=wh2060be.data;
wh2070be=importdata('Wood_Nov_25/FOHarvWoodEnergy2070A2.asc',' ',6);
wh2070be=wh2070be.data;
wh2080be=importdata('Wood_Nov_25/FOHarvWoodEnergy2080A2.asc',' ',6);
wh2080be=wh2080be.data;
wh2090be=importdata('Wood_Nov_25/FOHarvWoodEnergy2090A2.asc',' ',6);
wh2090be=wh2090be.data;
wh2100be=importdata('Wood_Nov_25/FOHarvWoodEnergy2100A2.asc',' ',6);
wh2100be=wh2100be.data;

wh2000sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2000A2.asc',' ',6);
wh2000sw=wh2000sw.data;
wh2010sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2010A2.asc',' ',6);
wh2010sw=wh2010sw.data;
wh2005sw=(wh2000sw+wh2010sw)/2;
wh2020sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2020A2.asc',' ',6);
wh2020sw=wh2020sw.data;
wh2030sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2030A2.asc',' ',6);
wh2030sw=wh2030sw.data;
wh2040sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2040A2.asc',' ',6);
wh2040sw=wh2040sw.data;
wh2050sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2050A2.asc',' ',6);
wh2050sw=wh2050sw.data;
wh2060sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2060A2.asc',' ',6);
wh2060sw=wh2060sw.data;
wh2070sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2070A2.asc',' ',6);
wh2070sw=wh2070sw.data;
wh2080sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2080A2.asc',' ',6);
wh2080sw=wh2080sw.data;
wh2090sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2090A2.asc',' ',6);
wh2090sw=wh2090sw.data;
wh2100sw=importdata('Wood_Nov_25/FOHarvWoodIndusty2100A2.asc',' ',6);
wh2100sw=wh2100sw.data;

years=[2005,2010:10:2100];

for ind=1:(length(years)-1)
   %interpolate to annual grids
    years(ind)
   
    eval(['wh_range(1,:,:)=(wh',num2str(years(ind)),'be+wh',num2str(years(ind)),'sw)*1.3;']);
    eval(['wh_range(2,:,:)=(wh',num2str(years(ind+1)),'be+wh',num2str(years(ind+1)),'sw)*1.3;']);
    
    annual=years(ind):years(ind+1);
    annual_wh = interp1([years(ind),years(ind+1)],wh_range,annual);

    clear wh_range

    for k=1:(length(annual)-1)
        eval(['wh_grid' num2str(annual(k)) '= reshape(annual_wh(k,:,:),[360,720]);']);
        dlmwrite(['processed/wh/gfwhd.', num2str(annual(k)), '.txt'], eval(['wh_grid',num2str(annual(k))]) ,'precision','%.0f','delimiter',' ');
    end;

end

dlmwrite('processed/wh/gfwhd.2100.txt',(wh2100sw+wh2100be)*1.3 ,'precision','%.0f','delimiter',' ');

