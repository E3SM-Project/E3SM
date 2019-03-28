clear all

wh2005=importdata('6w\wood_harvested/wh_2005.txt',' ',6);
wh2005=wh2005.data;
wh2010=importdata('6w\wood_harvested\wh_2010.txt',' ',6);
wh2010=wh2010.data;
wh2020=importdata('6w\wood_harvested/wh_2020.txt',' ',6);
wh2020=wh2020.data;
wh2030=importdata('6w\wood_harvested/wh_2030.txt',' ',6);
wh2030=wh2030.data;
wh2040=importdata('6w\wood_harvested/wh_2040.txt',' ',6);
wh2040=wh2040.data;
wh2050=importdata('6w\wood_harvested/wh_2050.txt',' ',6);
wh2050=wh2050.data;
wh2060=importdata('6w\wood_harvested/wh_2060.txt',' ',6);
wh2060=wh2060.data;
wh2070=importdata('6w\wood_harvested/wh_2070.txt',' ',6);
wh2070=wh2070.data;
wh2080=importdata('6w\wood_harvested/wh_2080.txt',' ',6);
wh2080=wh2080.data;
wh2090=importdata('6w\wood_harvested/wh_2090.txt',' ',6);
wh2090=wh2090.data;
wh2100=importdata('6w\wood_harvested/wh_2100.txt',' ',6);
wh2100=wh2100.data;

years=[2005,2010:10:2100];

for ind=1:(length(years)-1)
   %interpolate to annual grids
    years(ind)
   
    eval(['wh_range(1,:,:)=wh',num2str(years(ind)),'*700*0.8*0.48*1.3/1000;']);
    eval(['wh_range(2,:,:)=wh',num2str(years(ind+1)),'*700*0.8*0.48*1.3/1000;']);
    
    annual=years(ind):years(ind+1);
    annual_wh = interp1([years(ind),years(ind+1)],wh_range,annual);

    clear wh_range

    for k=1:(length(annual)-1)
        eval(['wh_grid' num2str(annual(k)) '= reshape(annual_wh(k,:,:),[360,720]);']);
        dlmwrite(['processed/gfwhd.', num2str(annual(k)), '.txt'], eval(['wh_grid',num2str(annual(k))]) ,'precision','%.0f','delimiter',' ');
    end;

end

dlmwrite('processed/gfwhd.2100.txt',wh2100*700*0.8*0.48*1.3/1000 ,'precision','%.0f','delimiter',' ');

