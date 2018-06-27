clear all

cellarea_half_deg=importdata('Z:\links\tarotdata\backup\projects\glm\inputs\cellarea\cellarea_halfdeg.txt');
ccodes = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\other\ccodes\add_ak\ccodes_half_deg.txt');

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


codelist=unique(ccodes);
codelist=codelist(2:end);
years=[2005,2010:10:2100];
t_size=length(years);
c_size=length(codelist);

wh_country=zeros(t_size,c_size);

for t=1:t_size
    for i=1:length(codelist)
        sites=ccodes==codelist(i);
        eval(['wh_country(t,i)=sum(sum(wh',num2str(years(t)),'be([sites])+wh',num2str(years(t)),'sw([sites])))*1.3;'])
        %eval(['wh_country(t,i)=sum(sum(wh',num2str(years(t)),'sw([sites])));'])
        %wh_country(t,i)=sum(sum(wh2000be([sites])+wh2000sw([sites])));
    end;
end;

rcodes = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\other\wood_harvest\future\image\codes2glm_halfdeg_new3.txt');    
r_size=24;
wh_region=zeros(t_size,r_size);

for t=1:t_size
    for r=1:r_size
        r_sites=find(rcodes(:,1)==r);
        wh_region(t,r)=sum(wh_country(t,[r_sites]));
    end;
end;
    
for t=1:(t_size-1)
    wh_annual((years(t)-2004):(years(t+1)-2004),:)=interp1([years(t),years(t+1)],wh_region(t:t+1,:),years(t):years(t+1));
end;

wh_annual(96,:) = wh_region(end,:);

wh_t6 = reshape(wh_annual(1:46,:)',24*46,1);
wh_t7 = reshape(wh_annual(46:96,:)',24*51,1);

dlmwrite('rcp_wh_iiasa4.tsix',wh_t6,'precision','%.0f');
dlmwrite('rcp_wh_iiasa4.tseven',wh_t7,'precision','%.0f');
