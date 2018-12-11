ccodes = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\other\ccodes\add_ak\ccodes_half_deg.txt');
rcodes = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\other\wood_harvest\future\image\codes2glm_minicam.txt',' ',0);

rmap=zeros(size(ccodes));

for i=1:max(size(rcodes));
    country=rcodes(i,2);
    csites=find(ccodes==country);
    rmap(csites)=rcodes(i,1);
end;

dlmwrite('rmap.txt',rmap ,'precision','%.0f','delimiter',' ');