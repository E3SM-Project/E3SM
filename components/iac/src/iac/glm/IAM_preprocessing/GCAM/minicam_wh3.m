% create MiniCAM wood harvest inputs for GLM land-use harmonization
clear all

load minicam_wh_data

% convert to carbon units
%minicam_wh=minicam_wh*700*0.8*0.48*1.3/1000;
minicam_wh=minicam_wh*1.3;
minicam_wh=minicam_wh';

% average the wh data for years 2000 and 2010 to get data for 2005
%minicam_wh(1,:) = (minicam_wh(1,:)+minicam_wh(2,:))/2;
%keyboard
% separate USA into USA + Alaska
%minicam

years = [2005,2010:10:2100];
% interpolate to annual time-steps
for i=1:(length(years)-1)
    wh_annual((years(i)-2004):(years(i+1)-2004),:)=interp1([years(i),years(i+1)],minicam_wh(i:i+1,:),years(i):years(i+1));
    %wh_annual((i-1)*10+1:(i-1)*10+10,:)=interp1([(i-1)*10,i*10],minicam_wh(i:i+1,:),(i-1)*10:(i-1)*10+9);
end;

wh_length = length(2005:2100);
wh_annual(wh_length,:) = minicam_wh(11,:);

wh_t6 = reshape(wh_annual(1:46,:)',14*46,1);
wh_t7 = reshape(wh_annual(46:wh_length,:)',14*51,1);

dlmwrite('ref_wh_minicam.tsix',wh_t6,'precision','%.0f');
dlmwrite('ref_wh_minicam.tseven',wh_t7,'precision','%.0f');