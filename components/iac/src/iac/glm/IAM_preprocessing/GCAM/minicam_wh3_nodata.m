% create MiniCAM wood harvest inputs for GLM land-use harmonization
clear all

load minicam_wh_data

% convert to carbon units
%minicam_wh=minicam_wh*700*0.8*0.48*1.3/1000;
minicam_wh=minicam_wh*1.3;
minicam_wh=minicam_wh';

[N,M] = size(minicam_wh);
wh_start = zeros(1,M);
wh_end = minicam_wh(N,:);
wh_annual = interp1([2005,2100],[wh_start;wh_end],2005:2100);
wh_length = length(2005:2100);

wh_t6 = reshape(wh_annual(1:46,:)',14*46,1);
wh_t7 = reshape(wh_annual(46:wh_length,:)',14*51,1);

dlmwrite('rcp_wh_minicam_nd.tsix',wh_t6,'precision','%.0f');
dlmwrite('rcp_wh_minicam_nd.tseven',wh_t7,'precision','%.0f');