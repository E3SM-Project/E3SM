% IMAGE WH processing
clear all

% Load data extracted from IMAGE WH spreadsheet into MATLAB. Call this
% matrix IMAGE_WH (81x25)
load IMAGE_WH_data
%keyboard

% Add Ukraine(14), Russia(16), and Asia-Stan(15) together to get USSR. Put USSR column
% at end of spreadhseet (replacing "World"). Delete Ukraine, Asia-Stan and
% Russia from spreadsheet. Add two columns of zeros after USSR for
% Greenland and Antarctica
IMAGE_WH(:,25)=IMAGE_WH(:,14)+IMAGE_WH(:,15)+IMAGE_WH(:,16);
IMAGE_WH(:,26)=zeros(81,1);
IMAGE_WH(:,27)=zeros(81,1);
IMAGE_WH = [IMAGE_WH(:,1:13) IMAGE_WH(:,17:27)];

% For each year, add all three WH types together
for i=1:(81/3)
    temp(i,:)=sum(IMAGE_WH(((i-1)*3+1):i*3,:));
end;
clear IMAGE_WH
IMAGE_WH = temp;
clear temp

% load GLM WH reconstruction (by country) in year 2000 (192 countries) 
%load wh_reg_data2 wh2000
% load matrix that assigns each country to an IMAGE region
%codes2glm_halfdeg_new3=importdata('Z:\links\tarotdata\backup\projects\glm\inputs\other\wood_harvest\future\image\codes2glm_halfdeg_new3.txt',' ');
% load codes for the 192 countries in GLM (in numerical order ie woodharvest order)
%load wh_reg_data2 ccodes

% compute the wood harvest (in carbon units) for each IMAGE region in the
% year 2000 by aggregating the year 2000 country-level wood harvest data
% into the IMAGE regions (24 of them)
%wh_reg_2000 = zeros(1,24);
%for i=1:192
%    ind=codes2glm_halfdeg_new3(i,:);
%    j=find(ccodes==ind(2));
%    wh_reg_2000(ind(1)) = wh_reg_2000(ind(1))+wh2000(j);
%end;

% load the reg_ordering vector that gives the region code associated with
% each column of the IMAGE WH data
load wh_reg_data2 reg_ordering2

% compute conversion factor for converting IMAGE data from m^3 to MgC by
% comparing regional WH data in year 2000 from GLM and IMAGE. Note that row
% 7 is year 2000 for IMAGE data (first row is 1970 and rows are separated
% by 5 years)
%conv_factor=zeros(1,24);
%for i=1:length(reg_ordering2)
%    conv_factor(i)=wh_reg_2000(reg_ordering2(i))/(IMAGE_WH(7,i)*1e3);
%end;

% compute IMAGE conversion factor plus slash fraction
conv_factor = 700*0.8*0.48*1.3/1000;

% multiply IMAGE WH data by conversion factor
IMAGE_WH = IMAGE_WH*1000*conv_factor;

% re-order IMAGE WH so region codes are in numerical order
for i=1:length(reg_ordering2)
    new_wh(:,reg_ordering2(i)) = IMAGE_WH(:,i);
end;

% truncate new_wh so it only includes WH data from 2005 onwards (ie index 8
% onwards)
new_wh=new_wh(8:27,:);

% separate USA into USA and Alaska (1%)
% new_wh(:,25) = new_wh(:,2)*0.01;
% new_wh(:,2) = new_wh(:,2)*0.99;
% keyboard
for i=1:19
    wh_annual((i-1)*5+1:(i-1)*5+5,:)=interp1([(i-1)*5,i*5],new_wh(i:i+1,:),(i-1)*5:(i-1)*5+4);
end;

wh_annual(101,:) = new_wh(20,:);

wh_t6 = reshape(wh_annual(1:46,:)',24*46,1);
wh_t7 = reshape(wh_annual(46:96,:)',24*51,1);

dlmwrite('rcp_wh_image4.tsix',wh_t6,'precision','%.0f');
dlmwrite('rcp_wh_image4.tseven',wh_t7,'precision','%.0f');



