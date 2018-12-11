% This script takes the MiniCAM crop and pasture data and computes new
% harmonized grids that allow smooth transitions from the HYDE grids.

clear all

crop_forest_abandon_percent = 0.9;
past_forest_abandon_percent = 0.9;

% load the HYDE half degree grids in the harmonization year (2005)
GLMcrop2005 = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\hyde_3.0\half_deg_grids\gcrop.2005.txt',' ',6);
GLMpast2005 = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\hyde_3.0\half_deg_grids\gpast.2005.txt',' ',6);
GLMothr2005 = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\hyde_3.0\half_deg_grids\gothr.2005.txt',' ',6);
GLMicew2005 = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\hyde_3.0\half_deg_grids\gicew.1700.txt',' ',6);

GLMcrop2005 = GLMcrop2005.data;
GLMpast2005 = GLMpast2005.data;
GLMothr2005 = GLMothr2005.data;
GLMicew2005 = GLMicew2005.data;

cellarea_half_deg=importdata('Z:\links\tarotdata\backup\projects\glm\inputs\cellarea\cellarea_halfdeg.txt');

landarea_halfdeg = GLMpast2005 + GLMcrop2005 + GLMothr2005;
pot_veg = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\other\miami_biomass_v3\miami_halfdeg_conform.txt',' ');
pot_veg = pot_veg*0.75;
fnf = (pot_veg>2);

% load the MiniCAM crop and pasture data (that has already been extracted
% from a spreadsheet and placed in this txt-file)
minicam_data=importdata('RCP_MiniCAM.txt',' ');
minicam_crop=minicam_data(2:3:41,:);
minicam_past=minicam_data(3:3:42,:);
% % average the year 2000 and 2010 data to get data for 2005
% minicam_crop(:,1) = (minicam_crop(:,1)+minicam_crop(:,2))/2;
% minicam_past(:,1) = (minicam_past(:,1)+minicam_past(:,2))/2;

% load the MiniCAM map of 14 world regions
rmap=importdata('rmap.txt',' ');
% find gridcells that have a region code associated with them
total_rsites = find(ismember(rmap,[1,2,3,4,5,6,7,8,9,10,11,12,13,14]));

years=[2005,2010:10:2100];

for ind=1:(length(years)-1)
    years(ind)

    % for gridcells that do not have a region code, assume that they have
    % the same crop and pasture fractions as HYDE 2005 (and that these do
    % not change from 2005-2100). For gridcells that do have a region code,
    % set crop and pasture to zero (to begin with)
    eval(['GLMpast',num2str(years(ind+1)),'=GLMpast2005;'])
    eval(['GLMcrop',num2str(years(ind+1)),'=GLMcrop2005;'])
    eval(['GLMpast',num2str(years(ind+1)),'(total_rsites)=0;'])
    eval(['GLMcrop',num2str(years(ind+1)),'(total_rsites)=0;'])   

    for r=1:14
        r
        % find locations with region code "r"
        r_sites=find(rmap==r);
        % get regional crop and pasture data from MiniCAM at time t and t+1
        crop1 = minicam_crop(r,ind);
        crop2 = minicam_crop(r,ind+1);
        past1 = minicam_past(r,ind);
        past2 = minicam_past(r,ind+1);

        % compute changes in regional crop and pasture between time t and 
        % t+1 (factor of 10 is due to MiniCAM raw data needing a factor of 
        % 10 to make it in same units as HYDE/GLM data)
        crop_d = (crop2-crop1)*10;
        past_d = (past2-past1)*10;

        % if the regional pasture change is zero, the new gridded pasture
        % values for all locations within the region "r" will be the same
        % as the previous timestep
        % Note: as of March 10, 2009 there are no pasture changes in
        % MiniCAM data, hence there are no checks for increases or
        % decreases to regional pasture
        if past_d==0
            eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites);'])
            if min(min(GLMpast2010))<0
                keyboard
            end;
        end;

        % if the regional crop change is negative, apply the regional
        % *percentage* change to cropland in all gridcells within the region
        % then apply the pasture change
        if crop_d<0
            disp('crop decrease')
            % try to abandon on naturally forested land
            eval(['forested_crop = sum(sum(GLMcrop',num2str(years(ind)),'.*(rmap==r).*fnf.*(GLMcrop',num2str(years(ind)),'>0).*cellarea_half_deg));'])
            eval(['nonforested_crop = sum(sum(GLMcrop',num2str(years(ind)),'.*(rmap==r).*(fnf==0).*(GLMcrop',num2str(years(ind)),'>0).*cellarea_half_deg));'])
            if (abs(crop_d)*crop_forest_abandon_percent<=forested_crop)&(abs(crop_d)*(1-crop_forest_abandon_percent)<=nonforested_crop)
                disp('case 1')
               % keyboard
                forested_crop_percent = crop_d*crop_forest_abandon_percent/forested_crop;
                nonforested_crop_percent = crop_d*(1-crop_forest_abandon_percent)/nonforested_crop;
                eval(['GLMcrop',num2str(years(ind+1)),'(r_sites)=GLMcrop',num2str(years(ind)),'(r_sites)+GLMcrop',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_crop_percent+GLMcrop',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_crop_percent;'])
                %eval(['GLMcrop',num2str(years(ind+1)),'(r_sites)=GLMcrop',num2str(years(ind)),'(r_sites)+GLMcrop',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_crop_percent;'])
            elseif (abs(crop_d)*crop_forest_abandon_percent>forested_crop)&(abs(crop_d)*(1-crop_forest_abandon_percent)<=nonforested_crop)
                disp('case 2')
                forested_crop_percent = -1;
                nonforested_crop_percent = -(abs(crop_d)-forested_crop)/nonforested_crop;
                eval(['GLMcrop',num2str(years(ind+1)),'(r_sites)=GLMcrop',num2str(years(ind)),'(r_sites)+GLMcrop',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_crop_percent+GLMcrop',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_crop_percent;'])
                %eval(['GLMcrop',num2str(years(ind+1)),'(r_sites)=GLMcrop',num2str(years(ind)),'(r_sites)+GLMcrop',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_crop_percent;'])
            elseif (abs(crop_d)*crop_forest_abandon_percent<=forested_crop)&(abs(crop_d)*(1-crop_forest_abandon_percent)>nonforested_crop)
                disp('case 3')
                forested_crop_percent = -(abs(crop_d)-nonforested_crop)/forested_crop;
                nonforested_crop_percent = -1;
                eval(['GLMcrop',num2str(years(ind+1)),'(r_sites)=GLMcrop',num2str(years(ind)),'(r_sites)+GLMcrop',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_crop_percent+GLMcrop',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_crop_percent;'])
                %eval(['GLMcrop',num2str(years(ind+1)),'(r_sites)=GLMcrop',num2str(years(ind)),'(r_sites)+GLMcrop',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_crop_percent;'])
            else
                disp('case 4')
              end;
                
              if past_d<0
                  disp('pasture decrease')
                  % try to abandon on naturally forested land
                  eval(['forested_past = sum(sum(GLMpast',num2str(years(ind)),'.*(rmap==r).*fnf.*(GLMpast',num2str(years(ind)),'>0).*cellarea_half_deg));'])
                  eval(['nonforested_past = sum(sum(GLMpast',num2str(years(ind)),'.*(rmap==r).*(fnf==0).*(GLMpast',num2str(years(ind)),'>0).*cellarea_half_deg));'])
                  if (abs(past_d)*past_forest_abandon_percent<=forested_past)&(abs(past_d)*(1-past_forest_abandon_percent)<=nonforested_past)
                      disp('case 1')
                      forested_past_percent = past_d*past_forest_abandon_percent/forested_past;
                      nonforested_past_percent = past_d*(1-past_forest_abandon_percent)/nonforested_past;
                      eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_past_percent+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                     % eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                  elseif (abs(past_d)*past_forest_abandon_percent>forested_past)&(abs(past_d)*(1-past_forest_abandon_percent)<=nonforested_past)
                      disp('case 2')
                      forested_past_percent = -1;
                      nonforested_past_percent = -(abs(past_d)-forested_past)/nonforested_past;
                      eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_past_percent+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                     % eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                  elseif (abs(past_d)*past_forest_abandon_percent<=forested_past)&(abs(past_d)*(1-past_forest_abandon_percent)>nonforested_past)
                      disp('case 3')
                      forested_past_percent = -(abs(past_d)-nonforested_past)/forested_past;
                      nonforested_past_percent = -1;
                      eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_past_percent+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                    %  eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                  else
                      disp('case 4')
                  end;
              
            elseif past_d>0
                avail_land0 = zeros(360,720);
                eval(['avail_land0(r_sites)=(landarea_halfdeg(r_sites)-GLMcrop',num2str(years(ind)),'(r_sites)-GLMpast',num2str(years(ind)),'(r_sites)).*cellarea_half_deg(r_sites).*(GLMpast',num2str(years(ind)),'(r_sites)>0);'])
                if sum(sum(avail_land0(r_sites)))>past_d
                    disp('pasture increase - land available')
                    eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=(GLMpast',num2str(years(ind)),'(r_sites).*cellarea_half_deg(r_sites)+avail_land0(r_sites)./sum(sum(avail_land0(r_sites)))*past_d)./cellarea_half_deg(r_sites);'])
                    if min(min(GLMpast2010))<0
                        keyboard
                    end;
                else
                    disp('pasture increase - land not available')
                    clear all
                end
            end;
        end;

        % if the regional crop change is positive, first calculate
        % available land for crop expansion. If there is not enough land
        % available for expansion, halt the simulation, if there *is*
        % enough land available for expansion, apply the crop increase to
        % all gridcells in the region, weighted by available land in each
        % cell
        if crop_d>0
            avail_land0 = zeros(360,720);
            eval(['avail_land0(r_sites)=(landarea_halfdeg(r_sites)-GLMcrop',num2str(years(ind)),'(r_sites)-GLMpast',num2str(years(ind)),'(r_sites)).*cellarea_half_deg(r_sites).*(GLMcrop',num2str(years(ind)),'(r_sites)>0);'])
            if sum(sum(avail_land0(r_sites)))>=crop_d
                disp('crop increase - land available')
                eval(['GLMcrop',num2str(years(ind+1)),'(r_sites)=(GLMcrop',num2str(years(ind)),'(r_sites).*cellarea_half_deg(r_sites)+avail_land0(r_sites)./sum(sum(avail_land0(r_sites)))*crop_d)./cellarea_half_deg(r_sites);'])
                if past_d<0
                    disp('pasture decrease')
                    % try to abandon on naturally forested land
                    eval(['forested_past = sum(sum(GLMpast',num2str(years(ind)),'.*(rmap==r).*fnf.*(GLMpast',num2str(years(ind)),'>0).*cellarea_half_deg));'])
                    eval(['nonforested_past = sum(sum(GLMpast',num2str(years(ind)),'.*(rmap==r).*(fnf==0).*(GLMpast',num2str(years(ind)),'>0).*cellarea_half_deg));'])
                    if (abs(past_d)*past_forest_abandon_percent<=forested_past)&(abs(past_d)*(1-past_forest_abandon_percent)<=nonforested_past)
                        disp('case 1')
                        forested_past_percent = past_d*past_forest_abandon_percent/forested_past;
                        nonforested_past_percent = past_d*(1-past_forest_abandon_percent)/nonforested_past;
                        eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_past_percent+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                    %    eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                    elseif (abs(past_d)*past_forest_abandon_percent>forested_past)&(abs(past_d)*(1-past_forest_abandon_percent)<=nonforested_past)
                        disp('case 2')
                        forested_past_percent = -1;
                        nonforested_past_percent = -(abs(past_d)-forested_past)/nonforested_past;
                        eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_past_percent+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                     %   eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                    elseif (abs(past_d)*past_forest_abandon_percent<=forested_past)&(abs(past_d)*(1-past_forest_abandon_percent)>nonforested_past)
                        disp('case 3')
                        forested_past_percent = -(abs(past_d)-nonforested_past)/forested_past;
                        nonforested_past_percent = -1;
                        eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_past_percent+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                      %  eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                    else
                        disp('case 4')
                    end;

                elseif past_d>0
                    avail_land0 = zeros(360,720);
                    eval(['avail_land0(r_sites)=(landarea_halfdeg(r_sites)-GLMcrop',num2str(years(ind)),'(r_sites)-GLMpast',num2str(years(ind)),'(r_sites)).*cellarea_half_deg(r_sites).*(GLMpast',num2str(years(ind)),'(r_sites)>0);'])
                    if sum(sum(avail_land0(r_sites)))>past_d
                        disp('pasture increase - land available')
                        eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=(GLMpast',num2str(years(ind)),'(r_sites).*cellarea_half_deg(r_sites)+avail_land0(r_sites)./sum(sum(avail_land0(r_sites)))*past_d)./cellarea_half_deg(r_sites);'])
                        if min(min(GLMpast2010))<0
                            keyboard
                        end;
                    else
                        disp('pasture increase - land not available')
                        clear all
                    end
                end;
            else
                if past_d<0
                    disp('pasture decrease')
                    % try to abandon on naturally forested land
                    eval(['forested_past = sum(sum(GLMpast',num2str(years(ind)),'.*(rmap==r).*fnf.*(GLMpast',num2str(years(ind)),'>0).*cellarea_half_deg));'])
                    eval(['nonforested_past = sum(sum(GLMpast',num2str(years(ind)),'.*(rmap==r).*(fnf==0).*(GLMpast',num2str(years(ind)),'>0).*cellarea_half_deg));'])
                    if (abs(past_d)*past_forest_abandon_percent<=forested_past)&(abs(past_d)*(1-past_forest_abandon_percent)<=nonforested_past)
                        disp('case 1')
                        forested_past_percent = past_d*past_forest_abandon_percent/forested_past;
                        nonforested_past_percent = past_d*(1-past_forest_abandon_percent)/nonforested_past;
                        eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_past_percent+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                     %   eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                    elseif (abs(past_d)*past_forest_abandon_percent>forested_past)&(abs(past_d)*(1-past_forest_abandon_percent)<=nonforested_past)
                        disp('case 2')
                        forested_past_percent = -1;
                        nonforested_past_percent = -(abs(past_d)-forested_past)/nonforested_past;
                        eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_past_percent+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                      %  eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                    elseif (abs(past_d)*past_forest_abandon_percent<=forested_past)&(abs(past_d)*(1-past_forest_abandon_percent)>nonforested_past)
                        disp('case 3')
                        forested_past_percent = -(abs(past_d)-nonforested_past)/forested_past;
                        nonforested_past_percent = -1;
                        eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*fnf(r_sites).*forested_past_percent+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                       % eval(['GLMpast',num2str(years(ind+1)),'(r_sites)=GLMpast',num2str(years(ind)),'(r_sites)+GLMpast',num2str(years(ind)),'(r_sites).*(fnf(r_sites)==0).*nonforested_past_percent;'])
                    else
                        disp('case 4')
                    end;
                  
                    avail_land0 = zeros(360,720);
                    eval(['avail_land0(r_sites)=(landarea_halfdeg(r_sites)-GLMcrop',num2str(years(ind)),'(r_sites)-GLMpast',num2str(years(ind)),'(r_sites)).*cellarea_half_deg(r_sites).*(GLMcrop',num2str(years(ind)),'(r_sites)>0);'])
                    if sum(sum(avail_land0(r_sites)))>=crop_d
                        disp('crop increase - land available')
                        eval(['GLMcrop',num2str(years(ind+1)),'(r_sites)=(GLMcrop',num2str(years(ind)),'(r_sites).*cellarea_half_deg(r_sites)+avail_land0(r_sites)./sum(sum(avail_land0(r_sites)))*crop_d)./cellarea_half_deg(r_sites);'])
                    else
                        disp('crop increase - land not available')
                        clear all
                    end
                else
                    disp('crop increase - land not available')
                    clear all
                end;
            end
        end;        
    end;
    
    %interpolate to annual grids
    eval(['crop_range(1,:,:) = GLMcrop' num2str(years(ind)) ';']);
    eval(['crop_range(2,:,:) = GLMcrop' num2str(years(ind+1)) ';']);
    eval(['past_range(1,:,:) = GLMpast' num2str(years(ind)) ';']);
    eval(['past_range(2,:,:) = GLMpast' num2str(years(ind+1)) ';']);

    annual=years(ind):years(ind+1);
    annual_crop = interp1([years(ind),years(ind+1)],crop_range,annual);
    annual_past = interp1([years(ind),years(ind+1)],past_range,annual);

    clear crop_range
    clear past_range

    for k=1:(length(annual)-1)
        eval(['GLMcrop' num2str(annual(k)) '= reshape(annual_crop(k,:,:),[360,720]);']);
        eval(['GLMpast' num2str(annual(k)) '= reshape(annual_past(k,:,:),[360,720]);']);
 
        % save new annual grids
        dlmwrite(['processed/45_forestpref/gcrop.', num2str(annual(k)), '.txt'], eval(['GLMcrop',num2str(annual(k))]) ,'precision','%.6f','delimiter',' ');
        dlmwrite(['processed/45_forestpref/gpast.', num2str(annual(k)), '.txt'], eval(['GLMpast',num2str(annual(k))]) ,'precision','%.6f','delimiter',' ');
        dlmwrite(['processed/45_forestpref/gothr.', num2str(annual(k)), '.txt'], eval(['landarea_halfdeg-GLMpast',num2str(annual(k)),'-GLMcrop',num2str(annual(k))]) ,'precision','%.6f','delimiter',' ');
    end;
    
end

% save new grids in the year 2100
dlmwrite('processed/45_forestpref/gcrop.2100.txt',GLMcrop2100 ,'precision','%.6f','delimiter',' ');
dlmwrite('processed/45_forestpref/gpast.2100.txt',GLMpast2100 ,'precision','%.6f','delimiter',' ');
dlmwrite('processed/45_forestpref/gothr.2100.txt',landarea_halfdeg-GLMpast2100-GLMcrop2100,'precision','%.6f','delimiter',' ');


%figure(1)
%plot(years,sum(minicam_crop*10),'r')
%hold on
%plot(years,[sum(sum(GLMcrop2000.*cellarea_half_deg)),sum(sum(GLMcrop2010.*cellarea_half_deg)),sum(sum(GLMcrop2020.*cellarea_half_deg)),sum(sum(GLMcrop2030.*cellarea_half_deg)),sum(sum(GLMcrop2040.*cellarea_half_deg)),sum(sum(GLMcrop2050.*cellarea_half_deg)),sum(sum(GLMcrop2060.*cellarea_half_deg)),sum(sum(GLMcrop2070.*cellarea_half_deg)),sum(sum(GLMcrop2080.*cellarea_half_deg)),sum(sum(GLMcrop2090.*cellarea_half_deg)),sum(sum(GLMcrop2100.*cellarea_half_deg))])

%figure(2)
%plot(years,sum(minicam_past*10),'r')
%hold on
%plot(years,[sum(sum(GLMpast2000.*cellarea_half_deg)),sum(sum(GLMpast2010.*cellarea_half_deg)),sum(sum(GLMpast2020.*cellarea_half_deg)),sum(sum(GLMpast2030.*cellarea_half_deg)),sum(sum(GLMpast2040.*cellarea_half_deg)),sum(sum(GLMpast2050.*cellarea_half_deg)),sum(sum(GLMpast2060.*cellarea_half_deg)),sum(sum(GLMpast2070.*cellarea_half_deg)),sum(sum(GLMpast2080.*cellarea_half_deg)),sum(sum(GLMpast2090.*cellarea_half_deg)),sum(sum(GLMpast2100.*cellarea_half_deg))])


