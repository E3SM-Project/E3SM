% This script takes the IMAGE crop and pasture grids and computes new
% harmonized grids that allow smooth transitions from the HYDE grids.

clear all

GLMcrop2005 = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\hyde_3.0\half_deg_grids\gcrop.2005.txt',' ',6);
GLMpast2005 = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\hyde_3.0\half_deg_grids\gpast.2005.txt',' ',6);
GLMothr2005 = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\hyde_3.0\half_deg_grids\gothr.2005.txt',' ',6);
GLMicew2005 = importdata('Z:\links\tarotdata\backup\projects\glm\inputs\hyde_3.0\half_deg_grids\gicew.1700.txt',' ',6);

GLMcrop2005 = GLMcrop2005.data;
GLMpast2005 = GLMpast2005.data;
GLMothr2005 = GLMothr2005.data;
GLMicew2005 = GLMicew2005.data;

cellarea_half_deg=importdata('Z:\links\tarotdata\backup\projects\glm\inputs\cellarea\cellarea_halfdeg.txt');
cellarea_1deg=importdata('Z:\links\tarotdata\backup\projects\glm\inputs\cellarea\cellarea_1deg.txt');

% create 2deg cellarea grid
[dataN,dataM]=size(cellarea_1deg);
cellarea_2deg = zeros(90,180);
ratio = dataN/90;

for n = 1:dataN
    for m = 1:dataM
        if rem(n,ratio) == 0
            i = n/ratio;
        else
            i = (n+(ratio-rem(n,ratio)))/ratio;
        end;
        if rem(m,ratio) == 0
            j = m/ratio;
        else
            j = (m+(ratio-rem(m,ratio)))/ratio;
        end;
        cellarea_2deg(i,j) = cellarea_2deg(i,j) + cellarea_1deg(n,m);
    end
end

[dataN,dataM]=size(cellarea_half_deg);
ratio = dataN/90;

% aggregate to 2deg grids
GLMcrop2005_2deg = zeros(90,180);
GLMpast2005_2deg = zeros(90,180);
GLMicew2005_2deg = zeros(90,180);
GLMothr2005_2deg = zeros(90,180);
for n = 1:dataN
    for m = 1:dataM
        if rem(n,ratio) == 0
            i = n/ratio;
        else
            i = (n+(ratio-rem(n,ratio)))/ratio;
        end;
        if rem(m,ratio) == 0
            j = m/ratio;
        else
            j = (m+(ratio-rem(m,ratio)))/ratio;
        end;

        GLMothr2005_2deg(i,j) = GLMothr2005_2deg(i,j) + GLMothr2005(n,m)*cellarea_half_deg(n,m)/cellarea_2deg(i,j);
        GLMcrop2005_2deg(i,j) = GLMcrop2005_2deg(i,j) + GLMcrop2005(n,m)*cellarea_half_deg(n,m)/cellarea_2deg(i,j);
        GLMpast2005_2deg(i,j) = GLMpast2005_2deg(i,j) + GLMpast2005(n,m)*cellarea_half_deg(n,m)/cellarea_2deg(i,j);
        GLMicew2005_2deg(i,j) = GLMicew2005_2deg(i,j) + GLMicew2005(n,m)*cellarea_half_deg(n,m)/cellarea_2deg(i,j);

    end;
end;

landarea = GLMpast2005_2deg + GLMcrop2005_2deg + GLMothr2005_2deg;
landarea_halfdeg = GLMpast2005 + GLMcrop2005 + GLMothr2005;

years=[2005,2010:10:2100];

for ind=1:(length(years)-1)
    %for ind=13:14
    years(ind)
    
    if years(ind)==2005
        crop1_2000=importdata('Cropland301/CR2000A2.asc',' ',5);
        crop1_2010=importdata('Cropland301/CR2010A2.asc',' ',5);
        crop1_2000=crop1_2000.data/100;
        crop1_2010=crop1_2010.data/100;
        crop1=(crop1_2000+crop1_2010)/2;
        
        crop2=importdata('Cropland301/CR2010A2.asc',' ',5);
        crop2=crop2.data/100;
        
        past1_2000=importdata('GrasslWoodland301/GR2000A2.asc',' ',5);
        past1_2010=importdata('GrasslWoodland301/GR2010A2.asc',' ',5);
        gcl3_1_2000 = importdata('grassland-classes_17Sep08\processed/erelcl3.2000.txt',' ');
        gcl4_1_2000 = importdata('grassland-classes_17Sep08\processed/erelcl4.2000.txt',' ');
        gcl3_1_2010 = importdata('grassland-classes_17Sep08\processed/erelcl3.2010.txt',' ');
        gcl4_1_2010 = importdata('grassland-classes_17Sep08\processed/erelcl4.2010.txt',' ');
        
        past1_2000=past1_2000.data/100.*(gcl3_1_2000+gcl4_1_2000);
        past1_2010=past1_2010.data/100.*(gcl3_1_2010+gcl4_1_2010);
        past1=(past1_2000+past1_2010)/2;
        
        past2=importdata(['GrasslWoodland301/GR',num2str(years(ind+1)),'A2.asc'],' ',5);
        gcl3_2 = importdata(['grassland-classes_17Sep08\processed/erelcl3.',num2str(years(ind+1)),'.txt'],' ');
        gcl4_2 = importdata(['grassland-classes_17Sep08\processed/erelcl4.',num2str(years(ind+1)),'.txt'],' ');
        past2=past2.data/100.*(gcl3_2+gcl4_2);
    else
        crop1=importdata(['Cropland301/CR',num2str(years(ind)),'A2.asc'],' ',5);
        crop2=importdata(['Cropland301/CR',num2str(years(ind+1)),'A2.asc'],' ',5);
        past1=importdata(['GrasslWoodland301/GR',num2str(years(ind)),'A2.asc'],' ',5);
        past2=importdata(['GrasslWoodland301/GR',num2str(years(ind+1)),'A2.asc'],' ',5);
        gcl3_1 = importdata(['grassland-classes_17Sep08\processed/erelcl3.',num2str(years(ind)),'.txt'],' ');
        gcl4_1 = importdata(['grassland-classes_17Sep08\processed/erelcl4.',num2str(years(ind)),'.txt'],' ');
        gcl3_2 = importdata(['grassland-classes_17Sep08\processed/erelcl3.',num2str(years(ind+1)),'.txt'],' ');
        gcl4_2 = importdata(['grassland-classes_17Sep08\processed/erelcl4.',num2str(years(ind+1)),'.txt'],' ');
        
        crop1=crop1.data/100;
        crop2=crop2.data/100;
        past1=past1.data/100.*(gcl3_1+gcl4_1);
        past2=past2.data/100.*(gcl3_2+gcl4_2);
    end;

    % aggregate to 2deg grids
    crop1_2deg = zeros(90,180);
    past1_2deg = zeros(90,180);
    crop2_2deg = zeros(90,180);
    past2_2deg = zeros(90,180);

    for n = 1:dataN
        for m = 1:dataM
            if rem(n,ratio) == 0
                i = n/ratio;
            else
                i = (n+(ratio-rem(n,ratio)))/ratio;
            end;
            if rem(m,ratio) == 0
                j = m/ratio;
            else
                j = (m+(ratio-rem(m,ratio)))/ratio;
            end;

            crop1_2deg(i,j) = crop1_2deg(i,j) + crop1(n,m)*cellarea_half_deg(n,m)/cellarea_2deg(i,j);
            crop2_2deg(i,j) = crop2_2deg(i,j) + crop2(n,m)*cellarea_half_deg(n,m)/cellarea_2deg(i,j);
            past1_2deg(i,j) = past1_2deg(i,j) + past1(n,m)*cellarea_half_deg(n,m)/cellarea_2deg(i,j);
            past2_2deg(i,j) = past2_2deg(i,j) + past2(n,m)*cellarea_half_deg(n,m)/cellarea_2deg(i,j);

        end;
    end;

    % compute changes
    crop_d = (crop2_2deg-crop1_2deg);
    past_d = (past2_2deg-past1_2deg);

    % apply changes to GLM grids (@2deg)
    eval(['GLMcrop',num2str(years(ind+1)),'_2deg=GLMcrop',num2str(years(ind)),'_2deg+crop_d;'])
    eval(['GLMpast',num2str(years(ind+1)),'_2deg=GLMpast',num2str(years(ind)),'_2deg+past_d;'])

    % check that neight crop nor pasture exceeds 1
    eval(['unmetcrop1_',num2str(years(ind+1)),'=(1-GLMcrop',num2str(years(ind+1)),'_2deg).*(GLMcrop',num2str(years(ind+1)),'_2deg>1);'])
    eval(['unmetpast1_',num2str(years(ind+1)),'=(1-GLMpast',num2str(years(ind+1)),'_2deg).*(GLMpast',num2str(years(ind+1)),'_2deg>1);'])
    eval(['GLMcrop',num2str(years(ind+1)),'_2deg=GLMcrop',num2str(years(ind+1)),'_2deg-unmetcrop1_',num2str(years(ind+1)),';'])
    eval(['GLMpast',num2str(years(ind+1)),'_2deg=GLMpast',num2str(years(ind+1)),'_2deg-unmetpast1_',num2str(years(ind+1)),';'])

    % check that neight crop nor pasture exceeds 0
    eval(['unmetcrop0_',num2str(years(ind+1)),'=GLMcrop',num2str(years(ind+1)),'_2deg.*(GLMcrop',num2str(years(ind+1)),'_2deg<0);'])
    eval(['unmetpast0_',num2str(years(ind+1)),'=GLMpast',num2str(years(ind+1)),'_2deg.*(GLMpast',num2str(years(ind+1)),'_2deg<0);'])
    eval(['GLMcrop',num2str(years(ind+1)),'_2deg=GLMcrop',num2str(years(ind+1)),'_2deg-unmetcrop0_',num2str(years(ind+1)),';'])
    eval(['GLMpast',num2str(years(ind+1)),'_2deg=GLMpast',num2str(years(ind+1)),'_2deg-unmetpast0_',num2str(years(ind+1)),';'])

    % check is land area is conserved
    eval(['exceed_land_',num2str(years(ind+1)),'=GLMcrop',num2str(years(ind+1)),'_2deg+GLMpast',num2str(years(ind+1)),'_2deg-landarea;'])
    eval(['unmetcrop',num2str(years(ind+1)),'=GLMcrop',num2str(years(ind+1)),'_2deg./(GLMcrop',num2str(years(ind+1)),'_2deg+GLMpast',num2str(years(ind+1)),'_2deg+1e-12).*exceed_land_',num2str(years(ind+1)),'.*(exceed_land_',num2str(years(ind+1)),'>0);'])
    eval(['unmetpast',num2str(years(ind+1)),'=GLMpast',num2str(years(ind+1)),'_2deg./(GLMcrop',num2str(years(ind+1)),'_2deg+GLMpast',num2str(years(ind+1)),'_2deg+1e-12).*exceed_land_',num2str(years(ind+1)),'.*(exceed_land_',num2str(years(ind+1)),'>0);'])
    eval(['GLMcrop',num2str(years(ind+1)),'_2deg=GLMcrop',num2str(years(ind+1)),'_2deg-unmetcrop',num2str(years(ind+1)),'.*(exceed_land_',num2str(years(ind+1)),'>0);'])
    eval(['GLMpast',num2str(years(ind+1)),'_2deg=GLMpast',num2str(years(ind+1)),'_2deg-unmetpast',num2str(years(ind+1)),'.*(exceed_land_',num2str(years(ind+1)),'>0);'])

    eval(['total_unmet_crop',num2str(years(ind+1)),'=unmetcrop0_',num2str(years(ind+1)),'+unmetcrop1_',num2str(years(ind+1)),'+unmetcrop',num2str(years(ind+1)),';'])
    eval(['total_unmet_past',num2str(years(ind+1)),'=unmetpast0_',num2str(years(ind+1)),'+unmetpast1_',num2str(years(ind+1)),'+unmetpast',num2str(years(ind+1)),';'])

    eval(['total_unmet_crop=total_unmet_crop',num2str(years(ind+1)),';'])
    eval(['total_unmet_past=total_unmet_past',num2str(years(ind+1)),';'])
    eval(['GLMcrop=GLMcrop',num2str(years(ind+1)),'_2deg;'])
    eval(['GLMpast=GLMpast',num2str(years(ind+1)),'_2deg;'])

    eval(['displaced_crop',num2str(years(ind+1)),'=zeros(90,180);'])
    eval(['displaced_past',num2str(years(ind+1)),'=zeros(90,180);'])

    crop_rings=zeros(90,180);
    past_rings=zeros(90,180);

    for k=0:89
        for m=0:179
            j=1;
            while abs(total_unmet_crop(k+1,m+1))>1e-8
                crop_rings(k+1,m+1)=j;
                if total_unmet_crop(k+1,m+1)>0
                    avail_space = landarea(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)-GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) - GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1);
                    total_avail_space = sum(sum(avail_space.*(avail_space>0).*cellarea_2deg(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)));
                    if total_unmet_crop(k+1,m+1)*cellarea_2deg(k+1,m+1) >= total_avail_space
                        GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) = GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) + avail_space.*(avail_space>0);
                        eval(['displaced_crop',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)=displaced_crop',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)+avail_space.*(avail_space>0);'])
                        total_unmet_crop(k+1,m+1) = (total_unmet_crop(k+1,m+1)*cellarea_2deg(k+1,m+1)-total_avail_space)/cellarea_2deg(k+1,m+1);
                    elseif total_unmet_crop(k+1,m+1)*cellarea_2deg(k+1,m+1) < total_avail_space
                        GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) = GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) + total_unmet_crop(k+1,m+1)*cellarea_2deg(k+1,m+1)*avail_space./total_avail_space.*(avail_space>0);
                        eval(['displaced_crop',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)=displaced_crop',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)+total_unmet_crop(k+1,m+1)*cellarea_2deg(k+1,m+1)*avail_space./total_avail_space.*(avail_space>0);'])
                        total_unmet_crop(k+1,m+1) = 0;
                    end
                elseif total_unmet_crop(k+1,m+1)<0
                    avail_space = GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1);
                    total_avail_space = sum(sum(avail_space.*(avail_space>0).*cellarea_2deg(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)));
                    if total_unmet_crop(k+1,m+1)*cellarea_2deg(k+1,m+1) <= -total_avail_space
                        GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) = GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) - avail_space.*(avail_space>0);
                        eval(['displaced_crop',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)=displaced_crop',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)-avail_space.*(avail_space>0);'])
                        total_unmet_crop(k+1,m+1) = (total_unmet_crop(k+1,m+1)*cellarea_2deg(k+1,m+1)+total_avail_space)/cellarea_2deg(k+1,m+1);
                    elseif total_unmet_crop(k+1,m+1)*cellarea_2deg(k+1,m+1) > -total_avail_space
                        GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) = GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) + total_unmet_crop(k+1,m+1)*cellarea_2deg(k+1,m+1)*avail_space./total_avail_space.*(avail_space>0);
                        eval(['displaced_crop',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)=displaced_crop',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)+total_unmet_crop(k+1,m+1)*cellarea_2deg(k+1,m+1)*avail_space./total_avail_space.*(avail_space>0);'])
                        total_unmet_crop(k+1,m+1) = 0;
                    end
                end;
                j=j+1;
            end

            while abs(total_unmet_past(k+1,m+1))>1e-8
                past_rings(k+1,m+1)=j;
                if total_unmet_past(k+1,m+1)>0
                    avail_space = landarea(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)-GLMcrop(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) - GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1);
                    total_avail_space = sum(sum(avail_space.*(avail_space>0).*cellarea_2deg(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)));
                    if total_unmet_past(k+1,m+1)*cellarea_2deg(k+1,m+1) >= total_avail_space
                        GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) = GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) + avail_space.*(avail_space>0);
                        eval(['displaced_past',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)=displaced_past',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)+avail_space.*(avail_space>0);'])
                        total_unmet_past(k+1,m+1) = (total_unmet_past(k+1,m+1)*cellarea_2deg(k+1,m+1)-total_avail_space)/cellarea_2deg(k+1,m+1);
                    elseif total_unmet_past(k+1,m+1)*cellarea_2deg(k+1,m+1) < total_avail_space
                        GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) = GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) + total_unmet_past(k+1,m+1)*cellarea_2deg(k+1,m+1)*avail_space./total_avail_space.*(avail_space>0);
                        eval(['displaced_past',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)=displaced_past',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)+total_unmet_past(k+1,m+1)*cellarea_2deg(k+1,m+1)*avail_space./total_avail_space.*(avail_space>0);'])
                        total_unmet_past(k+1,m+1) = 0;
                    end
                elseif total_unmet_past(k+1,m+1)<0
                    avail_space = GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1);
                    total_avail_space = sum(sum(avail_space.*(avail_space>0).*cellarea_2deg(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)));
                    if total_unmet_past(k+1,m+1)*cellarea_2deg(k+1,m+1) <= -total_avail_space
                        GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) = GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) - avail_space.*(avail_space>0);
                        eval(['displaced_past',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)=displaced_past',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)-avail_space.*(avail_space>0);'])
                        total_unmet_past(k+1,m+1) = (total_unmet_past(k+1,m+1)*cellarea_2deg(k+1,m+1)+total_avail_space)/cellarea_2deg(k+1,m+1);
                    elseif total_unmet_past(k+1,m+1)*cellarea_2deg(k+1,m+1) > -total_avail_space
                        GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) = GLMpast(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1) + total_unmet_past(k+1,m+1)*cellarea_2deg(k+1,m+1)*avail_space./total_avail_space.*(avail_space>0);
                        eval(['displaced_past',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)=displaced_past',num2str(years(ind+1)),'(mod(k+(-j:j),90)+1,mod(m+(-j:j),180)+1)+total_unmet_past(k+1,m+1)*cellarea_2deg(k+1,m+1)*avail_space./total_avail_space.*(avail_space>0);'])
                        total_unmet_past(k+1,m+1) = 0;
                    end
                end;
                j=j+1;
            end

        end
    end

    disp('Crop rings')
    sum(sum(crop_rings>0))
    sum(sum(crop_rings>3))/sum(sum(crop_rings>0))*100
    sum(sum(crop_rings>4))/sum(sum(crop_rings>0))*100
    sum(sum(crop_rings>5))/sum(sum(crop_rings>0))*100
    sum(sum(crop_rings>6))/sum(sum(crop_rings>0))*100

    disp('Pasture rings')
    sum(sum(past_rings>0))
    sum(sum(past_rings>3))/sum(sum(past_rings>0))*100
    sum(sum(past_rings>4))/sum(sum(past_rings>0))*100
    sum(sum(past_rings>5))/sum(sum(past_rings>0))*100
    sum(sum(past_rings>6))/sum(sum(past_rings>0))*100

    eval(['GLMcrop',num2str(years(ind+1)),'_2deg=GLMcrop;'])
    eval(['GLMpast',num2str(years(ind+1)),'_2deg=GLMpast;'])
    
    for i=1:90
        for j=1:180

            eval(['crop_d=GLMcrop',num2str(years(ind+1)),'_2deg(i,j)-GLMcrop',num2str(years(ind)),'_2deg(i,j);']);
            eval(['past_d=GLMpast',num2str(years(ind+1)),'_2deg(i,j)-GLMpast',num2str(years(ind)),'_2deg(i,j);']);
            eval(['avail_land = landarea_halfdeg(4*i-(0:3),4*j-3:4*j)-GLMcrop',num2str(years(ind)),'(4*i-(0:3),4*j-3:4*j)-GLMpast',num2str(years(ind)),'(4*i-(0:3),4*j-3:4*j);'])

            if crop_d<=0
                eval(['crop_p = crop_d/(GLMcrop',num2str(years(ind)),'_2deg(i,j)+1e-12);']);
                eval(['GLMcrop',num2str(years(ind+1)),'(4*i-(0:3),4*j-3:4*j) = GLMcrop',num2str(years(ind)),'(4*i-(0:3),4*j-3:4*j)*(1+crop_p);']);
            else
                if (sum(sum(avail_land.*cellarea_half_deg(4*i-(0:3),4*j-3:4*j)))<sum(sum(crop_d.*cellarea_2deg(i,j))))&(past_d<0)
                    eval(['past_p = past_d/(GLMpast',num2str(years(ind)),'_2deg(i,j)+1e-12);']);
                    eval(['GLMpast',num2str(years(ind+1)),'(4*i-(0:3),4*j-3:4*j) = GLMpast',num2str(years(ind)),'(4*i-(0:3),4*j-3:4*j)*(1+past_p);']);
                    eval(['avail_land_crop = landarea_halfdeg(4*i-(0:3),4*j-3:4*j)-GLMcrop',num2str(years(ind)),'(4*i-(0:3),4*j-3:4*j)-GLMpast',num2str(years(ind+1)),'(4*i-(0:3),4*j-3:4*j);'])
                    eval(['GLMcrop',num2str(years(ind+1)),'(4*i-(0:3),4*j-3:4*j) = GLMcrop',num2str(years(ind)),'(4*i-(0:3),4*j-3:4*j)+ avail_land_crop/(sum(sum(avail_land_crop.*cellarea_half_deg(4*i-(0:3),4*j-3:4*j)))+1e-12)*crop_d*cellarea_2deg(i,j);'])
                else
                    eval(['GLMcrop',num2str(years(ind+1)),'(4*i-(0:3),4*j-3:4*j) = GLMcrop',num2str(years(ind)),'(4*i-(0:3),4*j-3:4*j)+ avail_land/(sum(sum(avail_land.*cellarea_half_deg(4*i-(0:3),4*j-3:4*j)))+1e-12)*crop_d*cellarea_2deg(i,j);'])
                end;
            end;
            if (past_d<0)&(sum(sum(avail_land.*cellarea_half_deg(4*i-(0:3),4*j-3:4*j)))>sum(sum(crop_d.*cellarea_2deg(i,j))))
                eval(['past_p = past_d/(GLMpast',num2str(years(ind)),'_2deg(i,j)+1e-12);']);
                eval(['GLMpast',num2str(years(ind+1)),'(4*i-(0:3),4*j-3:4*j) = GLMpast',num2str(years(ind)),'(4*i-(0:3),4*j-3:4*j)*(1+past_p);']);
            elseif (past_d>=0)
                eval(['avail_land = landarea_halfdeg(4*i-(0:3),4*j-3:4*j)-GLMcrop',num2str(years(ind+1)),'(4*i-(0:3),4*j-3:4*j)-GLMpast',num2str(years(ind)),'(4*i-(0:3),4*j-3:4*j);'])
                eval(['GLMpast',num2str(years(ind+1)),'(4*i-(0:3),4*j-3:4*j) = GLMpast',num2str(years(ind)),'(4*i-(0:3),4*j-3:4*j)+ avail_land/(sum(sum(avail_land.*cellarea_half_deg(4*i-(0:3),4*j-3:4*j)))+1e-12)*past_d*cellarea_2deg(i,j);'])
            end;

            %             eval(['GLMcrop',num2str(years(ind+1)),'(4*i-(0:3),4*j-3:4*j) = landarea_halfdeg(4*i-(0:3),4*j-3:4*j)/(sum(sum(landarea_halfdeg(4*i-(0:3),4*j-3:4*j).*cellarea_half_deg(4*i-(0:3),4*j-3:4*j)))+1e-12)*GLMcrop',num2str(years(ind+1)),'_2deg(i,j)*cellarea_2deg(i,j);'])
            %             eval(['GLMpast',num2str(years(ind+1)),'(4*i-(0:3),4*j-3:4*j) = landarea_halfdeg(4*i-(0:3),4*j-3:4*j)/(sum(sum(landarea_halfdeg(4*i-(0:3),4*j-3:4*j).*cellarea_half_deg(4*i-(0:3),4*j-3:4*j)))+1e-12)*GLMpast',num2str(years(ind+1)),'_2deg(i,j)*cellarea_2deg(i,j);'])
        end
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

        dlmwrite(['processed/gcrop.', num2str(annual(k)), '.txt'], eval(['GLMcrop',num2str(annual(k))]) ,'precision','%.6f','delimiter',' ');
        dlmwrite(['processed/gpast.', num2str(annual(k)), '.txt'], eval(['GLMpast',num2str(annual(k))]) ,'precision','%.6f','delimiter',' ');
        dlmwrite(['processed/gothr.', num2str(annual(k)), '.txt'], eval(['landarea_halfdeg-GLMpast',num2str(annual(k)),'-GLMcrop',num2str(annual(k))]) ,'precision','%.6f','delimiter',' ');
    end;

end

dlmwrite('processed/gcrop.2100.txt',GLMcrop2100 ,'precision','%.6f','delimiter',' ');
dlmwrite('processed/gpast.2100.txt',GLMpast2100 ,'precision','%.6f','delimiter',' ');
dlmwrite('processed/gothr.2100.txt',landarea_halfdeg-GLMpast2100-GLMcrop2100 ,'precision','%.6f','delimiter',' ');
