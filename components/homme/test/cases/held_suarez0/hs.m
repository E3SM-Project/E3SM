function hs

% T85 NE=8
%vorname = '/scratch3/mataylo/preqx/hs-ne9t60l20-hnu3e15-3/zeta168480';
%psname = '/scratch3/mataylo/preqx/hs-ne9t60l20-hnu3e15-3/ps168480';

% T85  NE=4 

%vorname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hnu3e15-1/zeta084240';
%psname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hnu3e15-1/ps084240';

%vorname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hp10nu3e15-1/zeta-surf086400';
%psname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hp10nu3e15-1/ps086400';

%vorname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hp100nu3e15-1/zeta-surf086400';
%psname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hp100nu3e15-1/ps086400';

%vorname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hnu1.2e16-2/zeta084240';
%psname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hnu1.2e16-2/ps084240';

%vorname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hnu3e16-6/zeta084240';
%psname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hnu3e16-6/ps084240';

%psname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hp10nu3e16-6/ps086400';
%vorname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hp10nu3e16-6/zeta-surf086400';

%psname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hp100nu3e16-6/ps086400';
%vorname = '/scratch3/mataylo/preqx/hs-ne21-4t120l20-hp100nu3e16-6/zeta-surf086400';

% T160 NE=8
%psname = '/scratch3/mataylo/preqx/hs-ne18t30l20-hnu2e14-3/ps360000';
%vorname = '/scratch3/mataylo/preqx/hs-ne18t30l20-hnu2e14-3/zeta-surf360000';

% T160 NE=4
%psname = '/scratch3/mataylo/preqx/hs-ne42-4t60l20-hnu2e14-1/ps168480';
%vorname = '/scratch3/mataylo/preqx/hs-ne42-4t60l20-hnu2e14-1/zeta168480';

%vorname = '/scratch3/mataylo/preqx/hs-ne42-4t60l20-hnu6e14-3/zeta-surf172800';
%psname = '/scratch3/mataylo/preqx/hs-ne42-4t60l20-hnu6e14-3/ps172800';

%vorname = '/scratch3/mataylo/preqx/hs-ne42-4t60l20-hnu2e15-10/zeta-surf108000';
%psname = '/scratch3/mataylo/preqx/hs-ne42-4t60l20-hnu2e15-10/ps108000';

vorname = '/scratch3/mataylo/preqx/hs-ne42-4t60l20-hpnu2e15-10/zeta-surf100800';
psname = '/scratch3/mataylo/preqx/hs-ne42-4t60l20-hpnu2e15-10/ps100800';

%vorname = '/users/mataylo/data/preqx/hs-ne42-4t60l20-hpnu2e15-10/zeta-surf100800';
%psname = '/users/mataylo/data/preqx/hs-ne42-4t60l20-hpnu2e15-10/ps100800';


disp(vorname)
vor=textread(vorname);
ps=textread(psname);
disp(sprintf('ps min/max = %f %f',min(min(ps)),max(max(ps))))
disp(sprintf('vor min/max = %10.2e %10.2e',min(min(vor)),max(max(vor))))

l=length(vorname);
tname=vorname(25:l);

figure(1)
subplot(2,1,1)
plotzeta(vor,0)
title(tname);
subplot(2,1,2)
plotps(ps,0)

orient tall
print -depsc hs.ps


return




function plotzeta(data,region)
%
%  region = 0   sphere
%  region = 1   from paper
%  region = 2   larger version of 1    
%   
if (nargin==1) 
   region=1;       
end

nlat=size(data);
nlon=nlat(2);
nlat=nlat(1);

disp(sprintf('nlat x nlon = %d x %d ',nlat,nlon));

%contour(data,30)
%print -depsc interp.ps

bd=4e-4;
v=-bd : 2e-5 : bd;

lon=1:nlon;  lon=360*(lon-1)/nlon;
lat=1:nlat;  lat=180*(lat-1)/nlat-90;
contour(lon,lat,data,v)
caxis([-2e-4,2e-4])
colorbar

if (region==0)
   axis([0,360,-90,90]);
elseif (region==1) 
   axis([150,360,15,80]);
elseif (region==2)
   axis([0,360,0,90]);
end

return


function plotps(data,region)
%
nlat=size(data);
nlon=nlat(2);
nlat=nlat(1);

disp(sprintf('nlat x nlon = %d x %d ',nlat,nlon));

%contour(data,30)
%print -depsc interp.ps

v=950:5:1030;

lon=1:nlon;  lon=360*(lon-1)/nlon;
lat=1:nlat;  lat=180*(lat-1)/nlat-90;
contour(lon,lat,data,v)
caxis([950,1030])
colorbar

if (region==0)
   axis([0,360,-90,90]);
elseif (region==1) 
   axis([150,360,15,80]);
elseif (region==2)
   axis([0,360,0,90]);
end
return
