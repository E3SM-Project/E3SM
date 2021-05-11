comp_G=1;
if (exist('G'))
   comp_G=0;
end
dirname='/scratch3/mataylo/preqx/';
dirname='/users/mataylo/data/preqx/';

t160_1 = [dirname,'hs-ne42-4t60l20-hnu2e14-1/zeta168480'];
t160_1b= [dirname,'hs-ne18t30l20-hnu2e14-3/zeta-surf360000'];

t160_2 = [dirname,'hs-ne42-4t60l20-hnu6e14-3/zeta-surf172800'];   % BEST?
t160_3 = [dirname,'hs-ne42-4t60l20-hnu2e15-10/zeta-surf108000'];
t160_3b = [dirname,'hs-ne42-4t60l20-hpnu2e15-10/zeta-surf100800'];


t85_3 = [dirname,'hs-ne21-4t120l20-hnu3e15-1/zeta084240'];
t85_3b= [dirname,'hs-ne9t60l20-hnu3e15-3/zeta168480'];
t85_4 = [dirname,'hs-ne21-4t120l20-hnu1.2e16-2/zeta084240'];
t85_5 = [dirname,'hs-ne21-4t120l20-hnu3e16-6/zeta084240'];
t85_5b= [dirname,'hs-ne21-4t120l20-hp10nu3e16-6/zeta-surf086400'];

thick=2;
figure(2); clf
set(gca,'FontWeight','bold')
set(gca,'FontSize',15)
set(gca,'LineWidth',2.0)

if (1)

data=textread(t160_1);
data=data(1:2:288,1:2:576);
%data=data(1:256,1:512);
if (comp_G) 
   [k,e,G]=spec(data);
else
  [k,e]=spec(data,G);
end
size(k)
kminv = k/(6e3*pi);
loglog(kminv,e,'b','linewidth',thick); hold on;
%loglog(k,e,'b','linewidth',thick); hold on;
plot(1/75,1e-13,'o')   % spectral element max res
plot(1/112,1e-13,'o')  % 2/3 dealiasing cutoff
plot(1/225,1e-13,'o')  % element length: 6e3*2*pi / (4*42)

data=textread(t160_1b);
[k,e]=spec(data,G);
loglog(k,e,'y'); hold on

data=textread(t160_1b);
data=data(1:2:288,1:2:576);
[k,e]=spec(data,G);
loglog(kminv,e,'y','linewidth',thick); hold on




data=textread(t160_2);
data=data(1:2:288,1:2:576);
[k,e]=spec(data,G);
loglog(kminv,e,'r','linewidth',thick); hold on

data=textread(t160_3);
data=data(1:2:288,1:2:576);
[k,e]=spec(data,G);
loglog(kminv,e,'g','linewidth',thick); hold on


%data=textread(t160_3b);
%data=data(1:2:288,1:2:576);
%[k,e]=spec(data,G);
%loglog(kminv,e,'y','linewidth',thick); hold on


end

if (0)
data=textread(t85_3);
[k,e]=spec(data);
kminv = 1./(6e3*pi./k);
loglog(kminv,e,'b')
hold on

data=textread(t85_3b);
[k,e]=spec(data);
loglog(kminv,e,'y')
hold on

data=textread(t85_4);
[k,e]=spec(data);
loglog(kminv,e,'r')
hold on

data=textread(t85_5);
[k,e]=spec(data);
loglog(kminv,e,'g')
hold on

data=textread(t85_5b);
[k,e]=spec(data);
loglog(kminv,e,'y')
hold on



end


title('Enstrophy')
E=5e-10 * k.^-1;
plot(kminv,E,'k')



hold off;
axis([1/10000,1/100,1e-15,1e-10]);
xlabel('km^{-1}');
%xticks(['10000';' 1000';'  100'],'km');


