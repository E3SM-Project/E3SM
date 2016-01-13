function verifAE(nGlen, M0yr, Lkm, etamax, sector, tfyr, Nx, Mt);
%VERIFAE  Compares numerical to exact steady solution of shallow ice 
%    equation with basal sliding
%       0 = (M_s + M_b) - Div(q_f + H u_b),
%    incorporating compensatory accumulation.  Numerical method: type I 
%    explicit finite difference.  Takes advantage of radial symmetry to 
%    do 1/4 of work.  Basal sliding in sector.
%
%verifAE(n, M0, L, etamax, sector, tf, Nx, Mt);
%   n      = Glen exponent
%   M0     = accumulation in m/a
%   L      = margin radius (km)
%   etamax = maximum basal till thickness (mm)
%   sector = [r1 r2 theta1 theta2]  extent of sector; 
%            ri in km, thetai in degrees, 0<r1<r2<L, 0<=theta1<theta2<=90
%   tf     = final time (years; run is from t=0 to t=tf)
%   Nx     = number of grid intervals in both x and y
%   Mt     = number of time steps (dt=(tf-t0)/Mt=tadv/Mt)
%
%Notes: 
%   (1) Stability index  max(D)*dt/dx^2  of roughly 0.15 is stability limit
%       for test A.  Uses fixed time step but this index could be used to adapt.
%   (2) If etamax<0 then no basal sliding and sector is ignored.
%   (3) Displays in figures 1,2,3, and 4.
%   (4) Reference: Bueler et al (2004), "Exact solutions and the verification 
%       of numerical models for isothermal ice sheets", preprint.
%
%Examples:   (all near stability limit)
%TEST A:  (Bodvarsson-Vialov solution; no basal sliding)
%   >> verifAE(3,0.3,750,-1,[],25000,15,1800)  % 2.6 secs
%   >> verifAE(3,0.3,750,-1,[],25000,30,7000)  % 17 secs
%   >> verifAE(3,0.3,750,-1,[],25000,60,25000)  % 2.8 mins
%   >> verifAE(3,0.3,750,-1,[],25000,120,100000)  % 53 mins
%   >> verifAE(3,0.3,750,-1,[],25000,240,400000)  % (~ 16 hours)
%TEST E:  (with basal sliding)
%   >> verifAE(3,0.3,750,200,[200 700 10 40],25000,15,3000)  % 6.3 secs
%   >> verifAE(3,0.3,750,200,[200 700 10 40],25000,30,12000)  % 40 secs
%   >> verifAE(3,0.3,750,200,[200 700 10 40],25000,60,60000)  % 10.8 min
%   >> verifAE(3,0.3,750,200,[200 700 10 40],25000,120,300000)  % 142 min
%(ELB 4/24/04)

clear H
global H dx dy fx fy n2 nm Rx Ry 
global rho g Gam L n M0 Cs 
global etamonut bsflag r1 r2 theta1 theta2 mustgx mustgy

% physical constants, etc
SperA=31556926; % seconds per year (i.e. 365.2422 days)
A=1e-16/SperA;  %=3.17e-24  1/(Pa^3 s); (EISMINT value) flow law parameter
rho=910; % kg/m^3; density of ice
g=9.81; % m/s^2; gravity
n=nGlen; % Glen exponent
Gam=2*(rho*g)^n*A/(n+2); % overall constant in deformation discharge q_f
tf=tfyr*SperA; M0=M0yr/SperA; L=Lkm*1000; % convert to secs and meters
Cs=(2^(n-1)*M0/Gam)^(1/(2*n+2)); % constant in H_s(r)
nut=8.0e9; % Pa s; viscosity of till
etamonut=(etamax/1000)/nut; % constant in mu
H0=getH(0);
bsflag=(etamax>0);
if bsflag, r1=sector(1)*1000; r2=sector(2)*1000;
   theta1=sector(3)*pi/180; theta2=sector(4)*pi/180;
end
errcontours=[-500 -200 -100 -70 -50 -30 -20 -10 -5 -1 ...
             1 5 10 20 30 50 70 100 150 500];
        
% improve display
set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',1.0,...
'defaultlinelinewidth',1.5,'defaultpatchlinewidth',1.2)

% display exact; velocity fields first
figure(1), Nr=20; Nz=20; rvel=linspace(0,L,Nr);
[rrr,zzz]=meshgrid(rvel,linspace(0,H0*1.1,Nz));
HU=getH(rrr);  [dHdrU, discard]=getddr(rrr);
maskU = (zzz <= HU);   U=zeros(size(rrr)); 
U(maskU) = - 2*A*(rho*g)^n/(n+1) * (-dHdrU(maskU)).^n .* ...
    ( (HU(maskU)-zzz(maskU)).^(n+1) - HU(maskU).^(n+1) );
if bsflag % add basal sliding component to velocity
   center=((theta1+theta2)/2);
   Ub=- getmu(rrr,center*ones(size(rrr))).*HU.*dHdrU;
   maskB=maskU & (r1 < rrr) & (rrr < r2);
   U(maskB)=U(maskB)+Ub(maskB);
end
quiver(rrr/1000,zzz,U,zeros(size(U)),.3,'r')
axis([0 L/1000 0 H0*1.1])
text((L/10)/1000,H0*.1,...
   ['Max |U| = ' num2str(max(max(abs(U)))*SperA) ' m/a.'],'Color','r')
% now thickness and accumulation
ep=400;  r=linspace(0,L,ep);
if bsflag, M=M0+getMb(r,center*ones(size(r)));
else M=M0*ones(size(r)); end
hold on
[AX,H1,H2] = plotyy(r/1000,getH(r),r/1000,M*SperA); 
grid on, xlabel('r in km')
if bsflag, title('thickness and accumulation (along centerline of ice stream)'),
else, title('thickness and accumulation'); end
set(get(AX(1),'Ylabel'),'String','H in m (solid)'), set(H1,'LineStyle','-')
set(get(AX(2),'Ylabel'),'String','M in m/a (dotted)'), set(H2,'LineStyle',':')
hold off

% "overhead view" of region of basal sliding
if bsflag
   figure(4)
   [xxp,yyp]=meshgrid(linspace(0,1.05*L,100),linspace(0,1.05*L,100));
   [thp,rrp]=cart2pol(xxp,yyp);
   MMp=M0+getMb(rrp,thp);
   contour(xxp/1000,yyp/1000,MMp*SperA), colorbar
   %to show contours of till thickness:
   %  contour(xxp/1000,yyp/1000,nut*getmucart(xxp,yyp))
   xlabel('x in km'), ylabel('y in km')
   hold on, view(2)
   [xl,yl]=pol2cart([theta1 theta1],[r1 r2]); line(xl/1000,yl/1000)
   [xl,yl]=pol2cart([theta2 theta2],[r1 r2]); line(xl/1000,yl/1000)
   thl=linspace(theta1,theta2,100); 
   [xl,yl]=pol2cart(thl,r1*ones(1,100)); line(xl/1000,yl/1000)
   [xl,yl]=pol2cart(thl,r2*ones(1,100)); line(xl/1000,yl/1000)
   thl=linspace(0,pi/2,200); 
   [xl,yl]=pol2cart(thl,L*ones(1,200)); line(xl/1000,yl/1000,'LineStyle',':')
   text(1.01*xl(100)/1000,1.01*yl(100)/1000,'margin (dotted)')
   title('nonzero basal sliding in sector; contours are of accumulation in m/a')
   hold off
end

% start numeric comparison
if Mt==0, return, end
box=1.1*L;
dx=box/Nx; dy=dx; dt=tf/Mt; Rx=(dt/(dx)^2); Ry=(dt/(dy)^2);
[xx,yy]=ndgrid(linspace(0,box,Nx+1),linspace(0,box,Nx+1)); % grid in space
% ndgrid makes coord sys left-handed; better for computation
[ththeta,rr]=cart2pol(xx,yy);
H=getH(rr); % initial condition is exact profile
t=linspace(0,tf,Mt+1); 
in=2:Nx; fx=4*dx; fy=4*dy; n2=n+2; nm=(n-1)/2;
M=M0*ones(size(H));  % base accumulation rate
outice=(rr>=L); % outside of ice if true
disp(['dx   =   dy       = ' num2str(dx/1000) ' km'])
disp(['dt                = ' num2str(dt/SperA) ' years'])
% volume computation by 2 variable trapezoid, roughly
volc=4*ones(Nx+1,Nx+1); volc(1,:)=2; volc(:,1)=2; volc(Nx+1,:)=2; volc(:,Nx+1)=2;
volc(1,1)=1; volc(Nx+1,1)=1; volc(1,Nx+1)=1; volc(Nx+1,Nx+1)=1; 
V0=dx*dy*sum(sum(volc.*H))/4;
disp(['initial num vol   = ' num2str(V0/1e9) ' cubic km'])

% region of basal sliding; compensatory accumulation
if bsflag
   mustgx=getmucart(xx-dx/2,yy); % staggered grid values
   mustgy=getmucart(xx,yy-dy/2);
   msk=(rr>r1)&(rr<r2)&(ththeta>theta1)&(ththeta<theta2);  % inside sector if true
   % comment out to see effect of sliding w/o compensatory accum:
   M(msk)=M(msk)+getMb(rr(msk),ththeta(msk));
end

% time-stepping loop
wbhandle=waitbar(0,'COMPUTING NUMERICAL APPROXIMATION.  Ctrl-C halts.'); tic
for l=1:Mt
   Hn=zeros(size(H));
   % H(1,:) is edge with x=0; H(:,1) is edge with y=0; H(1,1) is corner (x,y)=(0,0)
   Hn(in,in)=H(in,in) + dt*M(in,in) - divQf(in,in);
   Hn(1,in)=H(1,in) + dt*M(1,in) - divQf(1,in);
   Hn(in,1)=H(in,1) + dt*M(in,1) - divQf(in,1);   
   Hn(1,1)=H(1,1) + dt*M(1,1) - divQf(1,1);

   % comment out to see effect of compensatory accum w/o sliding:
   if bsflag  % if sliding in sector
      Hn(in,in)=Hn(in,in) - divQb(in,in);
      Hn(1,in)=Hn(1,in) - divQb(1,in);
      Hn(in,1)=Hn(in,1) - divQb(in,1);   
      Hn(1,1)=Hn(1,1) - divQb(1,1);
   end

   Hn(outice)=0; % apply boundary condition
   H=Hn;

   % stability diagnostic (uses only interior points; may miss real max by a bit)
   if l==1 
      disp(['init stability    = ' num2str(stabindex(in,in)) ' (= max(D)*dt/dx^2)']), end
   if l==Mt
      disp(['final stability   = ' num2str(stabindex(in,in)) ' (= max(D)*dt/dx^2)']), end
   % try to check for disaster; waitbar
   if rem(l,30)==0
      if max(abs(H(floor(Nx/2),:))) > H0*2
         close(wbhandle), error(['instability (blowup) at step ' int2str(l)]), end
      if H(1,1)<.1*H0
         close(wbhandle), error(['instability (collapse) at step ' int2str(l)]), end
      waitbar(l/(Mt+1)), end
   % if lots of steps, estimate compute time
   if rem(l,500)==0
      remsecs=ceil( (Mt+1-l) * (toc/l) ); close(wbhandle)
      wbhandle=waitbar(0, ['COMPUTING ... ESTIMATED WAIT TIME ' ...
            int2str(remsecs) ' secs']);
      waitbar(l/(Mt+1),wbhandle), end
end;
disp(['actual comp. time = ' num2str(toc) ' secs'])
close(wbhandle)

% plot numerical final state
figure(2), clf
surf(xx/1000,yy/1000,H);
axis([0 box/1000 0 box/1000 0 H0*1.1]); view(90,0)
xlabel('x in km'); ylabel('y in km'); zlabel('h in m');
title(['Numerical final state at t = ' num2str(tf/SperA) ' yrs.  Rotatable 3D fig.'])

% contour of error at final time
figure(3), clf
HHexactf=getH(rr);
err=H-HHexactf;
[Cont,hand] = contour(xx/1000,yy/1000,err,errcontours);
clabel(Cont,hand), axis equal, axis square
xlabel('x in km'); ylabel('y in km');
disp(['max error         = ' num2str(max(max(abs(err)))) ' meters']);
disp(['interior error    = ' num2str(max(max(abs(err(rr<.8*L))))) ' meters ( r < .80 L)']);
disp(['center error      = ' num2str(abs(err(1,1))) ' meters']);
disp(['work (= N^2 M)    = ' int2str(4*Nx*Nx*Mt)]);
Vf=dx*dy*sum(sum(volc.*H))/4;
disp(['final num volume  = ' num2str(Vf/1e9) ' cubic km'])
disp(['volume differenc  = ' num2str((Vf-V0)/1e9) ' cubic km'])


%%%%%%%%%%% exact HELPER FUNCTIONS %%%%%%%%%%%
function H=getH(r)
global n L Cs
ind=(r<L);
H=zeros(size(r));
H(ind)=Cs *( L^(1+1/n) - r(ind).^(1+1/n) ).^(n/(2*n+2));

function [dHs, ddHs]=getddr(r)
% r must be inside ice region
global n L Cs
chi=L^(1+1/n) - r.^(1+1/n);
dHs=-(Cs/2)*r.^(1/n).*chi.^((-n-2)/(2*n+2));
ddHs=-(Cs/(2*n))*chi.^((-3*n-4)/(2*n+2)).*( r.^((1-n)/n).*chi + ((n+2)/2)*r.^(2/n) );

function [mu, dmudr]=getmu(r,theta)
% r,theta must be at grid points inside the sector
global etamonut r1 r2 theta1 theta2
rbot=(r2-r1)^2; thbot=(theta2-theta1)^2;
thfact=4*(theta-theta1).*(theta2-theta)/thbot;
mu=etamonut*(4*(r-r1).*(r2-r)/rbot).*thfact;
dmudr=etamonut*thfact.*(4/rbot).*(r1+r2-2*r);

function mu=getmucart(x,y)
% x,y must be same size and must correspond to same grid pts
global etamonut r1 r2 theta1 theta2
[theta,r]=cart2pol(x,y);
ind=(r>r1)&(r<r2)&(theta>theta1)&(theta<theta2);
rbot=(r2-r1)^2; thbot=(theta2-theta1)^2;
mu=zeros(size(x));
mu(ind)=etamonut*(4*(r(ind)-r1).*(r2-r(ind))/rbot).*...
         (4*(theta(ind)-theta1).*(theta2-theta(ind))/thbot);

function Mb=getMb(r,theta)
% r,theta must be same size and must correspond to same grid pts
global n M0 Gam L rho g r1 r2 theta1 theta2
ind=(r>r1)&(r<r2)&(theta>theta1)&(theta<theta2);
Hs=getH(r(ind)); [dHs, ddHs]=getddr(r(ind));
[mu, dmudr]=getmu(r(ind),theta(ind));
Mb=zeros(size(r));
Mb(ind)=-rho*g*( (Hs.^2).*dHs.*((mu./r(ind))+dmudr) + mu.*Hs.*(2*dHs.^2+Hs.*ddHs) );


%%%%%%%%%%% numerical HELPER FUNCTIONS %%%%%%%%%%%
function dQ=divQf(ix,iy)
% (numerical) Divergence of deformation (flow) flux
global H dx dy fx fy n2 nm Rx Ry Gam
if length(ix)==1, px=2; mx=2; else, px=ix+1; mx=ix-1; end
if length(iy)==1, py=2; my=2; else, py=iy+1; my=iy-1; end
Hin=H(ix,iy);
Hbr=(H(px,iy) + Hin)/2; Hbl=(Hin + H(mx,iy))/2;
Hbu=(H(ix,py) + Hin)/2; Hbd=(Hin + H(ix,my))/2;
dHsr=((H(px,iy)-Hin)/dx).^2 + ((H(px,py)+H(ix,py)-H(px,my)-H(ix,my))/fy).^2;
dHsl=((Hin-H(mx,iy))/dx).^2 + ((H(mx,py)+H(ix,py)-H(mx,my)-H(ix,my))/fy).^2;
dHsu=((H(px,py)+H(px,iy)-H(mx,py)-H(mx,iy))/fx).^2 + ((H(ix,py)-Hin)/dy).^2;
dHsd=((H(px,iy)+H(px,my)-H(mx,iy)-H(mx,my))/fx).^2 + ((Hin-H(ix,my))/dy).^2;
dQx=(Hbr.^n2.*dHsr.^nm).*(H(px,iy)-Hin) - (Hbl.^n2.*dHsl.^nm).*(Hin-H(mx,iy));
dQy=(Hbu.^n2.*dHsu.^nm).*(H(ix,py)-Hin) - (Hbd.^n2.*dHsd.^nm).*(Hin-H(ix,my));
dQ=-Gam*(Rx*dQx + Ry*dQy);

function dQ=divQb(ix,iy)
% (numerical) Divergence of basal sliding flux
global H Rx Ry rho g mustgx mustgy
if length(ix)==1, px=2; mx=2; else, px=ix+1; mx=ix-1; end
if length(iy)==1, py=2; my=2; else, py=iy+1; my=iy-1; end
Hin=H(ix,iy);
Hbr=(H(px,iy) + Hin)/2; Hbl=(Hin + H(mx,iy))/2;
Hbu=(H(ix,py) + Hin)/2; Hbd=(Hin + H(ix,my))/2;
dQx=mustgx(px,iy).*Hbr.^2.*(H(px,iy)-Hin) - mustgx(ix,iy).*Hbl.^2.*(Hin-H(mx,iy));
dQy=mustgy(ix,py).*Hbu.^2.*(H(ix,py)-Hin) - mustgy(ix,iy).*Hbd.^2.*(Hin-H(ix,my));
dQ=-rho*g*(Rx*dQx + Ry*dQy);

function RR=stabindex(ix,iy)
% compute certain diffusivities for max estimate
global H dx dy fx fy n2 nm Rx Ry Gam rho g bsflag mustgx
px=ix+1; py=iy+1; my=iy-1;
Hin=H(ix,iy); Hbr=(H(px,iy) + Hin)/2;
dHsr=((H(px,iy)-Hin)/dx).^2 + ((H(px,py)+H(ix,py)-H(px,my)-H(ix,my))/fy).^2;
Df=Gam*max(max( Hbr.^n2.*dHsr.^nm )); RR=Rx*Df;
if bsflag
   Db=rho*g*max(max( mustgx(ix,iy).*Hbr.^2 )); RR=RR+Rx*Db;
end

