function verifD(Cp, Tpyr, H0, Lkm, tfyr, Nx, Mt);
%VERIFD  Compares numerical to exact perturbed steady solution
%    of shallow ice equation
%       H_t = (Ms + Mc) - Div q_f,
%    incorporating compensatory accumulation.  Numerical method: type I explicit 
%    finite difference.  Takes advantage of radial symmetry to do 1/4 of work.
%    Special to  n=3.  Computation is on [0,1.2 L]^2.  Accumulation function 
%    extended by -0.1 m/a beyond margin; ice flow outside  r=L  indicates error.
%
%verifD(Cp, Tp, H0, L, tf, Nx, Mt);
%   Cp    = magnitude of perturbation (meters)
%   Tp    = period of perturbation (years)
%   H0    = thickness at center (meters)
%   L     = margin radius (km)
%   tf    = run from 0 to tf (years)
%   Nx    = number of grid intervals in both x and y
%   Mt    = number of time steps (dt=(tf-t0)/Mt=tadv/Mt)
%
%Notes: 
%   (1) Stability index  max(D)*dt/dx^2  of roughly 0.06 is stability limit.
%       Code uses fixed time step but this index could be used to adapt.
%   (2) Displays in figures 1,2,3,4 and 5.
%   (3) Reference: Bueler et al (2004), "Exact solutions and the verification 
%       of numerical models for isothermal ice sheets", preprint.
%
%Examples:   (all near stability limit)
%TEST D:  
%   >> verifD(200,5000,3600,750,25000,15,5000)  % 16 secs
%   >> verifD(200,5000,3600,750,25000,30,20000)  % 153 secs
%   >> verifD(200,5000,3600,750,25000,60,100000)  % 45 mins
%   >> verifD(200,5000,3600,750,25000,120,400000)  % 10.6 hours
%(ELB 4/24/04)

clear H
global H dx dy fx fy n2 nm Rx Ry Gam SperA
global L n H0G CpG Tp C

% physical constants, etc
SperA=31556926; % seconds per year (i.e. 365.2422 days)
A=1e-16/SperA;  %=3.17e-24  1/(Pa^3 s); (EISMINT value) flow law parameter
rho=910; % kg/m^3; density of ice
g=9.81; % m/s^2; gravity
n=3; % Glen exponent
Gam=2*(rho*g)^n*A/(n+2); % overall constant in deformation discharge q_f
L=Lkm*1000;
C=(Gam*H0^(2*n+2))/(2*L*(1-1/n))^n;
tf=tfyr*SperA; Tp=Tpyr*SperA;
t=linspace(0,tf,Mt+1);
H0G=H0; CpG=Cp; % necessary but irritating way to make variables global

% improve display
set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',1.0,...
'defaultlinelinewidth',1.5,'defaultpatchlinewidth',1.2)

% display exact
box=1.2*L;
ep=500;  r=linspace(0,box,ep);
figure(1), clf % thickness and accumulation at t0, tf
subplot(2,1,1)
plot(r/1000,getHs(r),r/1000,getHp(r,tf)), grid on
legend('thickness at t=0','thickness at t=tf'), ylabel('H in m')
title(['conditions at t=0 and t=' num2str(tfyr) ' years'])
subplot(2,1,2)
plot(r/1000,(getMs(r)+getMc(r,0))*SperA,r/1000,(getMs(r)+getMc(r,tf))*SperA), grid on
legend('M=M_s+M_c at t=0','M=M_s+M_c at t=t_f')
ylabel('M in m/a'), xlabel('r in km')

figure(2), clf % thickness and accumulation envelopes
subplot(2,1,1)
envelopes=8; profiles=zeros(envelopes,ep);
for k=1:envelopes,  profiles(k,:)=getHp(r,(k-1)*Tp/envelopes); end
plot(r/1000,profiles','k'), grid on
ylabel('H in m')
title(['envelopes (values at ' int2str(envelopes) ' times over one cycle)'])
subplot(2,1,2)
profiles=zeros(envelopes,ep); Mstemp=getMs(r);
for k=1:envelopes,  accums(k,:)=Mstemp+getMc(r,(k-1)*Tp/envelopes); end
plot(r/1000,SperA*accums','k'), grid on
ylabel('M in m/a'), xlabel('r in km')

% start numeric comparison
if Mt==0, return, end
dx=box/Nx; dy=dx;
[xx,yy]=ndgrid(linspace(0,box,Nx+1),linspace(0,box,Nx+1)); % grid in space
% ndgrid makes coord sys left-handed; better for computation
rr=sqrt(xx.^2+yy.^2);
H=getHs(rr); % initial condition  Hp(r,0)=Hs(r)
dt=tf/Mt; in=2:Nx; Rx=(dt/(dx)^2); Ry=(dt/(dy)^2);
fx=4*dx; fy=4*dy; n2=n+2; nm=(n-1)/2;
disp(['dx   =   dy       = ' num2str(dx/1000) ' km'])
disp(['dt                = ' num2str(dt/SperA) ' years'])
% volume computation by 2 variable trapezoid, roughly
volc=4*ones(Nx+1,Nx+1); volc(1,:)=2; volc(:,1)=2; volc(Nx+1,:)=2; volc(:,Nx+1)=2;
volc(1,1)=1; volc(Nx+1,1)=1; volc(1,Nx+1)=1; volc(Nx+1,Nx+1)=1; 
V0=dx*dy*sum(sum(volc.*H))/4;
disp(['initial num vol   = ' num2str(V0/1e9) ' cubic km'])

% time-stepping loop
wbhandle=waitbar(0,'COMPUTING NUMERICAL APPROXIMATION.  Ctrl-C halts.');
tic
Msteady=getMs(rr);
for l=1:Mt
   Hn=zeros(size(H));
   M=Msteady+getMc(rr,(l-1)*dt); 

   % H(1,:) is edge with x=0; H(:,1) is edge with y=0; H(1,1) is corner (x,y)=(0,0)
   Hn(in,in)=H(in,in) + dt*M(in,in) - divQf(in,in);
   Hn(1,in)=H(1,in) + dt*M(1,in) - divQf(1,in);
   Hn(in,1)=H(in,1) + dt*M(in,1) - divQf(in,1);   
   Hn(1,1)=H(1,1) + dt*M(1,1) - divQf(1,1);
   H=max(0,Hn); % apply boundary condition
   
   % stability diagnostic (uses only interior points; may miss real max by a bit)
   if l==1 
      disp(['init stability    = ' num2str(stabindex(in,in)) ' (= max(D)*dt/dx^2)']), end
   if l==Mt
      disp(['final stability   = ' num2str(stabindex(in,in)) ' (= max(D)*dt/dx^2)']), end
   % occasional check for disaster; waitbar
   if rem(l,30)==0
      if max(abs(H(floor(Nx/2),:))) > H0*2
         close(wbhandle), error(['instability (blowup) failure at step ' ...
               int2str(l)]), end
      waitbar(l/(Mt+1)), end
   % if lots of steps, estimate compute time
   if rem(l,500)==0
      remsecs=ceil( (Mt+1-l) * (toc/l) );
      waitbar(l/(Mt+1),wbhandle, ['COMPUTING ... ESTIMATED WAIT TIME ' ...
            int2str(remsecs) ' secs']), end
end;
disp(['actual comp. time = ' num2str(toc) ' secs'])
close(wbhandle)

% plot numerical final state
figure(3), clf
surf(xx/1000,yy/1000,H);
axis([0 box/1000 0 box/1000 0 H0*1.1]); view(90,0)
xlabel('x in km'); ylabel('y in km'); zlabel('h in m');
title(['Numerical final state at t = ' num2str(tf/SperA) ' yrs.  Rotatable 3D fig.'])

% contour of error at final time
figure(4), clf
HHexactf=getHp(rr,tf);
err=H-HHexactf;
[Cont,hand] = contour(xx/1000,yy/1000,err,...
   [-500 -100 -50 -20 -10 -5 -1 1 5 10 20 50 100 500]);
clabel(Cont,hand), axis square
xlabel('x in km'); ylabel('y in km');
disp(['max err (at tf)   = ' num2str(max(max(err))) ' meters']);
disp(['center err (tf)   = ' num2str(err(1,1)) ' meters']);
disp(['work (= N^2 M)    = ' int2str(4*Nx*Nx*Mt)]);
Vf=dx*dy*sum(sum(volc.*H))/4;
disp(['final num volume  = ' num2str(Vf/1e9) ' cubic km'])
disp(['volume differenc  = ' num2str((Vf-V0)/1e9) ' cubic km'])

% show margin "mislocation"
figure(5), spy(flipud(xor(H>0,HHexactf>0))), xlabel('')
title('margin mislocation grid points'), set(gca,'XTick',[],'YTick',[])


%%%%%%%%%%% exact HELPER FUNCTIONS %%%%%%%%%%%
function Hp=getHp(r,t)
global CpG Tp
Hp=getHs(r)+CpG*sin(2*pi*t/Tp)*getgp(r);

function Hs=getHs(r)
global n L H0G
ind=(r<L); chi=zeros(size(r));
chi(ind)=(1+1/n)*(r(ind)/L) - 1/n + (1-(r(ind)/L)).^(1+1/n) - (r(ind)/L).^(1+1/n);
Hs=( H0G/(1-1/n)^(n/(2*n+2)) ) * chi.^(n/(2*n+2));

function Mc=getMc(r,t)
global L Gam H0G CpG Tp
% special to n=3
ind=(r>=.3*L)&(r<=.9*L);  rr=r(ind); % needs to be the case that rr==0 is empty
Hp=getHp(rr,t);  [dHp,ddHp]=getddr(rr,t);
divterms=Gam*Hp.^4.*dHp.^2.* ( (1./rr).*Hp.*dHp + 5*dHp.^2 + 3*Hp.*ddHp );
Mc=zeros(size(r));
Mc(ind)=(2*pi*CpG/Tp)*cos(2*pi*t/Tp)*getgp(rr)-getMs(rr)-divterms;

function Ms=getMs(r)
global n L Gam C SperA
ind=(abs(r)>5*eps)&(r<L); Ms=zeros(size(r));
temp=(r(ind)/L).^(1/n)+(1-r(ind)/L).^(1/n)-1;
Ms(ind)=(C./r(ind)).*temp.^(n-1).*  ...
   ( 2*(r(ind)/L).^(1/n)+(1-r(ind)/L).^(1/n-1).*(1-2*r(ind)/L)-1 );
Ms((abs(r)<10*eps))=2*C/L; % Ms(r=0) found by limit argument
Ms(abs(r)>=L)=-0.1/SperA;

function gp=getgp(r)
global L
ind=(r>=.3*L)&(r<=.9*L);  gp=zeros(size(r));
gp(ind)=.5*cos(pi*(r(ind)-.6*L)/(.3*L))+.5;

function [dHp, ddHp]=getddr(r,t)
%special to n=3; assumes .3 L < r < .9 L
global L Gam H0G CpG Tp
gp=zeros(size(r)); dgp=gp; chi=gp;  dchi=gp; ddchi=gp;
gp=.5*cos(pi*(r-.6*L)/(.3*L))+.5;
dgp=(-pi/(.6*L))*sin(pi*(r-.6*L)/(.3*L));
ddgp=(-pi^2/(.18*L^2))*cos(pi*(r-.6*L)/(.3*L));
chi=4*r/(3*L) - 1/3 + (1-r/L).^(4/3) - (r/L).^(4/3);
dchi=(-4/(3*L))*( (r/L).^(1/3)+(1-r/L).^(1/3)-1 );
ddchi=(-4/(9*L^2))*( (r/L).^(-2/3) - (1-(r/L)).^(-2/3) );
c1=(3*H0G)/(8*(2/3)^(3/8));
dHs=c1 * chi.^(-5/8) .* dchi;
ddHs=c1 *( (-5/8)*chi.^(-13/8).*dchi.^2 + chi.^(-5/8).*ddchi );
dHp=dHs + CpG*sin(2*pi*t/Tp)*dgp;
ddHp=ddHs + CpG*sin(2*pi*t/Tp)*ddgp;


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

function RR=stabindex(ix,iy)
% compute certain diffusivities for max estimate
global H dx dy fx fy n2 nm Rx Ry Gam
px=ix+1; py=iy+1; my=iy-1;
Hin=H(ix,iy); Hbr=(H(px,iy) + Hin)/2;
dHsr=((H(px,iy)-Hin)/dx).^2 + ((H(px,py)+H(ix,py)-H(px,my)-H(ix,my))/fy).^2;
RR=Rx*Gam*max(max( Hbr.^n2.*dHsr.^nm ));;

