function verifBu(H0, R0km, tDelyr, Nx, Mt);
%VERIFBU  Compares numerical to exact radial similarity solution
%       H(r,t)=t^{-alpha} phi(t^{-beta} r), 
%    for lam=0 and n=3 of shallow ice equation
%       H_t = M - Div q_f,
%    incorporating nonzero accumulation.  Two numerical methods; both are
%    explicit finite difference, but first is in H and second in u=H^(8/3). 
%    Takes advantage of radial symmetry to do 1/4 of work.  Automatically 
%    determines horizontal domain.
%
%verifBu(H0, R0, tDel, Nx, Mt);
%   H0    = central thickness at t0 (meters)
%   R0    = margin radius at t0 (km)
%   tDel  = time to advance from t0 (tf = t0 + tDel) (years)
%   Nx    = number of grid intervals in both x and y
%   Mt    = number of time steps (dt=(tf-t0)/Mt=tadv/Mt)
%
%Notes: 
%   (1) Numerical computation if Mt>0; otherwise shows exact states only.
%   (2) Requires tDel>0; run is from t=t0 to tf=t0+tDel.
%   (3) Initial stability index  max(D)*dt/dx^2  of roughly 0.13 is limit.
%       Code uses fixed time step but this index could be used to adapt.
%   (4) Displays in figures 1,2,3,4,5 and 6.
%   (6) Reference: Bueler et al (2004), "Exact solutions and the verification
%       of numerical models for isothermal ice sheets", preprint.
%
%Examples:  (all near stability limit)
%TEST B done in H and in u:
%   >> verifBu(3600,750,25000,15,1000)  % 3 secs
%   >> verifBu(3600,750,25000,30,6000)  % 32 secs
%   >> verifBu(3600,750,25000,60,30000)  % 7.3 mins
%   >> verifBu(3600,750,25000,120,120000)  % 2.2 hours
%(ELB 4/24/04)

clear H u
global H u dx dy fx fy n2 nm Rx Ry Gam tilGam

% physical constants
SperA=31556926; % seconds per year (i.e. 365.2422 days)
A=1e-16/SperA;  %=3.17e-24  1/(Pa^3 s); (EISMINT value) flow law parameter
rho=910; % kg/m^3; density of ice
g=9.81; % m/s^2; gravity
n=3;
Gam=2*(rho*g)^n*A/(n+2); % overall constant in deformation discharge q_f
tilGam=(n/(2*n+2))^n*Gam;
errconts=[-500 -100 -20 -5 -1 0 1 5 20 100 500];

% improve display
set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',1.0,...
'defaultlinelinewidth',1.5,'defaultpatchlinewidth',1.2)

% constants in sim soln
lam=0;
alf=(2-(n+1)*lam)/(5*n+3);
bet=(1+(2*n+1)*lam)/(5*n+3); 

% time since creation (typically only thousands of years despite Darwin)
R0=R0km*1000;
t0 = (bet/Gam) * ( (2*n+1)/((n+1)) )^n * (R0^(n+1)/H0^(2*n+1));
if tDelyr<=0, error('tDel must be positive'), end
tDel=tDelyr*SperA;  tf=t0+tDel;  t=linspace(t0,tf,max(Mt,1)+1);

% internal constants; time grid; max dimensions
s0=t0^(-bet)*R0;
Rmax=tf^bet * s0; % margin at last time
L=Rmax*1.1; % domain: (x,y) in [0,L] x [0,L]

% draw exact initial and final profiles
Ner=1000; r=linspace(0,L,Ner); % for display; not numerics
Hexacti=getH(n,alf,bet,H0,R0,t0,t0,r);
Hexactf=getH(n,alf,bet,H0,R0,t0,tf,r);

figure(1); clf, set(gcf,'DefaultLineLineWidth',1.5)
plot(r/1000,Hexacti,r/1000,Hexactf);  
legend(['t_0 = ' num2str(t0/SperA) ' a'],['t_f = ' num2str(tf/SperA) ' a'])
hold on, ylabel('h in m'); grid on, axis([0 L/1000 0 H0*1.1]);
xlabel('r in km'); title('Exact profiles.');

figure(2); clf, set(gcf,'DefaultLineLineWidth',1.5)
plot(r/1000,Hexacti.^(8/3),'--',r/1000,Hexactf.^(8/3));  
legend(['t_0 = ' num2str(t0/SperA) ' a'],['t_f = ' num2str(tf/SperA) ' a'])
hold on, ylabel('u in m^{8/3}'); grid on
set(gca,'XLim',[0 L/1000])
xlabel('r in km'); title('Exact profiles *of u*.');

% start numeric comparison
if Mt==0, return, end
dx=L/Nx; dy=dx;
[xx,yy]=ndgrid(linspace(0,L,Nx+1),linspace(0,L,Nx+1)); % grid in space
% ndgrid makes coord sys left-handed; better for computation
rr=sqrt(xx.^2+yy.^2);
H=getH(n,alf,bet,H0,R0,t0,t0,rr); dt=tDel/Mt; % initial condition
u=H.^(8/3); % initial u
in=2:Nx; Rx=(dt/(dx)^2); Ry=(dt/(dy)^2);
fx=4*dx; fy=4*dy; n2=n+2; nm=(n-1)/2;
disp(['t0                = ' num2str(t0/SperA) ' years (time since delta mass)'])
disp(['Rmax              = ' num2str(Rmax/1000) ' km'])
disp(['dx   =   dy       = ' num2str(dx/1000) ' km'])
disp(['dt                = ' num2str(dt/SperA) ' years'])

% time-stepping loop
wbhandle=waitbar(0,'COMPUTING NUMERICAL APPROXIMATION.  Ctrl-C halts.'); tic
for l=1:Mt
   Hn=zeros(size(H)); 
   % H(1,:) is edge with x=0; H(:,1) is edge with y=0; H(1,1) is (x,y)=(0,0)
   Hn(in,in)=H(in,in) - divQf(in,in);
   Hn(1,in)=H(1,in) - divQf(1,in);
   Hn(in,1)=H(in,1) - divQf(in,1);   
   Hn(1,1)=H(1,1) - divQf(1,1);
   H=max(0,Hn); % apply boundary condition
   
   un=zeros(size(u)); 
   un(in,in)=u(in,in).^(3/8) - divuQf(in,in);
   un(1,in)=u(1,in).^(3/8) - divuQf(1,in);
   un(in,1)=u(in,1).^(3/8) - divuQf(in,1);   
   un(1,1)=u(1,1).^(3/8) - divuQf(1,1);
   un=max(0,un); % apply boundary condition
   u=un.^(8/3);
   
   % stability diagnostic (uses only interior points; may miss real max by a bit)
   if l==1 
      disp(['init stability    = ' num2str(stabindex(in,in))...
            ' (= max(D)*dt/dx^2)']), end
   if l==Mt
      disp(['final stability   = ' num2str(stabindex(in,in))...
            ' (= max(D)*dt/dx^2)']), end
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

% contour of error at final time
figure(3), clf
HHexactf=getH(n,alf,bet,H0,R0,t0,tf,rr);
err=HHexactf-H;
[Cont,hand] = contour(xx/1000,yy/1000,err,errconts);
clabel(Cont,hand), axis equal, axis square
xlabel('x in km'); ylabel('y in km');
title('Direct computation.')
disp('direct computation with H errors:')
disp(['max err (at tf)   = ' num2str(max(max(err))) ' meters']);
disp(['center err (tf)   = ' num2str(err(1,1)) ' meters']);

% contour of error at final time: via u
figure(4), clf
err=HHexactf-u.^(3/8);
[Cont,hand] = contour(xx/1000,yy/1000,err,errconts);
clabel(Cont,hand), axis equal, axis square
clabel(Cont,hand), axis equal, axis square
xlabel('x in km'); ylabel('y in km');
title('Computation via u');
disp('errors with indirect computation (in terms of u):')
disp(['max err (at tf)   = ' num2str(max(max(err))) ' meters']);
disp(['center err (tf)   = ' num2str(err(1,1)) ' meters']);

% contour of error in u at final time
figure(5), clf
err=(HHexactf.^(8/3)-u)/u(1,1);
[Cont,hand] = contour(xx/1000,yy/1000,abs(err),...
   [0 .0002 .0005 .001 .002 .005 .010 .020 .050 .100]);
clabel(Cont,hand), axis equal, axis square
xlabel('x in km'); ylabel('y in km');
title('Relative errors in u');
disp('relative error in u:')
disp(['max err (at tf)   = ' num2str(max(max(err)))]);
disp(['center err (tf)   = ' num2str(err(1,1))]);

% show margin "mislocation"
figure(6), spy(flipud(xor(u>0,HHexactf>0))), xlabel('')
title('margin mislocation grid points'), set(gca,'XTick',[],'YTick',[])


%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%
function HOUT=getH(n,alf,bet,H0,R0,t0,t,r)
% compute exact H
rscl=( t^(-bet)*r )/( t0^(-bet)*R0 );
temp=max(0, 1-rscl.^((n+1)/n) );
HOUT=H0*(t/t0)^(-alf)*temp.^(n/(2*n+1));

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

function dQ=divuQf(ix,iy)
% (numerical) Divergence of deformation (flow) flux *in terms of u*
global u dx dy fx fy nm Rx Ry tilGam
if length(ix)==1, px=2; mx=2; else, px=ix+1; mx=ix-1; end
if length(iy)==1, py=2; my=2; else, py=iy+1; my=iy-1; end
uin=u(ix,iy);
dusr=((u(px,iy)-uin)/dx).^2 + ((u(px,py)+u(ix,py)-u(px,my)-u(ix,my))/fy).^2;
dusl=((uin-u(mx,iy))/dx).^2 + ((u(mx,py)+u(ix,py)-u(mx,my)-u(ix,my))/fy).^2;
dusu=((u(px,py)+u(px,iy)-u(mx,py)-u(mx,iy))/fx).^2 + ((u(ix,py)-uin)/dy).^2;
dusd=((u(px,iy)+u(px,my)-u(mx,iy)-u(mx,my))/fx).^2 + ((uin-u(ix,my))/dy).^2;
dQx=(dusr.^nm).*(u(px,iy)-uin) - (dusl.^nm).*(uin-u(mx,iy));
dQy=(dusu.^nm).*(u(ix,py)-uin) - (dusd.^nm).*(uin-u(ix,my));
dQ=-tilGam*(Rx*dQx + Ry*dQy);

function RR=stabindex(ix,iy)
% compute certain diffusivities for max estimate
global H dx dy fx fy n2 nm Rx Ry Gam
px=ix+1; py=iy+1; my=iy-1;
Hin=H(ix,iy); Hbr=(H(px,iy) + Hin)/2;
dHsr=((H(px,iy)-Hin)/dx).^2 + ((H(px,py)+H(ix,py)-H(px,my)-H(ix,my))/fy).^2;
Df=Gam*max(max( Hbr.^n2.*dHsr.^nm ));
RR=Rx*Df;

