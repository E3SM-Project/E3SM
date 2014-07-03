function verifBC(n, lam, H0, R0km, tDelyr, Nx, Mt);
%VERIFBC  Compares numerical to exact radial similarity solution
%       H(r,t)=t^{-alpha} phi(t^{-beta} r), 
%    of shallow ice equation
%       H_t = M - Div q_f,
%    incorporating nonzero accumulation.  Numerical method: type I explicit 
%    finite difference.  Takes advantage of radial symmetry to do 1/4 of work.
%    Automatically determines horizontal domain.
%
%verifBC(n, lam, H0, R0, tDel, Nx, Mt);
%   n     = Glen exponent (n >= 1)
%   lam   = accumulation parameter (M_lam = lam t^-1 H)
%   H0    = central thickness at t0 (meters)
%   R0    = margin radius at t0 (km)
%   tDel  = time to advance from t0 (tf = t0 + tDel) (years)
%   Nx    = number of grid intervals in both x and y
%   Mt    = number of time steps (dt=(tf-t0)/Mt=tadv/Mt)
%
%Notes: 
%   (1) Numerical computation if Mt>0; otherwise shows exact states only.
%   (2) If tDel<0 then run from t=0 to t=t0 using zero initial condition; 
%       only advisable if lam>=2/(n+1).  If tDel>0 then run from t=t0 to tf=t0+tDel.
%   (3) Stability index  max(D)*dt/dx^2  of roughly 0.15 is stability limit.
%       Code uses fixed time step but this index could be used to adapt.
%   (4) Displays in figures 1,2,3, and 4.
%   (6) Reference: Bueler et al (2004), "Exact solutions and the verification 
%       of numerical models for isothermal ice sheets", preprint.
%
%Examples:   (all near stability limit)  
%TEST B:  (Halfar solution)
%   >> verifBC(3,0,3600,750,25000,15,1000)  % 1.4 sec
%   >> verifBC(3,0,3600,750,25000,30,6000)  % 16 secs
%   >> verifBC(3,0,3600,750,25000,60,30000)  % 232 secs
%   >> verifBC(3,0,3600,750,25000,120,120000)  % 85 min
%TEST C:  (grows from zero initial condition)
%   >> verifBC(3,5,3600,750,-1,15,3000)   % 5 secs
%   >> verifBC(3,5,3600,750,-1,30,20000)   % 77 secs
%   >> verifBC(3,5,3600,750,-1,60,150000)   % 24 minutes
%   >> verifBC(3,5,3600,750,-1,120,600000)   % 8.8 hours
%(ELB 4/24/04)

clear H
global H dx dy fx fy n2 nm Rx Ry Gam

% physical constants
SperA=31556926; % seconds per year (i.e. 365.2422 days)
A=1e-16/SperA;  %=3.17e-24  1/(Pa^3 s); (EISMINT value) flow law parameter
rho=910; % kg/m^3; density of ice
g=9.81; % m/s^2; gravity
Gam=2*(rho*g)^n*A/(n+2); % overall constant in deformation discharge q_f

% improve display
set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',1.0,...
'defaultlinelinewidth',1.5,'defaultpatchlinewidth',1.2)

% constants in sim soln
alf=(2-(n+1)*lam)/(5*n+3);
bet=(1+(2*n+1)*lam)/(5*n+3); 

% time since creation (typically only thousands of years despite Darwin)
R0=R0km*1000;
t0 = (bet/Gam) * ( (2*n+1)/((n+1)) )^n * (R0^(n+1)/H0^(2*n+1));
tDel=tDelyr*SperA;
if tDel<0
   tf=t0;
   t=linspace(0,t0,max(Mt,1)+1);
else
   if tDel==0, tDel=2*t0; Mt=0; end
   tf=t0+tDel;
   t=linspace(t0,tf,max(Mt,1)+1);
end

% internal constants; time grid; max dimensions
s0=t0^(-bet)*R0;
lamflag=(lam~=0);
Rmax=tf^bet * s0; % margin at last time
L=Rmax*1.1; % domain: (x,y) in [0,L] x [0,L]

% draw exact initial and final profiles
Ner=300;
r=linspace(0,L,Ner); % for display; not numerics
if tDel<0, Hexacti=zeros(size(r));
else, Hexacti=getH(n,alf,bet,H0,R0,t0,t0,r); end
Hexactf=getH(n,alf,bet,H0,R0,t0,tf,r);
figure(1); clf, set(gcf,'DefaultLineLineWidth',1.5)
if lamflag, subplot(3,1,1:2), end
plot(r/1000,Hexacti,r/1000,Hexactf);  
if tDel<0, legend(['t_0 = 0 a'], ['t_f = ' num2str(tf/SperA) ' a']),
else
   legend(['t_0 = ' num2str(t0/SperA) ' a'],['t_f = ' num2str(tf/SperA) ' a'])
end
hold on
ylabel('h in m'); grid on, axis([0 L/1000 0 H0*1.1]);
titletext='Exact profiles.  Horizontal velocities at t_f.';
if lamflag, titletext=[titletext '  Exact accumulations below.'];
else, xlabel('r in km'); end
title(titletext);

% fill in velocity fields; final only
Nr=20; Nz=20;
[rrr,zzz]=meshgrid(linspace(0,L,Nr),linspace(0,H0*1.1,Nz));
HU=getH(n,alf,bet,H0,R0,t0,tf,rrr);
sss=tf^(-bet)*rrr;
dHdrU=-H0*(tf/t0)^(-alf) * ((n+1)/(2*n+1)) * ( tf^(-bet)/(s0^((n+1)/n)) ) * ...
   sss.^(1/n) .* max(0, 1-(sss/s0).^((n+1)/n) ).^(-(n+1)/(2*n+1));
maskU = (zzz <= HU) & (zzz >= 0) & (HU > 0);
U=zeros(size(rrr)); 
U(maskU) = - 2*A*(rho*g)^n/(n+1) * (-dHdrU(maskU)).^n .* ...
    ( (HU(maskU)-zzz(maskU)).^(n+1) - HU(maskU).^(n+1) );
set(gcf,'DefaultLineLineWidth',.5)
quiver(rrr/1000,zzz,U,zeros(size(U)),.3,'r')
text((L/10)/1000,H0*.1,...
   ['(Max U at t_f) = ' num2str(max(max(abs(U)))*SperA) ' m/a.'],'Color','r')
hold off

%display accumulation if appropriate
if lamflag
   for tt=[t0 tf]
      HM=getH(n,alf,bet,H0,R0,t0,tt,r);
      MM=lam*tt^(-1)*HM;
      if tt==t0, Mp1=MM; end
      if tt==tf, Mp2=MM; end
   end
   subplot(3,1,3), set(gcf,'DefaultLineLineWidth',1.5)
   [axeshand,junk1,junk2] = plotyy(r/1000,Mp1*SperA,r/1000,Mp2*SperA); 
   if tDel>=0, set(get(axeshand(1),'Ylabel'),'String','M(t_0) in m/a'), end
   set(get(axeshand(2),'Ylabel'),'String','M(t_f) in m/a')
   grid on, xlabel('r in km')
end
drawnow

% start numeric comparison
if Mt==0, return, end
dx=L/Nx; dy=dx;
[xx,yy]=ndgrid(linspace(0,L,Nx+1),linspace(0,L,Nx+1)); % grid in space
% ndgrid makes coord sys left-handed; better for computation
rr=sqrt(xx.^2+yy.^2);
if tDel<0, H=zeros(size(rr)); dt=t0/Mt;
else,  H=getH(n,alf,bet,H0,R0,t0,t0,rr); dt=tDel/Mt; end % initial condition
in=2:Nx; clear temp
Rx=(dt/(dx)^2); Ry=(dt/(dy)^2);
fx=4*dx; fy=4*dy; n2=n+2; nm=(n-1)/2;
disp(['t0                = ' num2str(t0/SperA) ' years (time since creation as delta mass)'])
disp(['Rmax              = ' num2str(Rmax/1000) ' km'])
disp(['dx   =   dy       = ' num2str(dx/1000) ' km'])
disp(['dt                = ' num2str(dt/SperA) ' years'])
% volume computation by 2 variable trapezoid, roughly
volc=4*ones(Nx+1,Nx+1); volc(1,:)=2; volc(:,1)=2; volc(Nx+1,:)=2; volc(:,Nx+1)=2;
volc(1,1)=1; volc(Nx+1,1)=1; volc(1,Nx+1)=1; volc(Nx+1,Nx+1)=1; 
V0=dx*dy*sum(sum(volc.*H))/4;
disp(['initial num vol   = ' num2str(V0/1e9,30) ' cubic km'])

% time-stepping loop
wbhandle=waitbar(0,'COMPUTING NUMERICAL APPROXIMATION.  Ctrl-C halts.');
tic
for l=1:Mt
   Hn=zeros(size(H));  M=Hn; 
   if lamflag
      if (tDel<0)&(t(l)<1), M=zeros(size(H));
      else, M=lam*t(l)^(-1)*getH(n,alf,bet,H0,R0,t0,t(l),rr); end
   end
   % H(1,:) is edge with x=0; H(:,1) is edge with y=0; H(1,1) is corner (x,y)=(0,0)
   Hn(in,in)=H(in,in) - divQf(in,in);
   Hn(1,in)=H(1,in) - divQf(1,in);
   Hn(in,1)=H(in,1) - divQf(in,1);   
   Hn(1,1)=H(1,1) - divQf(1,1);
   if lamflag
      Hn(in,in)=Hn(in,in) + dt*M(in,in);
      Hn(1,in)=Hn(1,in) + dt*M(1,in);
      Hn(in,1)=Hn(in,1) + dt*M(in,1);   
      Hn(1,1)=Hn(1,1) + dt*M(1,1);
   end
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
figure(2), clf
surf(xx/1000,yy/1000,H);
axis([0 L/1000 0 L/1000 0 H0*1.1]); view(90,0)
xlabel('x in km'); ylabel('y in km'); zlabel('h in m');
title(['Numerical final state at t = ' num2str(tf/SperA) ' yrs.  Rotatable 3D fig.'])

% contour of error at final time
figure(3), clf
HHexactf=getH(n,alf,bet,H0,R0,t0,tf,rr);
err=abs(HHexactf-H);
[Cont,hand] = contour(xx/1000,yy/1000,err,[1 5 20 100 500]);
clabel(Cont,hand), axis equal, axis square
xlabel('x in km'); ylabel('y in km');
text(.35*Rmax/1000,1.05*Rmax/1000,'Thickness error (in m).');
disp(['max err (at tf)   = ' num2str(max(max(err))) ' meters']);
disp(['center err (tf)   = ' num2str(H(1,1)-HHexactf(1,1)) ' meters']);
disp(['work (= N^2 M)    = ' int2str(4*Nx*Nx*Mt)]);
Vf=dx*dy*sum(sum(volc.*H))/4;
disp(['final num volume  = ' num2str(Vf/1e9,30) ' cubic km'])
disp(['volume differenc  = ' num2str(Vf-V0,30) ' cubic m'])

% show margin "mislocation"
figure(4), spy(flipud(xor(H>0,HHexactf>0))), xlabel('')
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

function RR=stabindex(ix,iy)
% compute certain diffusivities for max estimate
global H dx dy fx fy n2 nm Rx Ry Gam
px=ix+1; py=iy+1; my=iy-1;
Hin=H(ix,iy); Hbr=(H(px,iy) + Hin)/2;
dHsr=((H(px,iy)-Hin)/dx).^2 + ((H(px,py)+H(ix,py)-H(px,my)-H(ix,my))/fy).^2;
Df=Gam*max(max( Hbr.^n2.*dHsr.^nm ));
RR=Rx*Df;

