function [r,W,P,h,vb] = radialsteady(dofigs)
% RADIALSTEADY Compute exact solution documented in dampnotes.pdf.

if nargin<1, dofigs=true; end

p = params();

% constants specific to exact solution
h0   = 500.0;            % m     center thickness
v0   = 100.0 / p.spera;  % m/s   sliding velocity at margin
R0   = 25.0e3;           % m     ideal ice cap radius
R1   = 5.0e3;            % m     onset of sliding
Phi0 = 0.2 / p.spera;    % m/s   water input rate is 20 cm/a

vphi0 = Phi0 / (2 * p.c0);
L = 0.9 * R0;            % m     actual margin location

% WcL is key constant in construction; is initial condition on W(r)
hL  = h0 * (1 - (L/R0).^2);
PoL = p.rhoi * p.g * hL;
sbL = ( (p.c1 * v0 / (p.c2 * p.A)) )^(1/3);
WcL = (sbL^3 * p.Wr - PoL^3 * p.Y0) / (sbL^3 + PoL^3);
fprintf('  W(L) should satisfy  %.6f = W_c(L) <= W(L) <= %.6f = W_r\n',WcL,p.Wr)

if dofigs
  set(0,'defaultaxesfontsize',14)
  set(0,'defaultlinelinewidth',3.0)

  % grid for showing O, N, U regions behind exact W(r) solution
  dr = L / 500;
  r = 0:dr:L;
  dW = p.Wr/500;
  W = dW:dW:1.2*p.Wr;
  [rr WW] = meshgrid(r,W);
  Po = p.rhoi * p.g * h0 * (1 - (rr/R0).^2);
  vb = v0 * (rr - R1).^5 / (L - R1)^5;
  vb(rr < R1) = 0.0;
  ratio = psteady(p,Po,vb,WW) ./ Po;
  figure(97), clf
  [cc ll] = contour(rr/1000.0,WW,ratio,[0.000001 0.999999],'k','linewidth',1.5);
  text(2,0.5,'O','fontsize',16.0)
  text(12,1.1,'O','fontsize',16.0)
  text(12,0.5,'N','fontsize',16.0)
  text(22,0.3,'U','fontsize',16.0)
  xlabel('r  (km)'), ylabel('W  (m)')
  axis([0 R0/1000.0 0 1.2*p.Wr]), grid on

  clear dr r rr h Po vb W WW
end

%return
% in octave this requires "odepkg", but then fails because of negative step direction

% solve the ODE
wopt = odeset('RelTol', 1e-12,'AbsTol', 1e-9);
% others: 'InitialStep', 'MaxStep', 'NormControl', 'OutputFcn'
%[r,W] = ode45(@WODE,[L 0.0],WcL,wopt);
[r,W] = ode15s(@WODE,[L 0.0],WcL,wopt);  % stiff solver gives same result but more efficiently
fprintf('  numerical ODE solution generated %d r values in 0 <= r <= L = %.3f km\n',...
        length(r),max(r)/1e3)

if dofigs
  % show ODE soln W(r) onto existing fig
  figure(97), hold on
  plot(r/1000.0,W,'k--');
  plot(r(1)/1000,W(1),'k.','markersize',25)
  hold off
end

% check bounds  W_c(r) <= W(r) <= W_r
vb = v0 * (r - R1).^5 / (L - R1)^5;
vb(r < R1) = 0.0;
sb = ((p.c1 * vb) ./ (p.c2 * p.A)).^(1/3);
h = h0 * (1 - (r/R0).^2);
Po = p.rhoi * p.g * h;
Wc = (sb.^3 * p.Wr - Po.^3 * p.Y0) ./ (sb.^3 + Po.^3);
if any(W < 0)
  error('points exist where solution W(r) < 0'), end
if any(W < Wc)
  error('points exist where solution W(r) < W_c(r)'), end
if any(W > p.Wr)
  error('points exist where solution W(r) > W_r'), end

% compute pressure solution
P = psteady(p,Po,vb,W);

if dofigs
  % show pressure solution
  figure(98), clf
  plot(r/1000.0,Po/1e5,'k',r/1000.0,P/1e5,'k--');
  hold on
  plot([L L]/1000.0,[Po(1)/1e5 0],'k')
  xlabel('r  (km)'), ylabel('pressure  (bar)')
end

  function dW = WODE(r,W)
  % this function defines the right side of the ODE for W:  W'(r) = WODE(r,W)
  if r < R1
    sb  = 0.0;
    dsb = 0.0;
  else
    CC  = p.c1 / (p.c2 * p.A);
    % vb = v0 * (r - R1).^5 / (L - R1)^5   and   sb = (CC * vb)^(1/3)
    CZ  = ( CC * v0 / (L-R1)^5 )^(1/3);
    sb  = CZ * (r - R1).^(5/3);
    dsb = (5/3) * CZ * (r - R1).^(2/3);
  end
  dPo   = - (2 * p.rhoi * p.g * h0 / R0^2) * r;
  numer = dsb .* (W + p.Y0) .* (p.Wr - W);
  tmp1  = (W + p.Y0).^(4/3) .* (p.Wr - W).^(2/3);
  numer = numer - ( vphi0 * r ./ W + dPo) .* tmp1;
  denom = (1/3) * (p.Wr + p.Y0) * sb + p.rhow * p.g * tmp1;
  dW    = numer ./ denom;
  end

end

