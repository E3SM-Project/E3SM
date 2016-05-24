
%% PLOT MODEL STREAM/PLASTIC BED RESULTS AND COMPARE W/ RAYMOND/SCHOOF
%% SOLUTIONS

clear all

rho = 910;
g = -9.81;

flag = 0;       %% USE RAYMOND PROFILE
flag = 1;       %% USE SCHOOF PROFILE 

kinflag = 1;    %% apply kinematic bc (analytic soln) at up/downstream ends
% kinflag = 0;    %% apply 0 vel bc at up/downstream ends

n = 3;
m = 1.55;

% r = 20;         % !! r needs to be even # divisible by 2 !!
% c = 20;
% % c = 200;

r = 40;
c = 40;
% c = 400;

levels = 5;
% levels  = 11;

L = 1.4e4;
A = 1e-16;
B = A^(-1/n);

H = 1000;
dsdx = -1.0e-3;

% ice stream width
w = 50e3;
W = w / 2;

yy = linspace( 0, W, (r-4)/2 );

dy = mean( diff( yy ) );
dx = dy;

taud = rho * g * H * dsdx;

% Schoof yield stress distribution
tau0s = taud * abs( yy / L ).^m;

% yield stress (for Raymond solution only)
% tau0r = mean( tau0s );
tau0r = 5.2e3;

% Raymond solution
ur = 2 * A / (n+1) * ( (taud - tau0r)/H )^n * ( W^(n+1) - yy.^(n+1) );

us = -2*taud^3*L^4/(B^3*H^3) * ( ((yy/L).^4 - (m+1)^(4/m))/4 - 3*( abs(yy/L).^(m+4) ...
    - (m+1)^(1+4/m) )/((m+1)*(m+4)) + 3*( abs(yy/L).^(2*m+4) - (m-1)^(2+4/m) )/((m+1)^2*(2*m+4)) ...
    - ( abs(yy/L).^(3*m+4) - (m+1)^(3+4/m) )/ ( (m+1)^3*(3*m+4)) );

sscale = 2*taud^3*L^4/(B^3*H^3);

ind = find( abs( yy ) >= W ); us(ind) = min( min( us ) );

us = us - min( us );

if( flag == 0 )
    figure(198), clf
    subplot(2,1,1), hold on
    plot( yy/1e3, ur - min(ur), 'r-', 'linewidth', 2.0 ), hold on
    xlabel( 'dist across flow (m)'), ylabel( 'velocity (m/a)'), title( 'Raymond solution' )
    box on
    subplot(2,1,2), hold on
    xlabel( 'dist across flow (m)'), ylabel( 'yield stress (kPa)')
    box on
else
    figure(199), clf
    subplot(2,1,1), hold on
    plot( yy/1e3, us - min(us), 'r-', 'linewidth', 2.0 ), hold on
    xlabel( 'dist across flow (m)'), ylabel( 'velocity (m/a)'), title( 'Schoof solution' )
    box on
    subplot(2,1,2), hold on
    xlabel( 'dist across flow (m)'), ylabel( 'yield stress (kPa)')
    box on
end


%% put data into .mat file, to be read by .py script w/ config file, to make input .nc files
thck = H * ones( r, c );
topg = repmat( 5e3 + dsdx*dx*[0:c-1], r, 1 );
usrf = topg + thck;


if( flag == 0 )         % assign Raymond profile
    tauf_profile = tau0r*ones(r-7,1);        
    tauf = 0.7e5 * ones( r-1, c-1 );
    if( kinflag ~= 1 )
        tauf(4:end-3,3:end-2) = repmat( tauf_profile, 1, c-5 );     %% no slip at/near up/downstream ends
    else
        tauf(4:end-3,:) = repmat( tauf_profile, 1, c-1 );     %% use periodic bcs for cont. along flow
    end
    u_profile = repmat( [ 0 0 fliplr(ur) ur(2:end) 0 0 ], levels, 1 );     
else                    % assign Schoof
    
    tauf_profile = [ fliplr(tau0s(2:end)), tau0s ]';        
    tauf = 0.7e5 * ones( r-1, c-1 );
    if( kinflag ~= 1 )
        tauf(3:end-2,3:end-2) = repmat( tauf_profile, 1, c-5 );     %% no slip at/near up/downstream ends
    else
        tauf(3:end-2,:) = repmat( tauf_profile, 1, c-1 );     %% use periodic bcs for cont. along flow
    end
    u_profile = repmat( [ 0 0 fliplr(us) us(2:end) 0 0 ], levels, 1 );
end


%% use shear strain rate at up/downstream ends to calculate across-flow velocity profile
%% necessary to balance the shear
x = [0:dx:dx*length(u_profile)-1];
dudy = gradient( u_profile(1,:), dy );
v_profile = -dudy*dy;

figure(999), clf
subplot( 2,1,1 ),hold on
plot( x, dudy, 'bo-' ), box on
xlabel( 'dist (m)' ), ylabel( 'shear strain rate (1/a)' )
subplot( 2,1,2 ),hold on
plot( x, v_profile, 'bo-' ), box on
xlabel( 'dist (m)' ), ylabel( 'across-flow component of vel (m/a)' )

v_profile = repmat( v_profile, levels, 1 );


% figure(200),clf
% subplot(2,2,1),hold on
% imagesc( thck ), axis xy, axis equal, axis tight, colorbar, title( 'thickness (m)' )
% subplot(2,2,2),hold on
% imagesc( usrf ), axis xy, axis equal, axis tight, colorbar, title( 'surface (m)' )
% subplot(2,2,3),hold on
% imagesc( topg ), axis xy, axis equal, axis tight, colorbar, title( 'bed (m)' )
% subplot(2,2,4),hold on
% imagesc( tauf/1e3 ), axis xy, axis equal, axis tight, colorbar, title( 'yield stress (kPa)' )


%% optional: add analytic solution at up/downstream ends as kin vel bc
kinbcmask = zeros(size(tauf));
uvelhom = zeros( levels, r-1, c-1 );
vvelhom = zeros( levels, r-1, c-1 );

% %% HACKS
% tauf(3:end-2,end-1) = 2.19e4; tauf(3:end-2,2) = 2.19e4;
% ind = find( tauf == max( max( tauf ) ) ); tauf(ind) = 4e4;
% kinbcmask(:,1:2) = 1; kinbcmask(:,end-1:end) = 1;
% kinbcmask(1,:) = 1; kinbcmask(end,:) = 1;

if( kinflag == 1)
    uvelhom( :, :, end ) = u_profile; uvelhom( :, :, 1 ) = u_profile; 
    vvelhom( :, :, end ) = v_profile; vvelhom( :, :, 1 ) = -v_profile; 
    kinbcmask(:,1) = 1; kinbcmask(:,end) = 1;
    ind = find( tauf <= 0 ); tauf(ind) = 1e-10;
end

%% for newer code, uvelhom = uvel, etc.
uvel = uvelhom; vvel = vvelhom;

beta = tauf;    %% for now, beta is being used as a proxy for tauf

save stream.mat usrf topg thck beta uvel vvel kinbcmask levels 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot output

%% open file
if( flag == 0 )
    filename = 'stream.raymond.nc'; 
else
    filename = 'stream.schoof.nc'; 
end

filename = 'stream.out.nc'; 

    
ncid = netcdf.open( filename, 'nowrite' );

%% get id for variable names (use ncview to find var names?)
id_uvel = netcdf.inqvarid( ncid, 'uvel' );
id_vvel = netcdf.inqvarid( ncid, 'vvel' );
vel_scale = netcdf.getAtt(ncid,id_uvel,'scale_factor' );

id_btractx = netcdf.inqvarid( ncid, 'btractx' );
id_btracty = netcdf.inqvarid( ncid, 'btracty' );
btract_scale = netcdf.getAtt(ncid,id_btractx,'scale_factor' );

uvel = permute( netcdf.getvar(ncid, id_uvel ), [ 2 1 3 ] ) * vel_scale;
vvel = permute( netcdf.getvar(ncid, id_vvel ), [ 2 1 3 ] ) * vel_scale;
btractx = netcdf.getvar(ncid, id_btractx )' * btract_scale;
btracty = netcdf.getvar(ncid, id_btracty )' * btract_scale;

yy2 = [ yy yy(end)+dx yy(end)+2*dx ];
yy2 = [ -fliplr(yy2(2:end)), yy2 ];

if( flag == 0 )
    figure(198)
    subplot(2,1,1), hold on
    plot( yy2/1e3, uvel(:,round(c/2),1), 'bo:' )
    plot( yy2/1e3, uvel(:,end,1), 'b*' )              %% boundary value
    legend( 'analytic', 'model', 'boundary' )
    subplot(2,1,2), hold on
    plot( yy2/1e3, tauf(:,round(c/2))/1e3, 'r-', 'linewidth', 2.0 )
    plot( yy2/1e3, -btractx(:,round(c/2),1)/1e3, 'bo:' )
    legend( 'specified', 'model' )
else
    figure(199)
    subplot(2,1,1), hold on
    plot( yy2/1e3, uvel(:,round(c/2),1), 'bo:' )
    plot( yy2/1e3, uvel(:,end,1), 'b*' )              %% boundary value
    legend( 'analytic', 'model', 'boundary' )
    subplot(2,1,2), hold on
    plot( yy2/1e3, tauf(:,round(c/2))/1e3, 'r-', 'linewidth', 2.0 )
    plot( yy2/1e3, -btractx(:,round(c/2),1)/1e3, 'bo:' )
    legend( 'specified', 'model' )
end



