% convergence plot with resolution
%

ne = [ 11 15 21  31 41 51 61 81 101];
deg = 360./(ne*3*4);
l1 = [ 0.167214E-01 0.117563E-01  0.659899E-02  0.167060E-02 0.304366E-03   ...
       0.7031610E-04 0.3060075E-04  0.8771899E-05  0.3054130E-05 ];
l2 = [   0.225783E-01    0.161304E-01  0.101460E-01 0.329885E-02  0.659757E-03 ...
         0.1286152E-03  0.5274167E-04  0.1644966E-04  0.5205277E-05  ];
l8 = [ 0.116495E+00  0.860934E-01  0.663686E-01 0.256036E-01 0.702307E-02 ...
       0.1664922E-02  0.5380025E-03  0.1677521E-03 0.8359257E-04];

% data with DT ~ NE^-2:
% $$$ ne = [  15   31 41 51 61 81];
% $$$ deg = 360./(ne*3*4);
% $$$ l1 = [0.1208423E-01   0.1692709E-02   0.3405964E-03   0.1028254E-03   0.5955022E-04   0.2812692E-04 ];
% $$$ l2 = [ 0.1650990E-01  0.3305126E-02  0.6896965E-03  0.1870320E-03 ...
% $$$        0.1129990E-03  0.5829328E-04 ];
% $$$ l8=[  0.9270795E-01     0.2544988E-01  0.7428411E-02  0.1870845E-02 ...
% $$$       0.6359561E-03  0.2822944E-03 ];
% $$$ 





last=length(ne);
last1=last-2;

disp(sprintf('l1 slope = %f',log ( l1(last1)/l1(last) )/log( ne(last)/ne(last1) ) ));
disp(sprintf('l2 slope = %f',log ( l2(last1)/l2(last) )/log( ne(last)/ne(last1) ) ));
disp(sprintf('l8 slope = %f',log ( l8(last1)/l8(last) )/log( ne(last)/ne(last1) ) ));


% bv filter results
bvl1 = [0.4844710E+00  0.4242726E-01  0.9967520E-02 0.3247885E-02 ];
bvl2 = [0.2567838E+00  0.2668581E-01  0.7614058E-02  0.2595111E-02 ];
bvl8 = [ 0.1692518E+00 0.2650571E-01 0.9347796E-02 0.3099508E-02];


thick=3.0;
figure(1); clf
set(gca,'FontWeight','bold')
set(gca,'FontSize',18)
set(gca,'LineWidth',2.0) 

loglog(ne,l2,'bo-','LineWidth',thick); hold on;
loglog(ne,l1,'go-','LineWidth',thick); 
loglog(ne,l8,'ro-','LineWidth',thick); 

%loglog(ne,bvl2,'co-','LineWidth',1); 
%loglog(ne,bvl1,'co-','LineWidth',1); 
%loglog(ne,bvl8,'co-','LineWidth',1); 

xref = [ne(1),ne(length(ne))];


yref2 = l2(6) * (xref/ne(6)).^(-4);
loglog(xref,yref2,'k-','LineWidth',2)

%yref2 = l2(3) * (xref/ne(3)).^(-1.7);
%loglog(xref,yref2,'k-','LineWidth',2)


axis([10.  110.    1e-6    .1]);

legend('L2','L1','L\infty','location','northeast')
title('CASE 1 DAY=12')
ylabel('ERROR')
xlabel('NE')

print -dpng converge.png




