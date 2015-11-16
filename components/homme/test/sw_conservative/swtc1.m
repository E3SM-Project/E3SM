% convergence plot with resolution
%

ne = [ 11 15 21  31 41 51 ];
deg = 360./(ne*3*4);
l1 = [   0.1873639E+00   0.7082156E-01   0.3180798E-01   0.1147578E-01 0.5490936E-02  0.3438183E-02];
l2 = [  0.2018837E+00   0.7761264E-01   0.3577571E-01   0.1393754E-01  0.6564592E-02  0.4183211E-02 ];
l8 = [ 0.2976752E+00  0.1460969E+00  0.7823920E-01  0.4278967E-01    0.2505294E-01  0.1881130E-01 ];

last=length(ne);
last1=last-1;

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


yref2 = l8(3) * (xref/ne(3)).^(-2);
loglog(xref,yref2,'k-','LineWidth',2)

yref2 = l8(3) * (xref/ne(3)).^(-1.7);
loglog(xref,yref2,'k-','LineWidth',2)

%yref1 = l8(3) * (xref/ne(3)).^(-1);
%loglog(xref,yref1,'k-','LineWidth',1)

axis([5.0000  100.0000    0.0010    1.0000]);

legend('L2','L1','L\infty','location','northeast')
title('CASE 1 DAY=12')
ylabel('ERROR')
xlabel('NE')

print -dpng converge.png




