% convergence plot with resolution
%

ne_mono = [ 11 15 21 31 41];

% monotone limiter + monotone hypervis
%           l2             l2           l8
data=[  0.1585986E-01  0.3322470E-01  0.9069445E-01 ; ... 
        0.5998192E-02  0.1474055E-01  0.5322885E-01 ; ...
        0.2219639E-02  0.5974962E-02  0.2945256E-01 ; ...
        0.6683278E-03  0.2117329E-02  0.1479235E-01 ; ...
        0.2865575E-03  0.1011187E-02  0.9072316E-02 ];

l1mono = data(:,1);
l2mono = data(:,2);
l8mono = data(:,3);



disp(sprintf('l1 slope = %f',log ( l1mono(2)/l1mono(4) )/log( ne_mono(4)/ne_mono(2) ) ));
disp(sprintf('l2 slope = %f',log ( l2mono(2)/l2mono(4) )/log( ne_mono(4)/ne_mono(2) ) ));
disp(sprintf('l8 slope = %f',log ( l8mono(2)/l8mono(4) )/log( ne_mono(4)/ne_mono(2) ) ));



% zero limiter + zero hypervis
data = [ 0.5244622E-02  0.5225726E-02  0.8395018E-02 ; ...
         0.2080151E-02  0.2071161E-02  0.3319189E-02 ; ...
         0.7425880E-03  0.7469148E-03  0.1211872E-02 ; ...
         0.2194419E-03  0.2234200E-03  0.3650040E-03 ; ...
         0.8434988E-04  0.8481369E-04  0.1397168E-03 ];

l1zero = data(:,1);
l2zero = data(:,2);
l8zero = data(:,3);
ne_zero = [ 11 15 21 31 41];

  
% no limiters, no viscosity
ne_sem = [ 9 15 21 31 41];
data = [   0.1647045E-01  0.8258288E-02  0.6885314E-02 ; ...
           0.3171664E-03  0.1689483E-03  0.2844085E-03 ; ...
           0.6453766E-04  0.3516009E-04  0.6467689E-04 ; ...
           0.1243756E-04  0.6461270E-05  0.1771740E-04 ; ...
           0.3995662E-05  0.2233213E-05  0.7222591E-05 ];
l1sem = data(:,1);
l2sem = data(:,2);
l8sem = data(:,3);
  

thick=3.0;
figure(1); clf
set(gca,'FontWeight','bold')
set(gca,'FontSize',18)
set(gca,'LineWidth',2.0) 

loglog(ne_mono,l8mono,'bo-','LineWidth',thick);  hold on;
loglog(ne_zero,l8zero,'go-','LineWidth',thick); 
loglog(ne_sem,l8sem,'ro-','LineWidth',thick); 

xref = [ne_mono(1)/1.2,1.2*ne_mono(length(ne_mono))];
yref = l8mono(3) * (xref/ne(3)).^(-2);
loglog(xref,yref,'k-','LineWidth',2)
text(xref(2)*.8,yref(2)*.60,'slope=-2','FontWeight','bold','FontSize',16)

xref = [ne_zero(1)/1.2,1.2*ne_zero(length(ne_zero))];
yref = l8zero(3) * (xref/ne_zero(3)).^(-3);
loglog(xref,yref,'k-','LineWidth',2)
text(xref(2)*.9,yref(2)*.70,'slope=-3','FontWeight','bold','FontSize',16)

xref = [ne_sem(1)/1.2,1.2*ne_sem(length(ne_sem))];
yref = l8sem(3) * (xref/ne_sem(3)).^(-4);
loglog(xref,yref,'k-','LineWidth',2)
text(xref(2)*.8,yref(2)*.70,'slope=-4','FontWeight','bold','FontSize',16)


axis([4.0000  100.    1e-6 1]);

legend('monotone','positive','oscillatory','location','northeast')
title('CASE 1    DAY=12')
ylabel('MAX ERROR')
xlabel('NE')
hold off;

print -r100 -dpng converge_max.png




loglog(ne,l2mono,'bo-','LineWidth',thick); hold on;
loglog(ne,l2zero,'go-','LineWidth',thick); 
loglog(ne_sem,l2sem,'ro-','LineWidth',thick); 

xref = [ne_mono(1)/1.2,1.2*ne_mono(length(ne_mono))];
yref = l2mono(3) * (xref/ne(3)).^(-2);
loglog(xref,yref,'k-','LineWidth',2)
text(xref(2)*1.0,yref(2)*.80,'slope=-2','FontWeight','bold','FontSize',16)

xref = [ne_zero(1)/1.2,1.2*ne_zero(length(ne_zero))];
yref = l2zero(3) * (xref/ne_zero(3)).^(-3);
loglog(xref,yref,'k-','LineWidth',2)
text(xref(2)*1.0,yref(2)*.80,'slope=-3','FontWeight','bold','FontSize',16)

xref = [ne_sem(1)/1.2,1.2*ne_sem(length(ne_sem))];
yref = l2sem(3) * (xref/ne_sem(3)).^(-4);
loglog(xref,yref,'k-','LineWidth',2)
text(xref(2)*1.0,yref(2)*.80,'slope=-4','FontWeight','bold','FontSize',16)


axis([4.0000  100.    1e-7 .1]);

legend('monotone','positive','oscillatory','location','northeast')
title('CASE 1    DAY=12')
ylabel('L2 ERROR')
xlabel('NE')
hold off;
print -r100 -dpng converge_l2.png






