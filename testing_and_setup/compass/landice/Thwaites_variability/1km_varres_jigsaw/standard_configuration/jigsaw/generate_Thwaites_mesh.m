function Thwaites_test()

% note need path to jigsaw-geo-matlab and jigsaw-geo-matlab/mesh-util in MATLAB path.

%------------------------------------ output from this script, input to JIGSAW    

    
    opts.geom_file = ...                % GEOM file
        ['./bounding_box.msh'];    

    opts.jcfg_file = ...                % JCFG file
        ['./config.jig'];

    opts.hfun_file = ...                % HFUN file
        ['./density2.msh'];    
    
%------------------------------------ output file from JIGSAW    

    opts.mesh_file = ...                % MESH file
        ['./thwaites_jigsaw_mesh.msh'];
    

 
%------------------------------------ define JIGSAW geometry
    
    geom.mshID = 'EUCLIDEAN-MESH';

    geom.point.coord = [    % list of xy "node" coordinates
        -1.63e6, -0.82e6, 0
        -1.63e6, -0.14e6, 0
        -1.05e6, -0.14e6, 0
        -1.05e6, -0.82e6, 0 ] ;
    
    geom.edge2.index = [    % list of "edges" between nodes
        1, 2, 0
        2, 3, 0
        3, 4, 0
        4, 1, 0  ] ;    
        
    savemsh(opts.geom_file, geom) ;
    
%------------------------------------ compute HFUN over GEOM    
       
% densFile='~/documents/mpas-git/mesh_generation/matlab_triangle_varres_meshgen/density.nc';
% xpos=ncread(densFile, 'x');
% ypos=ncread(densFile, 'y');
% dens = ncread(densFile, 'density');
% dens(1:10,:) = 0.2;  dens(end-10:end,:)=0.2; dens(:,1:10) = 0.2;  dens(:,end-10:end)=0.2; 
% dens=(dens');
% hfun = (1.0./dens).^1.0;
% hfun = hfun*10000.0;
    

% old method using data density function
densFile='~/documents/mpas-git/mesh_generation/matlab_triangle_varres_meshgen/bikedens.nc.for-jigsaw5.nc';
xpos=ncread(densFile, 'x');
ypos=ncread(densFile, 'y');
dens = ncread(densFile, 'spacing')';
hfun = dens*1000.0;




   [XPOS,YPOS] = meshgrid(xpos,ypos) ;
    

    hmat.mshID = 'EUCLIDEAN-GRID' ;
    hmat.point.coord{1} = xpos ;
    hmat.point.coord{2} = ypos ;
    hmat.value = hfun ;

    
%% use analytical density fn
pz1=[-1586136.27728525, -360943.403207161];
pz2=[-1525083.64347291, -536525.430375485];
m=(pz2(2)-pz1(2))/(pz2(1)-pz1(1));
b=pz1(2)-m*pz1(1);
mask1 = (YPOS > (m*XPOS+b));


pa1=[-1405819.85287125, -322799.544196507];
pa2=[-1307779.2345637, -359676.629507899];    
m=(pa2(2)-pa1(2))/(pa2(1)-pa1(1));
b=pa1(2)-m*pa1(1);
mask2 = (YPOS < (m*XPOS+b));

pb1=[-1316817.51517972, -507739.466672411];
pb2=[-1317026.23068749, -338870.888229361];    
m=(pb2(2)-pb1(2))/(pb2(1)-pb1(1));
b=pb1(2)-m*pb1(1);
mask3 = (YPOS < (m*XPOS+b));

pc1=[-1313349.89163329, -497336.596033141];
pc2=[-1406028.56837902, -542304.802952857];    
m=(pc2(2)-pc1(2))/(pc2(1)-pc1(1));
b=pc1(2)-m*pc1(1);
mask4 = (YPOS > (m*XPOS+b));

% circle location and radius
c=[-1413295.23341991, -440814.734214159];
r=120000.0;

mask = mask1 & mask2 & mask3 & mask4;

dist=zeros(size(XPOS));
basevalue=1000.0;
den=ones(size(XPOS)) * basevalue;
for i=1:length(xpos)
    for j=1:length(ypos)
        pt=[xpos(i), ypos(j)];
        if mask1(j,i)==0  % to west of initial GL
           d = point_to_line(pt, pz1, pz2);
           f = min(2, d * 1.0/100000.0 + 1);
        elseif xpos(i) > c(1);
           d = max(0.0, sqrt((pt(1)-c(1))^2+(pt(2)-c(2))^2) - r);
           f = min(8, d * 7.0/250000.0 + 1); 
        else
           f =1;
        end    
    %    d = point_to_line(pt, pz1, pz2);
     %   d = min(d, point_to_line(pt, pb1, pb2));
      %  d = min(d, point_to_line(pt, pc1, pc2));
        %f = d * 4.0/200000.0 + 1;
        %f=min(4,f);
        %f=max(1,f);
        den(j,i)=basevalue*f;

    end
end

%den(mask1 & mask2 & mask3 & mask4) = 10000.0;

hmat.value = den;
%%
    
    
    savemsh(opts.hfun_file,hmat) ;
    
%------------------------------------ build mesh via JIGSAW! 
  
    opts.hfun_scal = 'absolute';
    opts.hfun_hmax = +inf ;             % null HFUN limits
    opts.hfun_hmin = 0.00 ;
  
    opts.mesh_dims = +2 ;               % 2-dim. simplexes
    
    opts.optm_qlim = 0.9375 ;
   
    opts.mesh_top1 = true ;             % for sharp feat's
    opts.geom_feat = true ;
    
    mesh = jigsaw  (opts) ;
 
%------------------------------------ draw mesh/cost outputs

    ang2 = triang2( ...                 % calc. tri-angles
        mesh.point.coord(:,1:2), ...
        mesh.tria3.index(:,1:3)) ;
            
    t_90 = max(ang2,[],2) > 90.0 ;
    t_95 = max(ang2,[],2) > 95.0 ;
    
    figure(1); clf;
    patch ('faces',geom.edge2.index(:,1:2), ...
        'vertices',geom.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    hold on; axis image;
    title('JIGSAW GEOM data') ;

    figure(2); clf; hold all;
    pcolor(XPOS,YPOS,hmat.value);
    axis equal; 
    shading interp ;
    title('JIGSAW HFUN data') ;
        patch ('faces',geom.edge2.index(:,1:2), ...
        'vertices',geom.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[1,.1,.1], ...
        'linewidth',1.5) ;
    colorbar();


    figure(3); clf;
    patch ('faces',mesh.tria3.index(:,1:3), ...
        'vertices',mesh.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image;
    patch ('faces',mesh.tria3.index(t_90,1:3), ...
        'vertices',mesh.point.coord(:,1:2), ...
        'facecolor','y', ...
        'edgecolor',[.2,.2,.2]) ;
    patch ('faces',mesh.tria3.index(t_95,1:3), ...
        'vertices',mesh.point.coord(:,1:2), ...
        'facecolor','r', ...
        'edgecolor',[.2,.2,.2]) ;
    patch ('faces',mesh.edge2.index(:,1:2), ...
        'vertices',mesh.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    axis equal
    title('JIGSAW TRIA mesh') ;

    drawscr(mesh.point.coord (:,1:2), ...
            mesh.edge2.index (:,1:2), ...
            mesh.tria3.index (:,1:3)) ;
    
    drawnow ;        
%     set(figure(1),'units','normalized', ...
%         'position',[.05,.55,.30,.35]) ;
%     set(figure(2),'units','normalized', ...
%         'position',[.35,.55,.30,.35]) ;
%     set(figure(3),'units','normalized', ...
%         'position',[.35,.10,.30,.35]) ;
%     set(figure(4),'units','normalized', ...
%         'position',[.05,.10,.30,.35]) ;
    drawnow ;

end
function d = point_to_line(pt, v1, v2)

%       a = v1 - v2;
%       b = pt - v2;
%       d = norm(cross([a(1), a(2), 0], [b(1), b(2), 0])) / norm(a);
       
m=(v2(2)-v1(2))/(v2(1)-v1(1));
b=v1(2)-m*v1(1);
d=abs(m*pt(1) + -1*pt(2) + b)/sqrt(m^2+(-1)^2);

%       
% a=v1; b=v2; x=pt;
% d_ab = norm(a-b);
% d_ax = norm(a-x);
% d_bx = norm(b-x);
% 
% if dot(a-b,x-b)*dot(b-a,x-a)>=0
%     A = [a,1;b,1;x,1];
%     dist = abs(det(A))/d_ab;        
% else
%     dist = min(d_ax, d_bx);
% end      
%     d=dist;  
    
end