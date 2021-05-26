function varres_dome()

% note need path to jigsaw-geo-matlab and jigsaw-geo-matlab/mesh-util in MATLAB path.

%------------------------------------ output from this script, input to JIGSAW    

    
    opts.geom_file = ...                % GEOM file
        ['./bounding_box.msh'];    

    opts.jcfg_file = ...                % JCFG file
        ['./config.jig'];

    opts.hfun_file = ...                % HFUN file
        ['./density.msh'];    
    
%------------------------------------ output file from JIGSAW    

    opts.mesh_file = ...                % MESH file
        ['./varres_dome_jigsaw_mesh.msh'];
    

 
%------------------------------------ define JIGSAW geometry
    
    geom.mshID = 'EUCLIDEAN-MESH';

    geom.point.coord = [    % list of xy "node" coordinates
        -30.0e3, -30.0e3, 0
        -30.0e3, 30.0e3, 0
         30.0e3, 30.0e3, 0
         30.0e3, -30.0e3, 0 ] ;
    
    geom.edge2.index = [    % list of "edges" between nodes
        1, 2, 0
        2, 3, 0
        3, 4, 0
        4, 1, 0  ] ;    
        
    savemsh(opts.geom_file, geom) ;
    
%------------------------------------ compute HFUN over GEOM    
       

    xpos = [-30.0e3:1000.0:30.0e3];
    ypos = xpos;
    [XPOS,YPOS] = meshgrid(xpos,ypos) ;
    
    r = (XPOS.^2+YPOS.^2).^0.5/1000.0; % radius in km

    % min/max spacing of resulting mesh
    hmax=5.0; %km
    hmin=1.0; %km
    % location of min/max spacing in radial distance from center of dome
    rmax= 3.0;  %km
    rmin = 20.0; %Km
    
    slope = ((hmin-hmax)/(rmin-rmax));
    dens = r * slope + (hmax - slope*rmax);
    dens = max(dens, hmin);
    dens = min(dens, hmax);
    hfun = dens *1000.0; % convert km to m
    
    hmat.mshID = 'EUCLIDEAN-GRID' ;
    hmat.point.coord{1} = xpos ;
    hmat.point.coord{2} = ypos ;
    hmat.value = hfun ;

    

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