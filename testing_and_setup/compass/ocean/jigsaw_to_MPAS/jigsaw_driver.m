function jigsaw_driver
%-----------------------------------------------------------
%   Mark Petersen (mpetersen@lanl.gov)
%   Phillip Wolfram (pwolfram@lanl.gov)
%   07/01/2018
%-----------------------------------------------------------

   addpath('jigsaw-geo-matlab')
   
   %------------------------------------ Load cellWidth, lon, lat

   load('cellWidthVsLatLon.mat')

   %------------------------------------ setup files for JIGSAW
   
   opts.geom_file = [ 'mesh.msh'];      % GEOM file
   opts.jcfg_file = [ 'mesh.jig'];      % JCFG file
   opts.mesh_file = [ 'mesh-MESH.msh']; % MESH file
   opts.hfun_file = [ 'mesh-HFUN.msh']; % HFUN file
           
   %------------------------------------ save HFUN data to file
   
   hmat.mshID = 'ELLIPSOID-GRID';
   hmat.point.coord{1} = lon*pi/180 ;
   hmat.point.coord{2} = lat*pi/180 ;
   hmat.value = cellWidth ;
   
   savemsh(opts.hfun_file,hmat) ;
   
   %------------------------------------ define JIGSAW geometry
   
   geom.mshID = 'ELLIPSOID-MESH';
   geom.radii = 6371.*ones(3,1) ;
   
   savemsh(opts.geom_file,geom) ;
   
   %------------------------------------ build mesh via JIGSAW!
   
   opts.hfun_scal = 'absolute';
   opts.hfun_hmax = +inf ;
   opts.hfun_hmin = +0.0 ;
   opts.mesh_dims = +2 ;               % 2-dim. simplexes
   opts.optm_qlim = 0.9375 ;
   opts.verbosity = +1 ;
   
   mesh = jigsaw  (opts) ;
