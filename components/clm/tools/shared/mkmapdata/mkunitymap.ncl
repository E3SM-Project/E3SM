;
;  mkunitymap.ncl
;
;  Create a unity map file either between two identical grids or between two
;  grids that do NOT intersect at all.
;
;  Erik Kluzek
;  Dec/07/2011
;  $Id: getco2_historical.ncl 23741 2010-06-12 19:27:09Z erik $
;  $HeadURL;
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
  ; Set a few constants needed later
  cdate  = systemfunc( "date +%y%m%d" );
  ldate  = systemfunc( "date" );
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE or use ENV VARIABLE SETTINGS
  ; Edit the following as needed to interpolate to a new resolution.
  gridfile1 = getenv("GRIDFILE1");    ; Get name of the first  SCRIP grid file
  gridfile2 = getenv("GRIDFILE2");    ; Get name of the second SCRIP grid file

  outfilename = getenv("MAPFILE");    ; Get name of the output mapping file

  print_str  = getenv("PRINT");       ; Do Extra printing for debugging

  if ( ismissing(gridfile1) )then
     print( "ERROR: GRIDFILE1 is missing!" );
     exit
  end if
  if ( ismissing(gridfile2) )then
     print( "ERROR: GRIDFILE2 is missing!" );
     exit
  end if
  if ( ismissing(outfilename) )then
     print( "ERROR: MAPFILE is missing!" );
     exit
  end if
  if ( ismissing(print_str) )then
     printn = False;
  else
     if ( print_str .eq. "TRUE" )then
        printn = True;
     else
        printn = False;
     end if
  end if

  ;
  ; Open up the input grid files
  ;
  nca = addfile( gridfile1, "r" );
  ncb = addfile( gridfile2, "r" );

  system( "/bin/rm -f "+outfilename );
  if ( printn )then
     print( "output mapping file to create: "+outfilename );
  end if
  nc = addfile( outfilename, "c" );
  ;
  ; Define dimensions
  ;  
  n_a = dimsizes( nca->grid_center_lat );
  n_b = dimsizes( ncb->grid_center_lat );
  if ( n_a .ne. n_b )then
     print( "ERROR: dimensions of input SCRIP grid files is NOT the same!" );
     exit
  end if
  if ( any(ncb->grid_imask .ne. 1.0d00) )then
     print( "ERROR: the mask of the second file isn't identically 1!" );
     print( "(second file should be land grid file)");
     exit
  end if
  chkvars = (/ "grid_center_lat", "grid_center_lon", "grid_corner_lat", "grid_corner_lon" /);
  do i = 1, dimsizes(chkvars)-1
     if ( any(nca->$chkvars(i)$ .ne. ncb->$chkvars(i)$) )then
        print( "ERROR: the grid variables are different between the two files!: "+chkvars(i)  );
        exit
     end if
  end do
  n_s = n_a;
  dimnames = (/ "n_a", "n_b", "n_s", "nv_a", "nv_b", "num_wgts", "src_grid_rank", "dst_grid_rank" /);
  dsizes   = (/ n_a,     n_b,   n_a,      4,      4,          1,               2,        2/);
  is_unlim = (/ False, False, False,  False,  False,      False,           False,   False /);
  filedimdef( nc, dimnames, dsizes, is_unlim );

  ;
  ; Define grid dimensions
  ;
  filevardef( nc, "src_grid_dims", "integer", (/ "src_grid_rank" /))
  nc->src_grid_dims = (/nca->grid_dims/)
  filevardef( nc, "dst_grid_dims", "integer", (/ "dst_grid_rank" /))
  nc->dst_grid_dims = (/ncb->grid_dims/)

  ;
  ; Define variables
  ;
  cvars = (/ "yc",              "xc",              "yv",              "xv",              "mask"  /);
  gvars = (/ "grid_center_lat", "grid_center_lon", "grid_corner_lat", "grid_corner_lon", "grid_imask" /);

  do i = 0, dimsizes(cvars)-1
     var = cvars(i)+"_a";
     if ( cvars(i) .eq. "yv" .or. cvars(i) .eq. "xv" )then
        dnamesa = (/ "n_a", "nv_a" /);
        dnamesb = (/ "n_b", "nv_b" /);
     else
        dnamesa = (/ "n_a" /);
        dnamesb = (/ "n_b" /);
     end if
     filevardef ( nc, var, typeof(nca->$gvars(i)$), dnamesa );
     filevarattdef ( nc, var, nca->$gvars(i)$ );
     nc->$var$ = (/ nca->$gvars(i)$ /);
     var = cvars(i)+"_b";
     filevardef ( nc, var, typeof(nca->$gvars(i)$), dnamesb );
     filevarattdef ( nc, var, ncb->$gvars(i)$ );
     nc->$var$ = (/ ncb->$gvars(i)$ /);
     delete( dnamesa );
     delete( dnamesb );
  end do
  filevardef ( nc, "area_a", "double", (/ "n_a" /) );
  filevardef ( nc, "area_b", "double", (/ "n_b" /) );
  filevardef ( nc, "frac_a", "double", (/ "n_a" /) );
  filevardef ( nc, "frac_b", "double", (/ "n_b" /) );
  ;
  ; Attributes
  ;
  nc->area_a@units = "square radians";
  nc->frac_a@units = "unitless";
  nc->area_b@units = nc->area_a@units;
  nc->frac_b@units = nc->frac_a@units;
  nc@conventions   = "NCAR-CESM";
  nc@domain_a      = gridfile1;
  nc@domain_b      = gridfile2;
  nc@grid_file_src = gridfile1;
  nc@grid_file_dst = gridfile2;
  nc@title   = "SCRIP mapping file between identical grids without ocean";
  nc@history = ldate+": create using mkunitymap.ncl";
  nc@Version = "$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/branches/newgrid/models/lnd/clm/tools/ncl_scripts/ndeplintInterp.ncl $";
  nc@Revision = "$Id: ndeplintInterp.ncl 25175 2010-10-16 03:47:50Z erik $";

  ;
  ; Fraction
  ;
  nc->frac_a = int2dble( (/nc->mask_a/) );
  nc->frac_b = int2dble( (/nc->mask_b/) );
  ;
  ; Area
  ;
  nc->area_a = gc_qarea( nc->yv_a(:,:), nc->xv_a(:,:) );
  nc->area_b = gc_qarea( nc->yv_b(:,:), nc->xv_b(:,:) );
  ;
  ; Weights
  ;
  filevardef ( nc, "col", "integer", (/ "n_s" /) );
  filevardef ( nc, "row", "integer", (/ "n_s" /) );
  filevardef ( nc, "S",   "double",  (/ "n_s" /) );

  nc->col = ispan( 1, n_s, 1 );
  nc->row = nc->col;
  nc->S   = 1.0d00;

end
