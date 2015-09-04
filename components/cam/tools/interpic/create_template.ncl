;
;	create_template.ncl
;
;	Purpose: Provide a easier mechanism to create the
;		template file to interpolate a dataset. The
;		template file is used by interpic to interpolate a
;		given file. Use the same file both here and in
;		interpic.
;
;	Method: Set grid sizes and variables at top of
;		script. Copy variables and needed values from
;		file that will be interpolated. Use the template.ncl
;             subroutines to do most of the work.
;
;	Usage:	Edit the settings in the first section. Make sure
;		everything is set as it should. Then...
;
;			ncl < create_template.ncl
;
;		Settings for the standard configurations of 18, 26, and 
;		30 levels are given. If you need to change the standard
;		configuration, or if you need a different number of levels
;		you will have to edit the section that sets the hybrid levels.
;
;	Author: Erik Kluzek
;
;	$Id$
;
begin
;-----------------------------------------------------------------------------------
;  Here are all the important settings. You must edit these by hand
;  The following global variables must be set.
;
;  nlat                 Number of latitudes
;  nlon                 Number of longitudes
;  nlev                 Number of levels
;  grid                 Grid type (staggered, reduced or gaussian)
;  ntrm                 Number of wave truncation for M (Gaussian and reduced grid only)
;  ntrk                 Number of wave truncation for K (Gaussian and reduced grid only)
;  ntrn                 Number of wave truncation for N (Gaussian and reduced grid only)
;  interpfilename       Filename of file interpolating from
;  templatefilename     Output filename of template creating
;  caseid               case id (max 16 characters)
;  nlons                Array of number of longitudes for each latitude (reduced grid only)
;  var_type             Output type of fields creating
;  nstandard            Array of the number of levels for various standard level configurations
;  hyai_standard        Hybrid "A" interface levels for the standard levels.
;  hybi_standard        Hybrid "B" interface levels for the standard levels.
;  FillValue            Value to give the missing value
;  nchar                Number of characters to use for character data.
;  dimnames             Names of the dimensions for the file
;  dsizes               Dimension sizes.
;  is_unlim             Logical array to indicate if a dimension is unlimited or not.
;  rlon                 reduced grid longitudes (reduced grid only)
;  lat                  Latitudes
;  lon                  Longitudes
;  slat                 Staggered latitudes (staggered only)
;  slon                 Staggered longitudes (staggered only)
;  gw                   Gaussian weights
;  w_stag               Staggered weights (staggered only)
;
;-----------------------------------------------------------------------------------
  nlat = 19;   ; Number of latitudes
  nlon = 24;   ; Number of longitudes
  nlev = 26;   ; Number of vertical levels
  ntrm =  0;   ; Truncation in "M"
  ntrk =  0;   ; Truncation in "K"
  ntrn =  0;   ; Truncation in "N"
  ; Grid type ("gaussian", "staggered" or "reduced")
  grid = "staggered";
  ; Filename to interpolate from
  interpfilename = "/fs/cgd/csm/inputdata/atm/cam1/inic/gaus/SEP1.T42L26.fland.c020402.nc";
  caseid = "interp_newsstclim03Br";   ; Case id to use on output template file
  ; Name of output filename
  if ( grid .eq. "staggered" )then
    templatefilename = "SEP1."+nlat+"x"+nlon+"L"+nlev+"."+grid+".template.nc";
  else
    templatefilename = "SEP1.T"+ntrm+"L"+nlev+"."+grid+".template.nc";
  end if
  ; number of longitudes used For reduced grids only
  if ( grid .eq. "reduced" )then
    nlons = new( nlat, integer );
    nlons(:nlat-1) = (/10, 16, 30, 32, 40, 48, 50, 60, 64, 72, 80, 80, 90, \
    90, 96, 100, 108, 108, 120, 120, 120, 120, 128, 128, 128, 128, 128, 128, \
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 120, \
    120, 120, 120, 108, 108, 100, 96, 90, 90, 80, 80, 72, 64, 60, 50, 48, 40, \
    32, 30, 16, 10 /);
  end if
;-----------------------------------------------------------------------------------
; "Standard" definitions of hybrid coefficents for given level numbers.
; You only need to edit this if you need to go to a non-standard number of levels.
; Or if you need to change the levels for these cases.
;-----------------------------------------------------------------------------------
  var_type = "double";   ; Output type for variables creating
  hyai_standard= new( (/3,50/), var_type );
  hybi_standard = new( dimsizes(hyai_standard), var_type );
  nstandard = (/19,27,31/);  ; Sizes of hyai for each set
  hyai_standard(0,:nstandard(0)-1) =  \
  (/0.00251499470323319, 0.00710360519587988, 0.0190425943583249,  \
    0.046075608581305, 0.081818588078022, 0.0786980539560318,  \
    0.0746317282319069, 0.0695530474185939, 0.0633905678987499,  \
    0.0562177114188669, 0.0481529273092749, 0.039492335170507,  \
    0.0305845476686949, 0.0219334047287699, 0.01403667870909,  \
    0.0074586505070329, 0.0026468755677342, 0.0, 0.0/);
  hybi_standard(0,:nstandard(0)-1) = \
  (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0375697687268259, 0.0865262672305107,  \
    0.147671163082119, 0.22186444699764, 0.308222115039829, 0.40531820058823,  \
    0.509587526321408, 0.616832911968228, 0.720988512039177,  \
    0.81606125831604, 0.895257532596588, 0.953188955783837,  \
    0.985056042671197, 1./);
  ; Hybrid coefficients for standard 26 level case
  hyai_standard(1,:nstandard(1)-1) = \
  (/0.00219406700000001, 0.00489520900000001, 0.009882418, 0.01805201,  \
    0.02983724, 0.0446233400000002, 0.0616058700000002, 0.0785124300000004,  \
    0.0773127100000002, 0.0759013100000003, 0.0742408600000002,  \
    0.0722874400000002, 0.0699893299999998, 0.06728574, 0.06410509,  \
    0.0603632200000002, 0.0559611100000001, 0.0507822500000001,  \
    0.0446896000000001, 0.0375219099999999, 0.0290894900000001, 0.02084739,  \
    0.01334443, 0.00708499000000001, 0.00252136, 0.0, 0.0 /) ;

  hybi_standard(1,:nstandard(1)-1) = \
  (/0, 0, 0, 0, 0, 0, 0, 0, 0.01505309, 0.03276228, 0.05359622,  \
    0.0781062700000006, 0.1069411, 0.140863700000001, 0.180772, 0.227722,  \
    0.282956200000001, 0.347936400000002, 0.4243822, 0.514316800000003,  \
    0.620120200000002, 0.723535500000004, 0.817676800000001,  \
    0.896215300000001, 0.953476100000003, 0.9851122, 1.0 /) ;
  ; Hybrid coefficients for standard 30 level case
  hyai_standard(2,:nstandard(2)-1) = \
  (/0.00225523952394724, 0.00503169186413288, 0.0101579474285245,  \
    0.0185553170740604, 0.0306691229343414, 0.0458674766123295,  \
    0.0633234828710556, 0.0807014182209969, 0.0949410423636436,  \
    0.11169321089983, 0.131401270627975, 0.154586806893349,  \
    0.181863352656364, 0.17459799349308, 0.166050657629967,  \
    0.155995160341263, 0.14416541159153, 0.130248308181763,  \
    0.113875567913055, 0.0946138575673103, 0.0753444507718086,  \
    0.0576589405536652, 0.0427346378564835, 0.0316426791250706,  \
    0.0252212174236774, 0.0191967375576496, 0.0136180268600583,  \
    0.00853108894079924, 0.00397881818935275, 0, 0 /);
  hybi_standard(2,:nstandard(2)-1) = \
  (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0393548272550106,  \
    0.0856537595391273, 0.140122056007385, 0.204201176762581,  \
    0.279586911201477, 0.368274360895157, 0.47261056303978,  \
    0.576988518238068, 0.672786951065063, 0.753628432750702,  \
    0.813710987567902, 0.848494648933411, 0.881127893924713,  \
    0.911346435546875, 0.938901245594025, 0.963559806346893,  \
    0.985112190246582, 1 /);
;-----------------------------------------------------------------------------------
;  These are things that you won't have to worry about except with special grids
;  Set the dimensions, latitudes, longitudes, and gaussian weights.
;  All of this depends on the grid type. Right now it does the right thing
;  for "gaussian" , "reduced" and "staggered" grids.
;-----------------------------------------------------------------------------------
  FillValue = new( 1, var_type );
  FillValue = 1.e36;   ; Missing fill value to use
  ;
  ; Output dimension names and sizes (and the one that is unlimited
  ;
  nchar = 8;
  ntimes = 1;
  if ( (grid .eq. "gaussian") .or. (grid .eq. "reduced") )then
    dimnames = (/"lat", "lon", "lev", "ilev", "time", "chars" /);
    dsizes   = (/nlat,  nlon,  nlev,  nlev+1, ntimes,  nchar/);
    is_unlim = (/False, False, False, False,  True,  False/);
  else
    if ( grid .eq. "staggered" )then
      dimnames = (/"lat", "lon", "lev", "ilev", "slat", "slon", "time", "chars"  /);
      dsizes   = (/nlat,  nlon,  nlev,  nlev+1, nlat-1, nlon,   ntimes,  nchar /);
      is_unlim = (/False, False, False, False,  False,  False,  True,  False /);
    else
      print( "ERROR::Invalid grid type:"+grid );
      print( "Valid grids are: gaussian, reduced or staggered" );
      exit;
    end if
  end if
  ;
  ; Reduced grid
  ;
  if ( grid .eq. "reduced" )then
    rlon = new( (/nlat,nlon/), var_type );
  end if
  if ( grid .eq. "staggered" )then
    ;
    ; Staggered grid
    ;
    slat = new( nlat-1, var_type );
    slon = new( nlon, var_type );
    w_stag = new( nlat-1, var_type );
  end if
  gw   = new( nlat, var_type );
  lat  = new( nlat, var_type );
  lon = new( nlon, var_type );

  load "template.ncl";

  get_standard_lat_and_longs( );
  settings( );
  system( "/bin/rm " + templatefilename );
  print( "Create:"+templatefilename );
  nco = addfile( templatefilename, "c" );
  print( "From:"+interpfilename );
  nc = addfile( interpfilename, "r" );
  define_file( nco, nc );
  copy_vars_and_atts( nco, nc );
  set_vertical_levels( nco, nc );
  print( "Done!" );
end
