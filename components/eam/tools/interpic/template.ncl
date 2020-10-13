;
;	template.ncl
;
;	Purpose: Subroutines to help create templates for interpolation.
;
;	Method: Set grid sizes and variables at top of
;		script. Copy variables and needed values from
;		file that will be interpolated.
;
;	Usage: Main program must set the following global variables:
;
;  nlat                 Number of latitudes
;  nlon                 Number of longitudes
;  nlev                 Number of levels
;  grid                 Grid type (staggered, reduced or gaussian)
;  ntrm                 Number of wave truncation for M (Gaussian and reduced grid only)
;  ntrk                 Number of wave truncation for K (Gaussian and reduced grid only)
;  ntrn                 Number of wave truncation for N (Gaussian and reduced grid only)
;  ntimes               Number of times on file
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
; Example usage:
;
;  get_standard_lat_and_longs( );
;  settings( );
;  nco = addfile( templatefilename, "c" );
;  nc = addfile( interpfilename, "r" );
;  define_file( nco, nc );
;  copy_vars_and_atts( nco, nc );
;  set_vertical_levels( nco, nc );
;
;	Author: Erik Kluzek
;
;	$Id$
;
undef("copy_VarAtts")
procedure copy_VarAtts(var_from,var_to)    
;
; Procedure to copy attributes from one variable to another.
;
local att_names, i
begin                                       
  att_names =getvaratts(var_from);
    if(.not.all(ismissing(att_names)))
      do i = 0,dimsizes(att_names)-1
         var_to@$att_names(i)$ = var_from@$att_names(i)$
      end do
  end if
end

undef( "get_standard_lat_and_longs" );
procedure get_standard_lat_and_longs( )
;
; Compute the standard latitudes and longitudes for the standard grid types
; gaussian, reduced, and staggered
;
begin
  if ( grid .eq. "reduced" )then
    rlon(:,0) = 0.0;
    rlon(:,1) = 360.0 / nlons;
    do j = 2, nlons(j)-1
      rlon(:,j) = rlon(:,j-1) + rlon(:,1);
    end do
  end if
  ;
  ; weights and latitudes
  ;
  if ( grid .eq. "staggered" )then
    info = linrood_latwgt( nlat );
    lat  = info(:,0)
    gw   = info(:,1)
    slat = lat(:nlat-2) + (lat(1) - lat(0))*0.5
    slat@long_name = "Latitude";
    slat@units = "degrees_north";
    ; Currently don't have an easy way to get the weights on the staggered grid
  else
    ;
    ; Gaussian grid
    ;
    gau_info = gaus(nlat/2)   ; divide by 2 to get "per hemisphere"
    lat = gau_info(:,0)  ; gaussian latitudes ( 1st dimension of gau_info)
    gw  = gau_info(:,1)  ; gaussian weights   ( 2nd dimension of gau_info)
  end if
  ;
  ; Longitudes
  ;
  lon(0) = 0.0;
  lon(1) = 360.0 / nlon;
  do i = 2, nlon-1
    lon(i) = lon(i-1) + lon(1);
  end do
  if ( grid .eq. "staggered" )then
    slon(:) = lon(:) - lon(1)*0.5;
    slon@long_name = "Longitude";
    slon@units = "degrees_east";
  end if
  lon@long_name = "Longitude";
  lon@units = "degrees_east";
  lat@long_name = "Latitude";
  lat@units = "degrees_north";
end

undef( "settings" );
procedure settings( )
;
; Print what the important settings are:
;
begin
  print( "This is creating a template for a "+grid+" grid." );
  print( "Input filename is: "+interpfilename );
  print( "Output filename is: "+templatefilename );
  print( "Case id is: "+caseid );
  if ( grid .ne. "staggered" )then
    print( "# lats: "+nlat+" # lons: "+nlon+" # levs: "+nlev+" Trunc M:"+ntrm+\
           " # times: "+ntimes + " Trunc N:"+ntrn+" Trunc K:"+ntrk );
  else
    print( "# lats: "+nlat+" # lons: "+nlon+" # levs: "+nlev+" Trunc M:"+ntrm+\
           " # times: "+ntimes );
  end if
  if ( grid .eq. "reduced" )then
    print( "Number of longitudes at each latitude: "+nlons );
  end if
end

undef( "grid_dimlist" );
function grid_dimlist( nco:file, nc:file, varname:string )
;
; Return the list of names for the desired dimensions for the given variable.
; On a staggered grid, U and V are handled differently
;
begin
  if ( grid .eq. "staggered" )then
    if ( varname .eq. "U" )then
      print( "U will be interpolating to staggered grid -- be sure to rename U to US after interpic" );
      dimlist = (/"slat", "lev", "lon"/);
    else
      if ( varname .eq. "V" )then
        print( "V will be interpolating to staggered grid -- be sure to rename V to VS after interpic" );
        dimlist = (/"lat", "lev", "slon"/);
      else
        dimlist = getfilevardims( nc, varname )
      end if
    end if
  else
     dimlist = getfilevardims( nc, varname )
  end if
  return( dimlist );
end

undef( "define_file" );
procedure  define_file( nco:file, nc:file )
;
; Define the output netCDF file variables and dimensions
;
begin
  print( "Define the file according to the settings" );
  print( "Dimension names are: "+dimnames );
  ;
  ; Define dimensions
  ;
  filedimdef ( nco, dimnames, dsizes, is_unlim );
  ;
  ; Check that needed dimensions and variables are on file
  ; 
  if ( .not. isdim( nc, "lat" ) )then
    print( "ERROR:: lat not a dimension on this file" );
    exit;
  end if
  if ( .not. isdim( nc, "lon" ) )then
    print( "ERROR:: lon not a dimension on this file" );
    exit;
  end if
  ;
  ; Define staggered lat and lon
  ;
  if ( grid .eq. "staggered" )then
    filevardef ( nco, "slat", var_type, "slat" );
    if ( isfilevar( nc, "lat" ) )then
      filevarattdef ( nco, "slat", nc->lat );
    else
      filevarattdef ( nco, "slat", lat );
    end if
    filevardef ( nco, "slon", var_type, "slon" );
    if ( isfilevar( nc, "lat" ) )then
      filevarattdef ( nco, "slon", nc->lon );
    else
      filevarattdef ( nco, "slon", lon );
    end if
  end if

  print( "now define the variables" );
  varnames = getfilevarnames( nc );
  do j = 0, dimsizes(varnames)-1
    i = dimsizes(varnames)-1 - j;
    ; NCL can't seem to handle string data
    if ( typeof(nc->$varnames(i)$) .ne. "char" )then
      if ( (grid .ne. "reduced") .or. (varnames(i) .ne. "lon") )then
        ;
        ; Get the Dimension names same as from interp grid 
        ; if staggered grid U and V are handled differently
        ;
        dimlist = grid_dimlist( nco, nc, varnames(i) );
        ; Copy scalars
        if ( dimlist(0) .eq. "ncl_scalar" )then
          nco->$varnames(i)$ =  nc->$varnames(i)$;
        ; Define vectors
        else 
          if ( (varnames(i) .ne. "slat") .and. (varnames(i) .ne. "slon") )then
            filevardef ( nco, varnames(i), typeof(nc->$varnames(i)$), dimlist );
            filevarattdef ( nco, varnames(i), nc->$varnames(i)$ );
            if ( typeof(nc->$varnames(i)$) .eq. var_type )then
              nco->$varnames(i)$@_FillValue = FillValue;
            end if
          end if
        end if
        delete( dimlist );
      end if
    else
    ; Copy string variables
      nco->$varnames(i)$ =  nc->$varnames(i)$;
    end if
  end do
  ;
  ; Define gw
  ;
  filevardef ( nco, "gw", var_type, "lat" );
  if ( isfilevar( nc, "gw" ) )then
    filevarattdef ( nco, "gw", nc->gw );
  else
    nco->gw@long_name = "gauss weights";
  end if
  if ( grid .eq. "staggered" )then
    filevardef ( nco, "w_stag", var_type, "slat" );
    nco->w_stag@long_name = "staggered latitude weights";
  end if
  ;
  ; Define reduced grid "rlon" and "nlon"
  ;
  if ( grid .eq. "reduced" )then
    filevardef ( nco, "rlon", var_type, (/"lat","lon"/) );
    filevardef ( nco, "nlon", "integer", "lat" );
    nco->nlon@long_name = "number of longitudes";
    nco->rlon@_FillValue = FillValue;
    if ( isfilevar( nc, "lon" ) ) then
      nco->rlon@long_name = grid + nc->lon@long_name;
      nco->rlon@units = nc->lon@units;
    else
      nco->rlon@long_name = grid + " Longitude";
      nco->rlon@units = "degrees_east";
    end if
  end if
end

undef( "copy_vars_and_atts" );
procedure  copy_vars_and_atts( nco:file, nc:file )
;
; Copy variables and attributes from the interpolation file
;
begin
  print( "Copy variables and attributes from the interpolation file: " + interpfilename );
  ;
  ; Set date and time
  ;
  ;if ( isfilevar( nco, "date_written" ) )then
  ;  delete( nco->date_written );
  ;end if
  ;if ( isfilevar( nco, "time_written" ) )then
  ;  delete( nco->time_written );
  ;end if
  ;date_written =  stringtochar( systemfunc( "date +%m/%d/%y" ) );
  ;time_written =  stringtochar( systemfunc( "date +%H:%M:%S" ) );
  ;do i = 0, nchar-1
  ;  nco->date_written(0,i) =  date_written(i);
  ;  nco->time_written(0,i) =  time_written(i);
  ;end do
  ;
  ; Copy reference pressure and variables having to do with time
  ;
  vlist = (/"P0", "time", "ndbase", "nsbase", "nbdate", "nbsec", "mdt", "ndcur", \
            "nscur", "date", "datesec", "nsteph" /);
  do i = 0, dimsizes(vlist)-1
    print( "Copy:"+vlist(i) );
    if ( isfilevar( nc, vlist(i) ) )then
      nco->$vlist(i)$ = (/nc->$vlist(i)$/);
    else
      print( "WARNING::Variable "+vlist(i)+" does not exist on "+interpfilename );
    end if
  end do
  ;
  ; Set ntrm, ntrk, ntrn
  ;
  if ( grid .ne. "staggered" )then
    nco->ntrm = ntrm;
    nco->ntrn = ntrn;
    nco->ntrk = ntrk;
  end if
  ;
  ; Copy lats, lon, and gw
  ;
  print( "copy grid" );
  nco->lat = (/lat/);
  if ( grid .eq. "reduced" )then
    nco->rlon = (/rlon/);
    nco->nlon = (/nlons/);
  else
    nco->lon = (/lon/);
    if ( grid .eq. "staggered" )then
      nco->slon = (/slon/);
      nco->slat = (/slat/);
      filevardef ( nco, "w_stag", var_type, (/"slat"/) );
      nco->w_stag = (/w_stag/);
    end if
  end if
  nco->gw  = (/gw/);
  ;
  ; Copy attributes
  ;
  print( "copy attributes" );
  ;fileattdef( nco, nc );
  globalAtt = True                 ; temporary to attach attributes                  
  copy_VarAtts (nc, globalAtt)                                                       
  copy_VarAtts (globalAtt, nco)                                                      
  ; Modify source title and caseid
  if ( isatt( nco, "source" ) ) then
    delete( nco@source );
  end if
  if ( isatt( nco, "case" ) ) then
    delete( nco@case );
  end if
  if ( isatt( nco, "title" ) ) then
    delete( nco@title );
  end if
  if ( isatt( nc, "source" ) ) then
    nco@source = "Interpolated from:" + interpfilename + "::" + nc@source;
  else
    nco@source = "Interpolated from:" + interpfilename;
  end if
  nco@case = caseid;
  if ( isatt( nc, "title" ) ) then
    nco@title = "Interpolated from:" + interpfilename + "::" + nc@title
  else
    nco@title = "Interpolated from:" + interpfilename;
  end if
end

undef( "set_vertical_levels" );
procedure  set_vertical_levels( nco:file, nc:file )
;
; Set the vertical levels
;
begin
  print( "Set the vertical levels" );
  if ( .not. isfilevar(nc, "lev") )then
    print( "lev does not exist on this file do not create vertical levels" );
    return
  end if
  if ( nlev .eq. dimsizes(nc->lev) )then
    if ( .not. isfilevar( nc, "hyai" ) )then
      print( "ERROR:: Variable hyai does not exist on the input file" );
      exit;
    end if
    if ( .not. isfilevar( nc, "hybi" ) )then
      print( "ERROR:: Variable hybi does not exist on the input file" );
      exit;
    end if
    hyai = (/nc->hyai/);
    hybi = (/nc->hybi/);
  else
    set = False;
    do i = 0, dimsizes(hybi_standard(:,0))-1
      if ( (.not. set) .and. (nlev .eq. nstandard(i)-1) )then
        hyai = hyai_standard(i,:nstandard(i)-1);
        hybi = hybi_standard(i,:nstandard(i)-1);
        set = True;
      end if
    end do
    if ( .not. set )then
      print( "ERROR:: Must add the hybi and hyai interfaces levels for nlev = "+nlev );
      print( "Edit the section of the script that lists values for hybi_standard" );
      exit;
    end if
  end if
  hyam = (hyai(0:nlev-1) + hyai(1:nlev) )*0.5;
  hybm = (hybi(0:nlev-1) + hybi(1:nlev) )*0.5;
  nco->hyam = hyam;
  nco->hybm = hybm;
  nco->hyai = hyai;
  nco->hybi = hybi;
  nco->lev = 1000.*(hyam+hybm);
  nco->ilev = 1000.*(hyai+hybi);
  print( "lev  = "+nco->lev );
  print( "ilev = "+nco->ilev );
end

