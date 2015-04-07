;
; Extract out regional datasets needed to run clm from the global datasets.
; NOTE: Requires at least NCL version 5.1.0 or later...
;
;  Erik Kluzek
;  Aug/28/2009
;  $Id: getregional_datasets.ncl 65547 2014-11-18 20:34:35Z erik $
;  $HeadURL;
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl";

procedure  getfilecoord_namenlen( filenames[*]:string, dimnames[*]:string, dimlens[*]:integer, nlen:integer, name:string )
;
; get the name and size of either the latitude or longitude
;
  local d, l
begin
       if ( name .eq. "" )then
         do d = 0, dimsizes(filenames)-1
            if ( any(dimnames .eq. filenames(d) ) )then
              name = filenames(d)
              ; Get length of this dimension
              do l = 0, dimsizes(dimnames)-1
                 if ( dimnames(l) .eq. name )then
                    nlen = dimlens(l)
                 end if
              end do
            end if
         end do
       end if
end

begin
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE or use ENV VARIABLE SETTINGS
  ; Edit the following as needed to interpolate to a new resolution.
  ;
  ; Input resolution and position
  ;
  latS        = stringtodouble( getenv("S_LAT") );   ; Get south latitude from env variable
  latN        = stringtodouble( getenv("N_LAT") );   ; Get north latitude from env variable
  lonE        = stringtodouble( getenv("E_LON") );   ; Get east longitude from env variable
  lonW        = stringtodouble( getenv("W_LON") );   ; Get west longitude from env variable
  debug_str   = getenv("DEBUG");                     ; Don't run just -- debug
  print_str   = getenv("PRINT");                     ; Do Extra printing for debugging
  gridfile    = getenv("GRIDFILE");                  ; Input global grid file
  nfiles      = stringtointeger( getenv("NFILES") ); ; number of files to read in file lists
  filelistfil = getenv("INFILELIST");                ; filename of list of global files to work on
  regfilelistfil = getenv("OUTFILELIST");            ; filename of list of regional eiles to create

  if ( ismissing(nfiles) )then
     print( "NFILES is missing -- need to provide the number of files to process" );
     status_exit( -1 )
  end if
  if ( ismissing(filelistfil) .or. ismissing(regfilelistfil) )then
     print( "INFILELIST or OUTFILELIST is missing -- need to provide both" );
     status_exit( -1 )
  end if
  if ( ismissing(latS) )then
     latS = 52.0d00;
  end if
  if ( ismissing(latN) )then
     latN = 73.0d00;
  end if
  if ( ismissing(lonW) )then
     lonW = 190.0d00;
  end if
  if ( ismissing(lonE) )then
     lonE = 220.0d00;
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
  if ( ismissing(debug_str) )then
     debug = False;
  else
     if ( debug_str .eq. "TRUE" )then
        print( "DEBUG is TRUE do extra printing AND do NOT execute -- just print what WOULD happen" );
        debug  = True;
        printn = True;
     else
        debug = False;
     end if
  end if
  print( "Extract out regional datasets from global datasets" );
  if ( printn .eq. True )then
    print( "Regional:  Latitude="+latS+"-"+latN+" Longitude="+lonW+"-"+lonE );
  end if

  ;
  ; Setup the namelist query script
  ;
  ldate     = systemfunc( "date" );
  clmroot   = getenv("CLM_ROOT");

  ;
  ; list of latitude and longitude names
  ;
  filelatnames = (/ "lsmlat", "lat", "nj" /);
  filelonnames = (/ "lsmlon", "lon", "ni" /);

  ;
  ; Open file
  ;
  if ( systemfunc("test -f "+gridfile+"; echo $?" ) .ne. 0 )then
     print( "Input gridfile does not exist or not found: "+gridfile );
     status_exit( -1 )
  end if
  if ( printn .eq. True )then
     print( "gridfile:"+gridfile );
  end if
  ncg     = addfile( gridfile,  "r" );
  ;
  ; Get the names for latitude/longitude on the grid file
  ;
  varnames = getfilevarnames( ncg );
  gridlonnm = ""
  gridlatnm = ""
  glat      = 0
  glon      = 0
  varlens  = new( dimsizes(varnames), "integer" );
  getfilecoord_namenlen( (/ "yc", "LATIXY"/), varnames, varlens, glat, gridlatnm );
  getfilecoord_namenlen( (/ "xc", "LONGXY"/), varnames, varlens, glon, gridlonnm );
  delete( varnames );
  delete( varlens  );
  if ( gridlatnm .eq. "" )then
     print( "Could not find a recognizable latitude dimension name" )
     status_exit(-1);
  end if
  if ( printn .eq. True )then
     print( "gridlatname = "+gridlatnm )
     print( "gridlonname = "+gridlonnm )
  end if

  gridlon = ncg->$gridlonnm$;
  gridlon = where( gridlon < 0.0, 360.0 + gridlon, gridlon );

  indx = region_ind ( (/ncg->$gridlatnm$/), (/gridlon/), latS, latN, lonW, lonE );
  ; Indexes into indices
  ilat0 = 0;
  ilatN = 1;
  ilon0 = 2;
  ilonN = 3;

  latdim = dimsizes(ncg->$gridlatnm$(:,0))
  londim = dimsizes(gridlon(0,:))
  if ( any( ismissing(indx)) )then
     print( "Indices:"+indx );
     print( "Missing indices found" );
     print( "nlat: "+latdim );
     print( "nlon: "+londim );
     print( "yc: "+ncg->$gridlatnm$(:,0) );
     print( "xc: "+gridlon(0,:) );
     status_exit(-1);
  end if

  if ( debug .eq. True )then
     print( "Indices:"+indx );
  end if
  if ( printn .eq. True )then
     print( "Full grid size: nlat = "+latdim+" nlon = "+londim )
     loclatdim = indx(ilatN) - indx(ilat0) + 1;
     loclondim = indx(ilonN) - indx(ilon0) + 1;
     print( "Grid size:"+loclatdim+"x"+loclondim );
     LOLAT = ncg->$gridlatnm$(indx(ilat0),indx(ilon0));
     HILAT = ncg->$gridlatnm$(indx(ilatN),indx(ilonN));
     print( "Actual grid span: Latitude="+LOLAT+"-"+HILAT );
     LOLON = gridlon(indx(ilat0),indx(ilon0));
     HILON = gridlon(indx(ilatN),indx(ilonN));
     print( "Actual grid span: Longitude="+LOLON+"-"+HILON );
  end if

  ;
  ; Read in the list of files
  ;
  filelist    = asciiread(filelistfil(0),    (/ nfiles /), "string");
  regfilelist = asciiread(regfilelistfil(0), (/ nfiles /), "string");
  ;
  ; Loop over each of the files to process...
  ;
  do i = 0, nfiles-1
     ;
     ; Get the filename of the input global file and the output regional filename
     ;
     globalfile = filelist(i)
     if ( systemfunc("test -f "+globalfile+"; echo $?" ) .ne. 0 )then
        print( "Input global "+globalfile+" file does not exist or not found: "+globalfile );
        status_exit(-1);
     end if
     if ( debug .eq. True )then
        print( "Process file: "+globalfile );
     end if
     regfile = regfilelist(i)
     if ( ismissing(regfile) )then
        print( "Output regional filename was NOT found: "+regfile );
        status_exit(-1);
     end if
     
     nc = addfile( globalfile, "r" );
     varnames = getfilevarnames( nc );
     filelonnm = ""
     filelatnm = ""
     nlat      = 0
     nlon      = 0
     do v = 0, dimsizes(varnames)-1
       dimnames = getfilevardims(     nc, varnames(v) );
       dimlens  = getfilevardimsizes( nc, varnames(v) );
       getfilecoord_namenlen( filelatnames, dimnames, dimlens, nlat, filelatnm );
       getfilecoord_namenlen( filelonnames, dimnames, dimlens, nlon, filelonnm );
       delete( dimnames );
       delete( dimlens  );
     end do
     if ( filelatnm .eq. "" )then
        print( "Could not find a recognizable latitude dimension name" )
        status_exit(-1);
     end if
     if ( printn .eq. True )then
        print( "nlat = "+nlat+" nlon = "+nlon )
     end if
     ;
     ; Check to make sure number of latitudes and longitudes are the same as on the domain file
     ;
     if ( (latdim .ne. nlat) .or. (londim .ne. nlon) )then
        print( "Latitude or longitude dimensions do NOT match the grid file for file: "+globalfile );
        status_exit(-1);
     end if
     ;
     ; Run ncks on it over the region of interest
     ;
     do v = 0, dimsizes(varnames)-1
        cmd = "ncks -O -d "+filelatnm+","+indx(ilat0)+","+indx(ilatN)+" -d "+filelonnm+","+indx(ilon0)+","+indx(ilonN);
        cmd = cmd + " -v " + varnames(v) + " " + globalfile + " "+regfile+"_VAR"+varnames(v)+".nc"
        print( "Execute:"+cmd );
        if ( debug .eq. False )then
           if (  systemfunc( cmd+"; echo $?" ) .ne. 0 )then
              print( "Command did not complete successfully: " );
              status_exit( -1 )
           end if
        end if
        cmd = "ncks -A  "+regfile+"_VAR"+varnames(v)+".nc "+regfile
        print( "Execute:"+cmd );
        if ( debug .eq. False )then
           if (  systemfunc( cmd+"; echo $?" ) .ne. 0 )then
              print( "Command did not complete successfully: " );
              status_exit( -1 )
           end if
           system( "/bin/rm "+regfile+"_VAR"+varnames(v)+".nc" )
        end if
     end do
     delete( varnames );
     if ( debug .eq. False )then
        ;
        ; Open up resultant file for writing
        ;
        nco = addfile( regfile, "w" );
        nco@history = nco@history + ":"+ldate + ": ";
     end if
  end do

  print( "================================================================================================" );
  print( "Successfully created regional datasets from global datasets" );

end
