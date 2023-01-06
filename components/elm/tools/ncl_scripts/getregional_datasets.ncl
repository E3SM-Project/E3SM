;
; Extract out regional datasets needed to run clm from the global datasets.
; NOTE: Requires at least NCL version 5.1.0 or later...
;
;  Erik Kluzek
;  Aug/28/2009
;  $Id: getregional_datasets.ncl 62740 2014-08-16 20:58:21Z erik $
;  $HeadURL;
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl";
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
  filelist    = getenv("INFILELIST");                ; Comma  delimited list of global files to work on
  regfilelist = getenv("OUTFILELIST");               ; Comma  delimited list of regional eiles to create
  delim       = ","                                  ; Delimitor for file list

  if ( ismissing(filelist) .or. ismissing(regfilelist) )then
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

  indx = region_ind ( (/ncg->yc/), (/ncg->xc/), latS, latN, lonW, lonE );

  latdim = dimsizes(ncg->yc)
  londim = dimsizes(ncg->xc)
  if ( any( ismissing(indx)) )then
     print( "Indices:"+indx );
     print( "Missing indices found" );
     print( "yc: "+ncg->yc );
     print( "xc: "+ncg->xc );
     exit
  end if

  if ( debug .eq. True )then
     print( "Indices:"+indx );
  end if
  if ( printn .eq. True )then
     latdim = indx(3) - indx(2) + 1;
     londim = indx(1) - indx(0) + 1;
     print( "Grid size:"+latdim+"x"+londim );
     LOLAT = ncg->yc(indx(2),indx(0));
     HILAT = ncg->yc(indx(3),indx(1));
     print( "Actual grid span: Latitude="+LOLAT+"-"+HILAT );
     LOLON = ncg->xc(indx(2),indx(0));
     HILON = ncg->xc(indx(3),indx(1));
     print( "Actual grid span: Longitude="+LOLON+"-"+HILON );
  end if

  filelatnames = (/ "lsmlat", "lat", "nj" /);
  filelonnames = (/ "lsmlon", "lon", "ni" /);

  if ( str_fields_count( filelist, delim ) .ne. str_fields_count( regfilelist, delim ) )then
     print( "Number of files in global file list does NOT match the regional file list!" );
     print( "Global   file list: "+filelist );
     print( "Regional file list: "+regfilelist );
     status_exit(-1);
  end if
  ;
  ; Loop over each of the files to process...
  ;
  do i = 1, str_fields_count( filelist, delim )
     ;
     ; Get the filename of the input global file and the output regional filename
     ;
     globalfile = str_get_field( filelist, i, delim )
     if ( systemfunc("test -f "+globalfile+"; echo $?" ) .ne. 0 )then
        print( "Input global "+globalfile+" file does not exist or not found: "+globalfile );
        status_exit(-1);
     end if
     if ( debug .eq. True )then
        print( "Process file: "+globalfile );
     end if
     regfile = str_get_field( regfilelist, i, delim )
     if ( ismissing(regfile) )then
        print( "Output regional filename was NOT found: "+regfile );
        status_exit(-1);
     end if
     
     nc = addfile( globalfile, "r" );
     varnames = getfilevarnames( nc );
     filelonnm = ""
     filelatnm = ""
     do v = 0, dimsizes(varnames)-1
       if ( filelatnm .eq. "" )then
          dimnames = getfilevardims( nc, varnames(v) );
          do d = 0, dimsizes(filelatnames)-1
             if ( any(dimnames .eq. filelatnames(d) ) )then
               filelatnm = filelatnames(d)
               filelonnm = filelonnames(d)
             end if
          end do
          delete( dimnames );
       end if
     end do
     if ( filelatnm .eq. "" )then
        print( "Could not find a recognizable latitude dimension name" )
        status_exit(-1);
     end if
     ;
     ; Run ncks on it over the region of interest
     ;
     do v = 0, dimsizes(varnames)-1
        cmd = "ncks -O -d "+filelatnm+","+indx(0)+","+indx(1)+" -d "+filelonnm+","+indx(2)+","+indx(3);
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
