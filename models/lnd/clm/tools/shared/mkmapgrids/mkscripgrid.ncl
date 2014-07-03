;
;  mkscripgrid.ncl
;
;  Create SCRIP grid and mapping file for a land-only point or region.
;  Requires NCL 6.1.0 or later for the ESMF regridding functions
;
;  Erik Kluzek
;  Dec/07/2011
;  $Id: mkscripgrid.ncl 47629 2013-05-31 08:59:50Z erik $
;  $HeadURL;
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
begin
  ; ===========================================================================================================
  ; Set a few constants needed later
  cdate  = systemfunc( "date +%y%m%d" );
  ldate  = systemfunc( "date" );
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE or use ENV VARIABLE SETTINGS
  ; Edit the following as needed to interpolate to a new resolution.
  ;
  ; Input resolution and position
  ;
  name = getenv("PTNAME");                    ; Get name of this point

  latS = stringtodouble( getenv("S_LAT") );   ; Get south   latitude from env variable
  latN = stringtodouble( getenv("N_LAT") );   ; Get north   latitude from env variable
  lonE = stringtodouble( getenv("E_LON") );   ; Get east   longitude from env variable
  lonW = stringtodouble( getenv("W_LON") );   ; Get west   longitude from env variable

  nx   = stringtointeger( getenv("NX" )  );   ; Get number of grids along longitude lines
  ny   = stringtointeger( getenv("NY" )  );   ; Get number of grids along latitude lines

  imask = stringtointeger( getenv("IMASK") ); ; Get imask to use     from env variable

  print_str  = getenv("PRINT");               ; Do Extra printing for debugging

  outfilename = getenv("GRIDFILE");           ; Get filename from env variable

  if ( ismissing(nx) )then
     nx = 1;
  end if
  if ( ismissing(ny) )then
     ny = 1;
  end if
  if ( ismissing(imask) )then
     imask = 1;
  end if
  if ( ismissing(name) )then
     name = nx+"x"+ny+"pt_US-UMB";
  end if
  if ( ismissing(latS) )then
     latS = 45.5098;
  end if
  if ( ismissing(latN) )then
     latN = 45.6098;
  end if
  if ( ismissing(lonW) )then
     lonW = 275.2362;
  end if
  if ( ismissing(lonE) )then
     lonE = 275.3362;
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

  if ( ismissing(outfilename) )then
     if ( imask .eq. 1 )then
        outfilename = "SCRIPgrid_"+name+"_nomask_c"+cdate+".nc";
     else
        if ( imask .eq. 0 )then
           outfilename = "SCRIPgrid_"+name+"_noocean_c"+cdate+".nc";
        else
           outfilename = "SCRIPgrid_"+name+"_mask_c"+cdate+".nc";
        end if
     end if
  end if

  system( "/bin/rm -f "+outfilename );
  if ( printn )then
     print( "output file: "+outfilename );
  end if

  ;
  ; Compute derived quantities
  ;

  delX       = (lonE - lonW) / int2dble(nx);
  delY       = (latN - latS) / int2dble(ny);

  lonCenters = fspan( (lonW + delX/2.d0), (lonE - delX/2.d0), nx)
  latCenters = fspan( (latS + delY/2.d0), (latN - delY/2.d0), ny)
  lon = new( (/nx, ny/), "double" );
  lat = new( (/nx, ny/), "double" );
  lonCorners = new( (/nx, ny, 4/), "double" );
  latCorners = new( (/nx, ny, 4/), "double" );
  do i = 0, nx-1
    lat(i,:)          =  latCenters;
    latCorners(i,:,0) =  latCenters - delY/2.d0;
    latCorners(i,:,1) =  latCenters + delY/2.d0;
    latCorners(i,:,2) =  latCenters + delY/2.d0;
    latCorners(i,:,3) =  latCenters - delY/2.d0;
  end do
  do j = 0, ny-1
    lon(:,j)          =  lonCenters;
    lonCorners(:,j,0) =  lonCenters - delX/2.d0;
    lonCorners(:,j,1) =  lonCenters - delX/2.d0;
    lonCorners(:,j,2) =  lonCenters + delX/2.d0;
    lonCorners(:,j,3) =  lonCenters + delX/2.d0;
  end do

  ; for some reason, "No_FillValue" isn't working in the case where imask=1
  Mask2D = new( (/nx,ny/), "integer", "No_FillValue" )
  Mask2D(:,:) = imask
  gridSize = delX+"x"+delY

  ;
  ; Create SCRIP grid file
  ;
  
  Opt = True
  Opt@Mask2D   = Mask2D
  Opt@GridCornerLat = latCorners
  Opt@GridCornerLon = lonCorners
  Opt@Title = "SCRIP grid file for "+name
  if (printn) then
    Opt@Debug = True
  end if
  curvilinear_to_SCRIP(outfilename, lat, lon, Opt)

  ;
  ; Add global attributes to file
  ;

  nc = addfile( outfilename, "w" );
  nc@history = ldate+": create using mkscripgrid.ncl";
  nc@comment = "Ocean is assumed to non-existant at this point";
  nc@Version = "$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/branch_tags/cesm1_2_rel_tags/cesm1_2_0_n02_clm4_5_07/models/lnd/clm/tools/shared/mkmapgrids/mkscripgrid.ncl $";
  nc@Revision = "$Id: mkscripgrid.ncl 47629 2013-05-31 08:59:50Z erik $";
  if ( printn )then
    print( "================================================================================================" );
    print( "Successfully created SCRIP grid file: "+outfilename);
  end if

  ; ===========================================================================================================

end
