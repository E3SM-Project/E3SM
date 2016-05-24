;
; Take the greenhouse gas file used by CAM for historical (and future) representations of
; greenhouse gases, and convert it to a format that can be used by streams.
; So include domain data for a single point that covers the globe, as well
; as CO2 data over that single point. In the process we also discard the other
; greenhouse gases, as the datm can only pass CO2.
;
;  Erik Kluzek
;  Mar/03/2010
;  $Id: getco2_historical.ncl 46863 2013-05-08 03:25:21Z sacks $
;  $HeadURL;
;
begin
  ; ===========================================================================================================


  ; ===========================================================================================================
   ;
   ; Setup the namelist query script
   ;
   csmdata  = getenv("CSMDATA");
   clmroot  = getenv("CLM_ROOT");
   rcp      = getenv("RCP");          ; Get representative concentration pathway from env variable
   querynml = "bld/queryDefaultNamelist.pl -silent -justvalue ";
   if ( .not. ismissing(csmdata) )then
      querynml = querynml+" -csmdata "+csmdata;
   end if
   if ( ismissing(clmroot) )then
      querynml = "../../../"+querynml;
   else
      querynml = clmroot+"/models/lnd/clm/"+querynml;
   end if
   histrcp = -999.9
   if ( ismissing(rcp) )then
     rcp = histrcp
   end if
   ;
   ; Get input Greenhouse gas file and open it
   ;
   filetype  = "mkghg_bndtvghg";
   ghgfile  = systemfunc( querynml+" -namelist clmexp -var "+filetype+" -options rcp="+rcp );
   print( "Use "+filetype+" file: "+ghgfile );
   if ( systemfunc("test -f "+ghgfile+"; echo $?" ) .ne. 0 )then
      print( "Input "+filetype+" file does not exist or not found: "+ghgfile );
      exit
   end if
   ncg = addfile( ghgfile, "r" );

   ;
   ; Get date time-stamp to put on output CO2 file
   ;
   sdate     = systemfunc( "date +%y%m%d" );
   ldate     = systemfunc( "date" );

   sim_yr0 = ncg->date(0) / 10000;
   nyrs    = dimsizes( ncg->date );
   sim_yr2 = ncg->date(nyrs-1) / 10000;

   sim_yr_rng = sim_yr0 + "-" + sim_yr2;
   
   if ( rcp .eq. histrcp )then
      outco2filename = "fco2_datm_"+sim_yr_rng+"_c"+sdate+".nc";
   else
      outco2filename = "fco2_datm_rcp"+rcp+"_"+sim_yr_rng+"_c"+sdate+".nc";
   end if
   system( "/bin/rm -f "+outco2filename );
   print( "output file: "+outco2filename );
   nco = addfile( outco2filename, "c" );
   ;
   ; Define dimensions
   ;
   nlat = 1;
   nlon = 1;
   nv   = 4;
   dimnames = (/ "time", "lat", "lon", "nv" /);
   dsizes   = (/ nyrs, nlat,  nlon, nv /);
   is_unlim = (/ True, False, False, False /);
   filedimdef( nco, dimnames, dsizes, is_unlim );
   ;
   ; Define variables
   ;
   vars = (/ "lonc", "latc", "lonv", "latv", "mask", "frac", "area", "CO2" /);
   units= (/ "degrees_east", "degrees_north", "degree_east", "degrees_north", "unitless", "unitless",          "radians^2",        "ppmv" /);
   lname= (/ "Longitude of grid cell center", "Latitude of grid cell center", "Longitudesof grid cell vertices", "Latitudes of grid cell vertices", "Mask of active cells: 1=active", "Fraction of grid cell that is active", "Area of grid cell", "CO2 concentration" /);
   print( "Define variables: "+vars );
   do i= 0, dimsizes(vars)-1
      if ( vars(i) .eq. "lonv" .or. vars(i) .eq. "latv" )then
         filevardef ( nco, vars(i), "double",  (/ "lat", "lon", "nv" /) );
      else
         if ( vars(i) .eq. "CO2" )then
            filevardef ( nco, vars(i),  "float",   (/ "time", "lat", "lon" /) );
            nco->$vars(i)$@coordinate  = "latc lonc time";
         else
            filevardef ( nco, vars(i), "double",  (/ "lat", "lon" /) );
         end if
      end if
      nco->$vars(i)$@units = units(i);
      nco->$vars(i)$@lname = lname(i);
   end do
   varstatic = (/ "mask", "frac", "area" /);
   do i = 0, dimsizes(varstatic)-1
      nco->$varstatic(i)$@coordinate  = "latc lonc";
   end do
   nco->lonc@bounds      = "lonv";
   nco->latc@bounds      = "latv";
   ;
   ; Add attributes
   ;
   fileattdef ( nco, ncg );
   nco@history  = ldate+": Convert by getco2_historical.ncl";
   nco@source   = "Convert from:"+ghgfile;
   nco@Version  = "$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/trunk_tags/clm4_5_1_r085/models/lnd/clm/tools/shared/ncl_scripts/getco2_historical.ncl $";
   nco@Revision = "$Id: getco2_historical.ncl 46863 2013-05-08 03:25:21Z sacks $";
   ;
   ; Set static variables
   ;
   pi                   = 3.14159265358979323846d00;
   nco->area            = 4.0*pi;
   nco->mask            =   1;
   nco->frac            =   1.0;
   nco->latv(0,0,0:1)   =  90.0;
   nco->latc            =   0.0;
   nco->latv(0,0,2:3)   = -90.0;
   nco->lonv(0,0,0:3:3) =   0.0;
   nco->lonc            = 180.0;
   nco->lonv(0,0,1:2)   = 360.0;
   ;
   ; Time and date
   ;
   nco->time = ncg->time;
   nco->date = ncg->date;
   nco->date@comment = "This variable is NOT used when read by datm, the time coordinate is used";
   ;
   ; CO2
   ;
   print( "Copy CO2 for "+nyrs+" years of data" );
   nco->CO2(:,0,0) = (/ ncg->CO2(:) /);

   print( "================================================================================================" );
   print( "Successfully created output historical CO2 file: "+outco2filename);

end
