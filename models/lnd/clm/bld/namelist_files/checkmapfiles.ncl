;
; Check that the *_b values are the same between the mapping files
; at the same output resolution.
;
;  Erik Kluzek
;  Nov/18/2011
;  $Id: checkmapfiles.ncl 46006 2013-04-15 16:04:42Z sacks $
;  $HeadURL;
;

  print( "Check that datm mapping files are consistent" );
  resolutions = (/ "128x256", "64x128", "48x96", "32x64", "8x16", "94x192", "0.23x0.31", "0.47x0.63", "0.9x1.25", "1.9x2.5", "2.5x3.33", "4x5", "10x15", "5x5_amazon", "1x1_tropicAtl", "1x1_camdenNJ", "1x1_vancouverCAN", "1x1_mexicocityMEX", "1x1_asphaltjungleNJ", "1x1_brazil", "1x1_urbanc_alpha", "1x1_numaIA", "1x1_smallvilleIA", "ne4np4", "ne16np4", "ne30np4", "ne60np4", "ne120np4", "ne240np4" /);
 
  space  = "      ";
  badres = 0
  badresolutions = new( (/ 1000 /), string )
  chkres = 0
  chkresolutions = new( (/ 1000 /), string )

procedure checkit( desc:string, maxdiff:numeric, res:string, lmask:string, eps:numeric )
;
; check that difference is within reasonable tolerance...
;
begin
  reso = res+"_"+lmask;
  if ( maxdiff .gt. eps )then
     print( space+space+space+desc+" are off by more than tolerance for "+reso+" resolution" );
     print( space+space+space+"maximum difference = "+maxdiff );
     if ( .not. any(badresolutions .eq. reso ) )then
        badresolutions(badres) = reso;
        badres = badres + 1
     end if
  else
     print( space+space+space+"File OK for "+desc+"!" );
  end if
  if ( .not. any(chkresolutions .eq. reso ) )then
     chkresolutions(chkres) = reso;
     chkres = chkres + 1
  end if
end


function checkdims( desc:string, dsizefile1 [*]:integer, dsizefile2 [*]:integer, res:string, lmask:string )
;
; check that dimensions are the same between the file variables
;
begin
  reso = res+"_"+lmask;
  if ( any( dsizefile1 .ne. dsizefile2) )then
     print( space+space+space+desc+" dimensions are different for "+reso+" resolution" );
     print( space+space+space+"dim first file "+dsizefile1 );
     print( space+space+space+"dim second file "+dsizefile2 );
     if ( .not. any(badresolutions .eq. reso ) )then
        badresolutions(badres) = reso;
        badres = badres + 1
     end if
     return( False );
  else
     print( space+space+space+"File dims OK for "+desc+"!" );
     return( True );
  end if
  if ( .not. any(chkresolutions .eq. reso ) )then
     chkresolutions(chkres) = reso;
     chkres = chkres + 1
  end if
end

begin

  csmdata  = getenv("CSMDATA");
  clmroot  = getenv("CLM_ROOT");
  querynml = "bld/queryDefaultNamelist.pl -silent -justvalue -namelist clmexp";
  if ( .not. ismissing(csmdata) )then
     querynml = querynml+" -csmdata "+csmdata;
  end if
  if ( ismissing(clmroot) )then
     querynml = "../../"+querynml;
  else
     querynml = clmroot+"/models/lnd/clm*/"+querynml;
  end if

  print( "query string="+querynml )


  mapgrids = (/"0.5x0.5_MODIS", "0.5x0.5_AVHRR", "0.5x0.5_MODIS", "5x5min_nomask", "5x5min_IGBP-GSDP", "5x5min_ISRIC-WISE", "10x10min_nomask", "3x3min_MODIS", "3x3min_LandScan2004", "3x3min_GLOBE-Gardner", "3x3min_GLOBE-Gardner-mergeGIS", "0.9x1.25_GRDC", "360x720cru_cruncep", "1km-merge-10min_HYDRO1K-merge-nomask"/);
  do i = 0, dimsizes(resolutions)-1
     res = resolutions(i);
     print( "Go through maps for Resolution: "+res );
     do j = 0, dimsizes(mapgrids)-1
        grid = str_get_field( mapgrids(j), 1, "_" );
        lmask = str_get_field( mapgrids(j), 2, "_" );
        print( space+"Look for maps from Grid: "+grid+"_"+lmask);

        querynmlres = querynml+" -options frm_lmask="+lmask+",frm_hgrid="+grid+",to_hgrid="+res+",to_lmask=nomask";
        ;
        ; Get map filename and open it
        ;
        mapfile  = systemfunc( querynmlres+" -var map" );
        if ( systemfunc("test -f "+mapfile+"; echo $?" ) .ne. 0 )then
           delete( mapfile );
           continue;
        end if
        print( space+"Use mapfile:       "+mapfile );
        ncm     = addfile( mapfile,  "r" );
   
        if ( .not. isvar("ncm0") )then
           ncm0 = ncm;
        else
           vars = (/"yc_b", "xc_b", "area_b", "xv_b", "yv_b" /);
           k = 0;
           if ( checkdims( vars(k), dimsizes(ncm->$vars(k)$), dimsizes(ncm0->$vars(k)$), res, "nomask" ) )then
              do k = 0, dimsizes(vars)-1
                 maxdiff = max( abs(ncm->$vars(k)$ - ncm0->$vars(k)$) );
                 checkit( vars(k), maxdiff, res, "nomask", 1.e-12 );
                 delete( maxdiff );
              end do
              var = "mask_b"
              imaxdiff = max( abs(ncm->$var$ - ncm0->$var$) );
              checkit( var, imaxdiff, res, "nomask", 1.e-12 );
              delete( imaxdiff );
           end if
           delete( ncm );
        end if
        delete( mapfile );

     end do

     delete( grid  );
     delete( lmask );
     delete( res   );
     if ( isvar("ncm0") )then
        delete( ncm0  );
     end if

  end do
  ;
  ; go the other direction now check the _a variables
  ;
  mksrf_files = (/"mksrf_fvegtyp", "mksrf_fglacier", "mksrf_furbtopo", "mksrf_flndtopo", "mksrf_flai", "mksrf_fsoitex", "mksrf_fsoicol", "mksrf_ffrac", "mksrf_fmax", "mksrf_ftopo", "mksrf_firrig", "mksrf_forganic", "mksrf_flakwat", "mksrf_fwetlnd", "mksrf_furban", "mksrf_fvocef"/)
  do i = 0, dimsizes(mapgrids)-1
     grid = str_get_field( mapgrids(i), 1, "_" );
     lmask = str_get_field( mapgrids(i), 2, "_" );
     print( "Grid: "+grid);
     print( "Mask: "+lmask);
     do j = 0, dimsizes(resolutions)-1
        res = resolutions(j);
        print( "res: "+res );

        querynmlres = querynml+" -options frm_lmask="+lmask+",frm_hgrid="+grid+",to_hgrid="+res+",to_lmask=nomask";
        ;
        ; Get map filename and open it
        ;
        mapfile  = systemfunc( querynmlres+" -var map" );
        if ( systemfunc("test -f "+mapfile+"; echo $?" ) .ne. 0 )then
           delete( mapfile );
           continue;
        end if
        print( space+"Use mapfile:       "+mapfile );
        ncm     = addfile( mapfile,  "r" );
   
        if ( .not. isvar("ncm0") )then
           ncm0 = ncm;
        else
           vars  = (/"yc_a",   "xc_a",   "area_a", "xv_a", "yv_a" /);
           vars2 = (/"LATIXY", "LONGXY", "AREA" /);
           k = 0;
           if ( checkdims( vars(k), dimsizes(ncm->$vars(k)$), dimsizes(ncm0->$vars(k)$), res, "nomask" ) )then
              do k = 0, dimsizes(vars)-1
                 maxdiff = max( abs(ncm->$vars(k)$ - ncm0->$vars(k)$) );
                 checkit( vars(k), maxdiff, res, "nomask", 1.e-12 );
                 delete( maxdiff );
              end do
           end if
           var = "mask_a"
           imaxdiff = max( abs(ncm->$var$ - ncm0->$var$) );
           checkit( var, imaxdiff, res, "nomask", 1.e-12 );
           delete( imaxdiff );
           ;
           ; Get mksurfdata input datasets
           ;
           do k = 0, dimsizes(mksrf_files)-1
              srffile  = systemfunc( querynmlres+" -var "+mksrf_files(k) );
              if ( systemfunc("test -f "+srffile+"; echo $?" ) .ne. 0 )then
                 delete( srffile );
                 continue;
              end if
              print( space+"Use srffile:       "+srffile );
              ncs     = addfile( srffile,  "r" );
              n = 0;
              if ( checkdims( vars(n), dimsizes(ncm->$vars(n)$), ndtooned(dimsizes(ncs->$vars2(n)$)), res, "nomask" ) )then
                 do n = 0, dimsizes(vars2)-1
                    maxdiff = max( abs(ncm->$vars(n)$ - ndtooned(ncs->$vars2(n)$)) );
                    checkit( vars(n), maxdiff, res, "nomask", 1.e-12 );
                    delete( maxdiff );
                 end do
                 var  = "mask_a"
                 var2 = "LANDMASK"
                 imaxdiff = max( abs(ncm->$var$ - ndtooned(ncs->$var2$)) );
                 checkit( var, imaxdiff, res, "nomask", 1.e-12 );
              end if
              delete( ncs );
           end do
           delete( ncm );
        end if
        delete( mapfile );

     end do

     if ( isvar("vars") )then
        delete( vars  )
     end if
     if ( isvar("vars2") )then
        delete( vars2 )
     end if
     delete( grid  );
     delete( lmask );
     delete( res   );
     if ( isvar("ncm0") )then
        delete( ncm0  );
     end if

  end do
  if ( chkres .gt. 0 )then
     print( "resolutions checked = " );
     print( chkresolutions(0:chkres-1) );
  end if
  if ( badres .gt. 0 )then
     print( "badresolutions = " );
     print( badresolutions(0:badres-1) );
  end if

  print( "===============================" );
  print( "Successfully went through files" );

end

