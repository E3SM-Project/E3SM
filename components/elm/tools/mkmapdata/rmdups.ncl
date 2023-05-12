;
; Remove duplicate weights from a mapping file.
;
;  Mark Taylor (converted for use by CLM mkmapdata by Erik Kluzek)
;  Sep/01/2011
;  $Id: rmdups.ncl 47629 2013-05-31 08:59:50Z erik $
;  $HeadURL;
;
load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"
begin
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE or use ENV VARIABLE SETTINGS
  ; Edit the following as needed
  ;
  ; Input mapping file to remove duplicate weights from a mapping file
  ;
  mapfile    = getenv("MAPFILE")         ; Get the mapping file
  newmapfile = getenv("NEWMAPFILE")      ; The new mapping file to create
  logname    = getenv("LOGNAME")         ; Logname of user running the script

  if ( ismissing(mapfile) )then
     print( "You did NOT enter an input mapping file to convert" )
     status_exit( -1 )
  end if
  if ( ismissing(newmapfile) )then
     sdate      = systemfunc( "date +%y%m%d" );
     newmapfile = mapfile+"_c"+sdate+".nc";
  end if
  ; ===========================================================================================================

  if ( systemfunc("test -f "+mapfile+"; echo $?" ) .ne. 0 )then
     print( "Input file does not exist or not found: "+mapfile );
     status_exit( -1 )
  end if
  print("map file: "+mapfile)
  f      = addfile(mapfile,"r")       ; Open netCDF files.	


  n_s = dimsizes(f->col)
  if ( n_s .eq. 0 )then
     print( "n_s is size zero, so no overlap points just return: " );
     exit
  end if

  n_b = dimsizes(f->area_b)
  n_a = dimsizes(f->area_a)
  print("n_s = "+n_s+" max(row)="+max(f->row)+" max(col)="+max(f->col))



  row = f->row
  col = f->col


  print("checking for dups, sorting...")
  hash = new( n_s, double )
  hash = col
  hash= hash + row*n_b 
  index1d=dim_pqsort(hash,1)
  row2=row(index1d)
  col2=col(index1d)
  S=f->S
  print("zeroing out any dups...")
  ndups=0
  i0=0
  do i=1,n_s-1
    if ( (col2(i) .eq. col2(i0)) .and. (row2(i) .eq. row2(i0))) then
        iorig1 = index1d(i0)
        iorig2 = index1d(i)
        ;print("dup row: "+row2(i)+" "+row2(i0)+" "+row(iorig1)+" "+row(iorig2))
        ;print("dup col: "+col2(i)+" "+col2(i0)+" "+col(iorig1)+" "+col(iorig2))
        ;print("removing "+iorig2+"  keeping "+iorig1)
        S(iorig1)=S(iorig1)+S(iorig2)
        S(iorig2)=0
        ndups=ndups+1
        ; dont increment i0
    else
        i0=i 
    end if
  end do
  delete(row2)
  delete(col2)
  if ( ndups .gt. 0) then
    print("ndups = "+ndups)  
    print("compacting S...")
    ns2 = n_s-ndups
    S2 = new( ns2, double)
    row2= new( ns2, integer)
    col2 = new( ns2, integer)
    ns2=0
    do i=0,n_s-1
      if (S(i) .ne. 0) then
        S2(ns2)=S(i)
        row2(ns2)=row(i)
        col2(ns2)=col(i)
        ns2=ns2+1
      end if
    end do
    print("removed "+ndups+" dups")
    delete(S)
    delete(row)
    delete(col)
    S=S2
    row=row2
    col=col2
    n_s = ns2
    print("writing new netcdf file")
    cmdout = systemfunc("ncks -O -x -v S,row,col  "+mapfile+" "+newmapfile)
    nco  = addfile(newmapfile,"rw")       ; Open netCDF files.	
    nco->S = S
    nco->row = row
    nco->col = col
    ldate       = systemfunc( "date" );
    nco@history = nco@history + ":"+ldate + ": ";
    nco@history = nco@history + " Removed duplicate weights from mapping file with: rmdups.ncl "
    nco@rmdups_Logname       = logname;
    nco@rmdups_mod_date      = ldate;
    nco@rmdups_version       = "$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/trunk_tags/clm4_5_1_r085/models/lnd/clm/tools/shared/mkmapdata/rmdups.ncl $";
    nco@rmdups_revision_id   = "$Id: rmdups.ncl 47629 2013-05-31 08:59:50Z erik $";

    print("Successfully removed duplicate weights from mapping file" );

  else

    print("No duplicate weights to remove from mapping file" );

  end if



end
