;***************************************************************
; NCL script
; 
; Apply SCRIP format mapping file to regrid a list of netcdf files
;
; latest version available on github:  
;    git@github.com:mt5555/remap-ncl.git
; download directly via:
;    wget https://raw.github.com/mt5555/remap-ncl/master/regridfile.ncl
;
;
; based on regridclimo.ncl, John Truesdale, May 2012
; Mark Taylor, April 2014:  user specified map file, process wildcard filenames
; Mark Taylor, November 2014:  support remapping of lat/lon input files
; Mark Taylor, December 2014:  Added srcfile name and weightfile name to attributes
; Mark Taylor, January 2015:   Specify wgt and src files on NCL command line
; Mark Taylor, Februrary 2015: Large file support
; Mark Taylor, August 2015:    Process files without lev/ilev information
;
; NCL script to apply remap NETCDF files using an ESMF (SCRIP-format) mapping file
; 
; ncl 'wgtfile="map_file_name"' 'srcfile="*cam.h0.*.nc"'  regridfile.ncl
; 
; or edit this file below and set:
;   wgtFileName =  Name of ESMF mapping file to use
;   srcFileName =  file (or list of files with wildcards) to be remapped
;
; This script will then remap all appropriate variables in each file, creating
; new files with the original filename and a "-remap.nc" suffix.
;
;***************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"


;***********************************************************************;
; detect type of lat/lon grid
; used to determine proper latitude weights
;
; gauss_grid =-1    routine failed 
; gauss_grid = 1    Gauss grid (use shaeg routines)
; gauss_grid = 2    Equal angle grid with poles  (use shaec routines)
; gauss_grid = 3    Equal angle grid without poles  (must interpolate first)
;
;***********************************************************************;
function get_gridtype(lat)
local d1,d2
begin
  gauss_grid = -1
  print("lat(0) = "+lat(0))
  if (lat(0) .lt. -89.99999 ) then
    gauss_grid = 2    ; use shaec routines
  else
    d1=lat(1)-lat(0)
    d2=lat(2)-lat(1)
    if (abs(d1-d2) .lt. 1e-8 ) then
      gauss_grid = 3    ; need to interpolate using fo2fsh
    else
      print("detected Gauss grid")
      gauss_grid = 1  
    end if
  end if
  return(gauss_grid)
end


begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set the mapping file and the map_method attribute
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (isvar("wgtfile")) then    
  wgtFileName = wgtfile
  print("wgtfile = "+wgtFileName)
else
  ;wgtFileName = "~/cam/remap/map_ne240np4_to_fv0.1x0.1_US.nc"
  ;wgtFileName = "~/codes/mapping/maps/map_ne30np4_to_201x400_bilin.nc"
  ;wgtFileName = "~/codes/mapping/maps/map_ne240np4_to_fv0.35x0.35_aave.nc"
  wgtFileName = "~/codes/mapping/maps/map_fv1793x3584_to_fv1441x2880_aave.nc"
  print("default wgtfile = "+wgtFileName)
end if



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set files to remap.  wildcards are ok
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (isvar("srcfile")) then    
  srcFileName = systemfunc("ls "+srcfile)
else
  ;srcFileName = systemfunc("ls *climo.nc")
  ;srcFileName = systemfunc("ls h0-tavg.nc")
  ;srcFileName = systemfunc("ls *cam.h2.000[2345]-??-??-00000.VORDIV.nc")
  srcFileName = systemfunc("ls *cam.h2.0004-??-??-00000.VORDIV.nc")
end if


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set this if the mapping file has tweaked pole points so the are 
; not at exactly +/-90 deg.  Often done in older, SCRIP mapping files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
fixpoles = 0;


;--- Specify a list of 1D variables on the *source SE file* that should NOT be copied
var_in_exclude = (/"gw", "lat", "lon", "lat_d", "lon_d", "area_d", "area", "date_written", "time_written"/)  

;--- Specify a list of variables on the source SE file that should be directly copied 
var_in_copy = (/"time_bnds"/)  

;--- Specify a list of variables to be regridded
var_out        = "All_Variables"     ; to be regridded
;var_out        = (/"VOR850","U850"/)     ; to be regridded




; read input files
print("adding files: "+srcFileName)
src_file_list=addfiles(srcFileName,"r")
num_files=dimsizes(srcFileName)
print("number of files added: "+num_files)



;---Read from the weight file the method used to derive the remap weights
wgt_file       = addfile(wgtFileName, "r")

if ( dimsizes(wgt_file->dst_grid_dims) .ne. 2 ) then
   print("Error: destination grid is not of rank 2")
   print("This script currently requires destination grid to be of rank 2")
   exit
end if

dst_grid_dims  = wgt_file->dst_grid_dims
dst_mlon       = dst_grid_dims(0)
dst_nlat       = dst_grid_dims(1)

dst_lat        = wgt_file->yc_b(::dst_mlon)
dst_lon        = wgt_file->xc_b(:dst_mlon-1)


if (fixpoles .eq. 1 ) then
  dst_lat(0) = -90.0d0
  dst_lat(dst_nlat-1) = 90.0d0
end if



;---Use the destination (FV) grid info on the weight file to create lat/lon
lat       = dst_lat                  ; get from weight file
lat@long_name = "latitude"
lat!0     = "lat"
lat@units = "degrees_north"
lat&lat   =  lat
nlat      = dimsizes(lat)            ; same as dst_nlat

lon       = dst_lon
lon@long_name = "longitude"
lon!0     = "lon"
lon@units = "degrees_east"
lon&lon   =  lon
mlon      = dimsizes(lon)            ; same as dst_mlon
print("interpolating to output grid of size: "+nlat+" x "+mlon)
print("dst_lat range: "+dst_lat(0)+" "+dst_lat(nlat-1))
print("dst_lon range: "+dst_lon(0)+" "+dst_lon(mlon-1))
;===================================================================
; generate lat weights 
;===================================================================
gridtype= get_gridtype(dst_lat)
;gridtype = 2
;print("Overriding grid detection. hardcoding gridtype=2")
if ( gridtype .eq. 3  ) then
  print("Detected FV offset grid")
  gw  = NormCosWgtGlobe(dst_lat)
end if
if ( gridtype .eq. 2  ) then
  print("Detected FV-style grid: generating weights with latRegWgt() " + dst_lat(0))
  gw  = latRegWgt(dst_lat, "double", 0)
end if
if ( gridtype .eq. 1) then
  print("Detected Gauss grid - generating Gauss weights "+dst_lat(0))
  nlat_w = dimsizes(dst_lat)
  nlat_w@double = True
  gw  = latGauWgt(nlat_w, "lat", "gaussian weights", "")
end if
gw!0= "lat"






do fk=0,num_files-1
    src_file = src_file_list[fk]

    dstFileName = str_sub_str(srcFileName(fk),".nc","-remap.nc")
    if (dstFileName .eq. srcFileName(fk)) then
       print("Error forming dstFileName. exiting...")
       exit
    end if

    print(" ")
    print("Regridding: fk="+(fk+1)+"/"+num_files)
    print("SRC:"+srcFileName(fk))
    print("DST:"+dstFileName)
    if(isfilepresent(dstFileName)) then
      system("rm "+dstFileName)
    end if
    setfileoption("nc","Format","LargeFile") 
    dst_file=addfile(dstFileName,"c")


    ;---Get all variables on the SE file 
    var_in   = getfilevarnames( src_file )
    nvar_in  = dimsizes(var_in)
    Opt_RGRD  = True


;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;                    PREDEFINE MODE
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    setfileoption(dst_file,"DefineMode",True)

;===================================================================
; create global attributes of the netCDF file containing regridded data
;===================================================================
    global    = True
    global@separator2= "---------------------------------------"
    copy_VarAtts(src_file, global)
    global@separator1= "------- SOURCE FILE ATTRIBUTES --------"

    if (isatt(src_file,"title")) then
        global@TITLE = "REMAPPED: "+src_file@title
    end if

    global@remap         = "NCL: ESMF_regrid_with_weights"
    global@srcfile       = srcFileName(fk)
    global@weightfile    = wgtFileName
    global@creation_date = systemfunc("date")
    
    fileattdef( dst_file, global )        ; copy file attributes to netCDF file

;===================================================================
; predefine the coordinate dimension names and their sizes
;===================================================================
    dNames   = getvardims(src_file)       ; get level info from source file
    if (any(dNames.eq."ncol")) then 
      print("Found ncol in dimension list.  Assuming unstructured source file")
      unstructured_input=True   ; set for cubed-sphere native grid data
    else
      print("Did not find ncol in dimension list.  Assuming structured source file")
      unstructured_input=False   ; set for lat/lon input data
    end if
    dSizes   = getfiledimsizes(src_file)
    if (any(dNames.eq."lev")) then 
      ia       = ind(dNames.eq."lev" )   ; # of levels 
      ib       = ind(dNames.eq."ilev")
      klev     = dSizes(ia)
      klevi    = dSizes(ib)
    else
      klev     = 1
      klevi    = 2
    end if
                                       ; standard CAM dimensions
    dimNames = (/"time",    "lat",    "lon", "lev",  "ilev", "nbnd", "chars"/)  
    dimSizes = (/ -1   ,    nlat ,    mlon ,  klev,  klevi ,   2   ,    8   /) 
    dimUnlim = (/ True ,    False,    False, False,  False , False , False  /)   
    filedimdef(dst_file,dimNames,dimSizes,dimUnlim)



;--- The following are explicitly added because they are 'special'
   filevardef(dst_file, "lat", typeof(lat), getvardims(lat)) 
   filevarattdef(dst_file,"lat" ,lat)           

   filevardef(dst_file, "lon", typeof(lon), getvardims(lon)) 
   filevarattdef(dst_file,"lon" ,lon)           

   filevardef(dst_file, "gw", typeof(gw), getvardims(gw)) 
   filevarattdef(dst_file,"gw" ,gw)    


;--- Loop over all variables and predfine meta data
   regridvar = new(nvar_in, logical)

   do nv=0,nvar_in-1
     if (.not.any(var_in(nv).eq.var_in_exclude)) then 
       if(var_out(0).eq."All_Variables"          .or. \
         any(var_in(nv).eq.var_out)             .or. \
         any(var_in(nv).eq.var_in_copy)       ) then

         rank_in = dimsizes(getfilevardimsizes(src_file, var_in(nv)))
         vdims = getfilevardims(src_file, var_in(nv))
         
         regridvar(nv) = False
         if ( any(var_in(nv).eq.var_in_copy) ) then
           ; USER SPECIFIED COPY - dont regrid
         else
           if (unstructured_input) then
             if (any(vdims.eq."ncol")) then 
               regridvar(nv)=True
             end if
           else
             if (any(vdims.eq."lat") .and. any(vdims.eq."lon")) then 
               regridvar(nv)=True
             end if
           end if
         end if
           
         print("defining: var="+var_in(nv)+" rank="+rank_in+" remap="+regridvar(nv))

         if ( .not. regridvar(nv) ) then
           ; copy varaible, no interpolation  
           filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
           , getfilevardims(src_file,var_in(nv))  )
         else 
           if (unstructured_input) then
             ; UNSTRUCTURED DATA (remapping ncol dimension to lat/lon)
             ; special variables of size ncol (rare)
             if (rank_in .eq.1 ) then
               filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
               , (/ "lat", "lon" /)        )
             end if
             ; 2D variables, like PS(ncol,time)
             if (rank_in .eq.2 ) then
               filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
               , (/ "time", "lat", "lon" /)        )
             end if
              ; 3D variables, like T(ncol,nlev,time)
             if (rank_in .eq.3) then
               filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
               , (/ "time", vdims(1),  "lat", "lon" /)        )
             end if
             if (rank_in .ge.4) then
                print("Error: not coded for rank 4 ncol variables - please fix")
                exit
             end if
           else
             ; STRUCTURED DATA (remapping lat/lon dimension to lat/lon)
             ; 2D variables, no time
             if (rank_in .eq.2 ) then
               filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
               , (/ "lat", "lon" /)        )
             end if
                                ; 2D variables, like PS(nlat,nlon,time)
             if (rank_in .eq.3 ) then
               filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
               , (/ "time", "lat", "lon" /)        )
             end if
                                ; 3D variables, like T(nlat,nlon,nlev,time)
             if (rank_in .eq.4) then
               filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
               , (/ "time", vdims(1),  "lat", "lon" /)        )
             end if
             if (rank_in .ge.5) then
                print("Error: not coded for rank 5 lat/lon variables - please fix")
                exit
             end if
           end if  ; unstructured_input
         end if ; not regrdivar
         
         dumAtts = new( 1, getfilevartypes(src_file,var_in(nv)))
         varAtts = getfilevaratts(src_file, var_in(nv))
         if (.not.ismissing(varAtts(0))) then
           nAtts = dimsizes(varAtts)
           do na=0,nAtts-1
             dumAtts@$varAtts(na)$ = src_file->$var_in(nv)$@$varAtts(na)$
           end do
           filevarattdef(dst_file, var_in(nv) , dumAtts)
         end if
         delete([/varAtts, dumAtts, vdims/])   ; may change size next iteration
       end if
     end if
   end do       ; nv

;===================================================================
; explicitly exit file definition mode. **NOT REQUIRED in NCL**
;===================================================================
    setfileoption(dst_file,"DefineMode",False)

        
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;        Write the basic and regridded data values to the predefined structures
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    dst_file->lat = (/ lat /)
    dst_file->lon = (/ lon /)
    dst_file->gw  = (/ gw  /)

    do nv=0,nvar_in-1
       if (.not.any(var_in(nv).eq.var_in_exclude)) then 
        if(var_out(0).eq."All_Variables"          .or. \
           any(var_in(nv).eq.var_out)             .or. \
           any(var_in(nv).eq.var_in_copy)       ) then

           rank_in = dimsizes(getfilevardimsizes(src_file, var_in(nv)))

           if (regridvar(nv)) then
               print(nv+"/"+nvar_in+" interpolating:"+var_in(nv))
               dst_file->$var_in(nv)$ = (/ ESMF_regrid_with_weights(src_file->$var_in(nv)$,wgtFileName,Opt_RGRD) /)
           else
               print(nv+"/"+nvar_in+" copying:"+var_in(nv))
               dst_file->$var_in(nv)$ = (/ src_file->$var_in(nv)$ /)
           end if
        end if
       else
         print(nv+"/"+nvar_in+" skipping:"+var_in(nv))
       end if
    end do       ; nv

    delete(regridvar)
    delete(var_in)
    print("done with file")
end do ; loop over files
end
