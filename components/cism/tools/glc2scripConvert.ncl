; ***********************************************
; glc2scripConvert.ncl
; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
;************************************************
; read in glimmer grid file
;************************************************

; infile    = addfile ("gland20km.nc","r")
; infile    = addfile ("gland10km.nc","r")
infile    = addfile ("gland5km.nc","r")

lat = infile->lat
lon = infile->lon
dims2D = dimsizes(lat)
nx     = dims2D(2)
ny     = dims2D(1)

delete(infile)

; -------------------------------------
; calculate 1-D arrays of lats and lons
; -------------------------------------

grid_size    = (nx-1)*(ny-1)
grid_corners = 4
grid_rank    = 2
grid_dims    = new((/grid_rank/), "integer")
grid_imask   = new((/grid_size/), "integer")
lat1D_center = new((/grid_size/),typeof(lat))
lon1D_center = new((/grid_size/),typeof(lon))
lat1D_corner = new((/grid_size, grid_corners/),typeof(lat))
lon1D_corner = new((/grid_size, grid_corners/),typeof(lon))

grid_dims(0) = nx-1
grid_dims(1) = ny-1

do iy=0,ny-2
  do ix=0,nx-2
    npt = (iy*(nx-1)) + ix
    grid_imask(npt)   = 1
    lat1D_center(npt) = 0.25*(lat(0,iy  ,ix  ) + lat(0,iy+1,ix  ) + \
                              lat(0,iy  ,ix+1) + lat(0,iy+1,ix+1))
    lon1D_center(npt) = 0.25*(lon(0,iy  ,ix  ) + lon(0,iy+1,ix  ) + \
                              lon(0,iy  ,ix+1) + lon(0,iy+1,ix+1))
    lat1D_corner(npt,0) = lat(0,iy  ,ix  )
    lat1D_corner(npt,1) = lat(0,iy  ,ix+1)
    lat1D_corner(npt,2) = lat(0,iy+1,ix+1)
    lat1D_corner(npt,3) = lat(0,iy+1,ix  )
    lon1D_corner(npt,0) = lon(0,iy  ,ix  )
    lon1D_corner(npt,1) = lon(0,iy  ,ix+1)
    lon1D_corner(npt,2) = lon(0,iy+1,ix+1)
    lon1D_corner(npt,3) = lon(0,iy+1,ix  )
  end do
end do

; ---------------------------------------
; write out lats and lons in SCRIP format
; ---------------------------------------

; outfile = "gland20km_scrip.nc"
; outfile = "gland10km_scrip.nc"
outfile = "gland5km_scrip.nc"
system("/bin/rm -f " + outfile)

fout1 = addfile(outfile,"c")

globalAtt             = True
; globalAtt@title       = "Glimmer Greenland 20 km Grid"
; globalAtt@title       = "Glimmer Greenland 10 km Grid"
globalAtt@title       = "Glimmer Greenland 5 km Grid"
globalAtt@history     = "GLC_to_SCRIP conversion " + systemfunc("date")
fileattdef( fout1, globalAtt )

dimNames = (/"grid_size", "grid_corners", "grid_rank" /)  
dimSizes = (/ grid_size ,  4,  2 /)
dimUnlim = (/ False , False, False /)
filedimdef(fout1, dimNames  , dimSizes,  dimUnlim )

filevardef   (fout1, "grid_dims", "integer", "grid_rank" )

filevardef   (fout1, "grid_center_lat", typeof(lat1D_center), "grid_size" )
grid_center_latAtt=0
grid_center_latAtt@units = "degrees"
filevarattdef(fout1, "grid_center_lat", grid_center_latAtt)

filevardef   (fout1, "grid_center_lon", typeof(lon1D_center), "grid_size")
grid_center_lonAtt=0
grid_center_lonAtt@units = "degrees"
filevarattdef(fout1, "grid_center_lon", grid_center_lonAtt)

filevardef   (fout1, "grid_imask", "integer", "grid_size")
grid_imaskAtt=0
grid_imaskAtt@units = "unitless"
filevarattdef(fout1, "grid_imask", grid_imaskAtt)

filevardef   (fout1, "grid_corner_lat", typeof(lat1D_corner), (/ "grid_size", "grid_corners" /))
grid_corner_latAtt=0
grid_corner_latAtt@units = "degrees"
filevarattdef(fout1, "grid_corner_lat", grid_corner_latAtt)

filevardef   (fout1, "grid_corner_lon", typeof(lon1D_corner), (/ "grid_size", "grid_corners" /))
grid_corner_lonAtt=0
grid_corner_lonAtt@units = "degrees"
filevarattdef(fout1, "grid_corner_lon", grid_corner_lonAtt)

fout1->grid_dims       = (/grid_dims/)
fout1->grid_center_lat = (/lat1D_center/)
fout1->grid_center_lon = (/lon1D_center/)
fout1->grid_imask      = (/grid_imask/)
fout1->grid_corner_lat = (/lat1D_corner/)
fout1->grid_corner_lon = (/lon1D_corner/)

end

