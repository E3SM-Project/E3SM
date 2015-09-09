; ***********************************************
; xy_1.ncl
; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
; -------------------------------------------------------------
; read in scrip weight file --
; must contain regrid information from a CLM grid to a GLC grid
; -------------------------------------------------------------

; infileS    = addfile ("map_T31_to_GLC_Gland20km_bilin_da_100303.nc","r")
; infileS    = addfile ("map_T31_to_GLC_Gland10km_bilin_da_100304.nc","r")
; infileS    = addfile ("map_T31_to_GLC_Gland5km_bilin_da_100304.nc","r")
; infileS    = addfile ("map_fv1.9x2.5_to_GLC_Gland20km_bilin_da_100304.nc","r")
; infileS    = addfile ("map_fv1.9x2.5_to_GLC_Gland10km_bilin_da_100304.nc","r")
; infileS    = addfile ("map_fv1.9x2.5_to_GLC_Gland5km_bilin_da_100304.nc","r")
infileS    = addfile ("map_fv0.9x1.25_to_GLC_Gland20km_bilin_da_100304.nc","r")
; infileS    = addfile ("map_fv0.9x1.25_to_GLC_Gland10km_bilin_da_100304.nc","r")
; infileS    = addfile ("map_fv0.9x1.25_to_GLC_Gland5km_bilin_da_100304.nc","r")

; ---------------------------------------
; read in corresponding CLM fracdata file
; ---------------------------------------

; infileF    = addfile ("/fis01/cgd/cseg/csm/inputdata/lnd/clm2/griddata/fracdata_48x96_gx3v7_c090915.nc","r")
; infileF    = addfile ("/fis01/cgd/cseg/csm/inputdata/lnd/clm2/griddata/fracdata_1.9x2.5_gx1v6_c090206.nc","r")
infileF    = addfile ("/fis01/cgd/cseg/csm/inputdata/lnd/clm2/griddata/fracdata_0.9x1.25_gx1v6_c090317.nc","r")

; ------------------------------------------------
; read in necessary information from the two files
; ------------------------------------------------

src_grid_dims = infileS->src_grid_dims
nlon_map = src_grid_dims(0)
nlat_map = src_grid_dims(1)
pts_map  = infileS->col
dims1D   = dimsizes(pts_map)
nwgts    = dims1D(0)

longxy    = infileF->LONGXY
latixy    = infileF->LATIXY
landmask  = infileF->LANDMASK
dims2D    = dimsizes(landmask)
nlat_frac = dims2D(0)
nlon_frac = dims2D(1)

; ----------------------------------------------------------------
; first order check to make sure the two land grids are consistent
; ----------------------------------------------------------------

if (nlat_map .ne. nlat_frac) then
   print ("Number of latitudes not the same for both grids.  Stopping.")
   exit
end if

if (nlon_map .ne. nlon_frac) then
   print ("Number of longitudes not the same for both grids.  Stopping.")
   exit
end if

delete(infileS)
delete(infileF)

; ---------------------------------------------------
; calculate new 2-D array of GLCMASKs on the CLM grid
; ---------------------------------------------------

glcmask1D = new((/ nlon_map*nlat_map /),  "integer")
glcmask2D = new((/ nlat_frac, nlon_frac /), "integer")

do ij=0,nlat_map*nlon_map-1
   glcmask1D(ij) = 0
end do

do ij=0,nwgts-1
   npt = pts_map(ij)-1
   glcmask1D(npt) = 1
end do

do j=0,nlon_map-1
   do i=0,nlat_map-1
      ij = (i*nlon_map) + j
      glcmask2D(i,j) = 0
      if( (glcmask1D(ij) .eq. 1) .and. (landmask(i,j) .eq. 1)) then
         glcmask2D(i,j) = 1
      end if
   end do
end do

; ------------------------------------------
; write out GLCMASKs to new glcmaskdata file
; ------------------------------------------

; outfile = "glcmaskdata_48x96_Gland20km.nc"
; outfile = "glcmaskdata_48x96_Gland10km.nc"
; outfile = "glcmaskdata_48x96_Gland5km.nc"
; outfile = "glcmaskdata_1.9x2.5_Gland20km.nc"
; outfile = "glcmaskdata_1.9x2.5_Gland10km.nc"
; outfile = "glcmaskdata_1.9x2.5_Gland5km.nc"
outfile = "glcmaskdata_0.9x1.25_Gland20km.nc"
; outfile = "glcmaskdata_0.9x1.25_Gland10km.nc"
; outfile = "glcmaskdata_0.9x1.25_Gland5km.nc"
system("/bin/rm -f " + outfile)

fout1 = addfile(outfile,"c")

globalAtt             = True
; globalAtt@title       = "GLCMASK file for T31 to Glimmer Greenland 20 km Grid"
; globalAtt@title       = "GLCMASK file for T31 to Glimmer Greenland 10 km Grid"
; globalAtt@title       = "GLCMASK file for T31 to Glimmer Greenland 5 km Grid"
; globalAtt@title       = "GLCMASK file for fv1.9x2.5 to Glimmer Greenland 20 km Grid"
; globalAtt@title       = "GLCMASK file for fv1.9x2.5 to Glimmer Greenland 10 km Grid"
; globalAtt@title       = "GLCMASK file for fv1.9x2.5 to Glimmer Greenland 5 km Grid"
globalAtt@title       = "GLCMASK file for fv0.9x1.25 to Glimmer Greenland 20 km Grid"
; globalAtt@title       = "GLCMASK file for fv0.9x1.25 to Glimmer Greenland 10 km Grid"
; globalAtt@title       = "GLCMASK file for fv0.9x1.25 to Glimmer Greenland 5 km Grid"
globalAtt@history     = "SCRIP to GLCMASK conversion " + systemfunc("date")
fileattdef( fout1, globalAtt )

dimNames = (/"lsmlon", "lsmlat" /)  
dimSizes = (/ nlon_frac, nlat_frac /)
dimUnlim = (/ False , False /)
filedimdef(fout1, dimNames  , dimSizes,  dimUnlim )

filevardef   (fout1, "GLCMASK", "integer", (/ "lsmlat", "lsmlon" /))
GLCMASKAtt=0
GLCMASKAtt@longname = "land/glc mask"
GLCMASKAtt@units    = "0=no glc source and 1=glc source"
filevarattdef(fout1, "GLCMASK", GLCMASKAtt)

fout1->LONGXY    = longxy
fout1->LATIXY    = latixy
fout1->GLCMASK   = (/glcmask2D/)

end
