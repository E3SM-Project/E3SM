
load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"

begin

#
#  plot the cells in a HOMME grid template file:   *_latlon.nc 
#
#

rad2deg1 = 1.0
fnames = systemfunc("ls  ~/scratch1/preqx/template/arm_30_x8_lowconn.g_latlon.nc")


fg1     = addfile(fnames,"rw")       ; Open netCDF files.	
  field = "map"
  output_format = "ncgm"
;  output_format = "pdf"
  wks = gsn_open_wks(output_format,field)
  print ("Writing file "+field+"."+output_format)

  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;  gsn_reverse_colormap(wks)

;
  res                     = True

  
    res@mpOutlineOn          = True
    res@mpFillOn             = False
    res@mpShapeMode          = "FreeAspect"
    res@mpGridAndLimbOn        = False                ; no lat/lon lines
    res@mpPerimOn           = True
; increase viewport size (make map bigger)
; res@vpWidthF  = 0.90
; res@vpHeightF = 0.60
; res@vpXF      = 0.05
; res@vpYF      = 0.80
      
; res@gsnDraw  = False       ; don't draw the plots now
 res@gsnFrame = False       ; or advance the frame

 res@tiMainString        = "GLL control volumes"
;     contour = gsn_csm_contour_map(wks,fg1->area,res)     ; Draw a contour plot.

 contour = gsn_map(wks,"CylindricalEquidistant",res) 
 
; res@mpCenterLatF = 90
 res@mpCenterLatF = 45
 res@mpCenterLonF = -45
 res@gsnMaximize = True
; contour = gsn_map(wks,"Orthographic",res) 
 
 
 
; plot the grid points:
 resp=True
 resp@gsMarkerIndex = 1
;      resp@gsFillIndex = 0
 
 
 lon = fg1->lon(:)*rad2deg1
 lat = fg1->lat(:)*rad2deg1
 corners = fg1->element_corners
temp = dimsizes(corners)
ncorners = temp(0)
ncells = temp(1)



 print("ncells = "+ncells)
 print("lon min/max = "+min(lon)+" "+max(lon))
 print("lat min/max = "+min(lat)+" "+max(lat))


 irange1a = 0
; irange1a = n-100
 irange1b = ncells-1
if (irange1b .gt. 2000) then
   irange1b = 2000
end if

; size of grid points
rangesize = .001
;rangesize = .01
; size of CV verticies (0 to disable)
cvsize = .000
; size of special points in ilist array
listsize = .0003



       resp@gsMarkerSizeF = rangesize
       resp@gsMarkerColor = (/"green"/)
       resp@gsMarkerIndex = 1

    
       do i=irange1a,irange1b
   
         do j=0,ncorners-1
           j1 = j
           j2 = j+1
           if (j2 .gt. (ncorners-1) ) then
             j2 = 0
           end if
           j1_idx = corners(j1,i)-1
           j2_idx = corners(j2,i)-1
            
           xlon = (/ lon(j1_idx),lon(j2_idx) /)
           xlat = (/ lat(j1_idx),lat(j2_idx) /)
           if ( abs(xlon(1)-xlon(0)) .gt. 180.0) then
             if (xlon(1) .gt. xlon(0) ) then
               xlon(0)=xlon(0)+360.
             else
               xlon(1)=xlon(1)+360.
             end if
           end if
           resp@gsLineColor=(/"blue"/)
           gsn_polyline(wks,contour,xlon,xlat,resp)

         end do   
       end do



     frame(wks)
end

