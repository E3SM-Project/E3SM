
load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"

begin


rad2deg1 = 1.0
;fnames = systemfunc("ls  ~/scratch1/preqx/template/ne4np4_subelem_scrip.nc")
fnames = systemfunc("ls  ~/scratch1/preqx/template/ARM_30_x8.g_scrip.nc")
;fnames = systemfunc("ls  ~/scratch1/preqx/template/ne16np4_scrip.nc")


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
 
 
 lon = fg1->grid_center_lon(:)*rad2deg1
 lat = fg1->grid_center_lat(:)*rad2deg1
 
 cvdims = dimsizes(fg1->grid_corner_lon)     
 ncorners = cvdims(1)
 
 cv_lon = fg1->grid_corner_lon(:,:)*rad2deg1
 cv_lat = fg1->grid_corner_lat(:,:)*rad2deg1
 n = dimsizes(lon)
 print("n = "+n)
 print("lon min/max = "+min(lon)+" "+max(lon))
 print("lat min/max = "+min(lat)+" "+max(lat))


 irange1a = 0
; irange1a = n-100
 irange1b = n-1
if (irange1b .gt. 50000) then
   irange1b = 500
end if
; irange1b = 2000
 

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
       gsn_polymarker(wks,contour,lon(irange1a:irange1b),lat(irange1a:irange1b),resp)  

    
       do i=irange1a,irange1b
   
         do j=0,ncorners-1
           j1 = j
           j2 = j+1
           if (j2 .gt. (ncorners-1) ) then
             j2 = 0
           end if
           xlon = (/ cv_lon(i,j1),cv_lon(i,j2) /)
           xlat = (/ cv_lat(i,j1),cv_lat(i,j2) /)
           if ( abs(xlon(1)-xlon(0)) .gt. 180.0) then
             if (xlon(1) .gt. xlon(0) ) then
               xlon(0)=xlon(0)+360.
             else
               xlon(1)=xlon(1)+360.
             end if
           end if
           resp@gsLineColor=(/"blue"/)
           gsn_polyline(wks,contour,xlon,xlat,resp)

           if (j.eq.0) then
             resp@gsMarkerColor = (/"red"/)
           end if
           if (j.eq.1) then
             resp@gsMarkerColor = (/"green"/)
           end if
           if (j.eq.2) then
             resp@gsMarkerColor = (/"yellow"/)
           end if
           if (j.eq.3) then
             resp@gsMarkerColor = (/"orange"/)
           end if
           if (j.eq.4) then
             resp@gsMarkerColor = (/"cyan"/)
           end if
           if (j.eq.5) then
             resp@gsMarkerColor = (/"brown"/)
           end if
           resp@gsMarkerSizeF = cvsize
           if (cvsize .gt. 0 ) then
             gsn_polymarker(wks,contour,cv_lon(i,j1),cv_lat(i,j1),resp)  
           end if
         end do   
       end do



     frame(wks)
end

