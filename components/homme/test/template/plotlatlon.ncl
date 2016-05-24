
load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"

begin

;
;  plot the cells in a HOMME grid template file:   *_latlon.nc 
;  Mark Taylor 2008
;  used mostly for debuging the metadata
;

rad2deg1 = 1.0
;fnames = systemfunc("ls  ~/scratch1/preqx/template/arm_30_x8_lowconn.g_latlon.nc")
fnames = systemfunc("ls  ./ne8np4_latlon.nc")


  fg1     = addfile(fnames,"rw")       ; Open netCDF files.	
  field = "map"
  output_format = "ncgm"
;  output_format = "pdf"
  wks = gsn_open_wks(output_format,field)
  print ("Writing file "+field+"."+output_format)
  
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  
  res                     = True
  res@mpOutlineOn          = True
  res@mpFillOn             = False
  res@mpShapeMode          = "FreeAspect"
  res@mpGridAndLimbOn        = False                ; no lat/lon lines
  res@mpPerimOn           = True
  
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

 print("num subelement quads = "+ncells)
 print("lon min/max = "+min(lon)+" "+max(lon))
 print("lat min/max = "+min(lat)+" "+max(lat))
 
 
; size of CV verticies (0 to disable)
 cvsize = .000
; size of special points in ilist array
 listsize = .0003

 


; PLOT SUBELEM QUADS IN THIS RANGE:
 irange1a = 0
; irange1a = n-100
 irange1b = ncells-1
 if (irange1b .gt. 100) then
   irange1b = 100
 end if
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
; frame(wks)
   

 if (isfilevar(fg1,"GLLnodes")) then

   GLLnodes = fg1->GLLnodes
   dims = dimsizes(GLLnodes)
   np = dims(1)
   ncells = dims(0)
   ncorners=4
   print("number of spectral elemens = "+ncells)   

   lonloc=new( (/4/),"double")
   latloc=new( (/4/),"double")
   
; PLOT SUBELEM QUADS IN THIS RANGE:
   irange1a = 0
   irange1b = ncells-1
   if (irange1b .gt. 10) then
     irange1b = 10
   end if
   do ie=irange1a,irange1b
     
   do j=1,4
     if (j.eq.1) then
        j1_idx = GLLnodes(0,0,ie)
        j2_idx = GLLnodes(np-1,0,ie)
     end if     
     if (j.eq.2) then
        j1_idx = GLLnodes(np-1,0,ie)
        j2_idx = GLLnodes(np-1,np-1,ie)
     end if     
     if (j.eq.3) then
        j1_idx = GLLnodes(np-1,np-1,ie)
        j2_idx = GLLnodes(0,np-1,ie)
     end if     
     if (j.eq.4) then
        j1_idx = GLLnodes(0,np-1,ie)
        j2_idx = GLLnodes(0,0,ie)
     end if     

     ; 0 offset:
     j1_idx=j1_idx-1
     j2_idx=j2_idx-1

     xlon = (/ lon(j1_idx),lon(j2_idx) /)
     xlat = (/ lat(j1_idx),lat(j2_idx) /)
     if ( abs(xlon(1)-xlon(0)) .gt. 180.0) then
       if (xlon(1) .gt. xlon(0) ) then
         xlon(0)=xlon(0)+360.
       else
         xlon(1)=xlon(1)+360.
       end if
     end if
     resp@gsLineColor=(/"red"/)
     gsn_polyline(wks,contour,xlon,xlat,resp)
   end do   



   do i=0,np-1
     ii=0
     do j=0,np-1
       lonloc(ii) = lon( GLLnodes(i,j,ie)-1 )
       latloc(ii) = lat( GLLnodes(i,j,ie)-1 )
       ii=ii+1
     end do
     colors= (/"red","green","blue","yellow"/) 
     resp@gsMarkerSizeF = .005
     resp@gsMarkerColor = colors(i)
     resp@gsMarkerIndex = 1
     gsn_polymarker(wks,contour,lonloc,latloc,resp)  
   end do
   
   end do
   frame(wks)
 end if


end

