load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
;  read spectra.nc  (computed from history output with ke.ncl)
;  plot from spectra.nc with ke2.ncl
;
;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; smooth the spectra for plottign
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function spec_smooth(spc[*][*]:double,nlat:integer,kmax:integer)
local n1
;
begin
  ; SMOOTH wave numbers > 15
  spc_smooth=spc
  spc_smooth2=spc
  eps=.25

  do k=0,kmax
    do n1=1,nlat-2
      spc_smooth2(:,n1) = eps*spc_smooth(:,n1-1)+(1-2*eps)*spc_smooth(:,n1)+eps*spc_smooth(:,n1+1)
    end do
    spc_smooth = spc_smooth2
  end do
  return(spc_smooth)
end






begin

;fnames =  systemfunc("ls held_suarez0-nhu1.6e12-*.nc spectra_ape240.nc")
;fnames =  systemfunc("ls held_suarez0-nhu1.6e12-*.nc")
fnames =  systemfunc("ls held_suarez0-nhu*-750.nc")
avespc   = addfiles(fnames,"r")       ; Open netCDF files.	
print("fnames = "+fnames)
;output_format = "ncgm"
output_format = "pdf"
wks = gsn_open_wks(output_format,"spectra2")

; for each file opened above, give the name to use in the legend
;
legends=new( (/10/), "string")

;xi=(/0/)
;legends=(/":F21:0.125:F34:0  :F33:n:F21: = 2e12"/)

;xi=(/0,1,2/)
;  subtitle = "HOMME HS :F33:n:F21: = 1.6e12"
;legends(0)=(/":F21:0.125:F34:0  :F21:250mb "/)
;legends(1)=(/":F21:0.125:F34:0  :F21:500mb "/)
;legends(2)=(/":F21:0.125:F34:0  :F21:750mb "/)

xi=(/1,0,2/)
  subtitle = "HS 750mb"
legends(0)=(/":F21:0.125:F34:0  :F33:n:F21: = 5.0e11"/)
legends(1)=(/":F21:0.125:F34:0  :F33:n:F21: = 1.6e12"/)
legends(2)=(/":F21:0.125:F34:0  :F33:n:F21: = 5.0e12"/)



;xi=(/3,0/)
;  subtitle = ":F21:CAM3.1/HOMME APE  250mb"
;legends(0)=(/":F21:0.125:F34:0  :F21:APE"/)
;legends(1)=(/":F21:0.125:F34:0  :F21:HS"/)


fcnt = dimsizes(xi)
nlat=4096;  %dimsizes(avespc[0]->spc_avg);

deg = new( (/fcnt/), "double" )
spc =   new ( (/fcnt,nlat/), "double")
spc = spc@_FillValue
spc_vr=spc
spc_dv=spc
spc_conv=spc
spc_diff=spc
spc_f=spc
x   = new ( (/fcnt,nlat/), "double")
x = x@_FillValue

do i=0,fcnt-1
   fk=xi(i)
   nlat=dimsizes(avespc[fk]->spc_avg);
   if (isfilevar(avespc[fk],"spc_avg_conv")) then
   spc_conv(i,0:nlat-1)=avespc[fk]->spc_avg_conv
   spc_diff(i,0:nlat-1)=avespc[fk]->spc_avg_diff
   spc_f(i,0:nlat-1)=avespc[fk]->spc_avg_f
   end if

   spc(i,0:nlat-1)=avespc[fk]->spc_avg
   spc_vr(i,0:nlat-1)=avespc[fk]->spc_avg_vr
   spc_dv(i,0:nlat-1)=avespc[fk]->spc_avg_dv
   x(i,0:nlat-1) = avespc[fk]->x_avg
   deg(i) = avespc[fk]->deg
   fnames(i) = chartostring(avespc[fk]->fnames(0,:) )
end do




  do fk=0,fcnt-1

;    tmax=340
;     tmax = doubletoint(180/deg(fk)) ; linear grid
    tmax = doubletoint(120/deg(fk)) ; 2/3 dealiased grid
	

    ; dont plot wave numbers above tmax
    do i=tmax,nlat-1
       spc(fk,i)=spc@_FillValue
       spc_diff(fk,i)=spc@_FillValue
       spc_conv(fk,i)=spc@_FillValue
       spc_dv(fk,i)=spc@_FillValue
       spc_vr(fk,i)=spc@_FillValue
    end do         
    hdeg=sprintf("%4.2f", deg(fk))
    print(fnames(fk)+" degree="+hdeg)
  end do  


  athick=2
  lthick=2





  res=True
  ;res@gsLineThicknessF = 3.0
  res@xyXStyle ="Log"
  res@xyYStyle ="Log"
  res@tmXBMinorPerMajor = 8   
  res@tmXTMinorPerMajor = 8   
  res@tmYRMinorPerMajor = 8   
  res@tmYLMinorPerMajor = 8   
  res@tmXMajorGridThicknessF = athick
  res@tmYLLabelFont = 21    ;  21 = helvetica
  res@tmXBLabelFont = 21    ;  22 = helvetica-bold 

  res@tmXBMajorThicknessF = athick
  res@tmYLMajorThicknessF = athick
  res@tmXBMinorThicknessF = athick 
  res@tmYLMinorThicknessF = athick 
  res@tmBorderThicknessF = athick
  res@gsnFrame = False
  res@xyLineThicknessF = lthick

  colors = (/"white","black","OrangeRed4","Red3",\
"red1","Red3","Red4","Blue","lightblue","green",\
"OrangeRed1","OrangeRed2","OrangeRed3","OrangeRed4","DarkOrange",\
"red1","red2","red3","red4",\
"DarkTurquoise","Turquoise","Cyan","blue",\
"gold","greenyellow",\
"GoldenRod","DarkGoldenRod1","lightgoldenrod1","lightyellow","lightgoldenrodyellow","khaki1"\
/)
  gsn_define_colormap(wks, colors)

  res@xyLineColors = (/"blue","green","DarkGoldenRod1","red3","red3"/)
  res@xyDashPattern = 0





  resL=res
  resL@trXMinF = 1.       
  resL@trXMaxF = 1000.
  resL@trYMinF = -10e-2
  resL@trYMaxF =  10e-2
  resL@xyXStyle ="Log"
  resL@xyYStyle ="Linear"
  resL@tiXAxisString          = ":F21:spherical wavenumber"
  resL@tiYAxisString          = ":F21:KE/s"

  fk=0
  spc_diff = spc_diff*(100000/26.)/9.8
  spc_conv = spc_conv*(100000/26.)/9.8
  spc_f = spc_f*(100000/26.)/9.8
  print("sum spc_diff = "+sum(spc_diff(fk,:)))
  print("sum spc_conv = "+sum(spc_conv(fk,:)))
  print("sum spc_f = "+sum(spc_f(fk,:)))

  resL@tiMainString = "(IE->KE)/E"
  resL@tiYAxisString          = ":F21:1/s"
  resL@trYMinF = -.002
  resL@trYMaxF =  .002

  spc_conv = (/ spec_smooth(spc_conv,tmax,50) /)
  spc_diff = (/ spec_smooth(spc_diff,tmax,50) /)

;  spc4 = spc
;  do i=0,tmax-1
;    spc4(:,i) = -.08*spc(:,i)*((i+1)^4.0)/(tmax^4.0)
;  end do
;  print("sum spc4 = "+sum(spc4(fk,:)))
  
;  plot  = gsn_xy (wks,x,spc_diff/spc,resL) ; create plot
  plot  = gsn_xy (wks,x,spc_conv/spc,resL) ; create plot
;  plot  = gsn_xy (wks,x,spc4/spc,resL) ; create plot
;  plot  = gsn_xy (wks,x,spc_f/spc,resL) ; create plot
  frame(wks)                                     ; now advance



  res@tiMainString = "|F|, |C|, |D|"
  res@trYMinF =  1e-10
  res@trYMaxF =  1e-1

  res@tiYAxisString          = ":F21:KE/s"
  res@tiXAxisString          = ":F21:spherical wavenumber"

;  spc_f = (/ spec_smooth(spc_f,tmax,10) /)

;   spc_conv = (/ spec_smooth(spc_conv,tmax,40) /)
;   spc_conv = where(spc_conv.gt.0, spc_conv, spc_conv@_FillValue)

;   spc_diff = (/ spec_smooth(spc_diff,tmax,10) /)
;   spc_diff = where(spc_diff.lt.0, spc_diff, spc_diff@_FillValue)

;  plot  = gsn_xy (wks,x,abs(spc_f),res) ; create plot
  plot  = gsn_xy (wks,x,abs(spc_conv),res) ; create plot
  plot  = gsn_xy (wks,x,abs(spc_diff),res) ; create plot
  frame(wks)                                     ; now advance


  res@tiMainString           = ":F21:CLIMATOLOGY"
  res@tiXAxisString          = ":F21:spherical wavenumber"
  res@tiYAxisString          = ":F21:Kinetic Energy"

  res@trXMinF = 1       
  res@trXMaxF = 1e3
  res@trYMinF = 1e-6      
  res@trYMaxF = 1e4   


;spc=spc/2
  plot  = gsn_xy (wks,x,spc,res) ; create plot
;  res@xyDashPattern = 2
;  plot  = gsn_xy (wks,x,spc_dv,res) ; create plot
; VORTICITY COMPONENT
;  plot  = gsn_xy (wks,x,spc_vr,res) ; create plot
  res@xyDashPattern = 0

  xref = (/10.,300./)
  spc_ref = 100*(xref/10)^(-3) ; make it 6 @ k=15
  polyres                  = True
  polyres@gsLineThicknessF = 2.0
  polyres@gsLineColor = "black"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline

  ypos=6000.
  do fk=0,fcnt-1
    ypos = ypos/3
    xx = (/15,25/)
    yy = (/ypos,ypos/)
    polyres@gsLineThicknessF = lthick
    polyres@gsLineColor = res@xyLineColors(fk)
    gsn_polyline(wks,plot,xx,yy,polyres)              ; add polyline
    polytext = True
    polytext@txFontHeightF      = 0.020                 ; change text size
    polytext@txJust             = "CenterLeft"          ; text justification
    gsn_text(wks,plot,legends(fk),30,ypos,polytext); add text
    print(fk+" legend = "+legends(fk))
  end do 
  polytext@txJust             = "BottomCenter"          ; text justification
  gsn_text_ndc(wks,subtitle,.5,.805,polytext); add text

  frame(wks)                                     ; now advance

  if (False) then  
    spc2=spc
    do n1=0,tmax-1
      spc2(:,n1) = spc(:,n1)*(n1^3)
    end do
    res@tiYAxisString          = ":F21:Compensated Kinetic Energy"
    res@trYMinF = 1e3
;  res@trYMinF = 1e1
    res@trYMaxF = 2e5   
    plot  = gsn_xy (wks,x,spc2,res) ; create plot
    
    xref = (/10.,300./)
    spc_ref = (xref^3)*100*(xref/10)^(-3) ; make it 6 @ k=15
    polyres@gsLineThicknessF = 2.0
    polyres@gsLineColor = "black"
    gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline


    frame(wks)                                     ; now advance 
  end if

  spc2=spc
  do n1=0,tmax-1
    spc2(:,n1) = spc(:,n1)*(n1^(5./3.))
  end do
  res@tiYAxisString          = ":F21:Compensated Kinetic Energy"
  res@trYMinF = 1e0
  res@trYMaxF = 1e4
  plot  = gsn_xy (wks,x,spc2,res) ; create plot

  res@xyDashPattern = 2
  do n1=0,tmax-1
    spc2(:,n1) = spc_vr(:,n1)*(n1^(5./3.))
  end do
; VORTICITY COMPONENT
;  plot  = gsn_xy (wks,x,spc2,res) ; create plot

  do n1=0,tmax-1
    spc2(:,n1) = spc_dv(:,n1)*(n1^(5./3.))
  end do
; DIV COMPONENT
;  plot  = gsn_xy (wks,x,spc2,res) ; create plot
  res@xyDashPattern = 0

  ypos=10000.
  do fk=0,fcnt-1
    ypos = ypos/1.6
    xx = (/15,25/)
    yy = (/ypos,ypos/)
    polyres@gsLineThicknessF = lthick
    polyres@gsLineColor = res@xyLineColors(fk)
    gsn_polyline(wks,plot,xx,yy,polyres)              ; add polyline
    polytext = True
    polytext@txFontHeightF      = 0.020                 ; change text size
    polytext@txJust             = "CenterLeft"          ; text justification
    gsn_text(wks,plot,legends(fk),30,ypos,polytext); add text
    print(fk+" legend = "+legends(fk))
  end do 
  gsn_text(wks,plot,":F21:-3",25,800,polytext); add text
  gsn_text(wks,plot,":F21:-5/3",150,200,polytext); add text

  polytext@txJust             = "BottomCenter"          ; text justification
  gsn_text_ndc(wks,subtitle,.5,.805,polytext); add text


  xref = (/10.,70./)
  spc_ref = (xref^(5./3.))*50*(xref/10)^(-3) ; make it 6 @ k=15
  polyres@gsLineThicknessF = 2.0
  polyres@gsLineColor = "black"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline

  xref = (/90.,500./)
  spc_ref = (xref^(5./3.))*3*(xref/10)^(-5./3.) ; make it 6 @ k=15
  polyres@gsLineThicknessF = 2.0
  polyres@gsLineColor = "black"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
  frame(wks)                                     ; now advance 


end
