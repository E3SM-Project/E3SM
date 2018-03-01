load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
;  read spectra.nc  (computed from history output with ke.ncl)
;  plot from spectra.nc with ke2.ncl
;
;

begin

nedata = 30

if (nedata .eq. 30) then
  fnames =  systemfunc( \
  "ls ~/scratch1/preqx/jw/jw-ne30-4t90l26-hnu1.5e15-1/spectra.nc; " \
+  "ls ~/scratch1/preqx/jwnew/jw-ne30-4t90l26-hnu1e15-1/spectra.nc; " \
+  "ls ~/scratch1/preqx/jwnew/jw-ne30-4t90l26-hnu7e14-1/spectra.nc; " \
+  "ls ~/scratch1/preqx/jwnew/jw-ne30-4t90l26-hnu3e14-1/spectra.nc; " \
  )

  legends=(/\
  ":F21:1.0:F34:0  :F33:n:F21: = 1.5e15",\
  ":F21:1.0:F34:0  :F33:n:F21: = 1e15",\
  ":F21:1.0:F34:0  :F33:n:F21: = 7e14",\
  ":F21:1.0:F34:0  :F33:n:F21: = 3e14"/)
  
end if

if (nedata .eq. 60) then
  fnames =  systemfunc( \
;   "ls ~/scratch1/preqx/jw/jw-ne60-4t40l26-hnu4e14-2/spectra.nc; " \
  "ls ~/scratch1/preqx/jwnew/jw-ne60-4t40l26-hnu1.5e14-1/spectra.nc; " \
  + "ls ~/scratch1/preqx/jwnew/jw-ne60-4t40l26-hnu1e14-1/spectra.nc; " \
  + "ls ~/scratch1/preqx/jwnew/jw-ne60-4t40l26-hnu8e13-1/spectra.nc; " \
+  "ls ~/scratch1/preqx/jwnew/jw-ne60-4t40l26-hnu5e13-1/spectra.nc; " \
+  "ls ~/scratch1/preqx/jwnew/jw-ne60-4t40l26-hnu3e13-1/spectra.nc; " \
  )

  legends=(/\
;  ":F21:0.5:F34:0  :F33:n:F21: = 4e14",\
  ":F21:0.5:F34:0  :F33:n:F21: = 1.4e14",\
  ":F21:0.5:F34:0  :F33:n:F21: = 1e14",\
  ":F21:0.5:F34:0  :F33:n:F21: = 8e13",\
  ":F21:0.5:F34:0  :F33:n:F21: = 5e13",\
  ":F21:0.5:F34:0  :F33:n:F21: = 3e13"/)
  
end if


if (nedata .eq. 120) then
  fnames =  systemfunc( \
;   "ls ~/scratch1/preqx/jw/jw-ne120-4t20l26-hnu3e13-2/spectra.nc; " \
   "ls ~/scratch1/preqx/jwnew/jw-ne120-4t20l26-hnu1.5e13-1/spectra.nc; " \
  + "ls ~/scratch1/preqx/jwnew/jw-ne120-4t20l26-hnu1e13-1/spectra.nc; " \
  + "ls ~/scratch1/preqx/jwnew/jw-ne120-4t20l26-hnu8e12-1/spectra.nc; " \
+  "ls ~/scratch1/preqx/jwnew/jw-ne120-4t20l26-hnu5e12-1/spectra.nc; " \
  )
  
  legends=(/\
;  ":F21:0.25:F34:0  :F33:n:F21: = 3e13",\
  ":F21:0.25:F34:0  :F33:n:F21: = 1.5e13",\
  ":F21:0.25:F34:0  :F33:n:F21: = 1e13",\
  ":F21:0.25:F34:0  :F33:n:F21: = 8e12",\
  ":F21:0.25:F34:0  :F33:n:F21: = 5e12"\
/)
  
end if





print("fnames = "+fnames)
avespc = addfiles(fnames,"r")       ; Open netCDF files.	
fcnt = dimsizes(fnames)
;output_format = "ncgm"
output_format = "pdf"
wks = gsn_open_wks(output_format,"spectra")




ds=dimsizes(avespc[0]->spc);
nlat=ds(0);

spc = new ( (/fcnt,nlat/), "double", 1d-99 )
x   = new ( (/fcnt,nlat/), "double", 1d-99 )
do fk=0,fcnt-1
   spc(fk,:)=avespc[fk]->spc(:)
   x(fk,:) = avespc[fk]->x(:)
end do

  ; SMOOTH wave numbers > 15
  do k=0,1
  spc_smooth=0*spc
  do n1=1,nlat-2
    eps=.25
    spc_smooth(:,n1) = eps*spc(:,n1-1)+(1-2*eps)*spc(:,n1)+eps*spc(:,n1+1)
  end do
  end do



  lthick=4
  thick=4


  res=True
  res@xyXStyle ="Log"
  res@xyYStyle ="Log"
  res@tmXBMinorPerMajor = 8   
  res@tmXTMinorPerMajor = 8   
  res@tmYRMinorPerMajor = 8   
  res@tmYLMinorPerMajor = 8   
  res@tmXMajorGridThicknessF = lthick
  res@tmYLLabelFont = 21    ;  21 = helvetica
  res@tmXBLabelFont = 21    ;  22 = helvetica-bold 

  res@tmXBMajorThicknessF = thick
  res@tmYLMajorThicknessF = thick
  res@tmXBMinorThicknessF = thick 
  res@tmYLMinorThicknessF = thick 
  res@tmBorderThicknessF = thick
  res@tiMainString           = ":F21:DAY 30"
  res@tiXAxisString          = ":F21:spherical wavenumber"
  res@tiYAxisString          = ":F21:Kinetic Energy"
  res@gsnFrame = False
  res@xyLineThicknessF = lthick

  colors = (/"white","black","OrangeRed4","Red3",\
"red1","Red3","Red4","Blue","green",\
"orange","OrangeRed1","OrangeRed2","OrangeRed3","OrangeRed4","DarkOrange",\
"red1","red2","red3","red4",\
"DarkTurquoise","Turquoise","Cyan","blue",\
"gold","greenyellow",\
"GoldenRod","DarkGoldenRod1","lightgoldenrod1","lightyellow","lightgoldenrodyellow","khaki1"\
/)
  gsn_define_colormap(wks, colors)

  res@xyLineColors = (/"blue","green","DarkGoldenRod1","red3","turquoise","orange"/)
;  res@xyDashPatterns = (/0,1,0,1,0,1,0,1/)
  res@xyDashPattern = 0


  res@trXMinF = 10       
  res@trXMaxF = 1000
  res@trYMinF = 1e-8      
  res@trYMaxF = 1e3
; tighter settings:   
;  res@trXMinF = 1       
;  res@trXMaxF = 500
;  res@trYMinF = 1e-6      
;  res@trYMaxF = 1e3   



  plot  = gsn_xy (wks,x,spc_smooth,res) ; create plot

  xref = (/10.,300./)
  spc_ref = 40*(xref/10)^(-3) ; make it 6 @ k=15
  polyres                  = True
  polyres@gsLineThicknessF = 2.0
  polyres@gsLineColor = "black"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline

  ypos=300*2.5;
  do fk=0,fcnt-1
    ypos = ypos/2.5
    xx = (/25,45/)
    yy = (/ypos,ypos/)
    polyres@gsLineThicknessF = lthick
    polyres@gsLineColor = res@xyLineColors((fk))
    gsn_polyline(wks,plot,xx,yy,polyres)              ; add polyline
    polytext = True
    polytext@txFontHeightF      = 0.020                 ; change text size
    polytext@txJust             = "CenterLeft"          ; text justification
    gsn_text(wks,plot,legends(fk),50,ypos,polytext); add text
  end do 
  polytext@txJust             = "BottomCenter"          ; text justification
  gsn_text_ndc(wks,":F21:HOMME  850mb",.5,.805,polytext); add text


  frame(wks)                                     ; now advance


  do n1=0,nlat-1
    spc(:,n1) = spc(:,n1)*(n1^3)
  end do


  ; SMOOTH wave numbers > 15
  do k=0,1
  spc_smooth=0*spc
  do n1=1,nlat-2
    eps=.25
    spc_smooth(:,n1) = eps*spc(:,n1-1)+(1-2*eps)*spc(:,n1)+eps*spc(:,n1+1)
  end do
  end do



  res@tiYAxisString          = ":F21:Compensated Kinetic Energy"
  res@trYMinF = 9e3
  res@trYMaxF = 3e4   

;  res@xyLineThicknessF = 2
  plot  = gsn_xy (wks,x,spc_smooth,res) ; create plot
  xref = (/10.,300./)
  spc_ref = (xref^3)*25*(xref/10)^(-3) ; make it 6 @ k=15

  polyres@gsLineThicknessF = 2.0
  polyres@gsLineColor = "black"
  gsn_polyline(wks,plot,xref,spc_ref,polyres)              ; add polyline
  frame(wks)                                     ; now advance 


end
