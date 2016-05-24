
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
;  mergenc = systemfunc("cd movies; perl $HOMME_HOME/utils/mergenc.pl *.nc")
  fnames =  systemfunc("ls movies/baroclinic1.cube.nc")
;  fnames =  systemfunc("ls combined/*1_.nc")
  f      = addfiles(fnames,"r")       ; Open netCDF files.	

  ; change the name of the field to plot here
  ; current field names are Q, ps, ke, u, v, T, zeta - names are case sensitive
  ; prec, accum
  ; "wind" or "stream" to get a graph of barbs or streamlines only (without contour overlay)
  
  field = "T"
;  field = "zeta"
; JPE wind vectors not working yet.
;  wind = True
  wind = False
; scale the wind speed by this factor
  wscale = 100.
;  sample the winds by this factor
  wsample = 1

; plot streamlines over contours (choose this or wind)
;  stream = True
  stream = False

  docontour=True
  if(field.eq."wind") then
    wind=True  
    docontour=False
  end if
  if(field.eq."stream") then
    stream=True
    docontour=False
  end if
  if(stream .and. wind) then
    print("Cannot specify both streamlines and wind barbs "+field)
    exit
  end if
;
; Open a file to send output to. You can use "ps", "pdf", "x11",
; or "ncgm".
;
  output_format = "ncgm"


  lat1d=doubletofloat(f[:]->lat)*180.0/3.14

  lon1d=doubletofloat(f[:]->lon)*180.0/3.14

  


  fcnt = dimsizes(fnames)
  vcnt=new((/fcnt,5/),"integer")
  do i=0,fcnt-1
    print ("i=" + i + "f[i]=" + getfiledimsizes(f[i]) )
    vcnt(i,:) = getfiledimsizes(f[i])
  end do	   	 
  
  ntimes = vcnt(0,3)
  nlev = vcnt(0,1)
  nxy = dim_sum(vcnt(:,0))

; change the levels and or times to plot here
  firsttime=0
  lasttime=ntimes-1
;  lasttime=12
  timeinc=1


  firstlev=nlev-1
  lastlev=nlev-1
  levinc=2
  if(wind .or. stream) then	
    u = new((/ntimes,nlev,nxy/),"float")
    v = new((/ntimes,nlev,nxy/),"float")
  end if
  if(docontour) then
    var = new((/ntimes,nlev,nxy/),"float")
  end if
  nxe=0
 print ("ntimes="+ntimes+" nlev="+nlev+" nxy="+nxy)

;  hyam=doubletofloat(f[0]->hyam(firstlev:lastlev:levinc))
;  ps0 = doubletofloat(hyam@ps0)
;  hybm=doubletofloat(f[0]->hybm(firstlev:lastlev:levinc))
  fieldis2d = False
  if(field .eq. "ps" .or. field .eq. "prec" .or. field .eq. "accum") then
    fieldis2d = True
  end if
  if(field .eq. "udrag" .or. field .eq. "vdrag" .or. field .eq. "tsflx" .or. field .eq. "qsflx") then
    fieldis2d = True
  end if
  do n=firsttime,lasttime,timeinc
    if(fieldis2d) then
      varmin = 1.e9
      varmax = -1.
       var(n,0,:)= doubletofloat(f[:]->$field$(n,:))
       firstlev=0
       lastlev=0
;       do k=0,nxy-1
;         if(var(n,0,k) .lt. varmin) then
;           varmin = var(n,0,k)
;           kvarmin = k
;         end if
;         if(var(n,0,k) .gt. varmax) then
;           varmax = var(n,0,k)
;           kvarmax = k
;         end if
;        end do
;        print ("varmin("+kvarmin+")="+varmin+" lat,lon "+lat1d(kvarmin)+", "+lon1d(kvarmin))
;        print ("varmax("+kvarmax+")="+varmax+" lat,lon "+lat1d(kvarmax)+", "+lon1d(kvarmax))
    else
      if(docontour) then
        
        do k=firstlev,lastlev,levinc
	  var(n,k,:)=doubletofloat(f[:]->$field$(n,k,:))
        end do
      end if
    end if
    if(wind.or.stream) then
      do k=firstlev,lastlev,levinc
	u(n,k,:)=wscale*doubletofloat(f[:]->u(n,k,:))
	v(n,k,:)=wscale*doubletofloat(f[:]->v(n,k,:))
      end do
    end if


  end do


  if(stream) then
    nx = 90
    ny = 45
    plon = 360.*fspan(0.,1.,nx)
    plat =  90.*fspan(-1.,1.,ny)
  end if
;  print(min(plat)+" "+max(plat)+" "+min(plon)+" "+max(plon))
 
  time = f[0]->time
;
; It will also look for a "homme.res" file, where you can set 
; resources that you want set all the time, like the font, or
; the background/foreground colors.
;
 
  wks = gsn_open_wks(output_format,field)
  print ("Writing file "+field+"."+output_format)
;
; The first call below selects a colormap, and the second reverses it.
; The colors go from violet at the low valus to blue at the high
; values.
;
; To see the various colormaps, go to:
;
; http://ngwww.ucar.edu/ncl/coltable.html
;
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;  gsn_reverse_colormap(wks)

;
; Set some plot resources.  Order doesn't matter here.
;
; Any time you don't want a plot drawn, just set the resources:
;
;  res@gsnDraw  = False
;  res@gsnFrame = False
;
; Of course, you'll have to remember to delete these resource
; settings or change them to True if you want later plots to
; be drawn.
;
  res                     = True
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
;  res@gsnSpreadColorStart = 24           ; at color index 24.

  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = True        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels

  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = False        ; Turn of labelbar box lines.

  if(docontour) then  
    res@sfXArray            = lon1d        ; Required to tell NCL where to
    res@sfYArray            = lat1d        ; overlay data on globe.

    res@cnLevelSelectionMode = "ManualLevels"

;    res@cnMinLevelValF = 0.0025
;    res@cnMaxLevelValF = 0.0048 
;    res@cnLevelSpacingF = 0.25
  
    res@mpOutlineOn          = False
    res@mpFillOn             = False
    res@mpCenterLonF         = -90
    res@mpShapeMode          = "FreeAspect"
  end if
  if(wind) then
    res@gsnFrame=False     
    wmsetp("wdf", 1); barbs point toward wind direction
  end if
  if(stream) then
    res2=True
    res2@gsnMaximize         = True 
    res2@vfXArray = plon
    res2@vfYArray = plat



    if(docontour) then
      res2@gsnFrame=False     
      res2@gsnDraw = False
      res@gsnDraw = False
      res@gsnFrame=False     
    else
      res2@mpOutlineOn          = False
      res2@mpFillOn             = False
      res2@mpCenterLonF         = 90.0
      res2@mpShapeMode          = "FreeAspect"
    end if
  end if

    res@mpMinLatF = 20
    res@mpMaxLatF = 80
    ;res@vpWidthF            = 1           ; change aspect ratio of plot
    res@vpHeightF           = .3

  do k= firstlev,lastlev,levinc

    res@cnMinLevelValF = 250
    res@cnMaxLevelValF = 295
    res@cnLevelSpacingF = 2.5

    print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 
    do n = firsttime,lasttime,timeinc
;    res@cnMinLevelValF = min(var(n,k,:))
;    res@cnMaxLevelValF = max(var(n,k,:))
;    res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/32.
;    print(res@cnMinLevelValF+" "+res@cnMaxLevelValF) 
 
     if(docontour) then
        res@tiMainString        = field+" at level = "+k+" time="+time(n)+" days"    ;
        contour = gsn_csm_contour_map(wks,var(n,k,:),res)     ; Draw a contour plot.
        print("Contour complete: "+res@tiMainString)
      end if
      if(wind) then
        wmbarbmap(wks,lat1d(::wsample), lon1d(::wsample), u(n,k,::wsample), v(n,k,::wsample))
	frame(wks)
      end if
      if (stream) then
       	us = triple2grid(lon1d(::wsample),lat1d(::wsample),u(n,k,::wsample),plon,plat,False)
        vs = triple2grid(lon1d(::wsample),lat1d(::wsample),v(n,k,::wsample),plon,plat,False)
        if(docontour) then
          strm = gsn_streamline(wks,us,vs,res2)
	  overlay(contour,strm)
	  draw(contour)
	  frame(wks)
        else
          strm = gsn_streamline_map(wks,us,vs,res2)

        end if
      end if
    end do
  end do





end

