load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"


begin

  ; Check for required input variables
  ; name of SCRIP grid file
  if (.not.isdefined("gridfile")) then
    print((/"ERROR: you must specify the variable gridfile when you run this script."/))
    status_exit(1)
  end if
  ; short name of SCRIP grid file (used to generate title / filename of plot)
  if (.not.isdefined("gridname")) then
    print((/"ERROR: you must specify the variable gridname when you run this script."/))
    status_exit(2)
  end if
  ; plot_ortho and plot_stereo are required, but set by default in gridplot.sh
  if (.not.isdefined("plot_ortho")) then
    plot_ortho=False
  end if
  if (.not.isdefined("plot_stereo")) then
    plot_stereo=False
  end if

  if (.not.isfilepresent(gridfile)) then
    print((/"Can not find "+gridfile+"!"/))
    status_exit(1)
  end if

  if (.not.(plot_ortho.or.plot_stereo)) then
    print((/"WARNING: no plots specified!"/))
    exit
  end if

  ; Check for optional input variables
  ; out_type (format of plot: PDF, PS, X11, PNG)
  ; default = pdf
  if (.not.isdefined("out_type")) then
    out_type = "pdf"
  end if
  ; center_lat -- latitude of center of orthographic projection plot
  ; default = 20.
  if (.not.isdefined("center_lat")) then
    center_lat = 20.
  end if
  ; center_lon -- longitude of center of orthographic projection plot
  ; default = -10.
  if (.not.isdefined("center_lon")) then
    center_lon = -10.
  end if
  ; out_dir -- directory where plot[s] will be saved
  ; default is ./
  if (.not.isdefined("out_dir")) then
    out_dir = "./"
  end if

  refine_type = "poles"
  num_smth = 0

  print((/"Plotting mesh from "+gridfile/))

  pi4= atan(1.d)
  pi2 = acos(0.d)
  pi = pi2*2.d
  f = addfile(gridfile,"r")
  lat = f->grid_corner_lat
  lon = f->grid_corner_lon

  if (lat@units.eq."radians") then
    lat = lat*180.d/pi
  end if
  if (lon@units.eq."radians") then
    lon = lon*180.d/pi
  end if

  dims = dimsizes(lon)
  nCells = dims(0)
  nVerts = dims(1)
  xlon = new ( (/nVerts+1/), "double")
  xlat = new ( (/nVerts+1/), "double")

  print("number of cells = "+nCells)
  print("number of verticies = "+nVerts)
  print("lat min/max = "+min(lat)+" "+max(lat))

  ; polygon resources
  res_p             = True

  res_p@gsLineThicknessF = 1.0
  res_p@gsLineColor   = "black"

  if (plot_ortho) then
    ; Orthographic Projection
    if ((out_type.eq."pdf").and.(plot_stereo)) then
      wks = gsn_open_wks(out_type,out_dir+"/"+gridname)
    else
      wks = gsn_open_wks(out_type,out_dir+"/"+gridname+"_ortho")
    end if
    gsn_define_colormap(wks,(/"White","Black","Tan","LightBlue","Blue"/))

    res = True
    res@tiMainString = gridname
    res@mpProjection      = "Orthographic"

    ; Center more towards Equator
    res@mpCenterLatF      =  center_lat
    res@mpCenterLonF      =  center_lon

    res@vpXF      = 0.05
    res@vpYF      = 0.9
    res@vpWidthF  = 0.9
    res@vpHeightF = 0.8

    res@gsnDraw  = False       ; don't draw the plots now
    res@gsnFrame = False       ; or advance the frame

    res@mpOutlineOn            = False
    res@mpPerimOn              = False
    res@mpLandFillColor        = "tan"
    res@mpOceanFillColor       = "LightBlue"
    res@mpInlandWaterFillColor = "Blue"
    res@mpGreatCircleLinesOn = True

    plot = gsn_csm_map(wks,res) ; create the plot
    draw(plot)

    do i=0,nCells-1
      if ( mod(i,500).eq.0) then
        print("i = "+(i+1)+"/"+nCells)
      end if
      xlon(0:(nVerts-1)) = lon(i,:)
      xlat(0:(nVerts-1)) = lat(i,:)
      xlon(nVerts)=xlon(0)
      xlat(nVerts)=xlat(0)

      do j=0,nVerts-1
        if ( abs(xlon(j+1)-xlon(j)) .gt. 180.0) then
          if (xlon(j+1) .gt. xlon(j) ) then
            xlon(j)=xlon(j)+360.
          else
            xlon(j+1)=xlon(j+1)+360.
          end if
        end if
      end do

      gsn_polyline(wks, plot, xlon,xlat,res_p)
    end do
    frame(wks)

  end if

  if (plot_ortho.and.plot_stereo) then
    if (out_type.ne."pdf") then
      delete(wks)
      wks = gsn_open_wks(out_type,out_dir+"/"+gridname+"_stereo")
      gsn_define_colormap(wks,(/"White","Black","Tan","LightBlue","Blue"/))
    end if
  end if

  ; Stereographic
  if (plot_stereo) then
    if (.not.plot_ortho) then
      wks = gsn_open_wks(out_type,out_dir+"/"+gridname+"_stereo")
      gsn_define_colormap(wks,(/"White","Black","Tan","LightBlue","Blue"/))
    end if
    res2 = True
    res2@tiMainString = gridname

    res2@vpXF      = 0.05
    res2@vpYF      = 0.9
    res2@vpWidthF  = 0.9
    res2@vpHeightF = 0.8

    res2@gsnDraw = False
    res2@gsnFrame = False

    res2@mpLandFillColor        = "tan"
    res2@mpOceanFillColor       = "LightBlue"
    res2@mpInlandWaterFillColor = "Blue"
    res2@mpGreatCircleLinesOn = True

    res2@mpGridAndLimbOn = True
    res2@mpGridLineDashPattern = 2
    res2@mpGridLineColor = "DarkSlateGray"

    res2@gsnMajorLonSpacing = 60
    res2@mpGridLonSpacingF = 60
    res2@gsnMajorLatSpacing = 45
    res2@mpGridLatSpacingF = 45

    plot2 = gsn_csm_map(wks,res2)
    draw(plot2)

    do i=0,nCells-1
      if ( mod(i,500).eq.0) then
        print("i = "+(i+1)+"/"+nCells)
      end if
      xlon(0:(nVerts-1)) = lon(i,:)
      xlat(0:(nVerts-1)) = lat(i,:)
      xlon(nVerts)=xlon(0)
      xlat(nVerts)=xlat(0)

      do j=0,nVerts-1
        if ( abs(xlon(j+1)-xlon(j)) .gt. 180.0) then
          if (xlon(j+1) .gt. xlon(j) ) then
            xlon(j)=xlon(j)+360.
          else
            xlon(j+1)=xlon(j+1)+360.
          end if
        end if
      end do

      gsn_polyline(wks, plot2, xlon,xlat,res_p)
    end do
    frame(wks)

  end if

end
