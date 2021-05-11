begin

  if ((.not.isdefined("in_file")).or.(.not.isdefined("out_file"))) then
    print((/"Run this script from the MALI_to_SCRIP.sh driver!"/))
    status_exit(1)
  end if

  fin = addfile(in_file, "r")

  ; Read in dims / variables, set attributes
  dims  = filevardimsizes(fin, "edgesOnCell") ; nCells x maxEdges
  size  = filevardimsizes(fin, "indexToCellID") ; nCells
  maxV = fin->nEdgesOnCell ; number of vertices on each cell
  edges = max(maxV)

  ; Constants:
  ; r needed to get area on unit sphere
  ; r2d needed to convert lats and lons from radians to degrees
  r = 6371220.d
  r2d = 90.d/acos(0.d)

  area           = (/fin->areaCell/(r*r)/)
  area@units     = "radians^2"
  area@long_name = "area weights"

  lats       = (/fin->latCell*r2d/)
  lats@units = "degrees"
  lons       = (/fin->lonCell*r2d/)
  lons@units = "degrees"

  clats = new((/size,edges/), "double")
  clons = new((/size,edges/), "double")
  clats@_FillValue = -9999.
  clons@_FillValue = -9999.
  do i=0,size-1
    jmax = maxV(i)-1
    clats(i,0:jmax) = (/fin->latVertex(fin->verticesOnCell(i,0:jmax)-1)*r2d/)
    clons(i,0:jmax) = (/fin->lonVertex(fin->verticesOnCell(i,0:jmax)-1)*r2d/)
    if (jmax.ne.edges-1) then
      clats(i,jmax+1:edges-1) = clats(i,jmax)
      clons(i,jmax+1:edges-1) = clons(i,jmax)
    end if
  end do
  clats@units = "degrees"
  clons@units = "degrees"

  imask = new(size, "double")
  imask@_FillValue = -9999.
  imask = 1.d
  if (min(fin->maxLevelCell).lt.1) then
    imask(ind(fin->maxLevelCell.eq.0)) = 0.d
  end if


  ; Outputs
  system("rm -f "+out_file)
  fout = addfile(out_file, "c")

  setfileoption(fout, "DefineMode", True)

  ; Global Attributes
  fAtt = True
  fAtt@title = "mpas120km"
  fAtt@Created_by = "mlevy NCL script"
  fAtt@Creation_date = systemfunc("date")
  fileattdef(fout, fAtt)

  ; Dimensions
  dimNames = (/"grid_size", "grid_corners", "grid_rank"/)
  dimSizes = (/size,        edges,          1/)
  dimUnlim = (/False,       False,          False/)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  ; Variables
  filevardef(fout, "grid_area",       "double" , "grid_size")
  filevardef(fout, "grid_center_lat", "double" , "grid_size")
  filevardef(fout, "grid_center_lon", "double" , "grid_size")
  filevardef(fout, "grid_corner_lat", "double" , (/"grid_size","grid_corners"/))
  filevardef(fout, "grid_corner_lon", "double" , (/"grid_size","grid_corners"/))
  filevardef(fout, "grid_imask",      "double" , "grid_size")
  filevardef(fout, "grid_dims",       "integer", "grid_rank")

  ; Variable Attributes
  filevarattdef(fout, "grid_area",       area)
  filevarattdef(fout, "grid_center_lat", lats)
  filevarattdef(fout, "grid_center_lon", lons)
  filevarattdef(fout, "grid_corner_lat", clats)
  filevarattdef(fout, "grid_corner_lon", clons)
  filevarattdef(fout, "grid_imask",      imask)

  ; Exit file definition mode, can now write data
  setfileoption(fout, "DefineMode", False)

  fout->grid_area       = (/area/)
  fout->grid_center_lat = (/lats/)
  fout->grid_center_lon = (/lons/)
  fout->grid_corner_lat = (/clats/)
  fout->grid_corner_lon = (/clons/)
  fout->grid_imask      = (/imask/)
  fout->grid_dims       = (/size/)

end


