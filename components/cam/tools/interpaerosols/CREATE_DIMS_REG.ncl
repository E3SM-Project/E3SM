

begin

  RESLN =

  NLAT =  
  NLON = 

  BEGLAT =
  BEGLON =

  ENDLAT =
  ENDLON =

  lat = new ( (/NLAT/), double )
  lon = new ( (/NLON/), double )

  lat = fspan(BEGLAT,ENDLAT,NLAT)
  lon = fspan(BEGLON,ENDLON,NLON)


  lat@long_name  = "latitude"
  lat@short_name = "lat"

  lon@long_name  = "longitude"
  lon@short_name = "lon"

;------------------------
; Write to file
;------------------------

  nfn = "regrid."+RESLN+".nc"
  system("rm "+nfn)
  cdf = addfile(nfn, "c")

  filedimdef(cdf, (/ "lat","lon"  /), \
                  (/ NLAT, NLON   /), \
                  (/ False, False /) )


  filevardef(cdf, "lat", "double", (/"lat"/))
  filevarattdef(cdf, "lat", lat)

  filevardef(cdf, "lon", "double", (/"lon"/))
  filevarattdef(cdf, "lon", lon)

  cdf->lat = lat
  cdf->lon = lon


end
