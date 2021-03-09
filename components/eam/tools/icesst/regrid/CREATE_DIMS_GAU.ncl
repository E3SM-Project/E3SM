

begin

  RESLN = "T85"

  NLAT = 128 
  NLON = 256

  gau_info = doubletofloat(gaus(NLAT/2))
  lat     = gau_info(:,0)
  gwgt     = gau_info(:,1)
  del_lon  = (360./NLON)
  lon     = ispan(0,NLON-1,1)*del_lon

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


  filevardef(cdf, "lat", "float", (/"lat"/))
  filevarattdef(cdf, "lat", lat)

  filevardef(cdf, "lon", "float", (/"lon"/))
  filevarattdef(cdf, "lon", lon)

  cdf->lat = lat
  cdf->lon = lon


end
