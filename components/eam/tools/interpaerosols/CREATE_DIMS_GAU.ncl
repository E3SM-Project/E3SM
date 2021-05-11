

begin

  RESLN = 

  NLAT = 
  NLON = 

  del_lon = new( (/1/) , double)
  gau_info = gaus(NLAT/2)
  lat     = gau_info(:,0)
  gwgt     = gau_info(:,1)
  del_lon(0)  = (360./NLON)
  lon     = ispan(0,NLON-1,1)*del_lon(0)

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
