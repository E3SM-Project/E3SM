begin
;;
fil1="USGS-gtopo30_ne30np4pg2_16xdel2_forOroDrag.c20241029.nc"
fil2="USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
fil3="final-180.nc"
system("rm -r "+fil1)
system("cp -r "+fil2+" "+fil1)
system("ncks -A -v OA,OC,OL "+fil3+" "+fil1)
;;
end
