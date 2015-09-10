load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"

;
;  read a HOMME grid template file and convert to SCRIP format
;  Mark Taylor 2008
; 

begin

if (.not. isvar("np")) then
   np=4
   print("NPTS not specified on command line. Using default 4.")
end if

template = name + "_tmp.nc"
output = name + "_latlon.nc"
system("rm -f "+output)
print("input: "+template+" output:"+output)
 

fin   = addfile(template,"r")   ;
fout  = addfile(output,"c")   ;

grid_size = dimsizes(fin->lat)
print("grid_size = "+grid_size)
grid_rank = 1 

print("lon min/max = "+min(fin->lon)+" "+max(fin->lon))
print("lat min/max = "+min(fin->lat)+" "+max(fin->lat))

;grid_area = (/fin->area(0,:)/)
;grid_area!0 = "grid_size"
;grid_area@long_name = "area weights"
;grid_area@units = "radians^2"

if (isfilevar(fin,"corners")) then
  dims=dimsizes(fin->corners)
  corners = new ( (/4,dims(1)/), "integer")
  do k=0,3
;   corners(:,k) = (/fin->corners(k,:)/)
    corners(k,:) = (/fin->corners(k,:)/)
  end do 
  corners!0 = "ncorners"
  corners!1 = "ncells"
end if


; write the 'title' global attribute
fout@title = "ne"+ne+"np"+np
fout@Created_by = "HOMME2META.ncl mt 2011"
fout@history = "HOMME2META.ncl mt 2011"

; write all the data to the output field:
fout->lat = fin->lat
fout->lon = fin->lon
fout->area = fin->area
if (isfilevar(fin,"corners")) then
  fout->element_corners = corners
end if

end
