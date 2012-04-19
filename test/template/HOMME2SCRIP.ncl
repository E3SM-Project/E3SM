load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"

;
;  read a HOMME grid template file and convert to SCRIP format
;  Mark Taylor 2008
; 

begin

template = name + "_tmp.nc"
output = name + "_scrip.nc"
system("rm -f "+output)
print("input: "+template+" output:"+output)
 

fin   = addfile(template,"r")   ;
fout  = addfile(output,"c")   ;

grid_size = dimsizes(fin->lat)
print("grid_size = "+grid_size)
grid_rank = 1 

; see if there are any hexagons:  
grid_corners = 0
kmax = dimsizes(fin->cv_lon)
do k=1,kmax(0)-1
  maxlon = max(abs( fin->cv_lon(k,:) ))
  if ( maxlon .lt. 0.000000001  .and. grid_corners.eq.0)  then
     grid_corners = k
  end if
end do
if ( grid_corners .eq. 0 )  then
   print (" Failed to find max number of nodes per element.")
else
   if ( grid_corners .eq. 4 )  then
      print("all control volumes are at most squares")
   else
      if (grid_corners .eq. 5 ) then
         print("detected pentagon control volumes")
      else
         if (grid_corners .eq. 6 ) then
            print("detected hexagon control volumes ")
         else  
            print("detected unstructured control volumes. Max Corners: "+grid_corners)
         end if
      end if
   end if
end if


print("lon min/max = "+min(fin->lon)+" "+max(fin->lon))
print("lat min/max = "+min(fin->lat)+" "+max(fin->lat))

grid_area = (/fin->area(:)/)
grid_area!0 = "grid_size"
grid_area@long_name = "area weights"
grid_area@units = "radians^2"


grid_center_lat = (/ fin->lat /)
grid_center_lat!0 = "grid_size"
grid_center_lat@units = "degrees"

grid_center_lon = (/ fin->lon /)
grid_center_lon!0 = "grid_size"
grid_center_lon@units = "degrees"

grid_corner_lon = new ( (/grid_size,grid_corners/), "double")
grid_corner_lat = new ( (/grid_size,grid_corners/), "double")
print("dims of cv_lon: "+dimsizes(fin->cv_lon))
print("dims of grid_corner_lat: " + dimsizes(grid_corner_lat))

do k=0,grid_corners-1
   grid_corner_lon(:,k) = (/fin->cv_lon(k,:)/)
   grid_corner_lat(:,k) = (/fin->cv_lat(k,:)/)
end do 
grid_corner_lon!1 = "grid_corners"
grid_corner_lon!0 = "grid_size"
grid_corner_lon@units = "degrees"

grid_corner_lat!1 = "grid_corners"
grid_corner_lat!0 = "grid_size"
grid_corner_lat@units = "degrees"

; reverse orientation at pole point
do i=0,grid_size-1
  if ( abs( 90 - abs(grid_center_lat(i)) ) .lt. 1e-9 ) then
    print("i="+i+" pole point: "+grid_center_lat(i)+" "+grid_center_lon(i))
    print("cv: "+grid_corner_lat(i,:)+" "+grid_corner_lon(i,:))
    ; 1 2 3 4  ->  1 4 3 2    swap pos 1,3
  end if
end do

grid_imask = new ( (/grid_size/), "double")
grid_imask = 1
grid_imask!0 = "grid_size"

grid_rank=1
grid_dims=1
grid_dims!0 = "grid_rank"

; write the 'title' global attribute
fout@title = "ne"+ne+"np4"
fout@Created_by = "HOMME2SCRIP.ncl mt 2009"
fout@history = "HOMME2SCRIP.ncl mt 2009"

; write all the data to the output field:
fout->grid_area = grid_area
fout->grid_center_lat = grid_center_lat
fout->grid_center_lon = grid_center_lon
fout->grid_corner_lon = grid_corner_lon
fout->grid_corner_lat = grid_corner_lat
fout->grid_imask = grid_imask
fout->grid_dims = grid_dims

end
