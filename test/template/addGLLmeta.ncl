load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"

;
;  Starting with a OMME grid template file with subcell metadata (the *_latlon.nc file),
;  compute the spectral element metadata, in array GLLnodes(np,np,nelem) and
;  add this data to the template file
;  
;  Mark Taylor 2014
; 

begin

if (.not. isvar("name")) then
   name = "ne8np4"
   print("No filename specifed on command line. Using default.")
end if
if (.not. isvar("np")) then
   np=4
   print("NPTS not specified on command line. Using default 4.")
end if
print("name = "+name)

output = name + "_latlon.nc"
print("input/output:"+output)
fin  = addfile(output,"rw")   ;


grid_size = dimsizes(fin->lat)
print("grid_size = "+grid_size)
grid_rank = 1 

print("lon min/max = "+min(fin->lon)+" "+max(fin->lon))
print("lat min/max = "+min(fin->lat)+" "+max(fin->lat))

if (.not. isfilevar(fin,"element_corners")) then
   print("GLL subelement metadata, 'element_corners', is missing from input file")
   exit
end if

nc=np-1
ncol=dimsizes(fin->lat)
dims=dimsizes(fin->element_corners)
nsubelems = dims(1) ; 1st dim is 4 (corners), 2nd dim is nsubelems
nelem = nsubelems/nc/nc

print("np="+np+" ncol="+ncol)
print("nelem=      "+nelem)

; check that ncols = 2 + nsubelems
; check that nsubelems is divisable by nc^2
if (mod(nsubelems,nc*nc) .ne. 0) then
   print("np value incorrect - does not divide nelem")
   exit
end if 


; determine number of spectral elements
; we know that:
;   ncols = 2 + nelem * ( npts-1)^2   =  2 + nsubelems
;
; assume npts=4:
;  nelem = nsubelems/(npts-1)^2
;
;      (1,4) (2,4) (3,4) (4,4)
;      (1,3) (2,3) (3,3) (4,3)
;      (1,2) (2,2) (3,2) (4,2)
;      (1,1) (2,1) (3,1) (4,1)

element_corners_local = new ( (/4,nc,nc/), "integer")
GLLnodes = new ( (/np,np,nelem/), "integer")

isubs=0
do ie = 0,nelem-1
   ; note that the order of these loops (i changes fastest) has to
   ; match the fortran code that created element_corners:
  ;        ii=0
  ;        do j=1,np-1
  ;           do i=1,np-1
  ;              ii=ii+1
  ;              subelement_corners(base+ii,1) = -fdofp_local(i,j,ie)
  ;              subelement_corners(base+ii,2) = -fdofp_local(i,j+1,ie)
  ;              subelement_corners(base+ii,3) = -fdofp_local(i+1,j+1,ie)
  ;              subelement_corners(base+ii,4) = -fdofp_local(i+1,j,ie)
  ;           end do
  ;        end do
  
   do j=0,nc-1
   do i=0,nc-1
      element_corners_local(:,i,j)=fin->element_corners(:,isubs)
      GLLnodes(i,j,ie) = element_corners_local(0,i,j)
      isubs=isubs+1
   end do
   end do

   i=nc
   do j=0,nc-1
      GLLnodes(i,j,ie) = element_corners_local(3,nc-1,j)
   end do

   j=nc
   do i=0,nc-1
      GLLnodes(i,j,ie) = element_corners_local(1,i,nc-1)
   end do

   GLLnodes(nc,nc,ie) = element_corners_local(2,nc-1,nc-1)
end do



GLLnodes!2 = "nelem"
GLLnodes!1 = "np"
GLLnodes!0 = "np"

; write the 'title' global attribute
;fout@title = "ne"+ne+"np4"
;fout@Created_by = "HOMME2META.ncl mt 2011"
;fout@history = "HOMME2META.ncl mt 2011"

; write all the data to the output field:
print("writing GLLnodes metadata to file")
fin->GLLnodes = GLLnodes

end
