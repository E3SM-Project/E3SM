load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


;run like this
;ncl make-convergence-remap.ncl fff=\"Th\" tt=8 myname=\"cpstarhy\" mynameShort=\"P\"

begin

;;; SET fili, tt, shortc, ztop, fff
;;; as file to read, time frame, shortc name for the plot suffix
;;; ztop is that same as in the run, fff as field


gg=9.80616
P0=100000.0
numz=200
ztop=15000.0
nlev = 128
ncols = 4608
shift=3

;-- define file name
  diri  = "./"

;myname="cpstarnh"
;mynameShort="P"

farea      =  addfile("area.nc","r")
area = farea->area
f00001 = myname+"0.0001-planar_rising_bubble1.nc"
fref     =  addfile(diri+f00001, "r")  
zref=fref->geo(tt,:,:)/gg
vals_ref=new((/numz,ncols/),double)
zvals=new((/numz/),double) ; newgrid
vals_current=new((/numz,ncols/),double)
fieldref=fref->$fff$(tt,:,:)

;make z array
dz=(ztop - 0.0)/(numz-1)
do ii=0,numz-1
zvals(ii) = (ii)*dz
end do

aarea = new((/numz,ncols/),double)
do i=0,numz-1
 aarea(i,:) = area(:)
end do

;; remap ref solutions
do ii=0,ncols-1
  zvalues=zref(::-1,ii)
  values=fieldref(::-1,ii)
  remapped_z = linint1(zvalues,values,False,zvals,0)
  vals_ref(:,ii) = remapped_z
end do

;now cut 
newref = vals_ref(shift:numz-shift,:)

;l2 from ref soln
lref=sqrt(sum(newref*newref*aarea(shift:numz-shift,:)))

ddt=(/"0.2","0.1","0.05","0.02","0.01","0.005","0.002","0.001"/)
dim_ddt=dimsizes(ddt)
nn=dim_ddt(0)

errarray = new((/nn/),double)

do i=0,nn-1
  iii=ddt(i)
  f    = myname+iii+"-planar_rising_bubble1.nc"
  ff   =  addfile(diri+f, "r")  
  ;fref     =  addfile(diri+f0001, "r")  

  field   =ff   ->$fff$(tt,:,:)
  zfield       =ff   ->geo(tt,:,:)/gg

  do ii=0,ncols-1
    zvalues=zfield(::-1,ii)
    values=field(::-1,ii)
    remapped_z = linint1(zvalues,values,False,zvals,0)
    vals_current(:,ii) = remapped_z
  end do

  err = vals_current(shift:numz-shift,:)-vals_ref(shift:numz-shift,:)
  err = err*err*aarea(shift:numz-shift,:)
  nerr = sqrt(sum(err))/lref

  errarray(i) = nerr

end do


print("------------------ L2 errors normalized by L2 ref soln")

print("series"+mynameShort+"=[ "+errarray(0)+","+errarray(1)+","+errarray(2)+","+errarray(3)+","+errarray(4)+","+errarray(5)+","+errarray(6)+","+errarray(7)+" ]")


;print("------------------ Now compare rates with dt")
;print("error ratio for dt=0.2 and dt=0.1 (should be 2) "+nerr02/nerr01)
;print("error ratio for dt=0.1 and dt=0.05 (should be 2) "+nerr01/nerr005)
;print("error ratio for dt=0.05 and dt=0.02 (should be 2.5)  "+nerr005/nerr002)
;print("error ratio for dt=0.02 and dt=0.01 (should be 2)  "+nerr002/nerr001)
;print("error ratio for dt=0.01 and dt=0.005 (should be 2)  "+nerr001/nerr0005)
;print("error ratio for dt=0.005 and dt=0.002 (should be 2.5)  "+nerr0005/nerr0002)
;print("error ratio for dt=0.002 and dt=0.001 (should be 2)  "+nerr0002/nerr0001)



end
