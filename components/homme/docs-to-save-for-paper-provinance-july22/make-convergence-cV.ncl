load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


;run like this
;ncl make-convergence-cp.ncl fff=\"Th\" tt=5


begin

;;; SET fili, tt, shortc, ztop, fff
;;; as file to read, time frame, shortc name for the plot suffix
;;; ztop is that same as in the run, fff as field


gg=9.80616
P0=100000.0
;ztop=20000.0
numz=200
;fff="Th"


;-- define file name
  diri  = "./"

myname="cv"

  f02    = myname+"0.2-planar_rising_bubble1.nc"
  f01    = myname+"0.1-planar_rising_bubble1.nc"
  f005   = myname+"0.05-planar_rising_bubble1.nc"
  f002   = myname+"0.02-planar_rising_bubble1.nc"
  f001   = myname+"0.01-planar_rising_bubble1.nc"
  f0005  = myname+"0.005-planar_rising_bubble1.nc"
  f0002  = myname+"0.002-planar_rising_bubble1.nc"
  f0001  = myname+"0.001-planar_rising_bubble1.nc"
  f00001 = myname+"0.0001-planar_rising_bubble1.nc"

;-- open file and read variable 
  ff02     =  addfile(diri+f02, "r")  
  ff01     =  addfile(diri+f01, "r")  
  ff005    =  addfile(diri+f005, "r")  
  ff002    =  addfile(diri+f002, "r")  
  ff001    =  addfile(diri+f001, "r")  
  ff0005   =  addfile(diri+f0005, "r")  
  ff0002   =  addfile(diri+f0002, "r")  
  ff0001   =  addfile(diri+f0001, "r")  
  fref     =  addfile(diri+f00001, "r")  
  ;fref     =  addfile(diri+f0001, "r")  

;  var   = f->Th(0,:,{0},{-150000:147487})     ;-- first time step, latitude=40N, longitude=0-60E.
  lon_t = fref->lon            ;-- longitude=0-60E
  lat_t = fref->lat
  lev_t = fref->lev                    ;-- currently 17 levels

nlev = dimsizes(lev_t)
nlat = dimsizes(lat_t)
nlon = dimsizes(lon_t)

;print(lev_t)

field02   =ff02   ->$fff$(tt,:,:)
field01   =ff01   ->$fff$(tt,:,:)
field005  =ff005  ->$fff$(tt,:,:)
field002  =ff002  ->$fff$(tt,:,:)
field001  =ff001  ->$fff$(tt,:,:)
field0005 =ff0005 ->$fff$(tt,:,:)
field0002 =ff0002 ->$fff$(tt,:,:)
field0001 =ff0001 ->$fff$(tt,:,:)
fieldref=fref->$fff$(tt,:,:)


lref=sqrt(sum(fieldref*fieldref))

err02 = field02-fieldref
err02 = err02*err02
nerr02 = sqrt(sum(err02))/lref

err01 = field01-fieldref
err01 = err01*err01
nerr01 = sqrt(sum(err01))/lref

err005 = field005-fieldref
err005 = err005*err005
nerr005 = sqrt(sum(err005))/lref

err002 = field002-fieldref
err002 = err002*err002
nerr002 = sqrt(sum(err002))/lref

err001 = field001-fieldref
err001 = err001*err001
nerr001 = sqrt(sum(err001))/lref

err0005 = field0005-fieldref
err0005 = err0005*err0005
nerr0005 = sqrt(sum(err0005))/lref

err0002 = field0002-fieldref
err0002 = err0002*err0002
nerr0002 = sqrt(sum(err0002))/lref

err0001 = field0001-fieldref
err0001 = err0001*err0001
nerr0001 = sqrt(sum(err0001))/lref





print("------------------ L2 errors normalized by L2 ref soln")

print("nerr02 = "+nerr02)
print("nerr01 = "+nerr01)
print("nerr005 = "+nerr005)
print("nerr002 = "+nerr002)
print("nerr001 = "+nerr001)
print("nerr0005 = "+nerr0005)
print("nerr0002 = "+nerr0002)
print("nerr0001 = "+nerr0001)

print("seriesV=[ "+nerr02+","+nerr01+","+nerr005+","+nerr002+","+nerr001+","+nerr0005+","+nerr0002+","+nerr0001+" ]")


print("------------------ Now compare rates with dt")

print("error ratio for dt=0.2 and dt=0.1 (should be 2) "+nerr02/nerr01)
print("error ratio for dt=0.1 and dt=0.05 (should be 2) "+nerr01/nerr005)
print("error ratio for dt=0.05 and dt=0.02 (should be 2.5)  "+nerr005/nerr002)
print("error ratio for dt=0.02 and dt=0.01 (should be 2)  "+nerr002/nerr001)

print("error ratio for dt=0.01 and dt=0.005 (should be 2)  "+nerr001/nerr0005)
print("error ratio for dt=0.005 and dt=0.002 (should be 2.5)  "+nerr0005/nerr0002)
print("error ratio for dt=0.002 and dt=0.001 (should be 2)  "+nerr0002/nerr0001)



end
