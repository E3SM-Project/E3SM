load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


;run like this
;ncl compute-ref-solns-difference-using-remap.ncl fff=\"Th\" tt=8

begin

;;; SET fili, tt, shortc, ztop, fff
;;; as file to read, time frame, shortc name for the plot suffix
;;; ztop is that same as in the run, fff as field


gg=9.80616
P0=100000.0
ztop=15000.0
numz=200
;fff="Th"


;-- define file name
  diri  = "./"


dt=0.0001

f00001   = "cpstarnh"+dt+"-planar_rising_bubble1.nc"
f00001cv = "cpstarhy"+dt+"-planar_rising_bubble1.nc"
;-- open file and read variable 

print("comparing with "+f00001)

frefcp     =  addfile(diri+f00001, "r")  
frefcv     =  addfile(diri+f00001cv, "r")  
farea      =  addfile("area.nc","r")

fieldrefcp=frefcp->$fff$(tt,:,:)
fieldrefcv=frefcv->$fff$(tt,:,:)
zcp=frefcp->geo(tt,:,:)/gg
zcv=frefcv->geo(tt,:,:)/gg
area = farea->area

;nlev = dimsizes(frefcp->lev)
;ncol = dimsizes(frefcp->ps)

nlev = 128
ncols = 4608
shift=3

vals_cp_z=new((/numz,ncols/),double)
vals_cv_z=new((/numz,ncols/),double)
zvals=new((/numz/),double)

;make z array
dz=(ztop - 0.0)/(numz-1)
do ii=0,numz-1
zvals(ii) = (ii)*dz
end do

p0=100000.0

;printVarSummary(area)

;since area is 1d filed, let's make it 2d
aarea = new((/numz,ncols/),double)
do i=0,numz-1
 aarea(i,:) = area(:)
end do

do ii=0,ncols-1
zvalues=zcp(::-1,ii)
values=fieldrefcp(::-1,ii)
remapped_z = linint1(zvalues,values,False,zvals,0)
vals_cp_z(:,ii) = remapped_z

zvalues=zcv(::-1,ii)
values=fieldrefcv(::-1,ii)
remapped_z2 = linint1(zvalues,values,False,zvals,0)
vals_cv_z(:,ii) = remapped_z2

;print(remapped_z - remapped_z2)

end do

newcp = vals_cp_z(shift:numz-shift,:)
newcv = vals_cv_z(shift:numz-shift,:)
aaarea = aarea(shift:numz-shift,:)


lref=sqrt(sum(newcv*newcv*aaarea))

err = newcv - newcp
err = err*err*aaarea
nerr = sqrt(sum(err))/lref
NOerr = sqrt(sum(err))
print("------------------ L2 error, not normalized and normalized by CPHY frame ")
print("dt=" +dt+ " NOerr "+NOerr +" nerr = "+nerr)





end
