load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


;run like this
; ncl compute-ref-solns-difference.ncl fff=\"Th\" tt=5

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



dt=0.0001

f00001   = "cpstarnh"+dt+"-planar_rising_bubble1.nc"
f00001cv = "cv"+dt+"-planar_rising_bubble1.nc"
;-- open file and read variable 
fref       =  addfile(diri+f00001, "r")  
frefcv     =  addfile(diri+f00001cv, "r")  

fieldref=fref->$fff$(tt,:,:)
fieldrefcv=frefcv->$fff$(tt,:,:)
lref=sqrt(sum(fieldref*fieldref))

err = fieldrefcv-fieldref
err = err*err
nerr = sqrt(sum(err))/lref
NOerr = sqrt(sum(err))
print("------------------ L2 error, not normalized and normalized by P frame ")
print("dt=" +dt+ " NOerr "+NOerr +" nerr = "+nerr)

dt=0.001

f00001   = "cpstarnh"+dt+"-planar_rising_bubble1.nc"
f00001cv = "cv"+dt+"-planar_rising_bubble1.nc"
;-- open file and read variable
fref       =  addfile(diri+f00001, "r")
frefcv     =  addfile(diri+f00001cv, "r")

fieldref=fref->$fff$(tt,:,:)
fieldrefcv=frefcv->$fff$(tt,:,:)
lref=sqrt(sum(fieldref*fieldref))

err = fieldrefcv-fieldref
err = err*err
nerr = sqrt(sum(err))/lref
NOerr = sqrt(sum(err))
print("------------------ L2 error, not normalized and normalized by P frame ")
print("dt=" +dt+ " NOerr "+NOerr +" nerr = "+nerr)

dt=0.01

f00001   = "cpstarnh"+dt+"-planar_rising_bubble1.nc"
f00001cv = "cv"+dt+"-planar_rising_bubble1.nc"
;-- open file and read variable
fref       =  addfile(diri+f00001, "r")
frefcv     =  addfile(diri+f00001cv, "r")

fieldref=fref->$fff$(tt,:,:)
fieldrefcv=frefcv->$fff$(tt,:,:)
lref=sqrt(sum(fieldref*fieldref))

err = fieldrefcv-fieldref
err = err*err
nerr = sqrt(sum(err))/lref
NOerr = sqrt(sum(err))
print("------------------ L2 error, not normalized and normalized by P frame ")
print("dt=" +dt+ " NOerr "+NOerr +" nerr = "+nerr)




end
