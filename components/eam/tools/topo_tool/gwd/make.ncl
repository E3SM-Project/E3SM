load "/lcrc/group/e3sm/ac.xie7/Analysis/NCLep/self.ncl"
begin
vars=(/"PHIS","SGH","SGH30","LANDFRAC","LANDM_COSLAT"/)
;;
fil1="final-180-mod.nc"
fil2="USGS-gtopo30_ne30np4pg2_16xdel2.c20200108.nc"
fil3="final-180-ne30pg2.nc"
system("rm -r "+fil1)
system("cp -r "+fil3+" "+fil1)
;;
ff1=addfile(fil1,"w")
ff2=addfile(fil2,"r")
;;
do i=0,4
ff1->$vars(i)$=ff2->$vars(i)$
end do


end
