Module ogwd_sub
use shr_kind_mod, only: r8 => shr_kind_r8
!use transform

contains
#if 0
subroutine OAdir(terr,ntarget,ncube,n,nvar_dir,jall,weights_lgr_index_all,weights_eul_index_all1,weights_eul_index_all2,weights_eul_index_all3,weights_all,landfrac_target,lon_cen,lat_cen,lon_terr,lat_terr,area_target,oa_target)
!use shr_kind_mod, only: r8 => shr_kind_r8
IMPLICIT NONE
integer ,intent(in)  :: ncube,ntarget,n,nvar_dir,jall,weights_lgr_index_all(jall)
integer ,intent(in)  :: weights_eul_index_all1(jall),&
                        weights_eul_index_all2(jall),&
                        weights_eul_index_all3(jall)
real(r8),intent(in)  :: weights_all(jall,1),landfrac_target(ntarget)
real(r8),intent(in)  :: terr(n)
real(r8),intent(in)  :: lon_cen(ntarget),&
                        lat_cen(ntarget),&
                        area_target(ntarget)
real(r8),intent(in)  :: lon_terr(n),lat_terr(n)
real(r8),intent(out) :: oa_target(ntarget,nvar_dir)
!local
integer  :: count,i,ix,iy,ip,ii,ip2,ip3
real(r8) :: xxterr,yyterr,zzterr,ix2,iy2,xx3,yy3,zz3,ix3,iy3
real(r8) :: wt,xhds(ntarget),yhds(ntarget),zhds(ntarget),hds(ntarget),OA_var(ntarget)
real(r8) :: xbar(ntarget),ybar(ntarget),zbar(ntarget),lon_bar(ntarget),lat_bar(ntarget)
real(r8) :: rad,theta1
real(r8) :: OAlon(ntarget),OAlat(ntarget),OArad(ntarget),OAx1,OAy1,OAz1

xhds=0.0_r8
yhds=0.0_r8
zhds=0.0_r8
hds=0.0_r8

xbar=0.0_r8
ybar=0.0_r8
zbar=0.0_r8
lon_bar=0.0_r8
lat_bar=0.0_r8
OA_var=0.0_r8

!terr=1
!stop
do count=1,jall
        i   = weights_lgr_index_all(count)
        ix  = weights_eul_index_all1(count)!,1)
        iy  = weights_eul_index_all2(count)!,2)
        ip  = weights_eul_index_all3(count)!,3)
        !
        ! convert to 1D indexing of cubed-sphere
        ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!
        wt = weights_all(count,1)
        rad=4.0_r8*atan(1.0_r8)/180.0_r8
        call CubedSphereABPFromRLL(lon_terr(ii)*rad,lat_terr(ii)*rad,ix2,iy2,ip2,.true.)
        call CubedSphereXYZFromABP(ix2,iy2,ip2,xxterr,yyterr,zzterr)
!#if 0
        xhds(i)=xhds(i)+xxterr*terr(ii)*wt
        yhds(i)=yhds(i)+yyterr*terr(ii)*wt
        zhds(i)=zhds(i)+zzterr*terr(ii)*wt
        hds(i) =hds(i)+terr(ii)*wt

        !masscenter for every coarse grid
        !on Cartesian coord
        !looking the h as ro
        xbar(i)=xhds(i)/hds(i)
        ybar(i)=yhds(i)/hds(i)
        zbar(i)=zhds(i)/hds(i)

        call CubedSphereABPFromRLL(lon_cen(i)*rad,lat_cen(i)*rad,&
                                   ix3,iy3,ip3,.true.)
        call CubedSphereXYZFromABP(ix3,iy3,ip3,xx3,yy3,zz3)
        !under Cartesian, the variability of the scale in the wind
        !direction is the sqrt(x^2+y^2+z^2), the scale of the orthogonal
        !3 directions
        !then it is only a matter of using the original formula
        !in the single direction
        OA_var(i)=OA_var(i)+wt/area_target(i)&
                *((xxterr-xx3)**2+(yyterr-yy3)**2+(zzterr-zz3)**2)
        OAx1=(xx3-xbar(i))/sqrt(OA_var(i))
        OAy1=(yy3-ybar(i))/sqrt(OA_var(i))
        OAz1=(zz3-zbar(i))/sqrt(OA_var(i))
        !assuming a small change in lon_cen to lon_bar
        !so it does not matter whether lon_cen or lon_bar
        !thus we change onto lat-lon grid vector in target gridcell
#if 0
        OArad(i)= OAx1*sin(lat_cen(i)*rad)*cos(lon_cen(i)*rad)&
                 +OAy1*sin(lat_cen(i)*rad)*sin(lon_cen(i)*rad)&
                 +OAz1*cos(lat_cen(i)*rad)
        OAlat(i)= OAx1*cos(lat_cen(i)*rad)*cos(lon_cen(i)*rad)&
                 +OAy1*cos(lat_cen(i)*rad)*sin(lon_cen(i)*rad)&
                 -OAz1*sin(lat_cen(i)*rad)
        OAlon(i)=-OAx1*sin(lon_cen(i)*rad)&
                 +OAy1*cos(lon_cen(i)*rad)
#endif
        !all lat_cen must use (90-lat_cen) since we only have 
        !latitude rather than colatitude
        !this is equivalent to using induction formula sin(90-lat)=cos(lat)
        !latitude is opposite of colatitude, thus OAlat is negative
        OAlat(i)=-(OAx1*sin(lat_cen(i)*rad)*cos(lon_cen(i)*rad)&
                  +OAy1*sin(lat_cen(i)*rad)*sin(lon_cen(i)*rad)&
                  -OAz1*cos(lat_cen(i)*rad))
        OAlon(i)= -OAx1*sin(lon_cen(i)*rad)&
                  +OAy1*cos(lon_cen(i)*rad)
#if 0
        theta1=0.
        oa_target(i,1) = OAlon(i)*cos(theta1*rad)+OAlat(i)*sin(theta1*rad)
        theta1=90.
        oa_target(i,2) = OAlon(i)*cos(theta1*rad)+OAlat(i)*sin(theta1*rad)
        theta1=45.
        oa_target(i,3)=  OAlon(i)*cos(theta1*rad)+OAlat(i)*sin(theta1*rad)
        theta1=360.-45.
        oa_target(i,4)=  OAlon(i)*cos(theta1*rad)+OAlat(i)*sin(theta1*rad)
#endif
!#if 0
        !reverse in order to be 
        !(2,ntarget),OAx,OAy
        oa_target(i,1) = OAlon(i)
        oa_target(i,2) = OAlat(i)

!#endif
        !landfrac may cause coastal area par to diminish
        !oa_target(i,:) = oa_target(i,:)*landfrac_target(i)
enddo
        !takeout abnormal values
!#if 0
        where(abs(oa_target)<.001_r8.or.&
              abs(oa_target).gt.1e+7) oa_target=0.0_r8
        where(abs(oa_target).gt.1) oa_target=1.0_r8
        where(oa_target.ne.oa_target) oa_target=0.0_r8

!#endif
end subroutine OAdir
!============================================================
subroutine OAorig(terr,ntarget,ncube,n,jall,weights_lgr_index_all,weights_eul_index_all1,weights_eul_index_all2,weights_eul_index_all3,weights_all,landfrac_target,lon_terr,lat_terr,area_target,oa_target)
!use shr_kind_mod, only: r8 => shr_kind_r8
IMPLICIT NONE
integer ,intent(in)  :: ncube,ntarget,n,jall,weights_lgr_index_all(jall),weights_eul_index_all1(jall),weights_eul_index_all2(jall),weights_eul_index_all3(jall)
real(r8),intent(in)  :: weights_all(jall,1),terr(n)
real(r8),intent(in)  :: landfrac_target(ntarget),lon_terr(n),lat_terr(n),area_target(ntarget)
real(r8),intent(out) :: oa_target(ntarget,4)
!local
real(r8) :: xh(ntarget),yh(ntarget),height(ntarget),modexcoords(ntarget),modeycoords(ntarget),avgx(ntarget),avgy(ntarget),varx(ntarget),vary(ntarget),OAx,OAy,theta1,rad
integer(r8) :: count,i,ix,iy,ip,ii
real(r8)    :: wt

        xh=0.0_r8
        yh=0.0_r8
        height=0.0_r8
        modexcoords=0.0_r8
        modeycoords=0.0_r8
        avgx=0.0_r8
        avgy=0.0_r8
        varx=0.0_r8
        vary=0.0_r8
        OAx=0.0_r8
        OAy=0.0_r8
        theta1=0.0_r8
        rad=0.0_r8
        
do count=1,jall
        i    = weights_lgr_index_all(count)
        ix  = weights_eul_index_all1(count)!,1)
        iy  = weights_eul_index_all2(count)!,2)
        ip  = weights_eul_index_all3(count)!,3)
        !
        ! convert to 1D indexing of cubed-sphere
        !
        ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!
        wt = weights_all(count,1)
        !for OA
        avgx(i)=avgx(i)+wt/area_target(i)*lon_terr(ii)
        avgy(i)=avgy(i)+wt/area_target(i)*lat_terr(ii)
enddo


do count=1,jall
        i   = weights_lgr_index_all(count)
        ix  = weights_eul_index_all1(count)!,1)
        iy  = weights_eul_index_all2(count)!,2)
        ip  = weights_eul_index_all3(count)!,3)
        !
        ! convert to 1D indexing of cubed-sphere
        !
        ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!
        wt = weights_all(count,1)
        !mode for both dim      
        xh(i)=xh(i)+wt/area_target(i)*lon_terr(ii)*terr(ii)
        yh(i)=yh(i)+wt/area_target(i)*lat_terr(ii)*terr(ii)
        height(i)=height(i)+wt/area_target(i)*terr(ii)
        modexcoords(i)=xh(i)/(height(i))!+1e-14)
        modeycoords(i)=yh(i)/(height(i))!+1e-14)

        varx(i)=varx(i)+(wt/area_target(i))*(lon_terr(ii)-avgx(i))**2
        vary(i)=vary(i)+(wt/area_target(i))*(lat_terr(ii)-avgy(i))**2
        !OAx,OAy
        OAx=landfrac_target(i)*(avgx(i)-modexcoords(i))/sqrt(varx(i))
        OAy=landfrac_target(i)*(avgy(i)-modeycoords(i))/sqrt(vary(i))

        rad=4.0*atan(1.0)/180.0
        theta1=0.
        oa_target(i,1) = OAx*cos(theta1*rad)+OAy*sin(theta1*rad)
        theta1=90.
        oa_target(i,2) = OAx*cos(theta1*rad)+OAy*sin(theta1*rad)
        theta1=45.
        oa_target(i,3)=  OAx*cos(theta1*rad)+OAy*sin(theta1*rad)
        theta1=360.-45.
        oa_target(i,4)=  OAx*cos(theta1*rad)+OAy*sin(theta1*rad)
        oa_target(i,:)=  oa_target(i,:)*landfrac_target(i)
enddo
        !takeout abnormal values
        where(abs(oa_target)<.001_r8.or.abs(oa_target).gt.1e+7) oa_target=0.0
        where(abs(oa_target).gt.1) oa_target=0.0
        where(oa_target.ne.oa_target) oa_target=0.0
end subroutine OAorig
!#endif
!===================================
subroutine OC(terr,ntarget,ncube,n,jall,weights_lgr_index_all,weights_eul_index_all1,weights_eul_index_all2,weights_eul_index_all3,weights_all,landfrac_target,area_target,sgh_target,terr_target,oc_target)
!use shr_kind_mod, only: r8 => shr_kind_r8
IMPLICIT NONE
integer ,intent(in)  :: ncube,ntarget,n,jall,weights_lgr_index_all(jall),weights_eul_index_all1(jall),weights_eul_index_all2(jall),weights_eul_index_all3(jall)
real(r8),intent(in)  :: weights_all(jall,1)
real(r8),intent(in)  :: landfrac_target(ntarget),area_target(ntarget),sgh_target(ntarget),terr_target(ntarget),terr(n)
real(r8),intent(out) :: oc_target(ntarget)
!local 
integer  :: count,i,ix,iy,ip,ii
real(r8) :: wt

        oc_target=0.0_r8
        do count=1,jall
        i   = weights_lgr_index_all(count)
        ix  = weights_eul_index_all1(count)!,1)
        iy  = weights_eul_index_all2(count)!,2)
        ip  = weights_eul_index_all3(count)!,3)
        ! convert to 1D indexing of cubed-sphere
        ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!
        wt = weights_all(count,1)
        oc_target(i) = oc_target(i)+(wt/area_target(i))*((terr_target(i)-terr(ii))**4)/(sgh_target(i)**4)
        oc_target(i) = oc_target(i) * landfrac_target(i)
        enddo

      where(abs(oc_target)<.001_r8.or.abs(oc_target).gt.1e+7) oc_target=0.0_r8
      where(abs(sgh_target).eq.0.0_r8) oc_target=0.0_r8
end subroutine OC
!========================
subroutine OLorig(terr,ntarget,ncube,n,jall,weights_lgr_index_all,weights_eul_index_all1,weights_eul_index_all2,weights_eul_index_all3,weights_all,landfrac_target,lon_terr,lat_terr,area_target,sgh_target,target_center_lat,target_center_lon,target_corner_lat_deg,target_corner_lon_deg,ol_target)
!use shr_kind_mod, only: r8 => shr_kind_r8
IMPLICIT NONE
integer,intent(in)  :: ncube,ntarget,n,jall,weights_lgr_index_all(jall),weights_eul_index_all1(jall),weights_eul_index_all2(jall),weights_eul_index_all3(jall)
real(r8),intent(in)  :: weights_all(jall,1)
real(r8),intent(in)  :: landfrac_target(ntarget),area_target(ntarget),sgh_target(ntarget),terr(n),lon_terr(n),lat_terr(n)
real(r8),intent(in)  :: target_center_lat(ntarget),target_center_lon(ntarget),target_corner_lat_deg(4,ntarget),target_corner_lon_deg(4,ntarget)
real(r8),intent(out) :: ol_target(ntarget,4)
!local 
integer  :: count,i,ix,iy,ip,ii
real(r8) :: wt,terr_if,Nw(4,ntarget),area_target_par(4,ntarget),j


        ol_target=0.0_r8
        Nw=0.0_r8
        area_target_par=0.0_r8

        do count=1,jall
        i   = weights_lgr_index_all(count)
        ix  = weights_eul_index_all1(count)!,1)
        iy  = weights_eul_index_all2(count)!,2)
        ip  = weights_eul_index_all3(count)!,3)
        ! convert to 1D indexing of cubed-sphere
        ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!
        wt = weights_all(count,1)
        !determine terr_if
        terr_if=0._r8
        if (terr(ii).GT.(1116.2-0.878*sgh_target(i))) terr_if=1.
! (1):  the lower left corner
! (2):  the lower right corner
! (3):  the upper right corner
! (4):  the upper left corner
     !OL1
     if (lat_terr(ii) &!(ii)&
     .GT.(target_corner_lat_deg(1,i)+target_center_lat(i))/2..and. &
         lat_terr(ii) &!(ii)&
     .LT.(target_corner_lat_deg(4,i)+target_center_lat(i))/2.) then
        Nw(1,i)=Nw(1,i)+wt*terr_if
        area_target_par(1,i)=area_target_par(1,i)+wt
     endif
        
     !OL2
     if (lon_terr(ii) &!(ii)&
     .GT.(target_corner_lon_deg(1,i)+target_center_lon(i))/2..and. &
        lon_terr(ii) &!(ii)&
     .LT.(target_corner_lon_deg(3,i)+target_center_lon(i))/2.) then
        Nw(2,i)=Nw(2,i)+wt*terr_if
        area_target_par(2,i)=area_target_par(2,i)+wt
     end if


     !OL3
     if (lon_terr(ii) &!(ii)&
     .LT.target_center_lon(i).and. &
      lat_terr(ii) &!(ii)&
     .LT.target_center_lat(i).or.  &
      lon_terr(ii) &!(ii)&
     .GT.target_center_lon(i).and. &
      lat_terr(ii) &!(ii)&  
     .GT.target_center_lat(i)) then
        Nw(3,i)=Nw(3,i)+wt*terr_if
        area_target_par(3,i)=area_target_par(3,i)+wt
     end if


     !OL4
     if (lat_terr(ii) & !(ii)&
     .GT.target_center_lat(i).and. &
     lon_terr(ii) & !(ii)&
     .LT.target_center_lon(i).or.  &
     lat_terr(ii) & !(ii)&
     .LT.target_center_lat(i).and. &
     lon_terr(ii) & !(ii)&
     .GT.target_center_lon(i)) then
        Nw(4,i)=Nw(4,i)+wt*terr_if
        area_target_par(4,i)=area_target_par(4,i)+wt
     end if

         !Nw(4,i)=Nw(4,i)+wt*terr_if
        !area_target_par(4,i)=area_target_par(4,i)+wt



        do j=1,4
        ol_target(i,j)=Nw(j,i)/(area_target_par(j,i)+1e-14)!Nt(i)!/2.)
        enddo
        ol_target(i,:)=ol_target(i,:)*landfrac_target(i)
        end do
        where(abs(ol_target)<.001_r8.or.abs(ol_target).gt.1e+7) ol_target=0.0_r8
end subroutine OLorig
#endif
!=====================
!===================================================================
!=====================
!#if 0
subroutine OLgrid(terr,terrx,terry,wt,b,a,n,theta_in,hc,OLout)
!use physconst,      only: rh2o,zvir,pi,rearth
!use abortutils
!Xie add
IMPLICIT NONE
integer,intent(in)  :: n
real(r8),intent(in) :: hc,wt(n),terr(n),a,b,theta_in!a dy,b dx
real(r8),intent(in) :: terrx(n),terry(n)
real(r8),intent(out) :: OLout
real(r8) :: theta,theta1,theta2,rad,interval
real(r8) :: terr_count(n),terr_whole_count(n),cx(n),c1,c2
!local
integer  :: i
real(r8) :: j
        terr_count=0.0_r8
        terr_whole_count=0.0_r8
        c1=0.0_r8
        c2=0.0_r8
        cx=0.0_r8
        !determine an acute theta in the triangle
        !or minus 180 equilvalent acute angle
        !then turn into radian
        rad=4.0_r8*atan(1.0_r8)/180.0_r8
        interval=0.0_r8
        
        !initialize
        theta1=0.0_r8
        theta2=0.0_r8
        !set inside -360~360
        !this adds robustness of the scheme to different angle
        theta1=MOD(theta_in,360._r8)
        !set negative axis into 0~360
        if (theta1.ge.-360._r8.and.theta1.lt.0._r8) then
        theta1=theta1+360._r8
        endif
        !now we have only 0~360 angle
        if (theta1.ge.  0._r8.and.theta1.le. 90._r8) then
        theta=theta1*rad
        theta2=theta1
        else if (theta1.gt.  90._r8.and.theta1.le. 180._r8) then
        theta=(180._r8-theta1)*rad
        theta2=180._r8-theta1
        else if (theta1.gt. 180._r8.and.theta1.le. 270._r8) then
        theta=(theta1-180._r8)*rad
        theta2=theta1-180._r8
        !we only use 0~180, so this makes similar to 1st and 2nd quadrant
        else if (theta1.gt. 270._r8.and.theta1.le. 360._r8) then
        theta=(360._r8-theta1)*rad
        theta2=360._r8-theta1
        !we only use 0~180, so this makes similar to 1st and 2nd quadrant
        endif
        !we use theta2 to judge instead
        !theta2=theta1
        !theta1=theta2
        !we restrict the angle in the first and second quadrant
        !the third and fourth for OL are similar when theta is 
        !transformed by minus pi(180)
                !two parallel lines are included
                !xsin(theta)-ycos(theta)=c1
                !xsin(theta)-ycos(theta)=c2
                !xsin(theta)-ycos(theta)=cx,min(c1,c2)<cx<max(c1,c2)
                !so there are 6 conditions considered here
        !0<theta1 or theta2<atan(a/(2.*b))/rad
                !we use two points on the left to
                !determine two lines respectively
                !(-0.5*b,-0.5*a+l)
                !(-0.5*b,0.5*a-l-b*tan(theta))
                !l=(1./4.)*a-(1./2.)*b*tan(theta)
                !c1=-0.5*b*sin(theta)-(-0.5*a+l)*cos(theta)
                !c1=-0.5*b*sin(theta)-(-0.5*a+0.25*a-0.5*b*tan(theta))&
                !   *cos(theta)
                !c1=-0.5*b*sin(theta)-(-0.25*a-0.5*b*tan(theta))*cos(theta)
                !c1=-0.5*b*sin(theta)+0.25*a*cos(theta)+0.5*b*tan(theta)
                !c1= 0.25*a*cos(theta)
                !c2=-0.5*b*sin(theta)-(0.5*a-l-b*tan(theta))*cos(theta)
                !c2=-0.5*b*sin(theta)-(0.5*a-0.25*a+&
                !    0.5*b*tan(theta)-b*tan(theta))*cos(theta)
                !c2=-0.5*b*sin(theta)-(0.25*a-0.5*b*tan(theta))*cos(theta)
                !c2=-0.5*b*sin(theta)-0.25*a*cos(theta)+0.5*b*sin(theta)
                !c2=-0.25*a*cos(theta)
        !atan(a/(2.*b))/rad<theta1 or theta2<atan(2*a/b)/rad
                !l1=sqrt(a*b/(2*tan(theta)))*tan(theta)
                !right to the angle
                !l2=sqrt(a*b/(2*tan(theta)))           
                !length near the angle
                !we use two points on the left to
                !determine two lines respectively
                !(-0.5*b, 0.5*a-l1)
                !( 0.5*b-l2,-0.5*a)
                !c1=-0.5*b*sin(theta)-(0.5*a-l1)*cos(theta)
                !c1=-0.5*b*sin(theta)-(0.5*a-sqrt(a*b/(2*tan(theta)))&
                !         *tan(theta))*cos(theta)
                !c1=-0.5*b*sin(theta)-0.5*a*cos(theta)+&
                !               sqrt(a*b/(2*tan(theta)))*sin(theta))
                !c1=-0.5*b*sin(theta)-0.5*a*cos(theta)+&
                !               sqrt(a*b*sin(theta)*cos(theta)/2.)
                !c1=-0.5*b*sin(theta)-0.5*a*cos(theta)+&
                !                sqrt(a*b*sin(2*theta)/4.
                !c2=(0.5*b-l2)*sin(theta)-(-0.5*a)*cos(theta)
                !c2=(0.5*b-sqrt(a*b/(2*tan(theta))))*sin(theta)-&
                !  (-0.5*a)*cos(theta)
                !c2=0.5*b*sin(theta)-sqrt(a*b*sin(theta)*cos(theta)/2)+&
                !   0.5*a*cos(theta)
                !c2=0.5*b*sin(theta)+0.5*a*cos(theta)-&
                !      sqrt(a*b*sin(theta)*cos(theta)/2.)
                !c2=0.5*b*sin(theta)+0.5*a*cos(theta)-&
                !               sqrt(a*b*sin(2*theta)/4.)
        !atan(2*a/b)/rad<theta1 or theta2<=90                      
                !l=(1./4.)*b-(1./2.)*a/tan(theta)                          
                !use upper two points to check                             
                !(-0.5*b+l+a/tan(theta),0.5*a)                             
                !( 0.5*b-l,0.5*a)    
                !c1=(-0.5*b+l+a/tan(theta))*sin(theta)-0.5*a*cos(theta)    
                !c1=(-0.5*b+0.25*b-0.5*a/tan(theta)+&                      
                !          a/tan(theta))*sin(theta)-0.5*a*cos(theta)       
                !c1=(-0.25*b+0.5*a/tan(theta))*sin(theta)-0.5*a*cos(theta) 
                !c1=-0.25*b*sin(theta)+0.5*a*cos(theta)-0.5*a*cos(theta)   
                !c1=-0.25*b*sin(theta)
                !c2=(0.5*b-l)*sin(theta)-0.5*a*cos(theta)                  
                !c2=(0.5*b-0.25*b+0.5*a/tan(theta))*sin(theta)&            
                !                -0.5*a*cos(theta)                         
                !c2=(0.25*b+0.5*a/tan(theta))*sin(theta)&                  
                !                      -0.5*a*cos(theta)
                !c2= 0.25*b*sin(theta)+0.5*a*cos(theta)-0.5*a*cos(theta)
                !c2= 0.25*b*sin(theta)
        !90<theta1 or theta2<180-atan(2*a/b)/rad
                !l=(1./4.)*b-(1./2.)*a/tan(theta)
                !we use left two points
                !(-0.5*b+l,0.5*a)
                !( 0.5*b-l-a/tan(theta),0.5*a)
                !c1=(-0.5*b+l)*sin(theta)-(0.5*a)*cos(pi-theta)
                !c1=(-0.5*b+0.25*b-0.5*a/tan(theta))*sin(theta)+&
                !                  0.5*a*cos(theta)
                !c1=(-0.25*b-0.5*a/tan(theta))*sin(theta)+&
                !            0.5*a*cos(theta)
                !c1= -0.25*b*sin(theta)-0.5*a*cos(theta)+0.5*a*cos(theta)
                !c1= -0.25*b*sin(theta)
                !c2=(0.5*b-l-a/tan(theta))*sin(theta)-(0.5*a)*cos(pi-theta)
                !c2=(0.5*b-0.25*b+0.5*a/tan(theta)-a/tan(theta))*&         
                !                       sin(theta)+0.5*a*cos(theta)        
                !c2=(0.25*b-0.5*a/tan(theta))*sin(theta)+0.5*a*cos(theta)  
                !c2=0.25*b*sin(theta)-0.5*a*cos(theta)+0.5*a*cos(theta)    
                !c2=0.25*b*sin(theta)
       !180-atan(2*a/b)/rad<theta1 or theta2<180-atan(a/(2*b))/rad)
                !l1=sqrt(a*b/(2*tan(theta)))*tan(theta)
                !right to the angle
                !l2=sqrt(a*b/(2*tan(theta)))           
                !length near the angle
                !(-0.5*b,-0.5*a+l1)
                !( 0.5*b-l2,0.5*a)
                !c1=(-0.5*b)*sin(pi-theta)-(-0.5*a+l1)*cos(pi-theta)
                !c1=-0.5*b*sin(theta)+&
                !  (-0.5*a+sqrt(a*b/(2*tan(theta)))*tan(theta))*cos(theta)
                !c1=-0.5*b*sin(theta)-0.5*a*cos(theta)+&
                !  sqrt(a*b/(2*tan(theta)))*sin(theta)
                !c1=-0.5*b*sin(theta)-0.5*a*cos(theta)+&                   
                !       sqrt(a*b*sin(theta)*cos(theta)/2.))                
                !c1=-0.5*b*sin(theta)-0.5*a*cos(theta)+&                   
                !                sqrt(a*b*sin(2*theta)/4.) 
                !c2=(0.5*b-l2)*sin(pi-theta)-(0.5*a)*cos(pi-theta)
                !c2=(0.5*b-sqrt(a*b/(2*tan(theta))))*sin(theta)+&
                !                0.5*a*cos(theta)
                !c2=0.5*b*sin(theta)-sqrt(a*b*sin(theta)*cos(theta)/2.)+&
                !   0.5*a*cos(theta)
                !c2=0.5*b*sin(theta)+0.5*a*cos(theta)-&
                !      sqrt(a*b*sin(theta)*cos(theta)/2.)
                !c2=0.5*b*sin(theta)+0.5*a*cos(theta)-&
                !               sqrt(a*b*sin(2*theta)/4.)
        !180-atan(a/(2*b))/rad<theta1 or theta2<180 
                !l=(1./4.)*a-(1./2.)*b*tan(theta)                          
                !(-0.5*b, 0.5*a-l)   
                !(-0.5*b,-0.5*a+l+b*tan(theta))                            
                !c1=(-0.5*b)*sin(pi-theta)-(0.5*a-l)*cos(pi-theta)         
                !c1=(-0.5*b)*sin(theta)+(0.5*a-0.25*a+0.5*b*tan(theta))*&  
                !                                           cos(theta)     
                !c1=-0.5*b*sin(theta)+(0.25*a*cos(theta)+0.5*b*sin(theta)) 
                !c1=0.25*a*cos(theta)                                      
                !c2=(-0.5*b)*sin(theta)-(-0.5*a+l+b*tan(theta))&           
                !        *cos(pi-theta)                                    
                !c2=-0.5*b*sin(theta)+&                                    
                !(-0.5*a+0.25*a-0.5*b*tan(theta)+b*tan(theta))*cos(theta)  
                !c2=-0.5*b*sin(theta)+(-0.25*a+0.5b*tan(theta))*cos(theta)
                !c2=-0.5*b*sin(theta)-0.25*a*cos(theta)+0.5*b*sin(theta)
                !c2=-0.25*a*cos(theta)

!0<theta1 or theta2<atan(a/(2.*b))/rad
        !c1= 0.25*a*cos(theta)
        !c2=-0.25*a*cos(theta)
!atan(a/(2.*b))/rad<theta1 or theta2<atan(2*a/b)/rad
        !c1=-0.5*b*sin(theta)-0.5*a*cos(theta)+sqrt(a*b*sin(2*theta)/4.)
        !c2= 0.5*b*sin(theta)+0.5*a*cos(theta)-sqrt(a*b*sin(2*theta)/4.)
!atan(2*a/b)/rad<theta1 or theta2<=90 
        !c1=-0.25*b*sin(theta)
        !c2= 0.25*b*sin(theta)
!90<theta1 or theta2<180-atan(2*a/b)/rad
        !c1=-0.25*b*sin(theta)
        !c2= 0.25*b*sin(theta)
!180-atan(2*a/b)/rad<theta1 or theta2<180-atan(a/(2*b))/rad)
        !c1=-0.5*b*sin(theta)-0.5*a*cos(theta)+sqrt(a*b*sin(2*theta)/4.)
        !c2= 0.5*b*sin(theta)+0.5*a*cos(theta)-sqrt(a*b*sin(2*theta)/4.)
!180-atan(a/(2*b))/rad<theta1 or theta2<180
        !c1= 0.25*a*cos(theta)
        !c2=-0.25*a*cos(theta)
!6 merged into 3 conditions
        !0<theta1 or theta2<atan(a/(2.*b))/rad
        !180-atan(a/(2*b))/rad<theta1 or theta2<180
        !c1=-0.25*a*cos(theta)
        !c2= 0.25*a*cos(theta)
        !atan(2*a/b)/rad<theta1 or theta2<=90
        !90<theta1 or theta2<180-atan(2*a/b)/rad
        !c1=-0.25*b*sin(theta)
        !c2= 0.25*b*sin(theta)
        !atan(a/(2.*b))/rad<theta1 or theta2<atan(2*a/b)/rad
        !180-atan(2*a/b)/rad<theta1 or theta2<180-atan(a/(2*b))/rad)
        !c1=-0.5*b*sin(theta)-0.5*a*cos(theta)+sqrt(a*b*sin(2*theta)/4.)
        !c2= 0.5*b*sin(theta)+0.5*a*cos(theta)-sqrt(a*b*sin(2*theta)/4.)
                if      (theta1.ge.  0._r8.and.theta1.lt.atan2(a,2._r8*b)/rad.or.&
                         theta2.ge.  0._r8.and.theta2.lt.atan2(a,2._r8*b)/rad.or.&
                         theta1.ge.180._r8-atan2(a,2_r8*b)/rad.and.theta1.lt.180._r8.or.&
                         theta2.ge.180._r8-atan2(a,2_r8*b)/rad.and.theta2.lt.180._r8)&
                then
interval=1
                c1=-0.25_r8*a*cos(theta)
                c2= 0.25_r8*a*cos(theta)
                else if (theta1.ge.atan2(a,2_r8*b)/rad.and.&
                         theta1.lt.atan2(2_r8*a,b)/rad.or.&
                         theta2.ge.atan2(a,2_r8*b)/rad.and.&
                         theta2.lt.atan2(2_r8*a,b)/rad.or.&
                         theta1.ge.180._r8-atan2(2_r8*a,b)/rad.and.&
                         theta1.lt.180._r8-atan2(a,2_r8*b)/rad.or.&
                         theta2.ge.180._r8-atan2(2_r8*a,b)/rad.and.&
                         theta2.lt.180._r8-atan2(a,2_r8*b)/rad) then
interval=2
                c1=-0.5_r8*b*sin(theta)-0.5_r8*a*cos(theta)+&
                                sqrt(a*b*sin(2_r8*theta)/4._r8)
                c2= 0.5_r8*b*sin(theta)+0.5_r8*a*cos(theta)-&
                                sqrt(a*b*sin(2_r8*theta)/4._r8)
                else if (theta1.ge.atan2(2_r8*a,b)/rad.and.theta1.lt.90._r8.or.&
                         theta2.ge.atan2(2_r8*a,b)/rad.and.theta2.lt.90._r8.or.&
                         theta1.ge.90._r8.and.theta1.lt.180._r8-atan2(2_r8*a,b)/rad&
                     .or.theta2.ge.90._r8.and.theta2.lt.180._r8-atan2(2_r8*a,b)/rad)&
                then
interval=3
                c1=-0.25_r8*b*sin(theta)
                c2= 0.25_r8*b*sin(theta)
                endif
                !determine two line functions
                cx=terrx*sin(theta_in*rad)-terry*cos(theta_in*rad)
                !assuming rectangle grid or ladder-shape
                !would be pretty similar
                !since in 1.4,the max difference is of 0.02
                !although the above expression is actually ladder-shape
                !and use a rectangle with the same area
                !to derive c1 and c2
                where (cx.ge.min(c1,c2).and.cx.le.max(c1,c2)) &
                terr_whole_count=1._r8
                where (cx.ge.min(c1,c2).and.cx.le.max(c1,c2) &
                .and.terr.ge.hc) &
                terr_count=1._r8


                !deals with noise that affects OL
                !when there are no terrx center points
                !in the two lines in the center
                if (sum(wt*terr_whole_count).eq.0._r8) then
                !enlarge about 5 times interval
                !to include more points and avoid 
                !noise
                j=5._r8
                where ((cx+j*abs(c1-c2)).ge.min(c1,c2)&
                .and.  (cx-j*abs(c1-c2)).le.max(c1,c2)) &
                terr_whole_count=1._r8
                where ((cx+j*abs(c1-c2)).ge.min(c1,c2)&
                .and.  (cx-j*abs(c1-c2)).le.max(c1,c2) &
                .and.terr.ge.hc) &                      
                terr_count=1._r8
                endif


                !output
                OLout=sum(wt*terr_count)/sum(wt*terr_whole_count)!
                !when insufficient number of grids to support anisotropic
                !there may be a strong 0 or 1 jump between directions
                !set to istotropic instead
                if (n.le.20) then
                terr_whole_count=1._r8
                where(terr.gt.hc) terr_count=1._r8
                OLout=sum(wt*terr_count)/sum(wt*terr_whole_count)
                endif
                !take out pole point
                !where dx is usually smaller than 1m
                !if (a.lt.1.or.b.lt.1) OLout=0.0_r8
                !take out NaN
                if (OLout.ne.OLout) OLout=0.0_r8
end subroutine OLgrid
!===================================================================
#if 0
subroutine OLdir(terr,ntarget,ncube,n,jall,nlon,nlat,indexb,nvar_dir,weights_lgr_index_all,weights_eul_index_all1,weights_eul_index_all2,weights_eul_index_all3,weights_all,landfrac_target,lon_cen,lat_cen,lon_terr,lat_terr,sgh_target,ol_target,terrout,dxy)
IMPLICIT NONE
integer ,intent(in)  :: ncube,ntarget,n,jall,indexb
integer ,intent(in)  :: weights_lgr_index_all(jall) ,&
                        weights_eul_index_all1(jall),&
                        weights_eul_index_all2(jall),&
                        weights_eul_index_all3(jall),&
                        nlon,nlat,nvar_dir
real(r8),intent(in)  :: weights_all(jall,1),landfrac_target(ntarget)
real(r8),intent(in)  :: terr(n),lon_terr(n),lat_terr(n)
real(r8),intent(in)  :: lon_cen(ntarget),&
                        lat_cen(ntarget),&
                        sgh_target(ntarget)
real(r8),intent(out) :: ol_target(ntarget,nvar_dir)
real(r8),intent(out) :: terrout(4,ntarget,indexb)
real(r8),intent(out) :: dxy(ntarget,nvar_dir)


!local
!1 is lb,2 is ub
integer  :: index_b(3,ntarget),index_jall(jall)
integer,allocatable :: indexii_b(:,:)
integer  :: ix,iy,ip,i,count,alloc_error,j
real(r8) :: xterr(n),yterr(n),dx(ntarget),dy,hc(ntarget),theta1(nvar_dir)
real(r8) :: xterr_cen(ntarget),yterr_cen(ntarget),rad
real(r8) :: reflon_terr(n),reflat_terr(n)!,lon_terr2(n)
REAL(r8), PARAMETER :: pi    = 3.14159265358979323846264338327


        !initialize
        index_b=0
        indexii_b=0
        index_jall=0
        xterr=0.0_r8
        yterr=0.0_r8
        reflon_terr=0.0_r8
        reflat_terr=0.0_r8
        rad=4.0_r8*atan(1.0_r8)/180.0_r8
        !determine correspondent upper and 
        !lower bound of the large grid
        print*,"begin ol index_b"
        do count=1,jall
        i   = weights_lgr_index_all(count)
        index_b(3,i)=index_b(3,i)+1
        enddo
        !cumsum to form upper and lower bound of index_b
        !1 for lower bound, 2 for upper bound
        do i=1,ntarget
        index_b(2,i)=sum(index_b(3,1:i))
        enddo
        index_b(1,1)=1
        do i=2,ntarget
        index_b(1,i)=index_b(2,i-1)+1
        enddo
        print*,"after index_b"
        !leave largest possible dimension first
        allocate(indexii_b(maxval(index_b(3,:)),ntarget),stat=alloc_error)
        print*,"maxval(index_b(3,:)",maxval(index_b(3,:))
        indexii_b=0

        do i=1,ntarget               
                do j=1,index_b(3,i)  
                index_jall(index_b(1,i)+j-1)=j                             
                enddo                
        enddo
        print*,"after index_jall"
        !get correspondent ii for ub and lb                                
        do count=1,jall
        i   = weights_lgr_index_all(count)
        ix  = weights_eul_index_all1(count)!,1)
        iy  = weights_eul_index_all2(count)!,2)
        ip  = weights_eul_index_all3(count)!,3)
        ! convert to 1D indexing of cubed-sphere
        indexii_b(index_jall(count),i) = (ip-1)*ncube*ncube+(iy-1)*ncube+ix
        enddo
        print*,"convert indexii_b"
        !input small grids and make OL for individual large grid
        !critical height using empirical function
        hc=1116.2_r8-0.878_r8*sgh_target 
        !get terrx terry for the small grids
        !transform lon lat grid to distance to (0,0)
                !in case the grid spans 0 line
                !lon_terr2=lon_terr
                !where(lon_terr2.gt.180._r8) lon_terr2=lon_terr2-360._r8
!#if 0
        do i=1,ntarget
        reflon_terr(indexii_b(:index_b(3,i),i))=&
           lon_terr(indexii_b(:index_b(3,i),i))-lon_cen(i)
        reflat_terr(indexii_b(:index_b(3,i),i))=&
           lat_terr(indexii_b(:index_b(3,i),i))-lat_cen(i)
        enddo
        !take out boundary problem points
        !after that, the problem is left with
        !-45 where these points are in the pole
        !and are set out later elsewhere
        where(reflon_terr.gt.350) reflon_terr=reflon_terr-360
        where(reflon_terr.lt.-350) reflon_terr=reflon_terr+360

        print*,"get reflonlat"
        !determine real length on cartesian coordinate
        xterr=reflon_terr*cos(lat_terr*rad)
        yterr=reflat_terr
        !get dx dy for the large grid
        !dx is a function of latitude
        !at this time, we do not need real length
        !so R is set to normalized 1
!dx=rearth*cos(rlat)*(2._r8*pi/256._r8),dy=rearth*(pi/(128._r8-1._r8))
        dx=(2._r8*pi/real(nlon,kind=r8))*cos(lat_cen*rad)
        dy=(      pi/real(nlat-1,kind=r8))

        !test 4 direction for the time
        !only needs 0-180 half of the axis
        do j=1,ntarget
                do i=1,nvar_dir
                theta1(i)=(180._r8/nvar_dir)*real(i-1,kind=r8)
                call dxygrid(6.37100e6_r8*dx(j),6.37100e6_r8*dy,theta1(i),dxy(j,i))
                enddo
        enddo

        print*,"before into OLgrid"
        !input into every large grid
        do j=1,nvar_dir
                do i=1,ntarget
        call OLgrid(terr(indexii_b(:index_b(3,i),i)),&
                   xterr(indexii_b(:index_b(3,i),i)),&
                   yterr(indexii_b(:index_b(3,i),i)),&
            weights_all(index_b(1,i):index_b(2,i),1),&
                          dx(i),dy,&
                        index_b(3,i),&
                        theta1(j),hc(i),ol_target(i,j))
        !landfrac may cause coast area par diminish
        !ol_target(i,:)=ol_target(i,:)*landfrac_target(i)
                enddo
        enddo

        print*,"after OLgrid"

        !get correspondent relationship for terr,terrx,terry,wt
        terrout=1.d36
        do i=1,ntarget
        terrout(1,i,:index_b(3,i))= terr(indexii_b(:index_b(3,i),i))
        terrout(2,i,:index_b(3,i))=xterr(indexii_b(:index_b(3,i),i))
        terrout(3,i,:index_b(3,i))=yterr(indexii_b(:index_b(3,i),i))
        terrout(4,i,:index_b(3,i))=weights_all(index_b(1,i):index_b(2,i),1)
        enddo

!#if 0
        do j=1,nvar_dir
        where(abs(sgh_target)<.005_r8) ol_target(:,j)=0.0_r8
        enddo
        where(abs(ol_target)<.001_r8.or.&
              abs(ol_target).gt.1e+7) ol_target=0.0_r8
        where(abs(ol_target).gt.1) ol_target=1.0_r8
        where(ol_target.ne.ol_target) ol_target=0.0_r8
!#endif
end subroutine OLdir
#endif
!=====================
subroutine dxygrid(dx,dy,theta_in,dxy)
IMPLICIT NONE
real(r8),intent(in) :: dx,dy,theta_in
real(r8),intent(out):: dxy
real(r8) :: rad,theta,theta1
                rad=4.0_r8*atan(1.0_r8)/180.0_r8
                theta1=MOD(theta_in,360._r8)
                !set negative axis into 0~360
                if (theta1.ge.-360._r8.and.theta1.lt.0._r8) then
                theta1=theta1+360._r8
                endif   
                !in case the angle is not into the judgement
                theta=theta1
                !transform of angle into first quadrant
                if      (theta1.ge.  0._r8.and.theta1.lt. 90._r8) then
                theta=theta1
                else if (theta1.gt. 90._r8.and.theta1.lt.180._r8) then
                theta=(180._r8-theta1)
                else if (theta1.gt.180._r8.and.theta1.lt.270._r8) then
                theta=(theta1-180._r8)
                else if (theta1.gt.270._r8.and.theta1.lt.360._r8) then
                theta=(360._r8-theta1)
                else if (theta1.eq.90._r8.or.theta1.eq.270._r8) then
                theta=90._r8
                else if (theta1.eq.0._r8.or.theta1.eq.180._r8) then
                theta=0._r8
                endif

                !get dxy
                if   (theta.ge. 0._r8.and.theta.lt.atan2(dy,dx)/rad) then
                dxy=dx/cos(theta*rad)
                else if (theta.ge.atan2(dy,dx)/rad.and.theta.le.90._r8)then
                dxy=dy/sin(theta*rad)
                endif
!print*,"atan2(dy,dx)/rad,theta_in,theta,dy,sin(theta*rad),dxy",atan2(dy,dx)/rad,theta_in,theta,dy,sin(theta*rad),dxy
end subroutine dxygrid
!=======================
!=======================
end module ogwd_sub


