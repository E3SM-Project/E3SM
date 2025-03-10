Module orographic_shape_methods
use shr_kind_mod, only: r8 => shr_kind_r8


! updated methods following Xie et al. (2020)
public :: orographic_asymmetry_xie2020
public :: orographic_efflength_xie2020

! original methods following Kim and Doyle (2005)
public :: orographic_convexity_kim2005
public :: orographic_asymmetry_kim2005
public :: orographic_efflength_kim2005

! utility functions
private :: orographic_efflength_grid
private :: dxygrid

!===================================================================================================
! References

! Xie, J., Zhang, M., Xie, Z., Liu, H., Chai, Z., He, J. X., & Zhang, H. (2020).
! An orographic-drag parametrization scheme including orographic anisotropy for
! all flow directions. Journal of Advances in Modeling Earth Systems, 12, e2019MS001921.

! Kim, Y.-J., & Doyle, J. D. (2005). Extension of an orographic-drag parameterization 
! scheme to incorporate orographic anisotropy and flow blocking. Quarterly Journal 
! of the Royal Meteorological Society, 131, 1893–1921.

!===================================================================================================
contains
!===================================================================================================
subroutine orographic_asymmetry_xie2020( terr, ntarget, ncube, n, nvar_dir, jall, &
                                         weights_lgr_index_all, weights_eul_index_all1, &
                                         weights_eul_index_all2, weights_eul_index_all3, &
                                         weights_all, lon_cen, lat_cen, lon_terr, lat_terr, &
                                         area_target, landfrac_target, oa_target)
  IMPLICIT NONE
  integer ,intent(in)    :: ncube,ntarget
  integer ,intent(in)    :: n
  integer ,intent(in)    :: nvar_dir
  integer ,intent(in)    :: jall
  integer ,intent(in)    :: weights_lgr_index_all(jall)
  integer ,intent(in)    :: weights_eul_index_all1(jall)
  integer ,intent(in)    :: weights_eul_index_all2(jall)
  integer ,intent(in)    :: weights_eul_index_all3(jall)
  real(r8),intent(in)    :: weights_all(jall,1)
  real(r8),intent(in)    :: terr(n)
  real(r8),intent(inout) :: lon_cen(ntarget)
  real(r8),intent(inout) :: lat_cen(ntarget)
  real(r8),intent(inout) :: area_target(ntarget)
  real(r8),intent(in)    :: lon_terr(n)
  real(r8),intent(in)    :: lat_terr(n)
  real(r8),intent(in)    :: landfrac_target(ntarget)
  real(r8),intent(out)   :: oa_target(ntarget,nvar_dir)
  !local
  integer  :: count,i,ix,iy,ip,ii,ip2,ip3
  real(r8) :: xxterr
  real(r8) :: yyterr
  real(r8) :: zzterr
  real(r8) :: ix2
  real(r8) :: iy2
  real(r8) :: xx3
  real(r8) :: yy3
  real(r8) :: zz3
  real(r8) :: ix3
  real(r8) :: iy3
  real(r8) :: wt,xhds(ntarget)
  real(r8) :: yhds(ntarget)
  real(r8) :: zhds(ntarget)
  real(r8) :: hds(ntarget)
  real(r8) :: OAx_var(ntarget)
  real(r8) :: OAy_var(ntarget)
  real(r8) :: OAz_var(ntarget)
  real(r8) :: OA_var(ntarget)
  real(r8) :: xbar(ntarget)
  real(r8) :: ybar(ntarget)
  real(r8) :: zbar(ntarget)
  real(r8) :: lon_bar(ntarget)
  real(r8) :: lat_bar(ntarget)
  real(r8) :: rad
  real(r8) :: theta1
  real(r8) :: OAlon(ntarget)
  real(r8) :: OAlat(ntarget)
  real(r8) :: OArad(ntarget)
  real(r8) :: OAx1
  real(r8) :: OAy1
  real(r8) :: OAz1

  real(r8) :: terr_anom(n)
  real(r8) :: terr_avg(ntarget)
  real(r8) :: weights_ano(jall)
  real(r8) :: area_target_ano(ntarget)

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
  OAx_var=0.0_r8
  OAy_var=0.0_r8
  OAz_var=0.0_r8

  terr_anom=0.0_r8
  terr_avg=0.0_r8
  do count=1,jall
    i   = weights_lgr_index_all(count)
    ix  = weights_eul_index_all1(count)!,1)
    iy  = weights_eul_index_all2(count)!,2)
    ip  = weights_eul_index_all3(count)
    ! convert to 1D indexing of cubed-sphere
    ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!
    wt = weights_all(count,1)
    !
    terr_avg(i)=terr_avg(i)+(wt/area_target(i))*terr(ii)
    !terr(ii)*wt!(wt/area_target(i))*terr(ii)
  enddo

  do count=1,jall
    i   = weights_lgr_index_all(count)
    ix  = weights_eul_index_all1(count)!,1)
    iy  = weights_eul_index_all2(count)!,2)
    ip  = weights_eul_index_all3(count)
    ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix
    terr_anom(ii)=terr(ii)-terr_avg(i)
  enddo
  where(terr_anom.le.0) terr_anom=0.0_r8

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

    xhds(i)=xhds(i)+xxterr*terr_anom(ii)*wt
    yhds(i)=yhds(i)+yyterr*terr_anom(ii)*wt
    zhds(i)=zhds(i)+zzterr*terr_anom(ii)*wt
    hds(i) =hds(i)+terr_anom(ii)*wt

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
    OAx_var(i)=OAx_var(i)+(wt/area_target(i))*(xxterr-xx3)**2
    OAy_var(i)=OAy_var(i)+(wt/area_target(i))*(yyterr-yy3)**2
    OAz_var(i)=OAz_var(i)+(wt/area_target(i))*(zzterr-zz3)**2
    OAx1=(xx3-xbar(i))/sqrt(OAx_var(i))!OA_var(i))
    OAy1=(yy3-ybar(i))/sqrt(OAy_var(i))!OA_var(i))
    OAz1=(zz3-zbar(i))/sqrt(OAz_var(i))!OA_var(i))

    !all lat_cen must use (90-lat_cen) since we only have 
    !latitude rather than colatitude
    !this is equivalent to using induction formula sin(90-lat)=cos(lat)
    !latitude is opposite of colatitude, thus OAlat is negative
    OAlat(i)=-(OAx1*sin(lat_cen(i)*rad)*cos(lon_cen(i)*rad)&
              +OAy1*sin(lat_cen(i)*rad)*sin(lon_cen(i)*rad)&
              -OAz1*cos(lat_cen(i)*rad))
    OAlon(i)= -OAx1*sin(lon_cen(i)*rad)&
              +OAy1*cos(lon_cen(i)*rad)

    !reverse in order to be (2,ntarget),OAx,OAy
    oa_target(i,1) = OAlon(i)
    oa_target(i,2) = OAlat(i)

  enddo

  !takeout abnormal values
  where(abs(oa_target)<.001_r8.or.abs(oa_target).gt.1e+7) oa_target=0.0_r8
  !where(abs(oa_target).gt.1) oa_target=1.0_r8
  where(oa_target.ne.oa_target) oa_target=0.0_r8

end subroutine orographic_asymmetry_xie2020
!===================================================================================================
subroutine orographic_asymmetry_kim2005( terr, ntarget, ncube, n, jall, &
                                         weights_lgr_index_all, weights_eul_index_all1, &
                                         weights_eul_index_all2, weights_eul_index_all3, &
                                         weights_all, lon_terr, lat_terr, area_target, &
                                         landfrac_target, oa_target)
  IMPLICIT NONE
  integer ,intent(in)  :: ncube
  integer ,intent(in)  :: ntarget
  integer ,intent(in)  :: n
  integer ,intent(in)  :: jall
  integer ,intent(in)  :: weights_lgr_index_all(jall)
  integer ,intent(in)  :: weights_eul_index_all1(jall)
  integer ,intent(in)  :: weights_eul_index_all2(jall)
  integer ,intent(in)  :: weights_eul_index_all3(jall)
  real(r8),intent(in)  :: weights_all(jall,1)
  real(r8),intent(in)  :: terr(n)
  real(r8),intent(in)  :: lon_terr(n)
  real(r8),intent(in)  :: lat_terr(n)
  real(r8),intent(in)  :: area_target(ntarget)
  real(r8),intent(in)  :: landfrac_target(ntarget)
  real(r8),intent(out) :: oa_target(ntarget,4)
  !local
  real(r8) :: xh(ntarget)
  real(r8) :: yh(ntarget)
  real(r8) :: height(ntarget)
  real(r8) :: modexcoords(ntarget)
  real(r8) :: modeycoords(ntarget)
  real(r8) :: avgx(ntarget)
  real(r8) :: avgy(ntarget)
  real(r8) :: varx(ntarget)
  real(r8) :: vary(ntarget)
  real(r8) :: OAx
  real(r8) :: OAy
  real(r8) :: theta1
  real(r8) :: rad
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
    ! apply landfrac
    oa_target(i,:)=  oa_target(i,:)*landfrac_target(i)
  enddo
  !takeout abnormal values
  where(abs(oa_target)<.001_r8.or.abs(oa_target).gt.1e+7) oa_target=0.0
  where(abs(oa_target).gt.1) oa_target=0.0
  where(oa_target.ne.oa_target) oa_target=0.0
end subroutine orographic_asymmetry_kim2005
!===================================================================================================
subroutine orographic_convexity_kim2005( terr, ntarget, ncube, n, jall, &
                                         weights_lgr_index_all, weights_eul_index_all1, &
                                         weights_eul_index_all2, weights_eul_index_all3, &
                                         weights_all, area_target, sgh_target, terr_target, &
                                         landfrac_target, oc_target)
  IMPLICIT NONE
  integer ,intent(in)  :: ncube
  integer ,intent(in)  :: ntarget
  integer ,intent(in)  :: n
  integer ,intent(in)  :: jall
  integer ,intent(in)  :: weights_lgr_index_all(jall)
  integer ,intent(in)  :: weights_eul_index_all1(jall)
  integer ,intent(in)  :: weights_eul_index_all2(jall)
  integer ,intent(in)  :: weights_eul_index_all3(jall)
  real(r8),intent(in)  :: weights_all(jall,1)
  real(r8),intent(in)  :: area_target(ntarget)
  real(r8),intent(in)  :: sgh_target(ntarget)
  real(r8),intent(in)  :: terr_target(ntarget)
  real(r8),intent(in)  :: terr(n)
  real(r8),intent(in)  :: landfrac_target(ntarget)
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
    ! apply landfrac
    oc_target(i) = oc_target(i) * landfrac_target(i)
  enddo

  where(abs(oc_target)<.001_r8.or.abs(oc_target).gt.1e+7) oc_target=0.0_r8
  where(abs(sgh_target).eq.0.0_r8) oc_target=0.0_r8
  where(oc_target<0.0_r8) oc_target=0.0_r8
end subroutine orographic_convexity_kim2005
!===================================================================================================
subroutine orographic_efflength_kim2005( terr, ntarget, ncube, n, jall, &
                                         weights_lgr_index_all, weights_eul_index_all1, &
                                         weights_eul_index_all2, weights_eul_index_all3, &
                                         weights_all, lon_terr, lat_terr, area_target, sgh_target, &
                                         target_center_lat, target_center_lon, &
                                         target_corner_lat, target_corner_lon, &
                                         landfrac_target, ol_target)
  IMPLICIT NONE
  integer, intent(in)  :: ncube
  integer, intent(in)  :: ntarget
  integer, intent(in)  :: n
  integer, intent(in)  :: jall
  integer, intent(in)  :: weights_lgr_index_all(jall)
  integer, intent(in)  :: weights_eul_index_all1(jall)
  integer, intent(in)  :: weights_eul_index_all2(jall)
  integer, intent(in)  :: weights_eul_index_all3(jall)
  real(r8),intent(in)  :: weights_all(jall,1)
  real(r8),intent(in)  :: area_target(ntarget)
  real(r8),intent(in)  :: sgh_target(ntarget)
  real(r8),intent(in)  :: terr(n)
  real(r8),intent(in)  :: lon_terr(n)
  real(r8),intent(in)  :: lat_terr(n)
  real(r8),intent(in)  :: target_center_lat(ntarget)
  real(r8),intent(in)  :: target_center_lon(ntarget)
  real(r8),intent(in)  :: target_corner_lat(4,ntarget)
  real(r8),intent(in)  :: target_corner_lon(4,ntarget)
  real(r8),intent(in)  :: landfrac_target(ntarget)
  real(r8),intent(out) :: ol_target(ntarget,4)
  !local 
  integer  :: count,i,ix,iy,ip,ii,j
  real(r8) :: wt
  real(r8) :: terr_if
  real(r8) :: Nw(4,ntarget)
  real(r8) :: area_target_par(4,ntarget)

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
    .GT.(target_corner_lat(1,i)+target_center_lat(i))/2..and. &
    lat_terr(ii) &!(ii)&
    .LT.(target_corner_lat(4,i)+target_center_lat(i))/2.) then
      Nw(1,i)=Nw(1,i)+wt*terr_if
      area_target_par(1,i)=area_target_par(1,i)+wt
    endif
          
    !OL2
    if (lon_terr(ii) &!(ii)&
    .GT.(target_corner_lon(1,i)+target_center_lon(i))/2..and. &
      lon_terr(ii) &!(ii)&
    .LT.(target_corner_lon(3,i)+target_center_lon(i))/2.) then
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

    do j=1,4
      ol_target(i,j)=Nw(j,i)/(area_target_par(j,i)+1e-14)!Nt(i)!/2.)
    enddo
    ol_target(i,:)=ol_target(i,:)*landfrac_target(i)

  end do
  where(abs(ol_target)<.001_r8.or.abs(ol_target).gt.1e+7) ol_target=0.0_r8
end subroutine orographic_efflength_kim2005
!===================================================================================================
subroutine orographic_efflength_grid( terr, terrx, terry, wt, b, a, n, theta_in, hc, OLout)
  IMPLICIT NONE
  integer, intent(in)  :: n
  real(r8),intent(in)  :: hc
  real(r8),intent(in)  :: wt(n)
  real(r8),intent(in)  :: terr(n)
  real(r8),intent(in)  :: a
  real(r8),intent(in)  :: b
  real(r8),intent(in)  :: theta_in
  real(r8),intent(in)  :: terrx(n),terry(n)
  real(r8),intent(out) :: OLout
  real(r8) :: theta
  real(r8) :: theta1
  real(r8) :: theta2
  real(r8) :: rad
  real(r8) :: interval
  real(r8) :: terr_count(n)
  real(r8) :: terr_whole_count(n)
  real(r8) :: cx(n),c1,c2
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

  if      (theta1.ge.  0._r8.and.theta1.lt.atan2(a,2._r8*b)/rad.or.&
           theta2.ge.  0._r8.and.theta2.lt.atan2(a,2._r8*b)/rad.or.&
           theta1.ge.180._r8-atan2(a,2_r8*b)/rad.and.theta1.lt.180._r8.or.&
           theta2.ge.180._r8-atan2(a,2_r8*b)/rad.and.theta2.lt.180._r8) &
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
           theta2.lt.180._r8-atan2(a,2_r8*b)/rad) &
  then
    interval=2
    c1=-0.5_r8*b*sin(theta)-0.5_r8*a*cos(theta)+sqrt(a*b*sin(2_r8*theta)/4._r8)
    c2= 0.5_r8*b*sin(theta)+0.5_r8*a*cos(theta)-sqrt(a*b*sin(2_r8*theta)/4._r8)
  else if (theta1.ge.atan2(2_r8*a,b)/rad.and.theta1.lt.90._r8.or.&
           theta2.ge.atan2(2_r8*a,b)/rad.and.theta2.lt.90._r8.or.&
           theta1.ge.90._r8.and.theta1.lt.180._r8-atan2(2_r8*a,b)/rad.or.&
           theta2.ge.90._r8.and.theta2.lt.180._r8-atan2(2_r8*a,b)/rad)&
  then
    interval=3
    c1=-0.25_r8*b*sin(theta)
    c2= 0.25_r8*b*sin(theta)
  endif
  !determine two line functions
  cx=terrx*sin(theta_in*rad)-terry*cos(theta_in*rad)

  !assuming rectangle grid or ladder-shape would be pretty similar since in 1.4,
  !the max difference is of 0.02 although the above expression is actually
  !ladder-shape and use a rectangle with the same area to derive c1 and c2
  where ( cx.ge.min(c1,c2) .and. cx.le.max(c1,c2)                ) terr_whole_count=1._r8
  where ( cx.ge.min(c1,c2) .and. cx.le.max(c1,c2) .and.terr.ge.hc) terr_count      =1._r8

  !deals with noise that affects OL when there are no terrx center points
  !in the two lines in the center
  if (sum(wt*terr_whole_count).eq.0._r8) then
    !enlarge about 5 times interval
    !to include more points and avoid 
    !noise
    j=5._r8
    where ((cx+j*abs(c1-c2)).ge.min(c1,c2) .and. (cx-j*abs(c1-c2)).le.max(c1,c2)) &
    terr_whole_count=1._r8
    where ((cx+j*abs(c1-c2)).ge.min(c1,c2) .and. (cx-j*abs(c1-c2)).le.max(c1,c2) .and.terr.ge.hc) &                      
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
  !take out NaN
  if (OLout.ne.OLout) OLout=0.0_r8
end subroutine orographic_efflength_grid
!===================================================================================================
subroutine orographic_efflength_xie2020( terr, ntarget, ncube, n, jall, nlon, nlat, indexb, nvar_dir, &
                                         weights_lgr_index_all, weights_eul_index_all1, &
                                         weights_eul_index_all2, weights_eul_index_all3, &
                                         weights_all, lon_cen, lat_cen, lon_cor, lat_cor, &
                                         lon_terr, lat_terr, sgh_target, area_target, landfrac_target, &
                                         ol_target, terrout, dxy)
  IMPLICIT NONE
  integer ,intent(in)  :: ncube
  integer ,intent(in)  :: ntarget
  integer ,intent(in)  :: n
  integer ,intent(in)  :: jall
  integer ,intent(in)  :: indexb
  integer ,intent(in)  :: weights_lgr_index_all(jall)
  integer ,intent(in)  :: weights_eul_index_all1(jall)
  integer ,intent(in)  :: weights_eul_index_all2(jall)
  integer ,intent(in)  :: weights_eul_index_all3(jall)
  integer ,intent(in)  :: nlon
  integer ,intent(in)  :: nlat
  integer ,intent(in)  :: nvar_dir
  real(r8),intent(in)  :: weights_all(jall,1)
  real(r8),intent(in)  :: terr(n)
  real(r8),intent(in)  :: lon_terr(n)
  real(r8),intent(in)  :: lat_terr(n)
  real(r8),intent(in)  :: lon_cen(ntarget)
  real(r8),intent(in)  :: lat_cen(ntarget)
  real(r8),intent(in)  :: lon_cor(4,ntarget)
  real(r8),intent(in)  :: lat_cor(4,ntarget)
  real(r8),intent(in)  :: sgh_target(ntarget)
  real(r8),intent(in)  :: area_target(ntarget)
  real(r8),intent(in)  :: landfrac_target(ntarget)
  real(r8),intent(out) :: ol_target(ntarget,nvar_dir)
  real(r8),intent(out) :: terrout(4,ntarget,indexb)
  real(r8),intent(out) :: dxy(ntarget,nvar_dir)

  !local variables
  ! Note: 1 is lower bound,2 is upper bound
  REAL(r8), PARAMETER :: pi = 3.14159265358979323846264338327
  integer, allocatable :: indexii_b(:,:)
  integer  :: index_b(3,ntarget),index_jall(jall)
  integer  :: ix,iy,ip,i,count,alloc_error,j
  real(r8) :: xterr(n),yterr(n),dx(ntarget),dy(ntarget),hc(ntarget),theta1(nvar_dir)
  real(r8) :: xterr_cen(ntarget),yterr_cen(ntarget),rad
  real(r8) :: reflon_terr(n),reflat_terr(n)!,lon_terr2(n)

  real(r8) :: terr_avg(ntarget),wt
  integer  :: ii

  logical :: verbose = .true. ! verbosity flag for debugging

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
  if (verbose) print*,"  oro_efflen: begin ol index_b"
  do count=1,jall
    i   = weights_lgr_index_all(count)
    index_b(3,i)=index_b(3,i)+1
  enddo
  !cumulative sum to form upper and lower bound of index_b
  !1 for lower bound, 2 for upper bound
  do i=1,ntarget
    index_b(2,i)=sum(index_b(3,1:i))
  enddo
  index_b(1,1)=1
  do i=2,ntarget
    index_b(1,i)=index_b(2,i-1)+1
  enddo
  if (verbose) print*,"  oro_efflen: after index_b"
  !leave largest possible dimension first
  allocate(indexii_b(maxval(index_b(3,:)),ntarget),stat=alloc_error)
  if (verbose) print*,"  oro_efflen: maxval(index_b(3,:)",maxval(index_b(3,:))
  indexii_b=0

  do i=1,ntarget               
    do j=1,index_b(3,i)  
      index_jall(index_b(1,i)+j-1)=j                             
    enddo                
  enddo
  if (verbose) print*,"  oro_efflen: after index_jall"

  !get terr avg
  terr_avg=0.0_r8
  do count=1,jall
    i   = weights_lgr_index_all(count)
    ix  = weights_eul_index_all1(count)!,1)
    iy  = weights_eul_index_all2(count)!,2)
    ip  = weights_eul_index_all3(count)
    ! convert to 1D indexing of cubed-sphere
    ii = (ip-1)*ncube*ncube+(iy-1)*ncube+ix!
    wt = weights_all(count,1)
    terr_avg(i)=terr_avg(i)+(wt/area_target(i))*terr(ii)
  end do

  !get correspondent ii for ub and lb                                
  do count=1,jall
    i   = weights_lgr_index_all(count)
    ix  = weights_eul_index_all1(count)!,1)
    iy  = weights_eul_index_all2(count)!,2)
    ip  = weights_eul_index_all3(count)!,3)
    ! convert to 1D indexing of cubed-sphere
    indexii_b(index_jall(count),i) = (ip-1)*ncube*ncube+(iy-1)*ncube+ix
  enddo
  if (verbose) print*,"  oro_efflen: convert indexii_b"
  !input small grids and make OL for individual large grid
  !critical height using empirical function
  !set to avg of the height
  hc=terr_avg
  !hc(i)=1116.2_r8-0.878_r8*sgh_target(i)
  !hc=1116.2_r8-0.878_r8*sgh_target
  !get terrx terry for the small grids
  !transform lon lat grid to distance to (0,0)
  !in case the grid spans 0 line
  !lon_terr2=lon_terr
  !where(lon_terr2.gt.180._r8) lon_terr2=lon_terr2-360._r8

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
  !
  if (verbose) print*,"  oro_efflen: get reflonlat"
  !determine real length on cartesian coordinate
  xterr=reflon_terr*cos(lat_terr*rad)
  yterr=reflat_terr
  !get dx dy for the large grid
  !dx is a function of latitude
  !at this time, we do not need real length
  !so R is set to normalized 1
  !dx=rearth*cos(rlat)*(2._r8*pi/256._r8),dy=rearth*(pi/(128._r8-1._r8))
  !dx=(2._r8*pi/real(nlon,kind=r8))*cos(lat_cen*rad)
  !dy=(      pi/real(nlat-1,kind=r8))
  !! debug, get dx and dy for each grid
  ! (1):  the lower left corner
  ! (2):  the lower right corner
  ! (3):  the upper right corner
  ! (4):  the upper left corner
  !take half of upper and lower to approximate rectangular
  dx=0.5*(abs(lon_cor(3,:)-lon_cor(4,:))+abs(lon_cor(2,:)-lon_cor(1,:)))*cos(lat_cen*rad)
  dy=0.5*(abs(lat_cor(3,:)-lat_cor(2,:))+abs(lat_cor(4,:)-lat_cor(1,:)))
  !test 4 direction for the time
  !only needs 0-180 half of the axis
  do j=1,ntarget
    do i=1,nvar_dir
      theta1(i)=(180._r8/nvar_dir)*real(i-1,kind=r8)
      call dxygrid(6.37100e6_r8*dx(j),6.37100e6_r8*dy(i),theta1(i),dxy(j,i))
    enddo
  enddo

  if (verbose) print*,"  oro_efflen: before into orographic_efflength_grid"
  !input into every large grid
  do j=1,nvar_dir
    do i=1,ntarget
      call orographic_efflength_grid( terr(indexii_b(:index_b(3,i),i)), &
                                      xterr(indexii_b(:index_b(3,i),i)), &
                                      yterr(indexii_b(:index_b(3,i),i)), &
                                      weights_all(index_b(1,i):index_b(2,i),1), &
                                      dx(i), dy(i), index_b(3,i), &
                                      theta1(j), hc(i), ol_target(i,j) )
    enddo
  enddo
  if (verbose) print*,"  oro_efflen: after orographic_efflength_grid"

  !get correspondent relationship for terr,terrx,terry,wt
  terrout=1.d36
  do i=1,ntarget
    terrout(1,i,:index_b(3,i))= terr(indexii_b(:index_b(3,i),i))
    terrout(2,i,:index_b(3,i))=xterr(indexii_b(:index_b(3,i),i))
    terrout(3,i,:index_b(3,i))=yterr(indexii_b(:index_b(3,i),i))
    terrout(4,i,:index_b(3,i))=weights_all(index_b(1,i):index_b(2,i),1)
  enddo

  do j=1,nvar_dir
    where(abs(sgh_target)<.005_r8) ol_target(:,j)=0.0_r8
  enddo
  where(abs(ol_target)<.001_r8.or.abs(ol_target).gt.1e+7) ol_target=0.0_r8
  where(abs(ol_target).gt.1) ol_target=1.0_r8
  where(ol_target.lt.0) ol_target=0.0_r8
  where(ol_target.ne.ol_target) ol_target=0.0_r8

  if (verbose) print*,"  oro_efflen: done"

end subroutine orographic_efflength_xie2020
!===================================================================================================
! Utility function to determine effective grid length
subroutine dxygrid( dx, dy, theta_in, dxy)
  IMPLICIT NONE
  real(r8),intent(in) :: dx
  real(r8),intent(in) :: dy
  real(r8),intent(in) :: theta_in
  real(r8),intent(out):: dxy
  real(r8) :: rad
  real(r8) :: theta1
  real(r8) :: theta2
  rad = 4.0_r8*atan(1.0_r8)/180.0_r8
  theta1 = MOD(theta_in,360._r8)
  !set negative axis into 0~360
  if (theta1.ge.-360._r8 .and. theta1.lt.0._r8) theta1 = theta1 + 360._r8
  !transform of angle into first quadrant
  if (theta1.ge.  0._r8 .and. theta1.lt. 90._r8) theta2 = theta1
  if (theta1.gt. 90._r8 .and. theta1.lt.180._r8) theta2 = 180._r8 - theta1
  if (theta1.gt.180._r8 .and. theta1.lt.270._r8) theta2 = theta1 - 180._r8
  if (theta1.gt.270._r8 .and. theta1.lt.360._r8) theta2 = 360._r8 - theta1
  !get dxy
  if (theta2.ge. 0._r8 .and. theta2.lt.atan2(dy,dx)/rad) dxy = dx/cos(theta2*rad)
  if (theta2.le.90._r8 .and. theta2.ge.atan2(dy,dx)/rad) dxy = dy/sin(theta2*rad)

end subroutine dxygrid
!===================================================================================================
end module orographic_shape_methods
