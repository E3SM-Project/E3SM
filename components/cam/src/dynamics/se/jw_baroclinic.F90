#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!  This file contains the initial condititions for the baroclinic
!  instability probelm:
!
!  Jablonowski and Williamson, QJR (2006) 132 
!
!  SZ and HW 2018-12: Currently, no tracers are initialized. 
!
module jw_baroclinic
!
!  This module contains the initial condititions for the baroclinic
!  instability probelms:
!
!  Jablonowski and Williamson, QJR (2006) 132 

 ! ====================
  use cam_logfile,    only : iulog
  use element_mod,    only : element_t
  use kinds,          only : real_kind, int_kind
  use control_mod,    only : u_perturb
  use hycoef,         only : ps0, hyam, hybm

  implicit none
  private
  public :: init_jw_baroclinic


contains

  subroutine init_jw_baroclinic(elem,tl)

    use dimensions_mod,      only : nelemd, nlev, np, npsq
    use physconst,           only : omega, DD_PI=>pi, rearth, rgas =>rair, Cp =>cpair, g=>gravit

    implicit none

    type(element_t), intent(inout) :: elem(:)

    integer :: i, j, ie, k
    integer :: tl  ! time level to initialize

!=======================================================================================================!
    real (kind=real_kind)  :: tbar(nlev)
    real (kind=real_kind) ::  eta(nlev), etv(nlev)
    real (kind=real_kind) ::  latc, lonc, rr,rc,  aa, lat,lon, snlat, cslat , v1,v2
    real (kind=real_kind) ::  trm1,trm2,trm3,trm4, term
!=======================================================================================================!
    real (kind=real_kind),  parameter  :: u0        = 35.0_real_kind      ! Zonal Mean wind 
    real (kind=real_kind),  parameter  :: t0        = 288.0_real_kind     ! Mean temp
    real (kind=real_kind),  parameter  :: gama      = 0.005_real_kind     ! Lapse rate
    real (kind=real_kind),  parameter  :: ddt       = 4.8e05_real_kind    ! Temp  gradient 
    real (kind=real_kind),  parameter  :: eta_t     = 0.2_real_kind       ! eta level 
    real (kind=real_kind),  parameter  :: eta_s     = 1.0_real_kind       ! eta level 
    real (kind=real_kind),  parameter  :: eta_0     = 0.252_real_kind     ! eta level 
    real (kind=real_kind),  parameter  :: u_perturb = 1.0_real_kind       ! wind perturbation

    real (kind=real_kind) :: r_d,omg,grv,erad
    real (kind=real_kind) :: pmin, pmax

    pmin = ps0
    pmax = 0._real_kind

    r_d = Rgas
    omg=omega
    erad=rearth
    grv=g

    do k = 1, nlev
       eta(k) = hyam(k)+hybm(k)
       etv(k) = (eta(k) - eta_0) * DD_PI * 0.5_real_kind
    end do

    latc = DD_PI * (2.0_real_kind/9.0_real_kind)
    lonc = DD_PI * (1.0_real_kind/9.0_real_kind)

    do ie=1,nelemd
        ! initial velocity   
        elem(ie)%state%v = 0.0_real_kind
        do k=1, nlev
           do j=1,np
              do i=1,np
              lon = elem(ie)%spherep(i,j)%lon
              lat = elem(ie)%spherep(i,j)%lat
              snlat=SIN(lat)
              cslat=COS(lat)

              aa = SIN(latc)*snlat + COS(latc)*cslat*COS(lon - lonc)
              rc =  10.0_real_kind  * ACOS(aa)
              v1 = u0 * (cos(etv(k)))**1.5_real_kind * (sin(2.0_real_kind * lat))**2  +  u_perturb*exp(-rc*rc)
              elem(ie)%state%v(i,j,1,k,tl)=v1
              elem(ie)%state%v(i,j,2,k,tl)=0
              end do
           end do
        end do


        ! Layer-mean Temperature fields 
        elem(ie)%state%T = 0.0_real_kind
        do k = 1, nlev
         if (eta(k) <= eta_t) then
            tbar(k) = t0 * eta(k)**(r_d*gama/grv) + ddt * (eta_t - eta(k))**5
         else
            tbar(k) = t0 * eta(k)**(r_d*gama/grv)
         endif
        end do

        do k=1,nlev
           do j=1,np
             do i=1,np
             lon = elem(ie)%spherep(i,j)%lon
             lat = elem(ie)%spherep(i,j)%lat

             snlat=SIN(lat)
             cslat=COS(lat)

             trm1 = 0.75_real_kind * (eta(k) * DD_PI*u0 /r_d) * sin(etv(k)) *sqrt(cos(etv(k)))
             trm2 = -2.0_real_kind *snlat**6 *(cslat**2 + 1.0_real_kind/3.0_real_kind) + 10.0_real_kind/63.0_real_kind
             trm3 =  2.0_real_kind * u0* (cos(etv(k)))**1.5_real_kind
             trm4 = (1.60_real_kind *cslat**3 *(snlat**2 + 2.0_real_kind/3.0_real_kind) - DD_PI *0.25_real_kind)* erad*omg
        
             elem(ie)%state%T(i,j,k,tl) = tbar(k) + trm1 *(trm2 * trm3 + trm4 )

             end do
          end do
        end do

        !Surface geopotential
        elem(ie)%state%phis = 0.0_real_kind
        elem(ie)%state%ps_v = 0.0_real_kind
        do j=1,np
           do i=1,np
              lon = elem(ie)%spherep(i,j)%lon
              lat = elem(ie)%spherep(i,j)%lat
              snlat=SIN(lat)
              cslat=COS(lat)

              trm1 =  u0* ( cos((eta_s - eta_0)*DD_PI*0.5_real_kind) )**1.5_real_kind

              trm2 = -2.0_real_kind *snlat**6 *(cslat**2 + 1.0_real_kind/3.0_real_kind) + 10.0_real_kind/63.0_real_kind
              trm3 = (1.60_real_kind *cslat**3 *(snlat**2 + 2.0_real_kind/3.0_real_kind) - DD_PI *0.25_real_kind)* erad*omg

              elem(ie)%state%phis(i,j)      = trm1 *(trm2 * trm1 + trm3)
              elem(ie)%state%ps_v(i,j,tl)   = ps0
           end do
        end do
        pmin = min(minval(elem(ie)%state%ps_v(:,:,tl)),pmin)
        pmax = max(maxval(elem(ie)%state%ps_v(:,:,tl)),pmin)

    enddo !element loop

    write(iulog,'(a,e20.8,e20.8)')  "jw_baroclinic initiaization: Ps MinMax = ", pmin, pmax

!=======================================================================================================!
  end subroutine init_jw_baroclinic

end module jw_baroclinic
