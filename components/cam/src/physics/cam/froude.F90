MODULE froude

  use shr_kind_mod,   only: r8 => shr_kind_r8, SHR_KIND_CL
  use physics_types,  only: physics_state

  implicit none
  private
  save

  public :: calc_uovern

CONTAINS


  subroutine calc_uovern(state, cam_out)

!---- determine the ratio u/n, where u is the mean wind speed and n is the mean brunt-vaisalla frequency over the lowest 200 hpa of the atmosphere
!     call from phys_run1 (module physpkg)
!
!---------------------------Code history--------------------------------
!
! Author:            Steve Ghan, May 2018
! Applied to E3SM    Steve Ghan, Jun 2018
!
!-----------------------------------------------------------------------

    use physconst,      only: rair, gravit, cpair, latvap, cappa, zvir, cpwv, epsilo
    use dimensions_mod, only: nelemd
    use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
    use cam_abortutils, only: endrun
    use cam_logfile,    only: iulog
    use hycoef,         only: hyam, hybm, hyai, hybi, ps0
    use cam_history,    only: outfld
    use wv_saturation,  only: qsat
    use camsrfexch,     only: cam_out_t
!
!     i/o fields
!
    type(physics_state), intent(in) :: state
    type (cam_out_t), intent(inout) :: cam_out

!     local fields

    real(r8) :: uovern(pcols)  ! ratio of wind speed/brunt vaisalla frequency
    real(r8) :: thx(pverp) ! potential temperature
    real(r8) :: t(pver) ! temperature
    real(r8) :: tint(pverp) ! temperature at layer interfaces
    real(r8) :: q(pver) ! specific humidity
    real(r8) :: qint(pverp) ! specific humidity at layer interfaces
    real(r8) :: bvf(pver) ! brunt-vaisala frequency
    real(r8) :: pmid(pver) ! pressure at model levels
    real(r8) :: pint(pverp)   ! Interface pressures
    real(r8) :: lnpint(pverp) ! ln(pressure) at model interfaces
    real(r8) :: pdel(pver) ! pdel(k)   = pint  (k+1)-pint  (k)
    real(r8) :: zs ! surface elevation
    real(r8) :: zint(pverp) ! height at model interfaces
    real(r8) :: zmid(pver)  ! height at model levels
    real(r8) :: bvfmin ! minimum brunt-vaisala frequency
    real(r8) :: wspeed ! wind speed
    real(r8) :: cmass ! column mass
    real(r8) :: wspdbar, bvfbar
    real(r8) :: gam, gammadry, gammawet
    real(r8) :: es,qs,rh(pver)
    integer  :: k, kmax, kmin  ! layer index
    integer  :: i, n,m  ! constituent index
    real(r8) :: rovg  ! rair/gravit
    real(r8) :: d
    integer :: ncol, lchnk

    character(len=*), parameter      :: subname = "froude"

   ncol  = state%ncol


   do i=1,ncol
          do k = 1, pverp
            pint(k) = hyai(k)*ps0 + hybi(k)*state%ps(i)
            lnpint(k) = log(pint(k))
          end do
          do k = 1, pver
            pdel(k) = pint(k+1) - pint(k)
            pmid(k) = 0.5_r8*(pint(k) + pint(k+1))
            t(k)=state%t(i,k)
            q(k)=state%q(i,k,1)
          end do
          kmax=pver-1
          do k=kmax-1,1,-1
             if(pint(kmax)-pint(k).gt.2.e4)then
                kmin=k
                exit
             end if 
          end do

          do k = 1, pver
            call qsat(t(k),pmid(k),es,qs)
            rh(k)=q(k)/qs
          end do
!
!---- determine the brunt-vaisala frequency
!
          do k = 2, pver
            thx(k)=0.5*(state%t(i,k)*(state%ps(i)/pmid(k))**cappa+state%t(i,k-1)*(state%ps(i)/pmid(k-1))**cappa)
            tint(k)=0.5*(state%t(i,k)+state%t(i,k-1))
            qint(k)=0.5*(state%q(i,k,1)+state%q(i,k-1,1))
          end do
!
          thx(pver+1) = state%t(i,pver)
          tint(pver+1) = state%t(i,pver)
          qint(pver+1) = state%q(i,pver,1)
!
          zs = state%phis(i) / gravit
          zint(pver+1) = zs
!
          rovg=rair/gravit
          do k = pver, 1, -1
            zint(k) = zint(k+1)+                             &
                 rovg*state%t(i,k)*(1. + zvir * state%q(i,k,1))*(lnpint(k+1)-lnpint(k))
          end do
!
          do k = 1, pver
            zmid(k) = 0.5_r8*(zint(k)+zint(k+1))
          end do
        !  zmid(pver+1) = zint(pver+1) 
          zmid(pver) = zint(pver+1) !TKT
!
!---- maximum height rise of an air parcel
!
          bvfmin = 1.e5_r8
          do k = kmin,kmax
            if(rh(k).lt.0.9)then
              bvf(k) = gravit*(thx(k)-thx(k+1))/         &
                 (0.5_r8*(thx(k)+thx(k+1))*(zmid(k)-zmid(k+1)))
            else
!             moist bvf (Hughes et al., JAS, 2009)
              call qsat(t(k),pmid(k),es,qs,gam)
              gammadry=gravit/cpair
              gammawet=gammadry*(1+q(k))*(1+latvap*qs/(rair*t(k)))/  &
                 (1+cpwv*qs/cpair+latvap*latvap*epsilo*qs/(cpair*rair*t(k)*t(k))*(1+qs/epsilo))
              bvf(k) = gravit/t(k)*(gammawet+(tint(k+1)-tint(k))/(zint(k+1)-zint(k)))* &
                 (1+latvap*qs/(rair*t(k)))-gravit/(1+q(k))*(qint(k+1)-qint(k))/(zint(k+1)-zint(k))
            endif
            bvfmin=min(bvf(k),bvfmin)

          end do
          if(bvfmin <= 0.0_r8) then
            uovern(i) = 1.e8_r8
          else
            uovern(i) = 0.0_r8
            cmass = 0.0_r8
            wspdbar=0
            bvfbar=0
            do k = kmin,kmax
              wspeed=state%u(i,k)*state%u(i,k)+state%v(i,k)*state%v(i,k)
              wspdbar= wspdbar+ wspeed*pdel(k)
              bvfbar= bvfbar+bvf(k)*pdel(k)
              uovern(i)= uovern(i) + pdel(k)*sqrt(wspeed/bvf(k))
              cmass=cmass+pdel(k)
            end do
            uovern(i) = sqrt(wspdbar/bvfbar)
          end if
          cam_out%uovern(i) = uovern(i)

      end do   ! i

      lchnk=state%lchnk
    
       call outfld('UOVERN', uovern, pcols, lchnk   )



    return
  end subroutine calc_uovern

END MODULE froude


