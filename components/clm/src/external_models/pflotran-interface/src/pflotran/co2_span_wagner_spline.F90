  module co2_span_wagner_spline_module

  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"
  
  save

  private
  
  PetscInt :: nptab,nttab0,ncrit_pts
  PetscInt, allocatable :: nttab(:),ncrit(:)
  PetscReal, allocatable :: p_tab(:),t_tab(:,:),r_tab(:,:), &
                         h_tab(:,:),u_tab(:,:),s_tab(:,:),f_tab(:,:), &
                         rr(:,:),hh(:,:),uu(:,:),ss(:,:),ff(:,:)
  PetscReal, allocatable :: tcrit(:),pcrit(:),rhol(:),ul(:),hl(:),sl(:), &
                         rhov(:),hv(:),uv(:),sv(:),fv(:)
  
  public sw_spline_read, sw_prop
  
  contains
subroutine sw_spline_read

  use spline_module

  PetscInt :: i,ipx,j,n,iunit=9
  
  open (unit=iunit,file='co2_prop_TC.dat',status='old')
  
  read(iunit,*) nptab
  
  allocate(p_tab(nptab))
  allocate(nttab(nptab))
  allocate(ncrit(nptab))

  read(iunit,*) (p_tab(i),i=1,nptab)
  
  read(iunit,*) (nttab(i),i=1,nptab)
  
! read: t, rho, h,u,f,s
  nttab0 = 1
  do i=1,nptab
    nttab0 = max(nttab0,nttab(i))
  enddo
  print *,'nttab0=',nttab0
  
  allocate(t_tab(nptab,nttab0))
  allocate(r_tab(nptab,nttab0))
  allocate(h_tab(nptab,nttab0))
  allocate(u_tab(nptab,nttab0))
  allocate(s_tab(nptab,nttab0))
  allocate(f_tab(nptab,nttab0))
  allocate(rr(nptab,nttab0))
  allocate(hh(nptab,nttab0))
  allocate(uu(nptab,nttab0))
  allocate(ss(nptab,nttab0))
  allocate(ff(nptab,nttab0))

  do i=1,nptab
    read(iunit,*) (t_tab(i,j),j=1,nttab(i))
  enddo
! print *,'finished reading t'
  
  do i=1,nptab
    read(iunit,*) (r_tab(i,j),j=1,nttab(i))
  enddo
! print *,'finished reading r'

  do i=1,nptab
    read(iunit,*) (h_tab(i,j),j=1,nttab(i))
  enddo
! print *,'finished reading h'

  do i=1,nptab
    read(iunit,*) (u_tab(i,j),j=1,nttab(i))
  enddo
! print *,'finished reading u'

  do i=1,nptab
    read(iunit,*) (f_tab(i,j),j=1,nttab(i))
  enddo
! print *,'finished reading f'
  
  read(iunit,*) (ncrit(i),i=1,nptab)
  read(iunit,*) ncrit_pts
! print *,'ncrit_pts=',ncrit_pts
  
  allocate(tcrit(ncrit_pts))
  allocate(pcrit(ncrit_pts))
  allocate(rhol(ncrit_pts))
  allocate(ul(ncrit_pts))
  allocate(hl(ncrit_pts))
  allocate(sl(ncrit_pts))
  allocate(rhov(ncrit_pts))
  allocate(hv(ncrit_pts))
  allocate(uv(ncrit_pts))
  allocate(sv(ncrit_pts))
  allocate(fv(ncrit_pts))
  do i = 1, ncrit_pts
    read(iunit,*) tcrit(i),pcrit(i),rhol(i),hl(i),ul(i),rhov(i),hv(i),uv(i),fv(i)
  enddo

!---  Vapor temperature splines  ---
  rr = 0.d0
  do ipx = 1,nptab
    n = nttab(ipx)-ncrit(ipx)+1
    call spline(t_tab(ipx,ncrit(ipx):nttab(ipx)),r_tab(ipx,ncrit(ipx):nttab(ipx)),n,rr(ipx,ncrit(ipx):nttab(ipx)))
    call spline(t_tab(ipx,ncrit(ipx):nttab(ipx)),h_tab(ipx,ncrit(ipx):nttab(ipx)),n,hh(ipx,ncrit(ipx):nttab(ipx)))
    call spline(t_tab(ipx,ncrit(ipx):nttab(ipx)),u_tab(ipx,ncrit(ipx):nttab(ipx)),n,uu(ipx,ncrit(ipx):nttab(ipx)))
    call spline(t_tab(ipx,ncrit(ipx):nttab(ipx)),f_tab(ipx,ncrit(ipx):nttab(ipx)),n,ff(ipx,ncrit(ipx):nttab(ipx)))
!   call spline(t(ipx,ncrit(ipx)),s(ipx,ncrit(ipx)),n,ss(ipx,ncrit(ipx)))
!   print *,'p= ',n,ipx,p_tab(ipx),rr(ipx,ncrit(ipx)),ncrit(ipx)
!   if (ipx == 10 .or. ipx ==11) then
!     print *,ncrit(ipx), nttab(ipx), rr(ipx,ncrit(ipx):nttab(ipx))
!     print *,r_tab(ipx,ncrit(ipx):nttab(ipx))
!   endif
  enddo

  return
end subroutine sw_spline_read

! ************************************************************************** !

subroutine sw_prop(tx,px,rho,h,u,fg)
      
       use spline_module

!     density of liquid or vapor co2.

!     isrx liquid or vapor index: 1 - liquid 2 - vapor or supercritical

!     span, r., and w. wagner.  1996.  a new equation of state for
!     carbon dioxide covering the fluid region from the triple-point
!     to 1100 k at pressures up to 800 mpa.
!     J. phys. chem. ref. data 25(6):1509-1588.

      implicit none
      save

      PetscReal :: tx,px,rho, h,u,fg
      PetscInt :: ipx,jpx,n
      PetscReal :: tkx,pcx,ptx,tcx,ttx
      PetscReal :: rtab(nptab+1) !,rtab2(nptab+1)
      PetscReal :: htab(nptab+1) !,htab2(nptab+1)
      PetscReal :: utab(nptab+1) !,utab2(nptab+1)
      PetscReal :: fgtab(nptab+1) !,fgtab2(nptab+1)

      tkx = tx + 273.15d0
      pcx = 7.3773d0
      ptx = 0.51795d0
      tcx = 304.1282d0
      ttx = 216.592d0

          jpx = 0
          do ipx = 1,nptab
            n = nttab(ipx)-ncrit(ipx)+1
            if (tkx.gt.t_tab(ipx,1)) then
              jpx = jpx+1
              call splint(t_tab(ipx,ncrit(ipx):nttab(ipx)),r_tab(ipx,ncrit(ipx):nttab(ipx)), &
                rr(ipx,ncrit(ipx):nttab(ipx)),n,tkx,rtab(jpx))
 
              call splint(t_tab(ipx,ncrit(ipx):nttab(ipx)),h_tab(ipx,ncrit(ipx):nttab(ipx)), &
                hh(ipx,ncrit(ipx):nttab(ipx)),n,tkx,htab(jpx))

              call splint(t_tab(ipx,ncrit(ipx):nttab(ipx)),u_tab(ipx,ncrit(ipx):nttab(ipx)), &
                uu(ipx,ncrit(ipx):nttab(ipx)),n,tkx,utab(jpx))

              call splint(t_tab(ipx,ncrit(ipx):nttab(ipx)),f_tab(ipx,ncrit(ipx):nttab(ipx)), &
                ff(ipx,ncrit(ipx):nttab(ipx)),n,tkx,fgtab(jpx))

!             print *,ipx,jpx,t_tab(ipx,ncrit(ipx)),r_tab(ipx,ncrit(ipx)),rtab(jpx)
            endif
          enddo
          !call locate(p_tab,jpx,px,ipx)
          !ipx = min(max(1,ipx),jpx-1)

#if 0
! Density
          call spline(p_tab,rtab,nptab,rtab2)   
          call splint(p_tab,rtab,rtab2,nptab,px,rho)
! H
          call spline(p_tab,htab,nptab,htab2)   
          call splint(p_tab,htab,htab2,nptab,px,h)
! U
          call spline(p_tab,utab,nptab,utab2)   
          call splint(p_tab,utab,utab2,nptab,px,u)
! fg
          call spline(p_tab,fgtab,nptab,fgtab2)   
          call splint(p_tab,fgtab,fgtab2,nptab,px,fg)
#endif

! ************** linear interpolation in pressure *******************
          call locate(p_tab,jpx,px,ipx)
          ipx = min(max(1,ipx),jpx-1)
          rho = (rtab(ipx+1)-rtab(ipx))*(px-p_tab(ipx))/(p_tab(ipx+1)-p_tab(ipx)) + rtab(ipx)
          h   = (htab(ipx+1)-htab(ipx))*(px-p_tab(ipx))/(p_tab(ipx+1)-p_tab(ipx)) + htab(ipx)
          u   = (utab(ipx+1)-utab(ipx))*(px-p_tab(ipx))/(p_tab(ipx+1)-p_tab(ipx)) + utab(ipx)
          fg  = (fgtab(ipx+1)-fgtab(ipx))*(px-p_tab(ipx))/(p_tab(ipx+1)-p_tab(ipx)) + fgtab(ipx)

!         print *,'density: ',tx,px,ipx,jpx,p_tab(ipx+1),p_tab(ipx),rtab(ipx+1),rtab(ipx),rho
!************************************************************************************
    return
end subroutine sw_prop
    
end module co2_span_wagner_spline_module
