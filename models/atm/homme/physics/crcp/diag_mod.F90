module diag_mod
  use grid_init_mod, only : dt, dx, dz, dti, dxi, dzi 
  use gcrk_mod, only : rhsdiv_1 , eer, eem, niter, nitsm, icount
  implicit none
contains
  subroutine diagno_1(ux,uz,th,scr1,scr2,rho, nx, nz, time)
    integer, intent(in) :: nx, nz
    real, intent(in) :: time
    real :: ux (nx, nz), uz (nx, nz), th (nx, nz), scr1 (nx, nz),   &
         scr2 (nx, nz)                                                     
    real :: rho (nx, nz) 
    integer :: i, k, nitav
    real :: amn, amx, cour
!cc pressure solver diagnostics
      integer :: nxz

      nxz=nx*nz
      do i=1,nx
      do k=1,nz
      scr2(i,k)=rho(i,k)
      enddo
      enddo

      print 200, time
 200  format(1x,' ****** analysis for time (min): ',f8.2)

      call minmax(ux,nxz,amn,amx)
      print 201,amn,amx
 201  format(1x,' --> min, max ux: ',2e12.4)

      call minmax(uz,nxz,amn,amx)
      print 202,amn,amx
 202  format(1x,' --> min, max uz: ',2e12.4)

      cour=0.
      do i=1,nxz
      cour=max(cour,abs(ux(i,1))*dt/dx+abs(uz(i,1))*dt/dz)
      enddo
      print 302,cour
 302  format(1x,' --> max courant: ',e12.4)

      call minmax(th,nxz,amn,amx)
      print 203,amn,amx
 203  format(1x,' --> min, max th: ',2e12.4)
      call rhsdiv_1(ux,uz,scr2,scr1,1, nx, nz)

      call minmax(scr1,nxz,amn,amx)
      print 204,amn,amx
 204  format(1x,' --> min, max div: ',2e12.4)

      nitav=nitsm/max(icount,1)
      print 205, eer,eem,niter,nitav
  205 format(1x,'            eer,eem:',2e11.4/ &
           1x,'       niter, nitav:',2i4)

       if(cour.gt.1.) then
!      call clsgks
       stop 'courant'
       endif

       return
     end subroutine diagno_1

     subroutine diagno_2(ux,uz,th,scr1,scr2,rho,nx,nz)
       integer, intent(in) :: nx, nz
       real ux (nx, nz), uz (nx, nz), th (nx, nz), scr1 (nx, nz),   &
            scr2 (nx, nz)                                                     
       real rho (nx, nz) 
       integer :: i, k, nitav, nxz
       real :: amn, amx
       nxz=nx*nz
      call minmax(ux,nxz,amn,amx)
      print 201,amn,amx
 201  format(1x,' --> min, max qv: ',2e12.4)

      call minmax(uz,nxz,amn,amx)
      print 202,amn,amx
 202  format(1x,' --> min, max qc: ',2e12.4)

      call minmax(th,nxz,amn,amx)
      print 203,amn,amx
 203  format(1x,' --> min, max qr: ',2e12.4)

       return
     end subroutine diagno_2

     subroutine minmax(a,n,an,ax)    
       integer, intent(in) :: n
       real, intent(in):: a (n)
       real, intent(out) :: an, ax
       integer :: i

       an= 1.e15
       ax=-1.e15
       do i=1,n
       an=min(a(i),an)
       ax=max(a(i),ax)
       enddo
       return
     end subroutine minmax

   end module diag_mod
