module absor_mod
  use grid_init_mod, only :  dt
  use prof_init_mod, only : ux_e, th_e, qv_e
  real, allocatable :: tau(:)
contains                                                                        
  SUBROUTINE absor (ux, uz, th, qv, qc, qr, fth, fqv, fqc, fqr, tau, nx, nz) 
    implicit none
    integer, intent(in) :: nx, nz
    real, intent(in) :: tau(nz)
    real, intent(inout) :: ux(nx,nz), uz(nx,nz),  th(nx,nz), qv(nx,nz), &
         qc(nx,nz), qr(nx,nz), fqv(nx,nz), fqc(nx,nz), fqr(nx,nz) , fth(nx,nz)

    integer :: i, k
    real ::  coe1, coe2

    DO k = 1, nz 
       coe1 = .5 * dt * tau (k) 
       coe2 = 1. + coe1 
       DO i = 1, nx 
          ux (i, k) = (ux (i, k) + coe1 * ux_e (k) ) / coe2 
          uz (i, k) = (uz (i, k) ) / coe2 
          th (i, k) = (th (i, k) + coe1 * th_e (k) ) / coe2 
          qv (i, k) = (qv (i, k) + coe1 * qv_e (k) ) / coe2 
          qc (i, k) = (qc (i, k) ) / coe2 
          qr (i, k) = (qr (i, k) ) / coe2 

          fth (i, k) = fth (i, k) - tau (k) * (th (i, k) - th_e (k) ) 
          fqv (i, k) = fqv (i, k) - tau (k) * (qv (i, k) - qv_e (k) ) 
          fqc (i, k) = fqc (i, k) - tau (k) * (qc (i, k) - 0.) 
          fqr (i, k) = fqr (i, k) - tau (k) * (qr (i, k) - 0.) 
       enddo
    enddo

    RETURN 
  END SUBROUTINE absor
end module absor_mod
