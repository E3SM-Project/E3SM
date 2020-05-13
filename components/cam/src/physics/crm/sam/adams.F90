module adams_mod
  use params, only: asyncid
  implicit none

contains

  subroutine adams(ncrms)
    !       Adams-Bashforth scheme
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: dtdx, dtdy, dtdz, rhox, rhoy, rhoz
    real(crm_rknd) :: utemp, vtemp, wtemp
    real(crm_rknd) :: dudt_temp,dvdt_temp,dwdt_temp
    integer i,j,k,icrm

    dtdx = dtn/dx
    dtdy = dtn/dy
#if defined(_OPENACC)
    !$acc parallel loop collapse(4) async(asyncid)
#elif defined(_OPENMP)
    !$omp target teams distribute parallel do collapse(4) &
    !$omp private(dudt_temp, dvdt_temp, dwdt_temp, utemp, vtemp, wtemp) 
#endif
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            dtdz = dtn/dz(icrm)
            rhox = rho(icrm,k)*dtdx
            rhoy = rho(icrm,k)*dtdy
            rhoz = rhow(icrm,k)*dtdz

            dudt_temp = dt3(na)*(at*dudt(icrm,i,j,k,na)+bt*dudt(icrm,i,j,k,nb)+ct*dudt(icrm,i,j,k,nc))
            dvdt_temp = dt3(na)*(at*dvdt(icrm,i,j,k,na)+bt*dvdt(icrm,i,j,k,nb)+ct*dvdt(icrm,i,j,k,nc))
            dwdt_temp = dt3(na)*(at*dwdt(icrm,i,j,k,na)+bt*dwdt(icrm,i,j,k,nb)+ct*dwdt(icrm,i,j,k,nc))
            dudt(icrm,i,j,k,nc) = u(icrm,i,j,k) + dudt_temp
            dvdt(icrm,i,j,k,nc) = v(icrm,i,j,k) + dvdt_temp
            dwdt(icrm,i,j,k,nc) = w(icrm,i,j,k) + dwdt_temp

            utemp = 0.5*u(icrm,i,j,k)*rhox
            vtemp = 0.5*v(icrm,i,j,k)*rhoy
            wtemp = 0.5*w(icrm,i,j,k)*rhoz
            u(icrm,i,j,k) = utemp + 0.5*dudt(icrm,i,j,k,nc)*rhox
            v(icrm,i,j,k) = vtemp + 0.5*dvdt(icrm,i,j,k,nc)*rhoy
            w(icrm,i,j,k) = wtemp + 0.5*dwdt(icrm,i,j,k,nc)*rhoz
          
            misc(icrm,i,j,k) = 0.5*(w(icrm,i,j,k)+dwdt(icrm,i,j,k,nc))
          end do
        end do
      end do
    end do

  end subroutine adams

end module adams_mod
