module ice_fall_mod
  implicit none

contains

  subroutine ice_fall(ncrms)
    ! Sedimentation of ice:
    use vars
    use microphysics, only: micro_field, index_cloud_ice
    !use micro_params
    use params
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    integer, allocatable :: kmax(:)
    integer, allocatable :: kmin(:)
    real(crm_rknd), allocatable :: fz(:,:,:,:)
    integer :: i,j,k, kb, kc, ici,icrm
    real(crm_rknd) coef,dqi,lat_heat,vt_ice
    real(crm_rknd) omnu, omnc, omnd, qiu, qic, qid, tmp_theta, tmp_phi

    allocate( kmax(ncrms) )
    allocate( kmin(ncrms) )
    allocate( fz(ncrms,nx,ny,nz) )
    call prefetch( kmax )
    call prefetch( kmin )
    call prefetch( fz )

    !$acc parallel loop async(asyncid)
    do icrm = 1 , ncrms
      kmax(icrm)=0
      kmin(icrm)=nzm+1
    enddo
    !$acc parallel loop collapse(3) async(asyncid)
    do j = 1, ny
      do i = 1, nx
        do icrm = 1 , ncrms
          do k = 1,nzm
            if(qcl(icrm,i,j,k)+qci(icrm,i,j,k).gt.0.D0.and. tabs(icrm,i,j,k).lt.273.15D0) then
              !$acc atomic update
              kmin(icrm) = min(kmin(icrm),k)
              !$acc atomic update
              kmax(icrm) = max(kmax(icrm),k)
            end if
          end do
        end do
      end do
    end do
    !$acc parallel loop collapse(2) async(asyncid)
    do k = 1,nzm
      do icrm = 1 , ncrms
        qifall(icrm,k) = 0.D0
        tlatqi(icrm,k) = 0.D0
      end do
    end do

    if(index_cloud_ice.eq.-1) return

    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1,nz
      do j = 1, ny
        do i = 1, nx
          do icrm = 1 , ncrms
            fz(icrm,i,j,k) = 0.D0
          end do
        end do
      end do
    end do

    ! Compute cloud ice flux (using flux limited advection scheme, as in
    ! chapter 6 of Finite Volume Methods for Hyperbolic Problems by R.J.
    ! LeVeque, Cambridge University Press, 2002).
    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1 , nz
      do j = 1,ny
        do i = 1,nx
          do icrm = 1 , ncrms
            if (k >= max(1,kmin(icrm)-1) .and. k <= kmax(icrm) ) then
              ! Set up indices for x-y planes above and below current plane.
              kc = min(nzm,k+1)
              kb = max(1,k-1)
              ! CFL number based on grid spacing interpolated to interface i,j,k-1/2
              coef = dtn/(0.5D0*(adz(icrm,kb)+adz(icrm,k))*dz(icrm))

              ! Compute cloud ice density in this cell and the ones above/below.
              ! Since cloud ice is falling, the above cell is u(icrm,upwind),
              ! this cell is c (center) and the one below is d (downwind).
              qiu = rho(icrm,kc)*qci(icrm,i,j,kc)
              qic = rho(icrm,k) *qci(icrm,i,j,k)
              qid = rho(icrm,kb)*qci(icrm,i,j,kb)

              ! Ice sedimentation velocity depends on ice content. The fiting is
              ! based on the data by Heymsfield (JAS,2003). -Marat
              vt_ice = min(real(0.4D0,crm_rknd),8.66D0*(max(real(0.,crm_rknd),qic)+1.D-10)**0.24D0)   ! Heymsfield (JAS, 2003, p.2607)

              ! Use MC flux limiter in computation of flux correction.
              ! (MC = monotonized centered difference).
              !         if (qic.eq.qid) then
              if (abs(qic-qid).lt.1.0D-25) then  ! when qic, and qid is very small, qic_qid can still be zero
                ! even if qic is not equal to qid. so add a fix here +++mhwang
                tmp_phi = 0.
              else
                tmp_theta = (qiu-qic)/(qic-qid)
                tmp_phi = max(real(0.,crm_rknd),min(0.5D0*(1.+tmp_theta),real(2.D0,crm_rknd),2.D0*tmp_theta))
              end if

              ! Compute limited flux.
              ! Since falling cloud ice is a 1D advection problem, this
              ! flux-limited advection scheme is monotonic.
              fz(icrm,i,j,k) = -vt_ice*(qic - 0.5D0*(1.D0-coef*vt_ice)*tmp_phi*(qic-qid))
            endif
          end do
        end do
      end do
    enddo
    !$acc parallel loop collapse(3) async(asyncid)
    do j = 1, ny
      do i = 1, nx
        do icrm = 1 , ncrms
          fz(icrm,i,j,nz) = 0.
        end do
      end do
    end do

    ici = index_cloud_ice

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1, nz
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            if ( k >= max(1,kmin(icrm)-2) .and. k <= kmax(icrm) ) then
            coef=dtn/(dz(icrm)*adz(icrm,k)*rho(icrm,k))
            ! The cloud ice increment is the difference of the fluxes.
            dqi=coef*(fz(icrm,i,j,k)-fz(icrm,i,j,k+1))
            ! Add this increment to both non-precipitating and total water.
            micro_field(icrm,i,j,k,ici)  = micro_field(icrm,i,j,k,ici)  + dqi
            ! Include this effect in the total moisture budget.
            !$acc atomic update
            qifall(icrm,k) = qifall(icrm,k) + dqi

            ! The latent heat flux induced by the falling cloud ice enters
            ! the liquid-ice static energy budget in the same way as the
            ! precipitation.  Note: use latent heat of sublimation.
            lat_heat  = (fac_cond+fac_fus)*dqi
            ! Add divergence of latent heat flux to liquid-ice static energy.
            t(icrm,i,j,k)  = t(icrm,i,j,k)  - lat_heat
            ! Add divergence to liquid-ice static energy budget.
            !$acc atomic update
            tlatqi(icrm,k) = tlatqi(icrm,k) - lat_heat
            endif
          end do
        end do
      end do
    end do

    !$acc parallel loop collapse(3) async(asyncid)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          coef=dtn/dz(icrm)
          dqi=-coef*fz(icrm,i,j,1)
          precsfc(icrm,i,j) = precsfc(icrm,i,j)+dqi
          precssfc(icrm,i,j) = precssfc(icrm,i,j)+dqi
        end do
      end do
    end do

    deallocate( kmax )
    deallocate( kmin )
    deallocate( fz )

  end subroutine ice_fall

end module ice_fall_mod
