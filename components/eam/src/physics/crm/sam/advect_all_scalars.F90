module advect_all_scalars_mod
  use advect_scalar_mod
  implicit none

contains

  subroutine advect_all_scalars(ncrms)

    use vars
    use microphysics
    use sgs
    use crmtracers
    use params, only: dotracers
    use scalar_momentum_mod
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    integer k,icrm, i, j, kk
    real(crm_rknd), allocatable :: esmt_offset(:)    ! offset for advecting scalar momentum tracers
    real(crm_rknd), allocatable :: dummy(:,:)

    allocate( esmt_offset(ncrms) )
    allocate( dummy(ncrms,nz) )
    call prefetch( esmt_offset )
    call prefetch( dummy )

    !      advection of scalars :
    call advect_scalar(ncrms,t,dummy,dummy)

    !    Advection of microphysics prognostics:
    do k = 1,nmicro_fields
      if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
      .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
      .or. doprecip.and.flag_precip(k).eq.1 ) then
        call advect_scalar(ncrms,micro_field(1,dimx1_s,dimy1_s,1,k),mkadv(1,1,k),mkwle(1,1,k))
      endif
    end do

    !    Advection of sgs prognostics:
    if(dosgs.and.advect_sgs) then
      do k = 1,nsgs_fields
        call advect_scalar(ncrms,sgs_field(1,dimx1_s,dimy1_s,1,k),dummy,dummy)
      end do
    end if

    !   Precipitation fallout:
    if(doprecip) then
      !do icrm = 1 , ncrms
      !  total_water_prec(icrm) = total_water_prec(icrm) + total_water(ncrms,icrm)
      !enddo
      call micro_precip_fall(ncrms)
      !do icrm = 1 , ncrms
      !  total_water_prec(icrm) = total_water_prec(icrm) - total_water(ncrms,icrm)
      !enddo
    end if

    ! advection of tracers:
    !There aren't any of these. We need to delete crmtracers.F90 too at some point
    !if(dotracers) then
    !  do k = 1,ntracers
    !    call advect_scalar(ncrms,icrm,tracer(:,:,:,k,icrm),tradv(:,k,icrm),trwle(:,k,icrm))
    !  end do
    !end if

#if defined(MMF_ESMT)
    ! the esmt_offset simply ensures that the scalar momentum
    ! tracers are positive definite during the advection calculation
    do icrm = 1 , ncrms
      esmt_offset(icrm) = abs( minval( (/ minval(u_esmt(icrm,:,:,:)), minval(v_esmt(icrm,:,:,:)) /) ) ) + 50.
      u_esmt(icrm,:,:,:) = u_esmt(icrm,:,:,:) + esmt_offset(icrm)
      v_esmt(icrm,:,:,:) = v_esmt(icrm,:,:,:) + esmt_offset(icrm)
    enddo
    ! advection of scalar momentum tracers
    call advect_scalar(ncrms,u_esmt,dummy,dummy)
    call advect_scalar(ncrms,v_esmt,dummy,dummy)
    do icrm = 1 , ncrms
      u_esmt(icrm,:,:,:) = u_esmt(icrm,:,:,:) - esmt_offset(icrm)
      v_esmt(icrm,:,:,:) = v_esmt(icrm,:,:,:) - esmt_offset(icrm)
    enddo
#endif

    deallocate( esmt_offset )
    deallocate( dummy )

  end subroutine advect_all_scalars

end module advect_all_scalars_mod
