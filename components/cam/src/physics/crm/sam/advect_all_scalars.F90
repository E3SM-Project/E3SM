module advect_all_scalars_mod
  use advect_scalar_mod
#if defined(_OPENACC)
    use openacc_utils
#endif

  implicit none

contains

  subroutine advect_all_scalars(ncrms)
    use vars
    use microphysics
    use sgs
    use crmtracers
#ifdef CLUBB_CRM
    use params, only: dotracers, doclubb, doclubbnoninter
#else
    use params, only: dotracers
#endif
    use scalar_momentum_mod
    implicit none
    integer, intent(in) :: ncrms
    integer k,icrm, i, j, kk
    real(crm_rknd), allocatable :: esmt_offset(:)    ! whannah - offset for advecting scalar momentum tracers
    real(crm_rknd), allocatable :: fadv_dummy(:,:)
    real(crm_rknd), allocatable :: flux_dummy(:,:)
    integer, save :: icall = 0
    icall = icall + 1

    allocate( esmt_offset(ncrms) )
    allocate( fadv_dummy(ncrms,nz) )
    allocate( flux_dummy(ncrms,nz) )
#if defined(_OPENACC)
    call prefetch( esmt_offset )
    call prefetch( fadv_dummy )
    call prefetch( flux_dummy )
#elif defined(_OPENMP)
    !$omp target enter data map(alloc: esmt_offset)
    !$omp target enter data map(alloc: fadv_dummy)
    !$omp target enter data map(alloc: flux_dummy)
#endif
#if defined(_OPENMP)
    !$omp target update from (flag_precip)
#endif
!$omp target update from(micro_field)
!$omp target update from(qn)
!$omp target update from(t)

    !      advection of scalars :
    call advect_scalar(ncrms,t,fadv_dummy,flux_dummy)

!$omp target update from(micro_field)
!$omp target update from(qn)
!$omp target update from(t)

    !    Advection of microphysics prognostics:
    do k = 1,nmicro_fields
      if(k.eq.index_water_vapor             &! transport water-vapor variable no metter what
#ifdef CLUBB_CRM
      !Added preprocessor directives. - nielsenb UWM 30 July 2008
      .or. (docloud .or. doclubb .or. doclubbnoninter) .and.flag_precip(k).ne.1    & ! transport non-precipitation vars
#else
      .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
#endif
      .or. doprecip.and.flag_precip(k).eq.1 ) then
        call advect_scalar(ncrms,micro_field(1,dimx1_s,dimy1_s,1,k),mkadv(1,1,k),mkwle(1,1,k))
      endif
    end do

    !    Advection of sgs prognostics:
    if(dosgs.and.advect_sgs) then
      do k = 1,nsgs_fields
        call advect_scalar(ncrms,sgs_field(1,dimx1_s,dimy1_s,1,k),fadv_dummy,flux_dummy)
      end do
    end if

    !   Precipitation fallout:
    if(doprecip) then
      !do icrm = 1 , ncrms
      !  total_water_prec(icrm) = total_water_prec(icrm) + total_water(ncrms,icrm)
      !enddo
!$omp target update from(micro_field)
!$omp target update from(qn)
!$omp target update from(t)

      call micro_precip_fall(ncrms)

!$omp target update from(micro_field)
!$omp target update from(qn)
!$omp target update from(t)

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

#if defined(SP_ESMT)
    ! whannah - the esmt_offset simply ensures that the scalar momentum
    ! tracers are positive definite during the advection calculation
    do icrm = 1 , ncrms
      esmt_offset(icrm) = abs( minval( (/ minval(u_esmt(icrm,:,:,:)), minval(v_esmt(icrm,:,:,:)) /) ) ) + 50.
      u_esmt(icrm,:,:,:) = u_esmt(icrm,:,:,:) + esmt_offset(icrm)
      v_esmt(icrm,:,:,:) = v_esmt(icrm,:,:,:) + esmt_offset(icrm)
    enddo
    ! advection of scalar momentum tracers
    call advect_scalar(ncrms,u_esmt,fadv_dummy,flux_dummy)
    call advect_scalar(ncrms,v_esmt,fadv_dummy,flux_dummy)
    do icrm = 1 , ncrms
      u_esmt(icrm,:,:,:) = u_esmt(icrm,:,:,:) - esmt_offset(icrm)
      v_esmt(icrm,:,:,:) = v_esmt(icrm,:,:,:) - esmt_offset(icrm)
    enddo
#endif
#if defined(_OPENMP)
    !$omp target exit data map(delete: esmt_offset)
    !$omp target exit data map(delete: fadv_dummy)
    !$omp target exit data map(delete: flux_dummy)
#endif
    deallocate( esmt_offset )
    deallocate( fadv_dummy )
    deallocate( flux_dummy )

  end subroutine advect_all_scalars

end module advect_all_scalars_mod
