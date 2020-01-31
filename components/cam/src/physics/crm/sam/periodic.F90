module periodic_mod
  use bound_exchange_mod
  use task_util_mod
  implicit none

contains

  subroutine periodic(ncrms,flag)

    use vars
    use microphysics
    use sgs
    use params, only: dotracers, dosgs
    use crmtracers
    use scalar_momentum_mod
#ifdef CLUBB_CRM
    use params, only: doclubb, doclubbnoninter
#endif
    implicit none
    integer, intent(in) :: ncrms,flag
    integer :: i,icrm, j, ii, k

    !-------------------------------------------------
    ! Update velocity fields
    !-------------------------------------------------
    if(flag.eq.0) then
      call bound_exchange(ncrms,u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1,1,1,1,1)
      call bound_exchange(ncrms,v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1,1,1,1,2)
      ! use w at the top level  - 0s anyway - to exchange the sst boundaries (for
      ! surface fluxes call
      !$acc parallel loop collapse(3) async(asyncid)
      do j = 1 , ny
        do i = 1 , nx
          do icrm = 1 , ncrms
            w(icrm,i,j,nz) = sstxy(icrm,i,j)
          enddo
        enddo
      enddo
      call bound_exchange(ncrms,w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,1,1,1,1,3)
      !$acc parallel loop collapse(3) async(asyncid)
      do j = 1-YES3D , ny+YES3D
        do i = 0 , nx+1
          do icrm = 1 , ncrms
            if ( i >= 0 .and. i <= nx   .and. j >= 1-YES3D .and. j <= ny       ) sstxy(icrm,i,j) = w(icrm,i,j,nz)
          enddo
        enddo
      enddo
      !$acc parallel loop collapse(3) async(asyncid)
      do j = 1-YES3D , ny+YES3D
        do i = 0 , nx+1
          do icrm = 1 , ncrms
            if ( i >= 0 .and. i <= nx+1 .and. j >= 1-YES3D .and. j <= ny+YES3D ) w(icrm,i,j,nz) = 0.
          enddo
        enddo
      enddo
    endif

    !-------------------------------------------------
    ! update prognostic scalar fields for advection
    !-------------------------------------------------
    if(flag.eq.2) then
      call bound_exchange(ncrms,u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2+NADV,2+NADV,1)
      call bound_exchange(ncrms,v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2+NADV,2+NADV,2,3,2)
      call bound_exchange(ncrms,w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2+NADV,2+NADV,2+NADV,2+NADV,3)
      call bound_exchange(ncrms,t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4)
      do i = 1,nsgs_fields
        if(dosgs.and.advect_sgs) call bound_exchange(ncrms,sgs_field(:,:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+i)
      end do
      do i = 1,nmicro_fields
        if(   i.eq.index_water_vapor             &
#ifdef CLUBB_CRM
        ! Vince Larson (UWM) changed so that bound_exchange is called even if
        !     docloud = .false. and doclubb = .true.    11 Nov 2007
        .or. (docloud.or.doclubb.or.doclubbnoninter) .and.flag_precip(i).ne.1    &
#else
        .or. docloud.and.flag_precip(i).ne.1    &
#endif
        .or. doprecip.and.flag_precip(i).eq.1 ) then
          call bound_exchange(ncrms,micro_field(:,:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+nsgs_fields+nsgs_fields_diag+i)
        endif
      end do
      !if(dotracers) then
      !  do i=1,ntracers
      !    call bound_exchange(tracer(:,:,:,i,icrm),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+nsgs_fields+nsgs_fields_diag+nmicro_fields+i)
      !  end do
      !end if
#if defined(SP_ESMT)
      call bound_exchange(ncrms,u_esmt,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+nsgs_fields+nsgs_fields_diag+nmicro_fields+ntracers+1)
      call bound_exchange(ncrms,v_esmt,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4+nsgs_fields+nsgs_fields_diag+nmicro_fields+ntracers+2)
#endif
    endif

    !-------------------------------------------------
    ! Update all scalars before SGS
    !-------------------------------------------------
    if(flag.eq.3) then
      call bound_exchange(ncrms,t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4)
      do i = 1,nsgs_fields
        if(dosgs.and.advect_sgs) &
        call bound_exchange(ncrms,sgs_field(:,:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4+i)
      end do
      do i = 1,nmicro_fields
        if(   i.eq.index_water_vapor             &
#ifdef CLUBB_CRM
        ! Vince Larson (UWM) changed so that bound_exchange is called even if
        !     docloud = .false. and doclubb = .true.    11 Nov 2007
        .or. (docloud.or.doclubb.or.doclubbnoninter) .and.flag_precip(i).ne.1    &
#else
        .or. docloud.and.flag_precip(i).ne.1    &
#endif
        .or. doprecip.and.flag_precip(i).eq.1 ) then
          call bound_exchange(ncrms,micro_field(:,:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4+nsgs_fields+nsgs_fields_diag+i)
        endif
      end do
      !if(dotracers) then
      !  do i=1,ntracers
      !    call bound_exchange(tracer(:,:,:,i,icrm),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4+nsgs_fields+nsgs_fields_diag+nmicro_fields+i)
      !  end do
      !end if
#if defined(SP_ESMT)
      call bound_exchange(ncrms,u_esmt,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4+nsgs_fields+nsgs_fields_diag+nmicro_fields+ntracers+1)
      call bound_exchange(ncrms,v_esmt,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4+nsgs_fields+nsgs_fields_diag+nmicro_fields+ntracers+2)
#endif
    endif

    !-------------------------------------------------
    ! SGS diagnostic fields
    !-------------------------------------------------
    if(flag.eq.4) then
      do i = 1,nsgs_fields_diag
        if(dosgs.and.do_sgsdiag_bound) &
        call bound_exchange(ncrms,sgs_field_diag(:,:,:,:,i),dimx1_d,dimx2_d,dimy1_d,dimy2_d,nzm,1+dimx1_d,dimx2_d-nx,YES3D+dimy1_d,1-YES3D+dimy2_d-ny,4+nsgs_fields+i)
      end do
    end if

  end subroutine periodic

end module periodic_mod
