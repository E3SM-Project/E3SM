module conditional_diag_main

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  use conditional_diag, only: cnd_diag_info, cnd_diag_info_t

  implicit none

  private

  public cnd_diag_checkpoint

  ! types of comparison used in defining sampling condition

  integer, parameter :: GE  =  2
  integer, parameter :: GT  =  1
  integer, parameter :: EQ  =  0
  integer, parameter :: LT  = -1
  integer, parameter :: LE  = -2

  ! flags for grid cell/column masking

  real(r8),parameter :: ON  = 1._r8
  real(r8),parameter :: OFF = 0._r8

  ! start time step for increment calculation

  integer,parameter :: NS0INC = 3 ! start time step for increment calculation
  integer,parameter :: NS0SMP = 4 ! start time step for conditional sampling
contains

!======================================================
subroutine cnd_diag_checkpoint( diag, this_chkpt, state, pbuf, cam_in, cam_out )

  use perf_mod,            only: t_startf, t_stopf
  use time_manager,        only: get_nstep
  use ppgrid,              only: pcols, pver
  use cam_history_support, only: max_fieldname_len
  use cam_history,         only: outfld

  use physics_types,    only: physics_state
  use camsrfexch,       only: cam_in_t, cam_out_t
  use physics_buffer,   only: physics_buffer_desc

  use conditional_diag,    only: cnd_diag_t, FILLVALUE, NODP
  use conditional_diag_output_utils, only: get_metric_and_flag_names_for_output, &
                                           get_fld_name_for_output

  type(cnd_diag_t),    intent(inout), target :: diag
  character(len=*),    intent(in)            :: this_chkpt

  type(physics_state),     intent(in) :: state
  type(physics_buffer_desc),pointer   :: pbuf(:)
  type(cam_in_t),          intent(in) :: cam_in
  type(cam_out_t),optional,intent(in) :: cam_out

  integer :: ncnd, nchkpt, nqoi
  integer :: icnd, ichkpt, ii, iqoi
  integer :: ncol, lchnk
  integer :: nstep
  integer :: x_dp
  character(len=20) :: qoi_name

  real(r8) :: new3d(pcols,pver)

  real(r8),pointer     :: metric(:,:), flag(:,:), inc(:,:), old(:,:)
  real(r8),allocatable :: new(:,:)

  character(len=max_fieldname_len) :: outfldname, flag_name_out, metric_name_out

  !---------------------------------------------------------------------------
  call t_startf('cnd_diag_checkpoint')
  if (cnd_diag_info%ncnd == 0) then ! no conditional diagnostics requested 

     call t_stopf('cnd_diag_checkpoint')
     return

  else
  !---------------------------------------------------------------------------

  ncnd   = cnd_diag_info%ncnd
  nchkpt = cnd_diag_info%nchkpt
  nqoi   = cnd_diag_info%nqoi

  lchnk  = state%lchnk
  ncol   = state%ncol

  nstep  = get_nstep()  ! current time step. Later in this routine, we
                        ! - save QoI values starting from the first step.
                        ! - do increment calculation when nstep >= NS0INC.
                        ! - evaluate sampling metrics starting from the first step.
                        ! - do conditional sampling and outfld calls only when 
                        !   nstep >= NS0SMP, because the sampling time window might 
                        !   involve some checkpints from the previous nstep,
                        !   and valid increments are available only from nstep >= NS0INC.

  !=======================================
  ! Obtain QoI values and/or increments
  !=======================================
  ! First determine if this checkpoint is active for QoI monitoring

  ichkpt = 0  ! 0 = checkpoint inactive; this is the default

  do ii = 1,nchkpt
     if ( trim(cnd_diag_info%qoi_chkpt(ii)) == trim(this_chkpt) ) then 
        ichkpt = ii
        exit
     end if
  end do

  if (ichkpt>0) then 
  !---------------------------------------------------------------------------
  ! This checkpoint is active for QoI monitoring. Obtain the QoI values 
  ! and/or their increments if needed, and save to the "diag" data structure. 
  ! Note that here we only obtain and save the QoIs and/or increments.
  ! Conditional sampling won't be applied until the "cnd_end_chkpt" checkpoint
  !---------------------------------------------------------------------------
     do iqoi = 1,nqoi

        qoi_name = cnd_diag_info% qoi_name(iqoi)
        x_dp     = cnd_diag_info% x_dp(iqoi,ichkpt)

        ! Allocate memory for tmp array. This has to be done inside the iqoi
        ! loop as different QoIs might have different numbers of vertical levels

        allocate( new( pcols,cnd_diag_info%qoi_nver_save(iqoi) ))

        !------------------------------------------------------------------------
        ! Obtain the most up-to-date values of the QoI or its vertical integral
        !------------------------------------------------------------------------
        if (x_dp==NODP) then ! get the new QoI values

           call get_values( new, trim(qoi_name), state, pbuf, cam_in, cam_out ) !out, in, ... ,in

        else ! get the vertical integral of the new QoI values

           call get_values( new3d, trim(qoi_name), state, pbuf, cam_in, cam_out ) !out, in, ... ,in
           call mass_wtd_vert_intg( new, new3d, state, x_dp ) !out, in, in, in

        end if
        !-----------------------------------------------------------------------------------------
        ! Save the QoI and its increment to corresponding components of the "diag" data structure
        !-----------------------------------------------------------------------------------------
        ! The current implementation is such that the same set of checkpoints and QoIs are 
        ! monitored for all different sampling conditions. Below, the same QoI values 
        ! and/or increments are copied to all conditions. When conditional sampling is applied 
        ! later, the different copies will likely be sampled differently. In the future, if we 
        ! decide to allow for different QoIs and/or checkpoints under different sampling 
        ! conditions, then the loops in this subroutine will need to be reworked.
        !--------------------------------------------------------------------------
        ! Save the QoI 

        if (cnd_diag_info%l_output_state) then
           do icnd = 1,ncnd
              diag%cnd(icnd)%qoi(iqoi)% val(1:ncol,:,ichkpt) = new(1:ncol,:) 
           end do
        end if

        !-------------------------------
        ! Calculate and save inrements

        if (cnd_diag_info%l_output_incrm) then

           icnd = 1
           old => diag%cnd(icnd)%qoi(iqoi)% old

           if (nstep >= NS0INC) then 
              inc => diag%cnd(icnd)%qoi(iqoi)% inc(:,:,ichkpt)
              inc(1:ncol,:) = new(1:ncol,:) - old(1:ncol,:)
           end if

           ! Save the current value of QoI as "old" for next checkpoint

           old(1:ncol,:) = new(1:ncol,:)

           ! Save increments for other sampling conditions; update "old" value

           do icnd = 2,ncnd

              diag%cnd(icnd)%qoi(iqoi)% old(1:ncol,:) = new(1:ncol,:)
              if (nstep >= NS0INC) &
              diag%cnd(icnd)%qoi(iqoi)% inc(1:ncol,:,ichkpt) = inc(1:ncol,:)

           end do

        end if !l_output_incrm
        !----------------------

        ! Calculations done for this QoI. Clean up.
        deallocate(new)

     end do ! iqoi = 1,nqoi
  end if    ! ichkpt > 0

  !=======================================================
  ! Evaluate sampling condition if this is cnd_eval_chkpt 
  !=======================================================
  do icnd = 1,ncnd

     ! Check if sampling condition needs to be evaluated at this checkpoint.
     ! The answer could be .t. for multiple icnd values

     if (trim(cnd_diag_info% cnd_eval_chkpt(icnd)) .eq. trim(this_chkpt)) then 

        !---------------------------------
        ! Get metric values and set flags 
        !---------------------------------
        metric => diag%cnd(icnd)%metric
        flag   => diag%cnd(icnd)%flag

        !----------------------------------------------------------------
        ! The special metric name "ALL" is used to select all grid cells.

        if (trim(cnd_diag_info%metric_name(icnd)).eq.'ALL') then

           metric = 1._r8
           flag   = ON

        !------------------------------------------
        else ! If NOT selecting all grid cells...

           call get_values( metric,                                &! out
                            trim(cnd_diag_info%metric_name(icnd)), &! in
                            state, pbuf, cam_in, cam_out )          ! in

           call get_flags( metric, icnd, ncol, cnd_diag_info, flag )

           ! Apply conditional sampling to metric

           where(flag.eq.OFF)  metric = FILLVALUE

        end if

        !----------------------------------------------------
        ! Send both metric and flag values to history buffer
        !----------------------------------------------------
        call get_metric_and_flag_names_for_output( icnd, cnd_diag_info, metric_name_out, flag_name_out )

        call outfld( trim(metric_name_out), metric, pcols, lchnk )
        call outfld( trim(flag_name_out),     flag, pcols, lchnk )

     end if !right chkpt
  end do    !icnd

  !=======================================================================
  ! Apply conditional sampling, then send QoIs to history buffer.
  ! (Do this only when nstep >= NS0SMP, because the sampling time window might 
  ! involve some checkpints from the previous nstep,
  ! and valid increments are available only from nstep = NS0INC onwards)
  !=======================================================================
  if (nstep >= NS0SMP) then
  do icnd = 1,ncnd

     ! Check if conditional sampling needs to be completed at this checkpoint.
     ! The answer could be .t. for multiple icnd values

     if (trim(cnd_diag_info% cnd_end_chkpt(icnd)).eq.trim(this_chkpt)) then 

        ! Each sampling condition has its own flags that will be applied
        ! to all QoIs and checkpoints

        flag => diag%cnd(icnd)%flag

        !----------------------------------------------------------------
        ! Apply conditional sampling to QoIs and/or their increments,
        ! then do the outfld calls to send the values to history buffer. 
        ! In subroutine apply_masking,
        ! different actions are taken depending on the vertical 
        ! dimension sizes of the metric and the QoIs.
        !----------------------------------------------------------------
        if (cnd_diag_info%l_output_state) then        

           ! Apply conditional sampling to QoI values

           if (trim(cnd_diag_info%metric_name(icnd)).ne.'ALL') then
             do iqoi = 1,nqoi
                call apply_masking( flag, diag%cnd(icnd)%qoi(iqoi)%val, FILLVALUE )
             end do
           end if

           ! Send QoI values to history buffer

           do iqoi = 1,nqoi
             do ichkpt = 1,nchkpt
                call get_fld_name_for_output( '', cnd_diag_info, icnd, iqoi, ichkpt, outfldname)
                call outfld( trim(outfldname), diag%cnd(icnd)%qoi(iqoi)%val(:,:,ichkpt), pcols, lchnk )
             end do
           end do

        end if

        !------------
        if (cnd_diag_info%l_output_incrm) then

           ! Apply conditional sampling to QoI increments

           if (trim(cnd_diag_info%metric_name(icnd)).ne.'ALL') then
             do iqoi = 1,nqoi
                call apply_masking( flag, diag%cnd(icnd)%qoi(iqoi)%inc, FILLVALUE )
             end do
           end if

           ! Send QoI increments to history buffer

           do iqoi = 1,nqoi
             do ichkpt = 1,nchkpt
                call get_fld_name_for_output( '_inc', cnd_diag_info, icnd, iqoi, ichkpt, outfldname)
                call outfld( trim(outfldname), diag%cnd(icnd)%qoi(iqoi)%inc(:,:,ichkpt), pcols, lchnk )
             end do
           end do

        end if

     end if  !trim(this_chkpt).eq.trim(cnd_diag_info% cnd_end_chkpt(icnd))
  end do ! icnd = 1,ncnd
  end if ! nstep >= NS0SMP

  call t_stopf('cnd_diag_checkpoint')

  end if ! cnd_diag_info%ncnd == 0

end subroutine cnd_diag_checkpoint

!==================================================================
subroutine apply_masking( flag, array, fillvalue )

    real(r8),intent(in)    ::  flag(:,:)
    real(r8),intent(inout) :: array(:,:,:)
    real(r8),intent(in)    :: fillvalue

    integer :: kk             ! vertical level index for a loop
    integer :: flag_nver      ! # of vertical levels the flag array has
    integer :: array_nver     ! # of vertical levels the output array has
    integer :: nchkpt, ichkpt ! # of checkpoints to process, and the loop index
    integer :: pcols, icol    ! # of columns in grid chunk, and the loop index

         pcols = size(flag, 1)
     flag_nver = size(flag, 2) 
    array_nver = size(array,2) 
        nchkpt = size(array,3)

    if (flag_nver == array_nver) then 
    ! same vertical dimension size; simply apply masking - and do this 
    ! for all checkpoints

       do ichkpt = 1,nchkpt
          where(flag(:,:).eq.OFF) array(:,:,ichkpt) = fillvalue 
       end do

    elseif (flag_nver == 1 .and. array_nver > 1) then 
    ! apply the same masking to all vertical levels and checkpoints

       do icol = 1,pcols
          if (flag(icol,1).eq.OFF) array(icol,:,:) = fillvalue
       end do

    elseif (flag_nver > 1 .and. array_nver == 1) then
    ! if any level in a grid column is selected, select that column;
    ! In other words, mask out a column in the output array 
    ! only if cells on all levels in that column are masked out.

       do icol = 1,pcols
          if (all( flag(icol,:).eq.OFF )) then
             array(icol,1,:) = fillvalue
          end if
       end do

    else
    ! QoI and flag both have multiple levels but the number of 
    ! vertical levels are different. Do not apply masking.

       continue

    end if

end subroutine apply_masking 


!========================================================
subroutine get_values( arrayout, varname, state, pbuf, cam_in, cam_out )

  use ppgrid,         only: pcols,pver
  use physics_types,  only: physics_state
  use camsrfexch,     only: cam_in_t, cam_out_t
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use time_manager,   only: get_nstep
  use constituents,   only: cnst_get_ind, pcnst, sflxnam
  use misc_diagnostics
  use wv_saturation, only: qsat_water

  real(r8),           intent(out)   :: arrayout(:,:)
  character(len=*),   intent(in)    :: varname
  type(physics_state),intent(in)    :: state
  type(physics_buffer_desc), pointer:: pbuf(:)
  type(cam_in_t),     intent(in)    :: cam_in
  type(cam_out_t),    intent(in)    :: cam_out

  real(r8),pointer :: ptr2d(:,:)
  real(r8),pointer :: ptr1d(:)

  real(r8) :: tmp(pcols,pver)

  character(len=*),parameter :: subname = 'conditional_diag_main:get_values'

  integer :: ncol, idx, m

  ncol = state%ncol

  !--------------------------------------------------------------------------------
  ! If the requested variable is one of the advected tracers, get it from state%q
  !--------------------------------------------------------------------------------
  ! cnst_get_ind returns the index of a tracer in the host  model's advected 
  ! tracer array. If variable is not found on the tracer list, and index value 
  ! of -1 will be returned

  call cnst_get_ind(trim(adjustl(varname)),idx, abrtf=.false.)  !in, out, in

  if (idx /= -1) then ! This variable is a tracer field
     arrayout(1:ncol,:) = state%q(1:ncol,:,idx)

  else
  !--------------------------------------------------------------------------------
  ! If the requested variable is the surface flux of an advected tracer, get it
  ! from cam_in%cflx
  !--------------------------------------------------------------------------------
     do m = 1,pcnst
        if (trim(adjustl(varname)).eq.trim(sflxnam(m))) then
           idx = m ; exit
        end if
     end do

     if (idx /= -1) then ! This variable is the surface flux of an advected tracer
        arrayout(1:ncol,1) = cam_in%cflx(1:ncol,idx)

     else
     !============================================================================
     ! Non-tracer and non-sfc-flux variables
     !============================================================================

        select case (trim(adjustl(varname)))
        case('T')
           arrayout(1:ncol,:) = state%t(1:ncol,:)

        case('U')
           arrayout(1:ncol,:) = state%u(1:ncol,:)

        case('V')
           arrayout(1:ncol,:) = state%v(1:ncol,:)

        case('OMEGA')
           arrayout(1:ncol,:) = state%omega(1:ncol,:)

        case('PMID')
           arrayout(1:ncol,:) = state%pmid(1:ncol,:)

        case('PINT')
           arrayout(1:ncol,:) = state%pint(1:ncol,:)

        case('ZM')
           arrayout(1:ncol,:) = state%zm(1:ncol,:)

        case('ZI')
           arrayout(1:ncol,:) = state%zi(1:ncol,:)

        case('PS')
           arrayout(1:ncol,1) = state%ps(1:ncol)

        !---------
        ! cam_in 
        !---------

        case('LWUP')
           arrayout(1:ncol,1) = cam_in%lwup(1:ncol)

        case('LHF')
           arrayout(1:ncol,1) = cam_in%lhf(1:ncol)

        case('SHF')
           arrayout(1:ncol,1) = cam_in%shf(1:ncol)

        case('WSX')
           arrayout(1:ncol,1) = cam_in%wsx(1:ncol)

        case('WSY')
           arrayout(1:ncol,1) = cam_in%wsy(1:ncol)

        case('TREF')
           arrayout(1:ncol,1) = cam_in%tref(1:ncol)

        case('QREF')
           arrayout(1:ncol,1) = cam_in%qref(1:ncol)

        case('U10')
           arrayout(1:ncol,1) = cam_in%u10(1:ncol)

        case('TS')
           arrayout(1:ncol,1) = cam_in%ts(1:ncol)

        case('SST')
           arrayout(1:ncol,1) = cam_in%sst(1:ncol)

        !----------
        ! cam_out
        !----------

        case('FLWDS')
           arrayout(1:ncol,1) = cam_out%flwds(1:ncol)

        case('NETSW')
           arrayout(1:ncol,1) = cam_out%netsw(1:ncol)

        !----------
        ! pbuf
        !----------

        case('PBLH')
            idx = pbuf_get_index('pblh') ; call pbuf_get_field( pbuf, idx, ptr1d )
            arrayout(:,1) = ptr1d

        case('CLD')
            idx = pbuf_get_index('CLD')  ; call pbuf_get_field( pbuf, idx, ptr2d )
            arrayout(:,:) = ptr2d

        !-----------------------------------------------------------
        ! physical quantities that need to be calculated on the fly 
        !-----------------------------------------------------------

        case ('QSATW')
          call qsat_water( state%t(:ncol,:), state%pmid(:ncol,:), &! in
                               tmp(:ncol,:),   arrayout(:ncol,:)  )! out

        case ('QSATI')
          call qsat_ice( ncol, pver, state%t(:ncol,:), state%pmid(:ncol,:), &! in
                         arrayout(:ncol,:)                                  )! out

        case ('QSSATW')  ! supersaturation w.r.t. water in terms of mixing ratio

          call supersat_q_water( ncol, pver, state%t(:ncol,:),            &! in
                                 state%pmid(:ncol,:), state%q(:ncol,:,1), &! in
                                 arrayout(:ncol,:)  )                      ! out

        case ('QSSATI')  ! supersaturation w.r.t. ice in terms of mixing ratio

          call supersat_q_ice(   ncol, pver, state%t(:ncol,:),            &! in
                                 state%pmid(:ncol,:), state%q(:ncol,:,1), &! in
                                 arrayout(:ncol,:)  )                      ! out

        case ('RHW')
          call relhum_water_percent( ncol, pver, state%t(:ncol,:),            &! in
                                     state%pmid(:ncol,:), state%q(:ncol,:,1), &! in
                                     arrayout(:ncol,:)  )                      ! out
        case ('RHI')
          call relhum_ice_percent( ncol, pver, state%t(:ncol,:),            &! in
                                   state%pmid(:ncol,:), state%q(:ncol,:,1), &! in
                                   arrayout(:ncol,:)  )                      ! out

        case ('CAPE')
          call compute_cape( state, pbuf, pcols, pver, arrayout(:,1) ) ! in, in, in, out

        !-----------------------------------------------------------------------------------
        ! The following were added mostly for testing of the conditional diag functionality
        !-----------------------------------------------------------------------------------
        case('NSTEP')
           arrayout(1:ncol,:) = get_nstep()

        !case('LAT')
        !   arrayout(1:ncol,1) = state%lat(1:ncol)

        !case('LON')
        !   arrayout(1:ncol,1) = state%lon(1:ncol)

        case default
           call endrun(subname//': unknow varname - '//trim(varname))
        end select

     end if !whether the requested variable is the surface flux of a tracer
  end if    !whether the requested variable is a tracer field

end subroutine get_values

!==============================================================
subroutine mass_wtd_vert_intg( arrayout, arrayin, state, x_dp )

  use physics_types,   only: physics_state
  use conditional_diag,only: PDEL,PDELDRY,NODP
  use physconst,       only: gravit
  use ppgrid,          only: pcols, pver

  real(r8),           intent(out) :: arrayout(:,:)
  real(r8),           intent(in)  :: arrayin (:,:)
  type(physics_state),intent(in)  :: state
  integer,            intent(in)  :: x_dp 

  real(r8) :: tmp(pcols,pver)

  integer :: ncol
  character(len=*),parameter :: subname = 'conditional_diag_main: mass_wtd_vert_intg'

  ncol = state%ncol

  ! Multiply by the appropriate dp (dry or wet)

  select case (x_dp)
  case (PDEL)
     tmp(1:ncol,:) = arrayin(1:ncol,:) * state%pdel(1:ncol,:)

  case (PDELDRY)
     tmp(1:ncol,:) = arrayin(1:ncol,:) * state%pdeldry(1:ncol,:)

  case default
     call endrun(subname//': unexpected value of x_dp')
  end select

  ! Vertical sum divided by gravit

  arrayout(1:ncol,1) = sum( tmp(1:ncol,:), 2 )/gravit

end subroutine mass_wtd_vert_intg 

!==============================================================
subroutine get_flags( metric, icnd, ncol, cnd_diag_info, flag )

  real(r8),              intent(in) :: metric(:,:)
  integer,               intent(in) :: icnd
  integer,               intent(in) :: ncol
  type(cnd_diag_info_t), intent(in) :: cnd_diag_info

  real(r8), intent(inout) :: flag(:,:)

  character(len=*),parameter :: subname = 'get_flags'

  flag(:,:) = OFF

  select case (cnd_diag_info%metric_cmpr_type(icnd))
  case (GT)

    where (metric(1:ncol,:).gt.cnd_diag_info%metric_threshold(icnd))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (GE)

    where (metric(1:ncol,:).ge.cnd_diag_info%metric_threshold(icnd))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (LT)

    where (metric(1:ncol,:).lt.cnd_diag_info%metric_threshold(icnd))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (LE)

    where (metric(1:ncol,:).le.cnd_diag_info%metric_threshold(icnd))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (EQ)

    where ( abs(metric(1:ncol,:)-cnd_diag_info%metric_threshold(icnd)) &
            .le. cnd_diag_info%metric_tolerance(icnd)                  )
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case default
    call endrun(subname//': unknown cnd_diag_info%metric_cmpr_type')
  end select

end subroutine get_flags

end module conditional_diag_main
