module qneg_module

  use shr_kind_mod,        only: r8 => shr_kind_r8, CS => SHR_KIND_CS
  use perf_mod,            only: t_startf, t_stopf
  use cam_logfile,         only: iulog
  use cam_abortutils,      only: endrun
  use shr_sys_mod,         only: shr_sys_flush
  use cam_history_support, only: max_fieldname_len
  use ppgrid,              only: pcols
  use constituents,        only: pcnst, cnst_name

  implicit none
  private
  save

  ! Public interface.

  public :: qneg_readnl
  public :: qneg_init
  public :: qneg3_new
  public :: qneg4_new
  public :: qneg_print_summary

  ! Private module variables
  character(len=8) :: print_qneg_warn
  logical          :: log_warnings = .false.
  logical          :: collect_stats = .false.
  logical          :: timestep_reset = .false.

  real(r8), parameter :: tol = 1.e-12_r8
  real(r8), parameter :: worst_reset = 1.e35_r8

  ! Diagnostic field names
  integer, parameter               :: num_diag_fields = (2 * pcnst) + 1
  character(len=max_fieldname_len) :: diag_names(num_diag_fields) = ''
  logical             :: cnst_out_calc = .false.
  logical             :: cnst_outfld(num_diag_fields) = .false.

  ! Summary buffers
  integer, parameter  :: num3_bins = 24
  integer, parameter  :: num4_bins = 4
  character(len=CS)   :: qneg3_warn_labels(num3_bins) = ''
  character(len=CS)   :: qneg4_warn_labels(num4_bins) = ''
  integer             :: qneg3_warn_num(pcnst, num3_bins)   = 0
  integer             :: qneg4_warn_num(num4_bins)          = 0
  real(r8)            :: qneg3_warn_worst(pcnst, num3_bins) = worst_reset
  real(r8)            :: qneg4_warn_worst(num4_bins)        = worst_reset

  private             :: calc_cnst_out
  private             :: find_index3
  private             :: find_index4
  interface reset_stats
     module procedure reset_stats_scalar
     module procedure reset_stats_array
  end interface reset_stats

contains

  subroutine qneg_readnl(nlfile)
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_character, masterproc
    ! File containing namelist input.
    character(len=*), intent(in) :: nlfile

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: sub = 'qneg_readnl'

    namelist /qneg_nl/ print_qneg_warn

    print_qneg_warn = ''

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'qneg_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, qneg_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(sub // ':: ERROR reading namelist qneg_nl')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    call mpi_bcast(print_qneg_warn, len(print_qneg_warn), mpi_character, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: print_qneg_warn")

    select case(trim(print_qneg_warn))
    case('summary')
       collect_stats = .true.
       timestep_reset = .false.
    case('timestep')
       collect_stats = .true.
       timestep_reset = .true.
    case('off')
       collect_stats = .false.
       timestep_reset = .false.
    case default
       call endrun(sub//" FATAL: '"//trim(print_qneg_warn)//"' is not a valid value for print_qneg_warn")
    end select

    if (masterproc) then
       if (collect_stats) then
          if (timestep_reset) then
             write(iulog, *) sub, ": QNEG statistics will be collected and printed for each timestep"
          else
             write(iulog, *) sub, ": QNEG statistics will be collected and printed at the end of the run"
          end if
       else
          write(iulog, *) sub, ": QNEG statistics will not be collected"
       end if
    end if

  end subroutine qneg_readnl

  subroutine qneg_init()
    use cam_history,    only: addfld, horiz_only
    use constituents,   only: cnst_longname

    integer :: index

    do index = 1, pcnst
       diag_names(index) = trim(cnst_name(index))//'_qneg3'
       call addfld(diag_names(index), (/ 'lev' /), 'I', 'kg/kg',              &
            trim(cnst_longname(index))//' QNEG3 error (cell)')
       diag_names(pcnst+index) = trim(cnst_name(index))//'_qneg3_col'
       call addfld(diag_names(pcnst+index), horiz_only, 'I', 'kg/kg',         &
            trim(cnst_longname(index))//' QNEG3 error (column)')
    end do
    diag_names((2*pcnst) + 1) = 'qflux_exceeded'
    call addfld(diag_names((2*pcnst) + 1), horiz_only, 'I', 'kg/m^2/s',     &
         'qflux excess (QNEG4)')

  end subroutine qneg_init

  subroutine calc_cnst_out()
    use cam_history, only: hist_fld_active, history_initialized
    integer :: index

    if (history_initialized()) then
       ! to protect against routines that call qneg3 too early
       do index = 1, num_diag_fields
          cnst_outfld(index) = hist_fld_active(trim(diag_names(index)))
       end do
       cnst_out_calc = .true.
    end if

  end subroutine calc_cnst_out

  integer function find_index3(nam) result(index)
    ! Find a valid or new index for 'nam' entries
    character(len=*),  intent(in) :: nam

    integer                      :: i

    index = -1
    do i = 1, num3_bins
       if (trim(nam) == trim(qneg3_warn_labels(i))) then
          ! We found this entry, return its index
          index = i
          exit
       else if (len_trim(qneg3_warn_labels(i)) == 0) then
          ! We have run out of known entries, use a new one and reset its stats
          qneg3_warn_labels(i) = nam
          index = i
          call reset_stats(qneg3_warn_num(:, index), qneg3_warn_worst(:,index))
          exit
       end if
    end do
  end function find_index3

  integer function find_index4(nam) result(index)
    ! Find a valid or new index for 'nam' entries
    character(len=*),  intent(in) :: nam

    integer                      :: i

    index = -1
    do i = 1, num4_bins
       if (trim(nam) == trim(qneg4_warn_labels(i))) then
          ! We found this entry, return its index
          index = i
          exit
       else if (len_trim(qneg4_warn_labels(i)) == 0) then
          ! We have run out of known entries, use a new one and reset its stats
          qneg4_warn_labels(i) = nam
          index = i
          call reset_stats(qneg4_warn_num(index), qneg4_warn_worst(index))
          exit
       end if
    end do
  end function find_index4

  subroutine qneg3_new (subnam, idx, ncol, ncold, lver, lconst_beg, &
       lconst_end, qmin, q)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Check moisture and tracers for minimum value, reset any below
    ! minimum value to minimum value and return information to allow
    ! warning message to be printed. The global average is NOT preserved.
    !
    ! Method:
    ! <Describe the algorithm(s) used in the routine.>
    ! <Also include any applicable external references.>
    !
    ! Author: J. Rosinski
    !
    !-----------------------------------------------------------------------
    use cam_history, only: outfld

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    character(len=*), intent(in) :: subnam ! name of calling routine

    integer, intent(in) :: idx          ! chunk/latitude index
    integer, intent(in) :: ncol         ! number of atmospheric columns
    integer, intent(in) :: ncold        ! declared number of atmospheric columns
    integer, intent(in) :: lver         ! number of vertical levels in column
    integer, intent(in) :: lconst_beg   ! beginning constituent
    integer, intent(in) :: lconst_end   ! ending    constituent

    real(r8), intent(in) :: qmin(lconst_beg:lconst_end)      ! Global minimum constituent concentration

    !
    ! Input/Output arguments
    !
    real(r8), intent(inout) :: q(ncold,lver,lconst_beg:lconst_end) ! moisture/tracer field
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: nvals            ! number of values found < qmin
    integer  :: i, k             ! longitude, level indices
    integer  :: index            ! For storing stats
    integer  :: m                ! constituent index
    integer  :: iw,kw            ! i,k indices of worst violator

    logical  :: found            ! true => at least 1 minimum violator found

    real(r8) :: badvals(ncold, lver) ! Collector for outfld calls
    real(r8) :: badcols(ncold)  ! Column sum for outfld
    real(r8) :: worst           ! biggest violator
    !
    !-----------------------------------------------------------------------
    !

    call t_startf ('qneg3')
    ! The first time we call this, we need to determine whether to call outfld
    if (.not. cnst_out_calc) then
       call calc_cnst_out()
    end if

    if (collect_stats) then
       index = find_index3(trim(subnam))
    else
       index = -1
    end if

    do m = lconst_beg, lconst_end
       nvals = 0
       found = .false.
       worst = worst_reset
       badvals(:,:) = 0.0_r8
       iw = -1
       kw = -1
       !
       ! Test all field values for being less than minimum value. Set q = qmin
       ! for all such points. Trace offenders and identify worst one.
       !
       do k = 1, lver
          do i = 1, ncol
             if (q(i,k,m) < qmin(m)) then
                found = .true.
                nvals = nvals + 1
                badvals(i, k) = q(i, k, m)
                if (index > 0) then
                   qneg3_warn_num(m, index) = qneg3_warn_num(m, index) + 1
                end if
                if (q(i,k,m) < worst) then
                   worst = q(i,k,m)
                   iw = i
                   kw = k
                   if (index > 0) then
                      qneg3_warn_worst(m, index) = worst
                   end if
                end if
                q(i,k,m) = qmin(m)
             end if
          end do
       end do
       ! Maybe output bad values
       if ((cnst_outfld(m)) .and. (worst < worst_reset)) then
          call outfld(trim(diag_names(m)), badvals, pcols, idx)
       end if
       if ((cnst_outfld(pcnst+m)) .and. (worst < worst_reset)) then
          do i = 1, pcols
             badcols(i) = SUM(badvals(i,:))
          end do
          call outfld(trim(diag_names(pcnst+m)), badcols, pcols, idx)
       end if
    end do
    call t_stopf ('qneg3')

  end subroutine qneg3_new

  subroutine qneg4_new (subnam, lchnk, ncol, ztodt,                            &
       qbot, srfrpdel, shflx, lhflx, qflx)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Check if moisture flux into the ground is exceeding the total
    ! moisture content of the lowest model layer (creating negative moisture
    ! values).  If so, then subtract the excess from the moisture and
    ! latent heat fluxes and add it to the sensible heat flux.
    !
    ! Method:
    ! <Describe the algorithm(s) used in the routine.>
    ! <Also include any applicable external references.>
    !
    ! Author: J. Olson
    !
    !-----------------------------------------------------------------------
!    use phys_grid,    only: get_lat_p, get_lon_p
    use physconst,    only: gravit, latvap
    use constituents, only: qmin
    use cam_history,  only: outfld

    !
    ! Input arguments
    !
    character(len=*), intent(in) :: subnam   ! name of calling routine
    !
    integer, intent(in) :: lchnk             ! chunk index
    integer, intent(in) :: ncol              ! number of atmospheric columns
    !
    real(r8), intent(in) :: ztodt            ! two times model timestep (2 delta-t)
    real(r8), intent(in) :: qbot(ncol,pcnst) ! moisture at lowest model level
    real(r8), intent(in) :: srfrpdel(ncol)   ! 1./(pint(K+1)-pint(K))
    !
    ! Input/Output arguments
    !
    real(r8), intent(inout) :: shflx(ncol)   ! Surface sensible heat flux (J/m2/s)
    real(r8), intent(inout) :: lhflx(ncol)   ! Surface latent   heat flux (J/m2/s)
    real(r8), intent(inout) :: qflx (ncol,pcnst) ! surface water flux (kg/m^2/s)
    !
    !---------------------------Local workspace-----------------------------
    !
    integer :: i                 ! column index
    integer :: iw                ! i index of worst violator
    integer :: index             ! caller bin index
    !
    real(r8):: worst             ! biggest violator
    real(r8):: excess(ncol)     ! Excess downward sfc latent heat flux
    !
    !-----------------------------------------------------------------------

    call t_startf ('qneg4')
    ! The first time we call this, we need to determine whether to call outfld
    if (.not. cnst_out_calc) then
       call calc_cnst_out()
    end if

    if (collect_stats) then
       index = find_index4(trim(subnam))
    else
       index = -1
    end if

    !
    ! Compute excess downward (negative) q flux compared to a theoretical
    ! maximum downward q flux.  The theoretical max is based upon the
    ! given moisture content of lowest level of the model atmosphere.
    !
    worst = worst_reset
    do i = 1, ncol
       excess(i) = qflx(i,1) - (qmin(1) - qbot(i,1))/(ztodt*gravit*srfrpdel(i))
       !
       ! If there is an excess downward (negative) q flux, then subtract
       ! excess from "qflx" and "lhflx" and add to "shflx".
       !
       if (excess(i) < 0._r8) then
          if (excess(i) < worst) then
             iw = i
             worst = excess(i)
          end if
          qflx (i,1) = qflx (i,1) - excess(i)
          lhflx(i) = lhflx(i) - excess(i)*latvap
          shflx(i) = shflx(i) + excess(i)*latvap
          if (index > 0) then
             qneg4_warn_num(index) = qneg4_warn_num(index) + 1
          end if
       end if
    end do
    ! Maybe output bad values
    if ((cnst_outfld((2*pcnst)+1)) .and. (worst < worst_reset)) then
       do i = 1, ncol
          if (excess(i) > 0.0_r8) then
             excess(i) = 0.0_r8
          end if
       end do
       call outfld(trim(diag_names((2*pcnst)+1)), excess(1:ncol), ncol, lchnk)
    end if
    call t_stopf ('qneg4')

  end subroutine qneg4_new

  subroutine qneg_print_summary(end_of_run)
    use spmd_utils, only: mpicom, masterprocid, masterproc
    use spmd_utils, only: MPI_MIN, MPI_SUM, MPI_INTEGER, MPI_REAL8

    logical, intent(in) :: end_of_run

    integer             :: global_warn_num(pcnst)
    real(r8)            :: global_warn_worst(pcnst)
    integer             :: index, m
    integer             :: ierr

    if (collect_stats) then
       if (timestep_reset .or. end_of_run) then
          do index = 1, num3_bins
             ! QNEG3
             call reset_stats(global_warn_num(:), global_warn_worst(:))
             call MPI_REDUCE(qneg3_warn_num(:, index), global_warn_num(:),    &
                  pcnst, MPI_INTEGER, MPI_SUM, masterprocid, mpicom, ierr)
             call MPI_REDUCE(qneg3_warn_worst(:, index), global_warn_worst(:),&
                  pcnst, MPI_REAL8, MPI_MIN, masterprocid, mpicom, ierr)
             if (masterproc) then
                do m = 1, pcnst
                   if ( (global_warn_num(m) > 0) .and.                        &
                        (abs(global_warn_worst(m)) > tol)) then
                      write(iulog, 9100) trim(qneg3_warn_labels(index)),      &
                           trim(cnst_name(m)), global_warn_num(m),            &
                           global_warn_worst(m)
                   end if
                   call shr_sys_flush(iulog)
                end do
             end if
             call reset_stats(qneg3_warn_num(:,index), qneg3_warn_worst(:,index))
          end do
          do index = 1, num4_bins
             ! QNEG4
             call reset_stats(qneg4_warn_num(index), qneg4_warn_worst(index))
             call reset_stats(global_warn_num(1), global_warn_worst(1))
             call MPI_REDUCE(qneg4_warn_num(index), global_warn_num(1),       &
                  1, MPI_INTEGER, MPI_SUM, masterprocid, mpicom, ierr)
             call MPI_REDUCE(qneg4_warn_worst(index), global_warn_worst(1),   &
                  1, MPI_REAL8, MPI_MIN, masterprocid, mpicom, ierr)
             if (masterproc) then
                if ( (global_warn_num(1) > 0) .and.                           &
                     (abs(global_warn_worst(1)) > tol)) then
                   write(iulog, 9101) trim(qneg4_warn_labels(index)),          &
                        global_warn_num(1), global_warn_worst(1)
                end if
                call shr_sys_flush(iulog)
             end if
             call reset_stats(qneg4_warn_num(index), qneg4_warn_worst(index))
          end do
       end if
    end if

    return
9100 format(' QNEG3 from ', a, ':', a, &
         ' Min. mixing ratio violated at ', i9, ' points. Worst = ', e10.1)
9101 format(' QNEG4 from ',a,': moisture flux exceeded at', &
          i9, ' points. Worst = ', e10.1)
  end subroutine qneg_print_summary

  subroutine reset_stats_array(num_array, worst_array)
    ! Private routine to reset statistics
    integer,  intent(inout) :: num_array(:)
    real(r8), intent(inout) :: worst_array(:)

    num_array(:)    = 0
    worst_array(:)  = worst_reset
  end subroutine reset_stats_array

  subroutine reset_stats_scalar(num, worst)
    ! Private routine to reset statistics
    integer,  intent(inout) :: num
    real(r8), intent(inout) :: worst

    num    = 0
    worst  = worst_reset
  end subroutine reset_stats_scalar

end module qneg_module
