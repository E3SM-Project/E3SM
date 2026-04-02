module eam_derived
  !-------------------------------------------------------------------------------------------
  !
  ! EAM derived field wrapper.
  !
  ! Thin integration layer between EAM's physics state / history system and
  ! the shared shr_derived_mod routines. All expression parsing and evaluation
  ! math is delegated to the shared module.
  !
  ! Provides three capabilities:
  !   1. Derived fields: user-defined expressions combining state variables,
  !      constituents, physics buffer fields, and numeric constants.
  !      Supports both 2D (horizontal-only) and 3D (horizontal+vertical) fields.
  !      If all field operands are 2D, the result is 2D. If any operand is 3D,
  !      the result is 3D (2D operands are broadcast across all levels).
  !   2. Automatic tendencies: cache previous timestep values and output
  !      time derivatives for any field (state, constituent, pbuf, or derived).
  !   3. Stage-aware tendencies: split tendencies into physics and dynamics
  !      contributions. Physics tendency = change during all physics schemes.
  !      Dynamics tendency = total tendency - physics tendency.
  !
  ! Configuration via namelist (eam_derived_nl):
  !   derived_fld_defs  - expression definitions, e.g. "TOTAL_WATER=Q+CLDICE+CLDLIQ+RAINQM"
  !   tend_flds         - field names for automatic tendency output
  !   tend_stages       - which stage tendencies to output: 'total', 'phys', 'dyn'
  !
  ! Expression syntax:
  !   OUTPUT_NAME=OPERAND1 op OPERAND2 op ...
  !   where op is +, -, *, / and each operand is a field name or numeric constant.
  !   Evaluation is strict left-to-right (no operator precedence).
  !   Units are not inferred automatically; expressions must be dimensionally
  !   consistent and derived outputs are currently registered with unit string
  !   'derived'.
  !   Examples:
  !     TOTAL_WATER=Q+CLDICE+CLDLIQ+RAINQM  (3D: sum of constituents)
  !     HGTsfc=PHIS/9.80616                   (2D: surface geopotential / constant)
  !     TOTAL_WATER_g=TOTAL_WATER*1000.0       (chaining: uses earlier derived field)
  !
  ! Chaining: definitions are processed in order; earlier derived fields are
  ! cached and available as inputs to later definitions.
  !
  ! Tendencies: for each field in tend_flds, outputs d{NAME}_dt = (curr - prev) / dt
  ! each physics timestep. The first timestep stores the initial value and outputs zero.
  !
  ! Stage-aware tendencies: when tend_stages includes 'phys' or 'dyn':
  !   d{NAME}_dt_phys = (state_after_physics - state_before_physics) / dt
  !   d{NAME}_dt_dyn  = d{NAME}_dt - d{NAME}_dt_phys
  !   d{NAME}_dt      = total tendency (existing behavior)
  !
  ! Usage:
  !   call eam_derived_readnl(nlfile)            ! during namelist reading
  !   call eam_derived_register()                ! during phys_init / phys_register
  !   call eam_derived_stage(state, pbuf2d)      ! before physics (after dynamics output)
  !   call eam_derived_write(state, pbuf2d)      ! after all physics
  !
  !-------------------------------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use shr_derived_mod, only: shr_derived_expr_t, shr_derived_operand_t, &
                              shr_derived_parse, shr_derived_eval, &
                              shr_derived_tend, shr_derived_is_number, &
                              shr_derived_max_operands, shr_derived_max_namelen, &
                              shr_derived_max_deflen
  use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk
  use cam_logfile,     only: iulog
  use cam_abortutils,  only: endrun
  use spmd_utils,      only: masterproc

  implicit none
  private
  save

  public :: eam_derived_readnl
  public :: eam_derived_register
  public :: eam_derived_stage     ! snapshot state before physics (for stage tendencies)
  public :: eam_derived_write

  ! Parameters
  integer, parameter :: max_derived_flds = 50
  integer, parameter :: max_tend_flds    = 100
  integer, parameter :: max_name_len     = shr_derived_max_namelen
  integer, parameter :: max_def_len      = shr_derived_max_deflen

  ! Known 3D state variable names
  integer, parameter :: n_known_3d_state = 6
  character(len=max_name_len), parameter :: known_3d_state_vars(n_known_3d_state) = &
       (/ 'T     ', 'U     ', 'V     ', 'OMEGA ', 'Z3    ', 'Q     ' /)

  ! Known 2D state variable names
  integer, parameter :: n_known_2d_state = 2
  character(len=max_name_len), parameter :: known_2d_state_vars(n_known_2d_state) = &
       (/ 'PS   ', 'PHIS ' /)

  ! Stage tendency parameters
  integer, parameter :: max_stages = 10
  integer, parameter :: max_stage_len = 16

  ! Namelist variables
  character(len=max_def_len)   :: derived_fld_defs(max_derived_flds)
  character(len=max_name_len)  :: tend_flds(max_tend_flds)
  character(len=max_stage_len) :: tend_stages(max_stages)

  ! Parsed state
  integer :: n_derived = 0
  type(shr_derived_expr_t) :: expressions(max_derived_flds)
  logical :: expr_is_2d(max_derived_flds) = .false.  ! true if result is 2D (horiz_only)

  integer :: n_tend = 0
  character(len=max_name_len) :: tend_names(max_tend_flds)
  logical :: tend_is_2d(max_tend_flds) = .false.     ! true if tendency source is 2D

  ! Derived field cache for chaining (allocated during register)
  ! For 2D fields, only level 1 is used.
  real(r8), allocatable :: derived_cache(:,:,:,:)   ! (pcols, pver, n_derived, begchunk:endchunk)

  ! Tendency tracking (allocated during register)
  real(r8), allocatable :: tend_prev(:,:,:,:)        ! (pcols, pver, n_tend, begchunk:endchunk)
  logical :: tend_initialized = .false.

  ! Stage-aware tendency tracking
  integer :: n_stages = 0
  logical :: do_stage_phys = .false.   ! output d{NAME}_dt_phys
  logical :: do_stage_dyn  = .false.   ! output d{NAME}_dt_dyn
  real(r8), allocatable :: phys_snap(:,:,:,:)  ! (pcols, pver, n_tend, begchunk:endchunk)
  logical :: phys_snap_valid = .false.

  ! Flags
  logical :: has_derived = .false.
  logical :: has_tend    = .false.
  logical :: module_is_initialized = .false.

contains

  !============================================================================
  subroutine eam_derived_readnl(nlfile)
    !--------------------------------------------------------------------------
    ! Read the eam_derived_nl namelist group.
    !--------------------------------------------------------------------------
    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: mpicom, masterprocid, mpi_character

    character(len=*), intent(in) :: nlfile

    integer :: unitn, ierr, i

    namelist /eam_derived_nl/ derived_fld_defs, tend_flds, tend_stages

    ! Initialize defaults
    derived_fld_defs(:) = ''
    tend_flds(:) = ''
    tend_stages(:) = ''
    tend_stages(1) = 'total'   ! default: only output total tendency

    if (masterproc) then
      unitn = getunit()
      open(unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'eam_derived_nl', status=ierr)
      if (ierr == 0) then
        read(unitn, eam_derived_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun('eam_derived_readnl: ERROR reading namelist eam_derived_nl')
        end if
      end if
      close(unitn)
      call freeunit(unitn)
    end if

    ! Broadcast to all processors
    call mpi_bcast(derived_fld_defs, max_def_len*max_derived_flds, mpi_character, &
         masterprocid, mpicom, ierr)
    call mpi_bcast(tend_flds, max_name_len*max_tend_flds, mpi_character, &
         masterprocid, mpicom, ierr)
    call mpi_bcast(tend_stages, max_stage_len*max_stages, mpi_character, &
         masterprocid, mpicom, ierr)

    ! Parse derived field definitions
    n_derived = 0
    do i = 1, max_derived_flds
      if (len_trim(derived_fld_defs(i)) == 0) exit
      n_derived = n_derived + 1
      call shr_derived_parse(derived_fld_defs(i), expressions(n_derived), ierr)
      if (ierr /= 0) then
        call endrun('eam_derived_readnl: failed to parse definition: '// &
             trim(derived_fld_defs(i)))
      end if
    end do
    has_derived = (n_derived > 0)

    ! Parse tendency field names
    n_tend = 0
    do i = 1, max_tend_flds
      if (len_trim(tend_flds(i)) == 0) exit
      n_tend = n_tend + 1
      tend_names(n_tend) = adjustl(tend_flds(i))
    end do
    has_tend = (n_tend > 0)

    ! Parse tendency stages
    n_stages = 0
    do_stage_phys = .false.
    do_stage_dyn  = .false.
    do i = 1, max_stages
      if (len_trim(tend_stages(i)) == 0) exit
      n_stages = n_stages + 1
      select case (trim(tend_stages(i)))
      case ('phys')
        do_stage_phys = .true.
      case ('dyn')
        do_stage_dyn = .true.
      case ('total')
        ! Default behavior, always active when has_tend
      case default
        call endrun('eam_derived_readnl: unknown tend_stage "'// &
             trim(tend_stages(i))//'". Must be total, phys, or dyn.')
      end select
    end do

    if (masterproc) then
      if (has_derived) then
        write(iulog,*) 'eam_derived_readnl: ', n_derived, ' derived field(s) defined'
        do i = 1, n_derived
          write(iulog,*) '  ', trim(expressions(i)%output_name), ' = ', &
               trim(expressions(i)%long_name)
        end do
      end if
      if (has_tend) then
        write(iulog,*) 'eam_derived_readnl: ', n_tend, ' tendency field(s) requested'
        do i = 1, n_tend
          write(iulog,*) '  d', trim(tend_names(i)), '_dt'
        end do
        write(iulog,*) '  stages:', (trim(tend_stages(i))//' ', i=1,n_stages)
        if (do_stage_phys) write(iulog,*) '  physics stage tendencies enabled'
        if (do_stage_dyn)  write(iulog,*) '  dynamics stage tendencies enabled'
      end if
    end if

  end subroutine eam_derived_readnl

  !============================================================================
  subroutine eam_derived_register()
    !--------------------------------------------------------------------------
    ! Register derived output fields and tendency fields with the history
    ! system. Must be called during phys_register, after other addfld calls.
    !
    ! Determines 2D vs 3D dimensionality for each expression:
    !   - If ALL field operands are 2D -> result is 2D (horiz_only)
    !   - If ANY field operand is 3D  -> result is 3D (on 'lev')
    !--------------------------------------------------------------------------
    use cam_history, only: addfld, horiz_only

    integer :: i, n, ndims
    character(len=max_name_len) :: fname
    character(len=max_def_len)  :: lname
    logical :: all_2d

    if (.not. has_derived .and. .not. has_tend) return

    ! Validate and register derived fields, determine dimensionality
    do i = 1, n_derived
      all_2d = .true.

      do n = 1, expressions(i)%n_operands
        if (.not. expressions(i)%operands(n)%is_constant) then
          call validate_field_name(expressions(i)%operands(n)%field_name, i)
          ndims = get_field_ndims(expressions(i)%operands(n)%field_name, i)
          if (ndims > 1) all_2d = .false.
        end if
      end do

      expr_is_2d(i) = all_2d

      if (all_2d) then
        call addfld(trim(expressions(i)%output_name), horiz_only, 'A', &
             'derived', trim(expressions(i)%long_name))
      else
        call addfld(trim(expressions(i)%output_name), (/ 'lev' /), 'A', &
             'derived', trim(expressions(i)%long_name))
      end if
    end do

    ! Validate and register tendency fields
    do i = 1, n_tend
      call validate_tend_field(tend_names(i))
      tend_is_2d(i) = is_field_2d(tend_names(i))

      ! Total tendency (always registered when tend_flds is set)
      fname = 'd' // trim(tend_names(i)) // '_dt'
      lname = 'd(' // trim(tend_names(i)) // ')/dt'
      if (tend_is_2d(i)) then
        call addfld(trim(fname), horiz_only, 'A', '/s', trim(lname))
      else
        call addfld(trim(fname), (/ 'lev' /), 'A', '/s', trim(lname))
      end if

      ! Physics stage tendency
      if (do_stage_phys) then
        fname = 'd' // trim(tend_names(i)) // '_dt_phys'
        lname = 'd(' // trim(tend_names(i)) // ')/dt physics'
        if (tend_is_2d(i)) then
          call addfld(trim(fname), horiz_only, 'A', '/s', trim(lname))
        else
          call addfld(trim(fname), (/ 'lev' /), 'A', '/s', trim(lname))
        end if
      end if

      ! Dynamics stage tendency
      if (do_stage_dyn) then
        fname = 'd' // trim(tend_names(i)) // '_dt_dyn'
        lname = 'd(' // trim(tend_names(i)) // ')/dt dynamics'
        if (tend_is_2d(i)) then
          call addfld(trim(fname), horiz_only, 'A', '/s', trim(lname))
        else
          call addfld(trim(fname), (/ 'lev' /), 'A', '/s', trim(lname))
        end if
      end if
    end do

    ! Allocate derived field cache for chaining
    ! For 2D fields, only (:,1,:,:) is used, but allocating uniform shape is simpler
    if (has_derived) then
      allocate(derived_cache(pcols, pver, n_derived, begchunk:endchunk))
      derived_cache(:,:,:,:) = 0.0_r8
    end if

    ! Allocate tendency cache
    if (has_tend) then
      allocate(tend_prev(pcols, pver, n_tend, begchunk:endchunk))
      tend_prev(:,:,:,:) = 0.0_r8
      tend_initialized = .false.

      ! Allocate physics stage snapshot cache
      if (do_stage_phys .or. do_stage_dyn) then
        allocate(phys_snap(pcols, pver, n_tend, begchunk:endchunk))
        phys_snap(:,:,:,:) = 0.0_r8
        phys_snap_valid = .false.
      end if
    end if

    module_is_initialized = .true.

    if (masterproc) then
      write(iulog,*) 'eam_derived_register: registration complete'
      do i = 1, n_derived
        if (expr_is_2d(i)) then
          write(iulog,*) '  ', trim(expressions(i)%output_name), ' (2D)'
        else
          write(iulog,*) '  ', trim(expressions(i)%output_name), ' (3D)'
        end if
      end do
    end if

  end subroutine eam_derived_register

  !============================================================================
  subroutine eam_derived_stage(state, pbuf2d)
    !--------------------------------------------------------------------------
    ! Snapshot the current state of all tendency fields BEFORE physics begins.
    ! Called from cam_comp.F90 after stepon_run1 (dynamics output loaded into
    ! phys_state) and before phys_run1 (physics starts modifying state).
    !
    ! This snapshot is used later in eam_derived_write to compute:
    !   d{NAME}_dt_phys = (state_after_phys - phys_snap) / dt
    !   d{NAME}_dt_dyn  = d{NAME}_dt - d{NAME}_dt_phys
    !--------------------------------------------------------------------------
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc

    type(physics_state),       intent(in) :: state
    type(physics_buffer_desc), pointer    :: pbuf2d(:,:)

    real(r8) :: curr_field(pcols, pver)
    integer  :: i, ncol, lchnk, nlev

    if (.not. has_tend) return
    if (.not. do_stage_phys .and. .not. do_stage_dyn) return

    ncol  = state%ncol
    lchnk = state%lchnk

    ! Snapshot each tendency field at this stage boundary
    do i = 1, n_tend
      if (tend_is_2d(i)) then
        nlev = 1
      else
        nlev = pver
      end if

      ! Use get_field_or_derived so derived fields in tend_flds work
      ! (they read the cache from the previous timestep's write step).
      call get_field_or_derived(state, pbuf2d, tend_names(i), &
           curr_field, ncol, lchnk, nlev)
      phys_snap(1:ncol, 1:nlev, i, lchnk) = curr_field(1:ncol, 1:nlev)
    end do

    phys_snap_valid = .true.

  end subroutine eam_derived_stage

  !============================================================================
  subroutine eam_derived_write(state, pbuf2d)
    !--------------------------------------------------------------------------
    ! Compute and output derived fields and tendencies for the current chunk.
    ! Called each timestep from physpkg.F90.
    !--------------------------------------------------------------------------
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    use cam_history,    only: outfld
    use time_manager,   only: get_step_size

    type(physics_state),       intent(in) :: state
    type(physics_buffer_desc), pointer    :: pbuf2d(:,:)

    real(r8) :: tmp_fields(pcols, pver, shr_derived_max_operands)
    real(r8) :: result(pcols, pver)
    real(r8) :: curr_field(pcols, pver)
    real(r8) :: tend_field(pcols, pver)
    real(r8) :: phys_tend_field(pcols, pver)
    integer  :: i, n, ncol, lchnk, nlev
    real(r8) :: dt
    character(len=max_name_len+16) :: fname

    if (.not. has_derived .and. .not. has_tend) return

    ncol  = state%ncol
    lchnk = state%lchnk
    dt    = real(get_step_size(), r8)

    ! --- Step 1: Compute derived fields in definition order (enables chaining) ---
    do i = 1, n_derived
      if (expr_is_2d(i)) then
        nlev = 1
      else
        nlev = pver
      end if

      ! Load field data for each operand
      do n = 1, expressions(i)%n_operands
        if (.not. expressions(i)%operands(n)%is_constant) then
          call get_field_or_derived(state, pbuf2d, &
               expressions(i)%operands(n)%field_name, &
               tmp_fields(:,:,n), ncol, lchnk, nlev)
        end if
      end do

      ! Evaluate expression
      call shr_derived_eval(expressions(i), tmp_fields, ncol, nlev, pcols, result)

      ! Cache result for chaining and tendencies
      derived_cache(1:ncol, 1:nlev, i, lchnk) = result(1:ncol, 1:nlev)

      ! Output to history
      call outfld(trim(expressions(i)%output_name), result, pcols, lchnk)
    end do

    ! --- Step 2: Compute tendencies ---
    if (has_tend) then
      do i = 1, n_tend
        if (tend_is_2d(i)) then
          nlev = 1
        else
          nlev = pver
        end if

        ! Load current field value (after all physics)
        call get_field_or_derived(state, pbuf2d, tend_names(i), &
             curr_field, ncol, lchnk, nlev)

        ! --- Total tendency: d{NAME}_dt ---
        if (tend_initialized) then
          call shr_derived_tend(curr_field, tend_prev(:,:,i,lchnk), &
               ncol, nlev, pcols, dt, tend_field)
        else
          tend_field(1:ncol, 1:nlev) = 0.0_r8
        end if

        fname = 'd' // trim(tend_names(i)) // '_dt'
        call outfld(trim(fname), tend_field, pcols, lchnk)

        ! --- Physics tendency: d{NAME}_dt_phys ---
        if (do_stage_phys .and. phys_snap_valid) then
          call shr_derived_tend(curr_field, phys_snap(:,:,i,lchnk), &
               ncol, nlev, pcols, dt, phys_tend_field)

          fname = 'd' // trim(tend_names(i)) // '_dt_phys'
          call outfld(trim(fname), phys_tend_field, pcols, lchnk)
        end if

        ! --- Dynamics tendency: d{NAME}_dt_dyn = total - phys ---
        if (do_stage_dyn .and. tend_initialized .and. phys_snap_valid) then
          ! dyn = total - phys
          ! total = (curr - prev_timestep) / dt
          ! phys  = (curr - before_phys) / dt
          ! dyn   = (before_phys - prev_timestep) / dt
          call shr_derived_tend(phys_snap(:,:,i,lchnk), tend_prev(:,:,i,lchnk), &
               ncol, nlev, pcols, dt, tend_field)

          fname = 'd' // trim(tend_names(i)) // '_dt_dyn'
          call outfld(trim(fname), tend_field, pcols, lchnk)
        end if

        ! Store current value for next timestep
        tend_prev(1:ncol, 1:nlev, i, lchnk) = curr_field(1:ncol, 1:nlev)
      end do

      if (.not. tend_initialized) tend_initialized = .true.
    end if

  end subroutine eam_derived_write

  !============================================================================
  ! Private helper routines
  !============================================================================

  subroutine get_field_or_derived(state, pbuf2d, fname, field_out, ncol, lchnk, nlev)
    !--------------------------------------------------------------------------
    ! Look up a field by name. Checks derived field cache first (for chaining),
    ! then falls back to state variables, constituents, and physics buffer.
    !
    ! For 3D expressions (nlev=pver): if the source field is 2D, broadcasts
    ! it across all vertical levels.
    ! For 2D expressions (nlev=1): returns the 2D field in field_out(:,1).
    !--------------------------------------------------------------------------
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc

    type(physics_state),       intent(in)  :: state
    type(physics_buffer_desc), pointer     :: pbuf2d(:,:)
    character(len=*),          intent(in)  :: fname
    real(r8),                  intent(out) :: field_out(pcols, pver)
    integer,                   intent(in)  :: ncol
    integer,                   intent(in)  :: lchnk
    integer,                   intent(in)  :: nlev  ! target: 1 for 2D, pver for 3D

    integer :: i, k

    field_out(:,:) = 0.0_r8

    ! Check derived field cache first (for chaining)
    do i = 1, n_derived
      if (trim(fname) == trim(expressions(i)%output_name)) then
        if (expr_is_2d(i) .and. nlev > 1) then
          ! Source is 2D, target is 3D: broadcast across levels
          do k = 1, nlev
            field_out(1:ncol, k) = derived_cache(1:ncol, 1, i, lchnk)
          end do
        else
          field_out(1:ncol, 1:nlev) = derived_cache(1:ncol, 1:nlev, i, lchnk)
        end if
        return
      end if
    end do

    ! Fall back to state/constituent/pbuf lookup
    call get_field(state, pbuf2d, fname, field_out, ncol, nlev)

  end subroutine get_field_or_derived

  !============================================================================
  subroutine get_field(state, pbuf2d, fname, field_out, ncol, nlev)
    !--------------------------------------------------------------------------
    ! Look up a field name and extract field data.
    ! Supports 2D and 3D state variables, constituents, and physics buffer fields.
    !
    ! For 3D target (nlev=pver): if source is 2D, broadcasts across levels.
    ! For 2D target (nlev=1): returns the 2D value in field_out(:,1).
    !--------------------------------------------------------------------------
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field, &
                              pbuf_get_chunk, pbuf_get_field_ndims
    use constituents,   only: cnst_get_ind

    type(physics_state),       intent(in)  :: state
    type(physics_buffer_desc), pointer     :: pbuf2d(:,:)
    character(len=*),          intent(in)  :: fname
    real(r8),                  intent(out) :: field_out(pcols, pver)
    integer,                   intent(in)  :: ncol
    integer,                   intent(in)  :: nlev  ! target: 1 for 2D, pver for 3D

    integer :: idx, pbuf_idx, k, pbuf_ndims
    character(len=max_name_len) :: uname
    real(r8), pointer :: pbuf_fld_3d(:,:)
    real(r8), pointer :: pbuf_fld_2d(:)
    type(physics_buffer_desc), pointer :: pbuf_chunk(:)

    uname = adjustl(fname)
    field_out(:,:) = 0.0_r8

    ! Check 3D state variables
    select case (trim(uname))
    case ('T')
      field_out(1:ncol, 1:nlev) = state%t(1:ncol, 1:nlev)
      return
    case ('U')
      field_out(1:ncol, 1:nlev) = state%u(1:ncol, 1:nlev)
      return
    case ('V')
      field_out(1:ncol, 1:nlev) = state%v(1:ncol, 1:nlev)
      return
    case ('OMEGA')
      field_out(1:ncol, 1:nlev) = state%omega(1:ncol, 1:nlev)
      return
    case ('Z3')
      field_out(1:ncol, 1:nlev) = state%zm(1:ncol, 1:nlev)
      return
    case ('Q')
      field_out(1:ncol, 1:nlev) = state%q(1:ncol, 1:nlev, 1)
      return
    end select

    ! Check 2D state variables
    select case (trim(uname))
    case ('PS')
      if (nlev == 1) then
        field_out(1:ncol, 1) = state%ps(1:ncol)
      else
        do k = 1, nlev
          field_out(1:ncol, k) = state%ps(1:ncol)
        end do
      end if
      return
    case ('PHIS')
      if (nlev == 1) then
        field_out(1:ncol, 1) = state%phis(1:ncol)
      else
        do k = 1, nlev
          field_out(1:ncol, k) = state%phis(1:ncol)
        end do
      end if
      return
    end select

    ! Try constituent lookup (always 3D)
    call cnst_get_ind(trim(uname), idx, abrtf=.false.)
    if (idx > 0) then
      field_out(1:ncol, 1:nlev) = state%q(1:ncol, 1:nlev, idx)
      return
    end if

    ! Try physics buffer lookup (2D or 3D)
    pbuf_idx = pbuf_get_index(trim(uname), errcode=idx)
    if (pbuf_idx > 0) then
      pbuf_chunk => pbuf_get_chunk(pbuf2d, state%lchnk)
      pbuf_ndims = pbuf_get_field_ndims(pbuf_idx)

      if (pbuf_ndims == 1) then
        ! 2D pbuf field (horizontal only)
        call pbuf_get_field(pbuf_chunk, pbuf_idx, pbuf_fld_2d)
        if (nlev == 1) then
          field_out(1:ncol, 1) = pbuf_fld_2d(1:ncol)
        else
          ! Broadcast across levels
          do k = 1, nlev
            field_out(1:ncol, k) = pbuf_fld_2d(1:ncol)
          end do
        end if
      else
        ! 3D pbuf field (horizontal + vertical)
        call pbuf_get_field(pbuf_chunk, pbuf_idx, pbuf_fld_3d)
        field_out(1:ncol, 1:nlev) = pbuf_fld_3d(1:ncol, 1:nlev)
      end if
      return
    end if

    call endrun('eam_derived: get_field: unknown field "'//trim(uname)// &
         '". Must be a state variable (T,U,V,OMEGA,Z3,Q,PS,PHIS), '// &
         'constituent, or physics buffer field.')

  end subroutine get_field

  !============================================================================
  function get_field_ndims(fname, def_index) result(ndims)
    !--------------------------------------------------------------------------
    ! Return the number of spatial dimensions for a field:
    !   1 = horizontal only (2D: PS, PHIS, 1D pbuf fields)
    !   2 = horizontal + vertical (3D: T, U, V, Q, constituents, 2D+ pbuf fields)
    ! For derived fields, returns 1 if expr_is_2d, else 2.
    !--------------------------------------------------------------------------
    use constituents,   only: cnst_get_ind
    use physics_buffer, only: pbuf_get_index, pbuf_get_field_ndims

    character(len=*), intent(in) :: fname
    integer,          intent(in) :: def_index  ! current definition index (for chaining)
    integer :: ndims

    integer :: k, idx
    character(len=max_name_len) :: uname

    uname = adjustl(fname)

    ! Check 3D state variables
    do k = 1, n_known_3d_state
      if (trim(uname) == trim(known_3d_state_vars(k))) then
        ndims = 2
        return
      end if
    end do

    ! Check 2D state variables
    do k = 1, n_known_2d_state
      if (trim(uname) == trim(known_2d_state_vars(k))) then
        ndims = 1
        return
      end if
    end do

    ! Check previously defined derived fields
    do k = 1, def_index - 1
      if (trim(uname) == trim(expressions(k)%output_name)) then
        if (expr_is_2d(k)) then
          ndims = 1
        else
          ndims = 2
        end if
        return
      end if
    end do

    ! Constituents are always 3D
    call cnst_get_ind(trim(uname), idx, abrtf=.false.)
    if (idx > 0) then
      ndims = 2
      return
    end if

    ! Physics buffer: query dimension count
    idx = pbuf_get_index(trim(uname), errcode=k)
    if (idx > 0) then
      ndims = pbuf_get_field_ndims(idx)
      return
    end if

    ! Unknown field - will be caught by validate_field_name
    ndims = 2

  end function get_field_ndims

  !============================================================================
  function is_field_2d(fname) result(is_2d)
    !--------------------------------------------------------------------------
    ! Check if a field is 2D (horizontal only).
    ! Used for tendency field dimensionality determination.
    !--------------------------------------------------------------------------
    use constituents,   only: cnst_get_ind
    use physics_buffer, only: pbuf_get_index, pbuf_get_field_ndims

    character(len=*), intent(in) :: fname
    logical :: is_2d

    integer :: k, idx
    character(len=max_name_len) :: uname

    uname = adjustl(fname)
    is_2d = .false.

    ! Check 2D state variables
    do k = 1, n_known_2d_state
      if (trim(uname) == trim(known_2d_state_vars(k))) then
        is_2d = .true.
        return
      end if
    end do

    ! Check 3D state variables
    do k = 1, n_known_3d_state
      if (trim(uname) == trim(known_3d_state_vars(k))) then
        return  ! is_2d = .false.
      end if
    end do

    ! Check derived fields
    do k = 1, n_derived
      if (trim(uname) == trim(expressions(k)%output_name)) then
        is_2d = expr_is_2d(k)
        return
      end if
    end do

    ! Constituents are always 3D
    call cnst_get_ind(trim(uname), idx, abrtf=.false.)
    if (idx > 0) return  ! is_2d = .false.

    ! Physics buffer: check ndims
    idx = pbuf_get_index(trim(uname), errcode=k)
    if (idx > 0) then
      is_2d = (pbuf_get_field_ndims(idx) == 1)
      return
    end if

  end function is_field_2d

  !============================================================================
  subroutine validate_field_name(fname, def_index)
    !--------------------------------------------------------------------------
    ! Validate that a field name is a known state variable, constituent,
    ! physics buffer field, or a previously defined derived field.
    !--------------------------------------------------------------------------
    use constituents,   only: cnst_get_ind
    use physics_buffer, only: pbuf_get_index

    character(len=*), intent(in) :: fname
    integer,          intent(in) :: def_index   ! index of current definition (for chaining check)

    integer :: k, idx
    logical :: found
    character(len=max_name_len) :: uname

    uname = adjustl(fname)
    found = .false.

    ! Check known 3D state variables
    do k = 1, n_known_3d_state
      if (trim(uname) == trim(known_3d_state_vars(k))) then
        found = .true.
        exit
      end if
    end do

    ! Check known 2D state variables
    if (.not. found) then
      do k = 1, n_known_2d_state
        if (trim(uname) == trim(known_2d_state_vars(k))) then
          found = .true.
          exit
        end if
      end do
    end if

    ! Check previously defined derived fields (for chaining)
    if (.not. found) then
      do k = 1, def_index - 1
        if (trim(uname) == trim(expressions(k)%output_name)) then
          found = .true.
          exit
        end if
      end do
    end if

    ! Check constituents
    if (.not. found) then
      call cnst_get_ind(trim(uname), idx, abrtf=.false.)
      found = (idx > 0)
    end if

    ! Check physics buffer
    if (.not. found) then
      idx = pbuf_get_index(trim(uname), errcode=k)
      found = (idx > 0)
    end if

    if (.not. found) then
      call endrun('eam_derived: validate_field_name: unknown field "'//trim(uname)// &
           '" in derived definition "'//trim(expressions(def_index)%output_name)// &
           '". Must be a state variable, constituent, physics buffer field, '// &
           'or a previously defined derived field.')
    end if

  end subroutine validate_field_name

  !============================================================================
  subroutine validate_tend_field(fname)
    !--------------------------------------------------------------------------
    ! Validate that a tendency field name is a known field or derived field.
    !--------------------------------------------------------------------------
    use constituents,   only: cnst_get_ind
    use physics_buffer, only: pbuf_get_index

    character(len=*), intent(in) :: fname

    integer :: k, idx
    logical :: found
    character(len=max_name_len) :: uname

    uname = adjustl(fname)
    found = .false.

    ! Check known 3D state variables
    do k = 1, n_known_3d_state
      if (trim(uname) == trim(known_3d_state_vars(k))) then
        found = .true.
        exit
      end if
    end do

    ! Check known 2D state variables
    if (.not. found) then
      do k = 1, n_known_2d_state
        if (trim(uname) == trim(known_2d_state_vars(k))) then
          found = .true.
          exit
        end if
      end do
    end if

    ! Check derived fields
    if (.not. found) then
      do k = 1, n_derived
        if (trim(uname) == trim(expressions(k)%output_name)) then
          found = .true.
          exit
        end if
      end do
    end if

    ! Check constituents
    if (.not. found) then
      call cnst_get_ind(trim(uname), idx, abrtf=.false.)
      found = (idx > 0)
    end if

    ! Check physics buffer
    if (.not. found) then
      idx = pbuf_get_index(trim(uname), errcode=k)
      found = (idx > 0)
    end if

    if (.not. found) then
      call endrun('eam_derived: validate_tend_field: unknown field "'//trim(uname)// &
           '". Must be a state variable, constituent, physics buffer field, or derived field.')
    end if

  end subroutine validate_tend_field

end module eam_derived
