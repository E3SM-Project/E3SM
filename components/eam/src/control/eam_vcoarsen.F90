module eam_vcoarsen
  !-------------------------------------------------------------------------------------------
  !
  ! EAM vertical coarsening wrapper.
  !
  ! Thin integration layer between EAM's physics state / history system and
  ! the shared shr_vcoarsen_mod routines. All vertical coarsening math is
  ! delegated to the shared module.
  !
  ! Supports five modes of vertical coarsening:
  !   1. Overlap-weighted averaging onto coarser pressure layers (pdel-weighted)
  !   2. Level selection by index (e.g., U_at_L5)
  !   3. Level selection by nearest pressure value (e.g., U_at_P850)
  !   4. Column integration: sum(field * pdel / g) producing a 2D field (e.g., TOTAL_WATER_INT)
  !   5. Linear interpolation to height above surface (e.g., U_at_z10 for 10 m wind)
  !
  ! Configuration via namelist (eam_vcoarsen_nl):
  !   vcoarsen_pbounds         - pressure boundaries (Pa), top to surface
  !   vcoarsen_avg_flds        - fields to average onto coarsened layers
  !   vcoarsen_select_levs     - level indices for level selection (e.g., 1,5,10)
  !   vcoarsen_select_lev_flds - fields for level index selection
  !   vcoarsen_select_pres     - pressure values in hPa for nearest-level selection
  !   vcoarsen_select_pres_flds - fields for pressure value selection
  !   vcoarsen_select_heights   - heights (m) above surface for linear interpolation
  !   vcoarsen_select_height_flds - fields for height interpolation
  !   vcoarsen_int_flds        - fields to column-integrate (sum field*pdel/g)
  !
  ! Use 'all' as the sole entry in any field list to apply to all known 3D fields.
  !
  ! Field sources (checked in order):
  !   1. Physics state variables: T, U, V, OMEGA, Z3, Q
  !   2. Registered constituents (e.g., CLDLIQ, CLDICE)
  !   3. Physics buffer fields (e.g., CLD, CONCLD)
  !   4. Derived fields defined via eam_derived (e.g., TOTAL_WATER)
  !      -- requires eam_derived_write to run first for the same chunk
  !
  ! Usage from physpkg.F90:
  !   call eam_vcoarsen_readnl(nlfile)        ! during namelist reading
  !   call eam_vcoarsen_register()            ! during phys_register
  !   call eam_vcoarsen_write(state, pbuf)    ! after eam_derived_write in tphysac loop
  !
  !-------------------------------------------------------------------------------------------

  use shr_kind_mod,          only: r8 => shr_kind_r8
  use ppgrid,                only: pcols, pver, pverp
  use cam_logfile,           only: iulog
  use cam_abortutils,        only: endrun
  use spmd_utils,            only: masterproc
  use cam_history_support,   only: fillvalue

  implicit none
  private
  save

  public :: eam_vcoarsen_readnl
  public :: eam_vcoarsen_register
  public :: eam_vcoarsen_write

  ! Parameters
  integer, parameter :: max_pbounds       = 51    ! max pressure boundaries (=> max 50 layers)
  integer, parameter :: max_flds          = 100   ! max fields per mode
  integer, parameter :: max_select_vals   = 50    ! max level indices or pressure values
  integer, parameter :: max_name_len      = 34    ! matches fieldname_len in cam_history
  integer, parameter :: max_fname_len     = 48    ! output field names can be longer

  ! Known 3D state variable names
  integer, parameter :: n_known_state = 6
  character(len=max_name_len), parameter :: known_state_vars(n_known_state) = &
       (/ 'T     ', 'U     ', 'V     ', 'OMEGA ', 'Z3    ', 'Q     ' /)

  ! Namelist variables
  real(r8) :: vcoarsen_pbounds(max_pbounds)
  integer  :: vcoarsen_level_bounds(max_pbounds)
  character(len=max_name_len) :: vcoarsen_avg_flds(max_flds)
  integer  :: vcoarsen_select_levs(max_select_vals)
  character(len=max_name_len) :: vcoarsen_select_lev_flds(max_flds)
  real(r8) :: vcoarsen_select_pres(max_select_vals)
  character(len=max_name_len) :: vcoarsen_select_pres_flds(max_flds)
  real(r8) :: vcoarsen_select_heights(max_select_vals)
  character(len=max_name_len) :: vcoarsen_select_height_flds(max_flds)
  character(len=max_name_len) :: vcoarsen_int_flds(max_flds)

  ! Parsed counts
  integer :: n_avg_levs      = 0   ! number of coarsened output layers
  integer :: n_avg_flds      = 0   ! number of fields for averaging
  logical :: use_level_bounds = .false.  ! true = index-based, false = pressure-based
  integer :: n_sel_levs      = 0   ! number of level indices for selection
  integer :: n_sel_lev_flds  = 0   ! number of fields for level selection
  integer :: n_sel_pres      = 0   ! number of pressure values for selection
  integer :: n_sel_pres_flds = 0   ! number of fields for pressure selection
  integer :: n_sel_heights     = 0 ! number of heights for selection
  integer :: n_sel_height_flds = 0 ! number of fields for height selection
  integer :: n_int_flds      = 0   ! number of fields for column integration

  ! Flags
  logical :: has_avg        = .false.
  logical :: has_sel_lev    = .false.
  logical :: has_sel_pres   = .false.
  logical :: has_sel_height = .false.
  logical :: has_int        = .false.

  ! Column integration tendency: previous timestep values (pcols, max_flds, begchunk:endchunk)
  real(r8), allocatable :: int_prev(:,:,:)
  logical :: int_prev_valid = .false.

contains

  !============================================================================
  subroutine eam_vcoarsen_readnl(nlfile)
    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: mpicom, masterprocid, mpi_real8, mpi_character, mpi_integer

    character(len=*), intent(in) :: nlfile

    integer :: unitn, ierr, i

    namelist /eam_vcoarsen_nl/ vcoarsen_pbounds, vcoarsen_level_bounds, &
         vcoarsen_avg_flds, &
         vcoarsen_select_levs, vcoarsen_select_lev_flds, &
         vcoarsen_select_pres, vcoarsen_select_pres_flds, &
         vcoarsen_select_heights, vcoarsen_select_height_flds, &
         vcoarsen_int_flds

    ! Initialize defaults
    vcoarsen_pbounds(:)         = -1.0_r8
    vcoarsen_level_bounds(:)    = -1
    vcoarsen_avg_flds(:)        = ''
    vcoarsen_select_levs(:)     = -1
    vcoarsen_select_lev_flds(:) = ''
    vcoarsen_select_pres(:)     = -1.0_r8
    vcoarsen_select_pres_flds(:) = ''
    vcoarsen_select_heights(:)     = -1.0_r8
    vcoarsen_select_height_flds(:) = ''
    vcoarsen_int_flds(:)        = ''

    if (masterproc) then
      unitn = getunit()
      open(unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'eam_vcoarsen_nl', status=ierr)
      if (ierr == 0) then
        read(unitn, eam_vcoarsen_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun('eam_vcoarsen_readnl: ERROR reading namelist eam_vcoarsen_nl')
        end if
      end if
      close(unitn)
      call freeunit(unitn)
    end if

    ! Broadcast
    call mpi_bcast(vcoarsen_pbounds, max_pbounds, mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_level_bounds, max_pbounds, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_avg_flds, max_name_len*max_flds, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_select_levs, max_select_vals, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_select_lev_flds, max_name_len*max_flds, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_select_pres, max_select_vals, mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_select_pres_flds, max_name_len*max_flds, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_select_heights, max_select_vals, mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_select_height_flds, max_name_len*max_flds, mpi_character, masterprocid, mpicom, ierr)

    ! Determine averaging mode: level-index bounds take precedence over pressure bounds
    ! Parse level-index bounds first
    n_avg_levs = 0
    use_level_bounds = .false.
    do i = 2, max_pbounds
      if (vcoarsen_level_bounds(i) < 0) exit
      n_avg_levs = n_avg_levs + 1
    end do
    if (n_avg_levs > 0) then
      use_level_bounds = .true.
      has_avg = .true.
      ! Validate: bounds must be non-negative and strictly increasing
      if (vcoarsen_level_bounds(1) < 0) then
        call endrun('eam_vcoarsen_readnl: vcoarsen_level_bounds(1) must be non-negative')
      end if
      do i = 1, n_avg_levs
        if (vcoarsen_level_bounds(i+1) <= vcoarsen_level_bounds(i)) then
          call endrun('eam_vcoarsen_readnl: vcoarsen_level_bounds must be strictly increasing')
        end if
      end do
    else
      ! Fall back to pressure bounds
      do i = 2, max_pbounds
        if (vcoarsen_pbounds(i) < 0.0_r8) exit
        n_avg_levs = n_avg_levs + 1
      end do
      has_avg = (n_avg_levs > 0)

      ! Validate pressure bounds
      if (has_avg) then
        if (vcoarsen_pbounds(1) < 0.0_r8) then
          call endrun('eam_vcoarsen_readnl: vcoarsen_pbounds(1) must be non-negative')
        end if
        do i = 1, n_avg_levs
          if (vcoarsen_pbounds(i+1) <= vcoarsen_pbounds(i)) then
            call endrun('eam_vcoarsen_readnl: vcoarsen_pbounds must be strictly increasing')
          end if
        end do
      end if
    end if

    ! Count averaging fields; expand 'all'
    call count_and_expand_flds(vcoarsen_avg_flds, n_avg_flds)

    ! Count level selection indices
    n_sel_levs = 0
    do i = 1, max_select_vals
      if (vcoarsen_select_levs(i) < 0) exit
      n_sel_levs = n_sel_levs + 1
    end do
    has_sel_lev = (n_sel_levs > 0)

    ! Count level selection fields; expand 'all'
    call count_and_expand_flds(vcoarsen_select_lev_flds, n_sel_lev_flds)

    ! Count pressure selection values
    n_sel_pres = 0
    do i = 1, max_select_vals
      if (vcoarsen_select_pres(i) < 0.0_r8) exit
      n_sel_pres = n_sel_pres + 1
    end do
    has_sel_pres = (n_sel_pres > 0)

    ! Count pressure selection fields; expand 'all'
    call count_and_expand_flds(vcoarsen_select_pres_flds, n_sel_pres_flds)

    ! Count height selection values (heights >= 0 m above surface)
    n_sel_heights = 0
    do i = 1, max_select_vals
      if (vcoarsen_select_heights(i) < 0.0_r8) exit
      n_sel_heights = n_sel_heights + 1
    end do
    has_sel_height = (n_sel_heights > 0)

    ! Count height selection fields; expand 'all'
    call count_and_expand_flds(vcoarsen_select_height_flds, n_sel_height_flds)

    ! Column integration fields
    call mpi_bcast(vcoarsen_int_flds, max_name_len*max_flds, mpi_character, masterprocid, mpicom, ierr)
    call count_and_expand_flds(vcoarsen_int_flds, n_int_flds)
    has_int = (n_int_flds > 0)

    if (masterproc) then
      if (has_avg) then
        if (use_level_bounds) then
          write(iulog,*) 'eam_vcoarsen_readnl: index-based dp-weighted averaging, ', &
               n_avg_levs, ' layers, ', n_avg_flds, ' fields'
        else
          write(iulog,*) 'eam_vcoarsen_readnl: pressure-bounds averaging, ', &
               n_avg_levs, ' layers, ', n_avg_flds, ' fields'
        end if
      end if
      if (has_sel_lev) then
        write(iulog,*) 'eam_vcoarsen_readnl: level selection enabled, ', n_sel_levs, &
             ' levels, ', n_sel_lev_flds, ' fields'
      end if
      if (has_sel_pres) then
        write(iulog,*) 'eam_vcoarsen_readnl: pressure selection enabled, ', n_sel_pres, &
             ' pressures, ', n_sel_pres_flds, ' fields'
      end if
      if (has_sel_height) then
        write(iulog,*) 'eam_vcoarsen_readnl: height selection enabled, ', n_sel_heights, &
             ' heights, ', n_sel_height_flds, ' fields'
      end if
      if (has_int) then
        write(iulog,*) 'eam_vcoarsen_readnl: column integration enabled, ', n_int_flds, ' fields'
      end if
    end if

  end subroutine eam_vcoarsen_readnl

  !============================================================================
  subroutine eam_vcoarsen_register()
    use cam_history, only: addfld, horiz_only
    use constituents, only: cnst_get_ind
    use ppgrid,       only: begchunk, endchunk

    integer :: i, k, idx
    character(len=max_fname_len) :: fname
    character(len=256) :: lname
    logical :: found

    if (.not. has_avg .and. .not. has_sel_lev .and. .not. has_sel_pres &
        .and. .not. has_sel_height .and. .not. has_int) return

    ! Validate and register averaged fields
    if (has_avg) then
      do i = 1, n_avg_flds
        call validate_field_name(vcoarsen_avg_flds(i))
        do k = 1, n_avg_levs
          call make_avg_name(vcoarsen_avg_flds(i), k-1, fname)  ! 0-indexed for ACE
          if (use_level_bounds) then
            write(lname, '(A,A,I0,A,I0,A,I0,A)') &
                 trim(vcoarsen_avg_flds(i)), ' vcoarsen layer ', k-1, &
                 ' (levels ', vcoarsen_level_bounds(k), '-', &
                 vcoarsen_level_bounds(k+1), ')'
          else
            write(lname, '(A,A,I0,A,F0.1,A,F0.1,A)') &
                 trim(vcoarsen_avg_flds(i)), ' vcoarsen layer ', k-1, &
                 ' (', vcoarsen_pbounds(k)/100.0_r8, '-', &
                 vcoarsen_pbounds(k+1)/100.0_r8, ' hPa)'
          end if
          call addfld(trim(fname), horiz_only, 'A', 'varies', trim(lname), &
               flag_xyfill=.true.)
        end do
      end do
    end if

    ! Validate and register level-selected fields
    if (has_sel_lev) then
      do i = 1, n_sel_lev_flds
        call validate_field_name(vcoarsen_select_lev_flds(i))
        do k = 1, n_sel_levs
          call make_sel_lev_name(vcoarsen_select_lev_flds(i), vcoarsen_select_levs(k), fname)
          write(lname, '(A,A,I0)') trim(vcoarsen_select_lev_flds(i)), &
               ' at model level ', vcoarsen_select_levs(k)
          call addfld(trim(fname), horiz_only, 'A', 'varies', trim(lname), &
               flag_xyfill=.true.)
        end do
      end do
    end if

    ! Validate and register pressure-selected fields
    if (has_sel_pres) then
      do i = 1, n_sel_pres_flds
        call validate_field_name(vcoarsen_select_pres_flds(i))
        do k = 1, n_sel_pres
          call make_sel_pres_name(vcoarsen_select_pres_flds(i), vcoarsen_select_pres(k), fname)
          write(lname, '(A,A,F0.1,A)') trim(vcoarsen_select_pres_flds(i)), &
               ' at nearest level to ', vcoarsen_select_pres(k), ' hPa'
          call addfld(trim(fname), horiz_only, 'A', 'varies', trim(lname), &
               flag_xyfill=.true.)
        end do
      end do
    end if

    ! Validate and register height-interpolated fields
    if (has_sel_height) then
      do i = 1, n_sel_height_flds
        call validate_field_name(vcoarsen_select_height_flds(i))
        do k = 1, n_sel_heights
          call make_sel_height_name(vcoarsen_select_height_flds(i), &
               vcoarsen_select_heights(k), fname)
          write(lname, '(A,A,F0.1,A)') trim(vcoarsen_select_height_flds(i)), &
               ' linearly interpolated to ', vcoarsen_select_heights(k), &
               ' m above surface'
          call addfld(trim(fname), horiz_only, 'A', 'varies', trim(lname), &
               flag_xyfill=.true.)
        end do
      end do
    end if

    ! Validate and register column-integrated fields + their tendencies
    if (has_int) then
      ! Allocate previous-timestep storage for tendency computation
      allocate(int_prev(pcols, n_int_flds, begchunk:endchunk))
      int_prev(:,:,:) = 0.0_r8
      do i = 1, n_int_flds
        call validate_field_name(vcoarsen_int_flds(i))
        fname = trim(vcoarsen_int_flds(i)) // '_INT'
        write(lname, '(A,A)') 'Column-integrated ', trim(vcoarsen_int_flds(i))
        call addfld(trim(fname), horiz_only, 'A', 'varies', trim(lname), &
             flag_xyfill=.true.)
        ! Also register the tendency d{FIELD}_INT_dt
        fname = 'd' // trim(vcoarsen_int_flds(i)) // '_INT_dt'
        write(lname, '(A,A)') 'Tendency of column-integrated ', trim(vcoarsen_int_flds(i))
        call addfld(trim(fname), horiz_only, 'A', 'varies', trim(lname), &
             flag_xyfill=.true.)
      end do
    end if

    if (masterproc) then
      write(iulog,*) 'eam_vcoarsen_register: registration complete'
    end if

  end subroutine eam_vcoarsen_register

  !============================================================================
  subroutine eam_vcoarsen_write(state, pbuf_chunk)
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    use cam_history,    only: outfld
    use shr_vcoarsen_mod, only: shr_vcoarsen_avg_cols, shr_vcoarsen_select_index, &
                                shr_vcoarsen_select_nearest, &
                                shr_vcoarsen_avg_by_index_cols
    use physconst,      only: gravit
    use time_manager,   only: get_step_size

    type(physics_state), intent(in) :: state
    type(physics_buffer_desc), pointer :: pbuf_chunk(:)

    real(r8) :: src_field(pcols, pver)
    real(r8) :: coarsened(pcols, max(n_avg_levs, 1))
    real(r8) :: selected(pcols)
    real(r8) :: integrated(pcols)
    real(r8) :: coord_iface(pcols, pverp)
    real(r8) :: coord_mid(pcols, pver)
    real(r8) :: dt
    integer  :: nlev_max(pcols)
    integer  :: i, k, ncol, lchnk
    character(len=max_fname_len) :: fname

    if (.not. has_avg .and. .not. has_sel_lev .and. .not. has_sel_pres &
        .and. .not. has_sel_height .and. .not. has_int) return

    ncol  = state%ncol
    lchnk = state%lchnk

    ! All columns have pver valid levels in EAM
    nlev_max(1:ncol) = pver

    ! --- dp-weighted vertical averaging ---
    if (has_avg) then
      do i = 1, n_avg_flds
        call get_state_field(state, pbuf_chunk, vcoarsen_avg_flds(i), src_field, ncol)

        if (use_level_bounds) then
          ! Index-based: dp-weighted average over level index ranges.
          ! dp(k) = pint(k+1) - pint(k), the pressure thickness of each level.
          do k = 1, pver
            coord_mid(1:ncol, k) = state%pint(1:ncol, k+1) - state%pint(1:ncol, k)
          end do
          call shr_vcoarsen_avg_by_index_cols(src_field(1:ncol, :), &
               coord_mid(1:ncol, :), ncol, pver, &
               vcoarsen_level_bounds(1:n_avg_levs+1), n_avg_levs, &
               fillvalue, coarsened(1:ncol, :))
        else
          ! Pressure-bounds: overlap-weighted averaging between pressure interfaces.
          coord_iface(1:ncol, 1:pverp) = state%pint(1:ncol, 1:pverp)
          call shr_vcoarsen_avg_cols(src_field(1:ncol, :), coord_iface(1:ncol, :), &
               ncol, pver, vcoarsen_pbounds(1:n_avg_levs+1), n_avg_levs, &
               fillvalue, coarsened(1:ncol, :))
        end if

        do k = 1, n_avg_levs
          ! Fill inactive columns
          coarsened(ncol+1:pcols, k) = fillvalue

          call make_avg_name(vcoarsen_avg_flds(i), k-1, fname)  ! 0-indexed
          call outfld(trim(fname), coarsened(:, k), pcols, lchnk)
        end do
      end do
    end if

    ! --- Level index selection ---
    if (has_sel_lev) then
      do i = 1, n_sel_lev_flds
        call get_state_field(state, pbuf_chunk, vcoarsen_select_lev_flds(i), src_field, ncol)

        do k = 1, n_sel_levs
          call shr_vcoarsen_select_index(src_field, ncol, pver, &
               vcoarsen_select_levs(k), nlev_max, fillvalue, selected)

          ! Fill inactive columns
          selected(ncol+1:pcols) = fillvalue

          call make_sel_lev_name(vcoarsen_select_lev_flds(i), vcoarsen_select_levs(k), fname)
          call outfld(trim(fname), selected, pcols, lchnk)
        end do
      end do
    end if

    ! --- Nearest pressure selection ---
    if (has_sel_pres) then
      ! Build midpoint pressure array
      coord_mid(1:ncol, 1:pver) = state%pmid(1:ncol, 1:pver)

      do i = 1, n_sel_pres_flds
        call get_state_field(state, pbuf_chunk, vcoarsen_select_pres_flds(i), src_field, ncol)

        do k = 1, n_sel_pres
          ! Convert hPa to Pa for comparison with pmid
          call shr_vcoarsen_select_nearest(src_field, coord_mid, ncol, pver, &
               vcoarsen_select_pres(k) * 100.0_r8, nlev_max, fillvalue, selected)

          ! Fill inactive columns
          selected(ncol+1:pcols) = fillvalue

          call make_sel_pres_name(vcoarsen_select_pres_flds(i), vcoarsen_select_pres(k), fname)
          call outfld(trim(fname), selected, pcols, lchnk)
        end do
      end do
    end if

    ! --- Height-above-surface linear interpolation ---
    ! state%zm is geopotential height above surface at midpoints (m); it is
    ! monotonically decreasing with level index (level 1 is the model top).
    ! shr_vcoarsen_select_nearest handles either coordinate orientation.
    if (has_sel_height) then
      coord_mid(1:ncol, 1:pver) = state%zm(1:ncol, 1:pver)

      do i = 1, n_sel_height_flds
        call get_state_field(state, pbuf_chunk, vcoarsen_select_height_flds(i), &
             src_field, ncol)

        do k = 1, n_sel_heights
          call shr_vcoarsen_select_nearest(src_field, coord_mid, ncol, pver, &
               vcoarsen_select_heights(k), nlev_max, fillvalue, selected)

          ! Fill inactive columns
          selected(ncol+1:pcols) = fillvalue

          call make_sel_height_name(vcoarsen_select_height_flds(i), &
               vcoarsen_select_heights(k), fname)
          call outfld(trim(fname), selected, pcols, lchnk)
        end do
      end do
    end if

    ! --- Column integration: sum(field * pdel / g) ---
    if (has_int) then
      ! int_prev allocated in eam_vcoarsen_register

      do i = 1, n_int_flds
        call get_state_field(state, pbuf_chunk, vcoarsen_int_flds(i), src_field, ncol)

        integrated(:) = fillvalue
        integrated(1:ncol) = 0.0_r8
        do k = 1, pver
          integrated(1:ncol) = integrated(1:ncol) + &
               src_field(1:ncol, k) * state%pdel(1:ncol, k) / gravit
        end do

        fname = trim(vcoarsen_int_flds(i)) // '_INT'
        call outfld(trim(fname), integrated, pcols, lchnk)

        ! Compute tendency: (current - previous) / dt
        dt = real(get_step_size(), r8)
        selected(:) = fillvalue
        if (int_prev_valid) then
          selected(1:ncol) = (integrated(1:ncol) - int_prev(1:ncol, i, lchnk)) / dt
        else
          selected(1:ncol) = 0.0_r8
        end if
        fname = 'd' // trim(vcoarsen_int_flds(i)) // '_INT_dt'
        call outfld(trim(fname), selected, pcols, lchnk)

        ! Save current for next timestep
        int_prev(1:ncol, i, lchnk) = integrated(1:ncol)
      end do

      ! Mark previous values as valid after first complete pass through all chunks.
      ! This is approximate — it becomes valid after the first timestep completes.
      if (.not. int_prev_valid) int_prev_valid = .true.
    end if

  end subroutine eam_vcoarsen_write

  !============================================================================
  ! Private helper routines
  !============================================================================

  subroutine get_state_field(state, pbuf_chunk, fname, field_out, ncol)
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
    use constituents,   only: cnst_get_ind
    use eam_derived,    only: eam_derived_get_cache

    type(physics_state), intent(in)  :: state
    type(physics_buffer_desc), pointer :: pbuf_chunk(:)
    character(len=*),    intent(in)  :: fname
    real(r8),            intent(out) :: field_out(pcols, pver)
    integer,             intent(in)  :: ncol

    integer :: idx, pbuf_idx, errcode
    character(len=max_name_len) :: uname
    real(r8), pointer :: pbuf_fld(:,:)
    logical :: found

    uname = adjustl(fname)
    field_out(:,:) = 0.0_r8

    ! Check standard state variables first
    select case (trim(uname))
    case ('T')
      field_out(1:ncol, :) = state%t(1:ncol, :)
      return
    case ('U')
      field_out(1:ncol, :) = state%u(1:ncol, :)
      return
    case ('V')
      field_out(1:ncol, :) = state%v(1:ncol, :)
      return
    case ('OMEGA')
      field_out(1:ncol, :) = state%omega(1:ncol, :)
      return
    case ('Z3')
      field_out(1:ncol, :) = state%zm(1:ncol, :)
      return
    case ('Q')
      field_out(1:ncol, :) = state%q(1:ncol, :, 1)
      return
    end select

    ! Try constituent lookup
    call cnst_get_ind(trim(uname), idx, abrtf=.false.)
    if (idx > 0) then
      field_out(1:ncol, :) = state%q(1:ncol, :, idx)
      return
    end if

    ! Try physics buffer lookup
    ! pbuf_get_index returns index>0 if found, -1 if not (errcode=index)
    pbuf_idx = pbuf_get_index(trim(uname), errcode)
    if (pbuf_idx > 0) then
      call pbuf_get_field(pbuf_chunk, pbuf_idx, pbuf_fld)
      field_out(1:ncol, :) = pbuf_fld(1:ncol, :)
      return
    end if

    ! Try derived field cache (populated by eam_derived_write before vcoarsen runs)
    call eam_derived_get_cache(trim(uname), state%lchnk, field_out, ncol, pver, found)
    if (found) return

    call endrun('eam_vcoarsen: get_state_field: unknown field "'//trim(uname)// &
         '". Must be a state variable, constituent, physics buffer field, or derived field.')

  end subroutine get_state_field

  !============================================================================
  subroutine validate_field_name(fname)
    use constituents,   only: cnst_get_ind
    use physics_buffer, only: pbuf_get_index
    use eam_derived,    only: eam_derived_is_defined

    character(len=*), intent(in) :: fname
    integer :: k, idx, errcode
    logical :: found
    character(len=max_name_len) :: uname

    uname = adjustl(fname)
    found = .false.

    ! Check known state variables
    do k = 1, n_known_state
      if (trim(uname) == trim(known_state_vars(k))) then
        found = .true.
        exit
      end if
    end do

    ! Check constituents
    if (.not. found) then
      call cnst_get_ind(trim(uname), idx, abrtf=.false.)
      found = (idx > 0)
    end if

    ! Check physics buffer
    if (.not. found) then
      idx = pbuf_get_index(trim(uname), errcode)
      found = (idx > 0)
    end if

    ! Check derived fields (populated by eam_derived_readnl before phys_register)
    if (.not. found) found = eam_derived_is_defined(trim(uname))

    if (.not. found) then
      call endrun('eam_vcoarsen: unknown field "'//trim(uname)// &
           '". Must be a state variable, constituent, physics buffer field, or derived field.')
    end if

  end subroutine validate_field_name

  !============================================================================
  subroutine count_and_expand_flds(fld_list, n_flds)
    !--------------------------------------------------------------------------
    ! Count non-empty entries in fld_list. If the first entry is 'all',
    ! expand to all known state variables.
    !--------------------------------------------------------------------------
    character(len=max_name_len), intent(inout) :: fld_list(max_flds)
    integer, intent(out) :: n_flds

    integer :: i

    n_flds = 0

    ! Check for 'all' expansion
    if (trim(adjustl(fld_list(1))) == 'all') then
      do i = 1, n_known_state
        fld_list(i) = known_state_vars(i)
      end do
      fld_list(n_known_state + 1:) = ''
      n_flds = n_known_state
      return
    end if

    ! Count non-empty entries
    do i = 1, max_flds
      if (len_trim(fld_list(i)) == 0) exit
      n_flds = n_flds + 1
    end do

  end subroutine count_and_expand_flds

  !============================================================================
  subroutine make_avg_name(base_name, layer_idx, out_name)
    ! E.g., base_name="U", layer_idx=3 => out_name="U_3"
    character(len=*), intent(in)  :: base_name
    integer,          intent(in)  :: layer_idx
    character(len=*), intent(out) :: out_name

    character(len=4) :: idx_str

    write(idx_str, '(I0)') layer_idx
    out_name = trim(base_name) // '_' // trim(idx_str)

  end subroutine make_avg_name

  !============================================================================
  subroutine make_sel_lev_name(base_name, lev_idx, out_name)
    ! E.g., base_name="U", lev_idx=5 => out_name="U_at_L5"
    character(len=*), intent(in)  :: base_name
    integer,          intent(in)  :: lev_idx
    character(len=*), intent(out) :: out_name

    character(len=4) :: idx_str

    write(idx_str, '(I0)') lev_idx
    out_name = trim(base_name) // '_at_L' // trim(idx_str)

  end subroutine make_sel_lev_name

  !============================================================================
  subroutine make_sel_pres_name(base_name, pres_hpa, out_name)
    ! E.g., base_name="U", pres_hpa=850.0 => out_name="U_at_P850"
    ! For integer-like values, omit decimal. For fractional, include one decimal.
    character(len=*), intent(in)  :: base_name
    real(r8),         intent(in)  :: pres_hpa
    character(len=*), intent(out) :: out_name

    character(len=16) :: pstr
    integer :: int_pres

    int_pres = nint(pres_hpa)
    if (abs(pres_hpa - real(int_pres, r8)) < 0.01_r8) then
      write(pstr, '(I0)') int_pres
    else
      write(pstr, '(F0.1)') pres_hpa
    end if
    out_name = trim(base_name) // '_at_P' // trim(pstr)

  end subroutine make_sel_pres_name

  !============================================================================
  subroutine make_sel_height_name(base_name, height_m, out_name)
    ! E.g., base_name="U", height_m=10.0 => out_name="U_at_z10"
    ! For integer-like values, omit decimal. For fractional, include one decimal.
    character(len=*), intent(in)  :: base_name
    real(r8),         intent(in)  :: height_m
    character(len=*), intent(out) :: out_name

    character(len=16) :: hstr
    integer :: int_h

    int_h = nint(height_m)
    if (abs(height_m - real(int_h, r8)) < 0.01_r8) then
      write(hstr, '(I0)') int_h
    else
      write(hstr, '(F0.1)') height_m
    end if
    out_name = trim(base_name) // '_at_z' // trim(hstr)

  end subroutine make_sel_height_name

end module eam_vcoarsen
