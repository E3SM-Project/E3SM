module eam_vcoarsen
  !-------------------------------------------------------------------------------------------
  !
  ! EAM vertical coarsening wrapper.
  !
  ! Thin integration layer between EAM's physics state / history system and
  ! the shared shr_vcoarsen_mod routines. All vertical coarsening math is
  ! delegated to the shared module.
  !
  ! Supports three modes of vertical coarsening:
  !   1. Overlap-weighted averaging onto coarser pressure layers (pdel-weighted)
  !   2. Level selection by index (e.g., U_at_L5)
  !   3. Level selection by nearest pressure value (e.g., U_at_P850)
  !
  ! Configuration via namelist (eam_vcoarsen_nl):
  !   vcoarsen_pbounds         - pressure boundaries (Pa), top to surface
  !   vcoarsen_avg_flds        - fields to average onto coarsened layers
  !   vcoarsen_select_levs     - level indices for level selection (e.g., 1,5,10)
  !   vcoarsen_select_lev_flds - fields for level index selection
  !   vcoarsen_select_pres     - pressure values in hPa for nearest-level selection
  !   vcoarsen_select_pres_flds - fields for pressure value selection
  !
  ! Use 'all' as the sole entry in any field list to apply to all known 3D fields.
  !
  ! Usage from physpkg.F90:
  !   call eam_vcoarsen_readnl(nlfile)     ! during namelist reading
  !   call eam_vcoarsen_register()         ! during phys_init
  !   call eam_vcoarsen_write(state)       ! during physics timestep
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
  character(len=max_name_len) :: vcoarsen_avg_flds(max_flds)
  integer  :: vcoarsen_select_levs(max_select_vals)
  character(len=max_name_len) :: vcoarsen_select_lev_flds(max_flds)
  real(r8) :: vcoarsen_select_pres(max_select_vals)
  character(len=max_name_len) :: vcoarsen_select_pres_flds(max_flds)

  ! Parsed counts
  integer :: n_avg_levs      = 0   ! number of coarsened pressure layers
  integer :: n_avg_flds      = 0   ! number of fields for averaging
  integer :: n_sel_levs      = 0   ! number of level indices for selection
  integer :: n_sel_lev_flds  = 0   ! number of fields for level selection
  integer :: n_sel_pres      = 0   ! number of pressure values for selection
  integer :: n_sel_pres_flds = 0   ! number of fields for pressure selection

  ! Flags
  logical :: has_avg     = .false.
  logical :: has_sel_lev = .false.
  logical :: has_sel_pres = .false.
  logical :: module_is_initialized = .false.

contains

  !============================================================================
  subroutine eam_vcoarsen_readnl(nlfile)
    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: mpicom, masterprocid, mpi_real8, mpi_character, mpi_integer

    character(len=*), intent(in) :: nlfile

    integer :: unitn, ierr, i

    namelist /eam_vcoarsen_nl/ vcoarsen_pbounds, vcoarsen_avg_flds, &
         vcoarsen_select_levs, vcoarsen_select_lev_flds, &
         vcoarsen_select_pres, vcoarsen_select_pres_flds

    ! Initialize defaults
    vcoarsen_pbounds(:)         = -1.0_r8
    vcoarsen_avg_flds(:)        = ''
    vcoarsen_select_levs(:)     = -1
    vcoarsen_select_lev_flds(:) = ''
    vcoarsen_select_pres(:)     = -1.0_r8
    vcoarsen_select_pres_flds(:) = ''

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
    call mpi_bcast(vcoarsen_avg_flds, max_name_len*max_flds, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_select_levs, max_select_vals, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_select_lev_flds, max_name_len*max_flds, mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_select_pres, max_select_vals, mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast(vcoarsen_select_pres_flds, max_name_len*max_flds, mpi_character, masterprocid, mpicom, ierr)

    ! Parse pressure bounds for averaging
    n_avg_levs = 0
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

    if (masterproc) then
      if (has_avg) then
        write(iulog,*) 'eam_vcoarsen_readnl: averaging enabled, ', n_avg_levs, &
             ' layers, ', n_avg_flds, ' fields'
      end if
      if (has_sel_lev) then
        write(iulog,*) 'eam_vcoarsen_readnl: level selection enabled, ', n_sel_levs, &
             ' levels, ', n_sel_lev_flds, ' fields'
      end if
      if (has_sel_pres) then
        write(iulog,*) 'eam_vcoarsen_readnl: pressure selection enabled, ', n_sel_pres, &
             ' pressures, ', n_sel_pres_flds, ' fields'
      end if
    end if

  end subroutine eam_vcoarsen_readnl

  !============================================================================
  subroutine eam_vcoarsen_register()
    use cam_history, only: addfld, horiz_only
    use constituents, only: cnst_get_ind

    integer :: i, k, idx
    character(len=max_fname_len) :: fname
    character(len=256) :: lname
    logical :: found

    if (.not. has_avg .and. .not. has_sel_lev .and. .not. has_sel_pres) return

    ! Validate and register averaged fields
    if (has_avg) then
      do i = 1, n_avg_flds
        call validate_field_name(vcoarsen_avg_flds(i))
        do k = 1, n_avg_levs
          call make_avg_name(vcoarsen_avg_flds(i), k, fname)
          write(lname, '(A,A,I0,A,F0.1,A,F0.1,A)') &
               trim(vcoarsen_avg_flds(i)), ' vcoarsen layer ', k, &
               ' (', vcoarsen_pbounds(k)/100.0_r8, '-', &
               vcoarsen_pbounds(k+1)/100.0_r8, ' hPa)'
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

    module_is_initialized = .true.

    if (masterproc) then
      write(iulog,*) 'eam_vcoarsen_register: registration complete'
    end if

  end subroutine eam_vcoarsen_register

  !============================================================================
  subroutine eam_vcoarsen_write(state)
    use physics_types,  only: physics_state
    use cam_history,    only: outfld
    use shr_vcoarsen_mod, only: shr_vcoarsen_avg_cols, shr_vcoarsen_select_index, &
                                shr_vcoarsen_select_nearest

    type(physics_state), intent(in) :: state

    real(r8) :: src_field(pcols, pver)
    real(r8) :: coarsened(pcols, max(n_avg_levs, 1))
    real(r8) :: selected(pcols)
    real(r8) :: coord_iface(pcols, pverp)
    real(r8) :: coord_mid(pcols, pver)
    integer  :: nlev_max(pcols)
    integer  :: i, k, ncol, lchnk
    character(len=max_fname_len) :: fname

    if (.not. has_avg .and. .not. has_sel_lev .and. .not. has_sel_pres) return

    ncol  = state%ncol
    lchnk = state%lchnk

    ! All columns have pver valid levels in EAM
    nlev_max(1:ncol) = pver

    ! --- Overlap-weighted averaging ---
    if (has_avg) then
      ! Build interface pressure array
      coord_iface(1:ncol, 1:pverp) = state%pint(1:ncol, 1:pverp)

      do i = 1, n_avg_flds
        call get_state_field(state, vcoarsen_avg_flds(i), src_field, ncol)

        call shr_vcoarsen_avg_cols(src_field, coord_iface, ncol, pver, &
             vcoarsen_pbounds(1:n_avg_levs+1), n_avg_levs, fillvalue, coarsened)

        do k = 1, n_avg_levs
          call make_avg_name(vcoarsen_avg_flds(i), k, fname)
          call outfld(trim(fname), coarsened(:, k), pcols, lchnk)
        end do
      end do
    end if

    ! --- Level index selection ---
    if (has_sel_lev) then
      do i = 1, n_sel_lev_flds
        call get_state_field(state, vcoarsen_select_lev_flds(i), src_field, ncol)

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
        call get_state_field(state, vcoarsen_select_pres_flds(i), src_field, ncol)

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

  end subroutine eam_vcoarsen_write

  !============================================================================
  ! Private helper routines
  !============================================================================

  subroutine get_state_field(state, fname, field_out, ncol)
    use physics_types, only: physics_state
    use constituents,  only: cnst_get_ind

    type(physics_state), intent(in)  :: state
    character(len=*),    intent(in)  :: fname
    real(r8),            intent(out) :: field_out(pcols, pver)
    integer,             intent(in)  :: ncol

    integer :: idx
    character(len=max_name_len) :: uname

    uname = adjustl(fname)
    field_out(:,:) = 0.0_r8

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

    call endrun('eam_vcoarsen: get_state_field: unknown field: '//trim(uname))

  end subroutine get_state_field

  !============================================================================
  subroutine validate_field_name(fname)
    use constituents, only: cnst_get_ind

    character(len=*), intent(in) :: fname
    integer :: k, idx
    logical :: found
    character(len=max_name_len) :: uname

    uname = adjustl(fname)
    found = .false.
    do k = 1, n_known_state
      if (trim(uname) == trim(known_state_vars(k))) then
        found = .true.
        exit
      end if
    end do
    if (.not. found) then
      call cnst_get_ind(trim(uname), idx, abrtf=.false.)
      found = (idx > 0)
    end if
    if (.not. found) then
      call endrun('eam_vcoarsen: unknown field "'//trim(uname)// &
           '". Must be T, U, V, OMEGA, Z3, Q, or a registered constituent.')
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

end module eam_vcoarsen
