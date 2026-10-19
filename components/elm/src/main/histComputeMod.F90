module histComputeMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute side of the ELM history: accumulates each time step's model
  ! fields into the history buffers, and normalizes / zeroes the buffers at the
  ! end of a history interval. histFileMod is the API: it owns the master field
  ! list, the tapes, the namelist, I/O and restart, and describes the active
  ! fields to this module (hist_compute_build) as a table of
  ! hist_active_field.
  !
  ! Layout
  ! - The buffers of all active fields on all tapes are one contiguous pool,
  !   ordered by tape, then field; each field is stored [output point, level].
  !   The tapes' hbuf/nacs are pointer views into the pool (hist_compute_views).
  ! - Model arrays are reached through histFileMod's pointer lists
  !   (elmptr_rs/elmptr_ra, passed in as src1/src2). The first element of a
  !   model array is the first point of its subgrid level, whatever its lower
  !   bound (some are sections with the lower bound reset to 1).
  !
  ! Per time step (hist_compute_update)
  ! - Multi-layer snow fields are first copied into snow_values with absent
  !   layers set to the field's no-snow value (spval or 0).
  ! - One kernel accumulates every field: one gang per field, taking the fields
  !   in an order grouped by subgrid level (process_order). For each level of
  !   the field, the gang decides once what to do and calls vector routines
  !   that loop over the points: a field averaged to gridcells first gets its
  !   averages (the x2g_range routines of subgridAveViewMod, which subgridAveMod
  !   uses too, so the sums and results are the same) into gridcell_samples;
  !   then accumulate_samples applies the time operation, chosen once, to all
  !   points.
  ! - OpenACC directives only; a build without OpenACC runs serially.
  !
  ! Time operations (avgflag in histFileMod)
  !   'I' time_inst       the sample, spval if invalid
  !   'A' time_mean       mean of the valid samples
  !   'X' time_max        maximum, with the rules of the original ELM code:
  !   'M' time_min        minimum   - averaged to gridcells: a missing sample
  !                                   sets spval at any time in the interval;
  !                                 - native: a missing or inactive sample sets
  !                                   spval only if it is the first one;
  !                                 - nacs = 1 after every sample; once spval,
  !                                   max stays spval while min recovers.
  !   'U' time_max_valid  maximum of the valid samples (spval if none)
  !   'L' time_min_valid  minimum of the valid samples (spval if none)
  ! A sample is valid where its point is active (averaged fields: always) and
  ! its value isn't spval. For 'A', 'U' and 'L', nacs counts the valid samples.
  !
  ! The kernels need the model arrays and the subgrid view on the device
  ! before hist_compute_build (OpenACC builds); see subgrid_ave_view_attach.
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use subgridAveViewMod, only : subgrid_ave_view, spval
  use subgridAveViewMod, only : p2g_range, c2g_range, l2g_range, t2g_range
  use subgridAveViewMod, only : p2g_range_unity, c2g_range_unity, scale_unity
  !
  implicit none
  private
  !
  ! !PUBLIC DATA:
  ! subgrid level of a model array
  integer, parameter, public :: lev_grid = 1, lev_topo = 2, lev_lun = 3, lev_col = 4, lev_patch = 5
  ! time operations
  integer, parameter, public :: time_inst = 1, time_mean = 2, time_max = 3, time_min = 4
  integer, parameter, public :: time_max_valid = 5, time_min_valid = 6
  ! multi-layer snow fields: value used where a snow layer is absent
  integer, parameter, public :: no_snow_unset  = 0 ! not a snow field
  integer, parameter, public :: no_snow_normal = 1 ! spval (the layer doesn't count)
  integer, parameter, public :: no_snow_zero   = 2 ! 0

  ! model array registered for a history field (histFileMod: elmptr_rs/elmptr_ra).
  ! Contiguous, so they are passed to the averages without copies.
  type, public :: hist_src1d
     real(r8), pointer, contiguous :: ptr(:) => null()
  end type hist_src1d
  type, public :: hist_src2d
     real(r8), pointer, contiguous :: ptr(:,:) => null()
  end type hist_src2d

  ! one active field of a tape
  type, public :: hist_active_field
     ! set by histFileMod
     integer :: source_index = 0             ! index into src1 (rank 1) or src2 (rank 2)
     integer :: source_rank  = 1
     integer :: source_level = lev_grid      ! subgrid level of the model array
     integer :: num_levels   = 1             ! levels of the field (num2d)
     integer :: time_op      = time_mean
     logical :: average_to_gridcell = .false.
     integer :: c2l_scale = 1, l2g_scale = 1 ! scale types (subgridAveMod's numbering)
     integer :: no_snow_behavior = no_snow_unset
     ! set by hist_compute_build
     integer :: source_first = 1, source_last = 0 ! points of the source level
     integer :: first_point = 1, last_point = 0   ! output points (gridcells if averaged)
     integer :: buffer_offset = 0                 ! start of the field in hbuf_pool/nacs_pool
     integer :: snow_offset = -1                  ! start in snow_values (snow fields)
     integer :: sample_offset = -1                ! start in gridcell_samples (averaged fields)
  end type hist_active_field
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: hist_compute_build         ! set up buffers and the field table
  public :: hist_compute_views         ! pointer views of one field's buffers
  public :: hist_compute_update        ! accumulate the current time step
  public :: hist_compute_normalize     ! divide time averages by their sample counts (one tape)
  public :: hist_compute_zero          ! reset buffers and counters (one tape)
  public :: hist_compute_update_host   ! copy one tape's buffers device -> host (OpenACC)
  public :: hist_compute_update_device ! copy all buffers host -> device (OpenACC)
  !
  ! !PRIVATE DATA:
  real(r8), parameter :: big = 1.e50_r8   ! start of a max (-big) or min (+big)

  integer :: num_fields = 0
  integer :: max_levels = 0
  type(hist_active_field), allocatable :: fields(:)
  integer , allocatable :: tape_first_field(:)  ! fields of tape t: tape_first_field(t):tape_first_field(t+1)-1
  integer , allocatable :: process_order(:)     ! the fields grouped by (subgrid level, averaged)
  integer , allocatable :: snow_fields(:)       ! the fields that are snow fields
  real(r8), allocatable, target :: hbuf_pool(:)
  integer , allocatable, target :: nacs_pool(:)
  real(r8), allocatable :: snow_values(:)       ! snow fields with absent layers filled in
  real(r8), allocatable :: gridcell_samples(:)  ! averages of the averaged fields at this step
  logical , allocatable :: all_active(:)        ! .true.: points without an active check
  logical :: extrasnowlayers = .false.          ! snow levels are ordered top layer first
  !$acc declare create(fields, hbuf_pool, nacs_pool, snow_values, gridcell_samples, all_active, &
  !$acc                extrasnowlayers)

  ! model arrays attached on the device by hist_compute_build
  integer, allocatable :: attached1(:), attached2(:)
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine hist_compute_build(view, field_list, tape_first, src1, src2, extrasnow, ierr, errmsg)
    !
    ! !DESCRIPTION:
    ! Check the active fields, lay out their buffers (zeroed) and keep the
    ! field table. field_list(tape_first(t):tape_first(t+1)-1) are the fields
    ! of tape t. Existing buffers are discarded.
    !
    ! !ARGUMENTS:
    type(subgrid_ave_view) , intent(in)  :: view
    type(hist_active_field), intent(in)  :: field_list(:)
    integer                , intent(in)  :: tape_first(:)
    type(hist_src1d)       , intent(in)  :: src1(:)
    type(hist_src2d)       , intent(in)  :: src2(:)
    logical                , intent(in)  :: extrasnow  ! snow levels are ordered top layer first
    integer                , intent(out) :: ierr       ! field with an error, or 0
    character(len=*)       , intent(out) :: errmsg
    !
    ! !LOCAL VARIABLES:
    integer :: f, n, offset, snow_offset, sample_offset, num_points, num_snow, extent
    integer :: level, averaged, max_points
    !-----------------------------------------------------------------------

    ierr = 0
    errmsg = ' '
    call release(src1, src2)

    num_fields = size(field_list)
    allocate(fields(num_fields), tape_first_field(size(tape_first)))
    fields(:) = field_list(:)
    tape_first_field(:) = tape_first(:)
    extrasnowlayers = extrasnow

    offset = 0
    snow_offset = 0
    sample_offset = 0
    num_snow = 0
    max_levels = 0
    max_points = view%endg - view%begg + 1
    do f = 1, num_fields
       associate(fld => fields(f))
       call level_bounds(view, fld%source_level, fld%source_first, fld%source_last)
       num_points = fld%source_last - fld%source_first + 1

       ! the model array covers the level's points; a 2d array has levels
       ! 1:num_levels, a snow field's levels start at its lower bound
       if (fld%source_rank == 1) then
          if (.not. associated(src1(fld%source_index)%ptr)) then
             ierr = f; errmsg = 'model array not associated'; return
          end if
          extent = size(src1(fld%source_index)%ptr)
       else
          if (.not. associated(src2(fld%source_index)%ptr)) then
             ierr = f; errmsg = 'model array not associated'; return
          end if
          extent = size(src2(fld%source_index)%ptr, 1)
          if (fld%no_snow_behavior /= no_snow_unset) then
             if (size(src2(fld%source_index)%ptr, 2) /= fld%num_levels) then
                ierr = f; errmsg = 'snow field: model array levels differ from the field levels'; return
             end if
          else if (lbound(src2(fld%source_index)%ptr, 2) > 1 .or. &
               ubound(src2(fld%source_index)%ptr, 2) < fld%num_levels) then
             ierr = f; errmsg = 'model array does not have levels 1:num2d'; return
          end if
       end if
       if (extent < num_points) then
          ierr = f; errmsg = 'model array smaller than its subgrid level'; return
       end if
       if (fld%no_snow_behavior /= no_snow_unset .and. (fld%source_rank /= 2 .or. &
            (fld%source_level /= lev_col .and. fld%source_level /= lev_patch))) then
          ierr = f; errmsg = 'snow fields must be column or patch arrays with levels'; return
       end if
       if (fld%average_to_gridcell .and. fld%source_level == lev_grid) fld%average_to_gridcell = .false.

       if (fld%no_snow_behavior /= no_snow_unset) then
          fld%snow_offset = snow_offset
          snow_offset = snow_offset + num_points * fld%num_levels
          num_snow = num_snow + 1
       end if
       if (fld%average_to_gridcell) then
          fld%first_point = view%begg
          fld%last_point  = view%endg
          fld%sample_offset = sample_offset
          sample_offset = sample_offset + (view%endg - view%begg + 1) * fld%num_levels
       else
          fld%first_point = fld%source_first
          fld%last_point  = fld%source_last
       end if
       fld%buffer_offset = offset
       offset = offset + (fld%last_point - fld%first_point + 1) * fld%num_levels
       max_levels = max(max_levels, fld%num_levels)
       max_points = max(max_points, num_points)
       end associate
    end do

    allocate(hbuf_pool(offset), nacs_pool(offset), snow_values(snow_offset), snow_fields(num_snow))
    allocate(gridcell_samples(sample_offset), all_active(max_points), process_order(num_fields))
    hbuf_pool(:) = 0._r8
    nacs_pool(:) = 0
    all_active(:) = .true.
    num_snow = 0
    do f = 1, num_fields
       if (fields(f)%snow_offset >= 0) then
          num_snow = num_snow + 1
          snow_fields(num_snow) = f
       end if
    end do
    ! processing order: grouped by subgrid level, native fields before
    ! averaged ones, in tape order within a group
    n = 0
    do level = lev_grid, lev_patch
       do averaged = 0, 1
          do f = 1, num_fields
             if (fields(f)%source_level == level .and. &
                  (fields(f)%average_to_gridcell .eqv. (averaged == 1))) then
                n = n + 1
                process_order(n) = f
             end if
          end do
       end do
    end do

    !$acc update device(fields, hbuf_pool, nacs_pool, all_active, extrasnowlayers)
    !$acc enter data copyin(tape_first_field, snow_fields, process_order)
    call attach_sources(src1, src2)

  end subroutine hist_compute_build

  !-----------------------------------------------------------------------
  pure subroutine level_bounds(view, level, first, last)
    ! First and last point of a subgrid level.
    type(subgrid_ave_view), intent(in) :: view
    integer, intent(in)  :: level
    integer, intent(out) :: first, last

    select case (level)
    case (lev_topo)
       first = view%begt ; last = view%endt
    case (lev_lun)
       first = view%begl ; last = view%endl
    case (lev_col)
       first = view%begc ; last = view%endc
    case (lev_patch)
       first = view%begp ; last = view%endp
    case default
       first = view%begg ; last = view%endg
    end select
  end subroutine level_bounds

  !-----------------------------------------------------------------------
  subroutine attach_sources(src1, src2)
    !
    ! !DESCRIPTION:
    ! Put the pointer lists on the device and attach the model arrays the
    ! fields use (their data must already be on the device).
    !
    ! !ARGUMENTS:
    type(hist_src1d), intent(in) :: src1(:)
    type(hist_src2d), intent(in) :: src2(:)
    !
    ! !LOCAL VARIABLES:
    logical :: used1(size(src1)), used2(size(src2))
    integer :: f, i
    !-----------------------------------------------------------------------

    used1(:) = .false.
    used2(:) = .false.
    do f = 1, num_fields
       if (fields(f)%source_rank == 1) then
          used1(fields(f)%source_index) = .true.
       else
          used2(fields(f)%source_index) = .true.
       end if
    end do
    attached1 = pack([(i, i = 1, size(src1))], used1)
    attached2 = pack([(i, i = 1, size(src2))], used2)

    !$acc enter data copyin(src1, src2)
    do i = 1, size(attached1)
       !$acc enter data attach(src1(attached1(i))%ptr)
    end do
    do i = 1, size(attached2)
       !$acc enter data attach(src2(attached2(i))%ptr)
    end do

  end subroutine attach_sources

  !-----------------------------------------------------------------------
  subroutine release(src1, src2)
    ! Discard the buffers and the field table (and their device copies).
    type(hist_src1d), intent(in) :: src1(:)
    type(hist_src2d), intent(in) :: src2(:)
    integer :: i

    if (.not. allocated(fields)) return
    do i = 1, size(attached1)
       !$acc exit data detach(src1(attached1(i))%ptr)
    end do
    do i = 1, size(attached2)
       !$acc exit data detach(src2(attached2(i))%ptr)
    end do
    !$acc exit data delete(src1, src2)
    !$acc exit data delete(tape_first_field, snow_fields, process_order)
    deallocate(fields, tape_first_field, snow_fields, process_order, hbuf_pool, nacs_pool)
    deallocate(snow_values, gridcell_samples, all_active)
    deallocate(attached1, attached2)
    num_fields = 0

  end subroutine release

  !-----------------------------------------------------------------------
  subroutine hist_compute_views(f, hbuf, nacs)
    !
    ! !DESCRIPTION:
    ! Pointer views of field f's buffers: [output point, level].
    !
    ! !ARGUMENTS:
    integer , intent(in) :: f
    real(r8), pointer, intent(out) :: hbuf(:,:)
    integer , pointer, intent(out) :: nacs(:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: first, last, n
    !-----------------------------------------------------------------------

    first = fields(f)%first_point
    last  = fields(f)%last_point
    n = (last - first + 1) * fields(f)%num_levels
    hbuf(first:last, 1:fields(f)%num_levels) => hbuf_pool(fields(f)%buffer_offset+1 : fields(f)%buffer_offset+n)
    nacs(first:last, 1:fields(f)%num_levels) => nacs_pool(fields(f)%buffer_offset+1 : fields(f)%buffer_offset+n)

  end subroutine hist_compute_views

  !-----------------------------------------------------------------------
  subroutine hist_compute_update(view, src1, src2)
    !
    ! !DESCRIPTION:
    ! Accumulate the current values of all active fields into their buffers.
    !
    ! !ARGUMENTS:
    type(subgrid_ave_view), intent(in)  :: view
    type(hist_src1d)      , intent(in)  :: src1(:)
    type(hist_src2d)      , intent(in)  :: src2(:)
    !
    ! !LOCAL VARIABLES:
    integer :: n, f, s, lev, i
    !-----------------------------------------------------------------------

    ! snow fields: fill in the absent layers (the averages read plain value arrays)
    !$acc parallel loop gang collapse(2) default(present) private(f)
    do s = 1, size(snow_fields)
       do lev = 1, max_levels
          f = snow_fields(s)
          if (lev > fields(f)%num_levels) cycle
          !$acc loop vector
          do i = fields(f)%source_first, fields(f)%source_last
             snow_values(snow_index(fields(f), i, lev)) = &
                  snow_layer_value(view, fields(f), src2(fields(f)%source_index)%ptr, i, lev)
          end do
       end do
    end do

    ! every field: one gang per field, its levels in turn
    !$acc parallel loop gang default(present) private(f)
    do n = 1, num_fields
       f = process_order(n)
       !$acc loop seq
       do lev = 1, fields(f)%num_levels
          if (fields(f)%snow_offset >= 0) then
             call field_level(view, fields(f), lev, snow_values(snow_index(fields(f), fields(f)%source_first, lev)), &
                  all_active, gridcell_samples, hbuf_pool, nacs_pool)
          else if (fields(f)%source_rank == 1) then
             call field_level(view, fields(f), lev, src1(fields(f)%source_index)%ptr, &
                  all_active, gridcell_samples, hbuf_pool, nacs_pool)
          else
             call field_level(view, fields(f), lev, src2(fields(f)%source_index)%ptr(:, lev), &
                  all_active, gridcell_samples, hbuf_pool, nacs_pool)
          end if
       end do
    end do

  end subroutine hist_compute_update

  !-----------------------------------------------------------------------
  pure integer function snow_index(fld, i, lev)
    ! Position of point i, level lev of a snow field in snow_values.
    !$acc routine seq
    type(hist_active_field), intent(in) :: fld
    integer, intent(in) :: i, lev
    snow_index = fld%snow_offset + (lev - 1) * (fld%source_last - fld%source_first + 1) &
         + (i - fld%source_first + 1)
  end function snow_index

  !-----------------------------------------------------------------------
  pure real(r8) function snow_layer_value(view, fld, model_array, i, lev)
    !
    ! !DESCRIPTION:
    ! Value of a snow field at point i (column or patch), level lev, with the
    ! levels of absent snow layers set to the no-snow value. Levels are ordered
    ! bottom-up with the absent layers first, unless extrasnowlayers, where
    ! the top layer comes first and the absent layers last.
    !
    !$acc routine seq
    ! !ARGUMENTS:
    type(subgrid_ave_view) , intent(in) :: view
    type(hist_active_field), intent(in) :: fld
    real(r8)               , intent(in) :: model_array(:,:)  ! [point of the level, level from 1]
    integer                , intent(in) :: i, lev
    !
    ! !LOCAL VARIABLES:
    integer  :: c, point_in_level, num_absent
    real(r8) :: no_snow_value
    !-----------------------------------------------------------------------

    if (fld%source_level == lev_patch) then
       c = view%col_of_pft(i)
    else
       c = i
    end if
    num_absent = fld%num_levels - abs(view%snl(c))
    no_snow_value = spval
    if (fld%no_snow_behavior == no_snow_zero) no_snow_value = 0._r8
    point_in_level = i - fld%source_first + 1

    if (.not. extrasnowlayers) then
       if (lev <= num_absent) then
          snow_layer_value = no_snow_value
       else
          snow_layer_value = model_array(point_in_level, lev)
       end if
    else
       if (lev > fld%num_levels - num_absent) then
          snow_layer_value = no_snow_value
       else
          snow_layer_value = model_array(point_in_level, lev + num_absent)
       end if
    end if

  end function snow_layer_value

  !-----------------------------------------------------------------------
  pure subroutine field_level(view, fld, lev, values, all_active, gridcell_samples, hbuf_pool, nacs_pool)
    !
    ! !DESCRIPTION:
    ! Accumulate level lev of field fld: average it to gridcells first if the
    ! field is averaged, then apply its time operation to all output points.
    ! The scratch and buffer pools are passed in (pure: no module state is
    ! changed); the field's part of them is found from its offsets.
    !
    !$acc routine vector
    ! !ARGUMENTS:
    type(subgrid_ave_view) , intent(in)    :: view
    type(hist_active_field), intent(in)    :: fld
    integer                , intent(in)    :: lev
    real(r8)               , intent(in)    :: values(*)           ! the level's values; values(1) is the first point of the source level
    logical                , intent(in)    :: all_active(:)       ! .true. (points without an active check)
    real(r8)               , intent(inout) :: gridcell_samples(*) ! averages of the averaged fields
    real(r8)               , intent(inout) :: hbuf_pool(*)
    integer                , intent(inout) :: nacs_pool(*)
    !
    ! !LOCAL VARIABLES:
    integer :: num_points, buffer_first, sample_first
    !-----------------------------------------------------------------------

    num_points   = fld%last_point - fld%first_point + 1
    buffer_first = fld%buffer_offset + (lev - 1) * num_points + 1

    if (fld%average_to_gridcell) then
       sample_first = fld%sample_offset + (lev - 1) * num_points + 1
       select case (fld%source_level)
       case (lev_patch)
          if (fld%c2l_scale == scale_unity .and. fld%l2g_scale == scale_unity) then
             call p2g_range_unity(view, view%begg, view%endg, values, view%begp, gridcell_samples(sample_first))
          else
             call p2g_range(view, view%begg, view%endg, values, view%begp, fld%c2l_scale, fld%l2g_scale, &
                  gridcell_samples(sample_first))
          end if
       case (lev_col)
          if (fld%c2l_scale == scale_unity .and. fld%l2g_scale == scale_unity) then
             call c2g_range_unity(view, view%begg, view%endg, values, view%begc, gridcell_samples(sample_first))
          else
             call c2g_range(view, view%begg, view%endg, values, view%begc, fld%c2l_scale, fld%l2g_scale, &
                  gridcell_samples(sample_first))
          end if
       case (lev_lun)
          call l2g_range(view, view%begg, view%endg, values, view%begl, fld%l2g_scale, gridcell_samples(sample_first))
       case (lev_topo)
          call t2g_range(view, view%begg, view%endg, values, view%begt, gridcell_samples(sample_first))
       end select
       call accumulate_samples(fld%time_op, .true., num_points, gridcell_samples(sample_first), &
            all_active, hbuf_pool(buffer_first), nacs_pool(buffer_first))
    else
       select case (fld%source_level)
       case (lev_topo)
          call accumulate_samples(fld%time_op, .false., num_points, values, view%active_t, &
               hbuf_pool(buffer_first), nacs_pool(buffer_first))
       case (lev_lun)
          call accumulate_samples(fld%time_op, .false., num_points, values, view%active_l, &
               hbuf_pool(buffer_first), nacs_pool(buffer_first))
       case (lev_col)
          call accumulate_samples(fld%time_op, .false., num_points, values, view%active_c, &
               hbuf_pool(buffer_first), nacs_pool(buffer_first))
       case (lev_patch)
          call accumulate_samples(fld%time_op, .false., num_points, values, view%active_p, &
               hbuf_pool(buffer_first), nacs_pool(buffer_first))
       case default
          call accumulate_samples(fld%time_op, .false., num_points, values, all_active, &
               hbuf_pool(buffer_first), nacs_pool(buffer_first))
       end select
    end if

  end subroutine field_level

  !-----------------------------------------------------------------------
  pure subroutine accumulate_samples(time_op, gridcell_rule, n, samples, active, hbuf, nacs)
    !
    ! !DESCRIPTION:
    ! Accumulate n samples into their buffers (rules in the module
    ! description). gridcell_rule selects the 'X'/'M' rules of fields averaged
    ! to gridcells.
    !
    !$acc routine vector
    ! !ARGUMENTS:
    integer , intent(in)    :: time_op
    logical , intent(in)    :: gridcell_rule
    integer , intent(in)    :: n
    real(r8), intent(in)    :: samples(*)
    logical , intent(in)    :: active(:)   ! active(i): point i is active
    real(r8), intent(inout) :: hbuf(*)
    integer , intent(inout) :: nacs(*)
    !
    ! !LOCAL VARIABLES:
    integer :: i
    logical :: valid
    !-----------------------------------------------------------------------

    select case (time_op)
    case (time_inst)
       !$acc loop vector private(valid)
       do i = 1, n
          valid = active(i) .and. samples(i) /= spval
          if (valid) then
             hbuf(i) = samples(i)
          else
             hbuf(i) = spval
          end if
          nacs(i) = 1
       end do

    case (time_mean)
       !$acc loop vector private(valid)
       do i = 1, n
          valid = active(i) .and. samples(i) /= spval
          if (valid) then
             if (nacs(i) == 0) hbuf(i) = 0._r8
             hbuf(i) = hbuf(i) + samples(i)
             nacs(i) = nacs(i) + 1
          else if (nacs(i) == 0) then
             hbuf(i) = spval
          end if
       end do

    case (time_max)
       if (gridcell_rule) then
          !$acc loop vector private(valid)
          do i = 1, n
             valid = active(i) .and. samples(i) /= spval
             if (valid) then
                if (nacs(i) == 0) hbuf(i) = -big
                hbuf(i) = max(hbuf(i), samples(i))
             else
                hbuf(i) = spval
             end if
             nacs(i) = 1
          end do
       else
          !$acc loop vector private(valid)
          do i = 1, n
             valid = active(i) .and. samples(i) /= spval
             if (valid) then
                if (nacs(i) == 0) hbuf(i) = -big
                hbuf(i) = max(hbuf(i), samples(i))
             else if (nacs(i) == 0) then
                hbuf(i) = spval
             end if
             nacs(i) = 1
          end do
       end if

    case (time_min)
       if (gridcell_rule) then
          !$acc loop vector private(valid)
          do i = 1, n
             valid = active(i) .and. samples(i) /= spval
             if (valid) then
                if (nacs(i) == 0) hbuf(i) = big
                hbuf(i) = min(hbuf(i), samples(i))
             else
                hbuf(i) = spval
             end if
             nacs(i) = 1
          end do
       else
          !$acc loop vector private(valid)
          do i = 1, n
             valid = active(i) .and. samples(i) /= spval
             if (valid) then
                if (nacs(i) == 0) hbuf(i) = big
                hbuf(i) = min(hbuf(i), samples(i))
             else if (nacs(i) == 0) then
                hbuf(i) = spval
             end if
             nacs(i) = 1
          end do
       end if

    case (time_max_valid)
       !$acc loop vector private(valid)
       do i = 1, n
          valid = active(i) .and. samples(i) /= spval
          if (valid) then
             if (nacs(i) == 0) hbuf(i) = -big
             hbuf(i) = max(hbuf(i), samples(i))
             nacs(i) = nacs(i) + 1
          else if (nacs(i) == 0) then
             hbuf(i) = spval
          end if
       end do

    case (time_min_valid)
       !$acc loop vector private(valid)
       do i = 1, n
          valid = active(i) .and. samples(i) /= spval
          if (valid) then
             if (nacs(i) == 0) hbuf(i) = big
             hbuf(i) = min(hbuf(i), samples(i))
             nacs(i) = nacs(i) + 1
          else if (nacs(i) == 0) then
             hbuf(i) = spval
          end if
       end do
    end select

  end subroutine accumulate_samples

  !-----------------------------------------------------------------------
  subroutine hist_compute_normalize(t)
    !
    ! !DESCRIPTION:
    ! Divide the time means of tape t by their number of samples.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t
    !
    ! !LOCAL VARIABLES:
    integer :: f, n, first, last
    !-----------------------------------------------------------------------

    !$acc parallel loop gang default(present) private(first, last)
    do f = tape_first_field(t), tape_first_field(t+1) - 1
       if (fields(f)%time_op /= time_mean) cycle
       first = fields(f)%buffer_offset + 1
       last  = fields(f)%buffer_offset + (fields(f)%last_point - fields(f)%first_point + 1) * fields(f)%num_levels
       !$acc loop vector
       do n = first, last
          if (nacs_pool(n) /= 0) hbuf_pool(n) = hbuf_pool(n) / real(nacs_pool(n), r8)
       end do
    end do

  end subroutine hist_compute_normalize

  !-----------------------------------------------------------------------
  subroutine hist_compute_zero(t)
    !
    ! !DESCRIPTION:
    ! Zero the buffers and counters of tape t.
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t
    !
    ! !LOCAL VARIABLES:
    integer :: n, first, last
    !-----------------------------------------------------------------------

    call tape_range(t, first, last)
    !$acc parallel loop default(present)
    do n = first, last
       hbuf_pool(n) = 0._r8
       nacs_pool(n) = 0
    end do

  end subroutine hist_compute_zero

  !-----------------------------------------------------------------------
  subroutine hist_compute_update_host(t)
    !
    ! !DESCRIPTION:
    ! Copy tape t's buffers to the host before they are written (OpenACC).
    !
    ! !ARGUMENTS:
    integer, intent(in) :: t
    !
    ! !LOCAL VARIABLES:
    integer :: first, last
    !-----------------------------------------------------------------------

    call tape_range(t, first, last)
    if (last < first) return
    !$acc update self(hbuf_pool(first:last), nacs_pool(first:last))

  end subroutine hist_compute_update_host

  !-----------------------------------------------------------------------
  subroutine hist_compute_update_device()
    !
    ! !DESCRIPTION:
    ! Copy all buffers to the device after they were read on the host
    ! (restart; OpenACC).
    !-----------------------------------------------------------------------

    !$acc update device(hbuf_pool, nacs_pool)

  end subroutine hist_compute_update_device

  !-----------------------------------------------------------------------
  pure subroutine tape_range(t, first, last)
    ! The positions of tape t's buffers in the pool (last < first if none).
    integer, intent(in)  :: t
    integer, intent(out) :: first, last
    integer :: f_last

    first = 1
    last  = 0
    if (tape_first_field(t+1) <= tape_first_field(t)) return
    f_last = tape_first_field(t+1) - 1
    first = fields(tape_first_field(t))%buffer_offset + 1
    last  = fields(f_last)%buffer_offset + &
         (fields(f_last)%last_point - fields(f_last)%first_point + 1) * fields(f_last)%num_levels
  end subroutine tape_range

end module histComputeMod
