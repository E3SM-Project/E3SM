module shr_bounds_mod

  !---------------------------------------------------------------------
  ! Bounds checks on fields a component sends to or receives from the
  ! coupler.
  !
  ! Limits come from the limits dictionary (share/field_limits.yaml,
  ! generated into shr_field_limits_mod). A component calls these routines
  ! from its coupler interface layer, on the plain arrays it exchanges with
  ! the coupler, so the check works under both the MCT and MOAB drivers:
  !
  !   shr_bounds_init   : once, with the driver's field list and the
  !                       component's grid points
  !   shr_bounds_check  : every coupling step, on the export array after it
  !                       is filled or the import array before it is used;
  !                       keeps running statistics
  !   shr_bounds_report : once per model day and at the end of the run;
  !                       writes the statistics to the component log and
  !                       resets them
  !   shr_bounds_require: once, from the driver's seq_flds_set; stops the
  !                       model if a field a component sends has no entry
  !
  ! An export check applies the limits of entries the component itself
  ! sends. An import check applies the limits of every entry whose field
  ! the component receives, as set by the sending component. Entries that
  ! are unbounded or not yet reviewed have no limits and are not checked.
  !---------------------------------------------------------------------

  use shr_kind_mod   , only: r8 => SHR_KIND_R8, CL => SHR_KIND_CL, CS => SHR_KIND_CS
  use shr_sys_mod    , only: shr_sys_abort, shr_sys_flush
  use shr_infnan_mod , only: shr_infnan_isnan
  use shr_const_mod  , only: shr_const_isspval
  use shr_string_mod , only: shr_string_listGetNum, shr_string_listGetName, shr_string_listGetIndexF
  use shr_field_limits_mod

  implicit none
#include <mpif.h>
  private

  public :: shr_bounds_type
  public :: shr_bounds_init
  public :: shr_bounds_check
  public :: shr_bounds_report
  public :: shr_bounds_require

  type shr_bounds_type
     logical               :: initialized = .false.
     character(len=16)     :: label           ! e.g. 'lnd export'
     integer               :: nflds = 0       ! number of fields in the exchanged array
     ! per field in the exchanged array
     integer , allocatable :: ent(:)     ! dictionary entry of this field, 0 if none
     integer , allocatable :: ient(:)    ! entry whose limits are checked on this field, 0 if none
     integer , allocatable :: kwhere(:)  ! index of the where_positive field, 0 if none
     character(len=CS), allocatable :: name(:)
     real(r8), allocatable :: vmin(:)    ! smallest value since the last report
     real(r8), allocatable :: vmax(:)    ! largest value since the last report
     integer , allocatable :: nmin(:)    ! local point of vmin, 0 if no point checked
     integer , allocatable :: nmax(:)    ! local point of vmax, 0 if no point checked
     real(r8), allocatable :: nout(:)    ! number of out-of-bounds values since the last report
     ! per local point
     logical , allocatable :: active(:)  ! point is checked
     real(r8), allocatable :: gid(:)     ! global index
     real(r8), allocatable :: lat(:)
     real(r8), allocatable :: lon(:)
  end type shr_bounds_type

  character(len=*), parameter :: F01 = &
       "('bounds_check: ',a,1x,a16,' min ',es14.6,' (gid ',i9,' lat ',f7.2,' lon ',f7.2,')', &
       &' max ',es14.6,' (gid ',i9,' lat ',f7.2,' lon ',f7.2,')',' n_out ',i10,2x,a)"
  character(len=*), parameter :: F02 = "('bounds_check: ',a,1x,a16,' no points checked',2x,a)"

!===============================================================================
contains
!===============================================================================

  subroutine shr_bounds_init(bnd, comp, direction, fieldlist, gid, lat, lon, &
       mpicom, logunit, mask)

    ! Match the dictionary entries to the fields in fieldlist and store the
    ! grid points. Writes the fields that are not checked to logunit from
    ! the root of mpicom. Must be called by every task of mpicom.

    type(shr_bounds_type), intent(inout) :: bnd
    character(len=*)     , intent(in)    :: comp       ! component type, e.g. 'lnd'
    character(len=*)     , intent(in)    :: direction  ! 'export' or 'import'
    character(len=*)     , intent(in)    :: fieldlist  ! colon-separated field names, e.g. seq_flds_l2x_fields
    integer              , intent(in)    :: gid(:)     ! global index of each local point
    real(r8)             , intent(in)    :: lat(:)     ! latitude of each local point (degrees)
    real(r8)             , intent(in)    :: lon(:)     ! longitude of each local point (degrees)
    integer              , intent(in)    :: mpicom     ! component communicator
    integer              , intent(in)    :: logunit    ! component log
    logical, optional    , intent(in)    :: mask(:)    ! check only points where true (default: all)

    integer :: i, k, npts, rank, ierr
    logical :: export
    character(len=*), parameter :: subname = '(shr_bounds_init) '

    if (direction /= 'export' .and. direction /= 'import') then
       call shr_sys_abort(subname//"direction must be 'export' or 'import', not '"//trim(direction)//"'")
    end if
    export = (direction == 'export')
    bnd%label = trim(comp)//' '//trim(direction)

    bnd%nflds = shr_string_listGetNum(fieldlist)
    allocate(bnd%name(bnd%nflds), bnd%ent(bnd%nflds), bnd%ient(bnd%nflds), bnd%kwhere(bnd%nflds))
    allocate(bnd%vmin(bnd%nflds), bnd%vmax(bnd%nflds))
    allocate(bnd%nmin(bnd%nflds), bnd%nmax(bnd%nflds))
    allocate(bnd%nout(bnd%nflds))
    do k = 1, bnd%nflds
       call shr_string_listGetName(fieldlist, k, bnd%name(k))
    end do

    ! fields whose entry has limits (and whose where_positive field is in
    ! the list); an export checks only the entries the component itself sends
    bnd%ient = 0
    bnd%kwhere = 0
    do k = 1, bnd%nflds
       i = shr_bounds_entry(bnd%name(k))
       bnd%ent(k) = i
       if (i == 0) cycle
       if (shr_field_limits(i)%status /= 'limited') cycle
       if (export .and. shr_field_limits(i)%component /= comp) cycle
       if (len_trim(shr_field_limits(i)%where_positive) > 0) then
          bnd%kwhere(k) = shr_string_listGetIndexF(fieldlist, trim(shr_field_limits(i)%where_positive))
          if (bnd%kwhere(k) == 0) cycle
       end if
       bnd%ient(k) = i
    end do

    npts = size(gid)
    if (size(lat) /= npts .or. size(lon) /= npts) then
       call shr_sys_abort(subname//trim(bnd%label)//': gid, lat and lon differ in size')
    end if
    allocate(bnd%active(npts), bnd%gid(npts), bnd%lat(npts), bnd%lon(npts))
    bnd%active = .true.
    if (present(mask)) then
       if (size(mask) /= npts) call shr_sys_abort(subname//trim(bnd%label)//': mask has the wrong size')
       bnd%active = mask
    end if
    bnd%gid = real(gid, r8)
    bnd%lat = lat
    bnd%lon = lon

    ! list the fields that are not checked
    call mpi_comm_rank(mpicom, rank, ierr)
    if (rank == 0) then
       call write_names(bnd, logunit, 'unbounded' , ' unbounded:')
       call write_names(bnd, logunit, 'unreviewed', ' unreviewed:')
       call write_names(bnd, logunit, 'limited'   , ' limits not applied here:')
       call write_names(bnd, logunit, ''          , ' no entry:')
       write(logunit,'(a,a,i4,a,i4,a)') 'bounds_check: ', trim(bnd%label), &
            count(bnd%ient > 0), ' of ', bnd%nflds, ' fields are checked'
       call shr_sys_flush(logunit)
    end if

    call shr_bounds_reset(bnd)
    bnd%initialized = .true.

  end subroutine shr_bounds_init

  !===============================================================================

  subroutine shr_bounds_check(bnd, values, field_major, abort_on_violation)

    ! Check the exchanged array against its limits, update the statistics,
    ! and abort on the first out-of-bounds value if abort_on_violation is
    ! true. NaNs and special values are skipped; NaNs are check_fields' job.

    type(shr_bounds_type), intent(inout) :: bnd
    real(r8)             , intent(in)    :: values(:,:)
    logical              , intent(in)    :: field_major ! true: values(field,point); false: values(point,field)
    logical              , intent(in)    :: abort_on_violation

    integer  :: i, k, kw, n, npts, nflds
    real(r8) :: v, w
    logical  :: out
    character(len=CL) :: msg
    character(len=*), parameter :: subname = '(shr_bounds_check) '

    if (.not. bnd%initialized) call shr_sys_abort(subname//'called before shr_bounds_init')

    npts = size(bnd%active)
    if (field_major) then
       nflds = size(values,1)
       if (size(values,2) /= npts) nflds = -1
    else
       nflds = size(values,2)
       if (size(values,1) /= npts) nflds = -1
    end if
    if (nflds < bnd%nflds) then
       write(msg,'(a,a,a,2i8,a,2i8)') subname, trim(bnd%label), ': array shape ', &
            shape(values), ' does not match fields x points ', bnd%nflds, npts
       call shr_sys_abort(trim(msg))
    end if

    do k = 1, bnd%nflds
       i = bnd%ient(k)
       if (i == 0) cycle
       kw = bnd%kwhere(k)
       do n = 1, npts
          if (.not. bnd%active(n)) cycle
          if (field_major) then
             v = values(k,n)
          else
             v = values(n,k)
          end if
          if (shr_infnan_isnan(v)) cycle
          if (shr_const_isspval(v)) cycle
          if (kw > 0) then
             if (field_major) then
                w = values(kw,n)
             else
                w = values(n,kw)
             end if
             if (.not. (w > 0.0_r8)) cycle
          end if

          if (bnd%nmin(k) == 0 .or. v < bnd%vmin(k)) then
             bnd%vmin(k) = v
             bnd%nmin(k) = n
          end if
          if (bnd%nmax(k) == 0 .or. v > bnd%vmax(k)) then
             bnd%vmax(k) = v
             bnd%nmax(k) = n
          end if

          out = .false.
          if (shr_field_limits(i)%has_min) then
             if (shr_field_limits(i)%min_inclusive) then
                out = out .or. v < shr_field_limits(i)%min_value
             else
                out = out .or. v <= shr_field_limits(i)%min_value
             end if
          end if
          if (shr_field_limits(i)%has_max) then
             if (shr_field_limits(i)%max_inclusive) then
                out = out .or. v > shr_field_limits(i)%max_value
             else
                out = out .or. v >= shr_field_limits(i)%max_value
             end if
          end if

          if (out) then
             bnd%nout(k) = bnd%nout(k) + 1.0_r8
             if (abort_on_violation) then
                write(msg,'(a,a,1x,a,a,es14.6,a,i9,a,f7.2,a,f7.2,a,a)') subname, trim(bnd%label), &
                     trim(bnd%name(k)), ' = ', v, ' at gid ', nint(bnd%gid(n)), &
                     ' lat ', bnd%lat(n), ' lon ', bnd%lon(n), ' outside ', &
                     trim(limits_string(i))
                call shr_sys_abort(trim(msg))
             end if
          end if
       end do
    end do

  end subroutine shr_bounds_check

  !===============================================================================

  subroutine shr_bounds_report(bnd, ymd, tod, mpicom, logunit)

    ! Reduce the statistics over mpicom, write them to logunit from its
    ! root, and reset them. Must be called by every task of mpicom.

    type(shr_bounds_type), intent(inout) :: bnd
    integer              , intent(in)    :: ymd      ! model date (yyyymmdd)
    integer              , intent(in)    :: tod      ! model time of day (s)
    integer              , intent(in)    :: mpicom   ! component communicator
    integer              , intent(in)    :: logunit  ! component log

    integer  :: i, k, n, rank, ierr
    real(r8) :: gnout(bnd%nflds)
    real(r8) :: lmin(2,bnd%nflds), gmin(2,bnd%nflds), lmax(2,bnd%nflds), gmax(2,bnd%nflds)
    real(r8) :: lloc(6,bnd%nflds), gloc(6,bnd%nflds)

    if (.not. bnd%initialized) return
    ! the field list, and so the checked fields, are the same on every task
    if (all(bnd%ient == 0)) return

    call mpi_comm_rank(mpicom, rank, ierr)

    lmin(1,:) =  huge(1.0_r8)
    lmax(1,:) = -huge(1.0_r8)
    lmin(2,:) = real(rank, r8)
    lmax(2,:) = real(rank, r8)
    do k = 1, bnd%nflds
       if (bnd%nmin(k) > 0) lmin(1,k) = bnd%vmin(k)
       if (bnd%nmax(k) > 0) lmax(1,k) = bnd%vmax(k)
    end do

    call mpi_allreduce(bnd%nout, gnout, bnd%nflds, MPI_REAL8, MPI_SUM, mpicom, ierr)
    call mpi_allreduce(lmin, gmin, bnd%nflds, MPI_2DOUBLE_PRECISION, MPI_MINLOC, mpicom, ierr)
    call mpi_allreduce(lmax, gmax, bnd%nflds, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, mpicom, ierr)

    ! the task holding each extreme contributes its location
    lloc = 0.0_r8
    do k = 1, bnd%nflds
       n = bnd%nmin(k)
       if (n > 0 .and. nint(gmin(2,k)) == rank) &
            lloc(1:3,k) = [bnd%gid(n), bnd%lat(n), bnd%lon(n)]
       n = bnd%nmax(k)
       if (n > 0 .and. nint(gmax(2,k)) == rank) &
            lloc(4:6,k) = [bnd%gid(n), bnd%lat(n), bnd%lon(n)]
    end do
    call mpi_allreduce(lloc, gloc, 6*bnd%nflds, MPI_REAL8, MPI_SUM, mpicom, ierr)

    if (rank == 0) then
       write(logunit,'(a,a,a,i10,i6)') 'bounds_check: ', trim(bnd%label), ' model date = ', ymd, tod
       do k = 1, bnd%nflds
          i = bnd%ient(k)
          if (i == 0) cycle
          if (gmin(1,k) > gmax(1,k)) then
             write(logunit,F02) trim(bnd%label), bnd%name(k), trim(limits_string(i))
          else
             write(logunit,F01) trim(bnd%label), bnd%name(k), &
                  gmin(1,k), nint(gloc(1,k)), gloc(2,k), gloc(3,k), &
                  gmax(1,k), nint(gloc(4,k)), gloc(5,k), gloc(6,k), &
                  nint(gnout(k)), trim(limits_string(i))
          end if
       end do
       call shr_sys_flush(logunit)
    end if
    call shr_bounds_reset(bnd)

  end subroutine shr_bounds_report

  !===============================================================================

  subroutine shr_bounds_reset(bnd)

    type(shr_bounds_type), intent(inout) :: bnd

    bnd%vmin = 0.0_r8
    bnd%vmax = 0.0_r8
    bnd%nmin = 0
    bnd%nmax = 0
    bnd%nout = 0.0_r8

  end subroutine shr_bounds_reset

  !===============================================================================

  subroutine shr_bounds_require(comps, fieldlists, exempt, logunit, write_log)

    ! Stop the model if a field that a component sends has no entry in the
    ! limits dictionary, or an entry for another component. Fields in
    ! exempt (user-defined fields, e.g. cplflds_custom) are skipped. Writes
    ! each missing field, and a summary per component, to logunit if
    ! write_log.

    character(len=*), intent(in) :: comps(:)       ! component types, e.g. 'lnd'
    character(len=*), intent(in) :: fieldlists(:)  ! colon-separated fields each component sends
    character(len=*), intent(in) :: exempt         ! colon-separated fields to skip
    integer         , intent(in) :: logunit
    logical         , intent(in) :: write_log

    integer :: i, k, n, nfld, nlimited, nunbounded, nunreviewed, nexempt, nmissing, ntotal
    character(len=CS) :: name
    character(len=CL) :: msg
    character(len=*), parameter :: subname = '(shr_bounds_require) '

    ntotal = 0
    do n = 1, size(comps)
       nfld = shr_string_listGetNum(fieldlists(n))
       nlimited = 0
       nunbounded = 0
       nunreviewed = 0
       nexempt = 0
       nmissing = 0
       do k = 1, nfld
          call shr_string_listGetName(fieldlists(n), k, name)
          i = shr_bounds_entry(name)
          msg = ''
          if (shr_string_listGetIndexF(exempt, trim(name)) > 0) then
             nexempt = nexempt + 1
          else if (i == 0) then
             nmissing = nmissing + 1
             msg = trim(name)//' has no entry'
          else if (shr_field_limits(i)%component /= comps(n)) then
             nmissing = nmissing + 1
             msg = trim(name)//' is listed as sent by '//trim(shr_field_limits(i)%component)
          else if (shr_field_limits(i)%status == 'limited') then
             nlimited = nlimited + 1
          else if (shr_field_limits(i)%status == 'unbounded') then
             nunbounded = nunbounded + 1
          else
             nunreviewed = nunreviewed + 1
          end if
          if (write_log .and. len_trim(msg) > 0) &
               write(logunit,'(a,a,a,a)') subname, trim(comps(n)), ' sends ', trim(msg)
       end do
       ntotal = ntotal + nmissing
       if (write_log .and. nfld > 0) then
          write(logunit,'(a,a,a,i4,a,5(i4,a))') 'bounds_check: ', trim(comps(n)), ' sends', nfld, &
               ' fields:', nlimited, ' limited,', nunbounded, ' unbounded,', nunreviewed, &
               ' unreviewed,', nexempt, ' user-defined,', nmissing, ' missing'
       end if
    end do
    if (write_log) call shr_sys_flush(logunit)

    if (ntotal > 0) then
       call shr_sys_abort(subname//'fields sent to the coupler are missing from the limits '// &
            'dictionary (see log). Add each to share/field_limits.yaml with min/max or '// &
            '"unbounded: <reason>", then run share/gen_field_limits.py.')
    end if

  end subroutine shr_bounds_require

  !===============================================================================

  integer function shr_bounds_entry(name)

    ! Index of the dictionary entry for field name, 0 if none. An exact
    ! entry wins over a pattern; a pattern 'base*' matches base followed by
    ! one or more digits.

    character(len=*), intent(in) :: name

    integer :: i, nb

    shr_bounds_entry = 0
    do i = 1, shr_field_limits_nflds
       if (trim(shr_field_limits(i)%name) == trim(name)) then
          shr_bounds_entry = i
          return
       end if
    end do
    do i = 1, shr_field_limits_nflds
       nb = len_trim(shr_field_limits(i)%name) - 1
       if (shr_field_limits(i)%name(nb+1:nb+1) /= '*') cycle
       if (len_trim(name) <= nb) cycle
       if (name(1:nb) /= shr_field_limits(i)%name(1:nb)) cycle
       if (verify(trim(name(nb+1:)), '0123456789') /= 0) cycle
       shr_bounds_entry = i
       return
    end do

  end function shr_bounds_entry

  !===============================================================================

  subroutine write_names(bnd, logunit, status, label)

    ! Write the fields of bnd that are not checked and whose entry has the
    ! given status ('': no entry), wrapped at 100 characters

    type(shr_bounds_type), intent(in) :: bnd
    integer              , intent(in) :: logunit
    character(len=*)     , intent(in) :: status
    character(len=*)     , intent(in) :: label

    integer :: k
    character(len=CL) :: line
    character(len=16) :: kstatus

    line = ''
    do k = 1, bnd%nflds
       if (bnd%ient(k) > 0) cycle
       kstatus = ''
       if (bnd%ent(k) > 0) kstatus = shr_field_limits(bnd%ent(k))%status
       if (kstatus /= status) cycle
       if (len_trim(line) + len_trim(bnd%name(k)) + 1 > 100) then
          write(logunit,'(a,a,a,a)') 'bounds_check: ', trim(bnd%label), label, trim(line)
          line = ''
       end if
       line = trim(line)//' '//trim(bnd%name(k))
    end do
    if (len_trim(line) > 0) &
         write(logunit,'(a,a,a,a)') 'bounds_check: ', trim(bnd%label), label, trim(line)

  end subroutine write_names

  !===============================================================================

  function limits_string(i) result(str)

    ! Limits of entry i in interval notation, e.g. "(0, inf)" or "[0, 1]"

    integer, intent(in) :: i
    character(len=64)   :: str
    character(len=24)   :: lo, hi

    lo = '-inf'
    hi = 'inf'
    if (shr_field_limits(i)%has_min) lo = number_string(shr_field_limits(i)%min_value)
    if (shr_field_limits(i)%has_max) hi = number_string(shr_field_limits(i)%max_value)
    str = 'limits ('
    if (shr_field_limits(i)%has_min .and. shr_field_limits(i)%min_inclusive) str = 'limits ['
    str = trim(str)//trim(lo)//', '//trim(hi)
    if (shr_field_limits(i)%has_max .and. shr_field_limits(i)%max_inclusive) then
       str = trim(str)//']'
    else
       str = trim(str)//')'
    end if
    str = trim(str)//' '//trim(shr_field_limits(i)%units)

  end function limits_string

  !===============================================================================

  function number_string(v) result(str)

    real(r8), intent(in) :: v
    character(len=24)    :: str

    if (abs(v) < 1.0e9_r8 .and. v == real(nint(v), r8)) then
       write(str,'(i0)') nint(v)
    else
       write(str,'(es12.5)') v
       str = adjustl(str)
    end if

  end function number_string

end module shr_bounds_mod
