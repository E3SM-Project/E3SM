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
  !
  ! An export check applies the limits of entries the component itself
  ! sends. An import check applies the limits of every entry whose field
  ! the component receives, as set by the sending component.
  !---------------------------------------------------------------------

  use shr_kind_mod   , only: r8 => SHR_KIND_R8, CL => SHR_KIND_CL, CS => SHR_KIND_CS
  use shr_sys_mod    , only: shr_sys_abort, shr_sys_flush
  use shr_infnan_mod , only: shr_infnan_isnan
  use shr_const_mod  , only: shr_const_isspval
  use shr_string_mod , only: shr_string_listGetNum, shr_string_listGetName
  use shr_field_limits_mod

  implicit none
#include <mpif.h>
  private

  public :: shr_bounds_type
  public :: shr_bounds_init
  public :: shr_bounds_check
  public :: shr_bounds_report

  type shr_bounds_type
     logical               :: initialized = .false.
     character(len=16)     :: label           ! e.g. 'lnd export'
     integer               :: nflds = 0       ! number of fields in the exchanged array
     ! per dictionary entry
     integer , allocatable :: kfld(:)    ! index in the field list, 0 if the entry is not checked
     integer , allocatable :: kwhere(:)  ! index of the where_positive field, 0 if none
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
    ! grid points. Writes the fields without limits to logunit from the
    ! root of mpicom. Must be called by every task of mpicom.

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

    integer :: i, k, npts, rank, ierr, nnone
    logical :: export, listed
    character(len=CS), allocatable :: names(:)
    character(len=CL) :: line
    character(len=*), parameter :: subname = '(shr_bounds_init) '

    if (direction /= 'export' .and. direction /= 'import') then
       call shr_sys_abort(subname//"direction must be 'export' or 'import', not '"//trim(direction)//"'")
    end if
    export = (direction == 'export')
    bnd%label = trim(comp)//' '//trim(direction)

    bnd%nflds = shr_string_listGetNum(fieldlist)
    allocate(names(bnd%nflds))
    do k = 1, bnd%nflds
       call shr_string_listGetName(fieldlist, k, names(k))
    end do

    allocate(bnd%kfld(shr_field_limits_nflds), bnd%kwhere(shr_field_limits_nflds))
    allocate(bnd%vmin(shr_field_limits_nflds), bnd%vmax(shr_field_limits_nflds))
    allocate(bnd%nmin(shr_field_limits_nflds), bnd%nmax(shr_field_limits_nflds))
    allocate(bnd%nout(shr_field_limits_nflds))

    ! entries whose field (and where_positive field) is in the list; an
    ! export checks only the entries the component itself sends
    bnd%kfld = 0
    bnd%kwhere = 0
    do i = 1, shr_field_limits_nflds
       if (export .and. trim(shr_field_limits(i)%component) /= trim(comp)) cycle
       bnd%kfld(i) = list_index(names, shr_field_limits(i)%name)
       if (len_trim(shr_field_limits(i)%where_positive) > 0) then
          bnd%kwhere(i) = list_index(names, shr_field_limits(i)%where_positive)
          if (bnd%kwhere(i) == 0) bnd%kfld(i) = 0
       end if
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

    ! list the fields that have no limits
    call mpi_comm_rank(mpicom, rank, ierr)
    if (rank == 0) then
       nnone = 0
       line = ''
       do k = 1, bnd%nflds
          listed = .false.
          do i = 1, shr_field_limits_nflds
             if (bnd%kfld(i) == k) listed = .true.
          end do
          if (listed) cycle
          nnone = nnone + 1
          if (len_trim(line) + len_trim(names(k)) + 1 > 100) then
             write(logunit,'(a,a,a,a)') 'bounds_check: ', trim(bnd%label), ' no limits:', trim(line)
             line = ''
          end if
          line = trim(line)//' '//trim(names(k))
       end do
       if (len_trim(line) > 0) &
            write(logunit,'(a,a,a,a)') 'bounds_check: ', trim(bnd%label), ' no limits:', trim(line)
       write(logunit,'(a,a,i4,a,i4,a)') 'bounds_check: ', trim(bnd%label), &
            bnd%nflds - nnone, ' of ', bnd%nflds, ' fields have limits'
       call shr_sys_flush(logunit)
    end if

    deallocate(names)
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

    do i = 1, shr_field_limits_nflds
       k = bnd%kfld(i)
       if (k == 0) cycle
       kw = bnd%kwhere(i)
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

          if (bnd%nmin(i) == 0 .or. v < bnd%vmin(i)) then
             bnd%vmin(i) = v
             bnd%nmin(i) = n
          end if
          if (bnd%nmax(i) == 0 .or. v > bnd%vmax(i)) then
             bnd%vmax(i) = v
             bnd%nmax(i) = n
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
             bnd%nout(i) = bnd%nout(i) + 1.0_r8
             if (abort_on_violation) then
                write(msg,'(a,a,1x,a,a,es14.6,a,i9,a,f7.2,a,f7.2,a,a)') subname, trim(bnd%label), &
                     trim(shr_field_limits(i)%name), ' = ', v, ' at gid ', nint(bnd%gid(n)), &
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

    integer, parameter :: nf = shr_field_limits_nflds
    integer  :: i, n, rank, ierr
    real(r8) :: lpresent(nf), gpresent(nf), lnout(nf), gnout(nf)
    real(r8) :: lmin(2,nf), gmin(2,nf), lmax(2,nf), gmax(2,nf)
    real(r8) :: lloc(6,nf), gloc(6,nf)

    if (.not. bnd%initialized) return

    call mpi_comm_rank(mpicom, rank, ierr)

    lpresent = 0.0_r8
    lnout    = 0.0_r8
    lmin(1,:) =  huge(1.0_r8)
    lmax(1,:) = -huge(1.0_r8)
    lmin(2,:) = real(rank, r8)
    lmax(2,:) = real(rank, r8)
    do i = 1, nf
       if (bnd%kfld(i) == 0) cycle
       lpresent(i) = 1.0_r8
       lnout(i) = bnd%nout(i)
       if (bnd%nmin(i) > 0) lmin(1,i) = bnd%vmin(i)
       if (bnd%nmax(i) > 0) lmax(1,i) = bnd%vmax(i)
    end do

    call mpi_allreduce(lpresent, gpresent, nf, MPI_REAL8, MPI_MAX, mpicom, ierr)
    if (all(gpresent == 0.0_r8)) return

    call mpi_allreduce(lnout, gnout, nf, MPI_REAL8, MPI_SUM, mpicom, ierr)
    call mpi_allreduce(lmin, gmin, nf, MPI_2DOUBLE_PRECISION, MPI_MINLOC, mpicom, ierr)
    call mpi_allreduce(lmax, gmax, nf, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, mpicom, ierr)

    ! the task holding each extreme contributes its location
    lloc = 0.0_r8
    do i = 1, nf
       n = bnd%nmin(i)
       if (n > 0 .and. nint(gmin(2,i)) == rank) &
            lloc(1:3,i) = [bnd%gid(n), bnd%lat(n), bnd%lon(n)]
       n = bnd%nmax(i)
       if (n > 0 .and. nint(gmax(2,i)) == rank) &
            lloc(4:6,i) = [bnd%gid(n), bnd%lat(n), bnd%lon(n)]
    end do
    call mpi_allreduce(lloc, gloc, 6*nf, MPI_REAL8, MPI_SUM, mpicom, ierr)

    if (rank == 0) then
       write(logunit,'(a,a,a,i10,i6)') 'bounds_check: ', trim(bnd%label), ' model date = ', ymd, tod
       do i = 1, nf
          if (gpresent(i) == 0.0_r8) cycle
          if (gmin(1,i) > gmax(1,i)) then
             write(logunit,F02) trim(bnd%label), shr_field_limits(i)%name, trim(limits_string(i))
          else
             write(logunit,F01) trim(bnd%label), shr_field_limits(i)%name, &
                  gmin(1,i), nint(gloc(1,i)), gloc(2,i), gloc(3,i), &
                  gmax(1,i), nint(gloc(4,i)), gloc(5,i), gloc(6,i), &
                  nint(gnout(i)), trim(limits_string(i))
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

  integer function list_index(names, name)

    ! Index of name in names, 0 if absent

    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: name

    integer :: k

    list_index = 0
    do k = 1, size(names)
       if (trim(names(k)) == trim(name)) then
          list_index = k
          return
       end if
    end do

  end function list_index

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
