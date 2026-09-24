module shr_bounds_mod

  !---------------------------------------------------------------------
  ! Bounds checks on fields sent from components to the coupler (c2x).
  !
  ! Limits come from the limits dictionary (share/field_limits.yaml, generated
  ! into shr_field_limits_mod). shr_bounds_check runs on the component pes right
  ! after the component runs, on the component's own c2x_cc, and keeps
  ! running statistics. shr_bounds_report reduces those statistics onto
  ! the coupler root, writes them to the coupler log and resets them.
  !---------------------------------------------------------------------

  use shr_kind_mod   , only: r8 => SHR_KIND_R8, CL => SHR_KIND_CL
  use shr_sys_mod    , only: shr_sys_abort, shr_sys_flush
  use shr_infnan_mod , only: shr_infnan_isnan
  use shr_const_mod  , only: shr_const_isspval
  use shr_field_limits_mod
  use mct_mod

  implicit none
#include <mpif.h>
  private

  public :: shr_bounds_type
  public :: shr_bounds_check
  public :: shr_bounds_report

  type shr_bounds_type
     logical               :: initialized = .false.
     ! per dictionary entry
     integer , allocatable :: kfld(:)    ! index in c2x, 0 if the entry is not checked
     integer , allocatable :: kwhere(:)  ! index of the where_positive field in c2x, 0 if none
     real(r8), allocatable :: vmin(:)    ! smallest value since the last report
     real(r8), allocatable :: vmax(:)    ! largest value since the last report
     integer , allocatable :: nmin(:)    ! local point of vmin, 0 if no point checked
     integer , allocatable :: nmax(:)    ! local point of vmax, 0 if no point checked
     real(r8), allocatable :: nout(:)    ! number of out-of-bounds values since the last report
     ! per local point
     logical , allocatable :: active(:)  ! point is inside the component's mask
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

  subroutine shr_bounds_check(bnd, ntype, c2x, dom, gsmap, mpicom, abort_on_violation)

    ! Check c2x against its limits, update the statistics, and abort on
    ! the first out-of-bounds value if abort_on_violation is true.
    ! NaNs and special values are skipped; NaNs are check_fields' job.

    type(shr_bounds_type), intent(inout) :: bnd
    character(len=*)     , intent(in)    :: ntype  ! component type, e.g. 'atm'
    type(mct_aVect)      , intent(in)    :: c2x
    type(mct_gGrid)      , intent(in)    :: dom
    type(mct_gsMap)      , intent(in)    :: gsmap
    integer              , intent(in)    :: mpicom ! component communicator
    logical              , intent(in)    :: abort_on_violation

    integer  :: i, k, kw, n
    real(r8) :: v
    logical  :: out
    character(len=CL) :: msg
    character(len=*), parameter :: subname = '(shr_bounds_check) '

    if (.not. bnd%initialized) call shr_bounds_init(bnd, ntype, c2x, dom, gsmap, mpicom)

    do i = 1, shr_field_limits_nflds
       k = bnd%kfld(i)
       if (k == 0) cycle
       kw = bnd%kwhere(i)
       do n = 1, size(bnd%active)
          if (.not. bnd%active(n)) cycle
          v = c2x%rAttr(k,n)
          if (shr_infnan_isnan(v)) cycle
          if (shr_const_isspval(v)) cycle
          if (kw > 0) then
             if (.not. (c2x%rAttr(kw,n) > 0.0_r8)) cycle
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
                write(msg,'(a,a,1x,a,a,es14.6,a,i9,a,f7.2,a,f7.2,a,a)') subname, ntype, &
                     trim(shr_field_limits(i)%name), ' = ', v, ' at gid ', nint(bnd%gid(n)), &
                     ' lat ', bnd%lat(n), ' lon ', bnd%lon(n), ' outside limits ', &
                     trim(limits_string(i))
                call shr_sys_abort(trim(msg))
             end if
          end if
       end do
    end do

  end subroutine shr_bounds_check

  !===============================================================================

  subroutine shr_bounds_report(bnd, ntype, mpicom, iamroot_cpl, logunit)

    ! Reduce the statistics of one component instance onto all tasks of
    ! mpicom, write them from the coupler root, and reset them.
    ! Must be called by every task of mpicom (the union of coupler and
    ! component pes); tasks outside the component contribute nothing.

    type(shr_bounds_type), intent(inout) :: bnd
    character(len=*)     , intent(in)    :: ntype
    integer              , intent(in)    :: mpicom      ! coupler + component communicator
    logical              , intent(in)    :: iamroot_cpl ! writes the report
    integer              , intent(in)    :: logunit

    integer, parameter :: nf = shr_field_limits_nflds
    integer  :: i, n, rank, ierr
    real(r8) :: lpresent(nf), gpresent(nf), lnout(nf), gnout(nf)
    real(r8) :: lmin(2,nf), gmin(2,nf), lmax(2,nf), gmax(2,nf)
    real(r8) :: lloc(6,nf), gloc(6,nf)

    call mpi_comm_rank(mpicom, rank, ierr)

    lpresent = 0.0_r8
    lnout    = 0.0_r8
    lmin(1,:) =  huge(1.0_r8)
    lmax(1,:) = -huge(1.0_r8)
    lmin(2,:) = real(rank, r8)
    lmax(2,:) = real(rank, r8)
    if (bnd%initialized) then
       do i = 1, nf
          if (bnd%kfld(i) == 0) cycle
          lpresent(i) = 1.0_r8
          lnout(i) = bnd%nout(i)
          if (bnd%nmin(i) > 0) lmin(1,i) = bnd%vmin(i)
          if (bnd%nmax(i) > 0) lmax(1,i) = bnd%vmax(i)
       end do
    end if

    call mpi_allreduce(lpresent, gpresent, nf, MPI_REAL8, MPI_MAX, mpicom, ierr)
    if (all(gpresent == 0.0_r8)) return

    call mpi_allreduce(lnout, gnout, nf, MPI_REAL8, MPI_SUM, mpicom, ierr)
    call mpi_allreduce(lmin, gmin, nf, MPI_2DOUBLE_PRECISION, MPI_MINLOC, mpicom, ierr)
    call mpi_allreduce(lmax, gmax, nf, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, mpicom, ierr)

    ! the task holding each extreme contributes its location
    lloc = 0.0_r8
    if (bnd%initialized) then
       do i = 1, nf
          n = bnd%nmin(i)
          if (n > 0 .and. nint(gmin(2,i)) == rank) &
               lloc(1:3,i) = [bnd%gid(n), bnd%lat(n), bnd%lon(n)]
          n = bnd%nmax(i)
          if (n > 0 .and. nint(gmax(2,i)) == rank) &
               lloc(4:6,i) = [bnd%gid(n), bnd%lat(n), bnd%lon(n)]
       end do
    end if
    call mpi_allreduce(lloc, gloc, 6*nf, MPI_REAL8, MPI_SUM, mpicom, ierr)

    if (iamroot_cpl) then
       do i = 1, nf
          if (gpresent(i) == 0.0_r8) cycle
          if (gmin(1,i) > gmax(1,i)) then
             write(logunit,F02) ntype, shr_field_limits(i)%name, trim(limits_string(i))
          else
             write(logunit,F01) ntype, shr_field_limits(i)%name, &
                  gmin(1,i), nint(gloc(1,i)), gloc(2,i), gloc(3,i), &
                  gmax(1,i), nint(gloc(4,i)), gloc(5,i), gloc(6,i), &
                  nint(gnout(i)), trim(limits_string(i))
          end if
       end do
       call shr_sys_flush(logunit)
    end if

    if (bnd%initialized) call shr_bounds_reset(bnd)

  end subroutine shr_bounds_report

  !===============================================================================

  subroutine shr_bounds_init(bnd, ntype, c2x, dom, gsmap, mpicom)

    type(shr_bounds_type), intent(inout) :: bnd
    character(len=*)     , intent(in)    :: ntype
    type(mct_aVect)      , intent(in)    :: c2x
    type(mct_gGrid)      , intent(in)    :: dom
    type(mct_gsMap)      , intent(in)    :: gsmap
    integer              , intent(in)    :: mpicom

    integer :: i, n, lsize, rank, ierr, kmask, klat, klon
    integer, pointer :: gpts(:)

    allocate(bnd%kfld(shr_field_limits_nflds), bnd%kwhere(shr_field_limits_nflds))
    allocate(bnd%vmin(shr_field_limits_nflds), bnd%vmax(shr_field_limits_nflds))
    allocate(bnd%nmin(shr_field_limits_nflds), bnd%nmax(shr_field_limits_nflds))
    allocate(bnd%nout(shr_field_limits_nflds))

    ! entries for this component whose fields (and where_positive fields) are sent
    bnd%kfld = 0
    bnd%kwhere = 0
    do i = 1, shr_field_limits_nflds
       if (trim(shr_field_limits(i)%component) /= trim(ntype)) cycle
       bnd%kfld(i) = mct_aVect_indexRA(c2x, trim(shr_field_limits(i)%name), perrWith='quiet')
       if (len_trim(shr_field_limits(i)%where_positive) > 0) then
          bnd%kwhere(i) = mct_aVect_indexRA(c2x, trim(shr_field_limits(i)%where_positive), perrWith='quiet')
          if (bnd%kwhere(i) == 0) bnd%kfld(i) = 0
       end if
    end do

    lsize = mct_aVect_lsize(c2x)
    allocate(bnd%active(lsize), bnd%gid(lsize), bnd%lat(lsize), bnd%lon(lsize))

    kmask = mct_aVect_indexRA(dom%data, 'mask')
    klat  = mct_aVect_indexRA(dom%data, 'lat')
    klon  = mct_aVect_indexRA(dom%data, 'lon')
    do n = 1, lsize
       bnd%active(n) = dom%data%rAttr(kmask,n) /= 0.0_r8
       bnd%lat(n) = dom%data%rAttr(klat,n)
       bnd%lon(n) = dom%data%rAttr(klon,n)
    end do

    call mpi_comm_rank(mpicom, rank, ierr)
    call mct_gsMap_orderedPoints(gsmap, rank, gpts)
    bnd%gid(:) = real(gpts(1:lsize), r8)
    deallocate(gpts)

    call shr_bounds_reset(bnd)
    bnd%initialized = .true.

  end subroutine shr_bounds_init

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
