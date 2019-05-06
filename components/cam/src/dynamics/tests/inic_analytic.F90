module inic_analytic

  !-----------------------------------------------------------------------
  !
  ! Purpose: Set analytic initial conditions based on input coordinates
  !
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use shr_sys_mod,         only: shr_sys_flush
  use inic_analytic_utils, only: analytic_ic_active, analytic_ic_type

  implicit none
  private

  public :: analytic_ic_active ! forwarded from init_analytic_utils
  public :: analytic_ic_set_ic ! Set analytic initial conditions

  interface analytic_ic_set_ic
    module procedure dyn_set_inic_cblock
  end interface analytic_ic_set_ic

  ! Private module variables
  integer                              :: call_num = 0

  ! Private interface
#ifdef ANALYTIC_IC
  interface get_input_shape
    module procedure get_input_shape_2d
    module procedure get_input_shape_3d
  end interface get_input_shape
#endif

!==============================================================================
CONTAINS
!==============================================================================

  subroutine dyn_set_inic_col(vcoord, latvals, lonvals, glob_ind, U, V, T,    &
       PS, PHIS, Q, m_cnst, mask, verbose)
    use cam_initfiles,        only: pertlim
#ifdef ANALYTIC_IC
    use ic_held_suarez,       only: hs94_set_ic
    use ic_baroclinic,        only: bc_wav_set_ic
#endif
    use spmd_utils,           only: masterproc
    !-----------------------------------------------------------------------
    !
    ! Purpose: Set analytic initial values for dynamics state variables
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    integer           , intent(in)    :: vcoord      ! See dyn_tests_utils
    real(r8),           intent(in)    :: latvals(:)  ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:)  ! lon in degrees (ncol)
    integer,            intent(in)    :: glob_ind(:) ! global column index
    real(r8), optional, intent(inout) :: U(:,:)      ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)      ! meridional velocity
    real(r8), optional, intent(inout) :: T(:,:)      ! temperature
    real(r8), optional, intent(inout) :: PS(:)       ! surface pressure
    real(r8), optional, intent(inout) :: PHIS(:)     ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:)    ! tracer (ncol, lev, m)
    integer,  optional, intent(in)    :: m_cnst(:)   ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)     ! Only init where .true.
    logical,  optional, intent(in)    :: verbose     ! For internal use

    ! Local variables
    logical                           :: verbose_use
    logical, allocatable              :: mask_use(:)
    real(r8)                          :: pertval
    integer, allocatable              :: rndm_seed(:)
    integer                           :: rndm_seed_sz
    integer                           :: i, k
    integer                           :: ncol, nlev
    character(len=*), parameter       :: subname = 'DYN_SET_INIC_COL'

#ifdef ANALYTIC_IC
    allocate(mask_use(size(latvals)))
    if (present(mask)) then
      if (size(mask_use) /= size(mask)) then
        call endrun('cnst_init_default: input, mask, is wrong size')
      end if
      mask_use = mask
    else
      mask_use = .true.
    end if

    if (present(verbose)) then
      verbose_use = verbose
    else
      verbose_use = .true.
    end if

    ! Basic size sanity checks
    if (size(latvals) /= size(lonvals)) then
      call endrun(subname//'latvals and lonvals must have same size')
    end if
    if (present(U)) then
      if (size(U) > 0) then
        call check_array_size(U(:,1), 'U', latvals, subname)
      else
        return
      end if
    end if
    if (present(V)) then
      if (size(V) > 0) then
        call check_array_size(V(:,1), 'V', latvals, subname)
      else
        return
      end if
    end if
    if (present(T)) then
      if (size(T) > 0) then
        call check_array_size(T(:,1), 'T', latvals, subname)
      else
        return
      end if
    end if
    if (present(PS)) then
      if (size(PS) > 0) then
        call check_array_size(PS, 'PS', latvals, subname)
      else
        return
      end if
    end if
    if (present(PHIS)) then
      if (size(PHIS) > 0) then
        call check_array_size(PHIS, 'PHIS', latvals, subname)
      else
        return
      end if
    end if
    ! Some special checks on the tracer argument
    if (present(Q)) then
      if (.not. present(m_cnst)) then
        call endrun(subname//'m_cnst is required if Q is present')
      end if
      if (size(Q, 3) /= size(m_cnst, 1)) then
        call endrun(subname//': size of m_cnst must match last dimension of Q')
      end if
      if (size(Q) > 0) then
        call check_array_size(Q(:,1,1), 'Q', latvals, subname)
      else
        return
      end if
    end if

    select case(trim(analytic_ic_type))
    case('held_suarez_1994')
      call hs94_set_ic(latvals, lonvals, U=U, V=V, T=T, PS=PS, PHIS=PHIS,     &
           Q=Q, m_cnst=m_cnst, mask=mask_use, verbose=verbose_use)

    case('baroclinic_wave', 'dry_baroclinic_wave')


      call bc_wav_set_ic(vcoord, latvals, lonvals, U=U, V=V, T=T, PS=PS,      &
           PHIS=PHIS, Q=Q, m_cnst=m_cnst, mask=mask_use, verbose=verbose_use)

    case default
      call endrun(subname//': Unknown analytic_ic_type, "'//trim(analytic_ic_type)//'"')
    end select

    ! Maybe peturb T initial conditions
    if (present(T) .and. (pertlim /= 0.0_r8)) then

      ! Add random perturbation to temperature if required
      if(masterproc .and. verbose_use) then
        write(iulog,*) trim(subname), ': Adding random perturbation bounded by +/-', &
             pertlim,' to initial temperature field'
      end if
      call random_seed(size=rndm_seed_sz)
      allocate(rndm_seed(rndm_seed_sz))

      ncol = size(T, 1)
      nlev = size(T, 2)
      do i = 1, ncol
        if (mask_use(i)) then
          ! seed random_number generator based on global column index
          rndm_seed(:) = glob_ind(i)
          call random_seed(put=rndm_seed)
          do k = 1, nlev
            call random_number(pertval)
            pertval = 2.0_r8 * pertlim * (0.5_r8 - pertval)
            T(i,k) = T(i,k) * (1.0_r8 + pertval)
          end do
        end if
      end do

      deallocate(rndm_seed)
    end if

    ! To get different random seeds each time
    call_num = call_num + 1
#else
    call endrun(subname//': analytic initial conditions are not enabled')
#endif

  end subroutine dyn_set_inic_col

  subroutine dyn_set_inic_cblock(vcoord,latvals, lonvals, glob_ind, U, V, T,  &
       PS, PHIS, Q, m_cnst, mask)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Set analytic initial values for dynamics state variables
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    integer,            intent(in)    :: vcoord      ! See dyn_tests_utils
    real(r8),           intent(in)    :: latvals(:)  ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:)  ! lon in degrees (ncol)
    integer,            intent(in)    :: glob_ind(:) ! global column index
    real(r8), optional, intent(inout) :: U(:,:,:)    ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:,:)    ! meridional velocity
    real(r8), optional, intent(inout) :: T(:,:,:)    ! temperature
    real(r8), optional, intent(inout) :: PS(:,:)     ! surface pressure
    real(r8), optional, intent(inout) :: PHIS(:,:)   ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:,:)  ! tracer (ncol,lev,blk,m)
    integer,  optional, intent(in)    :: m_cnst(:)   ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)     ! Only init where .true.

    ! Local variables
    real(r8), allocatable          :: lat_use(:)
    integer                        :: i, bbeg, bend
    integer                        :: size1, size2, size3
    integer                        :: nblks, blksize
    logical                        :: verbose
    character(len=4)               :: mname
    character(len=*), parameter    :: subname = 'DYN_SET_INIC_CBLOCK'

#ifdef ANALYTIC_IC
    verbose = .true. ! So subroutines can report setting variables
    ! Figure out what sort of blocks we have, all variables should be the same
    size1 = -1
    mname = ''
    if (present(U)) then
      call get_input_shape(U, 'U', mname, size1, size2, size3, subname)
    end if
    if(present(V)) then
      call get_input_shape(V, 'V', mname, size1, size2, size3, subname)
    end if
    if(present(T)) then
      call get_input_shape(T, 'T', mname, size1, size2, size3, subname)
    end if
    if(present(Q)) then
      call get_input_shape(Q(:,:,:,1), 'Q', mname, size1, size2, size3, subname)
    end if
    ! Need to do all 3-D variables before any 2-D variables
    if(present(PS)) then
      call get_input_shape(PS, 'PS', mname, size1, size2, size3, subname)
    end if
    if(present(PHIS)) then
      call get_input_shape(PHIS, 'PHIS', mname, size1, size2, size3, subname)
    end if
    if (size1 < 0) then
      call endrun(subname//': No state variables to initialize')
    end if
    if ((size(latvals) == size1*size3) .and. (size(lonvals) == size1*size3)) then
      ! Case: unstructured with blocks in 3rd dim
      if (size(glob_ind) /= size(latvals)) then
        call endrun(subname//': there must be a global index for every column')
      end if
      nblks = size3
      blksize = size1
      bend = 0
      do i = 1, nblks
        bbeg = bend + 1
        bend = bbeg + blksize - 1
        if (present(mask)) then
          if (size(mask) /= size(latvals)) then
            call endrun(subname//': incorrect mask size')
          end if
          if (present(U)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), U=U(:,:,i), mask=mask(bbeg:bend), verbose=verbose)
          end if
          if (present(V)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), V=V(:,:,i), mask=mask(bbeg:bend), verbose=verbose)
          end if
          if (present(T)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), T=T(:,:,i), mask=mask(bbeg:bend), verbose=verbose)
          end if
          if (present(PS)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), PS=PS(:,i), mask=mask(bbeg:bend), verbose=verbose)
          end if
          if (present(PHIS)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), PHIS=PHIS(:,i), mask=mask(bbeg:bend), verbose=verbose)
          end if
          if (present(Q)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), Q=Q(:,:,i,:), m_cnst=m_cnst,            &
                 mask=mask(bbeg:bend), verbose=verbose)
          end if
        else
          if (present(U)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), U=U(:,:,i), verbose=verbose)
          end if
          if (present(V)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), V=V(:,:,i), verbose=verbose)
          end if
          if (present(T)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), T=T(:,:,i), verbose=verbose)
          end if
          if (present(PS)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), PS=PS(:,i), verbose=verbose)
          end if
          if (present(PHIS)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), PHIS=PHIS(:,i), verbose=verbose)
          end if
          if (present(Q)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), Q=Q(:,:,i,:), m_cnst=m_cnst,            &
                 verbose=verbose)
          end if
        end if
        verbose = .false.
      end do
    else if ((size(latvals) == size1*size2) .and. (size(lonvals) == size1*size2)) then
      ! Case: unstructured with blocks in 2nd dim
      if (size(glob_ind) /= size(latvals)) then
        call endrun(subname//': there must be a global index for every column')
      end if
      nblks = size2
      blksize = size1
      bend = 0
      do i = 1, nblks
        bbeg = bend + 1
        bend = bbeg + blksize - 1
        if (present(mask)) then
          if (size(mask) /= size(latvals)) then
            call endrun(subname//': incorrect mask size')
          end if
          if (present(U)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), U=U(:,i,:), mask=mask(bbeg:bend), verbose=verbose)
          end if
          if (present(V)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), V=V(:,i,:), mask=mask(bbeg:bend), verbose=verbose)
          end if
          if (present(T)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), T=T(:,i,:), mask=mask(bbeg:bend), verbose=verbose)
          end if
          if (present(PS)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), PS=PS(:,i), mask=mask(bbeg:bend), verbose=verbose)
          end if
          if (present(PHIS)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), PHIS=PHIS(:,i), mask=mask(bbeg:bend), verbose=verbose)
          end if
          if (present(Q)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), Q=Q(:,i,:,:), m_cnst=m_cnst,            &
                 mask=mask(bbeg:bend), verbose=verbose)
          end if
        else
          if (present(U)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), U=U(:,i,:), verbose=verbose)
          end if
          if (present(V)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), V=V(:,i,:), verbose=verbose)
          end if
          if (present(T)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), T=T(:,i,:), verbose=verbose)
          end if
          if (present(PS)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), PS=PS(:,i), verbose=verbose)
          end if
          if (present(PHIS)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), PHIS=PHIS(:,i), verbose=verbose)
          end if
          if (present(Q)) then
            call dyn_set_inic_col(vcoord,latvals(bbeg:bend), lonvals(bbeg:bend), &
                 glob_ind(bbeg:bend), Q=Q(:,i,:,:), m_cnst=m_cnst,            &
                 verbose=verbose)
          end if
        end if
        verbose = .false.
      end do
    else if ((size(latvals) == size2) .and. (size(lonvals) == size1)) then
      ! Case: lon,lat,lev
      if (size(glob_ind) /= (size2 * size1)) then
        call endrun(subname//': there must be a global index for every column')
      end if
      nblks = size2
      allocate(lat_use(size(lonvals)))
      if (present(mask)) then
        call endrun(subname//': mask not supported for lon/lat')
      else
        bend = 0
        do i = 1, nblks
          bbeg = bend + 1
          bend = bbeg + size1 - 1
          lat_use = latvals(i)
          if (present(U)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 U=U(:,i,:), verbose=verbose)
          end if
          if (present(V)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 V=V(:,i,:), verbose=verbose)
          end if
          if (present(T)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 T=T(:,i,:), verbose=verbose)
          end if
          if (present(PS)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 PS=PS(:,i), verbose=verbose)
          end if
          if (present(PHIS)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 PHIS=PHIS(:,i), verbose=verbose)
          end if
          if (present(Q)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 Q=Q(:,i,:,:), m_cnst=m_cnst, verbose=verbose)
          end if
          verbose = .false.
        end do
      end if
      deallocate(lat_use)
    else if ((size(latvals) == size3) .and. (size(lonvals) == size1)) then
      if (size(glob_ind) /= (size3 * size1)) then
        call endrun(subname//': there must be a global index for every column')
      end if
      ! Case: lon,lev,lat
      nblks = size3
      allocate(lat_use(size(lonvals)))
      if (present(mask)) then
        call endrun(subname//': mask not supported for lon/lat')
      else
        bend = 0
        do i = 1, nblks
          bbeg = bend + 1
          bend = bbeg + size1 - 1
          lat_use = latvals(i)
          if (present(U)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 U=U(:,:,i), verbose=verbose)
          end if
          if (present(V)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 V=V(:,:,i), verbose=verbose)
          end if
          if (present(T)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 T=T(:,:,i), verbose=verbose)
          end if
          if (present(PS)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 PS=PS(:,i), verbose=verbose)
          end if
          if (present(PHIS)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 PHIS=PHIS(:,i), verbose=verbose)
          end if
          if (present(Q)) then
            call dyn_set_inic_col(vcoord,lat_use, lonvals, glob_ind(bbeg:bend), &
                 Q=Q(:,:,i,:), m_cnst=m_cnst, verbose=verbose)
          end if
          verbose = .false.
        end do
      end if
      deallocate(lat_use)
    else
      call endrun(subname//': Unknown state variable layout')
    end if
#else
    call endrun(subname//': analytic initial conditions are not enabled')
#endif
  end subroutine dyn_set_inic_cblock

#ifdef ANALYTIC_IC
  subroutine get_input_shape_2d(array, aname, sname, size1, size2, size3, es)
    real(r8),         intent(in)    :: array(:,:)
    character(len=*), intent(in)    :: aname
    character(len=*), intent(inout) :: sname
    integer,          intent(inout) :: size1
    integer,          intent(inout) :: size2
    integer,          intent(inout) :: size3
    character(len=*), intent(in)    :: es

    if ((size1 < 0) .and. (size(array) == 0)) then
      ! The shape has not yet been set, set it to zero
      size1 = 0
      size2 = 0
      size3 = 0
      sname = trim(aname)
    else if (size1 < 0) then
      ! The shape has not yet been set, set it
      size1 = size(array, 1)
      size2 = size(array, 2)
      size3 = 1
      sname = trim(aname)
    else if ((size1 * size2 * size3) > 0) then
      ! For 2-D variables, the second dimension is always the block size
      ! However, since the shape may have been set by a 3-D variable, we
      ! need to pass either possibility
      if (  (size1 /= size(array, 1)) .or.                                    &
           ((size2 /= size(array, 2)) .and. (size3 /= size(array, 2)))) then
        call endrun(trim(es)//': shape of '//trim(aname)//' does not match shape of '//trim(sname))
      end if
    ! No else, we cannot compare to zero size master array
    end if

  end subroutine get_input_shape_2d

  subroutine get_input_shape_3d(array, aname, sname, size1, size2, size3, es)
    real(r8),         intent(in)    :: array(:,:,:)
    character(len=*), intent(in)    :: aname
    character(len=*), intent(inout) :: sname
    integer,          intent(inout) :: size1
    integer,          intent(inout) :: size2
    integer,          intent(inout) :: size3
    character(len=*), intent(in)    :: es

    if ((size1 < 0) .and. (size(array) == 0)) then
      ! The shape has not yet been set, set it to zero
      size1 = 0
      size2 = 0
      size3 = 0
      sname = trim(aname)
    else if (size1 < 0) then
      ! The shape has not yet been set, set it
      size1 = size(array, 1)
      size2 = size(array, 2)
      size3 = size(array, 3)
      sname = trim(aname)
    else if ((size1 * size2 * size3) > 0) then
      ! We have a shape, make sure array matches it
      if ((size1 /= size(array, 1)) .or. (size2 /= size(array, 2)) .or. (size3 /= size(array, 3))) then
        call endrun(trim(es)//': shape of '//trim(aname)//' does not match shape of '//trim(sname))
      end if
    end if
    ! No else, we cannot compare to zero size master array
  end subroutine get_input_shape_3d

  subroutine check_array_size(array, aname, check, subname)
    real(r8),         intent(in)    :: array(:)
    character(len=*), intent(in)    :: aname
    real(r8),         intent(in)    :: check(:)
    character(len=*), intent(in)    :: subname

    if (size(array, 1) /= size(check, 1)) then
      call endrun(trim(subname)//': '//trim(aname)//' has the wrong first dimension')
    end if

  end subroutine check_array_size
#endif

end module inic_analytic
