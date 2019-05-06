module ic_held_suarez

  !-----------------------------------------------------------------------
  !
  ! Purpose: Set Held-Suarez initial conditions based on input coordinates
  !
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc
  use shr_sys_mod,         only: shr_sys_flush

  implicit none
  private

  ! Public interface
  public :: hs94_set_ic

!==============================================================================
CONTAINS
!==============================================================================

  subroutine hs94_set_ic(latvals, lonvals, U, V, T, PS, PHIS,           &
       Q, m_cnst, mask, verbose)
    use const_init,    only: cnst_init_default
    use constituents,  only: cnst_name

    !-----------------------------------------------------------------------
    !
    ! Purpose: Set Held-Suarez initial values for dynamics state variables
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
    real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
    real(r8), optional, intent(inout) :: T(:,:)     ! temperature
    real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
    real(r8), optional, intent(inout) :: PHIS(:)    ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
    integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)    ! Only init where .true.
    logical,  optional, intent(in)    :: verbose    ! For internal use

    ! Local variables
    logical, allocatable              :: mask_use(:)
    logical                           :: verbose_use
    integer                           :: i, k, m
    integer                           :: ncol
    integer                           :: nlev
    integer                           :: ncnst
    character(len=*), parameter       :: subname = 'HS94_SET_IC'

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

    ncol = size(latvals, 1)
    nlev = -1
    if (present(U)) then
      nlev = size(U, 2)
      do k = 1, nlev
        where(mask_use)
          U(:,k) = 0.0_r8
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          U initialized by "',subname,'"'
      end if
    end if

    if (present(V)) then
      nlev = size(V, 2)
      do k = 1, nlev
        where(mask_use)
          V(:,k) = 0.0_r8
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          V initialized by "',subname,'"'
      end if
    end if

    if (present(T)) then
      nlev = size(T, 2)
      do k = 1, nlev
        where(mask_use)
          T(:,k) = 250.0_r8
        end where
      end do
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          T initialized by "',subname,'"'
      end if
    end if

    if (present(PS)) then
      where(mask_use)
        PS = 100000.0_r8
      end where
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PS initialized by "',subname,'"'
      end if
    end if

    if (present(PHIS)) then
      where(mask_use)
        PHIS = 0.0_r8
      end where
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PHIS initialized by "',subname,'"'
      end if
    end if

    if (present(Q)) then
      nlev = size(Q, 2)
      ncnst = size(m_cnst, 1)
      do m = 1, ncnst
        if (m_cnst(m) == 1) then
          ! No water vapor in Held-Suarez
          do k = 1, nlev
            where(mask_use)
              Q(:,k,m_cnst(m)) = 0.0_r8
            end where
          end do
          if(masterproc .and. verbose_use) then
            write(iulog,*) '          ', trim(cnst_name(m_cnst(m))), ' initialized by "',subname,'"'
          end if
        else
          call cnst_init_default(m_cnst(m), latvals, lonvals, Q(:,:,m_cnst(m)),&
               mask=mask_use, verbose=verbose_use, notfound=.false.)
        end if
      end do
    end if

    deallocate(mask_use)

  end subroutine hs94_set_ic

end module ic_held_suarez
