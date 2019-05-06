module const_init

! Initialize constituents to default values

use shr_kind_mod,     only: r8 => shr_kind_r8, max_chars=>shr_kind_cl
use spmd_utils,       only: masterproc
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog

implicit none
private
save

public :: cnst_init_default

interface cnst_init_default
  module procedure cnst_init_default_col
  module procedure cnst_init_default_cblock
end interface cnst_init_default

!==============================================================================
CONTAINS
!==============================================================================

  subroutine cnst_init_default_col(m_cnst, latvals, lonvals, q, mask,         &
       verbose, notfound)
    use constituents,  only: cnst_name
    use aoa_tracers,   only: aoa_tracers_implements_cnst,   aoa_tracers_init_cnst
    use carma_intr,    only: carma_implements_cnst,         carma_init_cnst
    use chemistry,     only: chem_implements_cnst,          chem_init_cnst
    use clubb_intr,    only: clubb_implements_cnst,         clubb_init_cnst
    use co2_cycle,     only: co2_implements_cnst,           co2_init_cnst
    use microp_driver, only: microp_driver_implements_cnst, microp_driver_init_cnst
    !!!use rk_stratiform, only: rk_stratiform_implements_cnst, rk_stratiform_init_cnst
    use stratiform, only: stratiform_implements_cnst, stratiform_init_cnst
    use tracers,       only: tracers_implements_cnst,       tracers_init_cnst
    use unicon_cam,    only: unicon_implements_cnst,        unicon_init_cnst

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize named tracer mixing ratio field
    !  This subroutine should be called ONLY at the beginning of an initial run
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    integer,           intent(in)  :: m_cnst     ! Constant index
    real(r8),          intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8),          intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    real(r8),          intent(out) :: q(:,:)     ! mixing ratio (ncol, plev)
    logical, optional, intent(in)  :: mask(:)    ! Only initialize where .true.
    logical, optional, intent(in)  :: verbose    ! For internal use
    logical, optional, intent(in)  :: notfound   ! Turn off initial dataset warn

    ! Local variables
    logical, allocatable           :: mask_use(:)
    character(len=max_chars)       :: name
    logical                        :: verbose_use
    logical                        :: notfound_use

    name = cnst_name(m_cnst)

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

    if (present(notfound)) then
      notfound_use = notfound
    else
      notfound_use = .true.
    end if

    q = 0.0_r8 ! Make sure we start fresh (insurance)

    if(masterproc .and. verbose_use .and. notfound_use) then
      write(iulog, *) 'Field ',trim(trim(name)),' not found on initial dataset'
    end if

    if (aoa_tracers_implements_cnst(trim(name))) then
      call aoa_tracers_init_cnst(trim(name), latvals, lonvals, mask_use, q)
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          ', trim(name), ' initialized by "aoa_tracers_init_cnst"'
      end if
    else if (carma_implements_cnst(trim(name))) then
      call carma_init_cnst(trim(name), latvals, lonvals, mask_use, q)
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          ', trim(name), ' initialized by "carma_init_cnst"'
      end if
    else if (chem_implements_cnst(trim(name))) then
      call chem_init_cnst(trim(name), latvals, lonvals, mask_use, q)
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          ', trim(name), ' initialized by "chem_init_cnst"'
      end if
    else if (clubb_implements_cnst(trim(name))) then
      call clubb_init_cnst(trim(name), latvals, lonvals, mask_use, q)
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          ', trim(name), ' initialized by "clubb_init_cnst"'
      end if
    else if (co2_implements_cnst(trim(name))) then
      call co2_init_cnst(trim(name), latvals, lonvals, mask_use, q)
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          ', trim(name), ' initialized by "co2_init_cnst"'
      end if
    else if (microp_driver_implements_cnst(trim(name))) then
      call microp_driver_init_cnst(trim(name), latvals, lonvals, mask_use, q)
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          ', trim(name), ' initialized by "microp_driver_init_cnst"'
      end if
    !!!else if (rk_stratiform_implements_cnst(trim(name))) then
    !!!  call rk_stratiform_init_cnst(trim(name), latvals, lonvals, mask_use, q)
    else if (stratiform_implements_cnst(trim(name))) then
      call stratiform_init_cnst(trim(name), latvals, lonvals, mask_use, q)
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          ', trim(name), ' initialized by "stratiform_init_cnst"'
      end if
    else if (tracers_implements_cnst(trim(name))) then
      call tracers_init_cnst(trim(name), latvals, lonvals, mask_use, q)
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          ', trim(name), ' initialized by "tracers_init_cnst"'
      end if
    else if (unicon_implements_cnst(trim(name))) then
      call unicon_init_cnst(trim(name), latvals, lonvals, mask_use, q)
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          ', trim(name), ' initialized by "unicon_init_cnst"'
      end if
    else
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          ', trim(name), ' set to minimum value'
      end if
      ! Q already set to zero
    end if

  end subroutine cnst_init_default_col

  subroutine cnst_init_default_cblock(m_cnst, latvals, lonvals, q, mask)

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize named tracer mixing ratio field
    !  This subroutine should be called ONLY at the beginning of an initial run
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    integer,           intent(in)  :: m_cnst     ! Constant index
    real(r8),          intent(in)  :: latvals(:) ! lat in degrees (ncol*blk)
    real(r8),          intent(in)  :: lonvals(:) ! lon in degrees (ncol*blk)
    real(r8),          intent(out) :: q(:,:,:)   ! mix ratio (ncol, plev, blk)
    logical, optional, intent(in)  :: mask(:)    ! Only initialize where .true.

    ! Local variables
    real(r8), allocatable         :: latblk(:)
    integer                       :: i, bbeg, bend
    integer                       :: size1, size2, size3
    integer                       :: nblks, blksize
    logical                       :: verbose

    verbose = .true.
    size1 = size(q, 1)
    size2 = size(q, 2)
    size3 = size(q, 3)
    if ((size(latvals) == size1*size3) .and. (size(lonvals) == size1*size3)) then
      ! Case: unstructured with blocks in 3rd dim
      nblks = size3
      blksize = size1
      bend = 0
      do i = 1, nblks
        bbeg = bend + 1
        bend = bbeg + blksize - 1
        if (present(mask)) then
          if (size(mask) /= size(latvals)) then
            call endrun('cnst_init_default_cblock: incorrect mask size')
          end if
          call cnst_init_default(m_cnst, latvals(bbeg:bend), lonvals(bbeg:bend), q(:,:,i), mask=mask(bbeg:bend), verbose=verbose)
        else
          call cnst_init_default(m_cnst, latvals(bbeg:bend), lonvals(bbeg:bend), q(:,:,i), verbose=verbose)
        end if
        verbose = .false.
      end do
    else if ((size(latvals) == size2) .and. (size(lonvals) == size1)) then
      ! Case: lon,lat,lev
      if (present(mask)) then
        call endrun('cnst_init_default_cblock: mask not supported for lon/lat')
      else
        nblks = size2
        allocate(latblk(size1))
        do i = 1, nblks
          latblk(:) = latvals(i)
          call cnst_init_default(m_cnst, latblk, lonvals, q(:,i,:), verbose=verbose)
          verbose = .false.
        end do
        deallocate(latblk)
      end if
    else if ((size(latvals) == size3) .and. (size(lonvals) == size1)) then
      ! Case: lon,lev,lat
      if (present(mask)) then
        call endrun('cnst_init_default_cblock: mask not supported for lon/lat')
      else
        nblks = size3
        allocate(latblk(size1))
        do i = 1, nblks
          latblk(:) = latvals(i)
          call cnst_init_default(m_cnst, latblk, lonvals, q(:,:,i), verbose=verbose)
          verbose = .false.
        end do
        deallocate(latblk)
      end if
    else
      call endrun('cnst_init_default_cblock: Unknown q layout')
    end if

  end subroutine cnst_init_default_cblock

end module const_init
