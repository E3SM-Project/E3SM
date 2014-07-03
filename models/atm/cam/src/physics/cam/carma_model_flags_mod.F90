!! This module handles reading the namelist and provides access to some other flags
!! that control a specific CARMA model's behavior.
!!
!! By default the specific CARMA model does not have any unique namelist values. If
!! a CARMA model wishes to have its own namelist, then this file needs to be copied
!! from physics/cam to physics/model/<model_name> and the code needed to read in the
!! namelist values added there. This file will take the place of the one in
!! physics/cam. 
!!
!! It needs to be in its own file to resolve some circular dependencies.
!!
!! @author  Chuck Bardeen
!! @version Mar-2011
module carma_model_flags_mod

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc

  ! Flags for integration with CAM Microphysics
  public carma_model_readnl                   ! read the carma model namelist
  

  ! Namelist flags
  !
  ! Create a public definition of any new namelist variables that you wish to have,
  ! and default them to an inital value.
  logical, public                :: carma_flag        = .false.   ! If .true. then turn on CARMA microphysics in CAM
  real(r8), public               :: carma_vf_const    = 0.0_r8    ! If specified and non-zero, constant fall velocity for all particles [cm/s]
  character(len=256), public     :: carma_reftfile    = 'carma_reft.nc'  ! path to the file containing the reference temperature profile

contains


  !! Read the CARMA model runtime options from the namelist
  !!
  !! @author  Chuck Bardeen
  !! @version Mar-2011
  subroutine carma_model_readnl(nlfile)
  
    ! Read carma namelist group.
  
    use abortutils,      only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
  
    ! args
  
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  
    ! local vars
  
    integer :: unitn, ierr
  
    ! read namelist for CARMA
!    namelist /carma_model_nl/ &
!      carma_flag, &
!      carma_maxretries, &
!      carma_conmax, &
!      carma_reftfile
  
!    if (masterproc) then
!       unitn = getunit()
!       open( unitn, file=trim(nlfile), status='old' )
!       call find_group_name(unitn, 'carma_model_nl', status=ierr)
!       if (ierr == 0) then
!          read(unitn, carma_model_nl, iostat=ierr)
!          if (ierr /= 0) then
!             call endrun('carma_model_readnl: ERROR reading namelist')
!          end if
!       end if
!       close(unitn)
!       call freeunit(unitn)
!    end if
  
#ifdef SPMD
!    call mpibcast (carma_flag,            1 ,mpilog, 0,mpicom)
!    call mpibcast (carma_maxretries,      1 ,mpiint, 0,mpicom)
!    call mpibcast (carma_conmax,          1 ,mpir8,  0,mpicom)
!    call mpibcast (carma_reftfile, len(carma_reftfile), mpichar, 0, mpicom)
#endif
  
  end subroutine carma_model_readnl

end module carma_model_flags_mod
