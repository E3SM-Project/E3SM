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
  character(len=256), public     ::  carma_soilerosion_file  = 'soil_erosion_factor_1x1_c120907.nc'   ! name of the dust erosion factor file

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
    namelist /carma_model_nl/ &
      carma_soilerosion_file
  
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'carma_model_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, carma_model_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('carma_model_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if
  
#ifdef SPMD
    call mpibcast(carma_soilerosion_file,      len(carma_soilerosion_file),       mpichar, 0, mpicom)
#endif
  
  end subroutine carma_model_readnl

end module carma_model_flags_mod
