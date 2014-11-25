subroutine init_hydrology( NLFilename )
!
!DESCRIPTION
! Initialize implementation methods for different hydrology sub-modules 
! This is created for unit-based sensitivity tests
! created by Jinyun Tang, Mar 22, 2014.

  ! !USES:
  use spmdMod       , only : masterproc, mpicom
  use fileutils     , only : getavu, relavu, opnfil
  use shr_nl_mod    , only : shr_nl_find_group_name
  use shr_mpi_mod   , only : shr_mpi_bcast
    
  use FuncPedotransferMod,  only : init_pedof
  use RootBiophysMod,       only : init_rootprof
  use SoilWaterMovementMod, only : init_soilwater_movement
  use SurfaceResistanceMod, only : init_soil_stress
  use SoilMoistStressMod,   only : init_root_moist_stress
implicit none

  character(len=*), intent(IN) :: NLFilename ! Namelist filename

  !In future versions, a namelist will be created here to
  !set up options for different sub-models, the namelist file
  !will also be passed into this different initializing methods
  !to read in their local parameters, Jinyun Tang, Mar 29, 2014
  
  call init_pedof

  call init_rootprof
  
  call init_soilwater_movement
  
  call init_soil_stress
  
end subroutine init_hydrology
