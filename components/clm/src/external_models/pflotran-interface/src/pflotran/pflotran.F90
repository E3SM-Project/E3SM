!=======================================================================
! PFLOTRAN v2.0 LA-CC-09-047
!=======================================================================

!Copyright 2009. Los Alamos National Security, LLC. This material was produced under U.S. 
!Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated 
!by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has 
!rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS 
!NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE 
!USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software 
!should be clearly marked, so as not to confuse it with the version available from LANL.
!Additionally, this library is free software; you can redistribute it and/or modify it under the 
!terms of the GNU Lesser General Public License as published by the Free Software Foundation; 
!either version 2.1 of the License, or (at your option) any later version. Accordingly, this 
!library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
!the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser 
!General Public License for more details.

! Send all bug reports/questions/comments to:
!
! Peter C. Lichtner
! OFM Research
! (505) 692-4029 (cell)
! peter.lichtner@gmail.com
! Santa Fe, New Mexico

! or

! Glenn E. Hammond
! Sandia National Laboratories
! Applied Systems Analysis & Research
! 413 Cherry Blossom Lp
! Richland, WA 99352
! (505) 235-0665
! gehammo@sandia.gov

!=======================================================================
program pflotran
  
  use Option_module
  use Simulation_Base_class
  use Multi_Simulation_module
  use Factory_PFLOTRAN_module
  use Factory_Subsurface_module

  use PFLOTRAN_Constants_module
  use PFLOTRAN_Provenance_module, only : PrintProvenanceToScreen
  
  implicit none

#include "petsc/finclude/petscsys.h"

  class(simulation_base_type), pointer :: simulation
  ! multisimulation enables multiple simulations to be run concurrently
  ! and/or one after another until a specified set of simulations has 
  ! completed.
  type(multi_simulation_type), pointer :: multisimulation
  type(option_type), pointer :: option
  
  nullify(simulation)
  nullify(multisimulation)
  option => OptionCreate()
  call OptionInitMPI(option)
  call PFLOTRANInitializePrePetsc(multisimulation,option)
  call OptionInitPetsc(option)
  if (option%myrank == option%io_rank .and. option%print_to_screen) then
    call PrintProvenanceToScreen()
  endif

  do ! multi-simulation loop
    call PFLOTRANInitializePostPetsc(simulation,multisimulation,option)

    call simulation%InitializeRun()

    if (option%status == PROCEED) then
      call simulation%ExecuteRun()
    endif

    call simulation%FinalizeRun()
    call simulation%Strip()
    deallocate(simulation)
    nullify(simulation)
    call PFLOTRANFinalize(option)
    if (MultiSimulationDone(multisimulation)) exit
  enddo ! multi-simulation loop
  call OptionFinalize(option)

end program pflotran
