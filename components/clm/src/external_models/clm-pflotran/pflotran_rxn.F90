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
! Los Alamos National Laboratory
! Earth and Environmental Sciences
! EES-16, MS: D469
! (505) 667-3420
! lichtner@lanl.gov
! Los Alamos, NM

! or

! Glenn E. Hammond
! Sandia National Laboratories
! Applied Systems Analysis & Research
! 413 Cherry Blossom Lp
! Richland, WA 99352
! (505) 235-0665
! gehammo@sandia.gov

!=======================================================================


module BatchChem

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: BatchChemInitializeReactions, &
            BatchChemProcessConstraints

contains

! ************************************************************************** !

subroutine BatchChemInitializeReactions(option, input, reaction)

  use Reaction_module
  use Reaction_Aux_module
  use Reaction_Database_module
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petsclog.h"

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(reaction_type), pointer :: reaction
  character(len=MAXSTRINGLENGTH) :: string

  ! check for a chemistry block in the  input file
  string = "CHEMISTRY"
  call InputFindStringInFile(input, option, string)
  if (.not.InputError(input)) then
    ! found a chemistry block, initialize the chemistry.

    ! NOTE(bja): ReactionInit() only does a first pass through the
    ! input file to check for a select items
    call ReactionInit(reaction, input, option)
    ! rewind the input file to prepare for the second pass
    call InputFindStringInFile(input, option, string)
    ! the second pass through the input file to read the remaining blocks
    call ReactionReadPass2(reaction, input, option)
  else
     ! TODO(bja): no chemistry block --> fatal error
  endif
    
  if (associated(reaction)) then
    if (reaction%use_full_geochemistry) then
       call DatabaseRead(reaction, option)
       call BasisInit(reaction, option)    
    else
      ! NOTE(bja): do we need this for the batch chemistry driver?

      ! turn off activity coefficients since the database has not been read
      reaction%act_coef_update_frequency = ACT_COEF_FREQUENCY_OFF
      allocate(reaction%primary_species_print(option%ntrandof))
      reaction%primary_species_print = PETSC_TRUE
    endif
  endif

end subroutine BatchChemInitializeReactions

! ************************************************************************** !

subroutine BatchChemProcessConstraints(option, input, reaction, &
     global_auxvars, rt_auxvars, material_auxvars, transport_constraints, &
     constraint_coupler)

  use Reaction_module
  use Reaction_Aux_module
  use Reaction_Database_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Transport_Constraint_module
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petsclog.h"

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(reaction_type), pointer :: reaction
  character(len=MAXSTRINGLENGTH) :: string
  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars
  class(material_auxvar_type), pointer :: material_auxvars

  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  type(tran_constraint_type), pointer :: tran_constraint
  type(tran_constraint_list_type), pointer :: transport_constraints
  type(tran_constraint_coupler_type), pointer :: constraint_coupler 
  PetscBool :: use_prev_soln_as_guess
  PetscInt :: num_iterations
  

  !
  ! read the constraints...
  !

  ! look through the input file
  rewind(input%fid)        
  do
    call InputReadPflotranString(input, option)
    if (InputError(input)) exit

    call InputReadWord(input, option, word, PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    option%io_buffer = 'pflotran card:: ' // trim(card)
    call printMsg(option)

    select case(trim(card))
      case('CONSTRAINT')
        if (.not.associated(reaction)) then
          option%io_buffer = 'CONSTRAINTs not supported without CHEMISTRY.'
          call printErrMsg(option)
        endif
        tran_constraint => TranConstraintCreate(option)
        call InputReadWord(input, option, tran_constraint%name, PETSC_TRUE)
        call InputErrorMsg(input, option, 'constraint', 'name') 
        call printMsg(option, tran_constraint%name)
        call TranConstraintRead(tran_constraint, reaction, input, option)
        call TranConstraintAddToList(tran_constraint, transport_constraints)
        nullify(tran_constraint)

      case default
         ! do nothing
    end select
  enddo

  !
  ! process constraints
  !
  num_iterations = 0
  use_prev_soln_as_guess = PETSC_FALSE
  tran_constraint => transport_constraints%first
  ! NOTE(bja): we only created one set of global and rt auxvars, so if
  ! there is more than one constratint in the input file, they will be
  ! over written.
  do 
     if (.not. associated(tran_constraint)) exit
     ! initialize constraints
     option%io_buffer = "initializing constraint : " // tran_constraint%name
     call printMsg(option)
     call ReactionProcessConstraint(reaction, &
                                    tran_constraint%name, &
                                    tran_constraint%aqueous_species, &
                                    tran_constraint%free_ion_guess, &
                                    tran_constraint%minerals, &
                                    tran_constraint%surface_complexes, &
                                    tran_constraint%colloids, &
                                    tran_constraint%immobile_species, &
                                    option)

     ! link the constraint to the constraint coupler
     constraint_coupler%constraint_name = tran_constraint%name
     constraint_coupler%aqueous_species => tran_constraint%aqueous_species
     constraint_coupler%free_ion_guess => tran_constraint%free_ion_guess
     constraint_coupler%minerals => tran_constraint%minerals
     constraint_coupler%surface_complexes => tran_constraint%surface_complexes
     constraint_coupler%colloids => tran_constraint%colloids
     constraint_coupler%global_auxvar => global_auxvars
     constraint_coupler%rt_auxvar => rt_auxvars
     
     ! equilibrate
     option%io_buffer = "equilibrate constraint : " // tran_constraint%name
     call printMsg(option)
     call ReactionEquilibrateConstraint(rt_auxvars, global_auxvars, &
                                        material_auxvars, reaction, &
                                        tran_constraint%name, &
                                        tran_constraint%aqueous_species, &
                                        tran_constraint%free_ion_guess, &
                                        tran_constraint%minerals, &
                                        tran_constraint%surface_complexes, &
                                        tran_constraint%colloids, &
                                        tran_constraint%immobile_species, &
                                        num_iterations, &
                                        use_prev_soln_as_guess, &
                                        option)
     call ReactionPrintConstraint(constraint_coupler, reaction, option)
     tran_constraint => tran_constraint%next
  enddo

end subroutine BatchChemProcessConstraints


end module BatchChem


! ************************************************************************** !
program pflotran_rxn
  
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Database_module
  use Option_module
  use Input_Aux_module
  use String_module
  
  use Transport_Constraint_module
  use PFLOTRAN_Constants_module

  use BatchChem

  implicit none

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petsclog.h"

  PetscErrorCode :: ierr
  PetscBool :: option_found  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: filename_out
  type(reaction_type), pointer :: reaction
  type(option_type), pointer :: option
  type(input_type), pointer :: input

  type(global_auxvar_type), pointer :: global_auxvars
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars
  class(material_auxvar_type), pointer :: material_auxvars

  character(len=MAXWORDLENGTH) :: card
  character(len=MAXWORDLENGTH) :: word
  type(tran_constraint_list_type), pointer :: transport_constraints
  type(tran_constraint_coupler_type), pointer :: constraint_coupler 

  option => OptionCreate()
  option%fid_out = OUT_UNIT

  call MPI_Init(ierr)
  option%global_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(MPI_COMM_WORLD, option%global_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, option%global_commsize, ierr)
  call MPI_Comm_group(MPI_COMM_WORLD, option%global_group, ierr)
  option%mycomm = option%global_comm
  option%myrank = option%global_rank
  option%mycommsize = option%global_commsize
  option%mygroup = option%global_group

  ! check for non-default input filename
  option%input_filename = "pflotran.in"
  string = '-pflotranin'
  call InputGetCommandLineString(string, option%input_filename, option_found, option)

  string = '-output_prefix'
  call InputGetCommandLineString(string, option%global_prefix, option_found, option)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)

  input => InputCreate(IN_UNIT, option%input_filename, option)

  filename_out = trim(option%global_prefix) // trim(option%group_prefix) // &
                 '.out'

  if (option%myrank == option%io_rank .and. option%print_to_file) then
    open(option%fid_out, file=filename_out, action="write", status="unknown")
  endif

  !
  ! manual initialization...
  !
  option%nphase = 1
  option%liquid_phase = 1
  option%reference_water_density = 998.2

  call BatchChemInitializeReactions(option, input, reaction)

  !
  ! create the storage containers
  !
  ! NOTE(bja) : batch chem --> one cell

  ! global_auxvars --> cell by cell temperature, pressure, saturation, density
  allocate(global_auxvars)
  call GlobalAuxVarInit(global_auxvars, option)

  ! rt_auxvars --> cell by cell chemistry data
  allocate(rt_auxvars)
  call RTAuxVarInit(rt_auxvars, reaction, option)

  ! material_auxvars --> cell by cell material property data
  allocate(material_auxvars)
  call MaterialAuxVarInit(material_auxvars, option)
  material_auxvars%porosity = option%reference_porosity

  ! assign default state values
  global_auxvars%pres = option%reference_pressure
  global_auxvars%temp = option%reference_temperature
  ! global_auxvars%den_kg = option%reference_water_density
  ! NOTE(bja): option%ref_density = 0.0, so we set it manually. This is a Bad Thing(TM)
  global_auxvars%den_kg = 998.2
  global_auxvars%sat = option%reference_saturation  

  ! create the constraint list
  allocate(transport_constraints)
  call TranConstraintInitList(transport_constraints)
  allocate(constraint_coupler)
  constraint_coupler => TranConstraintCouplerCreate(option)

  call BatchChemProcessConstraints(option, input, reaction, &
                                   global_auxvars, rt_auxvars, &
                                   material_auxvars, transport_constraints, &
                                   constraint_coupler)

  ! cleanup
  call TranConstraintCouplerDestroy(constraint_coupler)
  call TranConstraintDestroyList(transport_constraints)
  ! FIXME(bja) : causes error freeing memory.
  !call RTAuxVarDestroy(rt_auxvars)
  !call GlobalAuxVarDestroy(global_auxvars)
  call ReactionDestroy(reaction,option)
  call GlobalAuxVarStrip(global_auxvars)
  deallocate(global_auxvars)
  nullify(global_auxvars)
  call RTAuxVarStrip(rt_auxvars)
  deallocate(rt_auxvars)
  nullify(rt_auxvars)
  call MaterialAuxVarStrip(material_auxvars)
  deallocate(material_auxvars)
  nullify(material_auxvars)
  call InputDestroy(input)
  call OptionDestroy(option)
  call PetscFinalize(ierr);CHKERRQ(ierr)
  call MPI_Finalize(ierr)

end program pflotran_rxn

