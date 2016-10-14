
module MultiPhysicsProbThermal

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for thermal model
  !-----------------------------------------------------------------------

  ! !USES:
  use mpp_varctl                         , only : iulog
  use mpp_abortutils                     , only : endrun
  use mpp_shr_log_mod                    , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbBaseType           , only : multiphysicsprob_base_type
  use SystemOfEquationsThermalType       , only : sysofeqns_thermal_type
  use SystemOfEquationsBasePointerType   , only : sysofeqns_base_pointer_type

  implicit none
  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscts.h"
#include "finclude/petscts.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscsnes.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmda.h90"
#include "finclude/petscviewer.h"

  type, public, extends(multiphysicsprob_base_type) :: mpp_thermal_type
     class(sysofeqns_thermal_type),pointer          :: sysofeqns
     type(sysofeqns_base_pointer_type), pointer     :: sysofeqns_ptr
   contains
     procedure, public :: Init                        => ThermalMPPInit
     procedure, public :: AddGovEqn                   => ThermalMPPAddGovEqn
     procedure, public :: GovEqnAddCondition          => ThermalMPPGovEqnAddCondition
     procedure, public :: SetMeshesOfGoveqns          => ThermalMPPSetMeshesOfGoveqns
     procedure, public :: GovEqnAddCouplingCondition  => ThermalMPPGovEqnAddCouplingCondition
     procedure, public :: AllocateAuxVars             => ThermalMPPAllocateAuxVars
     procedure, public :: GovEqnSetCouplingVars       => ThermalMPPGovEqnSetCouplingVars
     procedure, public :: SetupProblem                => ThermalMPPSetupProblem
     procedure, public :: GovEqnUpdateBCConnectionSet => ThermalMPPGovEqnUpdateBCConnectionSet
     procedure, public :: SetMPIRank                  => ThermalMPPSetMPIRank

  end type mpp_thermal_type

  public :: MPPThermalSetSoils

  type(mpp_thermal_type), public, target :: thermal_mpp

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermalMPPInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize the thermal MPP
    !
    use MultiPhysicsProbBaseType , only : MPPBaseInit
    use MultiPhysicsProbConstants, only : MPP_THERMAL_TBASED_KSP_CLM
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this

    call MPPBaseInit(this)

    allocate(this%sysofeqns)
    call this%sysofeqns%Init()

    allocate(this%sysofeqns_ptr)
    nullify(this%sysofeqns_ptr%ptr)

  end subroutine ThermalMPPInit

  !------------------------------------------------------------------------
  subroutine ThermalMPPSetMPIRank(this, rank)
    !
    ! !DESCRIPTION:
    ! Sets MPI rank
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    PetscInt                :: rank

    if (associated(this%sysofeqns)) then
       this%sysofeqns%mpi_rank = rank
    endif

  end subroutine ThermalMPPSetMPIRank

  !------------------------------------------------------------------------
  subroutine MPPThermalSetSoils(therm_mpp, begc, endc, filter_thermal, &
       lun_type, watsat, csol, tkmg, tkdry)
    !
    ! !DESCRIPTION:
    ! Sets soil thermal properties
    !
    use GoverningEquationBaseType
    use GoveqnThermalKSPTemperatureSoilType
    use ThermalKSPTemperatureSoilAuxType
    use mpp_varpar      , only : nlevgrnd, nlevsoi, nlevsno
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                             :: therm_mpp
    integer              , intent(in)                   :: begc,endc
    integer, intent(in)                                 :: filter_thermal(:)
    PetscInt, pointer, intent(in)                       :: lun_type(:)
    PetscReal, pointer, intent(in)                      :: watsat(:,:)
    PetscReal, pointer, intent(in)                      :: csol(:,:)
    PetscReal, pointer, intent(in)                      :: tkmg(:,:)
    PetscReal, pointer, intent(in)                      :: tkdry(:,:)

    !
    ! !LOCAL VARIABLES:
    class (goveqn_thermal_ksp_temp_soil_type) , pointer :: goveq_soil
    class(sysofeqns_thermal_type)             , pointer :: therm_soe
    class(goveqn_base_type)                   , pointer :: cur_goveq
    type (therm_ksp_temp_soil_auxvar_type)    , pointer :: aux_vars_in(:)
    PetscInt                                            :: j,c,g,l
    PetscInt                                            :: icell
    PetscInt                                            :: col_id
    PetscInt                                            :: first_active_col_id
    PetscBool                                           :: found

    therm_soe => therm_mpp%sysofeqns

    found = PETSC_FALSE
    cur_goveq => therm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_soil_type)
          goveq_soil => cur_goveq
          found = PETSC_TRUE
       end select

       cur_goveq => cur_goveq%next
    enddo

    if (.not.found) then
       write(iulog,*)'Thermal governing equation for soil NOT found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    aux_vars_in => goveq_soil%aux_vars_in

    first_active_col_id = -1
    do c = begc, endc

       if (filter_thermal(c) == 1) then
          if (first_active_col_id == -1) then
             first_active_col_id = c
             exit
          endif
       endif
    enddo

    if (first_active_col_id == -1) then
       write(iulog,*)'No active soil column found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Set thermal properties

    icell = 0
    do c = begc, endc

       ! Soil layers
      do j = 1, nlevgrnd

         icell = icell + 1

         if (filter_thermal(c) == 1) then
            col_id                       = c
            aux_vars_in(icell)%is_active = PETSC_TRUE
         else
            col_id                       = first_active_col_id 
            aux_vars_in(icell)%is_active = PETSC_FALSE
         endif

         if (j > nlevsoi) then
            aux_vars_in(icell)%is_soil_shallow = PETSC_FALSE
         else
            aux_vars_in(icell)%is_soil_shallow = PETSC_TRUE
         endif         

         aux_vars_in(icell)%itype                 = lun_type(col_id)
         aux_vars_in(icell)%por                   = watsat(col_id,j)
         aux_vars_in(icell)%therm_cond_minerals   = tkmg(col_id,j)
         aux_vars_in(icell)%therm_cond_dry        = tkdry(col_id,j)
         aux_vars_in(icell)%heat_cap_minerals_puv = csol(col_id,j)

      enddo
   enddo

  end subroutine MPPThermalSetSoils

  !------------------------------------------------------------------------
  subroutine ThermalMPPAddGovEqn(this, geq_type, name, mesh_itype)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to the MPP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    PetscInt                :: geq_type
    character(len =*)       :: name
    PetscInt                :: mesh_itype

    call this%sysofeqns%AddGovEqn(geq_type, name, mesh_itype)

  end subroutine ThermalMPPAddGovEqn

  !------------------------------------------------------------------------
  subroutine ThermalMPPGovEqnAddCondition(this, igoveqn, ss_or_bc_type, name, unit, &
       cond_type, region_type, id_of_other_goveq)
    !
    ! !DESCRIPTION:
    ! Adds a boundary/source-sink condition to a governing equation
    !
    use GoverningEquationBaseType, only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    PetscInt                          :: igoveqn
    PetscInt                          :: ss_or_bc_type
    character(len =*)                 :: name
    character(len =*)                 :: unit
    PetscInt                          :: cond_type
    PetscInt                          :: region_type
    PetscInt, optional                :: id_of_other_goveq
    !
    class(goveqn_base_type),pointer   :: cur_goveq
    class(goveqn_base_type),pointer   :: other_goveq
    PetscInt                          :: ii

    if (igoveqn > this%sysofeqns%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%sysofeqns%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    if (.not.present(id_of_other_goveq)) then
       call cur_goveq%AddCondition(ss_or_bc_type, name, unit, &
            cond_type, region_type)
    else

       other_goveq => this%sysofeqns%goveqns
       do ii = 1,id_of_other_goveq-1
          other_goveq => other_goveq%next
       enddo

       call cur_goveq%AddCondition(ss_or_bc_type, name, unit, &
            cond_type, region_type, id_of_other_goveq, other_goveq%id )
    endif

  end subroutine ThermalMPPGovEqnAddCondition

  !------------------------------------------------------------------------
  subroutine ThermalMPPSetMeshesOfGoveqns(this)
    !
    ! !DESCRIPTION:
    ! Set association of governing equations and meshes
    !
    use GoverningEquationBaseType, only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this

    call this%sysofeqns%SetMeshesOfGoveqns(this%meshes, this%nmesh)

  end subroutine ThermalMPPSetMeshesOfGoveqns

  !------------------------------------------------------------------------
  subroutine ThermalMPPGovEqnAddCouplingCondition(this, ieqn_1, ieqn_2, &
       iregion_1, iregion_2)
    !
    ! !DESCRIPTION:
    ! Adds a boundary condition to couple ieqn_1 and ieqn_2
    !
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    PetscInt                :: ieqn_1
    PetscInt                :: ieqn_2
    PetscInt                :: iregion_1
    PetscInt                :: iregion_2
    !
    character(len=256)        :: name

    write(name,*) ieqn_2
    name = 'BC_for_coupling_with_equation_' // trim(adjustl(name))
    call this%GovEqnAddCondition(ieqn_1, COND_BC, &
         name, '[K]', COND_DIRICHLET_FRM_OTR_GOVEQ, iregion_1, ieqn_2)

    write(name,*) ieqn_1
    name = 'BC_for_coupling_with_equation_' // trim(adjustl(name))
    call this%GovEqnAddCondition(ieqn_2, COND_BC, &
         name, '[K]', COND_DIRICHLET_FRM_OTR_GOVEQ, iregion_2, ieqn_1)

  end subroutine ThermalMPPGovEqnAddCouplingCondition

  !------------------------------------------------------------------------
  subroutine ThermalMPPAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates auxvars for governing equations and system-of-governing-eqns
    !
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnThermalKSPTemperatureSnowType , only : goveqn_thermal_ksp_temp_snow_type
    use GoveqnThermalKSPTemperatureSSWType  , only : goveqn_thermal_ksp_temp_ssw_type
    use GoveqnThermalKSPTemperatureSoilType , only : goveqn_thermal_ksp_temp_soil_type
    use MultiPhysicsProbConstants           , only : COND_BC
    use MultiPhysicsProbConstants           , only : COND_SS
    use MultiPhysicsProbConstants           , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                :: this
    !
    class(sysofeqns_base_type), pointer    :: soe_base
    class(sysofeqns_thermal_type), pointer :: soe
    class(goveqn_base_type), pointer       :: cur_goveq
    PetscInt                               :: igoveqn
    PetscInt                               :: num_bc
    PetscInt                               :: num_ss
    PetscInt                               :: icond
    PetscInt                               :: iauxvar
    PetscInt                               :: iauxvar_beg, iauxvar_end
    PetscInt                               :: iauxvar_beg_bc, iauxvar_end_bc
    PetscInt                               :: iauxvar_beg_ss, iauxvar_end_ss
    PetscInt                               :: count_bc, count_ss
    PetscInt                               :: offset_bc, offset_ss
    PetscInt, pointer                      :: ncells_for_bc(:)
    PetscInt, pointer                      :: ncells_for_ss(:)
    PetscInt, pointer                      :: offsets_bc(:)
    PetscInt, pointer                      :: offsets_ss(:)

    soe_base => this%sysofeqns

    select type(soe_base)
    class is(sysofeqns_thermal_type)
       soe => this%sysofeqns
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    !
    ! Pass-1: Determine total number of BCs (excluding BCs
    !         needed for coupling various governing equations)
    !         and SSs for all governing equations
    !
    igoveqn = 0
    cur_goveq => soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_snow_type)
          call cur_goveq%AllocateAuxVars()
          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, &
               num_ss, ncells_for_ss)

       class is (goveqn_thermal_ksp_temp_ssw_type)
          call cur_goveq%AllocateAuxVars()
          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, &
               num_ss, ncells_for_ss)

       class is (goveqn_thermal_ksp_temp_soil_type)
          call cur_goveq%AllocateAuxVars()
          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, &
               num_ss, ncells_for_ss)

       end select

       igoveqn = igoveqn + 1

       soe%num_auxvars_in = soe%num_auxvars_in + &
            cur_goveq%mesh%ncells_all

       do icond = 1, num_bc
          soe%num_auxvars_bc = soe%num_auxvars_bc + &
               ncells_for_bc(icond)
       enddo

       do icond = 1, num_ss
          soe%num_auxvars_ss = soe%num_auxvars_ss + &
               ncells_for_ss(icond)
       enddo

       if (num_bc > 0) deallocate(ncells_for_bc)
       if (num_ss > 0) deallocate(ncells_for_ss)

       cur_goveq => cur_goveq%next
    enddo

    ! Allocate memory
    allocate(soe%aux_vars_in           (soe%num_auxvars_in))
    allocate(soe%aux_vars_bc           (soe%num_auxvars_bc))
    allocate(soe%aux_vars_ss           (soe%num_auxvars_ss))

    allocate(soe%soe_auxvars_bc_offset (soe%num_auxvars_bc))
    allocate(soe%soe_auxvars_ss_offset (soe%num_auxvars_bc))
    allocate(soe%soe_auxvars_bc_ncells (soe%num_auxvars_bc))
    allocate(soe%soe_auxvars_ss_ncells (soe%num_auxvars_ss))

    igoveqn        = 0
    iauxvar_beg    = 0
    iauxvar_end    = 0
    iauxvar_beg_bc = 0
    iauxvar_end_bc = 0
    iauxvar_beg_ss = 0
    iauxvar_end_ss = 0
    count_bc       = 0
    count_ss       = 0
    offset_bc      = 0
    offset_ss      = 0

    !
    ! Pass-2: Set values for auxvars of SoE
    !
    cur_goveq => soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_snow_type)
          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, &
               num_ss, ncells_for_ss)

       class is (goveqn_thermal_ksp_temp_ssw_type)
          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, &
               num_ss, ncells_for_ss)

       class is (goveqn_thermal_ksp_temp_soil_type)
          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, &
               num_ss, ncells_for_ss)

       end select

       igoveqn = igoveqn + 1

       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = iauxvar_end + cur_goveq%mesh%ncells_all

       do iauxvar = iauxvar_beg, iauxvar_end
          call soe%aux_vars_in(iauxvar)%Init()

          soe%aux_vars_in(iauxvar)%is_in     = PETSC_TRUE
          soe%aux_vars_in(iauxvar)%goveqn_id = igoveqn
       enddo

       allocate(offsets_bc(num_bc))
       allocate(offsets_ss(num_ss))

       do icond = 1, num_bc
          count_bc = count_bc + 1

          soe%soe_auxvars_bc_offset(count_bc) = offset_bc
          offsets_bc(icond)                   = offset_bc
          soe%soe_auxvars_bc_ncells(count_bc) = ncells_for_bc(icond)
          offset_bc                           = offset_bc + ncells_for_bc(icond)

          iauxvar_beg_bc = iauxvar_end_bc + 1
          iauxvar_end_bc = iauxvar_end_bc + ncells_for_bc(icond)

          do iauxvar = iauxvar_beg_bc, iauxvar_end_bc
             call soe%aux_vars_bc(iauxvar)%Init()

             soe%aux_vars_bc(iauxvar)%is_in        = PETSC_TRUE
             soe%aux_vars_bc(iauxvar)%goveqn_id    = igoveqn
             soe%aux_vars_bc(iauxvar)%condition_id = icond
          enddo
       enddo

       do icond = 1, num_ss
          count_ss = count_ss + 1
          soe%soe_auxvars_ss_offset(count_ss) = offset_ss
          offsets_ss(icond)                   = offset_ss
          soe%soe_auxvars_ss_ncells(count_ss) = ncells_for_ss(icond)
          offset_ss                           = offset_ss + ncells_for_ss(icond)

          iauxvar_beg_ss = iauxvar_end_ss + 1
          iauxvar_end_ss = iauxvar_end_ss + ncells_for_ss(icond)

          do iauxvar = iauxvar_beg_ss, iauxvar_end_ss
             call soe%aux_vars_ss(iauxvar)%Init()

             soe%aux_vars_ss(iauxvar)%is_in        = PETSC_TRUE
             soe%aux_vars_ss(iauxvar)%goveqn_id    = igoveqn
             soe%aux_vars_ss(iauxvar)%condition_id = icond
          enddo
       enddo

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_snow_type)
          call cur_goveq%SetSOEAuxVarOffsets(num_bc, offsets_bc, num_ss, offsets_ss)

       class is (goveqn_thermal_ksp_temp_ssw_type)
          call cur_goveq%SetSOEAuxVarOffsets(num_bc, offsets_bc, num_ss, offsets_ss)

       class is (goveqn_thermal_ksp_temp_soil_type)
          call cur_goveq%SetSOEAuxVarOffsets(num_bc, offsets_bc, num_ss, offsets_ss)

       end select

       if (num_bc > 0) deallocate(ncells_for_bc)
       if (num_ss > 0) deallocate(ncells_for_ss)

       deallocate(offsets_bc)
       deallocate(offsets_ss)

       cur_goveq => cur_goveq%next
    enddo

  end subroutine ThermalMPPAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine ThermalMPPGovEqnSetCouplingVars(this, igoveqn, nvars, &
       var_ids, goveqn_ids)
    !
    ! !DESCRIPTION:
    ! In order to couple the given governing equation, add:
    ! - ids of variables needed, and
    ! - ids of governing equations from which variables are needed.
    ! needed for coupling
    ! 
    ! !USES:
    use ConditionType             , only : condition_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)           :: this
    PetscInt                          :: igoveqn
    PetscInt                          :: nvars
    PetscInt, pointer                 :: var_ids(:)
    PetscInt, pointer                 :: goveqn_ids(:)
    !
    class(goveqn_base_type) , pointer :: cur_goveq_1
    class(goveqn_base_type) , pointer :: cur_goveq_2
    type(condition_type)    , pointer :: cur_cond_1
    type(condition_type)    , pointer :: cur_cond_2
    PetscInt                          :: ii
    PetscInt                          :: ieqn
    PetscInt                          :: ivar
    PetscInt                          :: bc_idx_1
    PetscInt                          :: bc_idx_2
    PetscInt                          :: bc_offset_1
    PetscBool                         :: bc_found

    if (igoveqn > this%sysofeqns%ngoveqns) then
       write(iulog,*) 'Attempting to set coupling vars for governing ' // &
            'equation that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq_1 => this%sysofeqns%goveqns
    do ii = 1, igoveqn-1
       cur_goveq_1 => cur_goveq_1%next
    end do

    call cur_goveq_1%AllocVarsFromOtherGEs(nvars)

    do ivar = 1, nvars

       if (goveqn_ids(ivar) > this%sysofeqns%ngoveqns) then
          write(iulog,*) 'Attempting to set coupling vars to a governing ' // &
               'equation that is not in the list'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       bc_found    = PETSC_FALSE
       bc_idx_1    = 1
       bc_offset_1 = 0

       cur_cond_1 => cur_goveq_1%boundary_conditions%first
       do
          if (.not.associated(cur_cond_1)) exit

          ! Is this the appropriate BC?
          if (cur_cond_1%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
             do ieqn = 1, cur_cond_1%num_other_goveqs
                if (cur_cond_1%list_id_of_other_goveqs(ieqn) == goveqn_ids(ivar) ) then
                   bc_found = PETSC_TRUE
                   exit
                endif
             enddo
          endif

          if (bc_found) exit

          bc_idx_1    = bc_idx_1    + 1
          bc_offset_1 = bc_offset_1 + cur_cond_1%conn_set%num_connections

          cur_cond_1 => cur_cond_1%next
       enddo

       if (.not.bc_found) then
          write(iulog,*)'For goveqn%name = ',trim(cur_goveq_1%name) // &
               ', no coupling boundary condition found to copule it with ' // &
               'equation_number = ', goveqn_ids(ivar)
       endif

       cur_goveq_2 => this%sysofeqns%goveqns
       do ii = 1, goveqn_ids(ivar)-1
          cur_goveq_2 => cur_goveq_2%next
       enddo

       bc_found    = PETSC_FALSE
       bc_idx_2    = 1

       cur_cond_2 => cur_goveq_2%boundary_conditions%first
       do
          if (.not.associated(cur_cond_2)) exit

          ! Is this the appropriate BC?
          if (cur_cond_2%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
             do ieqn = 1, cur_cond_2%num_other_goveqs
                if (cur_cond_2%list_id_of_other_goveqs(ieqn) == igoveqn ) then
                   bc_found = PETSC_TRUE
                   exit
                endif
             enddo
          endif

          if (bc_found) exit

          bc_idx_2    = bc_idx_2    + 1

          cur_cond_2 => cur_cond_2%next
       enddo

       if (.not.bc_found) then
          write(iulog,*)'For goveqn%name = ',trim(cur_goveq_2%name) // &
               ', no coupling boundary condition found to copule it with ' // &
               'equation_number = ', bc_idx_2
       endif

       cur_goveq_1%var_ids_needed_from_other_goveqns (ivar) = var_ids(ivar)
       cur_goveq_1%ids_of_other_goveqns              (ivar) = goveqn_ids(ivar)
       cur_goveq_1%is_bc_auxvar_type                 (ivar) = PETSC_TRUE
       cur_goveq_1%bc_auxvar_offset                  (ivar) = bc_offset_1
       cur_goveq_1%bc_auxvar_ncells                  (ivar) = cur_cond_1%conn_set%num_connections
       cur_goveq_1%bc_auxvar_idx                     (ivar) = bc_idx_1
       cur_goveq_1%bc_auxvar_idx_of_other_goveqn     (ivar) = bc_idx_2

    enddo

  end subroutine ThermalMPPGovEqnSetCouplingVars

  !------------------------------------------------------------------------
  subroutine ThermalMPPSetupProblem(this)
    !
    ! !DESCRIPTION:
    ! Sets up the thermal MPP
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : SOE_THERMAL_TBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)           :: this

    call ThermalMPPUpdatCouplingBCConnections(this)
    call ThermalMPPKSPSetup(this)

    this%sysofeqns%solver_type = this%solver_type
    this%sysofeqns%itype       = SOE_THERMAL_TBASED

  end subroutine ThermalMPPSetupProblem

  !------------------------------------------------------------------------
  subroutine ThermalMPPUpdatCouplingBCConnections(this)
    !
    ! !DESCRIPTION:
    ! For BCs used to coupling two governing equations, updates
    ! ID and distance of upwind grid cell.
    !
    use ConditionType             , only : condition_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                 :: this
    !
    class(sysofeqns_base_type)    , pointer :: soe_base
    class(sysofeqns_thermal_type) , pointer :: soe
    class(goveqn_base_type)       , pointer :: cur_goveq_1
    class(goveqn_base_type)       , pointer :: cur_goveq_2
    type(condition_type)          , pointer :: cur_cond_1
    type(condition_type)          , pointer :: cur_cond_2
    type(connection_set_type)     , pointer :: cur_conn_set_1
    type(connection_set_type)     , pointer :: cur_conn_set_2
    PetscInt                                :: igoveqn
    PetscInt                                :: ieqn
    PetscInt                                :: iconn
    PetscInt                                :: ii, jj
    PetscInt                                :: ivar
    PetscInt                                :: bc_idx_1
    PetscInt                                :: bc_idx_2
    PetscInt                                :: bc_offset_1
    PetscBool                               :: bc_found

    soe_base => this%sysofeqns

    select type(soe_base)
    class is(sysofeqns_thermal_type)
       soe => this%sysofeqns
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    do ii = 1, soe%ngoveqns

       cur_goveq_1 => soe%goveqns
       do igoveqn = 1, ii-1
          cur_goveq_1 => cur_goveq_1%next
       end do


       do jj = ii+1, soe%ngoveqns

          cur_goveq_2 => soe%goveqns
          do igoveqn = 1, jj-1
             cur_goveq_2 => cur_goveq_2%next
          end do

          bc_found    = PETSC_FALSE
          cur_cond_1 => cur_goveq_1%boundary_conditions%first
          do
             if (.not.associated(cur_cond_1)) exit

             ! Is this the appropriate BC?
             if (cur_cond_1%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
                do ieqn = 1, cur_cond_1%num_other_goveqs
                   if (cur_cond_1%list_id_of_other_goveqs(ieqn) == jj ) then
                      bc_found = PETSC_TRUE
                      exit
                   endif
                enddo
             endif

             if (bc_found) exit
             cur_cond_1 => cur_cond_1%next
          enddo

          if (bc_found) then

             bc_found    = PETSC_FALSE
             cur_cond_2 => cur_goveq_2%boundary_conditions%first
             do
                if (.not.associated(cur_cond_2)) exit

                ! Is this the appropriate BC?
                if (cur_cond_2%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
                   do ieqn = 1, cur_cond_2%num_other_goveqs
                      if (cur_cond_2%list_id_of_other_goveqs(ieqn) == ii ) then
                         bc_found = PETSC_TRUE
                         exit
                      endif
                   enddo
                endif

                if (bc_found) exit
                cur_cond_2 => cur_cond_2%next
             enddo

             if (.not.bc_found) then
                write(iulog,*)'cur_goveq_1: ',trim(cur_goveq_1%name)
                write(iulog,*)'cur_goveq_2: ',trim(cur_goveq_2%name)
                write(iulog,*)'cur_goveq_1 has a coupling BC for cur_eqn_2, but not vice-versa.'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             cur_conn_set_1 => cur_cond_1%conn_set
             cur_conn_set_2 => cur_cond_2%conn_set

             if (cur_conn_set_1%num_connections /= cur_conn_set_2%num_connections) then
                write(iulog,*) 'Number of connections in two equations are not same.'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             do iconn = 1, cur_conn_set_1%num_connections
                cur_conn_set_2%id_up(iconn)   = cur_conn_set_1%id_dn(iconn)
                cur_conn_set_2%dist_up(iconn) = cur_conn_set_1%dist_dn(iconn)

                cur_conn_set_1%id_up(iconn)   = cur_conn_set_2%id_dn(iconn)
                cur_conn_set_1%dist_up(iconn) = cur_conn_set_2%dist_dn(iconn)
             enddo

          endif

       enddo
    enddo

  end subroutine ThermalMPPUpdatCouplingBCConnections

  !------------------------------------------------------------------------
  subroutine ThermalMPPKSPSetup(therm_mpp)
    !
    ! !DESCRIPTION:
    ! Sets the PETSc KSP for the thermal mpp
    !
    use GoverningEquationBaseType , only : goveqn_base_type
    use SystemOfEquationsBasePointerType, only : SOEComputeRHS, SOEComputeOperators
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)               :: therm_mpp
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer       :: cur_goveq
    class(sysofeqns_thermal_type),pointer :: therm_soe
    PetscInt                              :: size
    PetscInt                              :: igoveq
    PetscErrorCode                        :: ierr
    DM, pointer                           :: dms(:)
    character(len=256)                    :: name

    therm_soe => therm_mpp%sysofeqns
    therm_mpp%sysofeqns_ptr%ptr => therm_mpp%sysofeqns

    ! Create PETSc DM for each governing equation

    allocate(dms(therm_soe%ngoveqns))

    igoveq = 0
    cur_goveq => therm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       igoveq = igoveq + 1
       size   = cur_goveq%mesh%ncells_local

       call DMDACreate1d(PETSC_COMM_SELF, &
            DM_BOUNDARY_NONE, size, 1, 0, &
            PETSC_NULL_INTEGER, dms(igoveq), ierr);
       CHKERRQ(ierr)

       write(name,*) igoveq
       name = 'fgoveq_' // trim(adjustl(name))
       call DMSetOptionsPrefix(dms(igoveq), name, ierr); CHKERRQ(ierr)

       write(name,*) igoveq
       name = 'goveq_' // trim(adjustl(name))
       call DMDASetFieldName(dms(igoveq), 0, name, ierr); CHKERRQ(ierr)

       call DMSetFromOptions(dms(igoveq), ierr); CHKERRQ(ierr)

       cur_goveq => cur_goveq%next
    enddo

    ! DM-Composite approach

    ! Create DMComposite: temperature
    call DMCompositeCreate(PETSC_COMM_SELF, therm_soe%dm, ierr); CHKERRQ(ierr)
    call DMSetOptionsPrefix(therm_soe%dm, "temperature_", ierr); CHKERRQ(ierr)

    ! Add DMs to DMComposite
    do igoveq = 1, therm_soe%ngoveqns
       call DMCompositeAddDM(therm_soe%dm, dms(igoveq), ierr); CHKERRQ(ierr)
    enddo

    ! Setup DM
    call DMSetUp(therm_soe%dm, ierr); CHKERRQ(ierr)

    ! Create matrix
    call DMCreateMatrix    (therm_soe%dm   , therm_soe%Amat, ierr); CHKERRQ(ierr)

    call MatSetOption      (therm_soe%Amat , MAT_NEW_NONZERO_LOCATION_ERR , &
         PETSC_FALSE, ierr); CHKERRQ(ierr)
    call MatSetOption      (therm_soe%Amat , MAT_NEW_NONZERO_ALLOCATION_ERR, &
         PETSC_FALSE, ierr); CHKERRQ(ierr)

    call MatSetFromOptions (therm_soe%Amat , ierr); CHKERRQ(ierr)

    ! Create vectors
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%soln         , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%rhs          , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%soln_prev    , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%soln_prev_clm, ierr); CHKERRQ(ierr)

    ! Initialize vectors
    call VecZeroEntries(therm_soe%soln          , ierr); CHKERRQ(ierr)
    call VecZeroEntries(therm_soe%rhs           ,  ierr); CHKERRQ(ierr)
    call VecZeroEntries(therm_soe%soln_prev     ,  ierr); CHKERRQ(ierr)
    call VecZeroEntries(therm_soe%soln_prev_clm ,  ierr); CHKERRQ(ierr)

    ! Create KSP
    call KSPCreate              (PETSC_COMM_SELF , therm_soe%ksp, ierr); CHKERRQ(ierr)
    call KSPSetOptionsPrefix    (therm_soe%ksp   , "temperature_", ierr); CHKERRQ(ierr)

    call KSPSetComputeRHS       (therm_soe%ksp   , SOEComputeRHS      , &
         therm_mpp%sysofeqns_ptr, ierr); CHKERRQ(ierr)
    call KSPSetComputeOperators (therm_soe%ksp   , SOEComputeOperators, &
         therm_mpp%sysofeqns_ptr, ierr); CHKERRQ(ierr)

    call KSPSetFromOptions      (therm_soe%ksp   , ierr); CHKERRQ(ierr)

    ! Cleanup
    do igoveq = 1, therm_soe%ngoveqns
       call DMDestroy(dms(igoveq), ierr); CHKERRQ(ierr)
    enddo
    deallocate(dms)

  end subroutine ThermalMPPKSPSetup

  !------------------------------------------------------------------------
  subroutine ThermalMPPGovEqnUpdateBCConnectionSet(this, igoveqn, icond, &
       var_type, nval, values)
    !
    ! !DESCRIPTION:
    ! For a boundary condition of a given governing equation, update distance
    ! for a downstream cell.
    !
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use GoverningEquationBaseType, only : goveqn_base_type
    use MultiPhysicsProbConstants, only : VAR_DIST_DN
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    PetscInt :: igoveqn
    PetscInt :: icond
    PetscInt :: nval
    PetscInt :: var_type
    PetscReal, pointer :: values (:)
    !
    class(goveqn_base_type),pointer   :: cur_goveq
    type(condition_type)    , pointer :: cur_cond
    type(connection_set_type)     , pointer :: cur_conn_set
    PetscInt :: ii
    PetscInt :: iconn
    PetscInt :: bc_idx
    PetscBool :: bc_found

    if (igoveqn > this%sysofeqns%ngoveqns) then
       write(iulog,*) 'Attempting to access governing equation that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%sysofeqns%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    bc_found = PETSC_FALSE
    bc_idx = 0
    cur_cond => cur_goveq%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit

       bc_idx = bc_idx + 1
       if (bc_idx == icond) then
          bc_found = PETSC_TRUE

          cur_conn_set => cur_cond%conn_set
          if (nval /= cur_conn_set%num_connections) then
             write(iulog,*) 'Number of values to update connections ' // &
                  'do not match number of connections.'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif

          do iconn = 1, cur_conn_set%num_connections

             select case(var_type)
             case (VAR_DIST_DN)
                cur_conn_set%dist_dn(iconn) = values(iconn)
             case default
                write(iulog,*) 'Unknown variable type'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select
          enddo

          exit

       end if

       cur_cond => cur_cond%next
    enddo

    if (.not.bc_found) then
       write(iulog,*) 'Failed to find icond = ',icond,' in the boundary condition list.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif


  end subroutine ThermalMPPGovEqnUpdateBCConnectionSet

#endif

end module MultiPhysicsProbThermal
