module GoveqnRichardsODEPressureType

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Govneqn data type allocation
  !-----------------------------------------------------------------------

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                    , only : iulog
  use mpp_abortutils                , only : endrun
  use mpp_shr_log_mod               , only : errMsg => shr_log_errMsg
  use GoverningEquationBaseType     , only : goveqn_base_type
  use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
  use RichardsODEPressureConnAuxType, only : rich_ode_pres_conn_auxvar_type
  use petscvec
  use petscmat
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(goveqn_base_type) :: goveqn_richards_ode_pressure_type
     Vec                                             :: accum_prev

     PetscReal                             , pointer :: internal_flux(:)          ! mass flux betwenn internal connections [kg/s]
     PetscReal                             , pointer :: boundary_flux(:)          ! mass flux between boundary connections [kg/s]
     PetscReal                             , pointer :: ss_flux(:)                ! mass flux for source-sink connections [kg/s]
     PetscReal                             , pointer :: lat_mass_exc(:)           ! lateral mass exchanged [kg]
     PetscReal                             , pointer :: bnd_mass_exc(:)           ! mass exchanged through boundary conditions [kg]

     type (rich_ode_pres_auxvar_type)      , pointer :: aux_vars_in(:)            ! Internal state.
     type (rich_ode_pres_auxvar_type)      , pointer :: aux_vars_bc(:)            ! Boundary conditions.
     type (rich_ode_pres_auxvar_type)      , pointer :: aux_vars_ss(:)            ! Source-sink.
     type (rich_ode_pres_conn_auxvar_type) , pointer :: aux_vars_conn_in(:)       ! Internal connection
     type (rich_ode_pres_conn_auxvar_type) , pointer :: aux_vars_conn_bc(:)       ! Boundary connection

     PetscInt                              , pointer :: soe_auxvars_bc_offset (:) ! SoE auxvar offset corresponding to BCs
     PetscInt                              , pointer :: soe_auxvars_ss_offset (:) ! SoE auxvar offset corresponding to SSs

   contains
     procedure, public :: AllocateAuxVars           => RichardsODEPressureAllocateAuxVars
     procedure, public :: SetDensityType            => RichardsODEPressureSetDensityType
     procedure, public :: Setup                     => RichardsODESetup
     procedure, public :: ComputeResidual           => RichardsODEComputeResidual
     procedure, public :: ComputeJacobian           => RichardsODEComputeJacobian
     procedure, public :: ComputeOffDiagJacobian    => RichardsODEComputeOffDiagJacobian

     procedure, public :: GetFromSOEAuxVarsIntrn    => RichardsODEPressureGetFromSOEAuxVarsIntrn
     procedure, public :: SetFromSOEAuxVarsIntrn    => RichardsODEPressureSetFromSOEAuxVarsIntrn
     procedure, public :: GetFromSOEAuxVarsBC       => RichardsODEPressureGetFromSOEAuxVarsBC
     procedure, public :: GetFromSOEAuxVarsSS       => RichardsODEPressureGetFromSOEAuxVarsSS
     procedure, public :: GetDataFromSOEAuxVar      => RichardsODEPressureGetDataFromSOEAuxVar

     procedure, public :: SetDataInSOEAuxVar        => RichardsODEPressureSetDataInSOEAuxVar
     procedure, public :: UpdateAuxVars             => RichardsODEPressureUpdateAuxVars
     procedure, public :: UpdateAuxVarsIntrn        => RichardsODEPressureUpdateAuxVarsIntrn
     procedure, public :: UpdateAuxVarsBC           => RichardsODEPressureUpdateAuxVarsBC
     procedure, public :: UpdateAuxVarsSS           => RichardsODEPressureUpdateAuxVarsSS
     procedure, public :: PreSolve                  => RichardsODEPressurePreSolve
     procedure, public :: PreStepDT                 => RichardsODEPressurePreStepDT
     procedure, public :: NumConditions             => RichardsODEPressureNumConditions
     procedure, public :: NumCellsInConditions      => RichardsODEPressureNumCellsInConditions
     procedure, public :: SetSOEAuxVarOffsets       => RichardsODEPressureSetSOEAuxVarOffsets
     procedure, public :: GetNumInternalConnections => RichardsODEGetNumInternalConnections
     procedure, public :: CreateVectors             => RichardsODEPressureCreateVectors

     procedure, public :: ComputeLateralFlux        => RichardsODEComputeLateralFlux

     procedure, public :: SetAuxVarConnRealValue    => RichardsODESetAuxVarConnRealValue
     procedure, public :: SetAuxVarConnIntValue     => RichardsODESetAuxVarConnIntValue
  end type

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine RichardsODESetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation for Richards equation.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_RE
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this

    call this%Create()

    this%name         = "Richards Equation ODE"
    this%id           = GE_RE
    this%mesh_itype   = MESH_CLM_SOIL_COL

    nullify(this%soe_auxvars_bc_offset)
    nullify(this%soe_auxvars_ss_offset)

  end subroutine RichardsODESetup

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates memory for storing auxiliary variables associated with:
    !   + Internal control volumes,
    !   + Boundary condtions,
    !   + Source-sink condition.
    !
    ! !USES:
    use ConnectionSetType , only : connection_set_type
    use ConditionType     , only : condition_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    !
    type(condition_type)      , pointer      :: cur_cond
    type(connection_set_type) , pointer      :: cur_conn_set
    PetscInt                                 :: ncells_cond
    PetscInt                                 :: sum_conn
    PetscInt                                 :: icond
    PetscInt                                 :: ncond

    ! Allocate memory and initialize aux vars: For internal connections
    allocate(this%aux_vars_in(this%mesh%ncells_all))
    do icond = 1,this%mesh%ncells_all
       call this%aux_vars_in(icond)%Init()
    enddo

    ! Find number of internal connections
    cur_conn_set => this%mesh%intrn_conn_set_list%first
    sum_conn = 0
    do
       if (.not.associated(cur_conn_set)) exit
       sum_conn = sum_conn + cur_conn_set%num_connections
       cur_conn_set => cur_conn_set%next
    enddo
    allocate(this%internal_flux(sum_conn))
    this%internal_flux(:) = 0.d0

    allocate(this%aux_vars_conn_in(sum_conn))
    do icond = 1, sum_conn
       call this%aux_vars_conn_in(icond)%Init()
    enddo

    allocate(this%lat_mass_exc(this%mesh%ncells_all))
    this%lat_mass_exc = 0.d0

    ! Allocate memory and initialize aux vars: For boundary connections
    ncells_cond = 0
    ncond       = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       ncells_cond = ncells_cond + cur_cond%ncells
       ncond       = ncond + 1
       cur_cond => cur_cond%next
    enddo
    allocate(this%aux_vars_bc      (ncells_cond))
    allocate(this%boundary_flux    (ncells_cond))
    allocate(this%bnd_mass_exc     (ncells_cond))
    allocate(this%aux_vars_conn_bc (ncells_cond))

    do icond = 1,ncells_cond
       call this%aux_vars_bc      (icond)%Init()
       call this%aux_vars_conn_bc (icond)%Init()
    enddo
    this%boundary_flux(:) = 0.d0
    this%bnd_mass_exc(:) = 0.d0
    allocate(this%soe_auxvars_bc_offset(ncond))
    this%soe_auxvars_bc_offset(:) = -1

    ! Allocate memory and initialize aux vars: For source sink connections
    ncells_cond = 0
    ncond       = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       ncells_cond = ncells_cond + cur_cond%ncells
       ncond       = ncond + 1
       cur_cond => cur_cond%next
    enddo
    allocate(this%aux_vars_ss    (ncells_cond))
    allocate(this%ss_flux        (ncells_cond))

    this%ss_flux(:) = 0.d0

    do icond = 1,ncells_cond
       call this%aux_vars_ss(icond)%Init()
    enddo
    allocate(this%soe_auxvars_ss_offset(ncond))
    this%soe_auxvars_ss_offset(:) = -1

  end subroutine RichardsODEPressureAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureSetDensityType(this, density_type)
    !
    ! !DESCRIPTION:
    ! Set type of density formulation to auxiliary variables associated with
    ! internal, boundary, and source-sink conditions.
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : VAR_DENSITY_TYPE
    use EOSWaterMod               , only : DENSITY_CONSTANT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use EOSWaterMod               , only : DENSITY_IFC67
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    PetscInt                                 :: density_type
    !
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: sum_conn
    PetscInt                                 :: icond

    if (density_type /= DENSITY_CONSTANT .and. &
        density_type /= DENSITY_TGDPB01  .and. &
        density_type /= DENSITY_IFC67          &
        ) then
       write(iulog,*) 'Unknown value for VAR_DENSITY_TYPE'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! For internal connections
    do icond = 1,this%mesh%ncells_all
       this%aux_vars_in(icond)%density_type = density_type
    enddo

    ! For boundary conditions
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       do icond = 1,cur_cond%ncells
          sum_conn = sum_conn + 1
          this%aux_vars_bc(sum_conn)%density_type = density_type
       enddo
       cur_cond => cur_cond%next
    enddo

    ! For source sink conditions
    sum_conn = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       do icond = 1,cur_cond%ncells
          sum_conn = sum_conn + 1
          this%aux_vars_ss(sum_conn)%density_type = density_type
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureSetDensityType

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureNumConditions(this, cond_type, &
              cond_type_to_exclude, num_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of conditions
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    PetscInt                                 :: cond_type
    PetscInt                                 :: cond_type_to_exclude
    PetscInt, intent(out)                    :: num_conds
    !
    type(condition_type),pointer             :: cur_cond
    character(len=256)                       :: string

    ! Choose the condition type
    select case (cond_type)
    case (COND_BC)
       cur_cond => this%boundary_conditions%first
    case (COND_SS)
      cur_cond => this%source_sinks%first
    case default
       write(string,*) cond_type
       write(iulog,*) 'Unknown cond_type = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    num_conds = 0
    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= cond_type_to_exclude) then
          num_conds = num_conds + 1
       endif
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureNumConditions

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureNumCellsInConditions(this, cond_type, &
                cond_type_to_exclude, num_conds, ncells_for_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of conditions (eg. boundary condition or
    ! source-sink) and number of control volumes associated with each condition
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_NULL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    PetscInt, intent(in)                     :: cond_type
    PetscInt, intent(in)                     :: cond_type_to_exclude
    PetscInt, intent(out)                    :: num_conds
    PetscInt, intent(out), pointer           :: ncells_for_conds(:)
    !
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: ncells_cond
    PetscInt                                 :: icond
    character(len=256)                       :: string

    ! Find number of BCs
    call this%NumConditions(cond_type, COND_NULL, num_conds)

    if (num_conds == 0) then
       nullify(ncells_for_conds)
       return
    endif

    allocate(ncells_for_conds(num_conds))

    ! Choose the condition type
    select case (cond_type)
    case (COND_BC)
       cur_cond => this%boundary_conditions%first
    case (COND_SS)
      cur_cond => this%source_sinks%first
    case default
       write(string,*) cond_type
       write(iulog,*) 'Unknown cond_type = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    num_conds = 0
    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= cond_type_to_exclude) then
          num_conds = num_conds + 1
          ncells_for_conds(num_conds) = cur_cond%ncells
       endif
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureNumCellsInConditions

  !------------------------------------------------------------------------
  subroutine RichardsODEComputeResidual(this, X, F, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the residual equation for the discretized Richards equation
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    Vec                                      :: X
    Vec                                      :: F
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    PetscReal, pointer                         :: f_p(:)
    PetscReal, pointer                         :: accum_prev_p(:)

    call VecGetArrayF90(F, f_p, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(this%accum_prev, accum_prev_p, ierr); CHKERRQ(ierr)

    ! Computes the following:
    ! F = \left( \frac{(\phi * s * \rho)_i V_i}{dt} \right)^{k+1}
    call RichardsODEPressureAccum(this, f_p)

    ! F += -\left( \frac{(\phi * s * \rho)_i V_i}{dt} \right)^{k}
    f_p(:) = f_p(:) - accum_prev_p(:)

    ! F += \sum (\rho \mathbf{q})_{i,j} \cdot \mathbf{n}_{i,j} A_{i,j}  - Q_i
    call RichardsODEPressureDivergence(this, f_p)

    call VecRestoreArrayF90(this%accum_prev, accum_prev_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(F, f_p, ierr); CHKERRQ(ierr)

  end subroutine RichardsODEComputeResidual


  !------------------------------------------------------------------------
  subroutine RichardsODEComputeJacobian(this, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the jacobian matrix for the discretized Richards equation
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    Vec                                      :: X
    Mat                                      :: A
    Mat                                      :: B
    PetscErrorCode                           :: ierr

    ! Computes the following:
    ! \sum (\rho \mathbf{q})_{i,j} \cdot \mathbf{n}_{i,j} A_{i,j}  - Q_i
    call RichardsODEPressureDivergenceDeriv(this, B, ierr); CHKERRQ(ierr)

    ! \left( \frac{(\phi * s * \rho)_i V_i}{dt} \right)^{k+1}
    call RichardsODEPressureAccumDeriv(this, B, ierr); CHKERRQ(ierr)

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    if ( A /= B) then
       call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
       call MatAssemblyEnd(  A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    endif

  end subroutine RichardsODEComputeJacobian

  !------------------------------------------------------------------------
  subroutine RichardsODEComputeOffDiagJacobian(this, X_1, X_2, A, B, &
       id_of_other_goveq, list_id_of_other_goveq,        &
       ierr)
    !
    ! !DESCRIPTION:
    ! Computes the off-diagonal jacobian sub matrix associated with the
    ! coupling the given governing equation with another governing equation.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_RE
    use MultiPhysicsProbConstants, only : GE_THERM_SOIL_EBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    Vec                                      :: X_1
    Vec                                      :: X_2
    Mat                                      :: A
    Mat                                      :: B
    PetscInt                                 :: id_of_other_goveq
    PetscInt                                 :: list_id_of_other_goveq
    PetscErrorCode                           :: ierr
    !
    ! LOCAL VARIABLES
    character(len=256)                       :: string

    select case(id_of_other_goveq)
    case (GE_RE)
       call RichardsODEPressureJacOffDiag_BC(this, list_id_of_other_goveq, B, ierr)
    case (GE_THERM_SOIL_EBASED)
       call RichardsODEPressureJacOffDiag_Temp(this, list_id_of_other_goveq, B, ierr)
    case default
       write(string,*) id_of_other_goveq
       write(iulog,*) 'Unknown id_of_other_goveq = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    if ( A /= B) then
       call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
       call MatAssemblyEnd(  A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    endif

  end subroutine RichardsODEComputeOffDiagJacobian


  !------------------------------------------------------------------------
  subroutine RichardsODEPressureGetFromSOEAuxVarsIntrn(this, soe_avars, offset)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for internal nodes
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants     , only : VAR_PRESSURE
    use SystemOfEquationsVSFMAuxType  , only : sysofeqns_vsfm_auxvar_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type), intent(inout)       :: this
    type(sysofeqns_vsfm_auxvar_type), dimension(:), intent(in)    :: soe_avars
    PetscInt, intent(in) :: offset
    !
    ! LOCAL VARIABLES
    PetscInt                                                      :: iauxvar
    PetscInt                                                      :: nauxvar
    type(rich_ode_pres_auxvar_type), dimension(:), pointer        :: ge_avars

    ge_avars => this%aux_vars_in

    nauxvar = size(ge_avars)
    if( nauxvar > size(soe_avars) ) then
       write(iulog,*) 'size(ge_avars) > size(soe_avars)'
       write(iulog,*) 'size(ge_avars)  ',size(ge_avars)
       write(iulog,*) 'size(soe_avars) ',size(soe_avars)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do iauxvar = 1, nauxvar

       if (this%mesh%is_active(iauxvar)) then
          ! Copy temperature.
          ge_avars(iauxvar)%temperature =  &
               soe_avars(iauxvar+offset)%temperature

          ! Copy frac_liq_sat.
          ge_avars(iauxvar)%frac_liq_sat =  &
               soe_avars(iauxvar+offset)%frac_liq_sat

          ! Copy pressure.
          ge_avars(iauxvar)%pressure =  &
               soe_avars(iauxvar+offset)%pressure
       endif
    enddo

  end subroutine RichardsODEPressureGetFromSOEAuxVarsIntrn

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureSetFromSOEAuxVarsIntrn(this, var_type, ndata, data)
    !
    ! !DESCRIPTION:
    ! Sets values in GE auxiliary variables of internal nodes
    !
    ! !USES:
    use MultiPhysicsProbConstants     , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants     , only : VAR_PRESSURE
    use MultiPhysicsProbConstants     , only : VAR_FRAC_LIQ_SAT
    use SystemOfEquationsVSFMAuxType  , only : sysofeqns_vsfm_auxvar_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type), intent(inout)       :: this
    PetscInt, intent(in) :: var_type
    PetscInt, intent(in) :: ndata
    PetscReal, intent(in), pointer :: data(:)
    !
    ! LOCAL VARIABLES
    PetscInt                                                      :: iauxvar
    PetscInt                                                      :: nauxvar
    type(rich_ode_pres_auxvar_type), dimension(:), pointer        :: ge_avars

    ge_avars => this%aux_vars_in

    nauxvar = size(ge_avars)
    if( nauxvar /= ndata ) then
       write(iulog,*) 'size(ge_avars) /= ndata'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case (var_type)
    case (VAR_TEMPERATURE)
       do iauxvar = 1, nauxvar
          if (this%mesh%is_active(iauxvar)) then
             ge_avars(iauxvar)%temperature = data(iauxvar)
          end if
       enddo

       case (VAR_FRAC_LIQ_SAT)
          do iauxvar = 1, nauxvar
             if (this%mesh%is_active(iauxvar)) then
                ge_avars(iauxvar)%frac_liq_sat = data(iauxvar)
             end if
          enddo

       case (VAR_PRESSURE)
          do iauxvar = 1, nauxvar
             if (this%mesh%is_active(iauxvar)) then
                ge_avars(iauxvar)%pressure = data(iauxvar)
             end if
          enddo

       case default
          write(iulog,*) 'Unknown var_type = ', var_type
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

  end subroutine RichardsODEPressureSetFromSOEAuxVarsIntrn

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureGetFromSOEAuxVarsBC(this, soe_avars)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for bondary conditions
    !
    ! !USES:
    use ConditionType               , only : condition_type
    use ConnectionSetType           , only : connection_set_type
    use MultiPhysicsProbConstants   , only : COND_DIRICHLET
    use MultiPhysicsProbConstants   , only : COND_MASS_FLUX
    use MultiPhysicsProbConstants   , only : COND_MASS_RATE
    use MultiPhysicsProbConstants   , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants   , only : COND_SEEPAGE_BC
    use MultiPhysicsProbConstants   , only : VAR_BC_SS_CONDITION
    use SystemOfEquationsVSFMAuxType, only : sysofeqns_vsfm_auxvar_type
    use MultiPhysicsProbConstants   , only : COND_BC
    use MultiPhysicsProbConstants   , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) , intent(inout)            :: this
    type(sysofeqns_vsfm_auxvar_type)         , dimension(:), intent(in) :: soe_avars
    !
    ! LOCAL VARIABLES
    integer                                                             :: iauxvar
    integer                                                             :: iauxvar_off
    integer                                                             :: iconn
    integer                                                             :: auxVarCt_ge
    integer                                                             :: auxVarCt_soe
    integer                                                             :: condition_id
    integer                                                             :: sum_conn
    type(rich_ode_pres_auxvar_type)          , dimension(:), pointer    :: ge_avars
    type(condition_type)                     , pointer                  :: cur_cond
    type(connection_set_type)                , pointer                  :: cur_conn_set
    character(len=256)                                                  :: string
    PetscInt                                                            :: num_bc
    PetscInt                                                            :: icond
    PetscInt                                                            :: cond_itype_to_exclude
    PetscInt                                 , pointer                  :: ncells_for_bc(:)

    ge_avars => this%aux_vars_bc

    cond_itype_to_exclude = COND_DIRICHLET_FRM_OTR_GOVEQ
    call this%GetNCellsInCondsExcptCondItype(COND_BC, &
         cond_itype_to_exclude, num_bc, ncells_for_bc)

    auxVarCt_ge = 0
    do icond = 1, num_bc
       auxVarCt_ge = auxVarCt_ge + ncells_for_bc(icond)
    enddo
    if (associated(ncells_for_bc)) deallocate(ncells_for_bc)

    if (auxVarCt_ge == 0) return

    auxVarCt_soe = size(soe_avars)
    if( auxVarCt_ge > auxVarCt_soe ) then
       write(iulog,*) 'size(ge_avars) > size(soe_avars)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    condition_id = 0
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       condition_id = condition_id + 1

       ! Find first soe-auxvar corresponding to goveqn-auxvar.
       iauxvar_off = -1
       do iauxvar = 1, auxVarCt_soe
          if(  &
               soe_avars(iauxvar)%is_bc  &
               .and.  &
               soe_avars(iauxvar)%goveqn_id == this%rank_in_soe_list  &
               .and.  &
               soe_avars(iauxvar)%condition_id == condition_id  &
               ) then
             iauxvar_off = iauxvar - 1
             exit
          end if
       end do
       if (iauxvar_off < 0) then
          write(iulog,*) 'RichardsODEPressureGetFromSOEAuxVarsBC: iauxvar_off < 0'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       cur_conn_set => cur_cond%conn_set
       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          select case(cur_cond%itype)
          case (COND_DIRICHLET, COND_MASS_RATE, COND_MASS_FLUX)
             ge_avars(sum_conn)%condition_value =  &
                  soe_avars(iconn + iauxvar_off)%condition_value
          case (COND_DIRICHLET_FRM_OTR_GOVEQ)
             ! Do nothing
          case (COND_SEEPAGE_BC)
             ge_avars(sum_conn)%condition_value = &
                  soe_avars(iconn + iauxvar_off)%condition_value
          case default
             write(string,*) cur_cond%itype
             write(iulog,*) 'Unknown cur_cond%itype = ' // trim(string)
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       enddo

       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureGetFromSOEAuxVarsBC


  !------------------------------------------------------------------------
  subroutine RichardsODEPressureGetFromSOEAuxVarsSS(this, soe_avars)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for bondary conditions
    !
    ! !USES:
    use ConditionType                , only : condition_type
    use ConnectionSetType            , only : connection_set_type
    use MultiPhysicsProbConstants    , only : COND_MASS_RATE
    use MultiPhysicsProbConstants    , only : COND_DOWNREG_MASS_RATE_CAMPBELL
    use MultiPhysicsProbConstants    , only : COND_DOWNREG_MASS_RATE_FETCH2
    use MultiPhysicsProbConstants    , only : VAR_BC_SS_CONDITION
    use SystemOfEquationsVSFMAuxType , only : sysofeqns_vsfm_auxvar_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type), intent(inout)    :: this
    type(sysofeqns_vsfm_auxvar_type), dimension(:), intent(in) :: soe_avars
    !
    ! LOCAL VARIABLES
    integer                                                    :: iauxvar
    integer                                                    :: iauxvar_off
    integer                                                    :: iconn
    integer                                                    :: auxVarCt_ge, auxVarCt_soe
    integer                                                    :: condition_id, sum_conn
    integer                                                    :: cell_id_dn
    PetscReal                                                  :: var_value
    type(rich_ode_pres_auxvar_type), dimension(:), pointer     :: ge_avars
    type(condition_type), pointer                              :: cur_cond
    type(connection_set_type), pointer                         :: cur_conn_set
    character(len=256)                                         :: string

    ge_avars => this%aux_vars_ss
    auxVarCt_ge = size(ge_avars)

    auxVarCt_soe = size(soe_avars)
    if( auxVarCt_ge > auxVarCt_soe ) then
       write(iulog,*) 'size(ge_avars) > size(soe_avars)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    condition_id = 0
    sum_conn = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       condition_id = condition_id + 1

       ! Find first soe-auxvar corresponding to goveqn-auxvar.
       iauxvar_off = -1
       do iauxvar = 1, auxVarCt_soe
          if(  &
               soe_avars(iauxvar)%is_ss  &
               .and.  &
               soe_avars(iauxvar)%goveqn_id == this%rank_in_soe_list  &
               .and.  &
               soe_avars(iauxvar)%condition_id == condition_id  &
               ) then
             iauxvar_off = iauxvar - 1
             exit
          end if
       end do
       if (iauxvar_off < 0) then
          write(iulog,*) 'RichardsODEPressureGetFromSOEAuxVarsSS: iauxvar_off < 0'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       if (trim(cur_cond%name) /= 'Lateral_flux') then
          !
          ! Do not update values associated with lateral flux because:
          !  - SOE auxvars associated with lateral flux source-sink
          !    have zero values, AND
          !  - GE auxvars already have pre-computed values of lateral flux.
          !
          cur_conn_set => cur_cond%conn_set
          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             cell_id_dn = cur_conn_set%conn(iconn)%GetIDDn()
             if ( (.not. this%mesh%is_active(cell_id_dn))) cycle

             select case(cur_cond%itype)
             case (COND_MASS_RATE, COND_DOWNREG_MASS_RATE_CAMPBELL, COND_DOWNREG_MASS_RATE_FETCH2)
                var_value = soe_avars(iconn + iauxvar_off)%condition_value
                ge_avars(sum_conn)%condition_value = var_value
                cur_cond%value(iconn) = var_value
             case default
                write(string,*) cur_cond%itype
                write(iulog,*) 'Unknown cur_cond%itype = ' // trim(string)
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select
          enddo
       endif

       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureGetFromSOEAuxVarsSS


  !------------------------------------------------------------------------
  subroutine RichardsODEPressureGetDataFromSOEAuxVar(this, soe_avar_type, soe_avars, &
       offset)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    !
    use MultiPhysicsProbConstants   , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants   , only : AUXVAR_BC
    use MultiPhysicsProbConstants   , only : AUXVAR_SS
    use SystemOfEquationsVSFMAuxType, only : sysofeqns_vsfm_auxvar_type
    use ConditionType               , only : condition_type
    use ConnectionSetType           , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type), intent(inout)        :: this
    PetscInt, intent(in)                                           :: soe_avar_type
    type (sysofeqns_vsfm_auxvar_type), dimension(:), intent(in)    :: soe_avars
    PetscInt, intent(in), optional :: offset
    !
    ! !LOCAL VARIABLES
    type (rich_ode_pres_auxvar_type), pointer       :: ge_avars(:)
    PetscInt                                        :: iauxvar_off

    select case(soe_avar_type)
    case (AUXVAR_INTERNAL)
       if (present(offset)) then
          iauxvar_off = offset
       else
          iauxvar_off = 0
       endif
       call RichardsODEPressureGetFromSOEAuxVarsIntrn(this, soe_avars, iauxvar_off)
    case (AUXVAR_BC)
       call RichardsODEPressureGetFromSOEAuxVarsBC(this, soe_avars)
    case (AUXVAR_SS)
       call RichardsODEPressureGetFromSOEAuxVarsSS(this, soe_avars)
    case default
       write(iulog,*) 'RichardsODEGetDataFromSOEAuxVar: soe_avar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine RichardsODEPressureGetDataFromSOEAuxVar


  !------------------------------------------------------------------------
  subroutine RichardsODEPressureSetDataInSOEAuxVar(this, soe_avar_type, soe_avars, &
       offset, num_avars_set)
    !
    ! !DESCRIPTION:
    ! Copies data GE auxiliary variable into SoE auxiliary variable.
    ! This is done as a part of post solve.
    !
    use MultiPhysicsProbConstants    , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants    , only : AUXVAR_BC
    use MultiPhysicsProbConstants    , only : AUXVAR_SS
    use MultiPhysicsProbConstants    , only : AUXVAR_CONN_INTERNAL
    use MultiPhysicsProbConstants    , only : VAR_LIQ_SAT
    use MultiPhysicsProbConstants    , only : VAR_MASS
    use MultiPhysicsProbConstants    , only : VAR_SOIL_MATRIX_POT
    use MultiPhysicsProbConstants    , only : FMWH2O
    use MultiPhysicsProbConstants    , only : PRESSURE_REF
    use MultiPhysicsProbConstants    , only : GRAVITY_CONSTANT
    use MultiPhysicsProbConstants    , only : COND_MASS_RATE
    use MultiPhysicsProbConstants    , only : CONN_HORIZONTAL
    use MultiPhysicsProbConstants    , only : COND_BC
    use MultiPhysicsProbConstants    , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants    , only : COND_DOWNREG_MASS_RATE_CAMPBELL
    use MultiPhysicsProbConstants    , only : COND_DOWNREG_MASS_RATE_FETCH2
    use SystemOfEquationsVSFMAuxType , only : sysofeqns_vsfm_auxvar_type
    use ConditionType                , only : condition_type
    use ConnectionSetType            , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type)                       :: this
    PetscInt                                                       :: soe_avar_type
    type (sysofeqns_vsfm_auxvar_type), dimension(:), intent(inout) :: soe_avars
    PetscInt                                                       :: offset
    PetscInt, intent(out)                                          :: num_avars_set
    !
    ! !LOCAL VARIABLES
    type (rich_ode_pres_auxvar_type), pointer                      :: ge_avars(:)
    PetscInt                                                       :: iauxvar
    PetscInt                                                       :: iauxvar_off
    PetscInt                                                       :: iconn
    PetscInt                                                       :: sum_conn
    PetscInt                                                       :: auxVarCt_ge
    PetscInt                                                       :: auxVarCt_soe
    PetscInt                                                       :: condition_id
    PetscReal                                                      :: mass
    PetscReal                                                      :: smp
    PetscReal                                                      :: Pa_to_Meters
    PetscInt                                                       :: cell_id_up
    PetscInt                                                       :: cell_id_dn
    PetscInt                                                       :: num_bc, icond
    PetscInt                                                       :: cond_itype_to_exclude
    PetscInt, pointer                                              :: ncells_for_bc(:)
    character(len=256)                                             :: string
    type(condition_type), pointer                                  :: cur_cond
    type(connection_set_type), pointer                             :: cur_conn_set

    iauxvar_off = offset

    select case(soe_avar_type)
    case (AUXVAR_INTERNAL)
       ge_avars => this%aux_vars_in

       if ( size(ge_avars) > size(soe_avars)) then
          write(iulog,*) 'size(ge_avars) > size(soe_avars)'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       ! Interior cells
       cur_conn_set => this%mesh%intrn_conn_set_list%first
       sum_conn = 0
       do
          if (.not.associated(cur_conn_set)) exit

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             cell_id_up = cur_conn_set%conn(sum_conn)%GetIDUp()
             cell_id_dn = cur_conn_set%conn(sum_conn)%GetIDDn()

             if (cur_conn_set%conn(sum_conn)%GetType() == CONN_HORIZONTAL) then

                this%lat_mass_exc(cell_id_up) = this%lat_mass_exc(cell_id_up) + &
                     this%internal_flux(sum_conn)*this%dtime

                this%lat_mass_exc(cell_id_dn) = this%lat_mass_exc(cell_id_dn) - &
                     this%internal_flux(sum_conn)*this%dtime
             endif

          enddo
          cur_conn_set => cur_conn_set%next
       enddo

       do iauxvar = 1, size(ge_avars)
          if (this%mesh%is_active(iauxvar)) then
             soe_avars(iauxvar+iauxvar_off)%liq_sat =  &
                  ge_avars(iauxvar)%sat

             mass =  &
                  this%aux_vars_in(iauxvar)%por*        & ! [-]
                  this%aux_vars_in(iauxvar)%den*FMWH2O* & ! [kg m^{-3}]
                  this%aux_vars_in(iauxvar)%sat*        & ! [-]
                  this%mesh%vol(iauxvar)                  ! [m^3]

             soe_avars(iauxvar+iauxvar_off)%mass = mass
             soe_avars(iauxvar+iauxvar_off)%lateral_mass_exchanged = &
                  this%lat_mass_exc(iauxvar)

             Pa_to_Meters = this%aux_vars_in(iauxvar)%den*FMWH2O*GRAVITY_CONSTANT
             smp          = (this%aux_vars_in(iauxvar)%pressure - &
                  PRESSURE_REF)/Pa_to_Meters

             soe_avars(iauxvar+iauxvar_off)%soil_matrix_pot = smp

          endif
       enddo

       num_avars_set = size(ge_avars)

    case (AUXVAR_CONN_INTERNAL)
       ge_avars => this%aux_vars_in

       ! Interior cells
       cur_conn_set => this%mesh%intrn_conn_set_list%first
       sum_conn = 0
       do
          if (.not.associated(cur_conn_set)) exit

          if ( sum_conn + cur_conn_set%num_connections + iauxvar_off > size(soe_avars)) then
             write(iulog,*) 'size(soe_avars) is smaller than number of internal connections'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             soe_avars(iconn + iauxvar_off)%mass_flux = this%internal_flux(sum_conn)

          enddo
          cur_conn_set => cur_conn_set%next
       enddo

       num_avars_set = sum_conn

    case (AUXVAR_SS)

       condition_id = 0
       sum_conn = 0

       auxVarCt_soe = size(soe_avars)

       cur_cond => this%source_sinks%first
       do
          if (.not.associated(cur_cond)) exit
          condition_id = condition_id + 1

          ! Find first soe-auxvar corresponding to goveqn-auxvar.
          iauxvar_off = -1
          do iauxvar = 1, auxVarCt_soe
             if(  &
                  soe_avars(iauxvar)%is_ss  &
                  .and.  &
                  soe_avars(iauxvar)%goveqn_id == this%rank_in_soe_list  &
                  .and.  &
                  soe_avars(iauxvar)%condition_id == condition_id  &
                  ) then
                iauxvar_off = iauxvar - 1
                exit
             end if
          end do
          if (iauxvar_off < 0) then
             write(iulog,*) 'RichardsODEPressureSetDataInSOEAuxVar: iauxvar_off < 0'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if

          if (trim(cur_cond%name) == 'Lateral_flux') then

             cur_conn_set => cur_cond%conn_set
             do iconn = 1, cur_conn_set%num_connections
                sum_conn = sum_conn + 1
                select case(cur_cond%itype)
                case (COND_MASS_RATE)
                   soe_avars(iconn + iauxvar_off)%condition_value = cur_cond%value(iconn)
                   soe_avars(iconn + iauxvar_off)%mass_flux       = this%ss_flux(sum_conn)

                case default
                   write(string,*) cur_cond%itype
                   write(iulog,*) 'Unknown cur_cond%itype = ' // trim(string)
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                end select
             enddo

          else

             cur_conn_set => cur_cond%conn_set
             do iconn = 1, cur_conn_set%num_connections

                sum_conn = sum_conn + 1
                select case(cur_cond%itype)
                case (COND_MASS_RATE)
                   soe_avars(iconn + iauxvar_off)%mass_flux = this%ss_flux(sum_conn)

                case (COND_DOWNREG_MASS_RATE_CAMPBELL, COND_DOWNREG_MASS_RATE_FETCH2)
                   soe_avars(iconn + iauxvar_off)%mass_flux = this%ss_flux(sum_conn)

                case default
                   write(string,*) cur_cond%itype
                   write(iulog,*) 'Unknown cur_cond%itype = ' // trim(string)
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                end select
             enddo

          endif

          cur_cond => cur_cond%next
       enddo

       num_avars_set = sum_conn

    case (AUXVAR_BC)

       cond_itype_to_exclude = COND_DIRICHLET_FRM_OTR_GOVEQ
       call this%GetNCellsInCondsExcptCondItype(COND_BC, &
            cond_itype_to_exclude, num_bc, ncells_for_bc)

       auxVarCt_ge = 0
       do icond = 1, num_bc
          auxVarCt_ge = auxVarCt_ge + ncells_for_bc(icond)
       enddo
       if (associated(ncells_for_bc)) deallocate(ncells_for_bc)

       if (auxVarCt_ge == 0) return

       ! Boundary condition cells
       cur_cond => this%boundary_conditions%first
       sum_conn = 0
       do
          if (.not.associated(cur_cond)) exit
          cur_conn_set => cur_cond%conn_set

          if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) cycle

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             this%bnd_mass_exc(sum_conn) = this%bnd_mass_exc(sum_conn) + &
                  this%boundary_flux(sum_conn)*this%dtime

             soe_avars(iconn + iauxvar_off)%boundary_mass_exchanged = &
                  this%bnd_mass_exc(sum_conn)

             soe_avars(iconn + iauxvar_off)%mass_flux = this%boundary_flux(sum_conn)
          enddo
          cur_cond => cur_cond%next
       enddo

       num_avars_set = sum_conn

    case default
       write(iulog,*) 'RichardsODESetDataInSOEAuxVar: soe_avar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine RichardsODEPressureSetDataInSOEAuxVar

  !------------------------------------------------------------------------

  subroutine RichardsODEPressureUpdateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Updates all auxiliary variables
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this

    call this%UpdateAuxVarsIntrn()
    call this%UpdateAuxVarsBC()
    call this%UpdateAuxVarsSS()

  end subroutine RichardsODEPressureUpdateAuxVars

  !------------------------------------------------------------------------

  subroutine RichardsODEPressureUpdateAuxVarsIntrn(this)
    !
    ! !DESCRIPTION:
    ! Updates auxiliary variable associated with internal control volumes
    !
    use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    type(connection_set_type), pointer       :: cur_conn_set
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: ghosted_id
    PetscInt                                 :: sum_conn
    PetscInt                                 :: iconn

    ! Update aux vars for internal cells
    do ghosted_id = 1, this%mesh%ncells_all
       if (this%mesh%is_active(ghosted_id)) then
          call this%aux_vars_in(ghosted_id)%AuxVarCompute()
       end if
    enddo

    ! Interior cells
    cur_conn_set => this%mesh%intrn_conn_set_list%first
    sum_conn = 0
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          select case(this%aux_vars_conn_in(sum_conn)%flux_type)
          case (DARCY_FLUX_TYPE)
          case (CONDUCTANCE_FLUX_TYPE)

             ! Set values for 'connection' auxvars
             ghosted_id = cur_conn_set%conn(iconn)%GetIDUp()

             this%aux_vars_conn_in(sum_conn)%pressure_up = &
                  this%aux_vars_in(ghosted_id)%pressure

             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

             this%aux_vars_conn_in(sum_conn)%pressure_dn = &
                  this%aux_vars_in(ghosted_id)%pressure

             call this%aux_vars_conn_in(sum_conn)%AuxVarCompute()
          end select
       end do
       cur_conn_set => cur_conn_set%next
    end do

  end subroutine RichardsODEPressureUpdateAuxVarsIntrn

  !------------------------------------------------------------------------

  subroutine RichardsODEPressureUpdateAuxVarsBC(this)
    !
    ! !DESCRIPTION:
    ! Updates auxiliary variable associated with boundary condition
    !
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_MASS_FLUX
    use MultiPhysicsProbConstants , only : COND_MASS_RATE
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_SEEPAGE_BC
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: ghosted_id
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscReal                                :: temperature
    type(condition_type),pointer             :: cur_cond
    type(connection_set_type), pointer       :: cur_conn_set
    character(len=256)                       :: string

    ! Update aux vars for boundary cells
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          select case(cur_cond%itype)
          case (COND_DIRICHLET, COND_SEEPAGE_BC)
             this%aux_vars_bc(sum_conn)%pressure = &
                  this%aux_vars_bc(sum_conn)%condition_value
          case (COND_MASS_RATE, COND_MASS_FLUX)
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()
             this%aux_vars_bc(sum_conn)%pressure =  &
                this%aux_vars_in(ghosted_id)%pressure
          case (COND_DIRICHLET_FRM_OTR_GOVEQ)
             ! Do nothing
          case default
             write(string,*) cur_cond%itype
             write(iulog,*) 'Unknown cur_cond%itype = ' // trim(string)
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select

          ! GB: Get rid of setting values for Temperature
          temperature = 273.15d0 + 25.d0
          this%aux_vars_bc(sum_conn)%temperature = temperature

          call this%aux_vars_bc(sum_conn)%AuxVarCompute()

          ! Set values for 'connection' auxvars
          this%aux_vars_conn_bc(sum_conn)%pressure_up = &
               this%aux_vars_bc(sum_conn)%pressure

          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          this%aux_vars_conn_bc(sum_conn)%pressure_dn = &
               this%aux_vars_in(ghosted_id)%pressure

          call this%aux_vars_conn_bc(sum_conn)%AuxVarCompute()

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureUpdateAuxVarsBC

  !------------------------------------------------------------------------

  subroutine RichardsODEPressureUpdateAuxVarsSS(this)
    !
    ! !DESCRIPTION:
    ! Updates auxiliary variable associated with soure sink
    !
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: ghosted_id
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    type(condition_type),pointer             :: cur_cond
    type(connection_set_type), pointer       :: cur_conn_set

    return

    ! Update aux vars for source/sink cells
    sum_conn = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()
          this%aux_vars_ss(sum_conn)%pressure =  &
              this%aux_vars_in(ghosted_id)%pressure
          call this%aux_vars_ss(sum_conn)%AuxVarCompute()
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureUpdateAuxVarsSS


  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAccum(this, ff)
    !
    ! !DESCRIPTION:
    ! Computes the accumulation term
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type), intent(inout)    :: this
    PetscReal, dimension(:), intent(out)                       :: ff

    ! !LOCAL VARIABLES
    PetscInt                                                  :: cell_id
    PetscReal                                                 :: dtInv

    dtInv = 1.d0 / this%dtime

    ff(:) = 0.d0

    ! Interior cells.
    do cell_id = 1, this%mesh%ncells_local

       if (this%mesh%is_active(cell_id)) then
          ff(cell_id) = this%aux_vars_in(cell_id)%por   * &
                        this%aux_vars_in(cell_id)%den   * &
                        this%aux_vars_in(cell_id)%sat   * &
                        this%mesh%vol(cell_id)          * &
                        dtInv
       end if
    enddo

  end subroutine RichardsODEPressureAccum


  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAccumDeriv(this, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes derivative of accumulation term
    ! \frac{\partial \xi}{\partial dt} dV
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    Mat                                      :: B
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: cell_id
    PetscInt                                 :: row
    PetscInt                                 :: col
    PetscReal                                :: por
    PetscReal                                :: den
    PetscReal                                :: sat
    PetscReal                                :: dpor_dp
    PetscReal                                :: dden_dp
    PetscReal                                :: dsat_dp
    PetscReal                                :: derivative
    PetscReal                                :: dtInv

    dtInv = 1.d0 / this%dtime

    ! Interior cells
    do cell_id = 1, this%mesh%ncells_local

       if (this%mesh%is_active(cell_id) ) then
          por     = this%aux_vars_in(cell_id)%por
          den     = this%aux_vars_in(cell_id)%den
          sat     = this%aux_vars_in(cell_id)%sat

          dpor_dp = this%aux_vars_in(cell_id)%dpor_dp
          dden_dp = this%aux_vars_in(cell_id)%dden_dp
          dsat_dp = this%aux_vars_in(cell_id)%dsat_dp

          derivative = (dpor_dp*den    *sat     + &
                        por    *dden_dp*sat     + &
                        por    *den    *dsat_dp    )*this%mesh%vol(cell_id) * dtInv

       else
          derivative = 1.d0
       endif

       row = cell_id - 1
       col = cell_id - 1

       call MatSetValuesLocal(B, 1, row, 1, col, derivative, ADD_VALUES, ierr); CHKERRQ(ierr)

    enddo

  end subroutine RichardsODEPressureAccumDeriv

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureDivergence(this, ff)
    !
    ! !DESCRIPTION:
    ! Computes the divergence associated with internal and boundary
    ! conditions for residual equation. Also computes the contribution of
    ! source-sink to residual equation.
    !
    ! \int \nabla \cdot (\rho \mathbf{q}) dV
    ! = \int (\rho \mathbf{q}) \cdot \mathbf{n} dA
    !
    ! !USES:
    use RichardsMod                 , only : RichardsFlux
    use RichardsMod                 , only : RichardsFluxConductanceModel
    use ConditionType               , only : condition_type
    use ConnectionSetType           , only : connection_set_type
    use MultiPhysicsProbConstants   , only : COND_NULL
    use MultiPhysicsProbConstants   , only : FMWH2O
    use MultiPhysicsProbConstants   , only : COND_MASS_RATE
    use MultiPhysicsProbConstants   , only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants   , only : CONDUCTANCE_FLUX_TYPE
    use MultiPhysicsProbConstants   , only : PRESSURE_REF
    use MultiPhysicsProbConstants   , only : COND_DOWNREG_MASS_RATE_CAMPBELL
    use MultiPhysicsProbConstants   , only : COND_DOWNREG_MASS_RATE_FETCH2
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    PetscReal, dimension(:), intent(inout)   :: ff
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id_dn
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id
    PetscReal                                :: flux
    PetscReal                                :: dP
    PetscReal                                :: factor
    PetscReal                                :: n
    PetscReal                                :: Pc
    PetscReal                                :: dummy_var1
    PetscReal                                :: dummy_var2
    PetscBool                                :: compute_deriv
    PetscBool                                :: internal_conn
    PetscInt                                 :: cond_type
    type(condition_type),pointer             :: cur_cond
    type(connection_set_type), pointer       :: cur_conn_set

    compute_deriv = PETSC_FALSE

    ! Interior cells
    cur_conn_set => this%mesh%intrn_conn_set_list%first
    sum_conn = 0
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id_up = cur_conn_set%conn(sum_conn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(sum_conn)%GetIDDn()

          internal_conn = PETSC_TRUE
          cond_type     = COND_NULL

          if ( (.not. this%mesh%is_active(cell_id_up)) .or. &
               (.not. this%mesh%is_active(cell_id_dn)) ) cycle

          select case(this%aux_vars_conn_in(sum_conn)%flux_type)

          case (DARCY_FLUX_TYPE)

             call RichardsFlux(                 &
                  this%aux_vars_in(cell_id_up), &
                  this%aux_vars_in(cell_id_dn), &
                  cur_conn_set%conn(iconn),     &
                  compute_deriv,                &
                  internal_conn,                &
                  cond_type,                    &
                  flux,                         &
                  dummy_var1,                   &
                  dummy_var2                    &
                  )

          case (CONDUCTANCE_FLUX_TYPE)

             call RichardsFluxConductanceModel(    &
                  this%aux_vars_in(cell_id_up),    &
                  this%aux_vars_in(cell_id_dn),    &
                  this%aux_vars_conn_in(sum_conn), &
                  cur_conn_set%conn(iconn),        &
                  compute_deriv,                   &
                  internal_conn,                   &
                  cond_type,                       &
                  flux,                            &
                  dummy_var1,                      &
                  dummy_var2                       &
                  )

          case default
             write(iulog,*) 'Unknown flux_type ', this%aux_vars_conn_in(sum_conn)%flux_type
             call endrun(msg=errMsg(__FILE__,__LINE__))
          end select

          ff(cell_id_up) = ff(cell_id_up) - flux;
          ff(cell_id_dn) = ff(cell_id_dn) + flux;

          this%internal_flux(sum_conn) = flux*FMWH2O
       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    ! Boundary cells
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id = cur_conn_set%conn(iconn)%GetIDDn()

          internal_conn = PETSC_FALSE
          cond_type     = cur_cond%itype

          if ( (.not. this%mesh%is_active(cell_id))) cycle

          if (.not.cur_cond%swap_order) then

             select case(this%aux_vars_conn_bc(sum_conn)%flux_type)

             case (DARCY_FLUX_TYPE)

                call RichardsFlux(                    &
                     this%aux_vars_bc(sum_conn),      &
                     this%aux_vars_in(cell_id ),      &
                     cur_conn_set%conn(iconn),        &
                     compute_deriv,                   &
                     internal_conn,                   &
                     cond_type,                       &
                     flux,                            &
                     dummy_var1,                      &
                     dummy_var2                       &
                     )

             case (CONDUCTANCE_FLUX_TYPE)
                call RichardsFluxConductanceModel(    &
                     this%aux_vars_bc(sum_conn ),     &
                     this%aux_vars_in(cell_id  ),     &
                     this%aux_vars_conn_bc(sum_conn), &
                     cur_conn_set%conn(iconn),        &
                     compute_deriv,                   &
                     internal_conn,                   &
                     cond_type,                       &
                     flux,                            &
                     dummy_var1,                      &
                     dummy_var2                       &
                     )

             case default
                write(iulog,*) 'Unknown flux_type ', this%aux_vars_conn_in(sum_conn)%flux_type
                call endrun(msg=errMsg(__FILE__,__LINE__))
             end select

          else

             select case(this%aux_vars_conn_bc(sum_conn)%flux_type)

             case (DARCY_FLUX_TYPE)

                call RichardsFlux(                    &
                     this%aux_vars_in(cell_id ),      &
                     this%aux_vars_bc(sum_conn),      &
                     cur_conn_set%conn(iconn),        &
                     compute_deriv,                   &
                     internal_conn,                   &
                     cond_type,                       &
                     flux,                            &
                     dummy_var1,                      &
                     dummy_var2                       &
                     )

             case (CONDUCTANCE_FLUX_TYPE)
                call RichardsFluxConductanceModel(    &
                     this%aux_vars_in(cell_id ),      &
                     this%aux_vars_bc(sum_conn ),     &
                     this%aux_vars_conn_bc(sum_conn), &
                     cur_conn_set%conn(iconn),        &
                     compute_deriv,                   &
                     internal_conn,                   &
                     cond_type,                       &
                     flux,                            &
                     dummy_var1,                      &
                     dummy_var2                       &
                     )

             case default
                write(iulog,*) 'Unknown flux_type ', this%aux_vars_conn_in(sum_conn)%flux_type
                call endrun(msg=errMsg(__FILE__,__LINE__))
             end select

             flux = -flux
          endif

          ff(cell_id) = ff(cell_id) + flux;

          this%boundary_flux(sum_conn) = flux*FMWH2O

       enddo
       cur_cond => cur_cond%next
    enddo

    ! Source-sink cells
    sum_conn = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          cell_id = cur_conn_set%conn(iconn)%GetIDDn()
          sum_conn = sum_conn + 1

          if ( (.not. this%mesh%is_active(cell_id))) cycle

          select case(cur_cond%itype)
          case (COND_MASS_RATE)
             ff(cell_id) = ff(cell_id) - cur_cond%value(iconn)/FMWH2O

             this%ss_flux(sum_conn) = cur_cond%value(iconn)

          case (COND_DOWNREG_MASS_RATE_CAMPBELL)
             dP          = this%aux_vars_in(cell_id)%pressure - PRESSURE_REF
             Pc          = this%aux_vars_ss(sum_conn)%pot_mass_sink_pressure
             n           = this%aux_vars_ss(sum_conn)%pot_mass_sink_exponent

             if (dP <= 0.d0) then
                factor = 1.d0 + (dP/Pc)**n
             else
                factor = 1.d0
             endif

             ff(cell_id) = ff(cell_id) - cur_cond%value(iconn)/factor/FMWH2O
             this%ss_flux(sum_conn) = cur_cond%value(iconn)/factor

          case (COND_DOWNREG_MASS_RATE_FETCH2)
             dP          = this%aux_vars_in(cell_id)%pressure - PRESSURE_REF
             Pc          = this%aux_vars_ss(sum_conn)%pot_mass_sink_pressure
             n           = this%aux_vars_ss(sum_conn)%pot_mass_sink_exponent

             if (dP <= 0.d0) then
                factor = exp( -(dP/Pc)**n )
             else
                factor = 1.d0
             endif

             ff(cell_id) = ff(cell_id) - cur_cond%value(iconn)*factor/FMWH2O
             this%ss_flux(sum_conn) = cur_cond%value(iconn)*factor

          case default
            write(iulog,*)'RichardsODEPressureDivergence: Unknown cond_type in SS.'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end select
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureDivergence

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureDivergenceDeriv(this, B, ierr)
    !
    ! !DESCRIPTION:
    !
    ! \int \nabla \cdot (\rho \mathbf{q}) dV
    ! = \int (\rho \mathbf{q}) \cdot \mathbf{n} dA
    !
    ! !USES:
    use RichardsMod               , only : RichardsFlux
    use RichardsMod               , only : RichardsFluxConductanceModel
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_NULL
    use MultiPhysicsProbConstants , only : COND_MASS_RATE
    use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
    use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_CAMPBELL
    use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_FETCH2
    use MultiPhysicsProbConstants , only : FMWH2O
    use MultiPhysicsProbConstants , only : PRESSURE_REF
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    Mat                                      :: B
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id_dn
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id
    PetscInt                                 :: row
    PetscInt                                 :: col
    PetscReal                                :: dummy_var
    PetscReal                                :: Jup
    PetscReal                                :: Jdn
    PetscReal                                :: val
    PetscReal                                :: factor
    PetscReal                                :: n
    PetscReal                                :: Pc
    PetscReal                                :: dP
    PetscBool                                :: compute_deriv
    PetscBool                                :: internal_conn
    PetscInt                                 :: cond_type
    type(condition_type),pointer             :: cur_cond
    type(connection_set_type), pointer       :: cur_conn_set

    compute_deriv = PETSC_TRUE

    ! Interior cells
    cur_conn_set => this%mesh%intrn_conn_set_list%first
    sum_conn = 0
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id_up = cur_conn_set%conn(sum_conn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(sum_conn)%GetIDDn()

          internal_conn = PETSC_TRUE
          cond_type     = COND_NULL

          if ( (.not. this%mesh%is_active(cell_id_up)) .or. &
               (.not. this%mesh%is_active(cell_id_dn)) ) cycle

          select case(this%aux_vars_conn_in(sum_conn)%flux_type)

          case (DARCY_FLUX_TYPE)

             call RichardsFlux(                    &
                  this%aux_vars_in(cell_id_up),    &
                  this%aux_vars_in(cell_id_dn),    &
                  cur_conn_set%conn(iconn),        &
                  compute_deriv,                   &
                  internal_conn,                   &
                  cond_type,                       &
                  dummy_var,                       &
                  Jup,                             &
                  Jdn                              &
                  )

          case (CONDUCTANCE_FLUX_TYPE)
             call RichardsFluxConductanceModel(    &
                  this%aux_vars_in(cell_id_up),    &
                  this%aux_vars_in(cell_id_dn),    &
                  this%aux_vars_conn_in(sum_conn), &
                  cur_conn_set%conn(iconn),        &
                  compute_deriv,                   &
                  internal_conn,                   &
                  cond_type,                       &
                  dummy_var,                       &
                  Jup,                             &
                  Jdn                              &
                  )

          case default
             write(iulog,*) 'Unknown flux_type ', this%aux_vars_conn_in(sum_conn)%flux_type
             call endrun(msg=errMsg(__FILE__,__LINE__))
          end select

          row = cell_id_up - 1
          col = cell_id_up - 1
          val = Jup
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_up - 1
          col = cell_id_dn - 1
          val = Jdn
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_dn - 1
          col = cell_id_up - 1
          val = -Jup
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_dn - 1
          col = cell_id_dn - 1
          val = -Jdn
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    ! Boundary cells
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id = cur_conn_set%conn(iconn)%GetIDDn()

          if ( (.not. this%mesh%is_active(cell_id))) cycle

          internal_conn = PETSC_FALSE
          cond_type     = cur_cond%itype

          if (.not.cur_cond%swap_order) then
             select case(this%aux_vars_conn_bc(sum_conn)%flux_type)

             case (DARCY_FLUX_TYPE)

                call RichardsFlux(                    &
                     this%aux_vars_bc(sum_conn),      &
                     this%aux_vars_in(cell_id),       &
                     cur_conn_set%conn(iconn),        &
                     compute_deriv,                   &
                     internal_conn,                   &
                     cond_type,                       &
                     dummy_var,                       &
                     Jup,                             &
                     Jdn                              &
                     )

             case (CONDUCTANCE_FLUX_TYPE)
                call RichardsFluxConductanceModel(    &
                     this%aux_vars_bc(sum_conn ),     &
                     this%aux_vars_in(cell_id ),      &
                     this%aux_vars_conn_bc(sum_conn), &
                     cur_conn_set%conn(iconn),        &
                     compute_deriv,                   &
                     internal_conn,                   &
                     cond_type,                       &
                     dummy_var,                       &
                     Jup,                             &
                     Jdn                              &
                     )

             case default
                write(iulog,*) 'Unknown flux_type ', this%aux_vars_conn_in(sum_conn)%flux_type
                call endrun(msg=errMsg(__FILE__,__LINE__))
             end select

             val = -Jdn
          else
             select case(this%aux_vars_conn_bc(sum_conn)%flux_type)

             case (DARCY_FLUX_TYPE)

                call RichardsFlux(                    &
                     this%aux_vars_in(cell_id ),      &
                     this%aux_vars_bc(sum_conn),      &
                     cur_conn_set%conn(iconn),        &
                     compute_deriv,                   &
                     internal_conn,                   &
                     cond_type,                       &
                     dummy_var,                       &
                     Jup,                             &
                     Jdn                              &
                     )

             case (CONDUCTANCE_FLUX_TYPE)
                call RichardsFluxConductanceModel(    &
                     this%aux_vars_in(cell_id ),      &
                     this%aux_vars_bc(sum_conn ),     &
                     this%aux_vars_conn_bc(sum_conn), &
                     cur_conn_set%conn(iconn),        &
                     compute_deriv,                   &
                     internal_conn,                   &
                     cond_type,                       &
                     dummy_var,                       &
                     Jup,                             &
                     Jdn                              &
                     )

             case default
                write(iulog,*) 'Unknown flux_type ', this%aux_vars_conn_in(sum_conn)%flux_type
                call endrun(msg=errMsg(__FILE__,__LINE__))
             end select

             val = Jup
          endif

          row = cell_id - 1
          col = cell_id - 1
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

       enddo
       cur_cond => cur_cond%next
    enddo

    ! Source-sink cells
    sum_conn = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id = cur_conn_set%conn(iconn)%GetIDDn()

          if ( (.not. this%mesh%is_active(cell_id))) cycle

          select case(cur_cond%itype)
          case (COND_MASS_RATE)

          case (COND_DOWNREG_MASS_RATE_CAMPBELL)
             dP          = this%aux_vars_in(cell_id)%pressure - PRESSURE_REF
             Pc          = this%aux_vars_ss(sum_conn)%pot_mass_sink_pressure
             n           = this%aux_vars_ss(sum_conn)%pot_mass_sink_exponent

             if (dP <= 0.d0) then
                cell_id = cur_conn_set%conn(iconn)%GetIDDn()

                factor = 1.d0 + (dP/Pc)**n
                val = (cur_cond%value(iconn)/FMWH2O) * ( n * (dP/Pc)**n) / ( dP * ( factor**2.d0))

                row = cell_id - 1
                col = cell_id - 1
                call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)
             endif

          case (COND_DOWNREG_MASS_RATE_FETCH2)
             dP          = this%aux_vars_in(cell_id)%pressure - PRESSURE_REF
             Pc          = this%aux_vars_ss(sum_conn)%pot_mass_sink_pressure
             n           = this%aux_vars_ss(sum_conn)%pot_mass_sink_exponent

             if (dP <= 0.d0) then
                cell_id = cur_conn_set%conn(iconn)%GetIDDn()

                factor = exp( -(dP/Pc)**n )
                val = (cur_cond%value(iconn)/FMWH2O) * (n * (dP/Pc)**n) * factor/dP

                row = cell_id - 1
                col = cell_id - 1
                call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)
             endif

          case default
            write(iulog,*)'RichardsDAEPressureDivergenceDeriv: Unknown cond_type in SS.'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end select

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureDivergenceDeriv

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureJacOffDiag_BC(this, list_id_of_other_goveq, &
                                              B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the jacobian submatrix associated with the coupling the given
    ! governing equation with another governing equation.
    !
    ! !USES:
    use RichardsMod               , only : RichardsFlux
    use RichardsMod               , only : RichardsFluxConductanceModel
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    PetscInt                                 :: list_id_of_other_goveq
    Mat                                      :: B
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: ieqn
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id_dn
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id
    PetscInt                                 :: row
    PetscInt                                 :: col
    PetscReal                                :: dummy_var
    PetscReal                                :: Jup
    PetscReal                                :: Jdn
    PetscReal                                :: val
    PetscBool                                :: compute_deriv
    PetscBool                                :: internal_conn
    PetscBool                                :: cur_cond_used
    PetscInt                                 :: cond_type
    type(condition_type),pointer             :: cur_cond
    type(connection_set_type), pointer       :: cur_conn_set

    compute_deriv = PETSC_TRUE

    ! Boundary cells
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       cur_cond_used = PETSC_FALSE

       if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then

          do ieqn = 1, cur_cond%num_other_goveqs
             if (cur_cond%list_id_of_other_goveqs(ieqn) == list_id_of_other_goveq) then

                cur_cond_used = PETSC_TRUE

                do iconn = 1, cur_conn_set%num_connections
                   sum_conn = sum_conn + 1

                   cell_id = cur_conn_set%conn(iconn)%GetIDDn()

                   internal_conn = PETSC_FALSE
                   cond_type     = cur_cond%itype

                   if (.not.cur_cond%swap_order) then
                      select case(this%aux_vars_conn_bc(sum_conn)%flux_type)

                      case (DARCY_FLUX_TYPE)

                         call RichardsFlux(                    &
                              this%aux_vars_bc(sum_conn),      &
                              this%aux_vars_in(cell_id),       &
                              cur_conn_set%conn(iconn),        &
                              compute_deriv,                   &
                              internal_conn,                   &
                              cond_type,                       &
                              dummy_var,                       &
                              Jup,                             &
                              Jdn                              &
                              )

                      case (CONDUCTANCE_FLUX_TYPE)
                         call RichardsFluxConductanceModel(    &
                              this%aux_vars_bc(sum_conn ),     &
                              this%aux_vars_in(cell_id ),      &
                              this%aux_vars_conn_bc(sum_conn), &
                              cur_conn_set%conn(iconn),        &
                              compute_deriv,                   &
                              internal_conn,                   &
                              cond_type,                       &
                              dummy_var,                       &
                              Jup,                             &
                              Jdn                              &
                              )

                      case default
                         write(iulog,*) 'Unknown flux_type ', this%aux_vars_conn_in(sum_conn)%flux_type
                         call endrun(msg=errMsg(__FILE__,__LINE__))
                      end select

                      val = -Jup

                   else

                      select case(this%aux_vars_conn_bc(sum_conn)%flux_type)

                      case (DARCY_FLUX_TYPE)

                         call RichardsFlux(                    &
                              this%aux_vars_in(cell_id ),      &
                              this%aux_vars_bc(sum_conn),      &
                              cur_conn_set%conn(iconn),        &
                              compute_deriv,                   &
                              internal_conn,                   &
                              cond_type,                       &
                              dummy_var,                       &
                              Jup,                             &
                              Jdn                              &
                              )

                      case (CONDUCTANCE_FLUX_TYPE)
                         call RichardsFluxConductanceModel(    &
                              this%aux_vars_in(cell_id ),      &
                              this%aux_vars_bc(sum_conn ),     &
                              this%aux_vars_conn_bc(sum_conn), &
                              cur_conn_set%conn(iconn),        &
                              compute_deriv,                   &
                              internal_conn,                   &
                              cond_type,                       &
                              dummy_var,                       &
                              Jup,                             &
                              Jdn                              &
                              )

                      case default
                         write(iulog,*) 'Unknown flux_type ', this%aux_vars_conn_in(sum_conn)%flux_type
                         call endrun(msg=errMsg(__FILE__,__LINE__))
                      end select
                      val = Jdn
                   endif

                   row = cell_id - 1
                   col = cur_conn_set%conn(iconn)%GetIDUp() - 1

                   call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

                enddo

             endif
          enddo

       endif

       if (.not. cur_cond_used) sum_conn = sum_conn + cur_cond%conn_set%num_connections

       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureJacOffDiag_BC

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureJacOffDiag_Temp(this, list_id_of_other_goveq, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the derivative of residual equation with respect to
    ! temperature
    !
    ! !USES:
    use RichardsMod               , only : RichardsFluxDerivativeWrtTemperature
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_NULL
    use CouplingVariableType      , only : coupling_variable_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    PetscInt                                 :: list_id_of_other_goveq
    Mat                                      :: B
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    type(coupling_variable_type), pointer    :: cpl_var
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id_dn
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id
    PetscInt                                 :: row
    PetscInt                                 :: col
    PetscReal                                :: dummy_var
    PetscReal                                :: Jup
    PetscReal                                :: Jdn
    PetscReal                                :: val
    PetscBool                                :: compute_deriv
    PetscBool                                :: internal_conn
    PetscInt                                 :: cond_type
    PetscInt                                 :: ieqn
    PetscReal                                :: por
    PetscReal                                :: den
    PetscReal                                :: sat
    PetscReal                                :: dpor_dT
    PetscReal                                :: dden_dT
    PetscReal                                :: dsat_dT
    PetscReal                                :: derivative
    PetscReal                                :: dtInv
    PetscBool                                :: coupling_via_BC
    PetscBool                                :: eqns_are_coupled
    PetscBool                                :: cur_cond_used
    PetscInt                                 :: ivar
    type(condition_type),pointer             :: cur_cond
    type(connection_set_type), pointer       :: cur_conn_set

    coupling_via_BC  = PETSC_FALSE
    eqns_are_coupled = PETSC_FALSE
    compute_deriv    = PETSC_TRUE

    ! Are the two equations coupled?
    cpl_var => this%coupling_vars%first
    do
       if (.not.associated(cpl_var)) exit
       
       if (cpl_var%rank_of_coupling_goveqn == list_id_of_other_goveq) then
          eqns_are_coupled = PETSC_TRUE
          exit
       endif

       cpl_var => cpl_var%next
    enddo

    if (.not.eqns_are_coupled) return

    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit

       cur_cond_used = PETSC_FALSE

       if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then

          do ieqn = 1, cur_cond%num_other_goveqs
             if (cur_cond%list_id_of_other_goveqs(ieqn) == list_id_of_other_goveq) then

                coupling_via_BC = PETSC_TRUE

                cur_conn_set => cur_cond%conn_set

                do iconn = 1, cur_conn_set%num_connections
                   sum_conn = sum_conn + 1

                   cell_id = cur_conn_set%conn(iconn)%GetIDDn()

                   internal_conn = PETSC_FALSE
                   cond_type     = cur_cond%itype

                   if (.not.cur_cond%swap_order) then

                      call RichardsFluxDerivativeWrtTemperature( &
                           this%aux_vars_bc(sum_conn),           &
                           this%aux_vars_in(cell_id),            &
                           cur_conn_set%conn(iconn),             &
                           compute_deriv,                        &
                           internal_conn,                        &
                           cond_type,                            &
                           dummy_var,                            &
                           Jup,                                  &
                           Jdn                                   &
                           )

                      val = Jup

                   else

                      call RichardsFluxDerivativeWrtTemperature( &
                           this%aux_vars_in(cell_id),            &
                           this%aux_vars_bc(sum_conn),           &
                           cur_conn_set%conn(iconn),             &
                           compute_deriv,                        &
                           internal_conn,                        &
                           cond_type,                            &
                           dummy_var,                            &
                           Jup,                                  &
                           Jdn                                   &
                           )

                      val = -Jdn

                   endif

                   row = cell_id - 1
                   col = cur_conn_set%conn(iconn)%GetIDUp() - 1

                   call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

                enddo

             endif
          enddo

       endif

       if (.not. cur_cond_used) sum_conn = sum_conn + cur_cond%conn_set%num_connections

       cur_cond => cur_cond%next
    enddo

    if (coupling_via_BC) return


    dtInv = 1.d0 / this%dtime

    ! Interior cells
    do cell_id = 1, this%mesh%ncells_local

       if (this%mesh%is_active(cell_id) ) then
          por     = this%aux_vars_in(cell_id)%por
          den     = this%aux_vars_in(cell_id)%den
          sat     = this%aux_vars_in(cell_id)%sat

          dpor_dT = 0.d0
          dden_dT = this%aux_vars_in(cell_id)%dden_dT
          dsat_dT = 0.d0

          derivative = (dpor_dT*den    *sat     + &
                        por    *dden_dT*sat     + &
                        por    *den    *dsat_dT   )*this%mesh%vol(cell_id) * dtInv

       else
          derivative = 1.d0
       endif

       row = cell_id - 1
       col = cell_id - 1

       call MatSetValuesLocal(B, 1, row, 1, col, derivative, ADD_VALUES, ierr); CHKERRQ(ierr)

    enddo

    ! Interior cells
    cur_conn_set => this%mesh%intrn_conn_set_list%first
    sum_conn = 0
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id_up = cur_conn_set%conn(sum_conn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(sum_conn)%GetIDDn()

          internal_conn = PETSC_TRUE
          cond_type     = COND_NULL

          if ( (.not. this%mesh%is_active(cell_id_up)) .or. &
               (.not. this%mesh%is_active(cell_id_dn)) ) cycle

          call RichardsFluxDerivativeWrtTemperature(      &
               this%aux_vars_in(cell_id_up),     &
               this%aux_vars_in(cell_id_dn),     &
               cur_conn_set%conn(iconn),                  &
               compute_deriv,                             &
               internal_conn,                             &
               cond_type,                                 &
               dummy_var,                                 &
               Jup,                                       &
               Jdn                                        &
               )

          row = cell_id_up - 1
          col = cell_id_up - 1
          val = -Jup
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_up - 1
          col = cell_id_dn - 1
          val = -Jdn
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_dn - 1
          col = cell_id_up - 1
          val = Jup
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_dn - 1
          col = cell_id_dn - 1
          val = Jdn
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

  end subroutine RichardsODEPressureJacOffDiag_Temp

  !------------------------------------------------------------------------
  subroutine RichardsODEComputeLateralFlux(this)
    !
    ! !DESCRIPTION:
    ! Computes the divergence associated with internal and boundary
    ! conditions for residual equation. Also computes the contribution of
    ! source-sink to residual equation.
    !
    ! \int \nabla \cdot (\rho \mathbf{q}) dV
    ! = \int (\rho \mathbf{q}) \cdot \mathbf{n} dA
    !
    ! !USES:
    use RichardsMod                 , only : RichardsFlux
    use ConditionType               , only : condition_type
    use ConnectionSetType           , only : connection_set_type
    use MultiPhysicsProbConstants   , only : COND_NULL
    use MultiPhysicsProbConstants   , only : FMWH2O
    use MultiPhysicsProbConstants   , only : COND_MASS_RATE
    use MultiPhysicsProbConstants   , only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants   , only : CONDUCTANCE_FLUX_TYPE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id_dn
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id
    PetscReal                                :: flux
    PetscReal                                :: dummy_var1
    PetscReal                                :: dummy_var2
    PetscBool                                :: compute_deriv
    PetscBool                                :: internal_conn
    PetscInt                                 :: cond_type
    type(condition_type),pointer             :: lateral_cond
    type(connection_set_type), pointer       :: cur_conn_set
    PetscBool                                :: lateral_cond_found
    PetscErrorCode                           :: ierr
    PetscViewer :: viewer
    character(len=256) :: string
    character(len=256) :: rank_string
    Vec :: xx
    PetscReal      , pointer       :: real_ptr(:)                   ! temporary

    compute_deriv = PETSC_FALSE

    call VecCreate(PETSC_COMM_SELF, xx, ierr); CHKERRQ(ierr)
    call VecSetSizes(xx, this%mesh%ncells_local, PETSC_DECIDE, ierr); CHKERRQ(ierr)
    call VecSetBlockSize(xx, 1, ierr); CHKERRQ(ierr)
    call VecSetFromOptions(xx, ierr); CHKERRQ(ierr)

    ! Source-sink cells
    lateral_cond_found = PETSC_FALSE
    lateral_cond => this%source_sinks%first
    do
       if (.not.associated(lateral_cond)) exit

       if (trim(lateral_cond%name) == 'Lateral_flux') then
          lateral_cond_found = PETSC_TRUE
          exit
       endif
       lateral_cond => lateral_cond%next
    enddo

    if (.not.lateral_cond_found) then
       write(iulog,*)'RichardsODEComputeLateralFlux: Lateral source-sink not found.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    lateral_cond%value(:) = 0.d0

    ! Lateral cells
    cur_conn_set => this%mesh%lateral_conn_set_list%first
    sum_conn = 0
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id_up = cur_conn_set%conn(sum_conn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(sum_conn)%GetIDDn()

          internal_conn = PETSC_TRUE
          cond_type     = COND_NULL

          if ( (.not. this%mesh%is_active(cell_id_up)) .or. &
               (.not. this%mesh%is_active(cell_id_dn)) ) cycle

          select case(this%aux_vars_conn_in(sum_conn)%flux_type)

          case (DARCY_FLUX_TYPE)

             call RichardsFlux( &
                  this%aux_vars_in(cell_id_up),     &
                  this%aux_vars_in(cell_id_dn),     &
                  cur_conn_set%conn(iconn),                  &
                  compute_deriv,                             &
                  internal_conn,                             &
                  cond_type,                                 &
                  flux,                                      &
                  dummy_var1,                                &
                  dummy_var2                                 &
                  )

          case (CONDUCTANCE_FLUX_TYPE)
             write(iulog,*) 'Add code for CONDUCTANCE_FLUX_TYPE'
             call endrun(msg=errMsg(__FILE__,__LINE__))

          case default
             write(iulog,*) 'Unknown flux_type ', this%aux_vars_conn_in(sum_conn)%flux_type
             call endrun(msg=errMsg(__FILE__,__LINE__))
          end select

          if (cell_id_up <= this%mesh%ncells_local) &
               lateral_cond%value(cell_id_up) = lateral_cond%value(cell_id_up) + flux*FMWH2O
          if (cell_id_dn <= this%mesh%ncells_local) &
               lateral_cond%value(cell_id_dn) = lateral_cond%value(cell_id_dn) - flux*FMWH2O
       enddo

       cur_conn_set => cur_conn_set%next
    enddo

  end subroutine RichardsODEComputeLateralFlux

  !------------------------------------------------------------------------

  subroutine RichardsODEPressurePreSolve(this)
    !
    ! !DESCRIPTION:
    ! Performs pre-solve computation
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this

    ! !LOCAL VARIABLES
    PetscErrorCode                           :: ierr
    PetscReal, dimension(:), pointer         :: f_p

    ! Computes contribution to residual equation (accumulation term) based on
    ! the pressure values from previous time step.

    call VecGetArrayF90(this%accum_prev, f_p, ierr); CHKERRQ(ierr)
    call RichardsODEPressureAccum(this, f_p)
    call VecRestoreArrayF90(this%accum_prev, f_p, ierr); CHKERRQ(ierr)

  end subroutine RichardsODEPressurePreSolve

  !------------------------------------------------------------------------

  subroutine RichardsODEPressurePreStepDT(this)
    !
    ! !DESCRIPTION:
    ! Performs pre-step computation
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this

    this%lat_mass_exc(:) = 0.d0
    this%bnd_mass_exc(:) = 0.d0

  end subroutine RichardsODEPressurePreStepDT

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureSetSOEAuxVarOffsets(this, &
       bc_offset_count, bc_offsets, &
       ss_offset_count, ss_offsets)
    !
    ! !DESCRIPTION:
    ! For auxvars corresponding to boundary and source-sink condition, set
    ! offsets corresponding to auxvars of SoE.
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_NULL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    PetscInt                                 :: bc_offset_count
    PetscInt, pointer                        :: bc_offsets(:)
    PetscInt                                 :: ss_offset_count
    PetscInt, pointer                        :: ss_offsets(:)
    !
    ! !LOCAL VARIABLES
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: cond_count
    PetscInt                                 :: cond_itype_to_exclude

    cond_itype_to_exclude = COND_NULL
    call this%GetNConditionsExcptCondItype(COND_BC, &
         cond_itype_to_exclude, cond_count)

    if (bc_offset_count > cond_count) then
       write(iulog,*) 'ERROR: bc_offset_count > cond_count'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cond_itype_to_exclude = COND_NULL
    call this%GetNConditionsExcptCondItype(COND_SS, &
         cond_itype_to_exclude, cond_count)

    if (ss_offset_count > cond_count) then
       write(iulog,*) 'ERROR: ss_offset_count > cond_count'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    allocate(this%soe_auxvars_bc_offset(bc_offset_count))
    cond_count = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       cond_count = cond_count + 1
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then
          this%soe_auxvars_bc_offset(cond_count) = bc_offsets(cond_count)
       endif
       cur_cond => cur_cond%next
    enddo

    allocate(this%soe_auxvars_ss_offset(ss_offset_count))
    cond_count = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       cond_count = cond_count + 1
       this%soe_auxvars_ss_offset(cond_count) = ss_offsets(cond_count)
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureSetSOEAuxVarOffsets

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureCreateVectors(this)
    !
    ! !DESCRIPTION:
    ! Creates required PETSc vectors.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    !
    PetscErrorCode                  :: ierr

    call VecCreateSeq(PETSC_COMM_SELF, this%mesh%ncells_local, &
         this%accum_prev, ierr)
    CHKERRQ(ierr)

  end subroutine RichardsODEPressureCreateVectors

  !------------------------------------------------------------------------
  subroutine RichardsODESetAuxVarConnRealValue(this, auxvar_type, var_type, &
       var_value)
    !
    ! !DESCRIPTION:
    ! Set real value in connection aux vars
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type)         :: this
    PetscInt                            , intent(in) :: auxvar_type
    PetscInt                            , intent(in) :: var_type
    PetscReal                 , pointer , intent(in) :: var_value(:)
    !
    ! !LOCAL VARAIABLES:
    type(connection_set_type) , pointer              :: cur_conn_set
    type(condition_type)      , pointer              :: cur_cond
    PetscInt                                         :: iconn
    PetscInt                                         :: nvars
    PetscInt                                         :: sum_conn

    nvars = size(var_value)

    select case(auxvar_type)
    case (AUXVAR_CONN_INTERNAL)
       ! Interior cells
       cur_conn_set => this%mesh%intrn_conn_set_list%first
       sum_conn = 0
       do
          if (.not.associated(cur_conn_set)) exit

          if (sum_conn + cur_conn_set%num_connections > nvars) then
             write(iulog,*) 'Number of values being set in aux_vars_conn_in is less than' // &
                  ' total number of aux_var_conn_in'
             write(iulog,*)'nvars    = ',nvars
             write(iulog,*)'sum_conn = ',sum_conn
             write(iulog,*)'num_conn = ',cur_conn_set%num_connections
             call endrun(msg=errMsg(__FILE__,__LINE__))
          endif

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             call this%aux_vars_conn_in(sum_conn)%SetRealValue(var_type, var_value(sum_conn))
          enddo
          cur_conn_set => cur_conn_set%next
       enddo

    case (AUXVAR_CONN_BC)
       ! Boundary cells
       sum_conn = 0
       cur_cond => this%boundary_conditions%first
       do
          if (.not.associated(cur_cond)) exit

          cur_conn_set => cur_cond%conn_set

          if (sum_conn + cur_conn_set%num_connections > nvars) then
             write(iulog,*) 'Number of values being set in aux_vars_conn_bc is less than' // &
                  ' total number of aux_var_conn_bc'
             write(iulog,*)'nvars    = ',nvars
             write(iulog,*)'sum_conn = ',sum_conn
             write(iulog,*)'num_conn = ',cur_conn_set%num_connections
             call endrun(msg=errMsg(__FILE__,__LINE__))
          endif

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             call this%aux_vars_conn_bc(sum_conn)%SetRealValue(var_type, var_value(sum_conn))
          end do

          cur_cond => cur_cond%next
       enddo

    case default
       write(iulog,*)'RichardsODESetAuxVarConnRealValue: Unkown auxvar_type'
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

  end subroutine RichardsODESetAuxVarConnRealValue

  !------------------------------------------------------------------------
  subroutine RichardsODESetAuxVarConnIntValue(this, auxvar_type, var_type, &
       var_value)
    !
    ! !DESCRIPTION:
    ! Set integer value in connection aux vars
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type)         :: this
    PetscInt                            , intent(in) :: auxvar_type
    PetscInt                            , intent(in) :: var_type
    PetscInt                  , pointer , intent(in) :: var_value(:)
    !
    !
    ! !LOCAL VARAIABLES:
    type(connection_set_type) , pointer              :: cur_conn_set
    type(condition_type)      , pointer              :: cur_cond
    PetscInt                                         :: iconn
    PetscInt                                         :: nvars
    PetscInt                                         :: sum_conn

    nvars = size(var_value)

    select case(auxvar_type)
    case (AUXVAR_CONN_INTERNAL)
       ! Interior cells
       cur_conn_set => this%mesh%intrn_conn_set_list%first
       sum_conn = 0
       do
          if (.not.associated(cur_conn_set)) exit

          if (sum_conn + cur_conn_set%num_connections > nvars) then
             write(iulog,*) 'Number of values being set in aux_vars_conn_in is less than' // &
                  ' total number of aux_var_conn_in'
             write(iulog,*)'nvars    = ',nvars
             write(iulog,*)'sum_conn = ',sum_conn
             write(iulog,*)'num_conn = ',cur_conn_set%num_connections
             call endrun(msg=errMsg(__FILE__,__LINE__))
          endif

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             call this%aux_vars_conn_in(sum_conn)%SetIntValue(var_type, var_value(sum_conn))
          enddo
          cur_conn_set => cur_conn_set%next
       enddo

    case (AUXVAR_CONN_BC)
       ! Boundary cells
       sum_conn = 0
       cur_cond => this%boundary_conditions%first
       do
          if (.not.associated(cur_cond)) exit

          cur_conn_set => cur_cond%conn_set

          if (sum_conn + cur_conn_set%num_connections > nvars) then
             write(iulog,*) 'Number of values being set in aux_vars_conn_bc is less than' // &
                  ' total number of aux_var_conn_bc'
             write(iulog,*)'nvars    = ',nvars
             write(iulog,*)'sum_conn = ',sum_conn
             write(iulog,*)'num_conn = ',cur_conn_set%num_connections
             call endrun(msg=errMsg(__FILE__,__LINE__))
          endif

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             call this%aux_vars_conn_bc(sum_conn)%SetIntValue(var_type, var_value(sum_conn))
          end do

          cur_cond => cur_cond%next
       enddo


    case default
       write(iulog,*)'RichardsODESetAuxVarConnIntValue: Unkown auxvar_type'
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

  end subroutine RichardsODESetAuxVarConnIntValue

  !------------------------------------------------------------------------
  subroutine RichardsODEGetNumInternalConnections(this, nconn_in)
    !
    ! !DESCRIPTION:
    ! Returns the number of internal connection
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_richards_ode_pressure_type) :: this
    PetscInt, intent(out)                    :: nconn_in
    !
    !
    ! !LOCAL VARAIABLES:
    type(connection_set_type) , pointer      :: cur_conn_set

    ! Interior cells
    cur_conn_set => this%mesh%intrn_conn_set_list%first
    nconn_in = 0
    do
       if (.not.associated(cur_conn_set)) exit

       nconn_in = nconn_in + cur_conn_set%num_connections

       cur_conn_set => cur_conn_set%next
    end do

  end subroutine RichardsODEGetNumInternalConnections

#endif

end module GoveqnRichardsODEPressureType
