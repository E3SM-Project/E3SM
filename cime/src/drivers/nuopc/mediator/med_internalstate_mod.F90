module med_internalstate_mod

  !-----------------------------------------------------------------------------
  ! Mediator Internal State Datatype.
  !-----------------------------------------------------------------------------

  use ESMF                  , only : ESMF_RouteHandle, ESMF_FieldBundle, ESMF_State
  use ESMF                  , only : ESMF_VM
  use esmFlds               , only : ncomps
  use shr_nuopc_fldList_mod , only : nmappers

  implicit none
  private

  integer, public :: logunit  ! logunit for mediator log output
  integer, public :: loglevel ! loglevel for mediator log output
  logical, public :: mastertask=.false. ! is this the mastertask

  ! Active coupling definitions
  ! This defines the med_mapping_allowed is a starting point for what is
  ! allowed in this coupled system.  It will be revised further after the system
  ! starts, but any coupling set to false will never be allowed.  As new connections
  ! are allowed, just update the table below.
  ! - the rows are the destination of coupling
  ! - the columns are the source of coupling
  ! - So, the second column indicates which models the atm is coupled to.
  ! - And the second row indicates which models are coupled to the atm.
  ! The mediator is not connected to any components because the mediator
  ! doesn't have it's own grid and only acts as a hub.

  ! tcraig, turned off glc2ocn and glc2ice for time being
  logical, public, parameter :: med_coupling_allowed(ncomps,ncomps) = &
   (/ .false., .false., .false., .false., .false., .false., .false., .false., &  ! med
      .false., .false., .true. , .true. , .true. , .false., .false., .false., &  ! atm
      .false., .true. , .false., .false., .false., .true. , .false., .true. , &  ! lnd
      .false., .true. , .false., .false., .true. , .true. , .true. , .false., &  ! ocn
      .false., .true. , .false., .true. , .false., .true. , .false., .false., &  ! ice
      .false., .false., .true. , .false., .false., .false., .false., .false., &  ! rof
      .false., .true. , .false., .true. , .true. , .false., .false., .false., &  ! wav
      .false., .false., .true. , .false., .false., .false., .false., .false.  /) ! glc
   !   med      atm      lnd      ocn      ice      rof      wav      glc

  ! private internal state to keep instance data
  type InternalStateStruct
    ! NState_Imp and NState_Exp are the standard NUOPC coupling datatypes
    ! FBImp and FBExp are the internal mediator datatypes
    ! NState_Exp(n) = FBExp(n), copied in the connector prep phase
    ! FBImp(n,n) = NState_Imp(n), copied in connector post phase
    ! FBImp(n,k) is the FBImp(n,n) interpolated to grid k
    ! RH(n,k,m) is a RH from grid n to grid k, map type m

    logical               :: comp_present(ncomps)               ! comp present flag
    logical               :: med_coupling_active(ncomps,ncomps) ! computes the active coupling
    type(ESMF_RouteHandle):: RH(ncomps,ncomps,nmappers)         ! Routehandles for pairs of components and different mappers
    type(ESMF_FieldBundle):: FBfrac(ncomps)                     ! Fraction data for various components, on their grid
    type(ESMF_FieldBundle):: FBNormOne(ncomps,ncomps,nmappers)  ! Unity static normalization
    type(ESMF_State)      :: NStateImp(ncomps)                  ! Import data from various component, on their grid
    type(ESMF_State)      :: NStateExp(ncomps)                  ! Export data to various component, on their grid
    type(ESMF_FieldBundle):: FBImp(ncomps,ncomps)               ! Import data from various components interpolated to various grids
    type(ESMF_FieldBundle):: FBImpAccum(ncomps)                 ! Accumulator for various components import on their grid
    integer               :: FBImpAccumcnt(ncomps)              ! Accumulator counter for each FBImpAccum
    type(ESMF_FieldBundle):: FBExp(ncomps)                      ! Export data for various components, on their grid
    type(ESMF_FieldBundle):: FBExpAccum(ncomps)                 ! Accumulator for various components export on their grid
    integer               :: FBExpAccumcnt(ncomps)              ! Accumulator counter for each FBExpAccum
    logical               :: FBExpAccumFlag(ncomps) = .false.   ! Accumulator flag, if true accumulation was done
    integer               :: conn_prep_cnt(ncomps)              ! Connector prep count
    integer               :: conn_post_cnt(ncomps)              ! Connector post count
    type(ESMF_VM)         :: vm

    ! CESM-specific internal state fields
    type(ESMF_FieldBundle):: FBMed_ocnalb_o                     ! Ocn albedo on ocn grid
    type(ESMF_FieldBundle):: FBMed_ocnalb_a                     ! Ocn albedo on atm grid
    type(ESMF_FieldBundle):: FBMed_aoflux_o                     ! Ocn/Atm flux fields on ocn grid
    type(ESMF_FieldBundle):: FBMed_aoflux_a                     ! Ocn/Atm flux fields on atm grid
    type(ESMF_FieldBundle):: FBMed_aoflux_diurnl_o              ! Ocn/Atm flux fields only needed for history
    type(ESMF_FieldBundle):: FBMed_aoflux_diurnl_a              ! Ocn/Atm flux fields only needed for history
    type(ESMF_FieldBundle):: FBMed_l2x_to_glc_l                 ! FB only in mediator- Land->glc on lnd grid
    type(ESMF_FieldBundle):: FBMed_l2x_to_glc_accum_l           ! FB only in mediator- Land->glc accumulator on lnd grid
 end type InternalStateStruct

  type, public :: InternalState
    type(InternalStateStruct), pointer :: wrap
 end type InternalState

  !-----------------------------------------------------------------------------

end module med_internalstate_mod
