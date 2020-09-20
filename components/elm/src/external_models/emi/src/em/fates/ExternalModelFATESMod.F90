module ExternalModelFATESMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides
  !
  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use decompMod                    , only : bounds_type
  use EMI_DataMod, only : emi_data_list, emi_data
  use ExternalModelBaseType        , only : em_base_type
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  !
  implicit none
  !

  type, public, extends(em_base_type) :: em_fates_type
     integer :: index_l2e_col_gridcell_index
     integer :: index_l2e_col_patch_index

     integer :: index_l2e_flux_solad
     integer :: index_l2e_flux_solai

     integer :: index_e2l_state_fsun
     integer :: index_e2l_state_laisun
     integer :: index_e2l_state_laisha
   contains
     procedure, public :: Populate_L2E_List       => EM_FATES_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_FATES_Populate_E2L_List
     procedure, public :: Solve                   => EM_FATES_Solve
  end type em_fates_type

contains

  !------------------------------------------------------------------------
  subroutine EM_FATES_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by FATES from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_fates_type)                :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !

    call EM_FATES_Populate_L2E_List_For_Sunfrac_Stage(this, l2e_list)

  end subroutine EM_FATES_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_FATES_Populate_E2L_List(this, e2l_list)
    !
    !
    ! !DESCRIPTION:
    ! Create a list of all variables to be returned by FATES from ALM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_fates_type)                 :: this
    class(emi_data_list) , intent(inout) :: e2l_list

    call EM_FATES_Populate_E2L_List_For_Surfac_Stage(this, e2l_list)

  end subroutine EM_FATES_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_FATES_Populate_L2E_List_For_Sunfrac_Stage(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by FATES from ALM for
    ! SUNFRAC computation
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_fates_type)                :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_FATES_SUNFRAC_STAGE

    id                                = L2E_FLUX_SOLAR_DIRECT_RADDIATION
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_solad         = index

    id                                = L2E_FLUX_SOLAR_DIFFUSE_RADDIATION
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_solai         = index

    id                                = L2E_COLUMN_GRIDCELL_INDEX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_gridcell_index = index

    id                                = L2E_COLUMN_PATCH_INDEX_BEGIN
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_patch_index    = index

    deallocate(em_stages)

  end subroutine EM_FATES_Populate_L2E_List_For_Sunfrac_Stage

  !------------------------------------------------------------------------
  subroutine EM_FATES_Populate_E2L_List_For_Surfac_Stage(this, e2l_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables to be returned by FATES from ALM
    ! after SUNFRAC computation
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_fates_type)                 :: this
    class(emi_data_list) , intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_FATES_SUNFRAC_STAGE

    id                          = E2L_STATE_FSUN
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_fsun   = index

    id                          = E2L_STATE_LAISUN
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_laisun = index

    id                          = E2L_STATE_LAISHA
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_laisha = index

    deallocate(em_stages)

  end subroutine EM_FATES_Populate_E2L_List_For_Surfac_Stage

    !------------------------------------------------------------------------
  subroutine EM_FATES_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use clm_varctl                , only : iulog
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_fates_type)                 :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent (in)   :: bounds_clump

    select case(em_stage)
    case (EM_FATES_SUNFRAC_STAGE)
       call EM_FATES_Sunfrac_Solve(this, clump_rank, l2e_list, e2l_list)
    case default
       write(iulog,*)'EM_FATES_Solve: Unknown em_stage.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine EM_FATES_Solve

    !------------------------------------------------------------------------
  subroutine EM_FATES_Sunfrac_Solve(this, clump_rank, l2e_list, e2l_list)
    !
    ! !DESCRIPTION:
    ! This interface function is a wrapper call on ED_SunShadeFracs. The only
    ! returned variable is a patch vector, fsun_patch, which describes the fraction
    ! of the canopy that is exposed to sun.
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use clm_varctl                , only : iulog
#ifdef FATES_VIA_EMI
    use elm_instMod               , only : clm_fates
    use EDSurfaceRadiationMod     , only : ED_SunShadeFracs
#endif
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_fates_type)                 :: this
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    !
    ! Local Variables
    real(r8)  , pointer                  :: l2e_solad(:,:)
    real(r8)  , pointer                  :: l2e_solai(:,:)
    real(r8)  , pointer                  :: e2l_fsun(:)
    real(r8)  , pointer                  :: e2l_laisun(:)
    real(r8)  , pointer                  :: e2l_laisha(:)
    integer   , pointer                  :: l2e_col_gridcell(:)
    integer   , pointer                  :: l2e_col_patchi(:)
    integer  :: p                           ! global index of the host patch
    integer  :: g                           ! global index of the host gridcell
    integer  :: c                           ! global index of the host column
    integer  :: s                           ! FATES site index
    integer  :: ifp                         ! FATEs patch index
                                            ! this is the order increment of patch
                                            ! on the site

    call l2e_list%GetPointerToInt1D(this%index_l2e_col_gridcell_index, l2e_col_gridcell)
    call l2e_list%GetPointerToInt1D(this%index_l2e_col_patch_index   , l2e_col_patchi)
    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_solad       , l2e_solad)
    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_solai       , l2e_solai)
    call e2l_list%GetPointerToReal1D(this%index_e2l_state_fsun       , e2l_fsun)
    call e2l_list%GetPointerToReal1D(this%index_e2l_state_laisun     , e2l_laisun)
    call e2l_list%GetPointerToReal1D(this%index_e2l_state_laisha     , e2l_laisha)

#ifdef FATES_VIA_EMI
        ! -------------------------------------------------------------------------------
        ! Convert input BC's
        ! The sun-shade calculations are performed only on FATES patches
        ! -------------------------------------------------------------------------------

    do s = 1, clm_fates%fates(clump_rank)%nsites
       c = clm_fates%f2hmap(clump_rank)%fcolumn(s)
       g = l2e_col_gridcell(c)

       do ifp = 1, clm_fates%fates(clump_rank)%sites(s)%youngest_patch%patchno

          p = ifp + l2e_col_patchi(c)
          clm_fates%fates(clump_rank)%bc_in(s)%solad_parb(ifp,:) = l2e_solad(g,:)
          clm_fates%fates(clump_rank)%bc_in(s)%solai_parb(ifp,:) = l2e_solai(g,:)

       end do
    end do

    ! -------------------------------------------------------------------------------
    ! Call FATES public function to calculate internal sun/shade structures
    ! as well as total patch sun/shade fraction output boundary condition
    ! -------------------------------------------------------------------------------

    call ED_SunShadeFracs(clm_fates%fates(clump_rank)%nsites, &
          clm_fates%fates(clump_rank)%sites,  &
          clm_fates%fates(clump_rank)%bc_in,  &
          clm_fates%fates(clump_rank)%bc_out)

    ! -------------------------------------------------------------------------------
    ! Transfer the FATES output boundary condition for canopy sun/shade fraction
    ! to the HLM
    ! -------------------------------------------------------------------------------

    do s = 1, clm_fates%fates(clump_rank)%nsites
       c = clm_fates%f2hmap(clump_rank)%fcolumn(s)
       do ifp = 1, clm_fates%fates(clump_rank)%sites(s)%youngest_patch%patchno
          p = ifp + l2e_col_patchi(c)
          e2l_fsun(p)   = clm_fates%fates(clump_rank)%bc_out(s)%fsun_pa(ifp)
          e2l_laisun(p) = clm_fates%fates(clump_rank)%bc_out(s)%laisun_pa(ifp)
          e2l_laisha(p) = clm_fates%fates(clump_rank)%bc_out(s)%laisha_pa(ifp)
       end do
    end do
#else
       call endrun('FATES is on but code was not compiled with -DFATES_VIA_EMI')
#endif

  end subroutine EM_FATES_Sunfrac_Solve

end module ExternalModelFATESMod
