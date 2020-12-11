module ExternalModelInterfaceMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides an interface to couple ALM with external model
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use spmdMod                               , only : masterproc, iam
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use decompMod                             , only : bounds_type, get_proc_clumps
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod         , only : emi_data_list, emi_data
  use EMI_DataDimensionMod , only : emi_data_dimension_list_type
#ifdef USE_PETSC_LIB
  use ExternalModelVSFMMod                  , only : em_vsfm_type
  use ExternalModelPTMMod                   , only : em_ptm_type
#endif
  use ExternalModelFATESMod                 , only : em_fates_type
  use ExternalModelStubMod                  , only : em_stub_type
  use EMI_TemperatureType_ExchangeMod       , only : EMI_Pack_TemperatureType_at_Column_Level_for_EM
  use EMI_TemperatureType_ExchangeMod       , only : EMI_Unpack_TemperatureType_at_Column_Level_from_EM
  use EMI_WaterStateType_ExchangeMod        , only : EMI_Pack_WaterStateType_at_Column_Level_for_EM
  use EMI_WaterStateType_ExchangeMod        , only : EMI_Unpack_WaterStateType_at_Column_Level_from_EM
  use EMI_SoilStateType_ExchangeMod         , only : EMI_Pack_SoilStateType_at_Column_Level_for_EM
  use EMI_SoilStateType_ExchangeMod         , only : EMI_Unpack_SoilStateType_at_Column_Level_from_EM
  use EMI_SoilHydrologyType_ExchangeMod     , only : EMI_Pack_SoilHydrologyType_at_Column_Level_for_EM
  use EMI_SoilHydrologyType_ExchangeMod     , only : EMI_Unpack_SoilHydrologyType_at_Column_Level_from_EM
  use EMI_WaterFluxType_ExchangeMod         , only : EMI_Pack_WaterFluxType_at_Column_Level_for_EM
  use EMI_WaterFluxType_ExchangeMod         , only : EMI_Unpack_WaterFluxType_at_Column_Level_from_EM
  use EMI_EnergyFluxType_ExchangeMod        , only : EMI_Pack_EnergyFluxType_at_Column_Level_for_EM
  use EMI_CanopyStateType_ExchangeMod       , only : EMI_Unpack_CanopyStateType_at_Patch_Level_from_EM
  use EMI_Atm2LndType_ExchangeMod           , only : EMI_Pack_Atm2LndType_at_Grid_Level_for_EM
  use EMI_ColumnType_Exchange               , only : EMI_Pack_ColumnType_for_EM
  use EMI_Filter_Exchange                   , only : EMI_Pack_Filter_for_EM
  use EMI_Landunit_Exchange                 , only : EMI_Pack_Landunit_for_EM
  use EMI_CNCarbonStateType_ExchangeMod
  !
  implicit none
  !
  private

  integer :: num_em              ! Number of external models
  integer :: nclumps

  ! Index of the various external models (EMs) in a simulation
  integer :: index_em_betr
  integer :: index_em_fates
  integer :: index_em_pflotran
  integer :: index_em_stub
  integer :: index_em_vsfm
  integer :: index_em_ptm

  class(emi_data_list)               , pointer :: l2e_driver_list(:)
  class(emi_data_list)               , pointer :: e2l_driver_list(:)
  class(emi_data_dimension_list_type), pointer :: emid_dim_list
#ifdef USE_PETSC_LIB
  class(em_vsfm_type)                , pointer :: em_vsfm(:)
  class(em_ptm_type)                 , pointer :: em_ptm(:)
#endif
  class(em_fates_type)               , pointer :: em_fates
  class(em_stub_type)                , pointer :: em_stub(:)

  public :: EMI_Determine_Active_EMs
  public :: EMI_Init_EM
  public :: EMI_Driver

contains

  !-----------------------------------------------------------------------
  subroutine EMI_Determine_Active_EMs()
    !
    ! !DESCRIPTION:
    ! Determine which EMs are active
    !
    ! !USES:
    use elm_varctl, only : use_fates
#ifndef FATES_VIA_EMI
    use elm_varctl, only : use_betr
    use elm_varctl, only : use_pflotran
    use elm_varctl, only : use_vsfm
#endif
    use elm_varctl, only : use_petsc_thermal_model
    use elm_varctl, only : use_em_stub
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer :: iem

    ! Initializes
    num_em               = 0
    index_em_betr        = 0
    index_em_fates       = 0
    index_em_pflotran    = 0
    index_em_stub        = 0
    index_em_vsfm        = 0

    nclumps = get_proc_clumps()

    ! Is FATES active?
    if (use_fates) then
       num_em            = num_em + 1
       index_em_fates    = num_em
       allocate(em_fates)
    endif

#ifndef FATES_VIA_EMI
    ! Is BeTR active?
    if (use_betr) then
       num_em            = num_em + 1
       index_em_betr     = num_em
    endif

    ! Is PFLOTRAN active?
    if (use_pflotran) then
       num_em            = num_em + 1
       index_em_pflotran = num_em
    endif

    ! Is VSFM active?
    if (use_vsfm) then
       num_em            = num_em + 1
       index_em_vsfm     = num_em
#ifdef USE_PETSC_LIB
       allocate(em_vsfm(nclumps))
#endif
    endif

    ! Is PETSc based Thermal Model active?
    if (use_petsc_thermal_model) then
       num_em            = num_em + 1
       index_em_ptm      = num_em
#ifdef USE_PETSC_LIB
       allocate(em_ptm(nclumps))
#endif
    endif

#endif

    ! Is Stub EM active?
    if (use_em_stub) then
       num_em            = num_em + 1
       index_em_stub     = num_em
       allocate(em_stub(nclumps))
    endif

    if ( masterproc ) then
       write(iulog,*) 'Number of External Models = ', num_em
       write(iulog,*) '  Is BeTR present?     ',(index_em_betr     >0)
       write(iulog,*) '  Is FATES present?    ',(index_em_fates    >0)
       write(iulog,*) '  Is PFLOTRAN present? ',(index_em_pflotran >0)
       write(iulog,*) '  Is PTM present?      ',(index_em_ptm      >0)
       write(iulog,*) '  Is Stub EM present?  ',(index_em_stub     >0)
       write(iulog,*) '  Is VSFM present?     ',(index_em_vsfm     >0)
    endif

    if (num_em > 1) then
       call endrun(msg='More than 1 external model is not supported.')
    endif

    allocate(l2e_driver_list(num_em*nclumps))
    allocate(e2l_driver_list(num_em*nclumps))

    do iem = 1, num_em*nclumps
       call l2e_driver_list(iem)%Init()
       call e2l_driver_list(iem)%Init()
    enddo

    allocate(emid_dim_list)
    call emid_dim_list%Init()

  end subroutine EMI_Determine_Active_EMs
  
  !-----------------------------------------------------------------------
  subroutine EMI_Init_EM(em_id)
    !
    ! !DESCRIPTION:
    ! Initialize EMI
    !
    ! !USES:
    use ExternalModelConstants, only : EM_INITIALIZATION_STAGE
    use ExternalModelConstants, only : EM_ID_BETR
    use ExternalModelConstants, only : EM_ID_FATES
    use ExternalModelConstants, only : EM_ID_PFLOTRAN
    use ExternalModelConstants, only : EM_ID_VSFM
    use ExternalModelConstants, only : EM_ID_PTM
    use ExternalModelConstants, only : EM_ID_STUB
#ifndef FATES_VIA_EMI
    use elm_instMod           , only : soilstate_vars
    use elm_instMod           , only : soilhydrology_vars
    use elm_instMod           , only : waterflux_vars
    use elm_instMod           , only : waterstate_vars
#else
    use elm_instMod           , only : soilstate_inst
    use elm_instMod           , only : soilhydrology_inst
    use elm_instMod           , only : waterflux_inst
    use elm_instMod           , only : waterstate_inst
#endif
    use ExternalModelBETRMod  , only : EM_BETR_Populate_L2E_List
    use ExternalModelBETRMod  , only : EM_BETR_Populate_E2L_List
    use decompMod             , only : get_clump_bounds
    use ColumnType            , only : col_pp
    use LandunitType          , only : lun_pp
    use landunit_varcon       , only : istsoil, istcrop,istice
    use column_varcon         , only : icol_road_perv
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer, intent(in)  :: em_id
    !
    ! !LOCAL VARIABLES:
    class(emi_data_list), pointer :: l2e_init_list(:)
    class(emi_data_list), pointer :: e2l_init_list(:)
    integer                       :: em_stage
    integer                       :: ii, c, l
    integer                       :: num_filter_col
    integer                       :: num_filter_lun
    integer, pointer              :: filter_col(:)
    integer, pointer              :: filter_lun(:)
    integer                       :: num_e2l_filter_col
    integer, pointer              :: e2l_filter_col(:)
    integer, pointer              :: tmp_col(:)
    type(bounds_type)             :: bounds_clump
    integer                       :: iem
    integer                       :: clump_rank

    em_stage = EM_INITIALIZATION_STAGE

    select case (em_id)
    case (EM_ID_BETR)

       ! -------------------------------------------------------------
       ! Data need during timestepping
       ! -------------------------------------------------------------

       ! Note: Each thread will exchange exactly the same data between
       !       ALM and FATES
       do clump_rank = 1, nclumps
          iem = (index_em_betr-1)*nclumps + clump_rank
          call EM_BETR_Populate_L2E_List(l2e_driver_list(iem))
          call EM_BETR_Populate_E2L_List(e2l_driver_list(iem))
       enddo

       !$OMP PARALLEL DO PRIVATE (clump_rank, iem, bounds_clump)
       do clump_rank = 1, nclumps

          call get_clump_bounds(clump_rank, bounds_clump)
          iem = (index_em_betr-1)*nclumps + clump_rank

          call EMI_Setup_Data_List(l2e_driver_list(iem), bounds_clump)
          call EMI_Setup_Data_List(e2l_driver_list(iem), bounds_clump)
       enddo
       !$OMP END PARALLEL DO

    case (EM_ID_FATES)

       ! -------------------------------------------------------------
       ! Data need during timestepping
       ! -------------------------------------------------------------

       ! Note: Each thread will exchange exactly the same data between
       !       ALM and FATES
       do clump_rank = 1, nclumps
          iem = (index_em_fates-1)*nclumps + clump_rank
          call em_fates%Populate_L2E_List(l2e_driver_list(iem))
          call em_fates%Populate_E2L_List(e2l_driver_list(iem))
       enddo


       !$OMP PARALLEL DO PRIVATE (clump_rank, iem, bounds_clump)
       do clump_rank = 1, nclumps

          call get_clump_bounds(clump_rank, bounds_clump)
          iem = (index_em_fates-1)*nclumps + clump_rank

          call EMI_Setup_Data_List(l2e_driver_list(iem), bounds_clump)
          call EMI_Setup_Data_List(e2l_driver_list(iem), bounds_clump)
       enddo
       !$OMP END PARALLEL DO

    case (EM_ID_PFLOTRAN)

    case (EM_ID_VSFM)

#ifdef USE_PETSC_LIB
       ! Initialize EM

       ! Initialize lists of data to be exchanged between ALM and VSFM
       ! during initialization step
       allocate(l2e_init_list(nclumps))
       allocate(e2l_init_list(nclumps))

       do clump_rank = 1, nclumps
          iem = (index_em_vsfm-1)*nclumps + clump_rank

          call l2e_init_list(clump_rank)%Init()
          call e2l_init_list(clump_rank)%Init()

          ! Fill the data list:
          !  - Data need during the initialization
          call em_vsfm(clump_rank)%Populate_L2E_Init_List(l2e_init_list(clump_rank))
          call em_vsfm(clump_rank)%Populate_E2L_Init_List(e2l_init_list(clump_rank))

          !  - Data need during timestepping
          call em_vsfm(clump_rank)%Populate_L2E_List(l2e_driver_list(iem))
          call em_vsfm(clump_rank)%Populate_E2L_List(e2l_driver_list(iem))
       enddo

       !$OMP PARALLEL DO PRIVATE (clump_rank, iem, bounds_clump)
       do clump_rank = 1, nclumps

          call get_clump_bounds(clump_rank, bounds_clump)
          iem = (index_em_vsfm-1)*nclumps + clump_rank

          ! Allocate memory for data
          call EMI_Setup_Data_List(l2e_init_list(clump_rank), bounds_clump)
          call EMI_Setup_Data_List(e2l_init_list(clump_rank), bounds_clump)
          call EMI_Setup_Data_List(l2e_driver_list(iem)     , bounds_clump)
          call EMI_Setup_Data_List(e2l_driver_list(iem)     , bounds_clump)

          ! GB_FIX_ME: Create a temporary filter
          num_filter_col = bounds_clump%endc - bounds_clump%begc + 1
          num_filter_lun = bounds_clump%endl - bounds_clump%begl + 1

          allocate(filter_col(num_filter_col))
          allocate(filter_lun(num_filter_lun))

          do ii = 1, num_filter_col
             filter_col(ii) = bounds_clump%begc + ii - 1
          enddo

          do ii = 1, num_filter_lun
             filter_lun(ii) = bounds_clump%begl + ii - 1
          enddo

          ! Reset values in the data list
          call EMID_Reset_Data_for_EM(l2e_init_list(clump_rank), em_stage)
          call EMID_Reset_Data_for_EM(e2l_init_list(clump_rank), em_stage)

          ! Pack all ALM data needed by the external model
          call EMI_Pack_WaterStateType_at_Column_Level_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, waterstate_vars)
          call EMI_Pack_WaterFluxType_at_Column_Level_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, waterflux_vars)
          call EMI_Pack_SoilHydrologyType_at_Column_Level_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, soilhydrology_vars)
          call EMI_Pack_SoilStateType_at_Column_Level_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, soilstate_vars)
          call EMI_Pack_ColumnType_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col)
          call EMI_Pack_Landunit_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_lun, filter_lun)

          ! Ensure all data needed by external model is packed
          call EMID_Verify_All_Data_Is_Set(l2e_init_list(clump_rank), em_stage)

          ! Initialize the external model
          call em_vsfm(clump_rank)%Init(l2e_init_list(clump_rank), e2l_init_list(clump_rank), &
               iam, bounds_clump)

          ! Build a column level filter on which VSFM is active.
          ! This new filter would be used during the initialization to
          ! unpack data from the EM into ALM's data structure.
          allocate(tmp_col(bounds_clump%begc:bounds_clump%endc))

          tmp_col(bounds_clump%begc:bounds_clump%endc) = 0

          num_e2l_filter_col = 0
          do c = bounds_clump%begc,bounds_clump%endc
             if (col_pp%active(c)) then
                l = col_pp%landunit(c)
                if (lun_pp%itype(l) == istsoil .or. col_pp%itype(c) == icol_road_perv .or. &
                    lun_pp%itype(l) == istcrop) then
                   num_e2l_filter_col = num_e2l_filter_col + 1
                   tmp_col(c) = 1
                end if
             end if
          end do

          allocate(e2l_filter_col(num_e2l_filter_col))

          num_e2l_filter_col = 0
          do c = bounds_clump%begc,bounds_clump%endc
             if (tmp_col(c) == 1) then
                num_e2l_filter_col = num_e2l_filter_col + 1
                e2l_filter_col(num_e2l_filter_col) = c
             endif
          enddo

          ! Unpack all data sent from the external model
          call EMI_Unpack_SoilStateType_at_Column_Level_from_EM(e2l_init_list(clump_rank), em_stage, &
               num_e2l_filter_col, e2l_filter_col, soilstate_vars)
          call EMI_Unpack_WaterStateType_at_Column_Level_from_EM(e2l_init_list(clump_rank), em_stage, &
               num_e2l_filter_col, e2l_filter_col, waterstate_vars)
          call EMI_Unpack_WaterFluxType_at_Column_Level_from_EM(e2l_init_list(clump_rank), em_stage, &
               num_e2l_filter_col, e2l_filter_col, waterflux_vars)
          call EMI_Unpack_SoilHydrologyType_at_Column_Level_from_EM(e2l_init_list(clump_rank), em_stage, &
               num_e2l_filter_col, e2l_filter_col, soilhydrology_vars)

          ! Ensure all data sent by external model is unpacked
          call EMID_Verify_All_Data_Is_Set(e2l_init_list(clump_rank), em_stage)

          ! Clean up memory
          call l2e_init_list(clump_rank)%Destroy()
          call e2l_init_list(clump_rank)%Destroy()

          deallocate(e2l_filter_col)
          deallocate(tmp_col)

       enddo
       !$OMP END PARALLEL DO

#else
       call endrun('VSFM is on but code was not compiled with -DUSE_PETSC_LIB')
#endif

    case (EM_ID_PTM)

#ifdef USE_PETSC_LIB

       ! Initialize lists of data to be exchanged between ALM and VSFM
       ! during initialization step
       allocate(l2e_init_list(nclumps))

       do clump_rank = 1, nclumps
          iem = (index_em_ptm - 1)*nclumps + clump_rank

          call l2e_init_list(clump_rank)%Init()

          ! Fill the data list:
          !  - Data need during the initialization
          call em_ptm(clump_rank)%Populate_L2E_Init_List(l2e_init_list(clump_rank))

          !  - Data need during timestepping
          call em_ptm(clump_rank)%Populate_L2E_List(l2e_driver_list(iem))
          call em_ptm(clump_rank)%Populate_E2L_List(e2l_driver_list(iem))
       enddo

       !$OMP PARALLEL DO PRIVATE (clump_rank, iem, bounds_clump)
       do clump_rank = 1, nclumps

          call get_clump_bounds(clump_rank, bounds_clump)
          iem = (index_em_ptm - 1)*nclumps + clump_rank

          ! Allocate memory for data
          call EMI_Setup_Data_List(l2e_init_list(iem), bounds_clump)
          call EMI_Setup_Data_List(l2e_driver_list(iem), bounds_clump)
          call EMI_Setup_Data_List(e2l_driver_list(iem), bounds_clump)

          ! Reset values in the data list
          call EMID_Reset_Data_for_EM(l2e_init_list(clump_rank), em_stage)

          ! GB_FIX_ME: Create a temporary filter
          num_filter_col = bounds_clump%endc - bounds_clump%begc + 1
          num_filter_lun = bounds_clump%endl - bounds_clump%begl + 1

          allocate(filter_col(num_filter_col))
          allocate(filter_lun(num_filter_lun))

          do ii = 1, num_filter_col
             filter_col(ii) = bounds_clump%begc + ii - 1
          enddo

          do ii = 1, num_filter_lun
             filter_lun(ii) = bounds_clump%begl + ii - 1
          enddo

          ! Reset values in the data list
          call EMID_Reset_Data_for_EM(l2e_init_list(clump_rank), em_stage)

          ! Pack all ALM data needed by the external model
          call EMI_Pack_SoilStateType_at_Column_Level_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, soilstate_vars)
          call EMI_Pack_ColumnType_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col)
          call EMI_Pack_Landunit_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_lun, filter_lun)

          ! Ensure all data needed by external model is packed
          call EMID_Verify_All_Data_Is_Set(l2e_init_list(clump_rank), em_stage)

          ! Initialize the external model
          call em_ptm(clump_rank)%Init(l2e_init_list(clump_rank), e2l_init_list(clump_rank), &
               iam, bounds_clump)

          ! Clean up memory
          call l2e_init_list(clump_rank)%Destroy()

       enddo
       !$OMP END PARALLEL DO

#else
       call endrun('PTM is on but code was not compiled with -DUSE_PETSC_LIB')
#endif

    case (EM_ID_STUB)

       !write(iulog,*)'*******************************************'
       write(iulog,*)'  In ELM: Initialization'
       write(iulog,*)'  1.1 Populate lists of variables that will be exchanged between ELM and EM'
       write(iulog,*)'      during initialization and timestepping.'

       ! Initialize lists of data to be exchanged between ALM and VSFM
       ! during initialization step
       allocate(l2e_init_list(nclumps))
       allocate(e2l_init_list(nclumps))

       do clump_rank = 1, nclumps

          iem = (index_em_stub-1)*nclumps + clump_rank

          call l2e_init_list(clump_rank)%Init()
          call e2l_init_list(clump_rank)%Init()

          ! Fill the data list:
          !  - Data need during the initialization
          call em_stub(clump_rank)%Populate_L2E_Init_List(l2e_init_list(clump_rank))
          call em_stub(clump_rank)%Populate_E2L_Init_List(e2l_init_list(clump_rank))

          !  - Data need during timestepping
          call em_stub(clump_rank)%Populate_L2E_List(l2e_driver_list(iem))
          call em_stub(clump_rank)%Populate_E2L_List(e2l_driver_list(iem))

       enddo

       write(iulog,*)'  1.2 Exchange variables between ELM and EM during initialization'

       !$OMP PARALLEL DO PRIVATE (clump_rank, iem, bounds_clump)
       do clump_rank = 1, nclumps

          call get_clump_bounds(clump_rank, bounds_clump)
          iem = (index_em_stub-1)*nclumps + clump_rank

          ! Allocate memory for data
          call EMI_Setup_Data_List(l2e_init_list(clump_rank), bounds_clump)
          call EMI_Setup_Data_List(e2l_init_list(clump_rank), bounds_clump)
          call EMI_Setup_Data_List(l2e_driver_list(iem)     , bounds_clump)
          call EMI_Setup_Data_List(e2l_driver_list(iem)     , bounds_clump)

          ! Reset values in the data list
          call EMID_Reset_Data_for_EM(l2e_init_list(clump_rank), em_stage)
          call EMID_Reset_Data_for_EM(e2l_init_list(clump_rank), em_stage)

          ! GB_FIX_ME: Create a temporary filter
          num_filter_col = bounds_clump%endc - bounds_clump%begc + 1
          num_filter_lun = bounds_clump%endl - bounds_clump%begl + 1

          allocate(filter_col(num_filter_col))
          allocate(filter_lun(num_filter_lun))

          do ii = 1, num_filter_col
             filter_col(ii) = bounds_clump%begc + ii - 1
          enddo

          do ii = 1, num_filter_lun
             filter_lun(ii) = bounds_clump%begl + ii - 1
          enddo

          ! Pack all ALM data needed by the external model
          call EMI_Pack_SoilStateType_at_Column_Level_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, soilstate_vars)

          ! Ensure all data needed by external model is packed
          write(iulog,*)'     1.2.1 Value of variables send by ELM'
          call EMID_Verify_All_Data_Is_Set(l2e_init_list(clump_rank), em_stage, print_data=.true.)

          ! Initialize the external model
          call em_stub(clump_rank)%Init(l2e_init_list(clump_rank), e2l_init_list(clump_rank), &
               iam, bounds_clump)

          ! Unpack all data sent from the external model
          call EMI_Unpack_WaterStateType_at_Column_Level_from_EM(e2l_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, waterstate_vars)

          ! Ensure all data sent by external model is unpacked
          write(iulog,*)'     1.2.4 Value of variables received by ELM'
          call EMID_Verify_All_Data_Is_Set(e2l_init_list(clump_rank), em_stage, print_data=.true.)

          call l2e_init_list(clump_rank)%Destroy()
          call e2l_init_list(clump_rank)%Destroy()

       enddo

    case default
       call endrun('Unknown External Model')
    end select

  end subroutine EMI_Init_EM

  !-----------------------------------------------------------------------
  subroutine EMI_Setup_Data_List(data_list, bounds_clump)
    !
    ! !DESCRIPTION:
    ! Setup the EMI data list
    !
    implicit none
    !
    ! !USES:
    class(emi_data_list) , intent(inout) :: data_list
    type(bounds_type)    , intent (in)   :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    class(emi_data)      , pointer       :: cur_data
    integer                              :: idata

    allocate(data_list%data_ptr(data_list%num_data))

    idata = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit

       idata = idata + 1

       data_list%data_ptr(idata)%data => cur_data
       cur_data => cur_data%next
    enddo

    do idata = 1, data_list%num_data
       cur_data => data_list%data_ptr(idata)%data
       call EMI_Setup_Data(cur_data, bounds_clump)
    enddo

    do idata = 1, data_list%num_data
       call data_list%data_ptr(idata)%data%AllocateMemory()
    enddo

  end subroutine EMI_Setup_Data_List

  !-----------------------------------------------------------------------
  subroutine EMI_Setup_Data(data, bounds_clump)
    !
    ! !DESCRIPTION:
    ! Setup a EMI data
    !
    implicit none
    !
    class(emi_data), pointer, intent(inout) :: data
    type(bounds_type)       , intent (in)   :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer                                 :: ndim
    integer                                 :: nc
    !
    integer                                 :: dim1_beg, dim1_end
    integer                                 :: dim2_beg, dim2_end
    integer                                 :: dim3_beg, dim3_end
    integer                                 :: dim4_beg, dim4_end

    dim1_beg  = 0
    dim2_beg  = 0
    dim3_beg  = 0
    dim4_beg  = 0
    dim1_end  = 0
    dim2_end  = 0
    dim3_end  = 0
    dim4_end  = 0

    ! Determine the dimension values
    select case (data%ndim)
    case (1)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim1_beg_name, dim1_beg)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim1_end_name, dim1_end)

    case (2)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim1_beg_name, dim1_beg)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim2_beg_name, dim2_beg)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim1_end_name, dim1_end)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim2_end_name, dim2_end)

    case (3)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim1_beg_name, dim1_beg)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim2_beg_name, dim2_beg)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim3_beg_name, dim3_beg)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim1_end_name, dim1_end)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim2_end_name, dim2_end)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim3_end_name, dim3_end)

    case (4)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim1_beg_name, dim1_beg)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim2_beg_name, dim2_beg)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim3_beg_name, dim3_beg)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim4_beg_name, dim4_beg)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim1_end_name, dim1_end)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim2_end_name, dim2_end)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim3_end_name, dim3_end)
       call emid_dim_list%GetDimValue( bounds_clump, data%dim4_end_name, dim4_end)

    end select

    call data%SetDimensions(                     &
         dim1_beg, dim1_end, dim2_beg, dim2_end, &
         dim3_beg, dim3_end, dim4_beg, dim4_end  )

  end subroutine EMI_Setup_Data

!-----------------------------------------------------------------------
  subroutine EMI_Driver(em_id, em_stage, dt, number_step,     &
       clump_rank, num_hydrologyc, filter_hydrologyc,         &
       num_nolakec, filter_nolakec,                           &
       num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc, &
       num_filter_lun, filter_lun,                            &
       soilhydrology_vars, soilstate_vars, waterflux_vars,    &
       waterstate_vars, temperature_vars,  atm2lnd_vars,      &
       canopystate_vars, energyflux_vars, carbonstate_vars)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ExternalModelConstants , only : EM_ID_BETR
    use ExternalModelConstants , only : EM_ID_FATES
    use ExternalModelConstants , only : EM_ID_PFLOTRAN
    use ExternalModelConstants , only : EM_ID_VSFM
    use ExternalModelConstants , only : EM_ID_PTM
    use ExternalModelConstants , only : EM_ID_STUB
    use SoilStateType          , only : soilstate_type
    use SoilHydrologyType      , only : soilhydrology_type
    use TemperatureType        , only : temperature_type
    use WaterFluxType          , only : waterflux_type
    use WaterStateType         , only : waterstate_type
    use atm2lndType            , only : atm2lnd_type
    use CanopyStateType        , only : canopystate_type
    use EnergyFluxType         , only : energyflux_type
    use CNCarbonStateType      , only : carbonstate_type
    use ExternalModelBETRMod   , only : EM_BETR_Solve
    use decompMod              , only : get_clump_bounds
    !
    implicit none
    !
    integer                             , intent(in)    :: em_id
    integer                             , intent(in)    :: em_stage
    real(r8)                 , optional , intent(in)    :: dt
    integer                  , optional , intent(in)    :: number_step
    integer                  , optional , intent(in)    :: clump_rank
    integer                  , optional , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , optional , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                  , optional , intent(in)    :: num_nolakec
    integer                  , optional , intent(in)    :: filter_nolakec(:)
    integer                  , optional , intent(in)    :: num_nolakec_and_nourbanc
    integer                  , optional , intent(in)    :: filter_nolakec_and_nourbanc(:)
    integer                  , optional , intent(in)    :: num_filter_lun
    integer                  , optional , intent(in)    :: filter_lun(:)
    type(soilhydrology_type) , optional , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , optional , intent(inout) :: soilstate_vars
    type(waterflux_type)     , optional , intent(inout) :: waterflux_vars
    type(waterstate_type)    , optional , intent(inout) :: waterstate_vars
    type(temperature_type)   , optional , intent(inout) :: temperature_vars
    type(atm2lnd_type)       , optional , intent(inout) :: atm2lnd_vars
    type(canopystate_type)   , optional , intent(inout) :: canopystate_vars
    type(energyflux_type)    , optional , intent(inout) :: energyflux_vars
    type(carbonstate_type)   , optional , intent(inout) :: carbonstate_vars
    !
    integer          :: index_em
    real(r8)         :: dtime
    integer          :: nstep
    type(bounds_type):: bounds_clump
    integer          :: num_filter_col
    integer          :: num_filter_patch
    integer          :: num_filter_grid
    integer, pointer :: filter_col(:)
    integer, pointer :: filter_patch(:)
    integer, pointer :: filter_grid(:)
    integer          :: ii
    integer          :: iem

    ! Find the index_em
    select case (em_id)
    case (EM_ID_BETR)
       index_em = index_em_betr
    case (EM_ID_FATES)
       index_em = index_em_fates
    case (EM_ID_PFLOTRAN)
       index_em = index_em_pflotran
    case (EM_ID_VSFM)
       index_em = index_em_vsfm
    case (EM_ID_PTM)
       index_em = index_em_ptm
    case (EM_ID_STUB)
       index_em = index_em_stub
       write(iulog,*)'     2.1 Value of variables send by ELM'
    case default
       call endrun('Unknown External Model')
    end select

    ! ------------------------------------------------------------------------
    ! Pack the data for EM
    ! ------------------------------------------------------------------------

    if (present(clump_rank)) then
       iem = (index_em-1)*nclumps + clump_rank
    else
       iem = (index_em-1)*nclumps + 1
    endif

    call EMID_Reset_Data_for_EM(l2e_driver_list(iem), em_stage)
    call EMID_Reset_Data_for_EM(e2l_driver_list(iem), em_stage)

    if ( present(temperature_vars) .and. &
         present(num_hydrologyc)   .and. &
         present(filter_hydrologyc)) then

       call EMI_Pack_TemperatureType_at_Column_Level_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, temperature_vars)

       elseif (present(num_nolakec_and_nourbanc)  .and. &
               present(filter_nolakec_and_nourbanc)) then

       call EMI_Pack_TemperatureType_at_Column_Level_for_EM(l2e_driver_list(iem), em_stage, &
            num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc, temperature_vars)
    endif

    if ( present(waterstate_vars)) then
       if (present(num_hydrologyc)  .and. &
           present(filter_hydrologyc)) then

          call EMI_Pack_WaterStateType_at_Column_Level_for_EM(l2e_driver_list(iem), em_stage, &
               num_hydrologyc, filter_hydrologyc, waterstate_vars)

       elseif (present(num_nolakec_and_nourbanc)  .and. &
               present(filter_nolakec_and_nourbanc)) then

          call EMI_Pack_WaterStateType_at_Column_Level_for_EM(l2e_driver_list(iem), em_stage, &
               num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc, waterstate_vars)
       else
          ! GB_FIX_ME: Create a temporary filter
          if (present(clump_rank)) then
             call get_clump_bounds(clump_rank, bounds_clump)
          else
             call get_clump_bounds(1, bounds_clump)
          endif
          num_filter_col = bounds_clump%endc - bounds_clump%begc + 1
          allocate(filter_col(num_filter_col))
          do ii = 1, num_filter_col
             filter_col(ii) = bounds_clump%begc + ii - 1
          enddo

          call EMI_Pack_WaterStateType_at_Column_Level_for_EM(l2e_driver_list(iem), em_stage, &
               num_filter_col, filter_col, waterstate_vars)
          deallocate(filter_col)
       endif
    endif

    if ( present(waterflux_vars) .and. &
         present(num_hydrologyc) .and. &
         present(filter_hydrologyc)) then

       call EMI_Pack_WaterFluxType_at_Column_Level_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, waterflux_vars)
    endif

    if ( present(num_nolakec_and_nourbanc) .and. &
         present(filter_nolakec_and_nourbanc)) then

       call EMI_Pack_EnergyFluxType_at_Column_Level_for_EM(l2e_driver_list(iem), em_stage, &
            num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc, energyflux_vars)

    endif

    if ( present(num_hydrologyc) .and. &
         present(filter_hydrologyc)) then

       call EMI_Pack_Filter_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc)

       call EMI_Pack_ColumnType_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc)

    endif

    if ( present(num_nolakec) .and. &
         present(filter_nolakec)) then

       call EMI_Pack_Filter_for_EM(l2e_driver_list(iem), em_stage, &
            num_nolakec, filter_nolakec)

       call EMI_Pack_ColumnType_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc)

    endif

    if ( present(num_nolakec_and_nourbanc) .and. &
         present(filter_nolakec_and_nourbanc)) then

       call EMI_Pack_Filter_for_EM(l2e_driver_list(iem), em_stage, &
            num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc)

       call EMI_Pack_ColumnType_for_EM(l2e_driver_list(iem), em_stage, &
            num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc)

    endif

    if ( present(num_filter_lun) .and. &
         present(filter_lun)) then

       call EMI_Pack_Landunit_for_EM(l2e_driver_list(iem), em_stage, &
            num_filter_lun, filter_lun)

    endif

    if (present(atm2lnd_vars)) then
       ! GB_FIX_ME: Create a temporary filter
       if (present(clump_rank)) then
          call get_clump_bounds(clump_rank, bounds_clump)
       else
          call get_clump_bounds(1, bounds_clump)
       endif

       num_filter_grid = bounds_clump%endg - bounds_clump%begg + 1

       allocate(filter_col(num_filter_grid))
       do ii = 1, num_filter_grid
          filter_grid(ii) = bounds_clump%begg + ii - 1
       enddo
       call EMI_Pack_Atm2LndType_at_Grid_Level_for_EM(l2e_driver_list(iem), em_stage, &
            num_filter_grid, filter_grid, atm2lnd_vars)
       deallocate(filter_grid)
    endif

    ! GB_FIX_ME: Create a temporary filter
    if (present(clump_rank)) then
       call get_clump_bounds(clump_rank, bounds_clump)
    else
       call get_clump_bounds(1, bounds_clump)
    endif

    num_filter_col = bounds_clump%endc - bounds_clump%begc + 1

    allocate(filter_col(num_filter_col))
    do ii = 1, num_filter_col
       filter_col(ii) = bounds_clump%begc + ii - 1
    enddo
    call EMI_Pack_ColumnType_for_EM(l2e_driver_list(iem), em_stage, &
            num_filter_col, filter_col)
    deallocate(filter_col)

    if (present(carbonstate_vars)  .and. &
         present(num_hydrologyc)   .and. &
         present(filter_hydrologyc)) then
       call EMI_Pack_CNCarbonStateType_at_Column_Level_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, carbonstate_vars)
    endif

    call EMID_Verify_All_Data_Is_Set(l2e_driver_list(iem), em_stage)

    ! ------------------------------------------------------------------------
    ! Solve EM
    ! ------------------------------------------------------------------------
    nstep = 0
    dtime = 0._r8
    if (present(number_step)) nstep = number_step
    if (present(dt))          dtime = dt

    select case (em_id)
    case (EM_ID_BETR)
       call EM_BETR_Solve(em_stage, dtime, nstep, bounds_clump, l2e_driver_list(iem), &
            e2l_driver_list(iem), bounds_clump)

    case (EM_ID_FATES)
       call em_fates%Solve(em_stage, dtime, nstep, clump_rank, l2e_driver_list(iem), &
            e2l_driver_list(iem), bounds_clump)

    case (EM_ID_PFLOTRAN)

    case (EM_ID_VSFM)
#ifdef USE_PETSC_LIB
       call em_vsfm(clump_rank)%Solve(em_stage, dtime, nstep, clump_rank, &
            l2e_driver_list(iem), e2l_driver_list(iem), bounds_clump)
#else
       call endrun('VSFM is on but code was not compiled with -DUSE_PETSC_LIB')
#endif

    case (EM_ID_PTM)
#ifdef USE_PETSC_LIB
       call em_ptm(clump_rank)%Solve(em_stage, dtime, nstep, clump_rank, &
            l2e_driver_list(iem), e2l_driver_list(iem), bounds_clump)
#else
       call endrun('PTM is on but code was not compiled with -DUSE_PETSC_LIB')
#endif

    case (EM_ID_STUB)
       call EMID_Verify_All_Data_Is_Set(l2e_driver_list(iem), em_stage, print_data=.true.)
       call em_stub(clump_rank)%Solve(em_stage, dtime, nstep, clump_rank, &
            l2e_driver_list(iem), e2l_driver_list(iem), bounds_clump)

    case default
       call endrun('Unknown External Model')
    end select

    ! ------------------------------------------------------------------------
    ! Unpack the data for EM
    ! ------------------------------------------------------------------------
    if ( present(waterstate_vars) .and. &
         present(num_hydrologyc)  .and. &
         present(filter_hydrologyc)) then

       call EMI_Unpack_WaterStateType_at_Column_Level_from_EM(e2l_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, waterstate_vars)
    endif

    if ( present(waterflux_vars) .and. &
         present(num_hydrologyc) .and. &
         present(filter_hydrologyc)) then

       call EMI_Unpack_WaterFluxType_at_Column_Level_from_EM(e2l_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, waterflux_vars)
    endif

    if ( present(soilstate_vars) .and. &
         present(num_hydrologyc) .and. &
         present(filter_hydrologyc)) then

       call EMI_Unpack_SoilStateType_at_Column_Level_from_EM(e2l_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, soilstate_vars)
    endif

    if ( present(soilhydrology_vars) .and. &
         present(num_hydrologyc)     .and. &
         present(filter_hydrologyc)) then

       call EMI_Unpack_SoilHydrologyType_at_Column_Level_from_EM(e2l_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, soilhydrology_vars)
    endif

    if (present(canopystate_vars)) then
          ! GB_FIX_ME: Create a temporary filter
          if (present(clump_rank)) then
             call get_clump_bounds(clump_rank, bounds_clump)
          else
             call get_clump_bounds(1, bounds_clump)
          endif
          num_filter_patch = bounds_clump%endp - bounds_clump%begp + 1
          allocate(filter_patch(num_filter_patch))
          do ii = 1, num_filter_patch
             filter_patch(ii) = bounds_clump%begp + ii - 1
          enddo
          call EMI_Unpack_CanopyStateType_at_Patch_Level_from_EM(e2l_driver_list(iem), em_stage, &
               num_filter_patch, filter_patch, canopystate_vars)
          deallocate(filter_patch)
    endif

    if ( present(temperature_vars) .and. &
         present(num_nolakec_and_nourbanc)     .and. &
         present(filter_nolakec_and_nourbanc)) then

       call EMI_Unpack_TemperatureType_at_Column_Level_from_EM(e2l_driver_list(iem), em_stage, &
            num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc, temperature_vars)
    endif

    if (present(carbonstate_vars)  .and. &
         present(num_hydrologyc)   .and. &
         present(filter_hydrologyc)) then
       call EMI_Unpack_CNCarbonStateType_at_Column_Level_from_EM(e2l_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, carbonstate_vars)
    endif

    if (em_id == EM_ID_STUB) then
       write(iulog,*)'     2.4 Value of variables received by ELM'
       call EMID_Verify_All_Data_Is_Set(e2l_driver_list(iem), em_stage, print_data=.true.)
    else
       call EMID_Verify_All_Data_Is_Set(e2l_driver_list(iem), em_stage)
    endif

  end subroutine EMI_Driver
  
!-----------------------------------------------------------------------
  subroutine EMID_Reset_Data_for_EM(data_list, em_stage)
    !
    ! !DESCRIPTION:
    ! Reset all EMI data that will be exchanged between ALM and external
    ! model for em_stage
    !
    implicit none
    !
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    !
    class(emi_data), pointer            :: cur_data
    integer                             :: istage

    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit

       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             call cur_data%Reset()
             exit
          endif
       enddo

       cur_data => cur_data%next
    enddo

  end subroutine EMID_Reset_Data_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Verify_All_Data_Is_Set(data_list, em_stage, print_data)
    !
    ! !DESCRIPTION:
    ! Verify that all EMI data that will be exchanged between ALM and external
    ! model for em_stage was set
    !
    implicit none
    !
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    logical, optional                   :: print_data
    !
    class(emi_data), pointer            :: cur_data
    integer                             :: istage
    integer                             :: rank

    rank = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit

       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             if (.not. cur_data%is_set) then
                write(iulog,*)'EMID%id   = ', cur_data%id
                write(iulog,*)'EMID%name = ', trim(cur_data%name)
                call endrun(msg='EMID is not set.')
             endif
             if (present(print_data)) then
                if (print_data) then
                   rank = rank + 1
                   call cur_data%PrintInfo(rank)
                   call cur_data%Print()
                endif
             endif
             exit
          endif
       enddo

       cur_data => cur_data%next
    enddo

  end subroutine EMID_Verify_All_Data_Is_Set

end module ExternalModelInterfaceMod
