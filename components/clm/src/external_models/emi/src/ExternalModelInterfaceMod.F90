module ExternalModelInterfaceMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides an interface to couple ALM with external model
  !
  use shr_kind_mod                         , only : r8 => shr_kind_r8
  use spmdMod                              , only : masterproc, iam
  use shr_log_mod                          , only : errMsg => shr_log_errMsg
  use decompMod                            , only : bounds_type, get_proc_clumps
  use abortutils                           , only : endrun
  use clm_varctl                           , only : iulog
  use ExternalModelInterfaceDataMod        , only : emi_data_list, emi_data
  use ExternalModelIntefaceDataDimensionMod, only : emi_data_dimension_list_type
#ifdef USE_PETSC_LIB
  use ExternalModelVSFMMod                 , only : em_vsfm_type
  use ExternalModelPTMMod                  , only : em_ptm_type
#endif
  use ExternalModelFATESMod                , only : em_fates_type
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
    use clm_varctl, only : use_ed
#ifndef FATES_VIA_EMI
    use clm_varctl, only : use_betr
    use clm_varctl, only : use_pflotran
    use clm_varctl, only : use_vsfm
#endif
    use clm_varctl, only : use_petsc_thermal_model
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
    index_em_vsfm        = 0

    nclumps = get_proc_clumps()

    ! Is FATES active?
    if (use_ed) then
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

    if ( masterproc ) then
       write(iulog,*) 'Number of Exteranl Models = ', num_em
       write(iulog,*) '  BeTR is present     ',(index_em_betr     >0)
       write(iulog,*) '  FATES is present    ',(index_em_fates    >0)
       write(iulog,*) '  PFLOTRAN is present ',(index_em_pflotran >0)
       write(iulog,*) '  VSFM is present     ',(index_em_vsfm     >0)
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
#ifndef FATES_VIA_EMI
    use clm_instMod           , only : soilstate_vars
    use clm_instMod           , only : soilhydrology_vars
    use clm_instMod           , only : waterflux_vars
    use clm_instMod           , only : waterstate_vars
#else
    use clm_instMod           , only : soilstate_inst
    use clm_instMod           , only : soilhydrology_inst
    use clm_instMod           , only : waterflux_inst
    use clm_instMod           , only : waterstate_inst
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
          call EMID_Pack_WaterState_Vars_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, waterstate_vars)
          call EMID_Pack_WaterFlux_Vars_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, waterflux_vars)
          call EMID_Pack_SoilHydrology_Vars_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, soilhydrology_vars)
          call EMID_Pack_SoilState_Vars_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, soilstate_vars)
          call EMID_Pack_Column_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col)
          call EMID_Pack_Landunit_for_EM(l2e_init_list(clump_rank), em_stage, &
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
          call EMID_Unpack_SoilState_Vars_for_EM(e2l_init_list(clump_rank), em_stage, &
               num_e2l_filter_col, e2l_filter_col, soilstate_vars)
          call EMID_Unpack_WaterState_Vars_for_EM(e2l_init_list(clump_rank), em_stage, &
               num_e2l_filter_col, e2l_filter_col, waterstate_vars)
          call EMID_Unpack_WaterFlux_Vars_for_EM(e2l_init_list(clump_rank), em_stage, &
               num_e2l_filter_col, e2l_filter_col, waterflux_vars)
          call EMID_Unpack_SoilHydrology_Vars_for_EM(e2l_init_list(clump_rank), em_stage, &
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
          call EMID_Pack_SoilState_Vars_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col, soilstate_vars)
          call EMID_Pack_Column_for_EM(l2e_init_list(clump_rank), em_stage, &
               num_filter_col, filter_col)
          call EMID_Pack_Landunit_for_EM(l2e_init_list(clump_rank), em_stage, &
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
       canopystate_vars, energyflux_vars)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ExternalModelConstants , only : EM_ID_BETR
    use ExternalModelConstants , only : EM_ID_FATES
    use ExternalModelConstants , only : EM_ID_PFLOTRAN
    use ExternalModelConstants , only : EM_ID_VSFM
    use ExternalModelConstants , only : EM_ID_PTM
    use SoilStateType          , only : soilstate_type
    use SoilHydrologyType      , only : soilhydrology_type
    use TemperatureType        , only : temperature_type
    use WaterFluxType          , only : waterflux_type
    use WaterStateType         , only : waterstate_type
    use atm2lndType            , only : atm2lnd_type
    use CanopyStateType        , only : canopystate_type
    use EnergyFluxType         , only : energyflux_type
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
    !
    integer          :: index_em
    real(r8)         :: dtime
    integer          :: nstep
    type(bounds_type):: bounds_clump
    integer          :: num_filter_col
    integer, pointer :: filter_col(:)
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

       call EMID_Pack_Temperature_Vars_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, temperature_vars)

       elseif (present(num_nolakec_and_nourbanc)  .and. &
               present(filter_nolakec_and_nourbanc)) then

       call EMID_Pack_Temperature_Vars_for_EM(l2e_driver_list(iem), em_stage, &
            num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc, temperature_vars)
    endif

    if ( present(waterstate_vars)) then
       if (present(num_hydrologyc)  .and. &
           present(filter_hydrologyc)) then

          call EMID_Pack_WaterState_Vars_for_EM(l2e_driver_list(iem), em_stage, &
               num_hydrologyc, filter_hydrologyc, waterstate_vars)

       elseif (present(num_nolakec_and_nourbanc)  .and. &
               present(filter_nolakec_and_nourbanc)) then

          call EMID_Pack_WaterState_Vars_for_EM(l2e_driver_list(iem), em_stage, &
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

          call EMID_Pack_WaterState_Vars_for_EM(l2e_driver_list(iem), em_stage, &
               num_filter_col, filter_col, waterstate_vars)
          deallocate(filter_col)
       endif
    endif

    if ( present(waterflux_vars) .and. &
         present(num_hydrologyc) .and. &
         present(filter_hydrologyc)) then

       call EMID_Pack_WaterFlux_Vars_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, waterflux_vars)
    endif

    if ( present(num_nolakec_and_nourbanc) .and. &
         present(filter_nolakec_and_nourbanc)) then

       call EMID_Pack_EnergyFlux_Vars_for_EM(l2e_driver_list(iem), em_stage, &
            num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc, energyflux_vars)

    endif

    if ( present(num_hydrologyc) .and. &
         present(filter_hydrologyc)) then

       call EMID_Pack_Filter_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc)

       call EMID_Pack_Column_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc)

    endif

    if ( present(num_nolakec) .and. &
         present(filter_nolakec)) then

       call EMID_Pack_Filter_for_EM(l2e_driver_list(iem), em_stage, &
            num_nolakec, filter_nolakec)

       call EMID_Pack_Column_for_EM(l2e_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc)

    endif

    if ( present(num_nolakec_and_nourbanc) .and. &
         present(filter_nolakec_and_nourbanc)) then

       call EMID_Pack_Filter_for_EM(l2e_driver_list(iem), em_stage, &
            num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc)

       call EMID_Pack_Column_for_EM(l2e_driver_list(iem), em_stage, &
            num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc)

    endif

    if ( present(num_filter_lun) .and. &
         present(filter_lun)) then

       call EMID_Pack_Landunit_for_EM(l2e_driver_list(iem), em_stage, &
            num_filter_lun, filter_lun)

    endif

    if (present(atm2lnd_vars)) then
       call EMID_Pack_Atm2Land_Forcings_for_EM(l2e_driver_list(iem), em_stage, atm2lnd_vars)
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
    call EMID_Pack_Column_for_EM(l2e_driver_list(iem), em_stage, &
            num_filter_col, filter_col)
    deallocate(filter_col)

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

    case default
       call endrun('Unknown External Model')
    end select

    ! ------------------------------------------------------------------------
    ! Unpack the data for EM
    ! ------------------------------------------------------------------------
    if ( present(waterstate_vars) .and. &
         present(num_hydrologyc)  .and. &
         present(filter_hydrologyc)) then

       call EMID_Unpack_WaterState_Vars_for_EM(e2l_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, waterstate_vars)
    endif

    if ( present(waterflux_vars) .and. &
         present(num_hydrologyc) .and. &
         present(filter_hydrologyc)) then

       call EMID_Unpack_WaterFlux_Vars_for_EM(e2l_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, waterflux_vars)
    endif

    if ( present(soilstate_vars) .and. &
         present(num_hydrologyc) .and. &
         present(filter_hydrologyc)) then

       call EMID_Unpack_SoilState_Vars_for_EM(e2l_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, soilstate_vars)
    endif

    if ( present(soilhydrology_vars) .and. &
         present(num_hydrologyc)     .and. &
         present(filter_hydrologyc)) then

       call EMID_Unpack_SoilHydrology_Vars_for_EM(e2l_driver_list(iem), em_stage, &
            num_hydrologyc, filter_hydrologyc, soilhydrology_vars)
    endif

    if (present(canopystate_vars)) then
       call EMID_Unpack_CanopyState_Vars_for_EM(e2l_driver_list(iem), em_stage, canopystate_vars)
    endif

    if ( present(temperature_vars) .and. &
         present(num_nolakec_and_nourbanc)     .and. &
         present(filter_nolakec_and_nourbanc)) then

       call EMID_Unpack_Temperature_Vars_for_EM(e2l_driver_list(iem), em_stage, &
            num_nolakec_and_nourbanc, filter_nolakec_and_nourbanc, temperature_vars)
    endif

    call EMID_Verify_All_Data_Is_Set(e2l_driver_list(iem), em_stage)

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
  subroutine EMID_Verify_All_Data_Is_Set(data_list, em_stage)
    !
    ! !DESCRIPTION:
    ! Verify that all EMI data that will be exchanged between ALM and external
    ! model for em_stage was set
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
             if (.not. cur_data%is_set) then
                write(iulog,*)'EMID%id   = ', cur_data%id
                write(iulog,*)'EMID%name = ', trim(cur_data%name)
                call endrun(msg='EMID is not set.')
             endif
             exit
          endif
       enddo

       cur_data => cur_data%next
    enddo

  end subroutine EMID_Verify_All_Data_Is_Set

!-----------------------------------------------------------------------
  subroutine EMID_Pack_WaterFlux_Vars_for_EM(data_list, em_stage, &
        num_hydrologyc, filter_hydrologyc, waterflux_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's waterflux_vars for EM
    !
    ! !USES:
    use ExternalModelConstants    , only : L2E_FLUX_INFIL_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_VERTICAL_ET_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_DEW_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_DRAINAGE_MASS_FLUX
    use WaterFluxType             , only : waterflux_type
    use clm_varpar                , only : nlevsoi, nlevgrnd
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer              , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(waterflux_type) , intent(in) :: waterflux_vars
    !
    integer                           :: c,fc,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

#ifdef USE_PETSC_LIB
    associate(&
         mflx_infl_col         => waterflux_vars%mflx_infl_col         , &
         mflx_dew_col          => waterflux_vars%mflx_dew_col          , &
         mflx_snowlyr_disp_col => waterflux_vars%mflx_snowlyr_disp_col , &
         mflx_snowlyr_col      => waterflux_vars%mflx_snowlyr_col      , &
         mflx_sub_snow_col     => waterflux_vars%mflx_sub_snow_col     , &
         mflx_et_col           => waterflux_vars%mflx_et_col           , &
         mflx_drain_col        => waterflux_vars%mflx_drain_col          &
         )
    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_FLUX_VERTICAL_ET_MASS_FLUX)

             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = mflx_et_col(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_INFIL_MASS_FLUX)

             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = mflx_infl_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_DEW_MASS_FLUX)

             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = mflx_dew_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX)

             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = mflx_sub_snow_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX)

             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = mflx_snowlyr_disp_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX)

             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = mflx_snowlyr_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_DRAINAGE_MASS_FLUX)

             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = mflx_drain_col(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate
#endif

  end subroutine EMID_Pack_WaterFlux_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Unpack_WaterFlux_Vars_for_EM(data_list, em_stage, &
        num_hydrologyc, filter_hydrologyc, waterflux_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from EM in ALM's waterflux_vars
    !
    ! !USES:
    use ExternalModelConstants    , only : E2L_FLUX_AQUIFER_RECHARGE
    use ExternalModelConstants    , only : E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    use WaterFluxType             , only : waterflux_type
    use clm_varpar                , only : nlevsoi, nlevgrnd
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer              , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(waterflux_type) , intent(in) :: waterflux_vars
    !
    integer                           :: c,fc,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_unpack
    integer                           :: istage
    integer                           :: count

#ifdef USE_PETSC_LIB
    associate( &
         mflx_recharge_col   => waterflux_vars%mflx_recharge_col, &
         mflx_snowlyr_col    => waterflux_vars%mflx_snowlyr_col   &
    )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_unpack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_unpack = .true.
             exit
          endif
       enddo

       if (need_to_unpack) then

          select case (cur_data%id)

          case (E2L_FLUX_AQUIFER_RECHARGE)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                mflx_recharge_col(c) = cur_data%data_real_1d(c)
             enddo
             cur_data%is_set = .true.

          case (E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                mflx_snowlyr_col(c) = cur_data%data_real_1d(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate
#endif

  end subroutine EMID_Unpack_WaterFlux_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Pack_Filter_for_EM(data_list, em_stage, &
        num_filter, filter)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's filter for EM
    !
    ! !USES:
    use ExternalModelConstants    , only : L2E_FILTER_HYDROLOGYC
    use ExternalModelConstants    , only : L2E_FILTER_NOLAKEC
    use ExternalModelConstants    , only : L2E_FILTER_NOLAKEC_AND_NOURBANC
    use ExternalModelConstants    , only : L2E_FILTER_NUM_HYDROLOGYC
    use ExternalModelConstants    , only : L2E_FILTER_NUM_NOLAKEC
    use ExternalModelConstants    , only : L2E_FILTER_NUM_NOLAKEC_AND_NOURBANC
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_filter ! number of column soil points in column filter
    integer              , intent(in) :: filter(:)  ! column filter for soil points
    !
    integer                           :: i
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage

    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_FILTER_HYDROLOGYC, L2E_FILTER_NOLAKEC, L2E_FILTER_NOLAKEC_AND_NOURBANC)
             do i = 1, num_filter
                cur_data%data_int_1d(i) = filter(i)
             enddo
             cur_data%is_set = .true.

          case (L2E_FILTER_NUM_HYDROLOGYC, L2E_FILTER_NUM_NOLAKEC, L2E_FILTER_NUM_NOLAKEC_AND_NOURBANC)

             cur_data%data_int_1d(1) = num_filter
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

  end subroutine EMID_Pack_Filter_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Pack_Column_for_EM(data_list, em_stage, &
        num_filter, filter)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's column type for EM
    !
    ! !USES:
    use ExternalModelConstants    , only : L2E_COLUMN_ACTIVE
    use ExternalModelConstants    , only : L2E_COLUMN_TYPE
    use ExternalModelConstants    , only : L2E_COLUMN_LANDUNIT_INDEX
    use ExternalModelConstants    , only : L2E_COLUMN_ZI
    use ExternalModelConstants    , only : L2E_COLUMN_DZ
    use ExternalModelConstants    , only : L2E_COLUMN_Z
    use ExternalModelConstants    , only : L2E_COLUMN_AREA
    use ExternalModelConstants    , only : L2E_COLUMN_GRIDCELL_INDEX
    use ExternalModelConstants    , only : L2E_COLUMN_PATCH_INDEX
    use ExternalModelConstants    , only : L2E_COLUMN_NUM_SNOW_LAYERS
    use ExternalModelConstants    , only : L2E_COLUMN_ZI_SNOW_AND_SOIL
    use ExternalModelConstants    , only : L2E_COLUMN_DZ_SNOW_AND_SOIL
    use ExternalModelConstants    , only : L2E_COLUMN_Z_SNOW_AND_SOIL
    use ColumnType                , only : col_pp
    use clm_varpar                , only : nlevgrnd, nlevsno
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_filter ! number of column soil points in column filter
    integer              , intent(in) :: filter(:)  ! column filter for soil points
    !
    integer                           :: fc,c,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_COLUMN_ACTIVE)
             do fc = 1, num_filter
                c = filter(fc)
                if (col_pp%active(c)) cur_data%data_int_1d(c) = 1
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_TYPE)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_int_1d(c) = col_pp%itype(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_LANDUNIT_INDEX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_int_1d(c) = col_pp%landunit(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_ZI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 0, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%zi(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_DZ)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%dz(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_Z)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%z(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_AREA)
             do fc = 1, num_filter
                c = filter(fc)
                !cur_data%data_real_1d(c) = col_pp%area(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_GRIDCELL_INDEX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_int_1d(c) = col_pp%gridcell(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_PATCH_INDEX)
             do fc = 1, num_filter
                c = filter(fc)
#ifndef FATES_VIA_EMI
                cur_data%data_int_1d(c) = col_pp%pfti(c)
#else
                cur_data%data_int_1d(c) = col_pp%patchi(c)
#endif
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_NUM_SNOW_LAYERS)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_int_1d(c) = col_pp%snl(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_ZI_SNOW_AND_SOIL)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%zi(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_DZ_SNOW_AND_SOIL)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno+1, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%dz(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_Z_SNOW_AND_SOIL)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno+1, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%z(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

  end subroutine EMID_Pack_Column_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Pack_Landunit_for_EM(data_list, em_stage, &
        num_filter, filter)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's landunit type for EM
    !
    ! !USES:
    use ExternalModelConstants    , only : L2E_LANDUNIT_TYPE
    use ExternalModelConstants    , only : L2E_LANDUNIT_LAKEPOINT
    use ExternalModelConstants    , only : L2E_LANDUNIT_URBANPOINT
    use LandunitType              , only : lun_pp
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_filter ! number of column soil points in column filter
    integer              , intent(in) :: filter(:)  ! column filter for soil points
    !
    integer                           :: fl,l
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_LANDUNIT_TYPE)
             do fl = 1, num_filter
                l = filter(fl)
                cur_data%data_int_1d(l) = lun_pp%itype(l)
             enddo
             cur_data%is_set = .true.

          case (L2E_LANDUNIT_LAKEPOINT)
             do fl = 1, num_filter
                l = filter(fl)
                if (lun_pp%lakpoi(l)) cur_data%data_int_1d(l) = 1
             enddo
             cur_data%is_set = .true.

          case (L2E_LANDUNIT_URBANPOINT)
             do fl = 1, num_filter
                l = filter(fl)
                if (lun_pp%urbpoi(l)) cur_data%data_int_1d(l) = 1
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

  end subroutine EMID_Pack_Landunit_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Pack_WaterState_Vars_for_EM(data_list, em_stage, &
        num_hydrologyc, filter_hydrologyc, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's waterstate_vars for EM
    !
    ! !USES:
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_LIQ_NLEVGRND
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_ICE_NLEVGRND
    use ExternalModelConstants    , only : L2E_STATE_VSFM_PROGNOSTIC_SOILP
    use ExternalModelConstants    , only : L2E_STATE_FRAC_H2OSFC
    use ExternalModelConstants    , only : L2E_STATE_FRAC_INUNDATED
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_AIR_VOL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_RHO_VAP_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_RHVAP_SOI_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_LIQ_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_ICE_NLEVSOI
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_LIQ_NLEVSNOW
    use ExternalModelConstants    , only : L2E_STATE_H2OSOI_ICE_NLEVSNOW
    use ExternalModelConstants    , only : L2E_STATE_H2OSNOW
    use ExternalModelConstants    , only : L2E_STATE_H2OSFC
    use ExternalModelConstants    , only : L2E_STATE_FRAC_SNOW_EFFECTIVE
    use WaterStateType            , only : waterstate_type
    use clm_varpar                , only : nlevgrnd
    use clm_varpar                , only : nlevsoi
    use clm_varpar                , only : nlevsno
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer              , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(waterstate_type), intent(in) :: waterstate_vars
    !
    integer                           :: c,fc,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

#ifndef FATES_VIA_EMI
    associate(&
         h2osoi_ice    => waterstate_vars%h2osoi_ice_col    , & ! Input:  [real(r8) (:,:) ] ice water (kg/m2)
         h2osoi_liq    => waterstate_vars%h2osoi_liq_col    , & ! Input:  [real(r8) (:,:) ] liquid water (kg/m2)
         soilp_col     => waterstate_vars%soilp_col         , & ! Input:  [real(r8) (:,:) ] soil water pressure (Pa)
         frac_h2osfc   => waterstate_vars%frac_h2osfc_col   , & ! Input:  [real(r8) (:)   ] col fractional area of surface water (-)
         finundated    => waterstate_vars%finundated_col    , & ! Input:  [real(r8) (:)   ] fraction of column that is inundated (-)
         h2osoi_liqvol => waterstate_vars%h2osoi_liqvol_col , & ! Input:  [real(r8) (:,:) ] volumetric liquid water content (m3/m3)
         h2osoi_icevol => waterstate_vars%h2osoi_icevol_col , & ! Input:  [real(r8) (:,:) ] volumetric ice content (m3/m3)
         h2osoi_vol    => waterstate_vars%h2osoi_vol_col    , & ! Input:  [real(r8) (:,:) ] volumetric soil water (m3/m3)
         air_vol       => waterstate_vars%air_vol_col       , & ! Input:  [real(r8) (:,:) ] air filled porosity (m3/m3)
         frac_sno_eff  => waterstate_vars%frac_sno_eff_col  , & ! Input:  [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)
         h2osno        => waterstate_vars%h2osno_col        , & ! Input:  [real(r8) (:)   ] snow water (mm H2O)
         h2osfc        => waterstate_vars%h2osfc_col        , & ! Input:  [real(r8) (:)   ] surface water (mm)
#ifdef BETR_VIA_EMI
         rho_vap       => waterstate_vars%rho_vap_col       , & ! Input:  [real(r8) (:,:) ] water vapor pressure (Pa)
         rhvap_soi     => waterstate_vars%rhvap_soi_col     , & ! Input:  [real(r8) (:,:) ] relative humidity (-)
#endif
         smp_l         => waterstate_vars%smp_l_col           & ! Input:  [real(r8) (:,:) ] soil matric potential (mm)
         )
#else
    associate(&
         h2osoi_ice =>    waterstate_vars%h2osoi_ice_col , & ! Input:  [real(r8) (:,:) ]  ice water (kg/m2)
         h2osoi_liq =>    waterstate_vars%h2osoi_liq_col   & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         )
#endif

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_STATE_H2OSOI_LIQ_NLEVGRND)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = h2osoi_liq(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_ICE_NLEVGRND)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = h2osoi_ice(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

#ifndef FATES_VIA_EMI
          case (L2E_STATE_VSFM_PROGNOSTIC_SOILP)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = soilp_col(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_FRAC_H2OSFC)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = frac_h2osfc(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_FRAC_INUNDATED)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = finundated(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = h2osoi_liqvol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = h2osoi_icevol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_VOL_NLEVSOI)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = h2osoi_vol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_AIR_VOL_NLEVSOI)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = air_vol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

#ifdef BETR_VIA_EMI
          case (L2E_STATE_RHO_VAP_NLEVSOI)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = rho_vap(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_RHVAP_SOI_NLEVSOI)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = rhvap_soi(c,j)
                enddo
             enddo
             cur_data%is_set = .true.
#endif

          case (L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = smp_l(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_LIQ_NLEVSOI)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = h2osoi_liq(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_ICE_NLEVSOI)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = h2osoi_ice(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSNOW)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = h2osno(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSFC)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = h2osfc(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_FRAC_SNOW_EFFECTIVE)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = frac_sno_eff(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_LIQ_NLEVSNOW)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = -nlevsno+1,0
                   cur_data%data_real_2d(c,j) = h2osoi_liq(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_ICE_NLEVSNOW)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = -nlevsno+1,0
                   cur_data%data_real_2d(c,j) = h2osoi_ice(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

#endif

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMID_Pack_WaterState_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Unpack_WaterState_Vars_for_EM(data_list, em_stage, &
        num_hydrologyc, filter_hydrologyc, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Save data from EM in ALM's waterflux_vars
    !
    ! !USES:
    use ExternalModelConstants    , only : E2L_STATE_H2OSOI_LIQ
    use ExternalModelConstants    , only : E2L_STATE_H2OSOI_ICE
    use ExternalModelConstants    , only : E2L_STATE_VSFM_PROGNOSTIC_SOILP
    use WaterStateType            , only : waterstate_type
    use clm_varpar                , only : nlevgrnd
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer              , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(waterstate_type), intent(in) :: waterstate_vars
    !
    integer                           :: c,fc,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_unpack
    integer                           :: istage
    integer                           :: count

#ifndef FATES_VIA_EMI
    associate(&
         h2osoi_ice =>    waterstate_vars%h2osoi_ice_col , & ! Output:  [real(r8) (:,:) ]  ice water (kg/m2)
         h2osoi_liq =>    waterstate_vars%h2osoi_liq_col , & ! Output:  [real(r8) (:,:) ]  liquid water (kg/m2)
         soilp_col  =>    waterstate_vars%soilp_col        & ! Output: [real(r8) (:,:) ]  soil water pressure (Pa)
         )
#else
    associate(&
         h2osoi_ice =>    waterstate_vars%h2osoi_ice_col , & ! Output:  [real(r8) (:,:) ]  ice water (kg/m2)
         h2osoi_liq =>    waterstate_vars%h2osoi_liq_col   & ! Output:  [real(r8) (:,:) ]  liquid water (kg/m2)
         )
#endif

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_unpack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_unpack = .true.
             exit
          endif
       enddo

       if (need_to_unpack) then

          select case (cur_data%id)

          case (E2L_STATE_H2OSOI_LIQ)

             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   h2osoi_liq(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_H2OSOI_ICE)

             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   h2osoi_ice(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

#ifndef FATES_VIA_EMI
          case (E2L_STATE_VSFM_PROGNOSTIC_SOILP)

             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   soilp_col(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.
#endif

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMID_Unpack_WaterState_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Pack_Temperature_Vars_for_EM(data_list, em_stage, &
        num_hydrologyc, filter_hydrologyc, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's temperature_vars for EM
    !
    ! !USES:
    use ExternalModelConstants , only : L2E_STATE_TSOIL_NLEVGRND
    use ExternalModelConstants , only : L2E_STATE_TSNOW
    use ExternalModelConstants , only : L2E_STATE_TH2OSFC
    use TemperatureType        , only : temperature_type
    use clm_varpar             , only : nlevgrnd, nlevsno
    !
    implicit none
    !
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    integer                , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer                , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(temperature_type) , intent(in) :: temperature_vars
    !
    integer                             :: c,fc,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         t_soisno => temperature_vars%t_soisno_col, & ! Input: [real(r8) (:,:) ]  soil temperature (Kelvin)
         t_h2osfc => temperature_vars%t_h2osfc_col  & ! Input: [real(r8) (:)   ]  surface water temperature
         )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_STATE_TSOIL_NLEVGRND)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = t_soisno(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_TSNOW)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = -nlevsno+1, 0
                   cur_data%data_real_2d(c,j) = t_soisno(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_TH2OSFC)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = t_h2osfc(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMID_Pack_Temperature_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Unpack_Temperature_Vars_for_EM(data_list, em_stage, &
        num_filter, filter, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Unpack data from EM into ALM's temperature_vars
    !
    ! !USES:
    use ExternalModelConstants , only : E2L_STATE_TSOIL_NLEVGRND
    use ExternalModelConstants , only : E2L_STATE_TSNOW_NLEVSNOW
    use ExternalModelConstants , only : E2L_STATE_TH2OSFC
    use TemperatureType        , only : temperature_type
    use clm_varpar             , only : nlevgrnd, nlevsno
    !
    implicit none
    !
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    integer                , intent(in) :: num_filter
    integer                , intent(in) :: filter(:)
    type(temperature_type) , intent(in) :: temperature_vars
    !
    integer                             :: c,fc,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_unpack
    integer                             :: istage
    integer                             :: count

    associate(& 
         t_soisno => temperature_vars%t_soisno_col, & ! Input: [real(r8) (:,:) ]  soil temperature (Kelvin)
         t_h2osfc => temperature_vars%t_h2osfc_col  & ! Input: [real(r8) (:)   ]  surface water temperature
         )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_unpack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_unpack = .true.
             exit
          endif
       enddo

       if (need_to_unpack) then

          select case (cur_data%id)

          case (E2L_STATE_TSOIL_NLEVGRND)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   t_soisno(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_TSNOW_NLEVSNOW)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno+1, 0
                   t_soisno(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_TH2OSFC)
             do fc = 1, num_filter
                c = filter(fc)
                t_h2osfc(c) = cur_data%data_real_1d(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMID_Unpack_Temperature_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Pack_SoilHydrology_Vars_for_EM(data_list, em_stage, &
        num_hydrologyc, filter_hydrologyc, soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! Save data from EM in ALM's soilhydrolgy_vars
    !
    ! !USES:
    use ExternalModelConstants    , only : L2E_STATE_WTD
    use SoilHydrologyType         , only : soilhydrology_type
    !
    implicit none
    !
    class(emi_data_list)     , intent(inout) :: data_list
    integer                  , intent(in)    :: em_stage
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    !
    integer                           :: c,fc,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

    associate( &
         zwt => soilhydrology_vars%zwt_col &
    )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_STATE_WTD)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                cur_data%data_real_1d(c) = zwt(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMID_Pack_SoilHydrology_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Unpack_SoilHydrology_Vars_for_EM(data_list, em_stage, &
        num_hydrologyc, filter_hydrologyc, soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! Save data from EM in ALM's soilhydrolgy_vars
    !
    ! !USES:
    use ExternalModelConstants    , only : E2L_STATE_WTD
    use ExternalModelConstants    , only : E2L_FLUX_AQUIFER_RECHARGE
    use SoilHydrologyType         , only : soilhydrology_type
    !
    implicit none
    !
    class(emi_data_list)     , intent(in) :: data_list
    integer                  , intent(in) :: em_stage
    integer                  , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type) , intent(in) :: soilhydrology_vars
    !
    integer                           :: c,fc,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_unpack
    integer                           :: istage
    integer                           :: count

    associate( &
         qcharge            =>    soilhydrology_vars%qcharge_col        , & ! Output:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)
         zwt => soilhydrology_vars%zwt_col &
    )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_unpack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_unpack = .true.
             exit
          endif
       enddo

       if (need_to_unpack) then

          select case (cur_data%id)

          case (E2L_STATE_WTD)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                zwt(c) = cur_data%data_real_1d(c)
             enddo
             cur_data%is_set = .true.

          case (E2L_FLUX_AQUIFER_RECHARGE)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                qcharge(c) = cur_data%data_real_1d(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMID_Unpack_SoilHydrology_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Pack_SoilState_Vars_for_EM(data_list, em_stage, &
        num_hydrologyc, filter_hydrologyc, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! Save data from EM in ALM's soilstate_vars
    !
    ! !USES:
    use ExternalModelConstants    , only : L2E_PARAMETER_WATSATC
    use ExternalModelConstants    , only : L2E_PARAMETER_HKSATC
    use ExternalModelConstants    , only : L2E_PARAMETER_BSWC
    use ExternalModelConstants    , only : L2E_PARAMETER_SUCSATC
    use ExternalModelConstants    , only : L2E_PARAMETER_EFFPOROSITYC
    use ExternalModelConstants    , only : L2E_PARAMETER_CSOL
    use ExternalModelConstants    , only : L2E_PARAMETER_TKMG
    use ExternalModelConstants    , only : L2E_PARAMETER_TKDRY
    use SoilStateType             , only : soilstate_type
    use clm_varpar                , only : nlevsoi, nlevgrnd
    !
    implicit none
    !
    class(emi_data_list) , intent(inout) :: data_list
    integer              , intent(in)    :: em_stage
    integer              , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer              , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilstate_type) , intent(in)    :: soilstate_vars
    !
    integer                           :: c,fc,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

    associate( &
         watsat       => soilstate_vars%watsat_col,       &
         hksat        => soilstate_vars%hksat_col,        &
         bsw          => soilstate_vars%bsw_col,          &
         sucsat       => soilstate_vars%sucsat_col,       &
         eff_porosity => soilstate_vars%eff_porosity_col, &
         csol         => soilstate_vars%csol_col,         &
         tkmg         => soilstate_vars%tkmg_col,         &
         tkdry        => soilstate_vars%tkdry_col         &
    )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_PARAMETER_WATSATC)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = watsat(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_HKSATC)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = hksat(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_BSWC)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = bsw(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_SUCSATC)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = sucsat(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_EFFPOROSITYC)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = eff_porosity(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_CSOL)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = csol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_TKMG)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = tkmg(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_TKDRY)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = tkdry(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMID_Pack_SoilState_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Unpack_SoilState_Vars_for_EM(data_list, em_stage, &
        num_hydrologyc, filter_hydrologyc, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! Save data from EM in ALM's soilstate_vars
    !
    ! !USES:
    use ExternalModelConstants    , only : E2L_STATE_SOIL_MATRIC_POTENTIAL
    use SoilStateType             , only : soilstate_type
    use clm_varpar                , only : nlevsoi, nlevgrnd
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer              , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(soilstate_type) , intent(in) :: soilstate_vars
    !
    integer                           :: c,fc,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_unpack
    integer                           :: istage
    integer                           :: count

    associate( &
         smp_l        => soilstate_vars%smp_l_col &
    )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_unpack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_unpack = .true.
             exit
          endif
       enddo

       if (need_to_unpack) then

          select case (cur_data%id)

          case (E2L_STATE_SOIL_MATRIC_POTENTIAL)
             do fc = 1, num_hydrologyc
                c = filter_hydrologyc(fc)
                do j = 1, nlevgrnd
                   smp_l(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMID_Unpack_SoilState_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Pack_Atm2Land_Forcings_for_EM(data_list, em_stage, atm2lnd_vars)
    !
    ! !DESCRIPTION:
    ! Save data for EM from ALM's atmospheric forcing
    !
    ! !USES:
    use ExternalModelConstants    , only : L2E_FLUX_SOLAR_DIRECT_RADDIATION
    use ExternalModelConstants    , only : L2E_FLUX_SOLAR_DIFFUSE_RADDIATION
    use atm2lndType               , only : atm2lnd_type
    use clm_varpar                , only : nlevsoi, nlevgrnd
    !
    implicit none
    !
    class(emi_data_list), intent(in) :: data_list
    integer             , intent(in) :: em_stage
    type(atm2lnd_type)  , intent(in) :: atm2lnd_vars
    !
    integer                           :: g
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

    associate( &
         forc_solad => atm2lnd_vars%forc_solad_grc, &
         forc_solai => atm2lnd_vars%forc_solai_grc  &
    )

    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_FLUX_SOLAR_DIRECT_RADDIATION)
             do g = cur_data%dim1_beg, cur_data%dim1_end
                cur_data%data_real_2d(g,1:2) = forc_solad(g,1:2)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SOLAR_DIFFUSE_RADDIATION)
             do g = cur_data%dim1_beg, cur_data%dim1_end
                cur_data%data_real_2d(g,1:2) = forc_solai(g,1:2)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMID_Pack_Atm2Land_Forcings_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Unpack_CanopyState_Vars_for_EM(data_list, em_stage, canopystate_vars)
    !
    ! !DESCRIPTION:
    ! Save data for EM from ALM's canopyState_vars
    !
    ! !USES:
    use ExternalModelConstants    , only : E2L_STATE_FSUN
    use ExternalModelConstants    , only : E2L_STATE_LAISUN
    use ExternalModelConstants    , only : E2L_STATE_LAISHA
    use CanopyStateType           , only : canopystate_type
    use clm_varpar                , only : nlevsoi, nlevgrnd
    !
    implicit none
    !
    class(emi_data_list), intent(in) :: data_list
    integer             , intent(in) :: em_stage
    type(canopystate_type), optional ,  intent(inout) :: canopystate_vars
    !
    integer                           :: p
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_unpack
    integer                           :: istage

    associate( &
         fsun       => canopystate_vars%fsun_patch,   &
         laisun     => canopystate_vars%laisun_patch, &
         laisha     => canopystate_vars%laisha_patch  &
    )

    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit

       need_to_unpack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_unpack = .true.
             exit
          endif
       enddo

       if (need_to_unpack) then

          select case (cur_data%id)

          case (E2L_STATE_FSUN)
             do p = cur_data%dim1_beg, cur_data%dim1_end
                fsun(p) = cur_data%data_real_1d(p)
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_LAISUN)
             do p = cur_data%dim1_beg, cur_data%dim1_end
                laisun(p) = cur_data%data_real_1d(p)
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_LAISHA)
             do p = cur_data%dim1_beg, cur_data%dim1_end
                laisha(p) = cur_data%data_real_1d(p)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMID_Unpack_CanopyState_Vars_for_EM

!-----------------------------------------------------------------------
  subroutine EMID_Pack_EnergyFlux_Vars_for_EM(data_list, em_stage, &
        num_filter_col, filter_col, energyflux_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's energyflux_vars for EM
    !
    ! !USES:
    use ExternalModelConstants    , only : L2E_FLUX_ABSORBED_SOLAR_RADIATION
    use ExternalModelConstants    , only : L2E_FLUX_SOIL_HEAT_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_SNOW_HEAT_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_H2OSFC_HEAT_FLUX
    use ExternalModelConstants    , only : L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX
    use EnergyFluxType            , only : energyflux_type
    use clm_varpar                , only : nlevsno
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_filter_col       ! number of points in column filter
    integer              , intent(in) :: filter_col(:)        ! column filter for soil points
    type(energyflux_type), intent(in) :: energyflux_vars
    !
    integer                           :: c,fc,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

#ifdef USE_PETSC_LIB
    associate(&
         sabg_lyr    => energyflux_vars%eflx_sabg_lyr_col    , &
         hs_soil     => energyflux_vars%eflx_hs_soil_col     , &
         hs_top_snow => energyflux_vars%eflx_hs_top_snow_col , &
         hs_h2osfc   => energyflux_vars%eflx_hs_h2osfc_col   , &
         dhsdT       => energyflux_vars%eflx_dhsdT_col         &
         )
    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_FLUX_ABSORBED_SOLAR_RADIATION)

             do fc = 1, num_filter_col
                c = filter_col(fc)
                do j = -nlevsno+1, 1
                   cur_data%data_real_2d(c,j) = sabg_lyr(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SOIL_HEAT_FLUX)

             do fc = 1, num_filter_col
                c = filter_col(fc)
                cur_data%data_real_1d(c) = hs_soil(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SNOW_HEAT_FLUX)

             do fc = 1, num_filter_col
                c = filter_col(fc)
                cur_data%data_real_1d(c) = hs_top_snow(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_H2OSFC_HEAT_FLUX)

             do fc = 1, num_filter_col
                c = filter_col(fc)
                cur_data%data_real_1d(c) = hs_h2osfc(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX)

             do fc = 1, num_filter_col
                c = filter_col(fc)
                cur_data%data_real_1d(c) = dhsdT(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate
#endif

  end subroutine EMID_Pack_EnergyFlux_Vars_for_EM

end module ExternalModelInterfaceMod
