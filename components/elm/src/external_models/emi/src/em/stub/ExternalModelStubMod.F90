module ExternalModelStubMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides wrapper for a stub external model
  !
  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use EMI_DataMod, only : emi_data_list, emi_data
  use decompMod                    , only : bounds_type
  use ExternalModelBaseType        , only : em_base_type
  use elm_varctl                            , only : iulog
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_CNCarbonStateType_Constants
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

  integer :: nlevgrnd        = 5 ! This data should also be exchanged
  integer :: nlevdecomp_full = 1 ! This data should also be exchanged
  integer :: ndecomp_pools   = 8 ! This data should also be exchanged

  type, public, extends(em_base_type) :: em_stub_type
     ! ----------------------------------------------------------------------
     ! Indicies required during the initialization
     ! ----------------------------------------------------------------------
     integer :: index_l2e_init_parameter_hksatc
     integer :: index_l2e_init_parameter_bswc

     integer :: index_e2l_init_state_h2osoi_liq
     integer :: index_e2l_init_state_h2osoi_ice

     ! ----------------------------------------------------------------------
     ! Indicies required during timestepping
     ! ----------------------------------------------------------------------
     integer :: index_l2e_flux_et
     integer :: index_l2e_state_h2osoi_liq
     integer :: index_l2e_state_h2osoi_ice
     integer :: index_l2e_flux_hs_soil
     integer :: index_l2e_state_decomp_cpools

     integer :: index_e2l_state_smp
     integer :: index_e2l_state_temperature_soil_nlevgrnd
     integer :: index_e2l_state_decomp_cpools

   contains

     procedure, public :: Populate_L2E_Init_List  => EM_Stub_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EM_Stub_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EM_Stub_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_Stub_Populate_E2L_List
     procedure, public :: Init                    => EM_Stub_Init
     procedure, public :: Solve                   => EM_Stub_Solve
  end type em_stub_type

contains

  !------------------------------------------------------------------------
  subroutine EM_Stub_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by stub EM from ELM during
    ! initialization stage
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_stub_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                         = L2E_PARAMETER_HKSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_hksatc       = index

    id                                         = L2E_PARAMETER_BSWC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_bswc         = index

    deallocate(em_stages)

    write(iulog,*)'    1.1.1 Stub EM will receive following variables to ELM during initalization: '
    call l2e_init_list%PrintInfo()

  end subroutine EM_Stub_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EM_Stub_Populate_E2L_Init_List(this, e2l_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables returned by stub EM to ELM during
    ! initialization stage
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_stub_type)                  :: this
    class(emi_data_list), intent(inout) :: e2l_init_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                        = E2L_STATE_H2OSOI_LIQ
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_liq      = index

    id                                        = E2L_STATE_H2OSOI_ICE
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_ice      = index

    deallocate(em_stages)

    write(iulog,*)'    1.1.2 Stub EM will return following variables to ELM during initalization: '
    call e2l_init_list%PrintInfo()

  end subroutine EM_Stub_Populate_E2L_Init_List

  !------------------------------------------------------------------------
  subroutine EM_Stub_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by stub EM from ELM during
    ! initialization stage
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_stub_type)                  :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index


    ! Request ELM send liq soil moisture and ice soil moisture during
    ! two EM solve stages
    number_em_stages = 2
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_STUB_SOIL_HYDRO_STAGE
    em_stages(2) = EM_STUB_SOIL_THERMAL_STAGE

    id                                   = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq      = index

    id                                   = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice      = index

    deallocate(em_stages)

    ! Request ELM send evapotranspiration during only one EM solve stage
    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_STUB_SOIL_HYDRO_STAGE
    id                                   = L2E_FLUX_VERTICAL_ET_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_et               = index

    ! Request ELM send soil heat flux during only one EM solve stage
    number_em_stages = 1
    em_stages(1) = EM_STUB_SOIL_THERMAL_STAGE
    id                                             = L2E_FLUX_SOIL_HEAT_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_hs_soil                    = index

    number_em_stages = 1
    em_stages(1) = EM_STUB_SOIL_THERMAL_STAGE
    id                                             = L2E_STATE_CARBON_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_decomp_cpools              = index

    deallocate(em_stages)

    write(iulog,*)'    1.1.3 Stub EM will receive following variables from ELM during timestepping: '
    call l2e_list%PrintInfo()

  end subroutine EM_Stub_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_Stub_Populate_E2L_List(this, e2l_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables returned by stub EM to ELM during
    ! initialization stage
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_stub_type)                  :: this
    class(emi_data_list), intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))

    ! Stub EM will send ELM soil matric potential during only one EM solve stage
    number_em_stages = 1
    em_stages(1) = EM_STUB_SOIL_HYDRO_STAGE
    id                              = E2L_STATE_SOIL_MATRIC_POTENTIAL
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_smp        = index

    ! Stub EM will send ELM soil temperature during only one EM solve stage
    number_em_stages = 1
    em_stages(1) = EM_STUB_SOIL_THERMAL_STAGE
    id                                             = E2L_STATE_TSOIL_NLEVGRND
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_temperature_soil_nlevgrnd = index

    number_em_stages = 1
    em_stages(1) = EM_STUB_SOIL_THERMAL_STAGE
    id                                             = E2L_STATE_CARBON_POOLS_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_decomp_cpools              = index

    deallocate(em_stages)

    write(iulog,*)'    1.1.4 Stub EM will return following variables to ELM during timestepping: '
    call e2l_list%PrintInfo()

  end subroutine EM_Stub_Populate_E2L_List

    !------------------------------------------------------------------------
  subroutine EM_Stub_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)

    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_stub_type)                  :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    real(r8)  , pointer                  :: l2e_hksatc(:,:)
    real(r8)  , pointer                  :: l2e_bswc(:,:)
    real(r8)  , pointer                  :: e2l_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2osoi_ice(:,:)
    integer                              :: begc, endc
    integer                              :: c,j

    begc     = bounds_clump%begc
    endc     = bounds_clump%endc

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_hksatc, l2e_hksatc)
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_bswc,   l2e_bswc)

    write(iulog,*)''
    write(iulog,*)'     1.2.2 Value of variables received by Stub EM'

    call Print2DData('l2e_hksatc', begc, endc, 1, nlevgrnd, l2e_hksatc)
    call Print2DData('l2e_bswc'  , begc, endc, 1, nlevgrnd, l2e_bswc)

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_liq, e2l_h2osoi_liq)
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_ice, e2l_h2osoi_ice)

    do c = begc, endc
       do j = 1, nlevgrnd
          e2l_h2osoi_liq(c,j) = (0.1_r8*(c**2._r8) + j)/100._r8
          e2l_h2osoi_ice(c,j) = 1._r8 - (0.1_r8*(c**2._r8) + j)/100._r8
       enddo
    enddo

    write(iulog,*)''
    write(iulog,*)'     1.2.3 Value of variables send by Stub EM'
    call Print2DData('e2l_h2osoi_liq', begc, endc, 1, nlevgrnd, e2l_h2osoi_liq)
    call Print2DData('e2l_h2osoi_ice', begc, endc, 1, nlevgrnd, e2l_h2osoi_ice)

  end subroutine EM_Stub_Init

    !------------------------------------------------------------------------
  subroutine EM_Stub_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)
    !
    ! !DESCRIPTION:
    ! 
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use elm_varctl                , only : iulog
    use ExternalModelConstants    , only : EM_STUB_SOIL_HYDRO_STAGE, EM_STUB_SOIL_THERMAL_STAGE

    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_stub_type)                  :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent (in)   :: bounds_clump

    select case(em_stage)

    case (EM_STUB_SOIL_HYDRO_STAGE)
       call EM_Stub_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, &
            bounds_clump)

    case (EM_STUB_SOIL_THERMAL_STAGE)
       call EM_Stub_Solve_Soil_Thermal(this, em_stage, dt, nstep, l2e_list, e2l_list, &
            bounds_clump)

    case default
       write(iulog,*)'EM_Stub_Solve: Unknown em_stage.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine EM_Stub_Solve

    !------------------------------------------------------------------------
  subroutine EM_Stub_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_stub_type)                  :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    real(r8)  , pointer                  :: l2e_mflux_et(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: e2l_smp_l(:,:)
    integer                              :: begc, endc
    integer                              :: c,j

    begc     = bounds_clump%begc
    endc     = bounds_clump%endc

    write(iulog,*)''
    write(iulog,*)'     2.2 Value of variables received by Stub EM for'
    write(iulog,*)'         solving the hydrology model'

    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_et, l2e_mflux_et)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq, l2e_h2osoi_liq)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice, l2e_h2osoi_ice)

    call Print2DData('l2e_mflux_et', begc, endc, 1, nlevgrnd, l2e_mflux_et)
    call Print2DData('l2e_h2osoi_liq', begc, endc, 1, nlevgrnd, l2e_h2osoi_liq)
    call Print2DData('l2e_h2osoi_ice', begc, endc, 1, nlevgrnd, l2e_h2osoi_ice)

    call e2l_list%GetPointerToReal2D(this%index_e2l_state_smp, e2l_smp_l)

    do c = begc, endc
       do j = 1, nlevgrnd
          e2l_smp_l(c,j) = (2._r8*(c**2._r8) + j)
       enddo
    enddo

    write(iulog,*)''
    write(iulog,*)'     2.3 Value of variables send by Stub EM after'
    write(iulog,*)'         solving the hydrology model'
    call Print2DData('e2l_smp_l', begc, endc, 1, nlevgrnd, e2l_smp_l)
    write(iulog,*)''

  end subroutine EM_Stub_Solve_Soil_Hydro

    !------------------------------------------------------------------------
  subroutine EM_Stub_Solve_Soil_Thermal(this, em_stage, dt, nstep, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_stub_type)                  :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    real(r8) , pointer                   :: l2e_hs_soil(:)
    real(r8)  , pointer                  :: l2e_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: l2e_decomp_cpools(:,:,:)
    real(r8) , pointer                   :: e2l_tsoil(:,:)
    real(r8)  , pointer                  :: e2l_decomp_cpools(:,:,:)
    integer                              :: begc, endc
    integer                              :: c,j,k

    begc     = bounds_clump%begc
    endc     = bounds_clump%endc

    write(iulog,*)''
    write(iulog,*)'     2.2 Value of variables received by Stub EM for'
    write(iulog,*)'         solving the thermal model'

    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_hs_soil, l2e_hs_soil)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq, l2e_h2osoi_liq)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice, l2e_h2osoi_ice)

    call l2e_list%GetPointerToReal3D(this%index_l2e_state_decomp_cpools, l2e_decomp_cpools)

    call Print2DData('l2e_h2osoi_liq', begc, endc, 1, nlevgrnd, l2e_h2osoi_liq)
    call Print2DData('l2e_h2osoi_ice', begc, endc, 1, nlevgrnd, l2e_h2osoi_ice)
    call Print1DData('l2e_hs_soil', begc, endc, l2e_hs_soil)
    call Print3DData('l2e_decomp_cpools', begc, endc, 1, nlevdecomp_full, 1, ndecomp_pools, l2e_decomp_cpools)

    call e2l_list%GetPointerToReal2D(this%index_e2l_state_temperature_soil_nlevgrnd, e2l_tsoil)
    do c = begc, endc
       do j = 1, nlevgrnd
          e2l_tsoil(c,j) = (2._r8*(c**2._r8) + j) + 273.15_r8
       enddo
    enddo

    call e2l_list%GetPointerToReal3D(this%index_e2l_state_decomp_cpools, e2l_decomp_cpools)
    do c = begc, endc
       do j = 1, nlevdecomp_full
          do k = 1, ndecomp_pools
             e2l_decomp_cpools(c,j,k) = l2e_decomp_cpools(c,j,k) + 10.d0
          end do
       enddo
    enddo
    write(iulog,*)''
    write(iulog,*)'     2.3 Value of variables send by Stub EM after'
    write(iulog,*)'         solving the thermal model'
    call Print2DData('e2l_tsoil', begc, endc, 1, nlevgrnd, e2l_tsoil)
    call Print3DData('e2l_decomp_cpools', begc, endc, 1, nlevdecomp_full, 1, ndecomp_pools, e2l_decomp_cpools)
    write(iulog,*)''

  end subroutine EM_Stub_Solve_Soil_Thermal

    !------------------------------------------------------------------------
  subroutine Print1DData(data_name, dim1_beg, dim1_end, data)
    !
    implicit none
    !
    character(len=*) :: data_name
    integer :: dim1_beg, dim1_end
    real(r8)  , pointer :: data(:)
    !
    integer :: i,j

    write(iulog,*)'      variable: ' // trim(data_name) //': '
    do i = dim1_beg, dim1_end
       write(iulog,*)'     ',data(i)
    enddo

  end subroutine Print1DData

    !------------------------------------------------------------------------
  subroutine Print2DData(data_name, dim1_beg, dim1_end, dim2_beg, dim2_end, data)
    !
    implicit none
    !
    character(len=*) :: data_name
    integer :: dim1_beg, dim1_end, dim2_beg, dim2_end
    real(r8)  , pointer :: data(:,:)
    !
    integer :: i,j

    write(iulog,*)'      variable: ' // trim(data_name) //': '
    do i = dim1_beg, dim1_end
       do j = dim2_beg, dim2_end
          write(iulog,*)'     ',data(i,j)
       enddo
    enddo

  end subroutine Print2DData

  !------------------------------------------------------------------------
  subroutine Print3DData(data_name, dim1_beg, dim1_end, dim2_beg, dim2_end, &
       dim3_beg, dim3_end, data)
    !
    implicit none
    !
    character(len=*) :: data_name
    integer :: dim1_beg, dim1_end, dim2_beg, dim2_end, dim3_beg, dim3_end
    real(r8)  , pointer :: data(:,:,:)
    !
    integer :: i,j,k

    write(iulog,*)'      variable: ' // trim(data_name) //': '
    do i = dim1_beg, dim1_end
       do j = dim2_beg, dim2_end
          do k = dim3_beg, dim3_end
             write(iulog,*)'     ',data(i,j,k)
          end do
       enddo
    enddo

  end subroutine Print3DData

end module ExternalModelStubMod
