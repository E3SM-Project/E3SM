module SFParamsMod
   !
   ! module that deals with reading the SF parameter file
   !
   use FatesConstantsMod , only: r8 => fates_r8
   use EDtypesMod        , only: NFSC,NCWD
   use FatesParametersInterface, only : param_string_length

   implicit none
   save
   ! private - if we allow this module to be private, it does not allow the protected values below to be
   ! seen outside of this module.

   !
   ! this is what the user can use for the actual values
   !
   real(r8),protected :: SF_val_fdi_a
   real(r8),protected :: SF_val_fdi_b
   real(r8),protected :: SF_val_fdi_alpha
   real(r8),protected :: SF_val_miner_total
   real(r8),protected :: SF_val_fuel_energy
   real(r8),protected :: SF_val_part_dens
   real(r8),protected :: SF_val_miner_damp
   real(r8),protected :: SF_val_max_durat
   real(r8),protected :: SF_val_durat_slope
   real(r8),protected :: SF_val_wind_max          ! Maximum wind speed expected by fire model (m/min)
   real(r8),protected :: SF_val_alpha_FMC(NFSC)
   real(r8),protected :: SF_val_CWD_frac(NCWD)
   real(r8),protected :: SF_val_max_decomp(NFSC)
   real(r8),protected :: SF_val_SAV(NFSC)
   real(r8),protected :: SF_val_FBD(NFSC)
   real(r8),protected :: SF_val_min_moisture(NFSC)
   real(r8),protected :: SF_val_mid_moisture(NFSC)
   real(r8),protected :: SF_val_low_moisture_Coeff(NFSC)
   real(r8),protected :: SF_val_low_moisture_Slope(NFSC)
   real(r8),protected :: SF_val_mid_moisture_Coeff(NFSC)
   real(r8),protected :: SF_val_mid_moisture_Slope(NFSC)

   character(len=param_string_length),parameter :: SF_name_fdi_a = "fates_fdi_a"
   character(len=param_string_length),parameter :: SF_name_fdi_b = "fates_fdi_b"
   character(len=param_string_length),parameter :: SF_name_fdi_alpha = "fates_fdi_alpha"
   character(len=param_string_length),parameter :: SF_name_miner_total = "fates_miner_total"
   character(len=param_string_length),parameter :: SF_name_fuel_energy = "fates_fuel_energy"
   character(len=param_string_length),parameter :: SF_name_part_dens = "fates_part_dens"
   character(len=param_string_length),parameter :: SF_name_miner_damp = "fates_miner_damp"
   character(len=param_string_length),parameter :: SF_name_max_durat = "fates_max_durat"
   character(len=param_string_length),parameter :: SF_name_durat_slope = "fates_durat_slope"
   character(len=param_string_length),parameter :: SF_name_alpha_FMC = "fates_alpha_FMC"
   character(len=param_string_length),parameter :: SF_name_CWD_frac = "fates_CWD_frac"
   character(len=param_string_length),parameter :: SF_name_max_decomp = "fates_max_decomp"
   character(len=param_string_length),parameter :: SF_name_SAV = "fates_SAV"
   character(len=param_string_length),parameter :: SF_name_FBD = "fates_FBD"
   character(len=param_string_length),parameter :: SF_name_min_moisture = "fates_min_moisture"
   character(len=param_string_length),parameter :: SF_name_mid_moisture = "fates_mid_moisture"
   character(len=param_string_length),parameter :: SF_name_low_moisture_Coeff = "fates_low_moisture_Coeff"
   character(len=param_string_length),parameter :: SF_name_low_moisture_Slope = "fates_low_moisture_Slope"
   character(len=param_string_length),parameter :: SF_name_mid_moisture_Coeff = "fates_mid_moisture_Coeff"
   character(len=param_string_length),parameter :: SF_name_mid_moisture_Slope = "fates_mid_moisture_Slope"
   character(len=param_string_length),parameter :: SF_name_wind_max = "fates_fire_wind_max"

   public :: SpitFireRegisterParams
   public :: SpitFireReceiveParams

   private :: SpitFireParamsInit
   private :: SpitFireRegisterScalars
   private :: SpitFireReceiveScalars
  
   private :: SpitFireRegisterNCWD
   private :: SpitFireReceiveNCWD
  
   private :: SpitFireRegisterNFSC
   private :: SpitFireReceiveNFSC
  
contains
  !-----------------------------------------------------------------------
  subroutine SpitFireParamsInit()
    ! Initialize all parameters to nan to ensure that we get valid
    ! values back from the host.
    
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)

    implicit none

    SF_val_fdi_a = nan
    SF_val_fdi_b = nan
    SF_val_fdi_alpha = nan
    SF_val_miner_total = nan
    SF_val_fuel_energy = nan
    SF_val_part_dens = nan
    SF_val_miner_damp = nan
    SF_val_max_durat = nan
    SF_val_durat_slope = nan
    SF_val_wind_max = nan

    SF_val_CWD_frac(:) = nan

    SF_val_alpha_FMC(:) = nan
    SF_val_max_decomp(:) = nan

    SF_val_SAV(:) = nan
    SF_val_FBD(:) = nan
    SF_val_min_moisture(:) = nan
    SF_val_mid_moisture(:) = nan
    SF_val_low_moisture_Coeff(:) = nan
    SF_val_low_moisture_Slope(:) = nan
    SF_val_mid_moisture_Coeff(:) = nan
    SF_val_mid_moisture_Slope(:) = nan

  end subroutine SpitFireParamsInit

  !-----------------------------------------------------------------------
  subroutine SpitFireRegisterParams(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar, dimension_shape_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    call SpitFireParamsInit()
    call SpitFireRegisterScalars(fates_params)
    call SpitFireRegisterNCWD(fates_params)
    call SpitFireRegisterNFSC(fates_params)

 end subroutine SpitFireRegisterParams

 !-----------------------------------------------------------------------
 subroutine SpitFireReceiveParams(fates_params)
   
   use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar
   
   implicit none
   
    class(fates_parameters_type), intent(inout) :: fates_params
    
    call SpitFireReceiveScalars(fates_params)
    call SpitFireReceiveNCWD(fates_params)
    call SpitFireReceiveNFSC(fates_params)
    
  end subroutine SpitFireReceiveParams

  !-----------------------------------------------------------------------
  subroutine SpitFireRegisterScalars(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar, dimension_shape_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names_scalar(1) = (/dimension_name_scalar/)
    
    call fates_params%RegisterParameter(name=SF_name_wind_max, dimension_shape=dimension_shape_scalar, &
          dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_fdi_a, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_fdi_b, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_fdi_alpha, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_miner_total, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_fuel_energy, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_part_dens, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_miner_damp, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_max_durat, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

    call fates_params%RegisterParameter(name=SF_name_durat_slope, dimension_shape=dimension_shape_scalar, &
         dimension_names=dim_names_scalar)

  end subroutine SpitFireRegisterScalars

 !-----------------------------------------------------------------------
 subroutine SpitFireReceiveScalars(fates_params)
   
    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    call fates_params%RetreiveParameter(name=SF_name_wind_max, &
          data=SF_val_wind_max)

    call fates_params%RetreiveParameter(name=SF_name_fdi_a, &
         data=SF_val_fdi_a)

    call fates_params%RetreiveParameter(name=SF_name_fdi_b, &
         data=SF_val_fdi_b)

    call fates_params%RetreiveParameter(name=SF_name_fdi_alpha, &
         data=SF_val_fdi_alpha)

    call fates_params%RetreiveParameter(name=SF_name_miner_total, &
         data=SF_val_miner_total)

    call fates_params%RetreiveParameter(name=SF_name_fuel_energy, &
         data=SF_val_fuel_energy)

    call fates_params%RetreiveParameter(name=SF_name_part_dens, &
         data=SF_val_part_dens)

    call fates_params%RetreiveParameter(name=SF_name_miner_damp, &
         data=SF_val_miner_damp)

    call fates_params%RetreiveParameter(name=SF_name_max_durat, &
         data=SF_val_max_durat)

    call fates_params%RetreiveParameter(name=SF_name_durat_slope, &
         data=SF_val_durat_slope)

  end subroutine SpitFireReceiveScalars

  !-----------------------------------------------------------------------
  subroutine SpitFireRegisterNCWD(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_cwd, dimension_shape_1d

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names_cwd(1) = (/dimension_name_cwd/)

    call fates_params%RegisterParameter(name=SF_name_CWD_frac, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_cwd)

  end subroutine SpitFireRegisterNCWD

 !-----------------------------------------------------------------------
 subroutine SpitFireReceiveNCWD(fates_params)
   
    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    call fates_params%RetreiveParameter(name=SF_name_CWD_frac, &
         data=SF_val_CWD_frac)

  end subroutine SpitFireReceiveNCWD

  !-----------------------------------------------------------------------
  subroutine SpitFireRegisterNFSC(fates_params)

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_fsc, dimension_shape_1d

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_fsc/)

    call fates_params%RegisterParameter(name=SF_name_SAV, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_FBD, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_min_moisture, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_mid_moisture, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_low_moisture_Coeff, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_low_moisture_Slope, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_mid_moisture_Coeff, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_mid_moisture_Slope, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_alpha_FMC, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=SF_name_max_decomp, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

  end subroutine SpitFireRegisterNFSC

 !-----------------------------------------------------------------------
 subroutine SpitFireReceiveNFSC(fates_params)
   
    use FatesParametersInterface, only : fates_parameters_type

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    
    call fates_params%RetreiveParameter(name=SF_name_SAV, &
         data=SF_val_SAV)

    call fates_params%RetreiveParameter(name=SF_name_FBD, &
         data=SF_val_FBD)

    call fates_params%RetreiveParameter(name=SF_name_min_moisture, &
         data=SF_val_min_moisture)

    call fates_params%RetreiveParameter(name=SF_name_mid_moisture, &
         data=SF_val_mid_moisture)

    call fates_params%RetreiveParameter(name=SF_name_low_moisture_Coeff, &
         data=SF_val_low_moisture_Coeff)

    call fates_params%RetreiveParameter(name=SF_name_low_moisture_Slope, &
         data=SF_val_low_moisture_Slope)

    call fates_params%RetreiveParameter(name=SF_name_mid_moisture_Coeff, &
         data=SF_val_mid_moisture_Coeff)

    call fates_params%RetreiveParameter(name=SF_name_mid_moisture_Slope, &
         data=SF_val_mid_moisture_Slope)

    call fates_params%RetreiveParameter(name=SF_name_alpha_FMC, &
         data=SF_val_alpha_FMC)

    call fates_params%RetreiveParameter(name=SF_name_max_decomp, &
         data=SF_val_max_decomp)

  end subroutine SpitFireReceiveNFSC
  !-----------------------------------------------------------------------


end module SFParamsMod
