module EDParamsMod
   !
   ! module that deals with reading the ED parameter file
   !
  
   use FatesParametersInterface, only : param_string_length
   use FatesGlobals        , only : fates_log
   use FatesGlobals        , only : endrun => fates_endrun

   ! CIME Globals
   use shr_log_mod         , only : errMsg => shr_log_errMsg
   use shr_kind_mod        , only: r8 => shr_kind_r8

   implicit none
   save
   ! private - if we allow this module to be private, it does not allow the protected values below to be 
   ! seen outside of this module.

   !
   ! this is what the user can use for the actual values
   !
   
   real(r8),protected :: ED_size_diagnostic_scale             ! Flag to switch between a linear and exponential
                                                              ! scale on the plant size axis in diagnostics (NOT USED YET)
   real(r8),protected :: fates_mortality_disturbance_fraction ! the fraction of canopy mortality that results in disturbance
   real(r8),protected :: ED_val_comp_excln
   real(r8),protected :: ED_val_stress_mort
   real(r8),protected :: ED_val_init_litter
   real(r8),protected :: ED_val_nignitions
   real(r8),protected :: ED_val_understorey_death
   real(r8),protected :: ED_val_cwd_fcel
   real(r8),protected :: ED_val_cwd_flig
   real(r8),protected :: ED_val_bbopt_c3
   real(r8),protected :: ED_val_bbopt_c4
   real(r8),protected :: ED_val_base_mr_20
   real(r8),protected :: ED_val_phen_drought_threshold
   real(r8),protected :: ED_val_phen_doff_time
   real(r8),protected :: ED_val_phen_a
   real(r8),protected :: ED_val_phen_b
   real(r8),protected :: ED_val_phen_c
   real(r8),protected :: ED_val_phen_chiltemp
   real(r8),protected :: ED_val_phen_mindayson
   real(r8),protected :: ED_val_phen_ncolddayslim
   real(r8),protected :: ED_val_phen_coldtemp
   real(r8),protected :: ED_val_cohort_fusion_tol
   real(r8),protected :: ED_val_patch_fusion_tol
   real(r8),protected :: ED_val_canopy_closure_thresh ! site-level canopy closure point where trees take on forest (narrow) versus savannah (wide) crown allometry

   ! two special parameters whose size is defined in the parameter file
   real(r8),protected,allocatable :: ED_val_history_sizeclass_bin_edges(:)
   real(r8),protected,allocatable :: ED_val_history_ageclass_bin_edges(:)

   character(len=param_string_length),parameter :: ED_name_size_diagnostic_scale = "fates_size_diagnostic_scale"
   character(len=param_string_length),parameter :: ED_name_mort_disturb_frac = "fates_mort_disturb_frac"
   character(len=param_string_length),parameter :: ED_name_comp_excln = "fates_comp_excln"
   character(len=param_string_length),parameter :: ED_name_stress_mort = "fates_stress_mort"
   character(len=param_string_length),parameter :: ED_name_init_litter = "fates_init_litter"
   character(len=param_string_length),parameter :: ED_name_nignitions = "fates_nignitions"
   character(len=param_string_length),parameter :: ED_name_understorey_death = "fates_understorey_death"
   character(len=param_string_length),parameter :: ED_name_cwd_fcel= "fates_cwd_fcel"   
   character(len=param_string_length),parameter :: ED_name_cwd_flig= "fates_cwd_flig"   
   character(len=param_string_length),parameter :: ED_name_bbopt_c3= "fates_bbopt_c3"   
   character(len=param_string_length),parameter :: ED_name_bbopt_c4= "fates_bbopt_c4"   
   character(len=param_string_length),parameter :: ED_name_base_mr_20= "fates_base_mr_20"   
   character(len=param_string_length),parameter :: ED_name_phen_drought_threshold= "fates_phen_drought_threshold"   
   character(len=param_string_length),parameter :: ED_name_phen_doff_time= "fates_phen_doff_time"   
   character(len=param_string_length),parameter :: ED_name_phen_a= "fates_phen_a"   
   character(len=param_string_length),parameter :: ED_name_phen_b= "fates_phen_b"   
   character(len=param_string_length),parameter :: ED_name_phen_c= "fates_phen_c"   
   character(len=param_string_length),parameter :: ED_name_phen_chiltemp= "fates_phen_chiltemp"   
   character(len=param_string_length),parameter :: ED_name_phen_mindayson= "fates_phen_mindayson"   
   character(len=param_string_length),parameter :: ED_name_phen_ncolddayslim= "fates_phen_ncolddayslim"   
   character(len=param_string_length),parameter :: ED_name_phen_coldtemp= "fates_phen_coldtemp"   
   character(len=param_string_length),parameter :: ED_name_cohort_fusion_tol= "fates_cohort_fusion_tol"   
   character(len=param_string_length),parameter :: ED_name_patch_fusion_tol= "fates_patch_fusion_tol"
   character(len=param_string_length),parameter :: ED_name_canopy_closure_thresh= "fates_canopy_closure_thresh"      

   ! non-scalar parameter names
   character(len=param_string_length),parameter :: ED_name_history_sizeclass_bin_edges= "fates_history_sizeclass_bin_edges"      
   character(len=param_string_length),parameter :: ED_name_history_ageclass_bin_edges= "fates_history_ageclass_bin_edges"      

   ! Hydraulics Control Parameters (ONLY RELEVANT WHEN USE_FATES_HYDR = TRUE)
   ! ----------------------------------------------------------------------------------------------
   real(r8),protected :: hydr_psi0          !  sapwood water potential at saturation (MPa)
   character(len=param_string_length),parameter :: hydr_name_psi0 = "fates_hydr_psi0"

   real(r8),protected :: hydr_psicap        !  sapwood water potential at which capillary reserves exhausted (MPa)
   character(len=param_string_length),parameter :: hydr_name_psicap = "fates_hydr_psicap"


   ! Logging Control Parameters (ONLY RELEVANT WHEN USE_FATES_LOGGING = TRUE)
   ! ----------------------------------------------------------------------------------------------

   real(r8),protected :: logging_dbhmin              ! Minimum dbh at which logging is applied (cm)
   character(len=param_string_length),parameter :: logging_name_dbhmin = "fates_logging_dbhmin"

   real(r8),protected :: logging_collateral_frac     ! Ratio of collateral mortality to direct logging mortality
   character(len=param_string_length),parameter :: logging_name_collateral_frac = "fates_logging_collateral_frac"

   real(r8),protected :: logging_coll_under_frac ! Fraction of understory plants that die when logging disturbance
                                                 ! is generated
   character(len=param_string_length),parameter :: logging_name_coll_under_frac = "fates_logging_coll_under_frac"
   
   real(r8),protected :: logging_direct_frac         ! Fraction of stems logged per event
   character(len=param_string_length),parameter :: logging_name_direct_frac = "fates_logging_direct_frac"

   real(r8),protected :: logging_mechanical_frac         ! Fraction of stems logged per event
   character(len=param_string_length),parameter :: logging_name_mechanical_frac = "fates_logging_mechanical_frac"

   real(r8),protected :: logging_event_code          ! Code that options how logging events are structured 
   character(len=param_string_length),parameter :: logging_name_event_code = "fates_logging_event_code"
   
   
   public :: FatesParamsInit
   public :: FatesRegisterParams
   public :: FatesReceiveParams
   public :: FatesReportParams
  
contains

  !-----------------------------------------------------------------------
  subroutine FatesParamsInit()
    ! Initialize all parameters to nan to ensure that we get valid
    ! values back from the host.
    
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)

    implicit none

    ED_size_diagnostic_scale              = nan
    fates_mortality_disturbance_fraction  = nan
    ED_val_comp_excln                     = nan
    ED_val_stress_mort                    = nan
    ED_val_init_litter                    = nan
    ED_val_nignitions                     = nan
    ED_val_understorey_death              = nan
    ED_val_cwd_fcel                       = nan
    ED_val_cwd_flig                       = nan
    ED_val_bbopt_c3                       = nan
    ED_val_bbopt_c4                       = nan
    ED_val_base_mr_20                     = nan
    ED_val_phen_drought_threshold         = nan
    ED_val_phen_doff_time                 = nan
    ED_val_phen_a                         = nan
    ED_val_phen_b                         = nan
    ED_val_phen_c                         = nan
    ED_val_phen_chiltemp                  = nan
    ED_val_phen_mindayson                 = nan
    ED_val_phen_ncolddayslim              = nan
    ED_val_phen_coldtemp                  = nan
    ED_val_cohort_fusion_tol              = nan
    ED_val_patch_fusion_tol               = nan
    ED_val_canopy_closure_thresh               = nan    

    hydr_psi0                             = nan
    hydr_psicap                           = nan

    logging_dbhmin                        = nan
    logging_collateral_frac               = nan
    logging_direct_frac                   = nan
    logging_mechanical_frac               = nan
    logging_event_code                    = nan

  end subroutine FatesParamsInit

  !-----------------------------------------------------------------------
  subroutine FatesRegisterParams(fates_params)
    ! Register the parameters we want the host to provide, and
    ! indicate whether they are fates parameters or host parameters
    ! that need to be synced with host values.

    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar1d, dimension_shape_1d
    use FatesParametersInterface, only : dimension_name_history_size_bins, dimension_name_history_age_bins

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    character(len=param_string_length), parameter :: dim_names(1) = (/dimension_name_scalar1d/)
    character(len=param_string_length), parameter :: dim_names_sizeclass(1) = (/dimension_name_history_size_bins/)
    character(len=param_string_length), parameter :: dim_names_ageclass(1) = (/dimension_name_history_age_bins/)

    call FatesParamsInit()

    call fates_params%RegisterParameter(name=ED_name_size_diagnostic_scale, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_mort_disturb_frac, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_comp_excln, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_stress_mort, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_init_litter, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_nignitions, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_understorey_death, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_cwd_fcel, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_cwd_flig, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_bbopt_c3, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_bbopt_c4, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_base_mr_20, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_drought_threshold, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_doff_time, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_a, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_b, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_c, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_chiltemp, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_mindayson, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_ncolddayslim, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_phen_coldtemp, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_cohort_fusion_tol, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_patch_fusion_tol, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=ED_name_canopy_closure_thresh, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=hydr_name_psi0, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=hydr_name_psicap, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=logging_name_dbhmin, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=logging_name_collateral_frac, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=logging_name_coll_under_frac, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=logging_name_direct_frac, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=logging_name_mechanical_frac, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)

    call fates_params%RegisterParameter(name=logging_name_event_code, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names)


    ! non-scalar parameters
    call fates_params%RegisterParameter(name=ED_name_history_sizeclass_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_sizeclass)

    call fates_params%RegisterParameter(name=ED_name_history_ageclass_bin_edges, dimension_shape=dimension_shape_1d, &
         dimension_names=dim_names_ageclass)

  end subroutine FatesRegisterParams

  
  !-----------------------------------------------------------------------
  subroutine FatesReceiveParams(fates_params)
    
    use FatesParametersInterface, only : fates_parameters_type, dimension_name_scalar

    implicit none

    class(fates_parameters_type), intent(inout) :: fates_params

    call fates_params%RetreiveParameter(name=ED_name_size_diagnostic_scale, &
         data=ED_size_diagnostic_scale)

    call fates_params%RetreiveParameter(name=ED_name_mort_disturb_frac, &
          data=fates_mortality_disturbance_fraction)

    call fates_params%RetreiveParameter(name=ED_name_comp_excln, &
         data=ED_val_comp_excln)

    call fates_params%RetreiveParameter(name=ED_name_stress_mort, &
         data=ED_val_stress_mort)

    call fates_params%RetreiveParameter(name=ED_name_init_litter, &
         data=ED_val_init_litter)

    call fates_params%RetreiveParameter(name=ED_name_nignitions, &
         data=ED_val_nignitions)

    call fates_params%RetreiveParameter(name=ED_name_understorey_death, &
         data=ED_val_understorey_death)

    call fates_params%RetreiveParameter(name=ED_name_cwd_fcel, &
         data=ED_val_cwd_fcel)

    call fates_params%RetreiveParameter(name=ED_name_cwd_flig, &
         data=ED_val_cwd_flig)

    call fates_params%RetreiveParameter(name=ED_name_bbopt_c3, &
         data=ED_val_bbopt_c3)

    call fates_params%RetreiveParameter(name=ED_name_bbopt_c4, &
         data=ED_val_bbopt_c4)

    call fates_params%RetreiveParameter(name=ED_name_base_mr_20, &
         data=ED_val_base_mr_20)

    call fates_params%RetreiveParameter(name=ED_name_phen_drought_threshold, &
         data=ED_val_phen_drought_threshold)

    call fates_params%RetreiveParameter(name=ED_name_phen_doff_time, &
         data=ED_val_phen_doff_time)

    call fates_params%RetreiveParameter(name=ED_name_phen_a, &
         data=ED_val_phen_a)

    call fates_params%RetreiveParameter(name=ED_name_phen_b, &
         data=ED_val_phen_b)

    call fates_params%RetreiveParameter(name=ED_name_phen_c, &
         data=ED_val_phen_c)

    call fates_params%RetreiveParameter(name=ED_name_phen_chiltemp, &
         data=ED_val_phen_chiltemp)

    call fates_params%RetreiveParameter(name=ED_name_phen_mindayson, &
         data=ED_val_phen_mindayson)

    call fates_params%RetreiveParameter(name=ED_name_phen_ncolddayslim, &
         data=ED_val_phen_ncolddayslim)

    call fates_params%RetreiveParameter(name=ED_name_phen_coldtemp, &
         data=ED_val_phen_coldtemp)

    call fates_params%RetreiveParameter(name=ED_name_cohort_fusion_tol, &
         data=ED_val_cohort_fusion_tol)

    call fates_params%RetreiveParameter(name=ED_name_patch_fusion_tol, &
         data=ED_val_patch_fusion_tol)
    
    call fates_params%RetreiveParameter(name=ED_name_canopy_closure_thresh, &
         data=ED_val_canopy_closure_thresh)
    
    call fates_params%RetreiveParameter(name=hydr_name_psi0, &
          data=hydr_psi0)

    call fates_params%RetreiveParameter(name=hydr_name_psicap, &
          data=hydr_psicap)

    call fates_params%RetreiveParameter(name=logging_name_dbhmin, &
          data=logging_dbhmin)
    
    call fates_params%RetreiveParameter(name=logging_name_collateral_frac, &
          data=logging_collateral_frac)
    
    call fates_params%RetreiveParameter(name=logging_name_coll_under_frac, &
          data=logging_coll_under_frac)

    call fates_params%RetreiveParameter(name=logging_name_direct_frac, &
          data=logging_direct_frac)

    call fates_params%RetreiveParameter(name=logging_name_mechanical_frac, &
          data=logging_mechanical_frac)
    
    call fates_params%RetreiveParameter(name=logging_name_event_code, &
          data=logging_event_code)

    ! parameters that are arrays of size defined within the params file and thus need allocating as well
    call fates_params%RetreiveParameterAllocate(name=ED_name_history_sizeclass_bin_edges, &
          data=ED_val_history_sizeclass_bin_edges)

    call fates_params%RetreiveParameterAllocate(name=ED_name_history_ageclass_bin_edges, &
          data=ED_val_history_ageclass_bin_edges)


  end subroutine FatesReceiveParams
  
  ! =====================================================================================

  subroutine FatesReportParams(is_master)

     logical,intent(in) :: is_master

     character(len=32),parameter :: fmt0 = '(a,(F12.4))'
     logical, parameter :: debug_report = .false.
     
     if(debug_report .and. is_master) then
        
        write(fates_log(),*) '-----------  FATES Scalar Parameters -----------------'
        write(fates_log(),fmt0) 'ED_size_diagnostic_scale = ',ED_size_diagnostic_scale
        write(fates_log(),fmt0) 'fates_mortality_disturbance_fraction = ',fates_mortality_disturbance_fraction
        write(fates_log(),fmt0) 'ED_val_comp_excln = ',ED_val_comp_excln
        write(fates_log(),fmt0) 'ED_val_stress_mort = ',ED_val_stress_mort
        write(fates_log(),fmt0) 'ED_val_init_litter = ',ED_val_init_litter
        write(fates_log(),fmt0) 'ED_val_nignitions = ',ED_val_nignitions
        write(fates_log(),fmt0) 'ED_val_understorey_death = ',ED_val_understorey_death
        write(fates_log(),fmt0) 'ED_val_cwd_fcel = ',ED_val_cwd_fcel
        write(fates_log(),fmt0) 'ED_val_cwd_flig = ',ED_val_cwd_flig
        write(fates_log(),fmt0) 'ED_val_bbopt_c3 = ',ED_val_bbopt_c3
        write(fates_log(),fmt0) 'ED_val_bbopt_c4 = ',ED_val_bbopt_c4
        write(fates_log(),fmt0) 'ED_val_base_mr_20 = ', ED_val_base_mr_20
        write(fates_log(),fmt0) 'ED_val_phen_drought_threshold = ',ED_val_phen_drought_threshold
        write(fates_log(),fmt0) 'ED_val_phen_doff_time = ',ED_val_phen_doff_time
        write(fates_log(),fmt0) 'ED_val_phen_a = ',ED_val_phen_a
        write(fates_log(),fmt0) 'ED_val_phen_b = ',ED_val_phen_b
        write(fates_log(),fmt0) 'ED_val_phen_c = ',ED_val_phen_c
        write(fates_log(),fmt0) 'ED_val_phen_chiltemp = ',ED_val_phen_chiltemp
        write(fates_log(),fmt0) 'ED_val_phen_mindayson = ',ED_val_phen_mindayson
        write(fates_log(),fmt0) 'ED_val_phen_ncolddayslim = ',ED_val_phen_ncolddayslim
        write(fates_log(),fmt0) 'ED_val_phen_coldtemp = ',ED_val_phen_coldtemp
        write(fates_log(),fmt0) 'ED_val_cohort_fusion_tol = ',ED_val_cohort_fusion_tol
        write(fates_log(),fmt0) 'ED_val_patch_fusion_tol = ',ED_val_patch_fusion_tol
        write(fates_log(),fmt0) 'ED_val_canopy_closure_thresh = ',ED_val_canopy_closure_thresh        
        write(fates_log(),fmt0) 'hydr_psi0 = ',hydr_psi0
        write(fates_log(),fmt0) 'hydr_psicap = ',hydr_psicap
        write(fates_log(),fmt0) 'logging_dbhmin = ',logging_dbhmin
        write(fates_log(),fmt0) 'logging_collateral_frac = ',logging_collateral_frac
        write(fates_log(),fmt0) 'logging_direct_frac = ',logging_direct_frac
        write(fates_log(),fmt0) 'logging_mechanical_frac = ',logging_mechanical_frac
        write(fates_log(),fmt0) 'logging_event_code = ',logging_event_code
        write(fates_log(),*) '------------------------------------------------------'

     end if

  end subroutine FatesReportParams

  
end module EDParamsMod
