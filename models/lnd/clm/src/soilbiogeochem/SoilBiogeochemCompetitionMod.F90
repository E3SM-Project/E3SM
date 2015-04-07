module SoilBiogeochemCompetitionMod
  
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Resolve plant/heterotroph competition for mineral N
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use clm_varcon                      , only : dzsoi_decomp
  use clm_varctl                      , only : use_nitrif_denitrif
  use abortutils                      , only : endrun
  use decompMod                       , only : bounds_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemNitrogenUptakeMod , only : SoilBiogeochemNitrogenUptake
  use ColumnType                      , only : col                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams
  public :: SoilBiogeochemCompetitionInit         ! Initialization
  public :: SoilBiogeochemCompetition             ! run method

  type :: params_type
     real(r8) :: bdnr              ! bulk denitrification rate (1/s)
     real(r8) :: compet_plant_no3  ! (unitless) relative compettiveness of plants for NO3
     real(r8) :: compet_plant_nh4  ! (unitless) relative compettiveness of plants for NH4
     real(r8) :: compet_decomp_no3 ! (unitless) relative competitiveness of immobilizers for NO3
     real(r8) :: compet_decomp_nh4 ! (unitless) relative competitiveness of immobilizers for NH4
     real(r8) :: compet_denit      ! (unitless) relative competitiveness of denitrifiers for NO3
     real(r8) :: compet_nit        ! (unitless) relative competitiveness of nitrifiers for NH4
  end type params_type
  !
  type(params_type), private :: params_inst  ! params_inst is populated in readParamsMod  
  !
  ! !PUBLIC DATA MEMBERS:
  character(len=* ), public, parameter :: suplnAll='ALL'       ! Supplemental Nitrogen for all PFT's
  character(len=* ), public, parameter :: suplnNon='NONE'      ! No supplemental Nitrogen
  character(len=15), public            :: suplnitro = suplnNon ! Supplemental Nitrogen mode
  !
  ! !PRIVATE DATA MEMBERS:
  real(r8) :: dt   ! decomp timestep (seconds)
  real(r8) :: bdnr ! bulk denitrification rate (1/s)
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io

    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNAllocParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! read in parameters

    tString='bdnr'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%bdnr=tempr

    tString='compet_plant_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%compet_plant_no3=tempr

    tString='compet_plant_nh4'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%compet_plant_nh4=tempr
   
    tString='compet_decomp_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%compet_decomp_no3=tempr

    tString='compet_decomp_nh4'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%compet_decomp_nh4=tempr
   
    tString='compet_denit'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%compet_denit=tempr

    tString='compet_nit'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%compet_nit=tempr   

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemCompetitionInit ( bounds)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use clm_varcon      , only: secspday
    use clm_time_manager, only: get_step_size
    use clm_varctl      , only: iulog, cnallocate_carbon_only_set
    use shr_infnan_mod  , only: nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'SoilBiogeochemCompetitionInit'
    logical :: carbon_only
    !-----------------------------------------------------------------------

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! set space-and-time parameters from parameter file
    bdnr = params_inst%bdnr * (dt/secspday)

    ! Change namelist settings into private logical variables
    select case(suplnitro)
    case(suplnNon)
       carbon_only = .false.
    case(suplnAll)
       carbon_only = .true.
    case default
       write(iulog,*) 'Supplemental Nitrogen flag (suplnitro) can only be: ', &
            suplnNon, ' or ', suplnAll
       call endrun(msg='ERROR: supplemental Nitrogen flag is not correct'//&
            errMsg(__FILE__, __LINE__))
    end select

    call cnallocate_carbon_only_set(carbon_only)

  end subroutine SoilBiogeochemCompetitionInit

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemCompetition (bounds, num_soilc, filter_soilc, &
       soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !USES:
    use clm_varctl , only: cnallocate_carbon_only
    use clm_varpar , only: nlevdecomp, ndecomp_cascade_transitions
    use clm_varcon , only: nitrif_n2o_loss_frac
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds
    integer                                 , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    type(soilbiogeochem_state_type)         , intent(inout) :: soilbiogeochem_state_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,l,pi,j                                            ! indices
    integer  :: fc                                                    ! filter column index
    real(r8) :: compet_plant_no3                                      ! (unitless) relative compettiveness of plants for NO3
    real(r8) :: compet_plant_nh4                                      ! (unitless) relative compettiveness of plants for NH4
    real(r8) :: compet_decomp_no3                                     ! (unitless) relative competitiveness of immobilizers for NO3
    real(r8) :: compet_decomp_nh4                                     ! (unitless) relative competitiveness of immobilizers for NH4
    real(r8) :: compet_denit                                          ! (unitless) relative competitiveness of denitrifiers for NO3
    real(r8) :: compet_nit                                            ! (unitless) relative competitiveness of nitrifiers for NH4
    real(r8) :: fpi_no3_vr(bounds%begc:bounds%endc,1:nlevdecomp)      ! fraction of potential immobilization supplied by no3(no units)
    real(r8) :: fpi_nh4_vr(bounds%begc:bounds%endc,1:nlevdecomp)      ! fraction of potential immobilization supplied by nh4 (no units)
    real(r8) :: sum_nh4_demand(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_nh4_demand_scaled(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_no3_demand(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_no3_demand_scaled(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_ndemand_vr(bounds%begc:bounds%endc, 1:nlevdecomp) !total column N demand (gN/m3/s) at a given level
    real(r8) :: nuptake_prof(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8) :: sminn_tot(bounds%begc:bounds%endc)
    integer  :: nlimit(bounds%begc:bounds%endc,0:nlevdecomp)          !flag for N limitation
    integer  :: nlimit_no3(bounds%begc:bounds%endc,0:nlevdecomp)      !flag for NO3 limitation
    integer  :: nlimit_nh4(bounds%begc:bounds%endc,0:nlevdecomp)      !flag for NH4 limitation
    real(r8) :: residual_sminn_vr(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8) :: residual_sminn(bounds%begc:bounds%endc)
    real(r8) :: residual_smin_nh4_vr(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8) :: residual_smin_no3_vr(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8) :: residual_smin_nh4(bounds%begc:bounds%endc)
    real(r8) :: residual_smin_no3(bounds%begc:bounds%endc)
    real(r8) :: residual_plant_ndemand(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------

    associate(                                                                                           &
         fpg                          => soilbiogeochem_state_inst%fpg_col                             , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)    
         fpi                          => soilbiogeochem_state_inst%fpi_col                             , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         fpi_vr                       => soilbiogeochem_state_inst%fpi_vr_col                          , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
         nfixation_prof               => soilbiogeochem_state_inst%nfixation_prof_col                  , & ! Output: [real(r8) (:,:) ]                                        
         plant_ndemand                => soilbiogeochem_state_inst%plant_ndemand_col                   , & ! Input:  [real(r8) (:)   ]  column-level plant N demand

         sminn_vr                     => soilbiogeochem_nitrogenstate_inst%sminn_vr_col                , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N                
         smin_nh4_vr                  => soilbiogeochem_nitrogenstate_inst%smin_nh4_vr_col             , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral NH4              
         smin_no3_vr                  => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col             , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral NO3              

         pot_f_nit_vr                 => soilbiogeochem_nitrogenflux_inst%pot_f_nit_vr_col             , & ! Input:  [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux
         pot_f_denit_vr               => soilbiogeochem_nitrogenflux_inst%pot_f_denit_vr_col           , & ! Input:  [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux
         f_nit_vr                     => soilbiogeochem_nitrogenflux_inst%f_nit_vr_col                 , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil nitrification flux     
         f_denit_vr                   => soilbiogeochem_nitrogenflux_inst%f_denit_vr_col               , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil denitrification flux   
         potential_immob              => soilbiogeochem_nitrogenflux_inst%potential_immob_col          , & ! Output: [real(r8) (:)   ]                                          
         actual_immob                 => soilbiogeochem_nitrogenflux_inst%actual_immob_col             , & ! Output: [real(r8) (:)   ]                                          
         sminn_to_plant               => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_col           , & ! Output: [real(r8) (:)   ]                                          
         sminn_to_denit_excess_vr     => soilbiogeochem_nitrogenflux_inst%sminn_to_denit_excess_vr_col , & ! Output: [real(r8) (:,:) ]                                        
         actual_immob_no3_vr          => soilbiogeochem_nitrogenflux_inst%actual_immob_no3_vr_col      , & ! Output: [real(r8) (:,:) ]                                        
         actual_immob_nh4_vr          => soilbiogeochem_nitrogenflux_inst%actual_immob_nh4_vr_col      , & ! Output: [real(r8) (:,:) ]                                        
         smin_no3_to_plant_vr         => soilbiogeochem_nitrogenflux_inst%smin_no3_to_plant_vr_col     , & ! Output: [real(r8) (:,:) ]                                        
         smin_nh4_to_plant_vr         => soilbiogeochem_nitrogenflux_inst%smin_nh4_to_plant_vr_col     , & ! Output: [real(r8) (:,:) ]                                        
         n2_n2o_ratio_denit_vr        => soilbiogeochem_nitrogenflux_inst%n2_n2o_ratio_denit_vr_col    , & ! Output: [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
         f_n2o_denit_vr               => soilbiogeochem_nitrogenflux_inst%f_n2o_denit_vr_col           , & ! Output: [real(r8) (:,:) ]  flux of N2O from denitrification [gN/m3/s]
         f_n2o_nit_vr                 => soilbiogeochem_nitrogenflux_inst%f_n2o_nit_vr_col             , & ! Output: [real(r8) (:,:) ]  flux of N2O from nitrification [gN/m3/s]
         supplement_to_sminn_vr       => soilbiogeochem_nitrogenflux_inst%supplement_to_sminn_vr_col   , & ! Output: [real(r8) (:,:) ]                                        
         sminn_to_plant_vr            => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_vr_col        , & ! Output: [real(r8) (:,:) ]                                        
         potential_immob_vr           => soilbiogeochem_nitrogenflux_inst%potential_immob_vr_col       , & ! Input:  [real(r8) (:,:) ]                                        
         actual_immob_vr              => soilbiogeochem_nitrogenflux_inst%actual_immob_vr_col            & ! Output: [real(r8) (:,:) ]                                        
         )

      ! calcualte nitrogen uptake profile
      ! nuptake_prof(:,:) = nan
      ! call SoilBiogelchemNitrogenUptakeProfile(bounds, &
      !     nlevdecomp, num_soilc, filter_soilc, &
      !     sminn_vr, dzsoi_decomp, nfixation_prof, nuptake_prof)

      ! column loops to resolve plant/heterotroph competition for mineral N

      if (.not. use_nitrif_denitrif) then

         ! init sminn_tot
         do fc=1,num_soilc
            c = filter_soilc(fc)
            sminn_tot(c) = 0.
         end do

         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_tot(c) = sminn_tot(c) + sminn_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)      
               if (sminn_tot(c)  >  0.) then
                  nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
               else
                  nuptake_prof(c,j) = nfixation_prof(c,j)
               endif
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)      
               sum_ndemand_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j) + potential_immob_vr(c,j)
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)      
               l = col%landunit(c)
               if (sum_ndemand_vr(c,j)*dt < sminn_vr(c,j)) then

                  ! N availability is not limiting immobilization or plant
                  ! uptake, and both can proceed at their potential rates
                  nlimit(c,j) = 0
                  fpi_vr(c,j) = 1.0_r8
                  actual_immob_vr(c,j) = potential_immob_vr(c,j)
                  sminn_to_plant_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j)
               else if ( cnallocate_carbon_only()) then !.or. &
                  ! this code block controls the addition of N to sminn pool
                  ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
                  ! model behave essentially as a carbon-only model, but with the
                  ! benefit of keeping track of the N additions needed to
                  ! eliminate N limitations, so there is still a diagnostic quantity
                  ! that describes the degree of N limitation at steady-state.

                  nlimit(c,j) = 1
                  fpi_vr(c,j) = 1.0_r8
                  actual_immob_vr(c,j) = potential_immob_vr(c,j)
                  sminn_to_plant_vr(c,j) =  plant_ndemand(c) * nuptake_prof(c,j)
                  supplement_to_sminn_vr(c,j) = sum_ndemand_vr(c,j) - (sminn_vr(c,j)/dt)
               else
                  ! N availability can not satisfy the sum of immobilization and
                  ! plant growth demands, so these two demands compete for available
                  ! soil mineral N resource.

                  nlimit(c,j) = 1
                  if (sum_ndemand_vr(c,j) > 0.0_r8) then
                     actual_immob_vr(c,j) = (sminn_vr(c,j)/dt)*(potential_immob_vr(c,j) / sum_ndemand_vr(c,j))
                  else
                     actual_immob_vr(c,j) = 0.0_r8
                  end if

                  if (potential_immob_vr(c,j) > 0.0_r8) then
                     fpi_vr(c,j) = actual_immob_vr(c,j) / potential_immob_vr(c,j)
                  else
                     fpi_vr(c,j) = 0.0_r8
                  end if

                  sminn_to_plant_vr(c,j) = (sminn_vr(c,j)/dt) - actual_immob_vr(c,j)
               end if
            end do
         end do

         ! sum up N fluxes to plant
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         ! give plants a second pass to see if there is any mineral N left over with which to satisfy residual N demand.
         do fc=1,num_soilc
            c = filter_soilc(fc)    
            residual_sminn(c) = 0._r8
         end do

         ! sum up total N left over after initial plant and immobilization fluxes
         do fc=1,num_soilc
            c = filter_soilc(fc)    
            residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
         end do
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               if (residual_plant_ndemand(c)  >  0._r8 ) then
                  if (nlimit(c,j) .eq. 0) then
                     residual_sminn_vr(c,j) = max(sminn_vr(c,j) - (actual_immob_vr(c,j) + sminn_to_plant_vr(c,j) ) * dt, 0._r8)
                     residual_sminn(c) = residual_sminn(c) + residual_sminn_vr(c,j) * dzsoi_decomp(j)
                  else
                     residual_sminn_vr(c,j)  = 0._r8
                  endif
               endif
            end do
         end do

         ! distribute residual N to plants
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               if ( residual_plant_ndemand(c)  >  0._r8 .and. residual_sminn(c)  >  0._r8 .and. nlimit(c,j) .eq. 0) then
                  sminn_to_plant_vr(c,j) = sminn_to_plant_vr(c,j) + residual_sminn_vr(c,j) * &
                       min(( residual_plant_ndemand(c) *  dt ) / residual_sminn(c), 1._r8) / dt
               endif
            end do
         end do

         ! re-sum up N fluxes to plant
         do fc=1,num_soilc
            c = filter_soilc(fc)    
            sminn_to_plant(c) = 0._r8
         end do
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
               sum_ndemand_vr(c,j) = potential_immob_vr(c,j) + sminn_to_plant_vr(c,j)
            end do
         end do

         ! under conditions of excess N, some proportion is assumed to
         ! be lost to denitrification, in addition to the constant
         ! proportion lost in the decomposition pathways
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               if ((sminn_to_plant_vr(c,j) + actual_immob_vr(c,j))*dt < sminn_vr(c,j)) then
                  sminn_to_denit_excess_vr(c,j) = max(bdnr*((sminn_vr(c,j)/dt) - sum_ndemand_vr(c,j)),0._r8)
               else
                  sminn_to_denit_excess_vr(c,j) = 0._r8
               endif
            end do
         end do

         ! sum up N fluxes to immobilization
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)    
               actual_immob(c) = actual_immob(c) + actual_immob_vr(c,j) * dzsoi_decomp(j)
               potential_immob(c) = potential_immob(c) + potential_immob_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         do fc=1,num_soilc
            c = filter_soilc(fc)    
            ! calculate the fraction of potential growth that can be
            ! acheived with the N available to plants      
            if (plant_ndemand(c) > 0.0_r8) then
               fpg(c) = sminn_to_plant(c) / plant_ndemand(c)
            else
               fpg(c) = 1.0_r8
            end if

            ! calculate the fraction of immobilization realized (for diagnostic purposes)
            if (potential_immob(c) > 0.0_r8) then
               fpi(c) = actual_immob(c) / potential_immob(c)
            else
               fpi(c) = 1.0_r8
            end if
         end do

      else  !----------NITRIF_DENITRIF-------------!

         ! column loops to resolve plant/heterotroph/nitrifier/denitrifier competition for mineral N
         !read constants from external netcdf file
         compet_plant_no3  = params_inst%compet_plant_no3
         compet_plant_nh4  = params_inst%compet_plant_nh4
         compet_decomp_no3 = params_inst%compet_decomp_no3
         compet_decomp_nh4 = params_inst%compet_decomp_nh4
         compet_denit      = params_inst%compet_denit
         compet_nit        = params_inst%compet_nit

         ! init total mineral N pools
         do fc=1,num_soilc
            c = filter_soilc(fc)
            sminn_tot(c) = 0.
         end do

         ! sum up total mineral N pools
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_tot(c) = sminn_tot(c) + (smin_no3_vr(c,j) + smin_nh4_vr(c,j)) * dzsoi_decomp(j)
            end do
         end do

         ! define N uptake profile for initial vertical distribution of plant N uptake, assuming plant seeks N from where it is most abundant
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)
               if (sminn_tot(c)  >  0.) then
                  nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
               else
                  nuptake_prof(c,j) = nfixation_prof(c,j)
               endif
            end do
         end do

         ! main column/vertical loop
         do j = 1, nlevdecomp  
            do fc=1,num_soilc
               c = filter_soilc(fc)
               l = col%landunit(c)

               !  first compete for nh4
               sum_nh4_demand(c,j) = plant_ndemand(c) * nuptake_prof(c,j) + potential_immob_vr(c,j) + pot_f_nit_vr(c,j)
               sum_nh4_demand_scaled(c,j) = plant_ndemand(c)* nuptake_prof(c,j) * compet_plant_nh4 + &
                    potential_immob_vr(c,j)*compet_decomp_nh4 + pot_f_nit_vr(c,j)*compet_nit

               if (sum_nh4_demand(c,j)*dt < smin_nh4_vr(c,j)) then

                  ! NH4 availability is not limiting immobilization or plant
                  ! uptake, and all can proceed at their potential rates
                  nlimit_nh4(c,j) = 0
                  fpi_nh4_vr(c,j) = 1.0_r8
                  actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j)
                  smin_nh4_to_plant_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j)

                  f_nit_vr(c,j) = pot_f_nit_vr(c,j)

               else

                  ! NH4 availability can not satisfy the sum of immobilization, nitrification, and
                  ! plant growth demands, so these three demands compete for available
                  ! soil mineral NH4 resource.
                  nlimit_nh4(c,j) = 1
                  if (sum_nh4_demand(c,j) > 0.0_r8) then
                     actual_immob_nh4_vr(c,j) = min((smin_nh4_vr(c,j)/dt)*(potential_immob_vr(c,j)* &
                          compet_decomp_nh4 / sum_nh4_demand_scaled(c,j)), potential_immob_vr(c,j))
                     smin_nh4_to_plant_vr(c,j) = min((smin_nh4_vr(c,j)/dt)*(plant_ndemand(c)* &
                          nuptake_prof(c,j)*compet_plant_nh4 / sum_nh4_demand_scaled(c,j)), plant_ndemand(c)*nuptake_prof(c,j))
                     f_nit_vr(c,j) =  min((smin_nh4_vr(c,j)/dt)*(pot_f_nit_vr(c,j)*compet_nit / &
                          sum_nh4_demand_scaled(c,j)), pot_f_nit_vr(c,j))
                  else
                     actual_immob_nh4_vr(c,j) = 0.0_r8
                     smin_nh4_to_plant_vr(c,j) = 0.0_r8
                     f_nit_vr(c,j) = 0.0_r8
                  end if

                  if (potential_immob_vr(c,j) > 0.0_r8) then
                     fpi_nh4_vr(c,j) = actual_immob_nh4_vr(c,j) / potential_immob_vr(c,j)
                  else
                     fpi_nh4_vr(c,j) = 0.0_r8
                  end if

               end if

               ! next compete for no3
               sum_no3_demand(c,j) = (plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j)) + &
                    (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j)) + pot_f_denit_vr(c,j)
               sum_no3_demand_scaled(c,j) = (plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))*compet_plant_no3 + &
                    (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))*compet_decomp_no3 + pot_f_denit_vr(c,j)*compet_denit

               if (sum_no3_demand(c,j)*dt < smin_no3_vr(c,j)) then

                  ! NO3 availability is not limiting immobilization or plant
                  ! uptake, and all can proceed at their potential rates
                  nlimit_no3(c,j) = 1
                  fpi_no3_vr(c,j) = 1.0_r8 -  fpi_nh4_vr(c,j)
                  actual_immob_no3_vr(c,j) = (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))
                  smin_no3_to_plant_vr(c,j) = (plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))

                  f_denit_vr(c,j) = pot_f_denit_vr(c,j)

               else 

                  ! NO3 availability can not satisfy the sum of immobilization, denitrification, and
                  ! plant growth demands, so these three demands compete for available
                  ! soil mineral NO3 resource.
                  nlimit_no3(c,j) = 1
                  if (sum_no3_demand(c,j) > 0.0_r8) then
                     actual_immob_no3_vr(c,j) = min((smin_no3_vr(c,j)/dt)*((potential_immob_vr(c,j)- &
                          actual_immob_nh4_vr(c,j))*compet_decomp_no3 / sum_no3_demand_scaled(c,j)), &
                          potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))
                     smin_no3_to_plant_vr(c,j) = min((smin_no3_vr(c,j)/dt)*((plant_ndemand(c)* &
                          nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))*compet_plant_no3 / sum_no3_demand_scaled(c,j)), &
                          plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))
                     f_denit_vr(c,j) =  min((smin_no3_vr(c,j)/dt)*(pot_f_denit_vr(c,j)*compet_denit / &
                          sum_no3_demand_scaled(c,j)), pot_f_denit_vr(c,j))
                  else
                     actual_immob_no3_vr(c,j) = 0.0_r8
                     smin_no3_to_plant_vr(c,j) = 0.0_r8
                     f_denit_vr(c,j) = 0.0_r8
                  end if

                  if (potential_immob_vr(c,j) > 0.0_r8) then
                     fpi_no3_vr(c,j) = actual_immob_no3_vr(c,j) / potential_immob_vr(c,j)
                  else
                     fpi_no3_vr(c,j) = 0.0_r8
                  end if

               end if

               ! n2o emissions: n2o from nitr is const fraction, n2o from denitr is calculated in nitrif_denitrif
               f_n2o_nit_vr(c,j) = f_nit_vr(c,j) * nitrif_n2o_loss_frac
               f_n2o_denit_vr(c,j) = f_denit_vr(c,j) / (1._r8 + n2_n2o_ratio_denit_vr(c,j))


               ! this code block controls the addition of N to sminn pool
               ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
               ! model behave essentially as a carbon-only model, but with the
               ! benefit of keeping track of the N additions needed to
               ! eliminate N limitations, so there is still a diagnostic quantity
               ! that describes the degree of N limitation at steady-state.

               if ( cnallocate_carbon_only()) then !.or. &
                  if ( fpi_no3_vr(c,j) + fpi_nh4_vr(c,j) < 1._r8 ) then
                     fpi_nh4_vr(c,j) = 1.0_r8 - fpi_no3_vr(c,j)
                     supplement_to_sminn_vr(c,j) = (potential_immob_vr(c,j) - actual_immob_no3_vr(c,j)) - actual_immob_nh4_vr(c,j)
                     ! update to new values that satisfy demand
                     actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j) -  actual_immob_no3_vr(c,j)   
                  end if
                  if ( smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j) < plant_ndemand(c)*nuptake_prof(c,j) ) then
                     supplement_to_sminn_vr(c,j) = supplement_to_sminn_vr(c,j) + &
                          (plant_ndemand(c)*nuptake_prof(c,j) - smin_no3_to_plant_vr(c,j)) - smin_nh4_to_plant_vr(c,j)  ! use old values
                     smin_nh4_to_plant_vr(c,j) = plant_ndemand(c)*nuptake_prof(c,j) - smin_no3_to_plant_vr(c,j)
                  end if
                  sminn_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)
               end if

               ! sum up no3 and nh4 fluxes
               fpi_vr(c,j) = fpi_no3_vr(c,j) + fpi_nh4_vr(c,j)
               sminn_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)
               actual_immob_vr(c,j) = actual_immob_no3_vr(c,j) + actual_immob_nh4_vr(c,j)
            end do
         end do

         do fc=1,num_soilc
            c = filter_soilc(fc)
            ! sum up N fluxes to plant after initial competition
            sminn_to_plant(c) = 0._r8
         end do
         do j = 1, nlevdecomp  
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         ! give plants a second pass to see if there is any mineral N left over with which to satisfy residual N demand.
         ! first take frm nh4 pool; then take from no3 pool
         do fc=1,num_soilc
            c = filter_soilc(fc)
            residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
            residual_smin_nh4(c) = 0._r8
         end do
         do j = 1, nlevdecomp  
            do fc=1,num_soilc
               c = filter_soilc(fc)
               if (residual_plant_ndemand(c)  >  0._r8 ) then
                  if (nlimit_nh4(c,j) .eq. 0) then
                     residual_smin_nh4_vr(c,j) = max(smin_nh4_vr(c,j) - (actual_immob_vr(c,j) + &
                          smin_nh4_to_plant_vr(c,j) ) * dt, 0._r8)
                     residual_smin_nh4(c) = residual_smin_nh4(c) + residual_smin_nh4_vr(c,j) * dzsoi_decomp(j)
                  else
                     residual_smin_nh4_vr(c,j)  = 0._r8
                  endif

                  if ( residual_smin_nh4(c) > 0._r8 .and. nlimit_nh4(c,j) .eq. 0 ) then
                     smin_nh4_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + residual_smin_nh4_vr(c,j) * &
                          min(( residual_plant_ndemand(c) *  dt ) / residual_smin_nh4(c), 1._r8) / dt
                  endif
               end if
            end do
         end do

         ! re-sum up N fluxes to plant after second pass for nh4
         do fc=1,num_soilc
            c = filter_soilc(fc)
            sminn_to_plant(c) = 0._r8
         end do
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)
               sminn_to_plant(c) = sminn_to_plant(c) + (sminn_to_plant_vr(c,j)) * dzsoi_decomp(j)
            end do
         end do

         !
         ! and now do second pass for no3
         do fc=1,num_soilc
            c = filter_soilc(fc)
            residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
            residual_smin_no3(c) = 0._r8
         end do

         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)
               if (residual_plant_ndemand(c) > 0._r8 ) then
                  if (nlimit_no3(c,j) .eq. 0) then
                     residual_smin_no3_vr(c,j) = max(smin_no3_vr(c,j) - (actual_immob_vr(c,j) + &
                          smin_no3_to_plant_vr(c,j) ) * dt, 0._r8)
                     residual_smin_no3(c) = residual_smin_no3(c) + residual_smin_no3_vr(c,j) * dzsoi_decomp(j)
                  else
                     residual_smin_no3_vr(c,j)  = 0._r8
                  endif

                  if ( residual_smin_no3(c) > 0._r8 .and. nlimit_no3(c,j) .eq. 0) then
                     smin_no3_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + residual_smin_no3_vr(c,j) * &
                          min(( residual_plant_ndemand(c) *  dt ) / residual_smin_no3(c), 1._r8) / dt
                  endif
               endif
            end do
         end do

         ! re-sum up N fluxes to plant after second passes of both no3 and nh4
         do fc=1,num_soilc
            c = filter_soilc(fc)
            sminn_to_plant(c) = 0._r8
         end do
         do j = 1, nlevdecomp
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)
               sminn_to_plant(c) = sminn_to_plant(c) + (sminn_to_plant_vr(c,j)) * dzsoi_decomp(j)
            end do
         end do

         ! sum up N fluxes to immobilization
         do fc=1,num_soilc
            c = filter_soilc(fc)
            actual_immob(c) = 0._r8
            potential_immob(c) = 0._r8
         end do
         do j = 1, nlevdecomp  
            do fc=1,num_soilc
               c = filter_soilc(fc)
               actual_immob(c) = actual_immob(c) + actual_immob_vr(c,j) * dzsoi_decomp(j)
               potential_immob(c) = potential_immob(c) + potential_immob_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         do fc=1,num_soilc
            c = filter_soilc(fc)   
            ! calculate the fraction of potential growth that can be
            ! acheived with the N available to plants
            if (plant_ndemand(c) > 0.0_r8) then
               fpg(c) = sminn_to_plant(c) / plant_ndemand(c)
            else
               fpg(c) = 1._r8
            end if

            ! calculate the fraction of immobilization realized (for diagnostic purposes)
            if (potential_immob(c) > 0.0_r8) then
               fpi(c) = actual_immob(c) / potential_immob(c)
            else
               fpi(c) = 1._r8
            end if
         end do ! end of column loops

      end if  !end of if_not_use_nitrif_denitrif

    end associate

  end subroutine SoilBiogeochemCompetition

end module SoilBiogeochemCompetitionMod
