module CNNDynamicsMod

  !-----------------------------------------------------------------------
  ! !MODULE: CNNDynamicsMod
  !
  ! !DESCRIPTION:
  ! Module for mineral nitrogen dynamics (deposition, fixation, leaching)
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use decompMod           , only : bounds_type
  use clm_varcon          , only : dzsoi_decomp, zisoi
  use clm_varctl          , only : use_nitrif_denitrif, use_vertsoilc
  use subgridAveMod       , only : p2c
  use atm2lndType         , only : atm2lnd_type
  use CNCarbonFluxType    , only : carbonflux_type, nfix_timeconst
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNStateType         , only : cnstate_type
  use WaterStateType      , only : waterstate_type
  use WaterFluxType       , only : waterflux_type
  use CropType            , only : crop_type
  use ColumnType          , only : col                
  use PatchType           , only : pft                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNNDynamicsInit
  public :: CNNDeposition
  public :: CNNFixation
  public :: CNNLeaching
  public :: CNNFert
  public :: CNSoyfix
  public :: readCNNDynamicsParams
  !
  ! !PRIVATE DATA:
  type, private :: CNNDynamicsParamsType
     real(r8):: sf        ! soluble fraction of mineral N (unitless)
     real(r8):: sf_no3    ! soluble fraction of NO3 (unitless)
  end type CNNDynamicsParamsType
  
  type(CNNDynamicsParamsType),private ::  CNNDynamicsParamsInst
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNNDynamicsInit ( )
    !
    ! !DESCRIPTION:
    ! Initialize module variables
    !-----------------------------------------------------------------------

    if (nfix_timeconst .eq. -1.2345_r8) then
       ! If nfix_timeconst is equal to the junk default value, then it
       ! wasn't specified by the user namelist and we need to assign
       ! it the correct default value. If the user specified it in the
       ! name list, we leave it alone.
       if (use_nitrif_denitrif) then
          nfix_timeconst = 10._r8
       else
          nfix_timeconst = 0._r8
       end if
    end if
   
  end subroutine CNNDynamicsInit

  !-----------------------------------------------------------------------
  subroutine readCNNDynamicsParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read in parameters
    !
    ! !USES:
    use ncdio_pio   , only : file_desc_t,ncd_io
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNNDynamicsParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    
    call CNNDynamicsInit()

    tString='sf_minn'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNDynamicsParamsInst%sf=tempr

    tString='sf_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNDynamicsParamsInst%sf_no3=tempr
   
  end subroutine readCNNDynamicsParams

  !-----------------------------------------------------------------------
  subroutine CNNDeposition( bounds, &
       atm2lnd_vars, nitrogenflux_vars )
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen deposition rate
    ! from atmospheric forcing. For now it is assumed that all the atmospheric
    ! N deposition goes to the soil mineral N pool.
    ! This could be updated later to divide the inputs between mineral N absorbed
    ! directly into the canopy and mineral N entering the soil pool.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: g,c                    ! indices
    !-----------------------------------------------------------------------
    
    associate(& 
         forc_ndep     =>  atm2lnd_vars%forc_ndep_grc           , & ! Input:  [real(r8) (:)]  nitrogen deposition rate (gN/m2/s)                
         ndep_to_sminn =>  nitrogenflux_vars%ndep_to_sminn_col   & ! Output: [real(r8) (:)]                                                    
         )
      
      ! Loop through columns
      do c = bounds%begc, bounds%endc
         g = col%gridcell(c)
         ndep_to_sminn(c) = forc_ndep(g)
      end do

    end associate

  end subroutine CNNDeposition

  !-----------------------------------------------------------------------
  subroutine CNNFixation(num_soilc, filter_soilc, &
       carbonflux_vars, nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen fixation rate
    ! as a function of annual total NPP. This rate gets updated once per year.
    ! All N fixation goes to the soil mineral N pool.
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(carbonflux_type)   , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars 
    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc                  ! indices
    real(r8) :: t                     ! temporary
    real(r8) :: dayspyr               ! days per year
    !-----------------------------------------------------------------------

    associate(& 
         cannsum_npp    => carbonflux_vars%annsum_npp_col      , & ! Input:  [real(r8) (:)]  nitrogen deposition rate (gN/m2/s)                
         col_lag_npp    => carbonflux_vars%lag_npp_col         , & ! Input: [real(r8) (:)]  (gC/m2/s) lagged net primary production           

         nfix_to_sminn  => nitrogenflux_vars%nfix_to_sminn_col   & ! Output: [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
         )

      dayspyr = get_days_per_year()

      if ( nfix_timeconst > 0._r8 .and. nfix_timeconst < 500._r8 ) then
         ! use exponential relaxation with time constant nfix_timeconst for NPP - NFIX relation
         ! Loop through columns
         do fc = 1,num_soilc
            c = filter_soilc(fc)         

            if (col_lag_npp(c) /= spval) then
               ! need to put npp in units of gC/m^2/year here first
               t = (1.8_r8 * (1._r8 - exp(-0.003_r8 * col_lag_npp(c)*(secspday * dayspyr))))/(secspday * dayspyr)  
               nfix_to_sminn(c) = max(0._r8,t)
            else
               nfix_to_sminn(c) = 0._r8
            endif
         end do
      else
         ! use annual-mean values for NPP-NFIX relation
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            t = (1.8_r8 * (1._r8 - exp(-0.003_r8 * cannsum_npp(c))))/(secspday * dayspyr)
            nfix_to_sminn(c) = max(0._r8,t)
         end do
      endif

    end associate

  end subroutine CNNFixation
 
  !-----------------------------------------------------------------------
  subroutine CNNLeaching(bounds, num_soilc, filter_soilc, &
       waterstate_vars, waterflux_vars, nitrogenstate_vars, nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen leaching rate
    ! as a function of soluble mineral N and total soil water outflow.
    !
    ! !USES:
    use clm_varpar       , only : nlevdecomp, nlevsoi
    use clm_time_manager , only : get_step_size
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars 
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,fc                                 ! indices
    real(r8) :: dt                                     ! radiation time step (seconds)
    real(r8) :: sf                                     ! soluble fraction of mineral N (unitless)
    real(r8) :: sf_no3                                 ! soluble fraction of NO3 (unitless)
    real(r8) :: disn_conc                              ! dissolved mineral N concentration (gN/kg water)
    real(r8) :: tot_water(bounds%begc:bounds%endc)     ! total column liquid water (kg water/m2)
    real(r8) :: surface_water(bounds%begc:bounds%endc) ! liquid water to shallow surface depth (kg water/m2)
    real(r8) :: drain_tot(bounds%begc:bounds%endc)     ! total drainage flux (mm H2O /s)
    real(r8), parameter :: depth_runoff_Nloss = 0.05   ! (m) depth over which runoff mixes with soil water for N loss to runoff
    !-----------------------------------------------------------------------

    associate(& 
         h2osoi_liq          => waterstate_vars%h2osoi_liq_col            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)

         qflx_drain          => waterflux_vars%qflx_drain_col             , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf           => waterflux_vars%qflx_surf_col              , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)                        
         
         sminn_vr            => nitrogenstate_vars%sminn_vr_col           , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N                          
         smin_no3_vr         => nitrogenstate_vars%smin_no3_vr_col        , & ! Input:  [real(r8) (:,:) ]                                                  
         sminn_leached_vr    => nitrogenflux_vars%sminn_leached_vr_col    , & ! Output: [real(r8) (:,:) ]  rate of mineral N leaching (gN/m3/s)            
         smin_no3_leached_vr => nitrogenflux_vars%smin_no3_leached_vr_col , & ! Output: [real(r8) (:,:) ]  rate of mineral NO3 leaching (gN/m3/s)          
         smin_no3_runoff_vr  => nitrogenflux_vars%smin_no3_runoff_vr_col    & ! Output: [real(r8) (:,:) ]  rate of mineral NO3 loss with runoff (gN/m3/s)  
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      if (.not. use_nitrif_denitrif) then
         ! set constant sf 
         sf = CNNDynamicsParamsInst%sf
      else
         ! Assume that 100% of the soil NO3 is in a soluble form
         sf_no3 =  CNNDynamicsParamsInst%sf_no3 
      end if

      ! calculate the total soil water
      tot_water(bounds%begc:bounds%endc) = 0._r8
      do j = 1,nlevsoi
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            tot_water(c) = tot_water(c) + h2osoi_liq(c,j)
         end do
      end do

      ! for runoff calculation; calculate total water to a given depth
      surface_water(bounds%begc:bounds%endc) = 0._r8
      do j = 1,nlevsoi
         if ( zisoi(j) <= depth_runoff_Nloss)  then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               surface_water(c) = surface_water(c) + h2osoi_liq(c,j)
            end do
         elseif ( zisoi(j-1) < depth_runoff_Nloss)  then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               surface_water(c) = surface_water(c) + h2osoi_liq(c,j) * ( (depth_runoff_Nloss - zisoi(j-1)) / col%dz(c,j))
            end do
         endif
      end do

      ! Loop through columns
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         drain_tot(c) = qflx_drain(c)
      end do


      if (.not. use_nitrif_denitrif) then

         !----------------------------------------
         ! --------- NITRIF_NITRIF OFF------------
         !----------------------------------------

         do j = 1,nlevdecomp
            ! Loop through columns
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (.not. use_vertsoilc) then
                  ! calculate the dissolved mineral N concentration (gN/kg water)
                  ! assumes that 10% of mineral nitrogen is soluble
                  disn_conc = 0._r8
                  if (tot_water(c) > 0._r8) then
                     disn_conc = (sf * sminn_vr(c,j) ) / tot_water(c)
                  end if

                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  sminn_leached_vr(c,j) = disn_conc * drain_tot(c)
               else
                  ! calculate the dissolved mineral N concentration (gN/kg water)
                  ! assumes that 10% of mineral nitrogen is soluble
                  disn_conc = 0._r8
                  if (h2osoi_liq(c,j) > 0._r8) then
                     disn_conc = (sf * sminn_vr(c,j) * col%dz(c,j) )/(h2osoi_liq(c,j) )
                  end if

                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  sminn_leached_vr(c,j) = disn_conc * drain_tot(c) * h2osoi_liq(c,j) / ( tot_water(c) * col%dz(c,j) )

               end if

               ! limit the flux based on current sminn state
               ! only let at most the assumed soluble fraction
               ! of sminn be leached on any given timestep
               sminn_leached_vr(c,j) = min(sminn_leached_vr(c,j), (sf * sminn_vr(c,j))/dt)

               ! limit the flux to a positive value
               sminn_leached_vr(c,j) = max(sminn_leached_vr(c,j), 0._r8)

            end do
         end do

      else     

         !----------------------------------------
         ! --------- NITRIF_NITRIF ON-------------
         !----------------------------------------

         do j = 1,nlevdecomp
            ! Loop through columns
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (.not. use_vertsoilc) then
                  ! calculate the dissolved mineral N concentration (gN/kg water)
                  ! assumes that 10% of mineral nitrogen is soluble
                  disn_conc = 0._r8
                  if (tot_water(c) > 0._r8) then
                     disn_conc = (sf_no3 * smin_no3_vr(c,j) )/tot_water(c)
                  end if

                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  smin_no3_leached_vr(c,j) = disn_conc * drain_tot(c)
               else
                  ! calculate the dissolved mineral N concentration (gN/kg water)
                  ! assumes that 10% of mineral nitrogen is soluble
                  disn_conc = 0._r8
                  if (h2osoi_liq(c,j) > 0._r8) then
                     disn_conc = (sf_no3 * smin_no3_vr(c,j) * col%dz(c,j) )/(h2osoi_liq(c,j) )
                  end if
                  !
                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  smin_no3_leached_vr(c,j) = disn_conc * drain_tot(c) * h2osoi_liq(c,j) / ( tot_water(c) * col%dz(c,j) )
                  !
                  ! ensure that leaching rate isn't larger than soil N pool
                  smin_no3_leached_vr(c,j) = min(smin_no3_leached_vr(c,j), smin_no3_vr(c,j) / dt )
                  !
                  ! limit the leaching flux to a positive value
                  smin_no3_leached_vr(c,j) = max(smin_no3_leached_vr(c,j), 0._r8)
                  !
                  !
                  ! calculate the N loss from surface runoff, assuming a shallow mixing of surface waters into soil and removal based on runoff
                  if ( zisoi(j) <= depth_runoff_Nloss )  then
                     smin_no3_runoff_vr(c,j) = disn_conc * qflx_surf(c) * &
                          h2osoi_liq(c,j) / ( surface_water(c) * col%dz(c,j) )
                  elseif ( zisoi(j-1) < depth_runoff_Nloss )  then
                     smin_no3_runoff_vr(c,j) = disn_conc * qflx_surf(c) * &
                          h2osoi_liq(c,j) * ((depth_runoff_Nloss - zisoi(j-1)) / &
                          col%dz(c,j)) / ( surface_water(c) * (depth_runoff_Nloss-zisoi(j-1) ))
                  else
                     smin_no3_runoff_vr(c,j) = 0._r8
                  endif
                  !
                  ! ensure that runoff rate isn't larger than soil N pool
                  smin_no3_runoff_vr(c,j) = min(smin_no3_runoff_vr(c,j), smin_no3_vr(c,j) / dt - smin_no3_leached_vr(c,j))
                  !
                  ! limit the flux to a positive value
                  smin_no3_runoff_vr(c,j) = max(smin_no3_runoff_vr(c,j), 0._r8)


               endif
               ! limit the flux based on current smin_no3 state
               ! only let at most the assumed soluble fraction
               ! of smin_no3 be leached on any given timestep
               smin_no3_leached_vr(c,j) = min(smin_no3_leached_vr(c,j), (sf_no3 * smin_no3_vr(c,j))/dt)

               ! limit the flux to a positive value
               smin_no3_leached_vr(c,j) = max(smin_no3_leached_vr(c,j), 0._r8)

            end do
         end do
      endif

    end associate
  end subroutine CNNLeaching

  !-----------------------------------------------------------------------
  subroutine CNNFert(bounds, num_soilc, filter_soilc, &
       nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen fertilizer for crops
    ! All fertilizer goes into the soil mineral N pool.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds  
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars 
    !
    ! !LOCAL VARIABLES:
    integer :: c,fc                 ! indices
    !-----------------------------------------------------------------------

    associate(&   
         fert          =>    nitrogenflux_vars%fert_patch          , & ! Input:  [real(r8) (:)]  nitrogen fertilizer rate (gN/m2/s)                
         fert_to_sminn =>    nitrogenflux_vars%fert_to_sminn_col   & ! Output: [real(r8) (:)]                                                    
         )
      
      call p2c(bounds, num_soilc, filter_soilc, &
           fert(bounds%begp:bounds%endp), &
           fert_to_sminn(bounds%begc:bounds%endc))

    end associate
  end subroutine CNNFert

  !-----------------------------------------------------------------------
  subroutine CNSoyfix (bounds, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       waterstate_vars, crop_vars, cnstate_vars, nitrogenstate_vars, nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    ! This routine handles the fixation of nitrogen for soybeans based on
    ! the EPICPHASE model M. Cabelguenne et al., Agricultural systems 60: 175-196, 1999
    ! N-fixation is based on soil moisture, plant growth phase, and availibility of
    ! nitrogen in the soil root zone.
    !
    ! !USES:
    use pftvarcon  , only : nsoybean
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars 
    !
    ! !LOCAL VARIABLES:
    integer :: fp,p,c
    real(r8):: fxw,fxn,fxg,fxr             ! soil water factor, nitrogen factor, growth stage factor
    real(r8):: soy_ndemand                 ! difference between nitrogen supply and demand
    real(r8):: GDDfrac
    real(r8):: sminnthreshold1, sminnthreshold2
    real(r8):: GDDfracthreshold1, GDDfracthreshold2
    real(r8):: GDDfracthreshold3, GDDfracthreshold4
    !-----------------------------------------------------------------------

    associate(                                                         & 
         wf               =>  waterstate_vars%wf_col                 , & ! Input:  [real(r8) (:) ]  soil water as frac. of whc for top 0.5 m          

         hui              =>  crop_vars%gddplant_patch               , & ! Input:  [real(r8) (:) ]  gdd since planting (gddplant)                    

         fpg              =>  cnstate_vars%fpg_col                   , & ! Input:  [real(r8) (:) ]  fraction of potential gpp (no units)              
         gddmaturity      =>  cnstate_vars%gddmaturity_patch         , & ! Input:  [real(r8) (:) ]  gdd needed to harvest                             
         croplive         =>  cnstate_vars%croplive_patch            , & ! Input:  [logical  (:) ]  true if planted and not harvested                  

         sminn            =>  nitrogenstate_vars%sminn_col           , & ! Input:  [real(r8) (:) ]  (kgN/m2) soil mineral N                           
         plant_ndemand    =>  nitrogenflux_vars%plant_ndemand_patch  , & ! Input:  [real(r8) (:) ]  N flux required to support initial GPP (gN/m2/s)  
         
         soyfixn          =>  nitrogenflux_vars%soyfixn_patch        , & ! Output: [real(r8) (:) ]  nitrogen fixed to each soybean crop               
         soyfixn_to_sminn =>  nitrogenflux_vars%soyfixn_to_sminn_col   & ! Output: [real(r8) (:) ]                                                    
         )

      sminnthreshold1 = 30._r8
      sminnthreshold2 = 10._r8
      GDDfracthreshold1 = 0.15_r8
      GDDfracthreshold2 = 0.30_r8
      GDDfracthreshold3 = 0.55_r8
      GDDfracthreshold4 = 0.75_r8

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = pft%column(p)

         ! if soybean currently growing then calculate fixation

         if (pft%itype(p) == nsoybean .and. croplive(p)) then

            ! difference between supply and demand

            if (fpg(c) < 1._r8) then
               soy_ndemand = 0._r8
               soy_ndemand = plant_ndemand(p) - plant_ndemand(p)*fpg(c)

               ! fixation depends on nitrogen, soil water, and growth stage

               ! soil water factor

               fxw = 0._r8
               fxw = wf(c)/0.85_r8

               ! soil nitrogen factor (Beth says: CHECK UNITS)

               if (sminn(c) > sminnthreshold1) then
                  fxn = 0._r8
               else if (sminn(c) > sminnthreshold2 .and. sminn(c) <= sminnthreshold1) then
                  fxn = 1.5_r8 - .005_r8 * (sminn(c) * 10._r8)
               else if (sminn(c) <= sminnthreshold2) then
                  fxn = 1._r8
               end if

               ! growth stage factor
               ! slevis: to replace GDDfrac, assume...
               ! Beth's crit_offset_gdd_def is similar to my gddmaturity
               ! Beth's ac_gdd (base 5C) similar to my hui=gddplant (base 10
               ! for soy) 
               ! Ranges below are not firm. Are they lit. based or tuning based?

               GDDfrac = hui(p) / gddmaturity(p)

               if (GDDfrac <= GDDfracthreshold1) then
                  fxg = 0._r8
               else if (GDDfrac > GDDfracthreshold1 .and. GDDfrac <= GDDfracthreshold2) then
                  fxg = 6.67_r8 * GDDfrac - 1._r8
               else if (GDDfrac > GDDfracthreshold2 .and. GDDfrac <= GDDfracthreshold3) then
                  fxg = 1._r8
               else if (GDDfrac > GDDfracthreshold3 .and. GDDfrac <= GDDfracthreshold4) then
                  fxg = 3.75_r8 - 5._r8 * GDDfrac
               else  ! GDDfrac > GDDfracthreshold4
                  fxg = 0._r8
               end if

               ! calculate the nitrogen fixed by the soybean

               fxr = min(1._r8, fxw, fxn) * fxg 
               fxr = max(0._r8, fxr)
               soyfixn(p) =  fxr * soy_ndemand
               soyfixn(p) = min(soyfixn(p), soy_ndemand)

            else ! if nitrogen demand met, no fixation

               soyfixn(p) = 0._r8

            end if

         else ! if not live soybean, no fixation

            soyfixn(p) = 0._r8

         end if
      end do

      call p2c(bounds, num_soilc, filter_soilc, &
           soyfixn(bounds%begp:bounds%endp), &
           soyfixn_to_sminn(bounds%begc:bounds%endc))

    end associate

  end subroutine CNSoyfix

end module CNNDynamicsMod
