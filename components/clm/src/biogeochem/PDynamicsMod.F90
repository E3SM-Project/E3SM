module PDynamicsMod

  !-----------------------------------------------------------------------
  ! !MODULE: PDynamicsMod
  !
  ! !DESCRIPTION:
  ! Module for inorganic phosphorus dynamics
  ! for coupled carbon-nitrogen-phosphorus code.
  ! X.YANG
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use decompMod           , only : bounds_type
  use clm_varcon          , only : dzsoi_decomp, zisoi
  use subgridAveMod       , only : p2c
  use atm2lndType         , only : atm2lnd_type
  use CNCarbonFluxType    , only : carbonflux_type, nfix_timeconst
  
  use clm_varpar          , only : nlevdecomp
  use clm_varctl          , only : use_vertsoilc
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use CNNitrogenStateType , only : nitrogenstate_type 

  use CNStateType         , only : cnstate_type
  use WaterStateType      , only : waterstate_type
  use WaterFluxType       , only : waterflux_type
  use CropType            , only : crop_type
  use ColumnType          , only : col
  use PatchType           , only : pft
  use EcophysConType      , only : ecophyscon
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: PDeposition
  public :: PWeathering
  public :: PAdsorption
  public :: PDesorption
  public :: POcclusion
  public :: PBiochemMin
  public :: PLeaching
  public :: PBiochemMin_balance

  !-----------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------
  subroutine PDeposition( bounds, &
       atm2lnd_vars, phosphorusflux_vars )
    ! BY X. SHI
    ! !DESCRIPTION:
    ! On the radiation time step, update the phosphorus deposition rate
    ! from atmospheric forcing. For now it is assumed that all the atmospheric
    ! P deposition goes to the soil mineral P pool.
    ! This could be updated later to divide the inputs between mineral P absorbed
    ! directly into the canopy and mineral P entering the soil pool.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(phosphorusflux_type) , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: g,c                    ! indices
    !-----------------------------------------------------------------------

    associate(&
         forc_pdep     =>  atm2lnd_vars%forc_pdep_grc           , & ! Input:  [real(r8) (:)]  Phosphorus deposition rate (gP/m2/s)                
         pdep_to_sminp =>  phosphorusflux_vars%pdep_to_sminp_col   & ! Output: [real(r8) (:)]                                                    
         )

      ! Loop through columns
      do c = bounds%begc, bounds%endc
         g = col%gridcell(c)
         pdep_to_sminp(c) = forc_pdep(g)
      end do

    end associate

  end subroutine PDeposition

  !-----------------------------------------------------------------------
  subroutine PWeathering(num_soilc, filter_soilc, &
       cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
    !
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
    use soilorder_varcon, only: r_weather
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(phosphorusstate_type), intent(in) ::  phosphorusstate_vars
    type(phosphorusflux_type) , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc                  ! indices
    real(r8) :: t                     ! temporary
    real(r8) :: dayspyr               ! days per year

!   !OTHER LOCAL VARIABLES
    real(r8)     :: r_weather_c
    real(r8)     :: rr
    real(r8):: dt           !decomp timestep (seconds)
    real(r8):: dtd          !decomp timestep (days)
    integer :: j

    !-----------------------------------------------------------------------

    associate(&

         isoilorder     => cnstate_vars%isoilorder                 ,&
         primp          => phosphorusstate_vars%primp_vr_col       ,& 
         primp_to_labilep => phosphorusflux_vars%primp_to_labilep_vr_col  &         

         )

      dayspyr = get_days_per_year()

      ! set time steps
      dt = real( get_step_size(), r8 )
      dtd = dt/(30._r8*secspday)
   
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
      
            !! read in monthly rate is converted to that in half hour
            r_weather_c = r_weather( isoilorder(c) )
            rr=-log(1._r8-r_weather_c)
            r_weather_c=1._r8-exp(-rr*dtd)
      
            primp_to_labilep(c,j) = primp(c,j)*r_weather_c/dt
      !     primp_to_labilep(c,j) = 0.005_r8/(365._r8*24._r8*3600._r8)
      !     primp_to_labilep(c,j) = 0._r8       
         end do
      enddo
    end associate

  end subroutine PWeathering
  !-----------------------------------------------------------------------



  !-----------------------------------------------------------------------
  subroutine PAdsorption(num_soilc, filter_soilc, &
       cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
    !
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
    use soilorder_varcon , only : r_adsorp
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(phosphorusstate_type), intent(in) ::  phosphorusstate_vars
    type(phosphorusflux_type) , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc                  ! indices
    real(r8) :: t                     ! temporary
    real(r8) :: dayspyr               ! days per year

!   !OTHER LOCAL VARIABLES
    real(r8)     :: r_adsorp_c
    real(r8)     :: rr
    real(r8):: dt           !decomp timestep (seconds)
    real(r8):: dtd          !decomp timestep (days)
    integer :: j

    !-----------------------------------------------------------------------

    associate(&

         isoilorder     => cnstate_vars%isoilorder                 ,&
         solutionp   => phosphorusstate_vars%solutionp_vr_col      ,&
         labilep     => phosphorusstate_vars%labilep_vr_col        ,&
         labilep_to_secondp => phosphorusflux_vars%labilep_to_secondp_vr_col &

         )

      dayspyr = get_days_per_year()

      ! set time steps
      dt = real( get_step_size(), r8 )
      dtd = dt/(30._r8*secspday)
   
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
   
            ! calculate rate at half-hour time step
            r_adsorp_c = r_adsorp( isoilorder(c) )
            rr=-log(1._r8-r_adsorp_c)
            r_adsorp_c = 1._r8-exp(-rr*dtd)

            if(labilep(c,j) > 0._r8)then
               labilep_to_secondp(c,j) = ( labilep(c,j) )*r_adsorp_c/dt
            else
               labilep_to_secondp(c,j) = 0._r8
            end if

         end do
       end do
    end associate

  end subroutine PAdsorption


  !-----------------------------------------------------------------------
  subroutine PDesorption(num_soilc, filter_soilc, &
       cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
    !
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
    use soilorder_varcon , only : r_desorp
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(phosphorusstate_type), intent(in) ::  phosphorusstate_vars
    type(phosphorusflux_type) , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc                  ! indices
    real(r8) :: t                     ! temporary
    real(r8) :: dayspyr               ! days per year

!   !OTHER LOCAL VARIABLES
    real(r8)     :: r_desorp_c
    real(r8)     :: rr
    real(r8):: dt           !decomp timestep (seconds)
    real(r8):: dtd          !decomp timestep (days)
    integer :: j

    !-----------------------------------------------------------------------

    associate(&

         isoilorder     => cnstate_vars%isoilorder              ,&
         secondp     => phosphorusstate_vars%secondp_vr_col     ,&
         secondp_to_labilep => phosphorusflux_vars%secondp_to_labilep_vr_col &

         )

      dayspyr = get_days_per_year()

      ! set time steps
      dt = real( get_step_size(), r8 )
      dtd = dt/(30._r8*secspday)
   
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
   
            ! calculate rate at half-hour time step
            r_desorp_c = r_desorp( isoilorder(c) )
            rr=-log(1._r8-r_desorp_c)
            r_desorp_c = 1._r8-exp(-rr*dtd)
    
            if(secondp(c,j) > 0._r8)then
              secondp_to_labilep(c,j) = secondp(c,j)*r_desorp_c/dt
            else
              secondp_to_labilep(c,j) = 0._r8
            endif

         end do
       end do
    end associate

  end subroutine PDesorption



  !-----------------------------------------------------------------------
  subroutine POcclusion(num_soilc, filter_soilc, &
       cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
    !
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
    use soilorder_varcon , only : r_occlude
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(phosphorusstate_type), intent(in) ::  phosphorusstate_vars
    type(phosphorusflux_type) , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc                  ! indices
    real(r8) :: t                     ! temporary
    real(r8) :: dayspyr               ! days per year

!   !OTHER LOCAL VARIABLES
    real(r8)     :: r_occlude_c
    real(r8)     :: rr
    real(r8):: dt           !decomp timestep (seconds)
    real(r8):: dtd          !decomp timestep (days)
    integer :: j
 
    !-----------------------------------------------------------------------

    associate(&

         isoilorder     => cnstate_vars%isoilorder                      ,&
         secondp     => phosphorusstate_vars%secondp_vr_col             ,&
         secondp_to_occlp => phosphorusflux_vars%secondp_to_occlp_vr_col &

         )

      dayspyr = get_days_per_year()

      ! set time steps
      dt = real( get_step_size(), r8 )
      dtd = dt/(30._r8*secspday)
   
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
   

            ! calculate rate at half-hour time step
            r_occlude_c = r_occlude( isoilorder(c) )
            rr=-log(1._r8-r_occlude_c)
            r_occlude_c = 1._r8-exp(-rr*dtd)
    
            if(secondp(c,j) > 0._r8)then
               secondp_to_occlp(c,j) = secondp(c,j)*r_occlude_c/dt
            else
               secondp_to_occlp(c,j) =0._r8
            endif

         end do
       end do
    end associate

  end subroutine POcclusion

  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  subroutine PLeaching(bounds, num_soilc, filter_soilc, &
       waterstate_vars, waterflux_vars, phosphorusstate_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the phosphorus leaching rate
    ! as a function of solution P and total soil water outflow.
    !
    ! !USES:
    use clm_varpar       , only : nlevsoi
    use clm_time_manager , only : get_step_size
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,fc                                 ! indices
    real(r8) :: dt                                     ! radiation time step(seconds)
    real(r8) :: disp_conc                              ! dissolved mineral N concentration (gP/kg water)
    real(r8) :: tot_water(bounds%begc:bounds%endc)     ! total column liquid water (kg water/m2)
    real(r8) :: surface_water(bounds%begc:bounds%endc) ! liquid water to shallow surface depth (kg water/m2)
    real(r8) :: drain_tot(bounds%begc:bounds%endc)     ! total drainage flux (mmH2O /s)
    real(r8), parameter :: depth_runoff_Ploss = 0.05   ! (m) depth over which runoff mixes with soil water for P loss to runoff
    !-----------------------------------------------------------------------

    associate(&
         h2osoi_liq          => waterstate_vars%h2osoi_liq_col            , & !Input:  [real(r8) (:,:) ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)

         qflx_drain          => waterflux_vars%qflx_drain_col             , & !Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf           => waterflux_vars%qflx_surf_col              , & !Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)                        

         solutionp_vr            => phosphorusstate_vars%solutionp_vr_col           , & !Input:  [real(r8) (:,:) ]  (gP/m3) soil mineral N                          
         sminp_leached_vr    => phosphorusflux_vars%sminp_leached_vr_col     & !Output: [real(r8) (:,:) ]  rate of mineral N leaching (gP/m3/s)            
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

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
         if ( zisoi(j) <= depth_runoff_Ploss)  then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               surface_water(c) = surface_water(c) + h2osoi_liq(c,j)
            end do
         elseif ( zisoi(j-1) < depth_runoff_Ploss)  then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               surface_water(c) = surface_water(c) + h2osoi_liq(c,j) * ((depth_runoff_Ploss - zisoi(j-1)) / col%dz(c,j))
            end do
         endif
      end do

      ! Loop through columns
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         drain_tot(c) = qflx_drain(c)
      end do

         do j = 1,nlevdecomp
            ! Loop through columns
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (.not. use_vertsoilc) then
                  disp_conc = 0._r8
                  if (tot_water(c) > 0._r8) then
                     disp_conc = ( solutionp_vr(c,j) ) / tot_water(c)
                  end if

                  ! calculate the P leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  sminp_leached_vr(c,j) = disp_conc * drain_tot(c)
               else
                  disp_conc = 0._r8
                  if (h2osoi_liq(c,j) > 0._r8) then
                     disp_conc = (solutionp_vr(c,j) * col%dz(c,j))/(h2osoi_liq(c,j) )
                  end if

                  ! calculate the P leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  sminp_leached_vr(c,j) = disp_conc * drain_tot(c) *h2osoi_liq(c,j) / ( tot_water(c) * col%dz(c,j) )

               end if
               ! limit the flux based on current sminp state
               ! only let at most the assumed soluble fraction
               ! of sminp be leached on any given timestep
               sminp_leached_vr(c,j) = min(sminp_leached_vr(c,j), (solutionp_vr(c,j))/dt)

               ! limit the flux to a positive value
               sminp_leached_vr(c,j) = max(sminp_leached_vr(c,j), 0._r8)
            end do
         end do

    end associate
  end subroutine PLeaching


  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------


  subroutine PBiochemMin(bounds,num_soilc, filter_soilc, &
       cnstate_vars, phosphorusstate_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the phosphorus leaching rate
    ! as a function of solution P and total soil water outflow.
    !
    ! !USES:
    use clm_varpar       , only : nlevsoi
    use clm_varpar       , only : ndecomp_pools
    use clm_time_manager , only : get_step_size
    use soilorder_varcon , only:k_s1_biochem,k_s2_biochem,k_s3_biochem,k_s4_biochem
    use clm_varcon       , only : secspday, spval

    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnstate_type)         , intent(in)    :: cnstate_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    integer  :: c,fc,j,l
    real(r8) :: rr
    real(r8):: dt           !decomp timestep (seconds)
    real(r8):: dtd          !decomp timestep (days)
    real(r8):: k_s1_biochem_c         !specfic biochemical rate constant SOM 1
    real(r8):: k_s2_biochem_c         !specfic biochemical rate constant SOM 1
    real(r8):: k_s3_biochem_c         !specfic biochemical rate constant SOM 1
    real(r8):: k_s4_biochem_c         !specfic biochemical rate constant SOM 1
    real(r8):: r_bc

    !-----------------------------------------------------------------------
    !!!!  decomp_ppools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)

    associate(&

         isoilorder     => cnstate_vars%isoilorder                            ,&
         decomp_ppools_vr_col => phosphorusstate_vars%decomp_ppools_vr_col    ,&
  
         biochem_pmin_ppools_vr_col  => phosphorusflux_vars%biochem_pmin_ppools_vr_col  ,&
         biochem_pmin_vr_col  => phosphorusflux_vars%biochem_pmin_vr_col      ,&
         biochem_pmin_col     => phosphorusflux_vars%biochem_pmin_col         , & 
         fpi_vr_col           => cnstate_vars%fpi_vr_col                      ,&
         fpi_p_vr_col           => cnstate_vars%fpi_p_vr_col                   &
         )

      dt = real( get_step_size(), r8 )
      dtd = dt/(30._r8*secspday)
      r_bc = -5._r8

      ! set initial values for potential C and N fluxes
      biochem_pmin_ppools_vr_col(bounds%begc : bounds%endc, :, :) = 0._r8

      do l = 1, ndecomp_pools
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               k_s1_biochem_c = k_s1_biochem( isoilorder(c) )
               k_s2_biochem_c = k_s2_biochem( isoilorder(c) )
               k_s3_biochem_c = k_s3_biochem( isoilorder(c) )
               k_s4_biochem_c = k_s4_biochem( isoilorder(c) )
         
               rr=-log(1._r8-k_s1_biochem_c)
               k_s1_biochem_c = 1-exp(-rr*dtd)
         
               rr=-log(1-k_s2_biochem_c)
               k_s2_biochem_c = 1-exp(-rr*dtd)
         
               rr=-log(1-k_s3_biochem_c)
               k_s3_biochem_c = 1-exp(-rr*dtd)
         
               rr=-log(1-k_s4_biochem_c)
               k_s4_biochem_c = 1-exp(-rr*dtd)

               if ( decomp_ppools_vr_col(c,j,l) > 0._r8 ) then

                 biochem_pmin_ppools_vr_col(c,j,l) = decomp_ppools_vr_col(c,j,l)* &
                                     k_s1_biochem_c * fpi_vr_col(c,j)*&
                                     (1._r8-exp(r_bc*(1-fpi_p_vr_col(c,j)) ) )/dt

               endif 

            end do
         end do
      end do

      
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            biochem_pmin_vr_col(c,j)=0._r8
            do l = 1, ndecomp_pools
               biochem_pmin_vr_col(c,j) = biochem_pmin_vr_col(c,j)+ &
                                          biochem_pmin_ppools_vr_col(c,j,l)
            enddo
         enddo
      enddo 


      
    end associate

  end subroutine PBiochemMin

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  
  subroutine PBiochemMin_balance(bounds,num_soilc, filter_soilc, &
       cnstate_vars,nitrogenstate_vars, phosphorusstate_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! created, Aug 2015 by Q. Zhu
    ! update the phosphatase activity induced P release based on Wang 2007
    !
    ! !USES:
    use pftvarcon              , only : noveg
    use clm_varpar             , only : ndecomp_pools
    use clm_time_manager       , only : get_step_size
    
    !
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)    :: bounds
    integer                    , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                    , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnstate_type)         , intent(in)    :: cnstate_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    integer  :: c,fc,p,j,l
    real(r8) :: lamda_up       ! nitrogen cost of phosphorus uptake
    real(r8) :: sop_profile(1:ndecomp_pools)
    real(r8) :: sop_tot
    integer  :: dt
    !-----------------------------------------------------------------------

    associate(                                                              &
         froot_prof           => cnstate_vars%froot_prof_patch            , & ! fine root vertical profile Zeng, X. 2001. Global vegetation root distribution for land modeling. J. Hydrometeor. 2:525-530
         biochem_pmin_vr      => phosphorusflux_vars%biochem_pmin_vr_col  , &
         biochem_pmin_ppools_vr_col  => phosphorusflux_vars%biochem_pmin_ppools_vr_col ,&
         npimbalance          => nitrogenstate_vars%npimbalance_patch     , &
         vmax_ptase_vr        => ecophyscon%vmax_ptase_vr                 , &
         km_ptase             => ecophyscon%km_ptase                      , &
         decomp_ppools_vr_col => phosphorusstate_vars%decomp_ppools_vr_col, &
         lamda_ptase          => ecophyscon%lamda_ptase                   ,  & ! critical value of nitrogen cost of phosphatase activity induced phosphorus uptake
         cn_scalar             => cnstate_vars%cn_scalar               , &
         cp_scalar             => cnstate_vars%cp_scalar                 &
         )

    dt = real( get_step_size(), r8 )

    ! set initial values for potential C and N fluxes
    biochem_pmin_ppools_vr_col(bounds%begc : bounds%endc, :, :) = 0._r8
      
    do j = 1,nlevdecomp
        do fc = 1,num_soilc
            c = filter_soilc(fc)
            biochem_pmin_vr(c,j) = 0.0_r8
            do p = col%pfti(c), col%pftf(c)
                if (pft%active(p).and. (pft%itype(p) .ne. noveg)) then
                    !lamda_up = npimbalance(p) ! partial_vcmax/partial_lpc / partial_vcmax/partial_lnc
                    lamda_up = cp_scalar(p)/max(cn_scalar(p),1e-20_r8)
                    lamda_up = min(max(lamda_up,0.0_r8), 150.0_r8)
                    biochem_pmin_vr(c,j) = biochem_pmin_vr(c,j) + &
                        vmax_ptase_vr(j) * max(lamda_up - lamda_ptase, 0.0_r8) / &
                        (km_ptase + max(lamda_up - lamda_ptase, 0.0_r8)) * froot_prof(p,j) * pft%wtcol(p)
                end if
            enddo
        enddo
    enddo 
    
    do j = 1,nlevdecomp
        do fc = 1,num_soilc
            c = filter_soilc(fc)
            sop_tot = 0._r8
            do l = 1,ndecomp_pools
                sop_tot = sop_tot + decomp_ppools_vr_col(c,j,l)
            end do
            do l = 1,ndecomp_pools
                if (sop_tot > 1e-12) then 
                    sop_profile(l) = decomp_ppools_vr_col(c,j,l)/sop_tot
                else
                    sop_profile(l) = 0._r8
                end if
            end do
            do l = 1,ndecomp_pools
                biochem_pmin_ppools_vr_col(c,j,l) = max(min(biochem_pmin_vr(c,j) * sop_profile(l), decomp_ppools_vr_col(c,j,l)),0._r8)
            end do
        end do
    end do

    do j = 1,nlevdecomp
        do fc = 1,num_soilc
            c = filter_soilc(fc)
            biochem_pmin_vr(c,j)=0._r8
            do l = 1, ndecomp_pools
               biochem_pmin_vr(c,j) = biochem_pmin_vr(c,j)+ &
                                          biochem_pmin_ppools_vr_col(c,j,l)*dt
            enddo
        enddo
    end do
    
    end associate

  end subroutine PBiochemMin_balance

end module PDynamicsMod
               
                
     

      








