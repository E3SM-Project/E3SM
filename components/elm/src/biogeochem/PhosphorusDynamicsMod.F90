module PhosphorusDynamicsMod

  !-----------------------------------------------------------------------
  ! !MODULE: PhosphorusDynamicsMod
  !
  ! !DESCRIPTION:
  ! Module for inorganic phosphorus dynamics
  ! for coupled carbon-nitrogen-phosphorus code.
  ! X.YANG
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use decompMod           , only : bounds_type
  use elm_varcon          , only : dzsoi_decomp, zisoi
  use subgridAveMod       , only : p2c
  use atm2lndType         , only : atm2lnd_type
  use CNCarbonFluxType    , only : carbonflux_type  
  use elm_varpar          , only : nlevdecomp
  use clm_varctl          , only : use_vertsoilc
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use CNNitrogenStateType , only : nitrogenstate_type 

  use CNStateType         , only : cnstate_type
  use WaterStateType      , only : waterstate_type
  use WaterFluxType       , only : waterflux_type
  use CropType            , only : crop_type
  use ColumnType          , only : col_pp
  use ColumnDataType      , only : col_ws, col_wf, nfix_timeconst, col_ps, col_pf
  use VegetationType      , only : veg_pp
  use VegetationDataType  , only : veg_ns, veg_pf
  use VegetationPropertiesType      , only : veg_vp
  use clm_varctl          , only : NFIX_PTASE_plant

  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: PhosphorusDeposition
  public :: PhosphorusWeathering
  public :: PhosphorusAdsportion
  public :: PhosphorusDesoprtion
  public :: PhosphorusOcclusion
  public :: PhosphorusBiochemMin
  public :: PhosphorusLeaching
  public :: PhosphorusBiochemMin_balance

  !-----------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------
  subroutine PhosphorusDeposition( bounds, &
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
         pdep_to_sminp =>  col_pf%pdep_to_sminp   & ! Output: [real(r8) (:)]                                                    
         )

      ! Loop through columns
      do c = bounds%begc, bounds%endc
         g = col_pp%gridcell(c)
         pdep_to_sminp(c) = forc_pdep(g)
      end do

    end associate

  end subroutine PhosphorusDeposition

  !-----------------------------------------------------------------------
  subroutine PhosphorusWeathering(num_soilc, filter_soilc, &
       cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
    !
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use elm_varcon       , only : secspday, spval
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
         primp          => col_ps%primp_vr       ,& 
         primp_to_labilep => col_pf%primp_to_labilep_vr  &         

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

  end subroutine PhosphorusWeathering
  !-----------------------------------------------------------------------



  !-----------------------------------------------------------------------
  subroutine PhosphorusAdsportion(num_soilc, filter_soilc, &
       cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
    !
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use elm_varcon       , only : secspday, spval
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
         solutionp   => col_ps%solutionp_vr      ,&
         labilep     => col_ps%labilep_vr        ,&
         labilep_to_secondp => col_pf%labilep_to_secondp_vr &

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

  end subroutine PhosphorusAdsportion


  !-----------------------------------------------------------------------
  subroutine PhosphorusDesoprtion(num_soilc, filter_soilc, &
       cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
    !
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use elm_varcon       , only : secspday, spval
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
         secondp     => col_ps%secondp_vr     ,&
         secondp_to_labilep => col_pf%secondp_to_labilep_vr &

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

  end subroutine PhosphorusDesoprtion



  !-----------------------------------------------------------------------
  subroutine PhosphorusOcclusion(num_soilc, filter_soilc, &
       cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
    !
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use elm_varcon       , only : secspday, spval
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
         secondp     => col_ps%secondp_vr             ,&
         secondp_to_occlp => col_pf%secondp_to_occlp_vr &

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

  end subroutine PhosphorusOcclusion

  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  subroutine PhosphorusLeaching(bounds, num_soilc, filter_soilc, &
       waterstate_vars, waterflux_vars, phosphorusstate_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the phosphorus leaching rate
    ! as a function of solution P and total soil water outflow.
    !
    ! !USES:
    use elm_varpar       , only : nlevsoi
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
         h2osoi_liq          => col_ws%h2osoi_liq            , & !Input:  [real(r8) (:,:) ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)

         qflx_drain          => col_wf%qflx_drain             , & !Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf           => col_wf%qflx_surf              , & !Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)                        

         solutionp_vr            => col_ps%solutionp_vr           , & !Input:  [real(r8) (:,:) ]  (gP/m3) soil mineral N                          
         sminp_leached_vr    => col_pf%sminp_leached_vr     & !Output: [real(r8) (:,:) ]  rate of mineral N leaching (gP/m3/s)            
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
               surface_water(c) = surface_water(c) + h2osoi_liq(c,j) * ((depth_runoff_Ploss - zisoi(j-1)) / col_pp%dz(c,j))
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
                     disp_conc = (solutionp_vr(c,j) * col_pp%dz(c,j))/(h2osoi_liq(c,j) )
                  end if

                  ! calculate the P leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  sminp_leached_vr(c,j) = disp_conc * drain_tot(c) *h2osoi_liq(c,j) / ( tot_water(c) * col_pp%dz(c,j) )

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
  end subroutine PhosphorusLeaching


  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------


  subroutine PhosphorusBiochemMin(bounds,num_soilc, filter_soilc, &
       cnstate_vars, phosphorusstate_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the phosphorus leaching rate
    ! as a function of solution P and total soil water outflow.
    !
    ! !USES:
    use elm_varpar       , only : nlevsoi
    use elm_varpar       , only : ndecomp_pools
    use clm_time_manager , only : get_step_size
    use soilorder_varcon , only:k_s1_biochem,k_s2_biochem,k_s3_biochem,k_s4_biochem
    use elm_varcon       , only : secspday, spval

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
         decomp_ppools_vr_col => col_ps%decomp_ppools_vr    ,&
  
         biochem_pmin_ppools_vr_col  => col_pf%biochem_pmin_ppools_vr  ,&
         biochem_pmin_vr_col  => col_pf%biochem_pmin_vr      ,&
         biochem_pmin_col     => col_pf%biochem_pmin         , & 
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
                                     (1._r8-exp(r_bc*(1._r8-fpi_p_vr_col(c,j)) ))/dt


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

  end subroutine PhosphorusBiochemMin

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  
  subroutine PhosphorusBiochemMin_balance(bounds,num_soilc, filter_soilc, &
       cnstate_vars,nitrogenstate_vars, phosphorusstate_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! created, Aug 2015 by Q. Zhu
    ! update the phosphatase activity induced P release based on Wang 2007
    !
    ! !USES:
    use pftvarcon              , only : noveg
    use elm_varpar             , only : ndecomp_pools
    use clm_time_manager       , only : get_step_size
    use CNDecompCascadeConType , only : decomp_cascade_con
 
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
    real(r8) :: biochem_pmin_to_ecosysp_vr_col_pot(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: biochem_pmin_to_plant_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp)
    real(r8) :: sop_tot
    integer  :: dt
    real(r8) :: ptase_tmp
    !-----------------------------------------------------------------------

    associate(                                                              &
         froot_prof           => cnstate_vars%froot_prof_patch            , & ! fine root vertical profile Zeng, X. 2001. Global vegetation root distribution for land modeling. J. Hydrometeor. 2:525-530
         biochem_pmin_vr      => col_pf%biochem_pmin_vr  , &
         biochem_pmin_ppools_vr_col  => col_pf%biochem_pmin_ppools_vr ,&
         biochem_pmin_to_ecosysp_vr_col => col_pf%biochem_pmin_to_ecosysp_vr , &
         biochem_pmin_to_plant_patch    => veg_pf%biochem_pmin_to_plant , &
         npimbalance          => veg_ns%npimbalance     , &
         vmax_ptase           => veg_vp%vmax_ptase                    , &
         km_ptase             => veg_vp%km_ptase                      , &
         alpha_ptase          => veg_vp%alpha_ptase                   , &
         decomp_ppools_vr_col => col_ps%decomp_ppools_vr, &
         lamda_ptase          => veg_vp%lamda_ptase                   ,  & ! critical value of nitrogen cost of phosphatase activity induced phosphorus uptake
         cn_scalar             => cnstate_vars%cn_scalar               , &
         cp_scalar             => cnstate_vars%cp_scalar               , &
         is_soil               => decomp_cascade_con%is_soil             &
         )

    dt = real( get_step_size(), r8 )

    ! set initial values for potential C and N fluxes
    biochem_pmin_ppools_vr_col(bounds%begc : bounds%endc, :, :) = 0._r8
    biochem_pmin_to_plant_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp) = 0._r8
    biochem_pmin_to_plant_patch(bounds%begp:bounds%endp) = 0._r8
      
    do j = 1,nlevdecomp
        do fc = 1,num_soilc
            c = filter_soilc(fc)
            biochem_pmin_vr(c,j) = 0.0_r8
            biochem_pmin_to_ecosysp_vr_col_pot(c,j) = 0._r8
            do p = col_pp%pfti(c), col_pp%pftf(c)
                if (veg_pp%active(p).and. (veg_pp%itype(p) .ne. noveg)) then
                    !lamda_up = npimbalance(p) ! partial_vcmax/partial_lpc / partial_vcmax/partial_lnc
                    lamda_up = cp_scalar(p)/max(cn_scalar(p),1e-20_r8)
                    lamda_up = min(max(lamda_up,0.0_r8), 150.0_r8)
                    ptase_tmp = vmax_ptase(veg_pp%itype(p)) * froot_prof(p,j) * max(lamda_up - lamda_ptase, 0.0_r8) / &
                        (km_ptase + max(lamda_up - lamda_ptase, 0.0_r8)) 
                    if (NFIX_PTASE_plant) then
                       biochem_pmin_to_plant_vr_patch(p,j) = ptase_tmp * alpha_ptase(veg_pp%itype(p)) 
                       biochem_pmin_vr(c,j) = biochem_pmin_vr(c,j) + ptase_tmp * veg_pp%wtcol(p) * &
                            (1._r8 - alpha_ptase(veg_pp%itype(p)))
                       biochem_pmin_to_ecosysp_vr_col_pot(c,j) = biochem_pmin_to_ecosysp_vr_col_pot(c,j) + ptase_tmp  * veg_pp%wtcol(p)
                    else
                       biochem_pmin_to_plant_vr_patch(p,j) = 0._r8
                       biochem_pmin_vr(c,j) = biochem_pmin_vr(c,j) + ptase_tmp * veg_pp%wtcol(p)
                       biochem_pmin_to_ecosysp_vr_col_pot(c,j) = biochem_pmin_to_ecosysp_vr_col_pot(c,j) + &
                            ptase_tmp  * veg_pp%wtcol(p)
                    endif
                end if
            enddo
        enddo
    enddo 
    
    do j = 1,nlevdecomp
        do fc = 1,num_soilc
            c = filter_soilc(fc)
            ! sum total
            sop_tot = 0._r8
            do l = 1,ndecomp_pools
              if (is_soil(l)) then
                sop_tot = sop_tot + decomp_ppools_vr_col(c,j,l)
              end if
            end do
            ! get profile
            do l = 1,ndecomp_pools
              if (is_soil(l)) then
                if (sop_tot > 1e-12) then 
                    sop_profile(l) = decomp_ppools_vr_col(c,j,l)/sop_tot
                else
                    sop_profile(l) = 0._r8
                end if
              end if
            end do
            ! cauculation actual biochem_pmin_ppool
            do l = 1,ndecomp_pools
                if (is_soil(l)) then
                   biochem_pmin_ppools_vr_col(c,j,l) = max(min(biochem_pmin_to_ecosysp_vr_col_pot(c,j) * sop_profile(l)&
                        , decomp_ppools_vr_col(c,j,l)/dt),0._r8)
                end if
            end do
        end do
    end do
    ! update biochem_pmin_to_ecosysp_vr_col,biochem_pmin_vr,biochem_pmin_to_plant_vr
    do j = 1,nlevdecomp
        do fc = 1,num_soilc
            c = filter_soilc(fc)
            biochem_pmin_to_ecosysp_vr_col(c,j)=0._r8
            do l = 1, ndecomp_pools
               if (is_soil(l)) then
                     biochem_pmin_to_ecosysp_vr_col(c,j) = biochem_pmin_to_ecosysp_vr_col(c,j)+ &
                                          biochem_pmin_ppools_vr_col(c,j,l)
               end if
            enddo
        enddo
    end do
    if (NFIX_PTASE_plant) then   
      ! rescale biochem_pmin_vr, biochem_pmin_to_plant_vr if necessary
      do j = 1,nlevdecomp
        do fc = 1,num_soilc
            c = filter_soilc(fc)
               if ( (biochem_pmin_to_ecosysp_vr_col_pot(c,j) > biochem_pmin_to_ecosysp_vr_col(c,j)) ) then
                  if ( biochem_pmin_to_ecosysp_vr_col_pot(c,j) > 0.0_r8 ) then
                     biochem_pmin_vr(c,j) = biochem_pmin_vr(c,j) * &
                        biochem_pmin_to_ecosysp_vr_col(c,j) / biochem_pmin_to_ecosysp_vr_col_pot(c,j)
                     do p = col_pp%pfti(c), col_pp%pftf(c)
                        biochem_pmin_to_plant_vr_patch(p,j) = biochem_pmin_to_plant_vr_patch(p,j) * &
                           biochem_pmin_to_ecosysp_vr_col(c,j) / biochem_pmin_to_ecosysp_vr_col_pot(c,j)
                     end do
                  else
                     biochem_pmin_vr(c,j) = 0.0_r8
                     do p = col_pp%pfti(c), col_pp%pftf(c)
                        biochem_pmin_to_plant_vr_patch(p,j) = 0.0_r8
                     end do
                  end if
               end if
        end do
      end do 
      ! sum up biochem_pmin_to_plant
      do fc = 1,num_soilc
       c = filter_soilc(fc)
       do p = col_pp%pfti(c), col_pp%pftf(c)
        if (veg_pp%active(p).and. (veg_pp%itype(p) .ne. noveg)) then
          !biochem_pmin_to_plant_patch(p) = 0._r8
          do j = 1,nlevdecomp
             biochem_pmin_to_plant_patch(p) = biochem_pmin_to_plant_patch(p) + &
                                              biochem_pmin_to_plant_vr_patch(p,j) * col_pp%dz(c,j)
          end do
        end if
       end do
      end do
    else
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             biochem_pmin_vr(c,j) = biochem_pmin_to_ecosysp_vr_col(c,j)
          enddo
       end do
    end if

    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          do l = 1, ndecomp_pools
             decomp_ppools_vr_col(c,j,l) = decomp_ppools_vr_col(c,j,l)- biochem_pmin_ppools_vr_col(c,j,l)*dt
          end do
       end do
    end do

    end associate

  end subroutine PhosphorusBiochemMin_balance

end module PhosphorusDynamicsMod
