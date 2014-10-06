module CanopyHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculation of
  ! (1) water storage of intercepted precipitation
  ! (2) direct throughfall and canopy drainage of precipitation
  ! (3) the fraction of foliage covered by water and the fraction
  !     of foliage that is dry and transpiring.
  ! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_sys_mod     , only : shr_sys_flush
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varctl      , only : iulog
  use LandunitType    , only : lun                
  use atm2lndType     , only : atm2lnd_type
  use AerosolType     , only : aerosol_type
  use CanopyStateType , only : canopystate_type
  use TemperatureType , only : temperature_type
  use WaterfluxType   , only : waterflux_type
  use WaterstateType  , only : waterstate_type
  use ColumnType      , only : col                
  use PatchType       , only : pft                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyHydrology_readnl ! Read namelist
  public :: CanopyHydrology        ! Run
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: FracWet    ! Determine fraction of vegetated surface that is wet
  private :: FracH2oSfc ! Determine fraction of land surfaces which are submerged  
  !
  ! !PRIVATE DATA MEMBERS:
  integer :: oldfflag=0  ! use old fsno parameterization (N&Y07) 
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CanopyHydrology_readnl( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for CanopyHydrology
    !
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'CanopyHydrology_readnl'  ! subroutine name
    !-----------------------------------------------------------------------

    namelist / clm_canopyhydrology_inparm / oldfflag

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input. 
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in clm_CanopyHydrology_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_CanopyHydrology_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_canopyhydrology_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading clm_canopyhydrology_inparm namelist"//errmsg(__FILE__, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if
    ! Broadcast namelist variables read in
    call shr_mpi_bcast(oldfflag, mpicom)

   end subroutine CanopyHydrology_readnl

   !-----------------------------------------------------------------------
   subroutine CanopyHydrology(bounds, &
        num_nolakec, filter_nolakec, num_nolakep, filter_nolakep, &
        atm2lnd_vars, canopystate_vars, temperature_vars, &
        aerosol_vars, waterstate_vars, waterflux_vars)
     !
     ! !DESCRIPTION:
     ! Calculation of
     ! (1) water storage of intercepted precipitation
     ! (2) direct throughfall and canopy drainage of precipitation
     ! (3) the fraction of foliage covered by water and the fraction
     !     of foliage that is dry and transpiring.
     ! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
     ! Note:  The evaporation loss is taken off after the calculation of leaf
     ! temperature in the subroutine clm\_leaftem.f90, not in this subroutine.
     !
     ! !USES:
     use clm_varcon         , only : hfus, denice, zlnd, rpi, spval, tfrz
     use column_varcon      , only : icol_roof, icol_sunwall, icol_shadewall
     use landunit_varcon    , only : istcrop, istice, istwet, istsoil, istice_mec 
     use clm_varctl         , only : subgridflag
     use clm_varpar         , only : nlevsoi,nlevsno
     use clm_time_manager   , only : get_step_size
     use subgridAveMod      , only : p2c
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds     
     integer                , intent(in)    :: num_nolakec          ! number of column non-lake points in column filter
     integer                , intent(in)    :: filter_nolakec(:)    ! column filter for non-lake points
     integer                , intent(in)    :: num_nolakep          ! number of pft non-lake points in pft filter
     integer                , intent(in)    :: filter_nolakep(:)    ! patch filter for non-lake points
     type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
     type(canopystate_type) , intent(in)    :: canopystate_vars
     type(temperature_type) , intent(inout) :: temperature_vars
     type(aerosol_type)     , intent(inout) :: aerosol_vars
     type(waterstate_type)  , intent(inout) :: waterstate_vars
     type(waterflux_type)   , intent(inout) :: waterflux_vars
     !
     ! !LOCAL VARIABLES:
     integer  :: f                                            ! filter index
     integer  :: pi                                           ! patch index
     integer  :: p                                            ! patch index
     integer  :: c                                            ! column index
     integer  :: l                                            ! landunit index
     integer  :: g                                            ! gridcell index
     integer  :: newnode                                      ! flag when new snow node is set, (1=yes, 0=no)
     real(r8) :: dtime                                        ! land model time step (sec)
     real(r8) :: h2ocanmx                                     ! maximum allowed water on canopy [mm]
     real(r8) :: fpi                                          ! coefficient of interception
     real(r8) :: xrun                                         ! excess water that exceeds the leaf capacity [mm/s]
     real(r8) :: dz_snowf                                     ! layer thickness rate change due to precipitation [mm/s]
     real(r8) :: bifall                                       ! bulk density of newly fallen dry snow [kg/m3]
     real(r8) :: fracsnow(bounds%begp:bounds%endp)            ! frac of precipitation that is snow
     real(r8) :: fracrain(bounds%begp:bounds%endp)            ! frac of precipitation that is rain
     real(r8) :: qflx_candrip(bounds%begp:bounds%endp)        ! rate of canopy runoff and snow falling off canopy [mm/s]
     real(r8) :: qflx_through_rain(bounds%begp:bounds%endp)   ! direct rain throughfall [mm/s]
     real(r8) :: qflx_through_snow(bounds%begp:bounds%endp)   ! direct snow throughfall [mm/s]
     real(r8) :: qflx_prec_grnd_snow(bounds%begp:bounds%endp) ! snow precipitation incident on ground [mm/s]
     real(r8) :: qflx_prec_grnd_rain(bounds%begp:bounds%endp) ! rain precipitation incident on ground [mm/s]
     real(r8) :: z_avg                                        ! grid cell average snow depth
     real(r8) :: rho_avg                                      ! avg density of snow column
     real(r8) :: temp_snow_depth,temp_intsnow                 ! temporary variables
     real(r8) :: fmelt
     real(r8) :: smr
     real(r8) :: delf_melt
     real(r8) :: fsno_new
     real(r8) :: accum_factor
     real(r8) :: newsnow(bounds%begc:bounds%endc)
     real(r8) :: snowmelt(bounds%begc:bounds%endc)
     integer  :: j
     !-----------------------------------------------------------------------

     associate(                                                     & 
          pgridcell            => pft%gridcell                             , & ! Input:  [integer  (:)   ]  pft's gridcell                           
          plandunit            => pft%landunit                             , & ! Input:  [integer  (:)   ]  pft's landunit                           
          pcolumn              => pft%column                               , & ! Input:  [integer  (:)   ]  pft's column                             
          ltype                => lun%itype                                , & ! Input:  [integer  (:)   ]  landunit type                            
          urbpoi               => lun%urbpoi                               , & ! Input:  [logical  (:)   ]  true => landunit is an urban point       
          cgridcell            => col%gridcell                             , & ! Input:  [integer  (:)   ]  columns's gridcell                       
          clandunit            => col%landunit                             , & ! Input:  [integer  (:)   ]  columns's landunit                       
          ctype                => col%itype                                , & ! Input:  [integer  (:)   ]  column type                              
          pfti                 => col%pfti                                 , & ! Input:  [integer  (:)   ]  column's beginning pft index             
          npfts                => col%npfts                                , & ! Input:  [integer  (:)   ]  number of patches in column                 
          snl                  => col%snl                                  , & ! Input:  [integer  (:)   ]  number of snow layers                    
          n_melt               => col%n_melt                               , & ! Input:  [real(r8) (:)   ]  SCA shape parameter                     
          zi                   => col%zi                                   , & ! Output: [real(r8) (:,:) ]  interface level below a "z" level (m) 
          dz                   => col%dz                                   , & ! Output: [real(r8) (:,:) ]  layer depth (m)                       
          z                    => col%z                                    , & ! Output: [real(r8) (:,:) ]  layer thickness (m)                   

          forc_rain            => atm2lnd_vars%forc_rain_downscaled_col    , & ! Input:  [real(r8) (:)   ]  rain rate [mm/s]                        
          forc_snow            => atm2lnd_vars%forc_snow_downscaled_col    , & ! Input:  [real(r8) (:)   ]  snow rate [mm/s]                        
          forc_t               => atm2lnd_vars%forc_t_downscaled_col       , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)        
          qflx_floodg          => atm2lnd_vars%forc_flood_grc              , & ! Input:  [real(r8) (:)   ]  gridcell flux of flood water from RTM   

          dewmx                => canopystate_vars%dewmx_patch             , & ! Input:  [real(r8) (:)   ]  Maximum allowed dew [mm]                
          frac_veg_nosno       => canopystate_vars%frac_veg_nosno_patch    , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
          elai                 => canopystate_vars%elai_patch              , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow
          esai                 => canopystate_vars%esai_patch              , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow

          t_grnd               => temperature_vars%t_grnd_col              , & ! Input:  [real(r8) (:)   ]  ground temperature (Kelvin)             
          t_soisno             => temperature_vars%t_soisno_col            , & ! Output: [real(r8) (:,:) ]  soil temperature (Kelvin)  

          do_capsnow           => waterstate_vars%do_capsnow_col           , & ! Output: [logical  (:)   ]  true => do snow capping                  
          h2ocan               => waterstate_vars%h2ocan_patch             , & ! Output: [real(r8) (:)   ]  total canopy water (mm H2O)             
          h2osfc               => waterstate_vars%h2osfc_col               , & ! Output: [real(r8) (:)   ]  surface water (mm)                      
          h2osno               => waterstate_vars%h2osno_col               , & ! Output: [real(r8) (:)   ]  snow water (mm H2O)                     
          snow_depth           => waterstate_vars%snow_depth_col           , & ! Output: [real(r8) (:)   ]  snow height (m)                         
          int_snow             => waterstate_vars%int_snow_col             , & ! Output: [real(r8) (:)   ]  integrated snowfall [mm]                
          frac_sno_eff         => waterstate_vars%frac_sno_eff_col         , & ! Output: [real(r8) (:)   ]  eff. fraction of ground covered by snow (0 to 1)
          frac_sno             => waterstate_vars%frac_sno_col             , & ! Output: [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
          frac_h2osfc          => waterstate_vars%frac_h2osfc_col          , & ! Output: [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
          frac_iceold          => waterstate_vars%frac_iceold_col          , & ! Output: [real(r8) (:,:) ]  fraction of ice relative to the tot water
          h2osoi_ice           => waterstate_vars%h2osoi_ice_col           , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                      
          h2osoi_liq           => waterstate_vars%h2osoi_liq_col           , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                  
          swe_old              => waterstate_vars%swe_old_col              , & ! Output: [real(r8) (:,:) ]  snow water before update              

          irrig_rate           => waterflux_vars%irrig_rate_patch          , & ! Input:  [real(r8) (:)   ]  current irrigation rate (applied if n_irrig_steps_left > 0) [mm/s]
          n_irrig_steps_left   => waterflux_vars%n_irrig_steps_left_patch  , & ! Output: [integer  (:)   ]  number of time steps for which we still need to irrigate today
          qflx_floodc          => waterflux_vars%qflx_floodc_col           , & ! Output: [real(r8) (:)   ]  column flux of flood water from RTM     
          qflx_snow_melt       => waterflux_vars%qflx_snow_melt_col        , & ! Output: [real(r8) (:)   ]  snow melt from previous time step       
          qflx_snow_h2osfc     => waterflux_vars%qflx_snow_h2osfc_col      , & ! Output: [real(r8) (:)   ]  snow falling on surface water (mm/s)     
          qflx_snwcp_liq       => waterflux_vars%qflx_snwcp_liq_patch      , & ! Output: [real(r8) (:)   ]  excess rainfall due to snow capping (mm H2O /s) [+]
          qflx_snwcp_ice       => waterflux_vars%qflx_snwcp_ice_patch      , & ! Output: [real(r8) (:)   ]  excess snowfall due to snow capping (mm H2O /s) [+]
          qflx_snow_grnd_col   => waterflux_vars%qflx_snow_grnd_col        , & ! Output: [real(r8) (:)   ]  snow on ground after interception (mm H2O/s) [+]
          qflx_snow_grnd_patch => waterflux_vars%qflx_snow_grnd_patch      , & ! Output: [real(r8) (:)   ]  snow on ground after interception (mm H2O/s) [+]
          qflx_prec_intr       => waterflux_vars%qflx_prec_intr_patch      , & ! Output: [real(r8) (:)   ]  interception of precipitation [mm/s]    
          qflx_prec_grnd       => waterflux_vars%qflx_prec_grnd_patch      , & ! Output: [real(r8) (:)   ]  water onto ground including canopy runoff [kg/(m2 s)]
          qflx_rain_grnd       => waterflux_vars%qflx_rain_grnd_patch      , & ! Output: [real(r8) (:)   ]  rain on ground after interception (mm H2O/s) [+]
          qflx_irrig           => waterflux_vars%qflx_irrig_patch            & ! Output: [real(r8) (:)   ]  irrigation amount (mm/s)                
          )

       ! Compute time step
       
       dtime = get_step_size()

       ! Start pft loop

       do f = 1, num_nolakep
          p = filter_nolakep(f)
          g = pgridcell(p)
          l = plandunit(p)
          c = pcolumn(p)

          ! Canopy interception and precipitation onto ground surface
          ! Add precipitation to leaf water

          if (ltype(l)==istsoil .or. ltype(l)==istwet .or. urbpoi(l) .or. &
               ltype(l)==istcrop) then

             qflx_candrip(p) = 0._r8      ! rate of canopy runoff
             qflx_through_snow(p) = 0._r8 ! rain precipitation direct through canopy
             qflx_through_rain(p) = 0._r8 ! snow precipitation direct through canopy
             qflx_prec_intr(p) = 0._r8    ! total intercepted precipitation
             fracsnow(p) = 0._r8          ! fraction of input precip that is snow
             fracrain(p) = 0._r8          ! fraction of input precip that is rain


             if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall) then

                if (frac_veg_nosno(p) == 1 .and. (forc_rain(c) + forc_snow(c)) > 0._r8) then

                   ! determine fraction of input precipitation that is snow and rain
                   fracsnow(p) = forc_snow(c)/(forc_snow(c) + forc_rain(c))
                   fracrain(p) = forc_rain(c)/(forc_snow(c) + forc_rain(c))

                   ! The leaf water capacities for solid and liquid are different,
                   ! generally double for snow, but these are of somewhat less
                   ! significance for the water budget because of lower evap. rate at
                   ! lower temperature.  Hence, it is reasonable to assume that
                   ! vegetation storage of solid water is the same as liquid water.
                   h2ocanmx = dewmx(p) * (elai(p) + esai(p))

                   ! Coefficient of interception
                   ! set fraction of potential interception to max 0.25
                   fpi = 0.25_r8*(1._r8 - exp(-0.5_r8*(elai(p) + esai(p))))

                   ! Direct throughfall
                   qflx_through_snow(p) = forc_snow(c) * (1._r8-fpi)
                   qflx_through_rain(p) = forc_rain(c) * (1._r8-fpi)

                   ! Intercepted precipitation [mm/s]
                   qflx_prec_intr(p) = (forc_snow(c) + forc_rain(c)) * fpi

                   ! Water storage of intercepted precipitation and dew
                   h2ocan(p) = max(0._r8, h2ocan(p) + dtime*qflx_prec_intr(p))

                   ! Initialize rate of canopy runoff and snow falling off canopy
                   qflx_candrip(p) = 0._r8

                   ! Excess water that exceeds the leaf capacity
                   xrun = (h2ocan(p) - h2ocanmx)/dtime

                   ! Test on maximum dew on leaf
                   ! Note if xrun > 0 then h2ocan must be at least h2ocanmx
                   if (xrun > 0._r8) then
                      qflx_candrip(p) = xrun
                      h2ocan(p) = h2ocanmx
                   end if

                end if
             end if

          else if (ltype(l)==istice .or. ltype(l)==istice_mec) then

             h2ocan(p)            = 0._r8
             qflx_candrip(p)      = 0._r8
             qflx_through_snow(p) = 0._r8
             qflx_through_rain(p) = 0._r8
             qflx_prec_intr(p)    = 0._r8
             fracsnow(p)          = 0._r8
             fracrain(p)          = 0._r8

          end if

          ! Precipitation onto ground (kg/(m2 s))

          if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall) then
             if (frac_veg_nosno(p) == 0) then
                qflx_prec_grnd_snow(p) = forc_snow(c)
                qflx_prec_grnd_rain(p) = forc_rain(c)
             else
                qflx_prec_grnd_snow(p) = qflx_through_snow(p) + (qflx_candrip(p) * fracsnow(p))
                qflx_prec_grnd_rain(p) = qflx_through_rain(p) + (qflx_candrip(p) * fracrain(p))
             end if
             ! Urban sunwall and shadewall have no intercepted precipitation
          else
             qflx_prec_grnd_snow(p) = 0.
             qflx_prec_grnd_rain(p) = 0.
          end if

          ! Determine whether we're irrigating here; set qflx_irrig appropriately
          if (n_irrig_steps_left(p) > 0) then
             qflx_irrig(p)         = irrig_rate(p)
             n_irrig_steps_left(p) = n_irrig_steps_left(p) - 1
          else
             qflx_irrig(p) = 0._r8
          end if

          ! Add irrigation water directly onto ground (bypassing canopy interception)
          ! Note that it's still possible that (some of) this irrigation water will runoff (as runoff is computed later)
          qflx_prec_grnd_rain(p) = qflx_prec_grnd_rain(p) + qflx_irrig(p)

          ! Done irrigation

          qflx_prec_grnd(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)

          if (do_capsnow(c)) then
             qflx_snwcp_liq(p) = qflx_prec_grnd_rain(p)
             qflx_snwcp_ice(p) = qflx_prec_grnd_snow(p)

             qflx_snow_grnd_patch(p) = 0._r8
             qflx_rain_grnd(p) = 0._r8
          else
             qflx_snwcp_liq(p) = 0._r8
             qflx_snwcp_ice(p) = 0._r8
             qflx_snow_grnd_patch(p) = qflx_prec_grnd_snow(p)           ! ice onto ground (mm/s)
             qflx_rain_grnd(p)     = qflx_prec_grnd_rain(p)           ! liquid water onto ground (mm/s)
          end if

       end do ! (end pft loop)

       ! Determine the fraction of foliage covered by water and the
       ! fraction of foliage that is dry and transpiring.

       call FracWet(num_nolakep, filter_nolakep, &
            canopystate_vars, waterstate_vars)

       ! Update column level state variables for snow.

       call p2c(bounds, num_nolakec, filter_nolakec, &
            qflx_snow_grnd_patch(bounds%begp:bounds%endp), &
            qflx_snow_grnd_col(bounds%begc:bounds%endc))

       ! apply gridcell flood water flux to non-lake columns
       do f = 1, num_nolakec
          c = filter_nolakec(f)
          g = cgridcell(c)
          if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall) then      
             qflx_floodc(c) = qflx_floodg(g)
          else
             qflx_floodc(c) = 0._r8
          endif
       enddo

       ! Determine snow height and snow water

       do f = 1, num_nolakec
          c = filter_nolakec(f)
          l = clandunit(c)
          g = cgridcell(c)

          ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
          ! U.S.Department of Agriculture Forest Service, Project F,
          ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

          qflx_snow_h2osfc(c) = 0._r8
          ! set temporary variables prior to updating
          temp_snow_depth=snow_depth(c)
          ! save initial snow content
          do j= -nlevsno+1,snl(c)
             swe_old(c,j) = 0.0_r8
          end do
          do j= snl(c)+1,0
             swe_old(c,j)=h2osoi_liq(c,j)+h2osoi_ice(c,j)
          enddo

          if (do_capsnow(c)) then
             dz_snowf = 0._r8
             newsnow(c) = (1._r8 - frac_h2osfc(c)) * qflx_snow_grnd_col(c) * dtime
             frac_sno(c)=1._r8
             int_snow(c) = 5.e2_r8
          else
             if (forc_t(c) > tfrz + 2._r8) then
                bifall=50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
             else if (forc_t(c) > tfrz - 15._r8) then
                bifall=50._r8 + 1.7_r8*(forc_t(c) - tfrz + 15._r8)**1.5_r8
             else
                bifall=50._r8
             end if

             ! newsnow is all snow that doesn't fall on h2osfc
             newsnow(c) = (1._r8 - frac_h2osfc(c)) * qflx_snow_grnd_col(c) * dtime

             ! update int_snow
             int_snow(c) = max(int_snow(c),h2osno(c)) !h2osno could be larger due to frost

             ! snowmelt from previous time step * dtime
             snowmelt(c) = qflx_snow_melt(c) * dtime

             ! set shape factor for accumulation of snow
             accum_factor=0.1

             if (h2osno(c) > 0.0) then

                !======================  FSCA PARAMETERIZATIONS  ======================
                ! fsca parameterization based on *changes* in swe
                ! first compute change from melt during previous time step
                if(snowmelt(c) > 0._r8) then

                   smr=min(1._r8,(h2osno(c))/(int_snow(c)))

                   frac_sno(c) = 1. - (acos(min(1._r8,(2.*smr - 1._r8)))/rpi)**(n_melt(c))

                endif

                ! update fsca by new snow event, add to previous fsca
                if (newsnow(c) > 0._r8) then
                   fsno_new = 1._r8 - (1._r8 - tanh(accum_factor*newsnow(c)))*(1._r8 - frac_sno(c))
                   frac_sno(c) = fsno_new

                   ! reset int_snow after accumulation events
                   temp_intsnow= (h2osno(c) + newsnow(c)) &
                        / (0.5*(cos(rpi*(1._r8-max(frac_sno(c),1e-6_r8))**(1./n_melt(c)))+1._r8))
                   int_snow(c) = min(1.e8_r8,temp_intsnow)
                endif

                !====================================================================

                ! for subgrid fluxes
                if (subgridflag ==1 .and. .not. urbpoi(l)) then
                   if (frac_sno(c) > 0._r8)then
                      snow_depth(c)=snow_depth(c) + newsnow(c)/(bifall * frac_sno(c))
                   else
                      snow_depth(c)=0._r8
                   end if
                else
                   ! for uniform snow cover
                   snow_depth(c)=snow_depth(c)+newsnow(c)/bifall
                endif

                ! use original fsca formulation (n&y 07)
                if (oldfflag == 1) then 
                   ! snow cover fraction in Niu et al. 2007
                   if(snow_depth(c) > 0.0_r8)  then
                      frac_sno(c) = tanh(snow_depth(c)/(2.5_r8*zlnd* &
                           (min(800._r8,(h2osno(c)+ newsnow(c))/snow_depth(c))/100._r8)**1._r8) )
                   endif
                   if(h2osno(c) < 1.0_r8)  then
                      frac_sno(c)=min(frac_sno(c),h2osno(c))
                   endif
                endif

             else !h2osno == 0
                ! initialize frac_sno and snow_depth when no snow present initially
                if (newsnow(c) > 0._r8) then 
                   z_avg = newsnow(c)/bifall
                   fmelt=newsnow(c)
                   frac_sno(c) = tanh(accum_factor*newsnow(c))

                   ! make int_snow consistent w/ new fsno, h2osno
                   int_snow(c) = 0. !reset prior to adding newsnow below
                   temp_intsnow= (h2osno(c) + newsnow(c)) &
                        / (0.5*(cos(rpi*(1._r8-max(frac_sno(c),1e-6_r8))**(1./n_melt(c)))+1._r8))
                   int_snow(c) = min(1.e8_r8,temp_intsnow)

                   ! update snow_depth and h2osno to be consistent with frac_sno, z_avg
                   if (subgridflag ==1 .and. .not. urbpoi(l)) then
                      snow_depth(c)=z_avg/frac_sno(c)
                   else
                      snow_depth(c)=newsnow(c)/bifall
                   endif
                   ! use n&y07 formulation
                   if (oldfflag == 1) then 
                      ! snow cover fraction in Niu et al. 2007
                      if(snow_depth(c) > 0.0_r8)  then
                         frac_sno(c) = tanh(snow_depth(c)/(2.5_r8*zlnd* &
                              (min(800._r8,newsnow(c)/snow_depth(c))/100._r8)**1._r8) )
                      endif
                   endif
                else
                   z_avg = 0._r8
                   snow_depth(c) = 0._r8
                   frac_sno(c) = 0._r8
                endif
             endif ! end of h2osno > 0

             ! snow directly falling on surface water melts, increases h2osfc
             qflx_snow_h2osfc(c) = frac_h2osfc(c)*qflx_snow_grnd_col(c)

             ! update h2osno for new snow
             h2osno(c) = h2osno(c) + newsnow(c) 
             int_snow(c) = int_snow(c) + newsnow(c)

             ! update change in snow depth
             dz_snowf = (snow_depth(c) - temp_snow_depth) / dtime

          end if !end of do_capsnow construct

          ! set frac_sno_eff variable
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
             if (subgridflag ==1) then 
                frac_sno_eff(c) = frac_sno(c)
             else
                frac_sno_eff(c) = 1._r8
             endif
          else
             frac_sno_eff(c) = 1._r8
          endif

          if (ltype(l)==istwet .and. t_grnd(c)>tfrz) then
             h2osno(c)=0._r8
             snow_depth(c)=0._r8
          end if

          ! When the snow accumulation exceeds 10 mm, initialize snow layer
          ! Currently, the water temperature for the precipitation is simply set
          ! as the surface air temperature

          newnode = 0    ! flag for when snow node will be initialized
          if (snl(c) == 0 .and. qflx_snow_grnd_col(c) > 0.0_r8 .and. frac_sno(c)*snow_depth(c) >= 0.01_r8) then
             newnode = 1
             snl(c) = -1
             dz(c,0) = snow_depth(c)                       ! meter
             z(c,0) = -0.5_r8*dz(c,0)
             zi(c,-1) = -dz(c,0)
             t_soisno(c,0) = min(tfrz, forc_t(c))      ! K
             h2osoi_ice(c,0) = h2osno(c)               ! kg/m2
             h2osoi_liq(c,0) = 0._r8                   ! kg/m2
             frac_iceold(c,0) = 1._r8

             ! intitialize SNICAR variables for fresh snow:
             call aerosol_vars%Reset(column=c)
             call waterstate_vars%Reset(column=c)
          end if

          ! The change of ice partial density of surface node due to precipitation.
          ! Only ice part of snowfall is added here, the liquid part will be added
          ! later.

          if (snl(c) < 0 .and. newnode == 0) then
             h2osoi_ice(c,snl(c)+1) = h2osoi_ice(c,snl(c)+1)+newsnow(c)
             dz(c,snl(c)+1) = dz(c,snl(c)+1)+dz_snowf*dtime
          end if

       end do

       ! update surface water fraction (this may modify frac_sno)
       call FracH2oSfc(bounds, num_nolakec, filter_nolakec, &
            waterstate_vars)

     end associate 

   end subroutine CanopyHydrology

   !-----------------------------------------------------------------------
   subroutine FracWet(numf, filter, canopystate_vars, waterstate_vars)
     !
     ! !DESCRIPTION:
     ! Determine fraction of vegetated surfaces which are wet and
     ! fraction of elai which is dry. The variable ``fwet'' is the
     ! fraction of all vegetation surfaces which are wet including
     ! stem area which contribute to evaporation. The variable ``fdry''
     ! is the fraction of elai which is dry because only leaves
     ! can transpire.  Adjusted for stem area which does not transpire.
     !
     ! !ARGUMENTS:
     integer                , intent(in)    :: numf                  ! number of filter non-lake points
     integer                , intent(in)    :: filter(numf)          ! patch filter for non-lake points
     type(canopystate_type) , intent(in)    :: canopystate_vars
     type(waterstate_type)  , intent(inout) :: waterstate_vars
     !
     ! !LOCAL VARIABLES:
     integer  :: fp,p             ! indices
     real(r8) :: vegt             ! lsai
     real(r8) :: dewmxi           ! inverse of maximum allowed dew [1/mm]
     !-----------------------------------------------------------------------

     associate(                                              & 
          frac_veg_nosno => canopystate_vars%frac_veg_nosno_patch , & ! Input:  [integer (:)]  fraction of veg not covered by snow (0/1 now) [-]
          dewmx          => canopystate_vars%dewmx_patch          , & ! Input:  [real(r8) (:) ]  Maximum allowed dew [mm]                
          elai           => canopystate_vars%elai_patch           , & ! Input:  [real(r8) (:) ]  one-sided leaf area index with burying by snow
          esai           => canopystate_vars%esai_patch           , & ! Input:  [real(r8) (:) ]  one-sided stem area index with burying by snow

          h2ocan         => waterstate_vars%h2ocan_patch          , & ! Input:  [real(r8) (:) ]  total canopy water (mm H2O)             
          
          fwet           => waterstate_vars%fwet_patch            , & ! Output: [real(r8) (:) ]  fraction of canopy that is wet (0 to 1) 
          fdry           => waterstate_vars%fdry_patch              & ! Output: [real(r8) (:) ]  fraction of foliage that is green and dry [-] (new)
          )

       do fp = 1,numf
          p = filter(fp)
          if (frac_veg_nosno(p) == 1) then
             if (h2ocan(p) > 0._r8) then
                vegt    = frac_veg_nosno(p)*(elai(p) + esai(p))
                dewmxi  = 1.0_r8/dewmx(p)
                fwet(p) = ((dewmxi/vegt)*h2ocan(p))**0.666666666666_r8
                fwet(p) = min (fwet(p),1.0_r8)   ! Check for maximum limit of fwet
             else
                fwet(p) = 0._r8
             end if
             fdry(p) = (1._r8-fwet(p))*elai(p)/(elai(p)+esai(p))
          else
             fwet(p) = 0._r8
             fdry(p) = 0._r8
          end if
       end do

     end associate
   end subroutine FracWet

   !-----------------------------------------------------------------------
   subroutine FracH2OSfc(bounds, num_h2osfc, filter_h2osfc, &
        waterstate_vars, no_update)
     !
     ! !DESCRIPTION:
     ! Determine fraction of land surfaces which are submerged  
     ! based on surface microtopography and surface water storage.
     !
     ! !USES:
     use shr_const_mod   , only : shr_const_pi
     use shr_spfn_mod    , only : erf => shr_spfn_erf
     use landunit_varcon , only : istsoil, istcrop
     !
     ! !ARGUMENTS:
     type(bounds_type)     , intent(in)           :: bounds           
     integer               , intent(in)           :: num_h2osfc       ! number of column points in column filter
     integer               , intent(in)           :: filter_h2osfc(:) ! column filter 
     type(waterstate_type) , intent(inout)        :: waterstate_vars
     integer               , intent(in), optional :: no_update        ! flag to make calculation w/o updating variables
     !
     ! !LOCAL VARIABLES:
     integer :: c,f,l          ! indices
     real(r8):: d,fd,dfdd      ! temporary variable for frac_h2oscs iteration
     real(r8):: sigma          ! microtopography pdf sigma in mm
     real(r8):: min_h2osfc
     !-----------------------------------------------------------------------

     associate(                                              & 
          micro_sigma  => col%micro_sigma                  , & ! Input:  [real(r8) (:)   ] microtopography pdf sigma (m)                     

          h2osno       => waterstate_vars%h2osno_col       , & ! Input:  [real(r8) (:)   ] snow water (mm H2O)                               
          
          h2osoi_liq   => waterstate_vars%h2osoi_liq_col   , & ! Output: [real(r8) (:,:) ] liquid water (col,lyr) [kg/m2]                  
          h2osfc       => waterstate_vars%h2osfc_col       , & ! Output: [real(r8) (:)   ] surface water (mm)                                
          frac_sno     => waterstate_vars%frac_sno_col     , & ! Output: [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)       
          frac_sno_eff => waterstate_vars%frac_sno_eff_col , & ! Output: [real(r8) (:)   ] eff. fraction of ground covered by snow (0 to 1)  
          frac_h2osfc  => waterstate_vars%frac_h2osfc_col    & ! Output: [real(r8) (:)   ] col fractional area with surface water greater than zero 
          )

       ! arbitrary lower limit on h2osfc for safer numerics...
       min_h2osfc=1.e-8_r8

       do f = 1, num_h2osfc
          c = filter_h2osfc(f)
          l = col%landunit(c)

          ! h2osfc only calculated for soil vegetated land units
          if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

             !  Use newton-raphson method to iteratively determine frac_h20sfc
             !  based on amount of surface water storage (h2osfc) and 
             !  microtopography variability (micro_sigma)

             if (h2osfc(c) > min_h2osfc) then
                ! a cutoff is needed for numerical reasons...(nonconvergence after 5 iterations)
                d=0.0

                sigma=1.0e3 * micro_sigma(c) ! convert to mm
                do l=1,10
                   fd = 0.5*d*(1.0_r8+erf(d/(sigma*sqrt(2.0)))) &
                        +sigma/sqrt(2.0*shr_const_pi)*exp(-d**2/(2.0*sigma**2)) &
                        -h2osfc(c)
                   dfdd = 0.5*(1.0_r8+erf(d/(sigma*sqrt(2.0))))

                   d = d - fd/dfdd
                enddo
                !--  update the submerged areal fraction using the new d value
                frac_h2osfc(c) = 0.5*(1.0_r8+erf(d/(sigma*sqrt(2.0))))

             else
                frac_h2osfc(c) = 0._r8
                h2osoi_liq(c,1) = h2osoi_liq(c,1) + h2osfc(c)
                h2osfc(c)=0._r8
             endif

             if (.not. present(no_update)) then

                ! adjust fh2o, fsno when sum is greater than zero
                if (frac_sno(c) > (1._r8 - frac_h2osfc(c)) .and. h2osno(c) > 0) then

                   if (frac_h2osfc(c) > 0.01_r8) then             
                      frac_h2osfc(c) = max(1.0_r8 - frac_sno(c),0.01_r8)
                      frac_sno(c) = 1.0_r8 - frac_h2osfc(c)
                   else
                      frac_sno(c) = 1.0_r8 - frac_h2osfc(c)
                   endif
                   frac_sno_eff(c)=frac_sno(c)

                endif

             endif ! end of no_update construct

          else !if landunit not istsoil/istcrop, set frac_h2osfc to zero

             frac_h2osfc(c) = 0._r8

          endif

       end do

     end associate

   end subroutine FracH2OSfc

end module CanopyHydrologyMod
