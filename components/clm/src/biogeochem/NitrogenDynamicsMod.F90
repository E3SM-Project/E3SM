module NitrogenDynamicsMod

  !-----------------------------------------------------------------------
  ! !MODULE: NitrogenDynamicsMod
  !
  ! !DESCRIPTION:
  ! Module for mineral nitrogen dynamics (deposition, fixation, leaching)
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use decompMod           , only : bounds_type
  use clm_varcon          , only : dzsoi_decomp, zisoi
  use clm_varctl          , only : use_nitrif_denitrif, use_vertsoilc, use_fan
  use subgridAveMod       , only : p2c
  use atm2lndType         , only : atm2lnd_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNStateType         , only : cnstate_type
  use WaterStateType      , only : waterstate_type
  use WaterFluxType       , only : waterflux_type
  use CropType            , only : crop_type
  use ColumnType          , only : col_pp
  use ColumnDataType      , only : col_es, col_ws, col_wf, col_cf, col_ns, col_nf 
  use ColumnDataType      , only : nfix_timeconst  
  use VegetationType      , only : veg_pp
  use VegetationDataType  , only : veg_cs, veg_ns, veg_nf  
  use VegetationPropertiesType  , only : veg_vp
  use CNCarbonStateType   , only : carbonstate_type
  use TemperatureType     , only : temperature_type
  use PhosphorusStateType , only : phosphorusstate_type
  use clm_varctl          , only : NFIX_PTASE_plant
  use FrictionVelocityType, only : frictionvel_type
  use SoilStateType       , only : soilstate_type

  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: NitrogenDynamicsInit
  public :: NitrogenDeposition
  public :: NitrogenFixation
  public :: NitrogenLeaching
  public :: NitrogenFert
  public :: CNSoyfix
  public :: readNitrogenDynamicsParams
  public :: NitrogenFixation_balance
  
  !
  ! !PRIVATE DATA:
  type, private :: CNNDynamicsParamsType
     real(r8):: sf        ! soluble fraction of mineral N (unitless)
     real(r8):: sf_no3    ! soluble fraction of NO3 (unitless)
  end type CNNDynamicsParamsType
  
  type(CNNDynamicsParamsType),private ::  CNNDynamicsParamsInst

  logical, private, parameter :: debug_fan = .false.

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine NitrogenDynamicsInit ( )
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
   
  end subroutine NitrogenDynamicsInit

  !-----------------------------------------------------------------------
  subroutine readNitrogenDynamicsParams ( ncid )
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
    
    call NitrogenDynamicsInit()

    tString='sf_minn'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNDynamicsParamsInst%sf=tempr

    tString='sf_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNDynamicsParamsInst%sf_no3=tempr
   
  end subroutine readNitrogenDynamicsParams

  !-----------------------------------------------------------------------
  subroutine NitrogenDeposition( bounds, &
       atm2lnd_vars, nitrogenflux_vars, nitrogenstate_vars, frictionvel_vars, waterstate_vars, waterflux_vars, temperature_vars, &
       soilstate_vars, filter_soilc, num_soilc)
    use GridcellType         , only: grc => grc_pp
    use LandunitType         , only: lun => lun_pp
    use ColumnType           , only : col => col_pp      
    use VegetationType       , only : patch => veg_pp
    use clm_varctl     , only : iulog
    use FanMod
    use abortutils     , only : endrun
    use pftvarcon, only : nc4_grass, nc3_nonarctic_grass
    use landunit_varcon,      only:  istsoil, istcrop
    use clm_varcon, only : spval, ispval
    use clm_time_manager     , only: get_step_size, get_curr_date, get_curr_calday, get_nstep

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
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(frictionvel_type), intent(in) :: frictionvel_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(temperature_type)  , intent(in) :: temperature_vars
    type(soilstate_type), intent(in) :: soilstate_vars
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter

    !
    ! !LOCAL VARIABLES:
    integer :: g,c                    ! indices
    real(r8) :: dt
    !-----------------------------------------------------------------------

    associate(& 
         forc_ndep     =>  atm2lnd_vars%forc_ndep_grc           , & ! Input:  [real(r8) (:)]  nitrogen deposition rate (gN/m2/s)                
         ndep_to_sminn =>  col_nf%ndep_to_sminn   & ! Output: [real(r8) (:)]                                                    
         )

      ! Loop through columns
      do c = bounds%begc, bounds%endc
         g = col_pp%gridcell(c)
         ndep_to_sminn(c) = forc_ndep(g)
      end do

    end associate

     if (use_fan) then
        dt = real( get_step_size(), r8 )
        call fan()
     end if

contains

  subroutine fan()
    integer, parameter :: num_substeps = 4, balance_check_freq = 1000
    integer :: c, g, patchcounter, p, status, c1, c2, l, fc, ind_substep
    real(r8) :: ndep_org(3), orgpools(3), tanprod(3), watertend, fluxes(6,3), tanpools3(3), ratm, tandep, &
         fluxes2(6,2), fluxes3(6,3), fluxes4(6,4), tanpools2(2), tanpools4(4), fluxes_tmp(6), garbage_total
    real(r8), parameter :: water_init_grz = 0.005_r8, cnc_nh3_air = 0.0_r8, depth_slurry = 0.005_r8

    !real(r8), parameter :: fract_resist=0.225_r8, fract_unavail=0.025_r8, fract_avail=0.25_r8, fract_tan=0.6_r8

    real(r8), parameter :: fract_tan=0.6_r8 ! of all N
    real(r8), parameter :: fract_resist=0.45_r8, fract_unavail=0.05_r8, fract_avail=0.5_r8 ! of organic N

    real(r8), parameter :: dz_layer_fert = 0.02_r8, dz_layer_grz = 0.02_r8
    !real(r8), parameter :: fract_resist=0._r8, fract_unavail=0._r8, fract_avail=0._r8, fract_tan=1.0_r8

    real(r8), parameter :: slurry_infiltr_time = 12*3600.0_r8, water_init_fert = 1e-6
    real(r8), parameter :: &
         poolranges_grz(3) = (/24*3600.0_r8, 10*24*3600.0_r8, 360*24*3600.0_r8/), &
         poolranges_fert(3) = (/2*24*3600.0_r8, 24*3600.0_r8, 360*24*3600.0_r8/), &
         poolranges_slr(4) = (/slurry_infiltr_time, 24*3600.0_r8, 10*24*3600.0_r8, 360*24*3600.0_r8/), &
                                !Hconc_grz(3) = (/10**(-8.5_r8), 10**(-8.0_r8), 10**(-7.0_r8)/), &
         Hconc_fert(3) = (/10**(-7.0_r8), 10**(-8.5_r8), 10**(-8.0_r8)/)

    real(r8) :: Hconc_grz(3), Hconc_slr(4), pH_soil, pH_crop

    !logical, parameter :: do_balance_checks = .false.
    logical :: do_balance_checks
    real(r8) :: tg, garbage, theta, thetasat, infiltr_m_s, evap_m_s, runoff_m_s, org_n_tot, &
         nstored_old, nsoilman_old, nsoilfert_old, fert_to_air, fert_to_soil, fert_total, fert_urea, fert_tan, &
         soilflux_org, urea_resid
    real(r8) :: tanprod_from_urea(3), ureapools(2), fert_no3, fert_generic
    !real(r8), parameter :: fract_urea=0.545, fract_no3=0.048
    real(r8) :: fract_urea, fract_no3, soilph_min, soilph_max, max_runoff
    integer, parameter :: ind_region = 1
    integer :: def_ph_count, bad_runoff_count

    Hconc_grz(1:2) = (/10**(-8.5_r8), 10**(-8.0_r8)/)
    Hconc_slr(1:3) = (/10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8), 10.0_r8**(-8.0_r8)/)
    soilph_min = 999
    soilph_max = -999
    def_ph_count = 0
    do_balance_checks = mod(get_nstep(), balance_check_freq) == 0
    bad_runoff_count = 0
    max_runoff = -999
    
    associate(&
         ngrz => nitrogenflux_vars%man_n_grz_col, &
         man_u_grz => nitrogenstate_vars%man_u_grz_col, &
         man_a_grz => nitrogenstate_vars%man_a_grz_col, &
         man_r_grz => nitrogenstate_vars%man_r_grz_col, &
         man_u_app => nitrogenstate_vars%man_u_app_col, &
         man_a_app => nitrogenstate_vars%man_a_app_col, &
         man_r_app => nitrogenstate_vars%man_r_app_col, &
         ns => nitrogenstate_vars, &
         nf => nitrogenflux_vars, &
         cnv_nf => nitrogenflux_vars, &
         ram1 => frictionvel_vars%ram1_patch, &
         rb1 => frictionvel_vars%rb1_patch)

    nf%fert_n_appl_col(bounds%begc:bounds%endc) = 0.0
    nf%man_n_appl_col(bounds%begc:bounds%endc) = 0.0
    nf%man_tan_appl_col(bounds%begc:bounds%endc) = 0.0

    call p2c(bounds, num_soilc, filter_soilc, &
         cnv_nf%fert_patch(bounds%begp:bounds%endp), &
         nf%fert_n_appl_col(bounds%begc:bounds%endc))
    !call p2c(bounds, num_soilc, filter_soilc, &
    !     cnv_nf%manu_patch(bounds%begp:bounds%endp), &
    !     nf%man_n_appl_col(bounds%begc:bounds%endc))

    nf%man_n_appl_col = 0.0

    if (any(nf%man_n_appl_col > 100)) then
       write(iulog, *) maxval(nf%man_n_appl_col)
       call endrun('bad man_n_appl_col')
    end if
    if (do_balance_checks) then
       nstored_old = get_total_n(ns, nf, 'pools_storage')
       nsoilman_old = get_total_n(ns, nf, 'pools_manure')
       nsoilfert_old = get_total_n(ns, nf, 'pools_fertilizer')
    end if

    ! Assign the "pastoral" manure entire to the natural vegetation column
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       l = col%landunit(c)
       if (.not. (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop)) cycle
       if (.not. col%active(c) .or. col%wtgcell(c) < 1e-6) cycle
       g = col%gridcell(c)
       if (lun%itype(l) == istsoil) then
          ngrz(c) = atm2lnd_vars%forc_ndep3_grc(g) / col%wtgcell(c) * 1e3 ! kg to g 
          if (debug_fan) then
             if (ngrz(c) > 1e12 .or. (isnan(ngrz(c)))) then
                write(iulog, *) 'bad ngrz', atm2lnd_vars%forc_ndep3_grc(g), col%wtgcell(c)
                call endrun('bad ngrz 1')
             end if
          end if
          if (nf%man_n_appl_col(c) > 0) then
             write(iulog, *) nf%man_n_appl_col(c)
             call endrun(msg='Found fertilizer in soil column')
          end if
       else
          ngrz(c) = 0.0
       end if

    end do
    !ngrz = 0
    
    if(debug_fan) then
       write(iulog, *) 'nan count of storage 1', count(isnan(ns%man_n_stored_col))
       if (any(isnan(nf%man_n_appl_col))) then
          call endrun('nan nh3 appl b')
       end if
    end if

    call handle_storage_v2(bounds, temperature_vars, frictionvel_vars, dt, &
         atm2lnd_vars%forc_ndep2_grc, &
         ns%man_n_stored_col, ns%man_tan_stored_col, &
         nf%man_n_appl_col, nf%man_tan_appl_col, &
         nf%man_n_grz_col, nf%man_n_mix_col, &
         nf%nh3_stores_col, nf%nh3_barns_col, &
         nf%man_n_transf_col, ns%fan_grz_fract_col, &
         nf%man_n_barns_col, &
         fract_tan, &
         filter_soilc, num_soilc)
    
    if (debug_fan) then
       if (any(isnan(nf%nh3_stores_col))) then
          call endrun('nan nh3 stores')
       end if
       if (any(isnan(nf%nh3_barns_col))) then
          call endrun('nan nh3 barns')
       end if
       if (any(isnan(nf%man_n_appl_col))) then
          call endrun('nan nh3 appl')
       end if
       if (any(isnan(nf%man_n_mix_col))) then
          call endrun('nan nh3 appl')
       end if
    end if

    do fc = 1, num_soilc
       c = filter_soilc(fc)
       l = col%landunit(c)
       g = col%gridcell(c)
       if (.not. (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop)) cycle
       if (.not. col%active(c) .or. col%wtgcell(c) < 1e-15) cycle

       if (nf%man_n_appl_col(c) > 1e12 .or. ngrz(c) > 1e12) then
          write(iulog, *) c, nf%man_n_appl_col(c), ngrz(c), cnv_nf%fert_patch(col%pfti(c):col%pftf(c))
          call endrun('nf%man_n_appl_col(c) is spval')
       end if

       ! Find and average the atmospheric resistances Rb and Ra.
       ! 
       if (lun%itype(col%landunit(c)) == istcrop) then
          ! for column, only one patch
          p = col%pfti(c)
          if (p /= col%pftf(c)) call endrun(msg='Strange patch for crop')
          ratm = ram1(p) + rb1(p)
       else
          ! if natural, find average over grasses
          ratm = 0.0
          patchcounter = 0
          do p = col%pfti(c), col%pftf(c)
             if (patch%itype(p) == nc4_grass .or. patch%itype(p) == nc3_nonarctic_grass) then
                if (.not. patch%active(p) .or. ram1(p) == spval .or. rb1(p) == spval) cycle
                ratm = ratm + ram1(p) + rb1(p)
                patchcounter = patchcounter + 1
             end if
          end do
          if (patchcounter > 0) then
             ratm = ratm / patchcounter
          else
             ! grass not found, take something.
             do p = col%pfti(c), col%pftf(c)
                if (.not. patch%active(p) .or. ram1(p) == spval .or. rb1(p) == spval) cycle
                ratm = ram1(p) + rb1(p)
                exit
             end do
             if (p == col%pftf(c) + 1) then
                ratm = 150.0_r8
             end if
          end if
          ns%fan_grz_fract_col(c) = 1.0_r8 ! for crops handled by handle_storage
       end if

       watertend = 0.0_r8 
       tg = temperature_vars%t_grnd_col(c)
       theta = waterstate_vars%h2osoi_vol_col(c,1)
       thetasat = soilstate_vars%watsat_col(c,1)
       theta = min(theta, 0.98_r8*thetasat)
       theta = max(theta, 1e-6) !?
       infiltr_m_s = max(waterflux_vars%qflx_infl_col(c), 0.0) * 1e-3 
       evap_m_s = waterflux_vars%qflx_evap_grnd_col(c) * 1e-3
       runoff_m_s = max(waterflux_vars%qflx_surf_col(c), 0.0) * 1e-3

       if (runoff_m_s > 1e-3) then
          bad_runoff_count = bad_runoff_count + 1
          runoff_m_s = 0.0
       end if
       if (runoff_m_s > max_runoff) max_runoff = runoff_m_s
       !
       ! grazing
       !

       ndep_org(ind_avail) = ngrz(c) * (1.0_r8-fract_tan) * fract_avail
       ndep_org(ind_resist) = ngrz(c) * (1.0_r8-fract_tan) * fract_resist
       ndep_org(ind_unavail) = ngrz(c) * (1.0_r8-fract_tan) * fract_unavail
       tandep = ngrz(c) * fract_tan

       orgpools(ind_avail) = man_a_grz(c)
       orgpools(ind_resist) = man_r_grz(c)
       orgpools(ind_unavail) = man_u_grz(c)
       call update_org_n(ndep_org, tg, orgpools, dt, tanprod, soilflux_org)
       man_a_grz(c) = orgpools(ind_avail)
       man_r_grz(c) = orgpools(ind_resist) 
       man_u_grz(c) = orgpools(ind_unavail)

       tanpools3(1) = ns%tan_g1_col(c)
       tanpools3(2) = ns%tan_g2_col(c)
       tanpools3(3) = ns%tan_g3_col(c)
       if (any(isnan(tanpools3))) then
          call endrun('nan1')
       end if

       ph_soil = atm2lnd_vars%forc_soilph_grc(g)
       if (ph_soil < 3.0) then
          ph_soil = 6.5_r8
          def_ph_count = def_ph_count + 1
       end if
       Hconc_grz(3) = 10**-(ph_soil)
       soilph_max = max(soilph_max, ph_soil)
       soilph_min = min(soilph_min, ph_soil)

       fluxes_tmp = 0.0
       garbage_total = 0.0
       fluxes3 = 0.0
       garbage = 0
       do ind_substep = 1, num_substeps
          call update_npool(tg, ratm, &
               theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_vars%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, tandep, (/0.0_r8, 0.0_r8, sum(tanprod)/), water_init_grz, &
               cnc_nh3_air, poolranges_grz, Hconc_grz, dz_layer_grz, tanpools3, &
               fluxes3(1:5,:), garbage, dt/num_substeps, status, 3)
          if (status /= 0) then
             write(iulog, *) 'status = ', status, tanpools2, ratm, theta, thetasat, tandep, tanprod
             call endrun(msg='update_npool status /= 0')
          end if
          if (debug_fan .and. any(isnan(tanpools2))) then
             call endrun('nan2')
          end if
          fluxes_tmp = fluxes_tmp + sum(fluxes3, dim=2)
          garbage_total = garbage_total + garbage
       end do
       fluxes_tmp = fluxes_tmp / num_substeps

       ns%tan_g1_col(c) = tanpools3(1)
       ns%tan_g2_col(c) = tanpools3(2)
       ns%tan_g3_col(c) = tanpools3(3)
       if (debug_fan .and. any(isnan(fluxes3))) then
          write(iulog, *) fluxes3
          call endrun('nan3')
       end if

       nf%nh3_grz_col(c) = fluxes_tmp(iflx_air)
       nf%manure_runoff_col(c) = fluxes_tmp(iflx_roff)
       nf%manure_no3_prod_col(c) = fluxes_tmp(iflx_no3)
       nf%manure_nh4_to_soil_col(c) &
            = fluxes_tmp(iflx_soild) + fluxes_tmp(iflx_soilq) + garbage_total / dt + soilflux_org

       !
       ! Manure application
       !

       org_n_tot = nf%man_n_appl_col(c) - nf%man_tan_appl_col(c)
       ! Use the the same fractionation of organic N as for grazing, after removing the
       ! "explicitly" calculated TAN.
       if (1-fract_tan > 1e-6) then
          ndep_org(ind_avail) = org_n_tot * fract_avail! / (1-fract_tan)
          ndep_org(ind_resist) = org_n_tot * fract_resist! / (1-fract_tan)
          ndep_org(ind_unavail) = org_n_tot * fract_unavail! / (1-fract_tan)
       else
          ndep_org = 0.0
       end if
       tandep = nf%man_tan_appl_col(c)

       orgpools(ind_avail) = man_a_app(c)
       orgpools(ind_resist) = man_r_app(c)
       orgpools(ind_unavail) = man_u_app(c)
       call update_org_n(ndep_org, tg, orgpools, dt, tanprod, soilflux_org)
       man_a_app(c) = orgpools(ind_avail)
       man_r_app(c) = orgpools(ind_resist)
       man_u_app(c) = orgpools(ind_unavail)
       tanpools4(1) = ns%tan_s0_col(c)
       tanpools4(2) = ns%tan_s1_col(c)
       tanpools4(3) = ns%tan_s2_col(c)
       tanpools4(4) = ns%tan_s3_col(c)

       ph_crop = min(max(ph_soil, 5.5_r8), 7.5_r8)
       Hconc_slr(4) = 10**-(ph_crop)

       if (debug_fan .and. any(isnan(tanpools4))) then
          call endrun('nan31')
       end if

       fluxes_tmp = 0.0
       garbage_total = 0.0
       fluxes4 = 0.0
       do ind_substep = 1, num_substeps
          if (debug_fan .and. any(abs(tanpools4) > 1e12)) then
             write(iulog, *) ind_substep, tanpools4, tandep, nf%fert_n_appl_col(c), &
                  nf%man_n_appl_col(c), ns%man_n_stored_col(c), ns%man_tan_stored_col(c)
             call endrun('bad tanpools (manure app)')
          end if

          call update_4pool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_vars%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, tandep, sum(tanprod), cnc_nh3_air, depth_slurry, &
               poolranges_slr, tanpools4, Hconc_slr, fluxes4(1:5,:), garbage, dt / num_substeps, status)
          if (status /= 0) then
             write(iulog, *) 'status = ', status, tanpools4, tg, ratm, 'th', theta, &
                  thetasat, tandep, 'tp', tanprod, 'fx', fluxes4(1:5,:), 'roff', runoff_m_s
             write(iulog, *) fluxes4(1:5,1)
             write(iulog, *) fluxes4(1:5,2)
             write(iulog, *) fluxes4(1:5,3)
             write(iulog, *) fluxes4(1:5,4)
             
             call endrun(msg='update_3pool status /= 0')
          end if
          if (debug_fan .and. any(isnan(fluxes4))) then
             write(iulog, *) fluxes4
             write(iulog, *), tanpools4
             write(iulog, * ) ratm, theta, thetasat, infiltr_m_s, tandep, tanprod
             call endrun('nan4')
          end if

          fluxes_tmp = fluxes_tmp + sum(fluxes4, dim=2)
          garbage_total = garbage_total + garbage
       end do
       fluxes_tmp = fluxes_tmp / num_substeps

       ns%tan_s0_col(c) = tanpools4(1)
       ns%tan_s1_col(c) = tanpools4(2)
       ns%tan_s2_col(c) = tanpools4(3)
       ns%tan_s3_col(c) = tanpools4(4)


       nf%nh3_man_app_col(c) = fluxes_tmp(iflx_air)
       nf%manure_runoff_col(c) = nf%manure_runoff_col(c) + fluxes_tmp(iflx_roff)
       nf%manure_no3_prod_col(c) = nf%manure_no3_prod_col(c) + fluxes_tmp(iflx_no3)
       nf%manure_nh4_to_soil_col(c) &
            = nf%manure_nh4_to_soil_col(c) + fluxes_tmp(iflx_soild) + fluxes_tmp(iflx_soilq) &
            + garbage_total / dt + soilflux_org

       !
       ! Fertilizer
       !

       fert_total = nf%fert_n_appl_col(c)
       fract_urea = atm2lnd_vars%forc_ndep_urea_grc(g)
       fract_no3 = atm2lnd_vars%forc_ndep_nitr_grc(g)

       if (fract_urea < 0 .or. fract_no3 < 0 .or. fract_urea + fract_no3 > 1) then
          call endrun('bad fertilizer fractions')
       end if

       fert_urea = fert_total * fract_urea
       fert_no3 = fert_total * fract_no3
       fert_generic = fert_total - fert_urea - fert_no3
       nf%otherfert_n_appl_col(c) = fert_no3 + fert_generic

       ! Urea decomposition 
       ! 
       ureapools(1) = ns%fert_u0_col(c)
       ureapools(2) = ns%fert_u1_col(c)
       fluxes2 = 0.0
       call update_urea(tg, theta, thetasat, infiltr_m_s, evap_m_s, watertend, &
            runoff_m_s, fert_urea, ureapools,  fluxes2, urea_resid, poolranges_fert(1:2), &
            dt, status, numpools=2)
       if (status /= 0) then
          call endrun(msg='Bad status after update_urea for fertilizer')
       end if
       ! Nitrogen fluxes from urea pool. Be sure to not zero below!
       fluxes_tmp = sum(fluxes2, dim=2)

       ns%fert_u0_col(c) = ureapools(1)
       ns%fert_u1_col(c) = ureapools(2)
       ! Collect the formed ammonia for updating the TAN pools
       tanprod_from_urea(1:2) = fluxes2(iflx_to_tan, 1:2)
       tanprod_from_urea(2) = tanprod_from_urea(2)
       ! There is no urea pool corresponding to tan_f2, because most of the urea will
       ! have decomposed. Here whatever remains gets sent to tan_f2. 
       tanprod_from_urea(3) = urea_resid / dt 

       tanpools3(1) = ns%tan_f0_col(c)
       tanpools3(2) = ns%tan_f1_col(c)
       tanpools3(3) = ns%tan_f2_col(c)         
       garbage_total = 0.0
       fluxes3 = 0.0
       nf%nh3_otherfert_col(c) = 0.0
       do ind_substep = 1, num_substeps
          ! Fertilizer pools f0...f2
          call update_npool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_vars%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, 0.0_r8, tanprod_from_urea, water_init_fert, cnc_nh3_air, &
               poolranges_fert, Hconc_fert, dz_layer_fert, tanpools3, fluxes3(1:5,:), &
               garbage, dt/num_substeps, status, numpools=3)
          if (status /= 0) then
             write(iulog, *) 'status:', status, tanpools3, nf%fert_n_appl_col(c)
             call endrun(msg='Bad status after npool for fertilizer')
          end if
          fluxes_tmp = fluxes_tmp + sum(fluxes3, dim=2) / num_substeps
          garbage_total = garbage_total + garbage

          ! Fertilizer pool f3
          call update_npool(tg, ratm, theta, thetasat, infiltr_m_s, evap_m_s, &
               atm2lnd_vars%forc_q_downscaled_col(c), watertend, &
               runoff_m_s, fert_generic, (/0.0_r8/), water_init_fert, cnc_nh3_air, &
                                !(/360*24*3600.0_r8/), (/10**(-6.0_r8)/), dz_layer_fert, ns%tan_f3_col(c:c), fluxes3(1:5,1:1), &
               (/360*24*3600.0_r8/), (/10**(-ph_crop)/), dz_layer_fert, ns%tan_f3_col(c:c), fluxes3(1:5,1:1), &
               garbage, dt/num_substeps, status, numpools=1)
          if (status /= 0) then
             write(iulog, *) 'status:', status, tanpools3, nf%fert_n_appl_col(c)
             call endrun(msg='Bad status after npool for generic')
          end if
          fluxes_tmp = fluxes_tmp + fluxes3(:, 1) / num_substeps
          garbage_total = garbage_total + garbage
          nf%nh3_otherfert_col(c) = nf%nh3_otherfert_col(c) + fluxes3(iflx_air, 1) / num_substeps
       end do

       ns%tan_f0_col(c) = tanpools3(1)
       ns%tan_f1_col(c) = tanpools3(2)
       ns%tan_f2_col(c) = tanpools3(3)
       ! !!tan_f3_col already updated above by update_npool!!

       nf%nh3_fert_col(c) = fluxes_tmp(iflx_air)
       nf%fert_runoff_col(c) = fluxes_tmp(iflx_roff)
       nf%fert_no3_prod_col(c) = fluxes_tmp(iflx_no3) + fert_no3
       nf%fert_nh4_to_soil_col(c) = fluxes_tmp(iflx_soild) + fluxes_tmp(iflx_soilq) + garbage_total/dt 

       ! Total flux
       ! 
       nf%nh3_total_col(c) = nf%nh3_fert_col(c) + nf%nh3_man_app_col(c) &
            + nf%nh3_grz_col(c) + nf%nh3_stores_col(c) +  nf%nh3_barns_col(c)
       if (nf%nh3_total_col(c) < -1e15) then
          call endrun(msg='ERROR: FAN, negative total emission')
       end if
    end do

    if (do_balance_checks) then
       call balance_check('Storage', nstored_old, &
            get_total_n(ns, nf, 'pools_storage'), get_total_n(ns, nf, 'fluxes_storage'))
       call balance_check('Manure', nsoilman_old, &
            get_total_n(ns, nf, 'pools_manure'), get_total_n(ns, nf, 'fluxes_manure'))
       call balance_check('Fertilizer', nsoilfert_old, &
            get_total_n(ns, nf, 'pools_fertilizer'), get_total_n(ns, nf, 'fluxes_fertilizer'))
       write(iulog, *) 'SoilPH check:', soilph_min, soilph_max, def_ph_count
       write(iulog, *) 'Runoff check:', bad_runoff_count, max_runoff
    end if

  end associate

end subroutine fan

real(r8) function get_total_n(ns, nf, which) result(total)
  type(nitrogenstate_type), intent(in) :: ns
  type(nitrogenflux_type), intent(in) :: nf
  character(len=*), intent(in) :: which

  total = 0

  associate(soilc => filter_soilc(1:num_soilc))

  select case(which)
  case('pools_storage')
     total = sum(ns%man_n_stored_col(soilc))

  case('fluxes_storage')
     total = sum(nf%man_n_mix_col(soilc))
     total = total - sum(nf%nh3_stores_col(soilc))
     total = total - sum(nf%nh3_barns_col(soilc)) - sum(nf%man_n_transf_col(soilc))

  case('pools_manure')
     total = total + sum(ns%tan_g1_col(soilc)) + sum(ns%tan_g2_col(soilc)) + sum(ns%tan_g3_col(soilc)) 
     total = total + sum(ns%man_u_grz_col(soilc)) &
          + sum(ns%man_a_grz_col(soilc)) + sum(ns%man_r_grz_col(soilc))
     total = total + sum(ns%tan_s0_col(soilc)) &
          + sum(ns%tan_s1_col(soilc)) + sum(ns%tan_s2_col(soilc)) + sum(ns%tan_s3_col(soilc))
     total = total + sum(ns%man_u_app_col(soilc)) &
          + sum(ns%man_a_app_col(soilc)) + sum(ns%man_r_app_col(soilc))

  case('fluxes_manure')
     total = sum(nf%man_n_grz_col(soilc)) + sum(nf%man_n_appl_col(soilc)) 
     total = total - sum(nf%nh3_man_app_col(soilc)) &
          - sum(nf%nh3_grz_col(soilc)) - sum(nf%manure_runoff_col(soilc))
     total = total - sum(nf%manure_no3_prod_col(soilc)) - sum(nf%manure_nh4_to_soil_col(soilc))

  case('pools_fertilizer')
     total = sum(ns%tan_f0_col((soilc))) + sum(ns%tan_f1_col((soilc))) + sum(ns%tan_f2_col(soilc)) &
          + sum(ns%tan_f3_col(soilc))
     total = total + sum(ns%fert_u0_col(soilc)) + sum(ns%fert_u1_col(soilc))

  case('fluxes_fertilizer')
     total = sum(nf%fert_n_appl_col(soilc))
     total = total - sum(nf%nh3_fert_col(soilc)) - sum(nf%fert_runoff_col(soilc))
     total = total - sum(nf%fert_no3_prod_col(soilc)) - sum(nf%fert_nh4_to_soil_col(soilc))

  case default
     call endrun(msg='Bad argument to get_total_n')

  end select

end associate

end function get_total_n

subroutine balance_check(label, total_old, total_new, flux)
  ! Check and report that the net flux equals the accumulated mass in pools. The
  ! total pools and fluxes can be evaluated by the function get_total_n.
character(len=*), intent(in) :: label
real(r8), intent(in) :: total_old, total_new, flux

real(r8) :: diff, accflux
real(r8) :: tol = 1e-6_r8

diff = total_new - total_old
accflux = flux*dt
write(iulog, *) 'Balance check:', label, diff, accflux

end subroutine balance_check


end subroutine NitrogenDeposition

  subroutine handle_storage_v2(bounds, temperature_inst, frictionvel_inst, dt,  &
       ndep_mixed_grc, n_stored_col, tan_stored_col, &
       n_manure_spread_col, tan_manure_spread_col, &
       n_manure_graze_col, n_manure_mixed_col, &
       nh3_flux_stores, nh3_flux_barns, man_n_transf, &
       grz_fract, man_n_barns, tan_fract_excr, &
       filter_soilc, num_soilc)
    use landunit_varcon, only : max_lunit
    use pftvarcon, only : nc4_grass, nc3_nonarctic_grass
    use clm_varcon, only : ispval
    use landunit_varcon,      only:  istsoil, istcrop
    use abortutils     , only : endrun
    use LandunitType   , only: lun => lun_pp
    use GridcellType   , only: grc => grc_pp
    use clm_varctl     , only : iulog
    use ColumnType     , only : col => col_pp
    use VegetationType       , only : patch => veg_pp
    use FanMod
    
    implicit none
    type(bounds_type), intent(in)    :: bounds
    type(temperature_type) , intent(in) :: temperature_inst
    type(frictionvel_type) , intent(in) :: frictionvel_inst
    real(r8), intent(in) :: dt
    
    ! N excreted in manure, mixed/pastoral systems, gN/m2:
    real(r8), intent(in) :: ndep_mixed_grc(bounds%begg:bounds%endg)
    real(r8), intent(inout) :: n_stored_col(bounds%begc:bounds%endc), tan_stored_col(bounds%begc:bounds%endc) ! N, TAN currently stored, gN/m2
    ! N, TAN spread on grasslands, gN/m2/s:
    real(r8), intent(inout) :: n_manure_spread_col(bounds%begc:bounds%endc) ! for crops, input, determined by crop model, otherwise output
    real(r8), intent(out) :: tan_manure_spread_col(bounds%begc:bounds%endc) ! output, calculated from the above and stored manure
    ! N excreted by animals allocated to mixed production systems temporarily grazing on grasslands:
    real(r8), intent(inout) :: n_manure_graze_col(bounds%begc:bounds%endc)
    ! N excreted by animals in mixed systems, total
    real(r8), intent(out) :: n_manure_mixed_col(bounds%begc:bounds%endc)
    ! NH3 emission fluxes from manure storage and housings, gN/m2/s
    real(r8), intent(out) :: nh3_flux_stores(bounds%begc:bounds%endc), nh3_flux_barns(bounds%begc:bounds%endc)
    ! total nitrogen flux transferred out of a crop column
    real(r8), intent(out) :: man_n_transf(bounds%begc:bounds%endc)
    real(r8), intent(out) :: man_n_barns(bounds%begc:bounds%endc)
    ! fraction of manure excreted when grazing
    real(r8), intent(out) :: grz_fract(bounds%begc:bounds%endc)
    ! TAN fraction in excreted N
    real(r8), intent(in) :: tan_fract_excr
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns

    integer :: begg, endg, g, l, c, il, counter, col_grass, status, p
    real(r8) :: flux_avail, flux_grazing
    real(r8) :: tempr_ave, windspeed_ave ! windspeed and temperature averaged over agricultural patches
    real(r8) :: tempr_barns, tempr_stores, vent_barns, flux_grass_crop, tempr_min_10day, &
         flux_grass_graze, flux_grass_spread, flux_grass_spread_tan, flux_grass_crop_tan
    real(r8) :: cumflux, totalinput
    real(r8) :: fluxes_nitr(4), fluxes_tan(4)
    ! The fraction of manure applied continuously on grasslands (if present in the gridcell)
    real(r8), parameter :: fract_continuous = 0.1_r8, kg_to_g = 1e3_r8, max_grazing_fract = 0.5_r8, &
         volat_coef_barns = 0.03_r8, volat_coef_stores = 0.025_r8, &
         tempr_min_grazing = 283.0_r8!!!!

    begg = bounds%begg; endg = bounds%endg
    nh3_flux_stores(bounds%begc:bounds%endc) = 0_r8
    nh3_flux_barns(bounds%begc:bounds%endc) = 0_r8
    man_n_barns(bounds%begc:bounds%endc) = 0.0_r8

    totalinput = 0.0
    cumflux = 0.0
    
    do g = begg, endg
       !totalinput = totalinput + ndep_mixed_grc(g)
       
       ! First find out if there are grasslands in this cell. If yes, a fraction of
       ! manure can be diverted to them before storage.
       col_grass = ispval
       do il = 1, max_lunit
          l = grc%landunit_indices(il, g)
          if (lun%itype(l) == istsoil) then
             do p = lun%pfti(l), lun%pftf(l)
                if (patch%itype(p) == nc4_grass .or. patch%itype(p) == nc3_nonarctic_grass) then
                   col_grass = patch%column(p)
                   exit
                end if
             end do
          end if
          if (col_grass /= ispval) exit
       end do
       if (col%wtgcell(col_grass) < 1e-6) col_grass = ispval
       ! Transfer of manure from all crop columns to the natural vegetation column:
       flux_grass_graze = 0_r8
       flux_grass_spread = 0_r8
       flux_grass_spread_tan = 0_r8

       do il = 1, max_lunit
          l = grc%landunit_indices(il, g)
          if (l == ispval) cycle
          if (lun%itype(l) == istcrop) then
             ! flux_avail = manure excreted per m2 of crops (ndep_mixed_grc = per m2 / all land units)
             do c = lun%coli(l), lun%colf(l)
                if (.not. col%active(c)) cycle
                if (col%wtgcell(c) < 1e-6) cycle

                if (col%landunit(c) /= l) then
                   write(iulog, *) g, il, c, col%landunit(c)
                   call endrun('something wrong')
                end if
                if (.not. any(c==filter_soilc(1:num_soilc))) then
                   write(iulog, *) c, n_manure_spread_col(c)
                   call endrun('column not in soilfilter')
                end if

                flux_avail = ndep_mixed_grc(g) * kg_to_g / lun%wtgcell(l)
                if (flux_avail > 1e12 .or. isnan(flux_avail)) then
                   write(iulog, *) 'bad flux_avail', ndep_mixed_grc(g), lun%wtgcell(l)
                   call endrun('bad flux_avail')
                end if
                n_manure_mixed_col(c) = flux_avail
                totalinput = totalinput + flux_avail

                counter = 0
                if (col_grass == c) call endrun('Something wrong with the indices')
                if (col%pfti(c) /= col%pftf(c)) then
                   call endrun(msg="ERROR crop column has multiple patches")
                end if

                tempr_ave = temperature_inst%t_ref2m_patch(col%pfti(c))
                windspeed_ave = frictionvel_inst%u10_patch(col%pfti(c))

                tempr_min_10day = temperature_inst%t_a10min_patch(col%pfti(c))
                if (tempr_min_10day > tempr_min_grazing) then
                   ! fraction of animals grazing -> allocate some manure to grasslands before barns
                   flux_grazing = max_grazing_fract * flux_avail
                   flux_avail = flux_avail - flux_grazing
                   grz_fract(c) = max_grazing_fract
                else
                   flux_grazing = 0
                   grz_fract(c) = 0
                end if
                flux_grass_graze = flux_grass_graze + flux_grazing*col%wtgcell(c)

                man_n_barns(c) = flux_avail
                
                call eval_fluxes_storage(flux_avail, tempr_ave, windspeed_ave, 0.0_r8, &
                     volat_coef_barns, volat_coef_stores, tan_fract_excr, fluxes_nitr, fluxes_tan, status)
                if (any(fluxes_nitr > 1e12)) then
                   write(iulog, *) 'bad fluxes', fluxes_nitr
                end if
                if (status /=0) then 
                   write(iulog, *) 'status = ', status
                   call endrun(msg='eval_fluxes_storage failed')
                end if
                cumflux = cumflux + sum(fluxes_nitr)
                
                if (fluxes_tan(iflx_to_store) < 0) then
                   call endrun(msg="ERROR too much manure lost")
                end if

                flux_grass_spread = flux_grass_spread + fluxes_nitr(iflx_to_store)*col%wtgcell(c)
                flux_grass_spread_tan = flux_grass_spread_tan + fluxes_tan(iflx_to_store)*col%wtgcell(c)

                man_n_transf(c) = flux_grazing + fluxes_nitr(iflx_to_store)
                
                nh3_flux_stores(c) = fluxes_nitr(iflx_air_stores)
                nh3_flux_barns(c) = fluxes_nitr(iflx_air_barns)
                
             end do ! column
          end if ! crop land unit
       end do ! landunit

       if (col_grass /= ispval) then
          if (tan_manure_spread_col(col_grass) > 1) then
             write(iulog, *) 'bad tan_manure col_grass before adding', n_manure_spread_col(col_grass), &
                  tan_manure_spread_col(col_grass)
          end if
          n_manure_spread_col(col_grass) = n_manure_spread_col(col_grass) &
               + flux_grass_spread / col%wtgcell(col_grass)
          tan_manure_spread_col(col_grass) = tan_manure_spread_col(col_grass) &
               + flux_grass_spread_tan / col%wtgcell(col_grass)
          n_manure_graze_col(col_grass) = n_manure_graze_col(col_grass) + flux_grass_graze / col%wtgcell(col_grass)
          !write(iulog, *) 'to grass:', n_manure_spread(col_grass), col_grass
          if (tan_manure_spread_col(col_grass) > 1) then
             write(iulog, *) 'bad tan_manure col_grass', flux_grass_spread_tan, col%wtgcell(col_grass)
          end if
       else if (flux_grass_spread > 0) then
          call endrun('Cannot spread manure')
       end if

    end do ! grid

  end subroutine handle_storage_v2


  !-----------------------------------------------------------------------
  subroutine NitrogenFixation(num_soilc, filter_soilc, waterflux_vars, &
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
    type(waterflux_type)     , intent(in)    :: waterflux_vars    
    type(carbonflux_type)   , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars 
    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc                  ! indices
    real(r8) :: t                     ! temporary
    real(r8) :: dayspyr               ! days per year
    real(r8) :: secspyr              ! seconds per yr
    logical  :: do_et_bnf = .false.
    !-----------------------------------------------------------------------

    associate(& 
         cannsum_npp    => col_cf%annsum_npp      , & ! Input:  [real(r8) (:)]  nitrogen deposition rate (gN/m2/s)                
         col_lag_npp    => col_cf%lag_npp         , & ! Input: [real(r8) (:)]  (gC/m2/s) lagged net primary production           

         qflx_tran_veg  => col_wf%qflx_tran_veg    , & ! col vegetation transpiration (mm H2O/s) (+ = to atm)
         
         qflx_evap_veg  => col_wf%qflx_evap_veg    , & ! col vegetation evaporation (mm H2O/s) (+ = to atm)
         nfix_to_sminn  => col_nf%nfix_to_sminn   & ! Output: [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
         )

      dayspyr = get_days_per_year()

      if (do_et_bnf) then
         secspyr = dayspyr * 86400._r8
         do fc = 1, num_soilc
            c =filter_soilc(fc)
            !use the cleveland equation
            t = 0.00102_r8*(qflx_evap_veg(c)+qflx_tran_veg(c))+0.0524_r8/secspyr
            nfix_to_sminn(c) = max(0._r8, t)
         enddo
      else
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
      endif

    end associate

  end subroutine NitrogenFixation
 
  !-----------------------------------------------------------------------
  subroutine NitrogenLeaching(bounds, num_soilc, filter_soilc, &
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
         h2osoi_liq          => col_ws%h2osoi_liq            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)

         qflx_drain          => col_wf%qflx_drain             , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf           => col_wf%qflx_surf              , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)                        
         
         sminn_vr            => col_ns%sminn_vr           , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N                          
         smin_no3_vr         => col_ns%smin_no3_vr        , & ! Input:  [real(r8) (:,:) ]                                                  
         sminn_leached_vr    => col_nf%sminn_leached_vr    , & ! Output: [real(r8) (:,:) ]  rate of mineral N leaching (gN/m3/s)            
         smin_no3_leached_vr => col_nf%smin_no3_leached_vr , & ! Output: [real(r8) (:,:) ]  rate of mineral NO3 leaching (gN/m3/s)          
         smin_no3_runoff_vr  => col_nf%smin_no3_runoff_vr    & ! Output: [real(r8) (:,:) ]  rate of mineral NO3 loss with runoff (gN/m3/s)  
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
               surface_water(c) = surface_water(c) + h2osoi_liq(c,j) * ( (depth_runoff_Nloss - zisoi(j-1)) / col_pp%dz(c,j))
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
                     disn_conc = (sf * sminn_vr(c,j) * col_pp%dz(c,j) )/(h2osoi_liq(c,j) )
                  end if

                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  sminn_leached_vr(c,j) = disn_conc * drain_tot(c) * h2osoi_liq(c,j) / ( tot_water(c) * col_pp%dz(c,j) )

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
                     disn_conc = (sf_no3 * smin_no3_vr(c,j) * col_pp%dz(c,j) )/(h2osoi_liq(c,j) )
                  end if
                  !
                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  smin_no3_leached_vr(c,j) = disn_conc * drain_tot(c) * h2osoi_liq(c,j) / ( tot_water(c) * col_pp%dz(c,j) )
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
                          h2osoi_liq(c,j) / ( surface_water(c) * col_pp%dz(c,j) )
                  elseif ( zisoi(j-1) < depth_runoff_Nloss )  then
                     smin_no3_runoff_vr(c,j) = disn_conc * qflx_surf(c) * &
                          h2osoi_liq(c,j) * ((depth_runoff_Nloss - zisoi(j-1)) / &
                          col_pp%dz(c,j)) / ( surface_water(c) * (depth_runoff_Nloss-zisoi(j-1) ))
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
  end subroutine NitrogenLeaching

  !-----------------------------------------------------------------------
  subroutine NitrogenFert(bounds, num_soilc, filter_soilc, &
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
         fert          =>    veg_nf%fert          , & ! Input:  [real(r8) (:)]  nitrogen fertilizer rate (gN/m2/s)                
         fert_to_sminn =>    col_nf%fert_to_sminn   & ! Output: [real(r8) (:)]                                                    
         )
      
      call p2c(bounds, num_soilc, filter_soilc, &
           fert(bounds%begp:bounds%endp), &
           fert_to_sminn(bounds%begc:bounds%endc))

    end associate
  end subroutine NitrogenFert

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
         wf               =>  col_ws%wf                 , & ! Input:  [real(r8) (:) ]  soil water as frac. of whc for top 0.5 m          

         hui              =>  crop_vars%gddplant_patch               , & ! Input:  [real(r8) (:) ]  gdd since planting (gddplant)                    

         fpg              =>  cnstate_vars%fpg_col                   , & ! Input:  [real(r8) (:) ]  fraction of potential gpp (no units)              
         gddmaturity      =>  cnstate_vars%gddmaturity_patch         , & ! Input:  [real(r8) (:) ]  gdd needed to harvest                             
         croplive         =>  crop_vars%croplive_patch            , & ! Input:  [logical  (:) ]  true if planted and not harvested                  

         sminn            =>  col_ns%sminn           , & ! Input:  [real(r8) (:) ]  (kgN/m2) soil mineral N                           
         plant_ndemand    =>  veg_nf%plant_ndemand  , & ! Input:  [real(r8) (:) ]  N flux required to support initial GPP (gN/m2/s)  
         
         soyfixn          =>  veg_nf%soyfixn        , & ! Output: [real(r8) (:) ]  nitrogen fixed to each soybean crop               
         soyfixn_to_sminn =>  col_nf%soyfixn_to_sminn   & ! Output: [real(r8) (:) ]                                                    
         )

      sminnthreshold1 = 30._r8
      sminnthreshold2 = 10._r8
      GDDfracthreshold1 = 0.15_r8
      GDDfracthreshold2 = 0.30_r8
      GDDfracthreshold3 = 0.55_r8
      GDDfracthreshold4 = 0.75_r8

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = veg_pp%column(p)

         ! if soybean currently growing then calculate fixation

         if (veg_pp%itype(p) == nsoybean .and. croplive(p)) then

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

  !-----------------------------------------------------------------------
  subroutine NitrogenFixation_balance(num_soilc, filter_soilc,cnstate_vars, carbonflux_vars, &
             nitrogenstate_vars, nitrogenflux_vars, temperature_vars, waterstate_vars, carbonstate_vars, phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! created, Aug 2015 by Q. Zhu
    ! On the radiation time step, update the nitrogen fixation rate
    ! as a function of (1) root NP status, (2) fraction of root that is nodulated, (3) carbon cost of root nitrogen uptake
    ! N2 fixation is based on Fisher 2010 GBC doi:10.1029/2009GB003621; Wang 2007 GBC doi:10.1029/2006GB002797; and Grand 2012 ecosys model
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
    use pftvarcon        , only : noveg
        
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnstate_type)      , intent(inout) :: cnstate_vars
    type(carbonflux_type)   , intent(inout) :: carbonflux_vars
    type(nitrogenstate_type), intent(in)    :: nitrogenstate_vars
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars
    type(temperature_type)  , intent(inout) :: temperature_vars
    type(waterstate_type)   , intent(in)    :: waterstate_vars
    type(carbonstate_type)  , intent(inout) :: carbonstate_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc,p                     ! indices
    real(r8) :: r_fix                      ! carbon cost of N2 fixation, gC/gN
    real(r8) :: r_nup                      ! carbon cost of root N uptake, gC/gN
    real(r8) :: f_nodule                   ! empirical, fraction of root that is nodulated
    real(r8) :: N2_aq                      ! aqueous N2 bulk concentration gN/m3 soil
    real(r8) :: nfix_tmp
    !-----------------------------------------------------------------------

    associate(& 
         ivt                   => veg_pp%itype                            , & ! input:  [integer  (:) ]  pft vegetation type  
         cn_scalar             => cnstate_vars%cn_scalar               , &
         cp_scalar             => cnstate_vars%cp_scalar               , &
         vmax_nfix             => veg_vp%vmax_nfix                 , &
         km_nfix               => veg_vp%km_nfix                   , &
         frootc                => veg_cs%frootc        , &
         nfix_to_sminn         => col_nf%nfix_to_sminn  , & ! output: [real(r8) (:)]  symbiotic/asymbiotic n fixation to soil mineral n (gn/m2/s)
         nfix_to_plantn        => veg_nf%nfix_to_plantn , &
         nfix_to_ecosysn       => col_nf%nfix_to_ecosysn, &
         pnup_pfrootc          => veg_ns%pnup_pfrootc, &
         benefit_pgpp_pleafc   => veg_ns%benefit_pgpp_pleafc , &
         t_soi10cm_col         => col_es%t_soi10cm       , &
         h2osoi_vol            => col_ws%h2osoi_vol       , &
         t_scalar              => col_cf%t_scalar           &
         )

      do fc=1,num_soilc
          c = filter_soilc(fc)
          nfix_to_sminn(c) = 0.0_r8
          nfix_to_ecosysn(c) = 0._r8
          do p = col_pp%pfti(c), col_pp%pftf(c)
              if (veg_pp%active(p).and. (veg_pp%itype(p) .ne. noveg)) then
                  ! calculate c cost of n2 fixation: fisher 2010 gbc doi:10.1029/2009gb003621
                  r_fix = -6.25_r8*(exp(-3.62_r8 + 0.27_r8*(t_soi10cm_col(c)-273.15_r8)*(1.0_r8-0.5_r8&
                       *(t_soi10cm_col(c)-273.15_r8)/25.15_r8))-2.0_r8) 
                  ! calculate c cost of root n uptake: rastetter 2001, ecosystems, 4(4), 369-388.
                  r_nup = benefit_pgpp_pleafc(p) / max(pnup_pfrootc(p),1e-20_r8)
                  ! calculate fraction of root that is nodulated: wang 2007 gbc doi:10.1029/2006gb002797
                  f_nodule = 1 - min(1.0_r8,r_fix / max(r_nup, 1e-20_r8))
                  ! np limitation factor of n2 fixation (not considered now)
                  ! calculate aqueous N2 concentration and bulk aqueous N2 concentration
                  ! aqueous N2 concentration under pure nitrogen is 6.1e-4 mol/L/atm (based on Hery's law)
                  ! 78% atm * 6.1e-4 mol/L/atm * 28 g/mol * 1e3L/m3 * water content m3/m3 at 10 cm
                  N2_aq = 0.78_r8 * 6.1e-4_r8 *28._r8 *1.e3_r8 * h2osoi_vol(c,4)
                  ! calculate n2 fixation rate for each pft and add it to column total
                  nfix_tmp = vmax_nfix(veg_pp%itype(p)) * frootc(p) * cn_scalar(p) *f_nodule * t_scalar(c,1) * &
                             N2_aq/ (N2_aq + km_nfix(veg_pp%itype(p))) 
                  if (NFIX_PTASE_plant) then
                     nfix_to_sminn(c) = nfix_to_sminn(c) + nfix_tmp  * veg_pp%wtcol(p) * (1._r8-veg_vp%alpha_nfix(veg_pp%itype(p)))
                     nfix_to_plantn(p) = nfix_tmp * veg_vp%alpha_nfix(veg_pp%itype(p))
                     nfix_to_ecosysn(c) = nfix_to_ecosysn(c) + nfix_tmp  * veg_pp%wtcol(p)
                  else
                     nfix_to_sminn(c) = nfix_to_sminn(c) + nfix_tmp  * veg_pp%wtcol(p)
                     nfix_to_plantn(p) = 0.0_r8
                     nfix_to_ecosysn(c) = nfix_to_ecosysn(c) + nfix_tmp  * veg_pp%wtcol(p)
                  end if
              end if
          end do
      end do

    end associate

  end subroutine NitrogenFixation_balance
  
end module NitrogenDynamicsMod
