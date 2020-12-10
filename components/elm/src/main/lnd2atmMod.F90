module lnd2atmMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use abortutils           , only : endrun
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use shr_megan_mod        , only : shr_megan_mechcomps_n
  use elm_varpar           , only : numrad, ndst, nlevgrnd, nlevsno, nlevsoi !ndst = number of dust bins.
  use elm_varcon           , only : rair, grav, cpair, hfus, tfrz, spval
  use elm_varctl           , only : iulog, use_c13, use_cn, use_lch4, use_voc
  use tracer_varcon        , only : is_active_betr_bgc
  use seq_drydep_mod       , only : n_drydep, drydep_method, DD_XLND
  use decompMod            , only : bounds_type
  use subgridAveMod        , only : p2g, c2g 
  use lnd2atmType          , only : lnd2atm_type
  use atm2lndType          , only : atm2lnd_type
  use CH4Mod               , only : ch4_type
  use CNCarbonFluxType     , only : carbonflux_type
  use DUSTMod              , only : dust_type
  use DryDepVelocity       , only : drydepvel_type
  use VocEmissionMod       , only : vocemis_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityType , only : frictionvel_type
  use SolarAbsorbedType    , only : solarabs_type
  use SurfaceAlbedoType    , only : surfalb_type
  use TemperatureType      , only : temperature_type
  use WaterFluxType        , only : waterflux_type
  use WaterstateType       , only : waterstate_type
  use GridcellType         , only : grc_pp
  use GridcellDataType     , only : grc_ef, grc_ws, grc_wf
  use ColumnDataType       , only : col_ws, col_wf, col_cf, col_es  
  use VegetationDataType   , only : veg_es, veg_ef, veg_ws, veg_wf
  use SoilHydrologyType    , only : soilhydrology_type 
  
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: lnd2atm
  public :: lnd2atm_minimal
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine lnd2atm_minimal(bounds, &
      waterstate_vars, surfalb_vars, energyflux_vars, lnd2atm_vars)
    !
    ! !DESCRIPTION:
    ! Compute clm_l2a_vars component of gridcell derived type. This routine computes
    ! the bare minimum of components necessary to get the first step of a
    ! run started.
    !
    ! !USES:
    use elm_varcon, only : sb
    
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds  
    type(waterstate_type) , intent(in)    :: waterstate_vars
    type(surfalb_type)    , intent(in)    :: surfalb_vars
    type(energyflux_type) , intent(in)    :: energyflux_vars
    type(lnd2atm_type)    , intent(inout) :: lnd2atm_vars 
    !
    ! !LOCAL VARIABLES:
    integer :: g                                    ! index
    
    !------------------------------------------------------------------------

    call c2g(bounds, &
         col_ws%h2osno (bounds%begc:bounds%endc), &
         lnd2atm_vars%h2osno_grc    (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       lnd2atm_vars%h2osno_grc(g) = lnd2atm_vars%h2osno_grc(g)/1000._r8
    end do

    call c2g(bounds, nlevgrnd, &
         col_ws%h2osoi_vol (bounds%begc:bounds%endc, :), &
         lnd2atm_vars%h2osoi_vol_grc    (bounds%begg:bounds%endg, :), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         surfalb_vars%albd_patch (bounds%begp:bounds%endp, :), &
         lnd2atm_vars%albd_grc   (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         surfalb_vars%albi_patch (bounds%begp:bounds%endp, :), &
         lnd2atm_vars%albi_grc   (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         veg_ef%eflx_lwrad_out (bounds%begp:bounds%endp), &
         lnd2atm_vars%eflx_lwrad_out_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       lnd2atm_vars%t_rad_grc(g) = sqrt(sqrt(lnd2atm_vars%eflx_lwrad_out_grc(g)/sb))
    end do

  end subroutine lnd2atm_minimal

  !------------------------------------------------------------------------
  subroutine lnd2atm(bounds, &
       atm2lnd_vars, surfalb_vars, frictionvel_vars, &
       waterstate_vars, waterflux_vars, energyflux_vars, &
       solarabs_vars, carbonflux_vars, drydepvel_vars, &
       vocemis_vars, dust_vars, ch4_vars, soilhydrology_vars, lnd2atm_vars) 
    !
    ! !DESCRIPTION:
    ! Compute lnd2atm_vars component of gridcell derived type
    !
    ! !USES:
    use CH4varcon  , only : ch4offline
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)     :: bounds  
    type(atm2lnd_type)     , intent(in)     :: atm2lnd_vars
    type(surfalb_type)     , intent(in)     :: surfalb_vars
    type(frictionvel_type) , intent(in)     :: frictionvel_vars
    type(waterstate_type)  , intent(inout)  :: waterstate_vars
    type(waterflux_type)   , intent(in)     :: waterflux_vars
    type(energyflux_type)  , intent(in)     :: energyflux_vars
    type(solarabs_type)    , intent(in)     :: solarabs_vars
    type(carbonflux_type)  , intent(in)     :: carbonflux_vars
    type(drydepvel_type)   , intent(in)     :: drydepvel_vars
    type(vocemis_type)     , intent(in)     :: vocemis_vars
    type(dust_type)        , intent(in)     :: dust_vars
    type(ch4_type)         , intent(in)     :: ch4_vars
    type(soilhydrology_type), intent(in)    :: soilhydrology_vars
    type(lnd2atm_type)     , intent(inout)  :: lnd2atm_vars 
    !
    ! !LOCAL VARIABLES:
    integer :: g, lvl             ! index
    real(r8), parameter :: amC   = 12.0_r8 ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8 ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    !----------------------------------------------------
    ! lnd -> atm
    !----------------------------------------------------
    
    ! First, compute the "minimal" set of fields.
    call lnd2atm_minimal(bounds, &
         waterstate_vars, surfalb_vars, energyflux_vars, lnd2atm_vars)

    call p2g(bounds, &
         veg_es%t_ref2m (bounds%begp:bounds%endp), &
         lnd2atm_vars%t_ref2m_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         veg_ws%q_ref2m (bounds%begp:bounds%endp), &
         lnd2atm_vars%q_ref2m_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_vars%u10_elm_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%u_ref10m_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         veg_ef%taux (bounds%begp:bounds%endp), &
         lnd2atm_vars%taux_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         veg_ef%tauy (bounds%begp:bounds%endp), &
         lnd2atm_vars%tauy_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         veg_wf%qflx_evap_tot (bounds%begp:bounds%endp), &
         lnd2atm_vars%qflx_evap_tot_grc     (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         solarabs_vars%fsa_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%fsa_grc    (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_vars%fv_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%fv_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g(bounds, &
         frictionvel_vars%ram1_patch (bounds%begp:bounds%endp), &
         lnd2atm_vars%ram1_grc       (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

    call p2g( bounds, &
         veg_ef%eflx_sh_tot (bounds%begp:bounds%endp), &
         lnd2atm_vars%eflx_sh_tot_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity',c2l_scale_type='urbanf',l2g_scale_type='unity')
    do g = bounds%begg, bounds%endg
       lnd2atm_vars%eflx_sh_tot_grc(g) =  lnd2atm_vars%eflx_sh_tot_grc(g) - &
            grc_ef%eflx_dynbal(g) 
    enddo

    call p2g(bounds, &
         veg_ef%eflx_lh_tot (bounds%begp:bounds%endp), &
         lnd2atm_vars%eflx_lh_tot_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    if (use_cn) then
       call c2g(bounds, &
            col_cf%nee(bounds%begc:bounds%endc), &
            lnd2atm_vars%nee_grc   (bounds%begg:bounds%endg), &
            c2l_scale_type= 'unity', l2g_scale_type='unity')

       if (use_lch4) then
          if (.not. ch4offline) then
             ! Adjust flux of CO2 by the net conversion of mineralizing C to CH4
             do g = bounds%begg,bounds%endg
                ! nem is in g C/m2/sec
                lnd2atm_vars%nee_grc(g) = lnd2atm_vars%nee_grc(g) + lnd2atm_vars%nem_grc(g) 
             end do
          end if
       end if

       ! Convert from gC/m2/s to kgCO2/m2/s
       do g = bounds%begg,bounds%endg
          lnd2atm_vars%nee_grc(g) = lnd2atm_vars%nee_grc(g)*convertgC2kgCO2
       end do
    else
       do g = bounds%begg,bounds%endg
          lnd2atm_vars%nee_grc(g) = 0._r8
       end do
    end if

    ! drydepvel
    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
       call p2g(bounds, n_drydep, &
            drydepvel_vars%velocity_patch (bounds%begp:bounds%endp, :), &
            lnd2atm_vars%ddvel_grc        (bounds%begg:bounds%endg, :), &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
    endif

    ! voc emission flux
    if (use_voc .and. shr_megan_mechcomps_n>0) then
       call p2g(bounds, shr_megan_mechcomps_n, &
            vocemis_vars%vocflx_patch(bounds%begp:bounds%endp,:), &
            lnd2atm_vars%flxvoc_grc  (bounds%begg:bounds%endg,:), &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
    end if

    ! dust emission flux
    call p2g(bounds, ndst, &
         dust_vars%flx_mss_vrt_dst_patch(bounds%begp:bounds%endp, :), &
         lnd2atm_vars%flxdst_grc        (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')


    ! ch4 flux
    if (use_lch4 .and. (.not. is_active_betr_bgc)) then
       call c2g( bounds,     &
            ch4_vars%ch4_surf_flux_tot_col (bounds%begc:bounds%endc), &
            lnd2atm_vars%flux_ch4_grc      (bounds%begg:bounds%endg), &
            c2l_scale_type= 'unity', l2g_scale_type='unity' )
    end if

    !----------------------------------------------------
    ! lnd -> rof
    !----------------------------------------------------

    call c2g( bounds, &
         col_wf%qflx_runoff (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofliq_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       lnd2atm_vars%qflx_rofliq_grc(g) = lnd2atm_vars%qflx_rofliq_grc(g) - grc_wf%qflx_liq_dynbal(g)
    enddo

    call c2g( bounds, &
         col_wf%qflx_surf (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofliq_qsur_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         col_wf%qflx_h2osfc_surf (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofliq_qsurp_grc  (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
		 
    call c2g( bounds, &
         col_wf%qflx_irr_demand (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_irr_demand_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         col_wf%qflx_drain (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofliq_qsub_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         col_wf%qflx_drain_perched (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofliq_qsubp_grc  (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         col_wf%qflx_qrgwl (bounds%begc:bounds%endc), &
         lnd2atm_vars%qflx_rofliq_qgwl_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )

    call c2g( bounds, &
         col_wf%qflx_snwcp_ice(bounds%begc:bounds%endc),  &
         lnd2atm_vars%qflx_rofice_grc     (bounds%begg:bounds%endg),  & 
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
	call c2g(bounds,  &
         surfalb_vars%coszen_col (bounds%begc:bounds%endc), &
         lnd2atm_vars%coszen_str (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')
		 
	do g = bounds%begg, bounds%endg
       lnd2atm_vars%coszen_str(g) = lnd2atm_vars%coszen_str(g)        
    enddo
	
    do g = bounds%begg, bounds%endg
       lnd2atm_vars%qflx_rofice_grc(g) = lnd2atm_vars%qflx_rofice_grc(g) - grc_wf%qflx_ice_dynbal(g) 
    enddo

    ! calculate total water storage for history files
    ! first set tws to gridcell total endwb
    ! second add river storage as gridcell average depth (1.e-3 converts [m3/km2] to [mm])
    ! TODO - this was in BalanceCheckMod - not sure where it belongs?

    call c2g( bounds, &
         col_ws%endwb(bounds%begc:bounds%endc), &
         grc_ws%tws  (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = bounds%begg, bounds%endg
       grc_ws%tws(g) = grc_ws%tws(g) + atm2lnd_vars%volr_grc(g) / grc_pp%area(g) * 1.e-3_r8
    enddo


    call c2g( bounds, &
         col_es%t_grnd (bounds%begc:bounds%endc), &
         lnd2atm_vars%t_grnd_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbans', l2g_scale_type='unity' )

    do lvl = -nlevsno+1, nlevgrnd
        call c2g( bounds, &
            col_es%t_soisno(bounds%begc:bounds%endc, lvl), & 
            lnd2atm_vars%t_soisno_grc(bounds%begg:bounds%endg, lvl), &
            c2l_scale_type= 'urbans', l2g_scale_type='unity' )
    enddo
    
    call c2g( bounds, &
         soilhydrology_vars%zwt_col (bounds%begc:bounds%endc), &
         lnd2atm_vars%zwt_grc   (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbans', l2g_scale_type='unity' )    
        
    do g = bounds%begg,bounds%endg
       ! TODO temperary treatment in case weird values after c2g
       if(lnd2atm_vars%t_soisno_grc(g, 1) > 400._r8) then
           write(iulog,*)'lnd2atm_vars%t_soisno_grc(g, 1) is',lnd2atm_vars%t_soisno_grc(g, 1)
           call endrun( msg=' lnd2atm ERROR: lnd2atm_vars%t_soisno_grc >  400 Kelvin degree.'//errMsg(__FILE__, __LINE__))
       end if
       lnd2atm_vars%Tqsur_grc(g) = avg_tsoil_surf(lnd2atm_vars%t_soisno_grc(g,-nlevsno+1:nlevgrnd))
       lnd2atm_vars%Tqsub_grc(g) = avg_tsoil(lnd2atm_vars%zwt_grc(g),lnd2atm_vars%t_soisno_grc(g,-nlevsno+1:nlevgrnd))

    end do

  end subroutine lnd2atm


    function avg_tsoil_surf(Tsoil_) result(avgT_)
    ! Function for estimating average soil temperature within the top few layers (which closely interacts with surface runoff)
        implicit none
        real(r8), intent(in) :: Tsoil_(-nlevsno+1:nlevgrnd)       ! water table depth, soil temperature
        real(r8) :: avgT_             ! average soil temperature within the saturated layers
        
        integer :: ilvl   !local index
        real(r8) :: depth_(nlevsoi), h_layer(nlevsoi), sum_h, sum_ht
        
        ! calculate the thickness of each 15 soil layer, refer to Eqn. (6.5) and (6.6) in CLM4.0 tech note
        do ilvl = 1, nlevsoi
            depth_(ilvl) = 0.025_r8*(EXP(0.5_r8*(REAL(ilvl)-0.5_r8))-1._r8)
        enddo
        h_layer(1) = 0.5_r8*(depth_(1)+depth_(2))
        do ilvl = 2, nlevsoi-1
            h_layer(ilvl) = 0.5_r8*(depth_(ilvl+1)-depth_(ilvl-1))
        end do
        h_layer(nlevsoi) = depth_(nlevsoi) - depth_(nlevsoi-1)
        
        sum_h = 0._r8
        sum_ht = 0._r8 
        do ilvl = 1, 3
            sum_h = sum_h + h_layer(ilvl)
            sum_ht = sum_ht + h_layer(ilvl)*Tsoil_(ilvl)
        enddo
        avgT_ = sum_ht/sum_h
        
        return
    end function avg_tsoil_surf  

    function avg_tsoil(zwt_, Tsoil_) result(avgT_)
    ! Function for estimating average soil temperature within the saturated soil zone (which produces subsurface runoff)
        implicit none
        real(r8), intent(in) :: zwt_, Tsoil_(-nlevsno+1:nlevgrnd)       ! water table depth, soil temperature
        real(r8) :: avgT_             ! average soil temperature within the saturated layers
        
        integer :: ilvl,izwt   !local index
        real(r8) :: depth_(nlevsoi), h_layer(nlevsoi), sum_h, sum_ht
        
        ! calculate the thickness of each 15 soil layer, refer to Eqn. (6.5) and (6.6) in CLM4.0 tech note
        do ilvl = 1, nlevsoi
            depth_(ilvl) = 0.025_r8*(EXP(0.5_r8*(REAL(ilvl)-0.5_r8))-1._r8)
        enddo
        h_layer(1) = 0.5_r8*(depth_(1)+depth_(2))
        do ilvl = 2, nlevsoi-1
            h_layer(ilvl) = 0.5_r8*(depth_(ilvl+1)-depth_(ilvl-1))
        end do
        h_layer(nlevsoi) = depth_(nlevsoi) - depth_(nlevsoi-1)
        
        if(zwt_ <= 0._r8) then ! water table close to ground surface, average over the whole soil column
            sum_h = 0._r8
            sum_ht = 0._r8
            do ilvl = 1, nlevsoi
                sum_h = sum_h + h_layer(ilvl)
                sum_ht = sum_ht + h_layer(ilvl)*Tsoil_(ilvl)
            enddo
            avgT_ = sum_ht/sum_h
        else if(zwt_ >= depth_(nlevsoi)) then ! water table deeper than the total soil depth, taking the temperature of the deepes
            avgT_ = Tsoil_(nlevsoi)
        else
            sum_h = 0._r8
            sum_ht = 0._r8
            ! find out which soil layer the water table is located
            do ilvl = 1, nlevsoi
                if(zwt_ <= depth_(ilvl)) then
                    izwt = ilvl
                    sum_h = depth_(ilvl) - zwt_
                    sum_ht = (depth_(ilvl) - zwt_)*Tsoil_(ilvl)
                    exit
                end if
            enddo
            ! calculate mean soil temperature of the total saturated soil zone
            do ilvl = izwt + 1, nlevsoi
                sum_h = sum_h + h_layer(ilvl)
                sum_ht = sum_ht + h_layer(ilvl)*Tsoil_(ilvl)
            enddo
            avgT_ = sum_ht/sum_h
        end if
        
        return
    end function avg_tsoil

end module lnd2atmMod
