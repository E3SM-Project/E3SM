module SedYieldMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate the sediment flux caused by soil erosion based on the 
  ! improved Morgan model (Tan et al., 2017 & 2018) 
  !
  use shr_const_mod     , only : T0 => SHR_CONST_TKFRZ
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_varcon        , only : grav, denh2o, rpi
  use clm_varpar        , only : mxpft, nlevsno, max_patch_per_col
  use atm2lndType       , only : atm2lnd_type
  use CanopyStateType   , only : CanopyState_type
  use EnergyFluxType    , only : energyflux_type
  use SoilHydrologyType , only : soilhydrology_type
  use SoilStateType     , only : soilstate_type
  use WaterfluxType     , only : waterflux_type
  use WaterStateType    , only : waterstate_type
  use TemperatureType   , only : temperature_type
  use ColumnType        , only : col_pp
  use LandunitType      , only : lun_pp
  use VegetationType    , only : veg_pp
  use SedFluxType       , only : sedflux_type 
  use TopounitDataType  , only : top_as, top_af ! Atmospheric state and flux variables
  use ColumnDataType    , only : col_ws, col_wf
  use VegetationDataType, only : veg_wf
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilErosion          ! Calculate hillslope sediment flux
  ! !MODULE CONSTANTS
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilErosion (bounds, num_hydrologyc, filter_hydrologyc, &
    atm2lnd_vars, canopystate_vars, soilstate_vars, waterstate_vars, &
    waterflux_vars, sedflux_vars)
    !
    ! !DESCRIPTION:
    ! Calculate rainfall and runoff driven erosion 
    !
    ! !USES:
    use clm_time_manager, only : get_step_size
    use clm_varctl      , only : iulog
    use landunit_varcon , only : istcrop, istsoil
    use pftvarcon       , only : gcpsi, pftcc
    use spmdMod         , only : masterproc
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(CanopyState_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(sedflux_type)       , intent(inout) :: sedflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c, fc, p, pi, l, j                         ! indices
    integer  :: dtime                                      ! timestep size [seconds]
    integer  :: nslp                                       ! slope level #
    real(r8) :: sinslp, fslp, factor_slp                   ! topo gradient
    real(r8) :: gndbare, vegcc, fbare                      ! veg cover factors
    real(r8) :: fungrvl                                    ! ground uncovered by gravel
    real(r8) :: K, COH                                     ! soil erodibility
    real(r8) :: Qs, Ptot, Ie, Dl                           ! water fluxes
    real(r8) :: Tc, Es_Q, Es_P, KE_DT, KE_LD               ! temporaries 
    real(r8) :: stxt(4)                                    ! soil texture including gravel
    character(len=32) :: subname = 'SoilErosion'           ! subroutine name
    !-----------------------------------------------------------------------

    associate(                                                        &
         forc_rain        =>    top_af%rain                         , & ! Input: [real(r8) (:) ] rain rate [mm/s]
         forc_t           =>    top_as%tbot                         , & ! Input: [real(r8) (:) ] atmospheric temperature (Kelvin) 

         tlai             =>    canopystate_vars%tlai_patch         , & ! Input: [real(r8) (:) ] LAI
         hbot             =>    canopystate_vars%hbot_patch         , & ! Input: [real(r8) (:) ] canopy bottom (m)
         htop             =>    canopystate_vars%htop_patch         , & ! Input: [real(r8) (:) ] canopy top (m) 

         hslp_p10         =>    col_pp%hslp_p10                     , & ! Input: [real(r8) (:,:) ] hillslope gradient percentiles
    
         bd               =>    soilstate_vars%bd_col               , & ! Input: [real(r8) (:,:) ] soil bulk density (kg/m3)
         fsand            =>    soilstate_vars%cellsand_col         , & ! Input: [real(r8) (:,:) ] sand percentage
         fclay            =>    soilstate_vars%cellclay_col         , & ! Input: [real(r8) (:,:) ] clay percentage
         fgrvl            =>    soilstate_vars%cellgrvl_col         , & ! Input: [real(r8) (:,:) ] gravel percentage

         hsnow            =>    col_ws%snow_depth                   , & ! Input: [real(r8) (:) ] snow depth (m)

         pfactor          =>    sedflux_vars%pfactor_col            , & ! Input: [real(r8) (:) ] rainfall-driven erosion scaling factor
         qfactor          =>    sedflux_vars%qfactor_col            , & ! Input: [real(r8) (:) ] runoff-driven erosion scaling factor
         tfactor          =>    sedflux_vars%tfactor_col            , & ! Input: [real(r8) (:) ] transport capacity scaling factor 

         qflx_surf        =>    col_wf%qflx_surf                    , & ! Input: [real(r8) (:) ] surface runoff (mm/s)
         qflx_dirct_rain  =>    veg_wf%qflx_dirct_rain              , & ! Input: [real(r8) (:) ] direct throughfall rain (mm/s)
         qflx_leafdrip    =>    veg_wf%qflx_leafdrip                , & ! Input: [real(r8) (:) ] leaf rain drip (mm/s)

         flx_p_ero        =>    sedflux_vars%sed_p_ero_col          , & ! Output: [real(r8) (:) ] sed detached by rainfall (kg/m2/s)
         flx_q_ero        =>    sedflux_vars%sed_q_ero_col          , & ! Output: [real(r8) (:) ] sed detached by runoff (kg/m2/s)
         flx_sed_ero      =>    sedflux_vars%sed_ero_col            , & ! Output: [real(r8) (:) ] total detachment (kg/m2/s)
         flx_sed_yld      =>    sedflux_vars%sed_yld_col              & ! Output: [real(r8) (:) ] sed flux to inland waters (kg/m2/s)
         )

         ! Get time step
         dtime = get_step_size()

         ! nolakec or other col filters
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            l = col_pp%landunit(c)

            ! initialization
            flx_p_ero(c)          = 0._r8
            flx_q_ero(c)          = 0._r8
            flx_sed_ero(c)        = 0._r8
            flx_sed_yld(c)        = 0._r8

            ! check landunit type and ground covered by snow/ice
            if ( lun_pp%itype(l)/=istsoil .and. lun_pp%itype(l)/=istcrop ) then
               cycle
            else if (hsnow(c)>0) then
               cycle
            end if

            ! soil detachment by rainfall
            stxt = (/fclay(c,1), 100.-fclay(c,1)-fsand(c,1), fsand(c,1), &
                fgrvl(c,1)/)  
            K = SoilDetachability(stxt)
            COH = SoilCohesion(stxt)
            Es_P = 0._r8    ! detachment by throughfall + leap drip
            if (forc_t(c)>T0 .and. forc_rain(c)>0._r8) then
               fungrvl = 1. - 0.01 * fgrvl(c,1)
               do pi = 1, max_patch_per_col
                  if ( pi<=col_pp%npfts(c) ) then
                     p = col_pp%pfti(c) + pi - 1
                     if (veg_pp%active(p) .and. veg_pp%wtcol(p)>0) then
                        ! throughfall power
                        Ptot = qflx_dirct_rain(p) * dtime   ! mm
                        Ie = 3.6d3 * forc_rain(c)           ! mm/hr
                        KE_DT = Ptot * fungrvl * max(0._r8,8.95+8.44*log10(Ie))
                        ! leaf drip power
                        Dl = max(0._r8, qflx_leafdrip(p)*dtime)      ! mm
                        KE_LD = max(0._r8,15.8*sqrt(0.5*(htop(p)+hbot(p)))-5.87) * &
                           fungrvl * Dl
                        Es_P = Es_P + pfactor(c) * veg_pp%wtcol(p) * K * (KE_DT+KE_LD)
                     end if
                  end if
               end do
            end if
            Es_P = 1d-7 / 8.64 * Es_P   ! kg/m2/s

            ! soil detachment by runoff
            gndbare = 0._r8
            vegcc = 0._r8
            do pi = 1, max_patch_per_col
               if ( pi<=col_pp%npfts(c) ) then
                  p = col_pp%pfti(c) + pi - 1
                  if (veg_pp%active(p)) then
                     fbare = exp( -gcpsi(veg_pp%itype(p)) * tlai(p) )
                     gndbare = gndbare + veg_pp%wtcol(p) * fbare
                     vegcc = vegcc + veg_pp%wtcol(p) * pftcc(veg_pp%itype(p))
                  end if
               end if
            end do
            
            Es_Q = 0._r8
            Tc = 0._r8
            if (qflx_surf(c)>0) then
               Qs = 8.64d4 * qflx_surf(c)  ! mm/d
               nslp = size(hslp_p10(c,:)) - 1
               fslp = 1.0 / DBLE(nslp)
               factor_slp = 0._r8
               do j = 1, nslp
                  sinslp = sin(atan(max(0.5*(hslp_p10(c,j)+hslp_p10(c,j+1)),1d-4)))
                  factor_slp = factor_slp + fslp * sinslp
               end do
               Es_Q = 19.1 * qfactor(c) * 2.0 / COH * factor_slp * &
                  gndbare * Qs**1.5
               Tc = 19.1 * tfactor(c) * factor_slp * vegcc * Qs**2.0
            end if
            Es_Q = 1d-7 / 8.64 * Es_Q
            Tc = 1d-7 / 8.64 * Tc

            ! assign flux values
            flx_p_ero(c) = flx_p_ero(c) + Es_P
            flx_q_ero(c) = flx_q_ero(c) + Es_Q
            flx_sed_ero(c) = flx_sed_ero(c) + Es_P + Es_Q 
            flx_sed_yld(c) = flx_sed_yld(c) + min(Es_P+Es_Q, Tc)
         end do
    
    end associate

  end subroutine SoilErosion

  !------------------------------------------------------------------------------
  character(len=32) function SoilTextureType(stxt)
    !
    ! !DESCRIPTION:
    ! soil texture 
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: stxt(4)  ! [clay, silt, sand, gravel] in percentage
    !
    ! !LOCAL VARIABLES:
    real(r8) :: clay, silt, sand, tsum 
    !------------------------------------------------------------------------------

    tsum = sum(stxt(1:3))
    clay = stxt(1) / tsum 
    silt = stxt(2) / tsum
    sand = stxt(3) / tsum
    if (silt + 1.5*clay < 15) then
       SoilTextureType = 'sand'
    else if (silt + 2*clay < 30) then
       SoilTextureType = 'loamy sand'
    else if (((clay>=7 .and. clay<20) .and. sand>52) .or. &
            (clay<7 .and. silt<50)) then
       SoilTextureType = 'sandy loam'
    else if ((clay>=7 .and. clay<27) .and. (silt>=28 .and. silt<50) .and. &
            sand<=52) then
       SoilTextureType = 'loam'
    else if ((silt>=50 .and. (clay>=12 .and. clay<27)) .or. &
            ((silt>=50 .and. silt<80) .and. clay<12)) then
       SoilTextureType = 'silt loam'
    else if (silt>=80 .and. clay<12) then
       SoilTextureType = 'silt'
    else if ((clay>=20 .and. clay<35) .and. silt<28 .and. sand>45) then
       SoilTextureType = 'sandy clay loam'
    else if ((clay>=27 .and. clay<40) .and. (sand>20 .and. sand<=45)) then
       SoilTextureType = 'clay loam'
    else if ((clay>=27 .and. clay<40) .and. sand<=20) then
       SoilTextureType = 'silty clay loam'
    else if (clay>=35 .and. sand>45) then
       SoilTextureType = 'sandy clay'
    else if (clay>=40 .and. silt>=40) then
       SoilTextureType = 'silty clay'
    else if (clay>=40 .and. sand<=45 .and. silt<40) then
       SoilTextureType = 'clay'     
    end if

  end function SoilTextureType

  !------------------------------------------------------------------------------
  real(r8) function SoilDetachability(stxt)
    !
    ! !DESCRIPTION:
    ! soil detachability by rain drops (g J-1)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: stxt(4)  ! [clay, silt, sand, gravel] in fraction 
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: txttype
    !------------------------------------------------------------------------------

    txttype = SoilTextureType(stxt) 
    if (trim(txttype)=='clay') then
       SoilDetachability = 0.05
    else if (trim(txttype)=='silty clay') then
       SoilDetachability = 0.5
    else if (trim(txttype)=='sandy clay') then
       SoilDetachability = 0.3
    else if (trim(txttype)=='silty clay loam') then
       SoilDetachability = 0.8
    else if (trim(txttype)=='clay loam') then
       SoilDetachability = 0.7
    else if (trim(txttype)=='sandy clay loam') then
       SoilDetachability = 0.1
    else if (trim(txttype)=='silt') then
       SoilDetachability = 1.0
    else if (trim(txttype)=='silt loam') then
       SoilDetachability = 0.9
    else if (trim(txttype)=='loam') then
       SoilDetachability = 0.8
    else if (trim(txttype)=='sandy loam') then
       SoilDetachability = 0.7
    else if (trim(txttype)=='loamy sand') then
       SoilDetachability = 0.3
    else if (trim(txttype)=='sand') then
       SoilDetachability = 1.9
    end if

  end function SoilDetachability

  !------------------------------------------------------------------------------
  real(r8) function SoilCohesion(stxt)
    !
    ! !DESCRIPTION:
    ! soil cohesion against concentrated flow (kPa)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: stxt(4)  ! [clay, silt, sand, gravel] in fraction 
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: txttype
    !------------------------------------------------------------------------------

    txttype = SoilTextureType(stxt)
    if (trim(txttype)=='clay') then
       SoilCohesion = 12.0
    else if (trim(txttype)=='silty clay') then
       SoilCohesion = 10.0
    else if (trim(txttype)=='sandy clay') then
       SoilCohesion = 10.0
    else if (trim(txttype)=='silty clay loam') then
       SoilCohesion = 9.0
    else if (trim(txttype)=='clay loam') then
       SoilCohesion = 10.0
    else if (trim(txttype)=='sandy clay loam') then
       SoilCohesion = 3.0
    else if (trim(txttype)=='silt') then
       SoilCohesion = 3.0
    else if (trim(txttype)=='silt loam') then
       SoilCohesion = 3.0
    else if (trim(txttype)=='loam') then
       SoilCohesion = 3.0
    else if (trim(txttype)=='sandy loam') then
       SoilCohesion = 2.0
    else if (trim(txttype)=='loamy sand') then
       SoilCohesion = 2.0
    else if (trim(txttype)=='sand') then
       SoilCohesion = 2.0
    end if

  end function SoilCohesion

end module SedYieldMod
