module SedYieldMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate the sediment flux caused by soil erosion documented in
  ! Tan, Z., et al. (2022), Representing global soil erosion and sediment 
  ! flux in Earth System Models, J. Adv. Model. Earth Sy., 14, e2021MS002756. 
  !
  use shr_const_mod     , only : T0 => SHR_CONST_TKFRZ
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun
  use decompMod         , only : bounds_type
  use elm_varcon        , only : dzsoi_decomp
  use elm_varcon        , only : grav, denh2o, rpi
  use elm_varcon        , only : ispval
  use elm_varpar        , only : mxpft, nlevsno, max_patch_per_col
  use elm_varpar        , only : nlevslp
  use elm_varpar        , only : ndecomp_pools
  use atm2lndType       , only : atm2lnd_type
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CanopyStateType   , only : CanopyState_type
  use CNStateType       , only : cnstate_type
  use EnergyFluxType    , only : energyflux_type
  use SoilHydrologyType , only : soilhydrology_type
  use SoilStateType     , only : soilstate_type
  use TemperatureType   , only : temperature_type
  use GridcellType      , only : grc_pp
  use TopounitType      , only : top_pp
  use ColumnType        , only : col_pp
  use LandunitType      , only : lun_pp
  use VegetationType    , only : veg_pp
  use SedFluxType       , only : sedflux_type 
  use TopounitDataType  , only : top_as, top_af ! Atmospheric state and flux variables
  use ColumnDataType    , only : col_ws, col_wf, col_cs
  use VegetationDataType, only : veg_wf, veg_cs

  use timeinfoMod
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
  subroutine SoilErosion (bounds, num_soilc, filter_soilc, &
    canopystate_vars, cnstate_vars, soilstate_vars, sedflux_vars)
    !
    ! !DESCRIPTION:
    ! Calculate rainfall and runoff driven erosion 
    ! It should be noted that because glacier runoff is now inactive for
    ! inland grid cells, the glacier factor is implemented on surface
    ! runoff.
    !
    ! !USES:
    use clm_time_manager, only : get_step_size
    use landunit_varcon , only : istcrop, istsoil, istice
    use pftvarcon       , only : gcbc_p, gcbc_q, gcbr_p, gcbr_q 
    use pftvarcon       , only : nc4_grass
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_soilc(:) ! column filter for soil points
    type(CanopyState_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(sedflux_type)       , intent(inout) :: sedflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c, fc, p, t, l, g, j                       ! indices
    integer  :: lg, cg, lt                                 ! indices of glaciers
    integer  :: lp                                         ! indices of soil pools
    integer  :: dtime                                      ! timestep size [seconds]
    logical  :: found                                      ! flags
    real(r8) :: sinslp, frac_slp                           ! topo gradient
    real(r8) :: fslp, fslp_tc                              ! slope factor
    real(r8) :: fgndcov, fsr, ftillage_tc                  ! factors of ground cover and roughness
    real(r8) :: fungrvl                                    ! ground uncovered by gravel
    real(r8) :: ftillage, flitho, fglacier                 ! factors of tillage, lithology and glacier
    real(r8) :: wglc                                       ! weight of glacier land unit
    real(r8) :: Brsd, Broot                                ! biomass of residue and roots
    real(r8) :: Crsd, Clai, PCT_gnd                        ! ground cover calculated from residue and LAI
    real(r8) :: nh                                         ! Manning's coefficient 
    real(r8) :: K, COH                                     ! soil erodibility
    real(r8) :: Qs, Qss, Qg, Ptot, Ie, Dl                  ! water fluxes
    real(r8) :: Tc, Es_Q, Es_P, KE_DT, KE_LD, cheight      ! temporaries 
    real(r8) :: Es_Qcrp, Es_Pcrp                           ! cropland temporaries
    real(r8) :: stxt(4)                                    ! soil texture including gravel
    character(len=32) :: subname = 'SoilErosion'           ! subroutine name
    !-----------------------------------------------------------------------

    associate(                                                        &
         forc_rain        =>    top_af%rain                         , & ! Input: [real(r8) (:) ] rain rate [mm/s]
         forc_t           =>    top_as%tbot                         , & ! Input: [real(r8) (:) ] atmospheric temperature (Kelvin)

         tlai             =>    canopystate_vars%tlai_patch         , & ! Input: [real(r8) (:) ] LAI
         hbot             =>    canopystate_vars%hbot_patch         , & ! Input: [real(r8) (:) ] canopy bottom (m)
         htop             =>    canopystate_vars%htop_patch         , & ! Input: [real(r8) (:) ] canopy top (m)

         froot_prof       =>    cnstate_vars%froot_prof_patch       , & ! Input: [real(r8) (:,:) ] fine root vertical profile (1/m)
         croot_prof       =>    cnstate_vars%croot_prof_patch       , & ! Input: [real(r8) (:,:) ] coarse root vertical profile (1/m)

         hslp_p10         =>    col_pp%hslp_p10                     , & ! Input: [real(r8) (:,:) ] hillslope gradient percentiles

         bd               =>    soilstate_vars%bd_col               , & ! Input: [real(r8) (:,:) ] soil bulk density (kg/m3)
         fsand            =>    soilstate_vars%cellsand_col         , & ! Input: [real(r8) (:,:) ] sand percentage
         fclay            =>    soilstate_vars%cellclay_col         , & ! Input: [real(r8) (:,:) ] clay percentage
         fgrvl            =>    soilstate_vars%cellgrvl_col         , & ! Input: [real(r8) (:,:) ] gravel percentage

         tillage          =>    soilstate_vars%tillage_col          , & ! Input: [real(r8) (:) ] conserved tillage fraction 
         litho            =>    soilstate_vars%litho_col            , & ! Input: [real(r8) (:) ] lithology erodiblity index

         decomp_cpools_vr =>    col_cs%decomp_cpools_vr             , & ! Input: [real(r8) (:,:,:) ] soil carbon pools [gC/m3]

         frac_sno         =>    col_ws%frac_sno                     , & ! Input: [real(r8) (:) ] fraction of ground covered by snow (0 to 1)

         pfactor          =>    sedflux_vars%pfactor_col            , & ! Input: [real(r8) (:) ] rainfall-driven erosion scaling factor
         qfactor          =>    sedflux_vars%qfactor_col            , & ! Input: [real(r8) (:) ] runoff-driven erosion scaling factor
         tfactor          =>    sedflux_vars%tfactor_col            , & ! Input: [real(r8) (:) ] transport capacity scaling factor

         qflx_surf        =>    col_wf%qflx_surf                    , & ! Input: [real(r8) (:) ] surface runoff (mm/s)
         qflx_qrgwl       =>    col_wf%qflx_qrgwl                   , & ! Input: [real(r8) (:) ] glacier runoff (mm/s) 
         qflx_dirct_rain  =>    veg_wf%qflx_dirct_rain              , & ! Input: [real(r8) (:) ] direct throughfall rain (mm/s)
         qflx_leafdrip    =>    veg_wf%qflx_leafdrip                , & ! Input: [real(r8) (:) ] leaf rain drip (mm/s)
         qflx_real_irrig  =>    veg_wf%qflx_real_irrig_patch        , & ! Input: [real(r8) (:) ] actual irrigation amount (mm/s)

         flx_p_ero        =>    sedflux_vars%sed_p_ero_col          , & ! Output: [real(r8) (:) ] sed detached by rainfall (kg/m2/s)
         flx_q_ero        =>    sedflux_vars%sed_q_ero_col          , & ! Output: [real(r8) (:) ] sed detached by runoff (kg/m2/s)
         flx_sed_ero      =>    sedflux_vars%sed_ero_col            , & ! Output: [real(r8) (:) ] total detachment (kg/m2/s)
         flx_sed_crop_ero =>    sedflux_vars%sed_crop_ero_col       , & ! Output: [real(r8) (:) ] cropland detachment (kg/m2/s)
         flx_sed_yld      =>    sedflux_vars%sed_yld_col              & ! Output: [real(r8) (:) ] sed flux to inland waters (kg/m2/s)
         )

         ! Get time step
         dtime = dtime_mod

         ! nolakec or other col filters
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            l = col_pp%landunit(c)
            t = col_pp%topounit(c)
            g = col_pp%gridcell(c)

            ! initialization
            flx_p_ero(c)          = 0._r8
            flx_q_ero(c)          = 0._r8
            flx_sed_ero(c)        = 0._r8
            flx_sed_crop_ero(c)   = 0._r8
            flx_sed_yld(c)        = 0._r8

            ! check landunit type and ground covered by snow/ice
            if ( lun_pp%itype(l)/=istsoil .and. lun_pp%itype(l)/=istcrop ) then
               cycle
            end if

            ! check the glacier landunit in the topounit
            lt = top_pp%landunit_indices(istice, t) 
            fglacier = 1._r8
            if (lt /= ispval) then
               found = .false.
               cg = lun_pp%coli(lt)
               do while (.not. found .and. cg <= lun_pp%colf(lt))
                  if (col_pp%active(cg)) then
                     found = .true.
                  else
                     cg = cg + 1
                  end if
               end do
               if (found .and. qflx_qrgwl(cg)<1.e10_r8) then
                  fglacier = 1._r8 + 9._r8 * lun_pp%wttopounit(lt)
               end if
            end if

            ! tillage, lithology factors
            ftillage = 2.7_r8 - 1.7_r8 * tillage(c)
            flitho = litho(c)

            ! ground cover calculated by plant residue
            Brsd = 0.0_r8  ! kg/m2
            do lp = 1, ndecomp_pools
               if ( decomp_cascade_con%is_litter(lp) .and. decomp_cpools_vr(c,1,lp)>0._r8 ) then
                  Brsd = Brsd + 1.e-3_r8 * decomp_cpools_vr(c,1,lp) * dzsoi_decomp(1) 
               end if
            end do
            Crsd = 1._r8 - exp(-6.68_r8 * Brsd)

            stxt = (/fclay(c,1), 100._r8-fclay(c,1)-fsand(c,1), fsand(c,1), &
                fgrvl(c,1)/)
            K = SoilDetachability(stxt)
            COH = SoilCohesion(stxt)
            
            Es_P = 0._r8    ! detachment by throughfall + leap drip
            Es_Pcrp = 0._r8 ! cropland detachment by throughfall + leap drip
            if (forc_t(t)>T0 .and. forc_rain(t)>0._r8) then
               fungrvl = 1._r8 - 0.01_r8 * fgrvl(c,1)
               do p = col_pp%pfti(c), col_pp%pftf(c)
                  if (veg_pp%active(p) .and. veg_pp%wtcol(p)>0._r8) then
                     ! throughfall power
                     Ptot = (qflx_dirct_rain(p) + qflx_real_irrig(p)) * dtime    ! mm
                     Ie = 3.6e3_r8 * (forc_rain(t) + qflx_real_irrig(p))         ! mm/hr
                     KE_DT = Ptot * fungrvl * max(8.95_r8+8.44_r8*log10(Ie), 0._r8)
                     ! leaf drip power
                     Dl = max(qflx_leafdrip(p)*dtime, 0._r8)      ! mm
                     cheight = 0.5_r8*min(htop(p)+hbot(p),40._r8) ! m
                     KE_LD = max(15.8_r8*sqrt(cheight)-5.87_r8, 0._r8) * fungrvl * Dl
                     ! LAI and root biomass (kgC/m3): OM/OC = 1.72 
                     Clai = 1._r8 - exp(-tlai(p))
                     PCT_gnd = 100._r8 * max(Crsd,Clai) ! ground cover in percentage
                     Broot = 1.e-3_r8 * ( veg_cs%frootc(p)*froot_prof(p,1) + &
                        (veg_cs%livecrootc(p)+veg_cs%deadcrootc(p))*croot_prof(p,1) )
                     fgndcov = exp( -gcbc_p(veg_pp%itype(p))*PCT_gnd - &
                        gcbr_p(veg_pp%itype(p))*Broot ) 
                     if( veg_pp%itype(p) > nc4_grass )then
                        Es_Pcrp = Es_Pcrp + pfactor(c) * ftillage * flitho * &
                           fgndcov * veg_pp%wtcol(p) * K * (KE_DT+KE_LD)
                        
                        Es_P = Es_P + pfactor(c) * ftillage * flitho * &
                           fgndcov * veg_pp%wtcol(p) * K * (KE_DT+KE_LD)
                     else
                        Es_P = Es_P + pfactor(c) * flitho * fgndcov * &
                           veg_pp%wtcol(p) * K * (KE_DT+KE_LD)
                     end if 
                  end if
               end do
            end if
            Es_P = 1.e-3_r8 / dtime * (1._r8 - frac_sno(c)) * Es_P        ! kg/m2/s 
            Es_Pcrp = 1.e-3_r8 / dtime * (1._r8 - frac_sno(c)) * Es_Pcrp  ! kg/m2/s

            ! snow scaling factor from T factor of BQART
            Qs = 8.64e4_r8 * qflx_surf(c)                ! mm/d
            Qss = (1._r8 - 0.7846_r8*frac_sno(c)) * Qs   ! mm/d

            Es_Q = 0._r8
            Es_Qcrp = 0._r8
            Tc = 0._r8
            if (Qs>0._r8) then
               frac_slp = 1.0_r8 / DBLE(nlevslp-1)
               fslp = 0._r8
               fslp_tc = 0._r8
               do j = 1, nlevslp-1
                  sinslp = sin(atan(max(0.5_r8*(hslp_p10(c,j)+hslp_p10(c,j+1)),1.e-4_r8)))
                  fslp = fslp + frac_slp * sinslp
                  fslp_tc = fslp_tc + frac_slp * sinslp**1.25_r8
               end do
               
               fsr = 0._r8
               ftillage_tc = 0._r8 
               do p = col_pp%pfti(c), col_pp%pftf(c)
                  if (veg_pp%active(p) .and. veg_pp%wtcol(p)>0._r8) then
                     ! LAI and root biomass (kgC/m3): OM/OC = 1.72
                     Clai = 1._r8 - exp(-tlai(p))
                     PCT_gnd = 100._r8 * max(Crsd,Clai) ! ground cover in percentage
                     Broot = 1.e-3_r8 * ( veg_cs%frootc(p)*froot_prof(p,1) + &
                        (veg_cs%livecrootc(p)+veg_cs%deadcrootc(p))*croot_prof(p,1) )
                     fgndcov = exp( -gcbc_q(veg_pp%itype(p))*PCT_gnd - &
                        gcbr_q(veg_pp%itype(p))*Broot )

                     nh = 0.03_r8 + 0.05_r8*max(Crsd,Clai)
                     fsr = fsr + veg_pp%wtcol(p) * (0.03_r8/nh)**0.6_r8
                  
                     if ( veg_pp%itype(p) > nc4_grass ) then
                        ftillage_tc = ftillage_tc + ftillage * veg_pp%wtcol(p)

                        Es_Q = Es_Q + 19.1_r8 * qfactor(c) * 2._r8/COH * flitho * fslp * &
                           ftillage * fgndcov * Qss**1.5_r8 * veg_pp%wtcol(p)

                        Es_Qcrp = Es_Qcrp + 19.1_r8 * qfactor(c) * 2._r8/COH * flitho * fslp * &
                           ftillage * fgndcov * Qss**1.5_r8 * veg_pp%wtcol(p)
                     else
                        ftillage_tc = ftillage_tc + veg_pp%wtcol(p)

                        Es_Q = Es_Q + 19.1_r8 * qfactor(c) * 2._r8/COH * flitho * fslp * &
                           fgndcov * fglacier * Qss**1.5_r8 * veg_pp%wtcol(p)
                     end if
                  end if
               end do

               Tc = 19.1_r8 * tfactor(c) * fslp_tc * fsr * ftillage_tc * &
                  flitho * fglacier * Qs**2._r8
            end if
            Es_Q = 1.e-7_r8 / 8.64_r8 * Es_Q        ! kg/m2/s
            Es_Qcrp = 1.e-7_r8 / 8.64_r8 * Es_Qcrp  ! kg/m2/s
            Tc = 1.e-7_r8 / 8.64_r8 * Tc            ! kg/m2/s

            ! assign flux values
            flx_p_ero(c) = flx_p_ero(c) + Es_P
            flx_q_ero(c) = flx_q_ero(c) + Es_Q
            flx_sed_ero(c) = flx_sed_ero(c) + Es_P + Es_Q
            flx_sed_crop_ero(c) = flx_sed_crop_ero(c) + Es_Pcrp + Es_Qcrp
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

    ! calculate the percentage of clay, silt and sand in surface soil
    ! assuming the total is 100% 
    tsum = sum(stxt(1:3))
    clay = 100._r8 * stxt(1) / tsum 
    silt = 100._r8 * stxt(2) / tsum
    sand = 100._r8 * stxt(3) / tsum
    if ( silt + 1.5_r8*clay < 15._r8 ) then
       SoilTextureType = 'sand'
    else if ( silt + 2.0_r8*clay < 30._r8 ) then
       SoilTextureType = 'loamy sand'
    else if ( ((clay>=7._r8 .and. clay<20._r8) .and. sand>52._r8) .or. &
            (clay<7._r8 .and. silt<50._r8) ) then
       SoilTextureType = 'sandy loam'
    else if ( (clay>=7._r8 .and. clay<27._r8) .and. &
            (silt>=28._r8 .and. silt<50._r8) .and. sand<=52._r8 ) then
       SoilTextureType = 'loam'
    else if ( (silt>=50._r8 .and. (clay>=12._r8 .and. clay<27._r8)) .or. &
            ((silt>=50._r8 .and. silt<80._r8) .and. clay<12._r8) ) then
       SoilTextureType = 'silt loam'
    else if ( silt>=80._r8 .and. clay<12._r8 ) then
       SoilTextureType = 'silt'
    else if ( (clay>=20._r8 .and. clay<35._r8) .and. &
            silt<28._r8 .and. sand>45._r8 ) then
       SoilTextureType = 'sandy clay loam'
    else if ( (clay>=27._r8 .and. clay<40._r8) .and. &
            (sand>20._r8 .and. sand<=45._r8) ) then
       SoilTextureType = 'clay loam'
    else if ( (clay>=27._r8 .and. clay<40._r8) .and. sand<=20._r8 ) then
       SoilTextureType = 'silty clay loam'
    else if ( clay>=35._r8 .and. sand>45._r8 ) then
       SoilTextureType = 'sandy clay'
    else if ( clay>=40._r8 .and. silt>=40._r8 ) then
       SoilTextureType = 'silty clay'
    else if ( clay>=40._r8 .and. sand<=45._r8 .and. silt<40._r8 ) then
       SoilTextureType = 'clay'
    else
       call endrun(msg=' ERROR: soil clay, silt and sand are out of bounds.'//&
            errMsg(__FILE__, __LINE__))
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
       SoilDetachability = 0.05_r8
    else if (trim(txttype)=='silty clay') then
       SoilDetachability = 0.5_r8
    else if (trim(txttype)=='sandy clay') then
       SoilDetachability = 0.3_r8
    else if (trim(txttype)=='silty clay loam') then
       SoilDetachability = 0.8_r8
    else if (trim(txttype)=='clay loam') then
       SoilDetachability = 0.7_r8
    else if (trim(txttype)=='sandy clay loam') then
       SoilDetachability = 0.1_r8
    else if (trim(txttype)=='silt') then
       SoilDetachability = 1.0_r8
    else if (trim(txttype)=='silt loam') then
       SoilDetachability = 0.9_r8
    else if (trim(txttype)=='loam') then
       SoilDetachability = 0.8_r8
    else if (trim(txttype)=='sandy loam') then
       SoilDetachability = 0.7_r8
    else if (trim(txttype)=='loamy sand') then
       SoilDetachability = 0.3_r8
    else if (trim(txttype)=='sand') then
       SoilDetachability = 1.9_r8
    else
       call endrun(msg=' ERROR: no soil texture type is found.'//&
            errMsg(__FILE__, __LINE__))
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
       SoilCohesion = 12.0_r8
    else if (trim(txttype)=='silty clay') then
       SoilCohesion = 10.0_r8
    else if (trim(txttype)=='sandy clay') then
       SoilCohesion = 10.0_r8
    else if (trim(txttype)=='silty clay loam') then
       SoilCohesion = 9.0_r8
    else if (trim(txttype)=='clay loam') then
       SoilCohesion = 10.0_r8
    else if (trim(txttype)=='sandy clay loam') then
       SoilCohesion = 3.0_r8
    else if (trim(txttype)=='silt') then
       SoilCohesion = 3.0_r8
    else if (trim(txttype)=='silt loam') then
       SoilCohesion = 3.0_r8
    else if (trim(txttype)=='loam') then
       SoilCohesion = 3.0_r8
    else if (trim(txttype)=='sandy loam') then
       SoilCohesion = 2.0_r8
    else if (trim(txttype)=='loamy sand') then
       SoilCohesion = 2.0_r8
    else if (trim(txttype)=='sand') then
       SoilCohesion = 2.0_r8
    else
       call endrun(msg=' ERROR: no soil texture type is found.'//&
            errMsg(__FILE__, __LINE__))
    end if

  end function SoilCohesion

end module SedYieldMod
