module SnowHydrologyMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SnowHydrologyMod
!
! !DESCRIPTION:
! Calculate snow hydrology.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar  , only: nlevsno
  use clm_varctl  , only: iulog
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SnowWater         ! Change of snow mass and the snow water onto soil
  public :: SnowCompaction    ! Change in snow layer thickness due to compaction
  public :: CombineSnowLayers ! Combine snow layers less than a min thickness
  public :: DivideSnowLayers  ! Subdivide snow layers if they exceed maximum thickness
  public :: DivideSnowLayers_Lake ! Adjusted so that snow layer thicknesses are larger over lakes
  public :: BuildSnowFilter   ! Construct snow/no-snow filters
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: Combo            ! Returns the combined variables: dz, t, wliq, wice.
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SnowWater
!
! !INTERFACE:
  subroutine SnowWater(lbc, ubc, num_snowc, filter_snowc, &
                       num_nosnowc, filter_nosnowc)
!
! !DESCRIPTION:
! Evaluate the change of snow mass and the snow water onto soil.
! Water flow within snow is computed by an explicit and non-physical
! based scheme, which permits a part of liquid water over the holding
! capacity (a tentative value is used, i.e. equal to 0.033*porosity) to
! percolate into the underlying layer.  Except for cases where the
! porosity of one of the two neighboring layers is less than 0.05, zero
! flow is assumed. The water flow out of the bottom of the snow pack will
! participate as the input of the soil water and runoff.  This subroutine
! uses a filter for columns containing snow which must be constructed prior
! to being called.
!
! !USES:
    use clmtype
    use clm_varcon        , only : denh2o, denice, wimp, ssi, isturb,istsoil,istdlak
    use clm_time_manager  , only : get_step_size
    use clm_atmlnd        , only : clm_a2l
    use SNICARMod         , only : scvng_fct_mlt_bcphi, scvng_fct_mlt_bcpho, &
                                   scvng_fct_mlt_ocphi, scvng_fct_mlt_ocpho, &
                                   scvng_fct_mlt_dst1,  scvng_fct_mlt_dst2,  &
                                   scvng_fct_mlt_dst3,  scvng_fct_mlt_dst4
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column bounds
    integer, intent(in) :: num_snowc                   ! number of snow points in column filter
    integer, intent(in) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
    integer, intent(in) :: num_nosnowc                 ! number of non-snow points in column filter
    integer, intent(in) :: filter_nosnowc(ubc-lbc+1)   ! column filter for non-snow points
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 15 November 2000: Mariana Vertenstein
! 2/26/02, Peter Thornton: Migrated to new data structures.
! 03/28/08, Mark Flanner: Added aerosol deposition and flushing with meltwater
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: qflx_snow_melt(:)   ! net snow melt
    integer, pointer :: clandunit(:)         ! columns's landunit
    integer, pointer :: ltype(:)             ! landunit type
    real(r8), pointer :: frac_sno_eff(:)     ! eff. fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_sno(:)         ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: int_snow(:)         ! integrated snowfall [mm]
    real(r8), pointer :: h2osno(:)           ! snow water (mm H2O)
    real(r8), pointer :: qflx_ev_snow(:)     ! evaporation flux from snow (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_ev_soil(:)     ! evaporation flux from soil (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi(:)    ! evaporation flux from soil (W/m**2) [+ to atm]
    integer , pointer :: snl(:)              !number of snow layers
    logical , pointer :: do_capsnow(:)       !true => do snow capping
    real(r8), pointer :: qflx_snomelt(:)     !snow melt (mm H2O /s)
    real(r8), pointer :: qflx_rain_grnd(:)   !rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_sub_snow(:)    !sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_evap_grnd(:)   !ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_dew_snow(:)    !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd(:)    !ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: dz(:,:)             !layer depth (m)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: qflx_top_soil(:)     !net water input into soil from top (mm/s)
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: h2osoi_ice(:,:)     !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)     !liquid water (kg/m2)
    integer , pointer :: cgridcell(:)        ! columns's gridcell (col) 
    real(r8), pointer :: mss_bcphi(:,:)      ! hydrophillic BC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bcpho(:,:)      ! hydrophobic BC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocphi(:,:)      ! hydrophillic OC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocpho(:,:)      ! hydrophobic OC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst1(:,:)       ! mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst2(:,:)       ! mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst3(:,:)       ! mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst4(:,:)       ! mass of dust species 4 in snow (col,lyr) [kg]
    real(r8), pointer :: flx_bc_dep_dry(:)   ! dry BC deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_bc_dep_wet(:)   ! wet BC deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_bc_dep(:)       ! total BC deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_bc_dep_pho(:)   ! hydrophobic BC deposition (col) [kg m-1 s-1]
    real(r8), pointer :: flx_bc_dep_phi(:)   ! hydrophillic BC deposition (col) [kg m-1 s-1]
    real(r8), pointer :: flx_oc_dep_dry(:)   ! dry OC deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_oc_dep_wet(:)   ! wet OC deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_oc_dep(:)       ! total OC deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_oc_dep_pho(:)   ! hydrophobic OC deposition (col) [kg m-1 s-1]
    real(r8), pointer :: flx_oc_dep_phi(:)   ! hydrophillic OC deposition (col) [kg m-1 s-1]
    real(r8), pointer :: flx_dst_dep_dry1(:) ! dry dust (species 1) deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_dst_dep_wet1(:) ! wet dust (species 1) deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_dst_dep_dry2(:) ! dry dust (species 2) deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_dst_dep_wet2(:) ! wet dust (species 2) deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_dst_dep_dry3(:) ! dry dust (species 3) deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_dst_dep_wet3(:) ! wet dust (species 3) deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_dst_dep_dry4(:) ! dry dust (species 4) deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_dst_dep_wet4(:) ! wet dust (species 4) deposition (col) [kg m-2 s-1]
    real(r8), pointer :: flx_dst_dep(:)      ! total dust deposition (col) [kg m-2 s-1]
    real(r8), pointer :: forc_aer(:,:)       ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: c, j, fc, l                        !do loop/array indices
    real(r8) :: dtime                              !land model time step (sec)
    real(r8) :: qin(lbc:ubc)                       !water flow into the elmement (mm/s) 
    real(r8) :: qout(lbc:ubc)                      !water flow out of the elmement (mm/s)
    real(r8) :: wgdif                              !ice mass after minus sublimation
    real(r8) :: vol_liq(lbc:ubc,-nlevsno+1:0)      !partial volume of liquid water in layer
    real(r8) :: vol_ice(lbc:ubc,-nlevsno+1:0)      !partial volume of ice lens in layer
    real(r8) :: eff_porosity(lbc:ubc,-nlevsno+1:0) !effective porosity = porosity - vol_ice
    integer  :: g                                 ! gridcell loop index
    real(r8) :: qin_bc_phi(lbc:ubc)               ! flux of hydrophilic BC into layer [kg]
    real(r8) :: qout_bc_phi(lbc:ubc)              ! flux of hydrophilic BC out of layer [kg]
    real(r8) :: qin_bc_pho(lbc:ubc)               ! flux of hydrophobic BC into layer [kg]
    real(r8) :: qout_bc_pho(lbc:ubc)              ! flux of hydrophobic BC out of layer [kg]
    real(r8) :: qin_oc_phi(lbc:ubc)               ! flux of hydrophilic OC into layer [kg]
    real(r8) :: qout_oc_phi(lbc:ubc)              ! flux of hydrophilic OC out of layer [kg]
    real(r8) :: qin_oc_pho(lbc:ubc)               ! flux of hydrophobic OC into layer [kg]
    real(r8) :: qout_oc_pho(lbc:ubc)              ! flux of hydrophobic OC out of layer [kg]
    real(r8) :: qin_dst1(lbc:ubc)                 ! flux of dust species 1 into layer [kg]
    real(r8) :: qout_dst1(lbc:ubc)                ! flux of dust species 1 out of layer [kg]
    real(r8) :: qin_dst2(lbc:ubc)                 ! flux of dust species 2 into layer [kg]
    real(r8) :: qout_dst2(lbc:ubc)                ! flux of dust species 2 out of layer [kg]
    real(r8) :: qin_dst3(lbc:ubc)                 ! flux of dust species 3 into layer [kg]
    real(r8) :: qout_dst3(lbc:ubc)                ! flux of dust species 3 out of layer [kg]
    real(r8) :: qin_dst4(lbc:ubc)                 ! flux of dust species 4 into layer [kg]
    real(r8) :: qout_dst4(lbc:ubc)                ! flux of dust species 4 out of layer [kg]
    real(r8) :: mss_liqice                        ! mass of liquid+ice in a layer
 
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtype components (column-level)

    frac_sno_eff     => clm3%g%l%c%cps%frac_sno_eff 
    frac_sno         => clm3%g%l%c%cps%frac_sno 
    clandunit        => clm3%g%l%c%landunit
    ltype            => clm3%g%l%itype
    int_snow         => clm3%g%l%c%cws%int_snow
    h2osno           => clm3%g%l%c%cws%h2osno
    qflx_ev_snow     => clm3%g%l%c%cwf%pwf_a%qflx_ev_snow
    qflx_ev_soil     => clm3%g%l%c%cwf%pwf_a%qflx_ev_soil
    qflx_evap_soi    => clm3%g%l%c%cwf%pwf_a%qflx_evap_soi
    qflx_snow_melt   => clm3%g%l%c%cwf%qflx_snow_melt
    snl              => clm3%g%l%c%cps%snl
    do_capsnow       => clm3%g%l%c%cps%do_capsnow
    qflx_snomelt     => clm3%g%l%c%cwf%qflx_snomelt
    qflx_rain_grnd   => clm3%g%l%c%cwf%pwf_a%qflx_rain_grnd
    qflx_sub_snow    => clm3%g%l%c%cwf%pwf_a%qflx_sub_snow
    qflx_evap_grnd   => clm3%g%l%c%cwf%pwf_a%qflx_evap_grnd
    qflx_dew_snow    => clm3%g%l%c%cwf%pwf_a%qflx_dew_snow
    qflx_dew_grnd    => clm3%g%l%c%cwf%pwf_a%qflx_dew_grnd
    qflx_top_soil    => clm3%g%l%c%cwf%qflx_top_soil
    dz               => clm3%g%l%c%cps%dz
    h2osoi_ice       => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq       => clm3%g%l%c%cws%h2osoi_liq
    cgridcell        => clm3%g%l%c%gridcell
    mss_bcphi        => clm3%g%l%c%cps%mss_bcphi
    mss_bcpho        => clm3%g%l%c%cps%mss_bcpho
    mss_ocphi        => clm3%g%l%c%cps%mss_ocphi
    mss_ocpho        => clm3%g%l%c%cps%mss_ocpho
    mss_dst1         => clm3%g%l%c%cps%mss_dst1
    mss_dst2         => clm3%g%l%c%cps%mss_dst2
    mss_dst3         => clm3%g%l%c%cps%mss_dst3
    mss_dst4         => clm3%g%l%c%cps%mss_dst4
    flx_bc_dep       => clm3%g%l%c%cwf%flx_bc_dep
    flx_bc_dep_wet   => clm3%g%l%c%cwf%flx_bc_dep_wet
    flx_bc_dep_dry   => clm3%g%l%c%cwf%flx_bc_dep_dry
    flx_bc_dep_phi   => clm3%g%l%c%cwf%flx_bc_dep_phi
    flx_bc_dep_pho   => clm3%g%l%c%cwf%flx_bc_dep_pho
    flx_oc_dep       => clm3%g%l%c%cwf%flx_oc_dep
    flx_oc_dep_wet   => clm3%g%l%c%cwf%flx_oc_dep_wet
    flx_oc_dep_dry   => clm3%g%l%c%cwf%flx_oc_dep_dry
    flx_oc_dep_phi   => clm3%g%l%c%cwf%flx_oc_dep_phi
    flx_oc_dep_pho   => clm3%g%l%c%cwf%flx_oc_dep_pho
    flx_dst_dep      => clm3%g%l%c%cwf%flx_dst_dep
    flx_dst_dep_wet1 => clm3%g%l%c%cwf%flx_dst_dep_wet1
    flx_dst_dep_dry1 => clm3%g%l%c%cwf%flx_dst_dep_dry1
    flx_dst_dep_wet2 => clm3%g%l%c%cwf%flx_dst_dep_wet2
    flx_dst_dep_dry2 => clm3%g%l%c%cwf%flx_dst_dep_dry2
    flx_dst_dep_wet3 => clm3%g%l%c%cwf%flx_dst_dep_wet3
    flx_dst_dep_dry3 => clm3%g%l%c%cwf%flx_dst_dep_dry3
    flx_dst_dep_wet4 => clm3%g%l%c%cwf%flx_dst_dep_wet4
    flx_dst_dep_dry4 => clm3%g%l%c%cwf%flx_dst_dep_dry4
    forc_aer         => clm_a2l%forc_aer

    ! Determine model time step

    dtime = get_step_size()

    ! Renew the mass of ice lens (h2osoi_ice) and liquid (h2osoi_liq) in the
    ! surface snow layer resulting from sublimation (frost) / evaporation (condense)

    do fc = 1,num_snowc
       c = filter_snowc(fc)
       l=clandunit(c)

       if (do_capsnow(c)) then
          wgdif = h2osoi_ice(c,snl(c)+1) - frac_sno_eff(c)*qflx_sub_snow(c)*dtime
          h2osoi_ice(c,snl(c)+1) = wgdif
          if (wgdif < 0._r8) then
             h2osoi_ice(c,snl(c)+1) = 0._r8
             h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
          end if
             h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) &
                  - frac_sno_eff(c)*qflx_evap_grnd(c) * dtime
       else
          wgdif = h2osoi_ice(c,snl(c)+1) &
                  + frac_sno_eff(c) * (qflx_dew_snow(c) - qflx_sub_snow(c)) * dtime
          h2osoi_ice(c,snl(c)+1) = wgdif
          if (wgdif < 0._r8) then
             h2osoi_ice(c,snl(c)+1) = 0._r8
             h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
          end if
          h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) +  &
                  frac_sno_eff(c) * (qflx_rain_grnd(c) + qflx_dew_grnd(c) &
                  - qflx_evap_grnd(c)) * dtime
          end if
          ! if negative, reduce deeper layer's liquid water content sequentially
          if(h2osoi_liq(c,snl(c)+1) < 0._r8) then
             do j = snl(c)+1, 1
                wgdif=h2osoi_liq(c,j)
                if (wgdif >= 0._r8) exit
                h2osoi_liq(c,j) = 0._r8
                h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + wgdif
             enddo
       end if

    end do

    ! Porosity and partial volume

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             ! need to scale dz by frac_sno to convert to grid cell average depth
             vol_ice(c,j)      = min(1._r8, h2osoi_ice(c,j)/(dz(c,j)*frac_sno_eff(c)*denice))
             eff_porosity(c,j) = 1._r8 - vol_ice(c,j)
             vol_liq(c,j)      = min(eff_porosity(c,j),h2osoi_liq(c,j)/(dz(c,j)*frac_sno_eff(c)*denh2o))
          end if
       end do
    end do

    ! Capillary forces within snow are usually two or more orders of magnitude
    ! less than those of gravity. Only gravity terms are considered.
    ! the genernal expression for water flow is "K * ss**3", however,
    ! no effective parameterization for "K".  Thus, a very simple consideration
    ! (not physically based) is introduced:
    ! when the liquid water of layer exceeds the layer's holding
    ! capacity, the excess meltwater adds to the underlying neighbor layer.

    ! Also compute aerosol fluxes through snowpack in this loop:
    ! 1) compute aerosol mass in each layer
    ! 2) add aerosol mass flux from above layer to mass of this layer
    ! 3) qout_xxx is mass flux of aerosol species xxx out bottom of
    !    layer in water flow, proportional to (current) concentration
    !    of aerosol in layer multiplied by a scavenging ratio.
    ! 4) update mass of aerosol in top layer, accordingly
    ! 5) update mass concentration of aerosol accordingly

    qin(:) = 0._r8
    qin_bc_phi(:) = 0._r8
    qin_bc_pho(:) = 0._r8
    qin_oc_phi(:) = 0._r8
    qin_oc_pho(:) = 0._r8
    qin_dst1(:)   = 0._r8
    qin_dst2(:)   = 0._r8
    qin_dst3(:)   = 0._r8
    qin_dst4(:)   = 0._r8

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             h2osoi_liq(c,j) = h2osoi_liq(c,j) + qin(c)
             
             mss_bcphi(c,j) = mss_bcphi(c,j) + qin_bc_phi(c)
             mss_bcpho(c,j) = mss_bcpho(c,j) + qin_bc_pho(c)
             mss_ocphi(c,j) = mss_ocphi(c,j) + qin_oc_phi(c)
             mss_ocpho(c,j) = mss_ocpho(c,j) + qin_oc_pho(c)
             mss_dst1(c,j)  = mss_dst1(c,j) + qin_dst1(c)
             mss_dst2(c,j)  = mss_dst2(c,j) + qin_dst2(c)
             mss_dst3(c,j)  = mss_dst3(c,j) + qin_dst3(c)
             mss_dst4(c,j)  = mss_dst4(c,j) + qin_dst4(c)

             if (j <= -1) then
                ! No runoff over snow surface, just ponding on surface
                if (eff_porosity(c,j) < wimp .OR. eff_porosity(c,j+1) < wimp) then
                   qout(c) = 0._r8
                else
                   ! dz must be scaled by frac_sno to obtain gridcell average value
                   qout(c) = max(0._r8,(vol_liq(c,j) &
                        - ssi*eff_porosity(c,j))*dz(c,j)*frac_sno_eff(c))
                   qout(c) = min(qout(c),(1._r8-vol_ice(c,j+1) &
                        - vol_liq(c,j+1))*dz(c,j+1)*frac_sno_eff(c))
                end if
             else
                qout(c) = max(0._r8,(vol_liq(c,j) &
                     - ssi*eff_porosity(c,j))*dz(c,j)*frac_sno_eff(c))
             end if
             qout(c) = qout(c)*1000._r8
             h2osoi_liq(c,j) = h2osoi_liq(c,j) - qout(c)
             qin(c) = qout(c)

             ! mass of ice+water: in extremely rare circumstances, this can
             ! be zero, even though there is a snow layer defined. In
             ! this case, set the mass to a very small value to
             ! prevent division by zero.
             mss_liqice = h2osoi_liq(c,j)+h2osoi_ice(c,j)
             if (mss_liqice < 1E-30_r8) then
                mss_liqice = 1E-30_r8
             endif
         
             ! BCPHI:
             ! 1. flux with meltwater:
             qout_bc_phi(c) = qout(c)*scvng_fct_mlt_bcphi*(mss_bcphi(c,j)/mss_liqice)
             if (qout_bc_phi(c) > mss_bcphi(c,j)) then
                qout_bc_phi(c) = mss_bcphi(c,j)
             endif
             mss_bcphi(c,j) = mss_bcphi(c,j) - qout_bc_phi(c)
             qin_bc_phi(c) = qout_bc_phi(c)

             ! BCPHO:
             ! 1. flux with meltwater:
             qout_bc_pho(c) = qout(c)*scvng_fct_mlt_bcpho*(mss_bcpho(c,j)/mss_liqice)
             if (qout_bc_pho(c) > mss_bcpho(c,j)) then
                qout_bc_pho(c) = mss_bcpho(c,j)
             endif
             mss_bcpho(c,j) = mss_bcpho(c,j) - qout_bc_pho(c)
             qin_bc_pho(c) = qout_bc_pho(c)

             ! OCPHI:
             ! 1. flux with meltwater:
             qout_oc_phi(c) = qout(c)*scvng_fct_mlt_ocphi*(mss_ocphi(c,j)/mss_liqice)
             if (qout_oc_phi(c) > mss_ocphi(c,j)) then
                qout_oc_phi(c) = mss_ocphi(c,j)
             endif
             mss_ocphi(c,j) = mss_ocphi(c,j) - qout_oc_phi(c)
             qin_oc_phi(c) = qout_oc_phi(c)

             ! OCPHO:
             ! 1. flux with meltwater:
             qout_oc_pho(c) = qout(c)*scvng_fct_mlt_ocpho*(mss_ocpho(c,j)/mss_liqice)
             if (qout_oc_pho(c) > mss_ocpho(c,j)) then
                qout_oc_pho(c) = mss_ocpho(c,j)
             endif
             mss_ocpho(c,j) = mss_ocpho(c,j) - qout_oc_pho(c)
             qin_oc_pho(c) = qout_oc_pho(c)

             ! DUST 1:
             ! 1. flux with meltwater:
             qout_dst1(c) = qout(c)*scvng_fct_mlt_dst1*(mss_dst1(c,j)/mss_liqice)
             if (qout_dst1(c) > mss_dst1(c,j)) then
                qout_dst1(c) = mss_dst1(c,j)
             endif
             mss_dst1(c,j) = mss_dst1(c,j) - qout_dst1(c)
             qin_dst1(c) = qout_dst1(c)

             ! DUST 2:
             ! 1. flux with meltwater:
             qout_dst2(c) = qout(c)*scvng_fct_mlt_dst2*(mss_dst2(c,j)/mss_liqice)
             if (qout_dst2(c) > mss_dst2(c,j)) then
                qout_dst2(c) = mss_dst2(c,j)
             endif
             mss_dst2(c,j) = mss_dst2(c,j) - qout_dst2(c)
             qin_dst2(c) = qout_dst2(c)

             ! DUST 3:
             ! 1. flux with meltwater:
             qout_dst3(c) = qout(c)*scvng_fct_mlt_dst3*(mss_dst3(c,j)/mss_liqice)
             if (qout_dst3(c) > mss_dst3(c,j)) then
                qout_dst3(c) = mss_dst3(c,j)
             endif
             mss_dst3(c,j) = mss_dst3(c,j) - qout_dst3(c)
             qin_dst3(c) = qout_dst3(c)

             ! DUST 4:
             ! 1. flux with meltwater:
             qout_dst4(c) = qout(c)*scvng_fct_mlt_dst4*(mss_dst4(c,j)/mss_liqice)
             if (qout_dst4(c) > mss_dst4(c,j)) then
                qout_dst4(c) = mss_dst4(c,j)
             endif
             mss_dst4(c,j) = mss_dst4(c,j) - qout_dst4(c)
             qin_dst4(c) = qout_dst4(c)
             
          end if
       end do
    end do

    ! Adjust layer thickness for any water+ice content changes in excess of previous 
    ! layer thickness. Strictly speaking, only necessary for top snow layer, but doing
    ! it for all snow layers will catch problems with older initial files.
    ! Layer interfaces (zi) and node depths (z) do not need adjustment here because they
    ! are adjusted in CombineSnowLayers and are not used up to that point.

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
              dz(c,j) = max(dz(c,j),h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)
          end if
       end do
    end do
 
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       ! Qout from snow bottom
       qflx_snow_melt(c) = qflx_snow_melt(c) + (qout(c) / dtime)

       qflx_top_soil(c) = (qout(c) / dtime) &
            + (1.0_r8 - frac_sno_eff(c)) * qflx_rain_grnd(c)
       int_snow(c) = int_snow(c) + frac_sno_eff(c) * qflx_dew_snow(c) * dtime &
            + frac_sno_eff(c) * qflx_rain_grnd(c) * dtime
    end do

    do fc = 1, num_nosnowc
       c = filter_nosnowc(fc)
       qflx_snow_melt(c) = qflx_snomelt(c)

       qflx_top_soil(c) = qflx_rain_grnd(c) + qflx_snomelt(c)
       ! reset accumulated snow when no snow present
       if (h2osno(c) <= 0) int_snow(c) = 0.
       if (h2osno(c) <= 0) frac_sno(c) = 0.
    end do

    
    !  set aerosol deposition fluxes from forcing array
    !  The forcing array is either set from an external file 
    !  or from fluxes received from the atmosphere model
    do c = lbc,ubc
       g = cgridcell(c)
       
       flx_bc_dep_dry(c)   = forc_aer(g,1) + forc_aer(g,2)
       flx_bc_dep_wet(c)   = forc_aer(g,3)
       flx_bc_dep_phi(c)   = forc_aer(g,1) + forc_aer(g,3)
       flx_bc_dep_pho(c)   = forc_aer(g,2)
       flx_bc_dep(c)       = forc_aer(g,1) + forc_aer(g,2) + forc_aer(g,3)
       
       flx_oc_dep_dry(c)   = forc_aer(g,4) + forc_aer(g,5)
       flx_oc_dep_wet(c)   = forc_aer(g,6)
       flx_oc_dep_phi(c)   = forc_aer(g,4) + forc_aer(g,6)
       flx_oc_dep_pho(c)   = forc_aer(g,5)
       flx_oc_dep(c)       = forc_aer(g,4) + forc_aer(g,5) + forc_aer(g,6)
       
       flx_dst_dep_wet1(c) = forc_aer(g,7)
       flx_dst_dep_dry1(c) = forc_aer(g,8)
       flx_dst_dep_wet2(c) = forc_aer(g,9)
       flx_dst_dep_dry2(c) = forc_aer(g,10)
       flx_dst_dep_wet3(c) = forc_aer(g,11)
       flx_dst_dep_dry3(c) = forc_aer(g,12)
       flx_dst_dep_wet4(c) = forc_aer(g,13)
       flx_dst_dep_dry4(c) = forc_aer(g,14)
       flx_dst_dep(c)      = forc_aer(g,7) + forc_aer(g,8) + forc_aer(g,9) + &
                             forc_aer(g,10) + forc_aer(g,11) + forc_aer(g,12) + &
                             forc_aer(g,13) + forc_aer(g,14)
    
    end do

    ! aerosol deposition fluxes into top layer
    ! This is done after the inter-layer fluxes so that some aerosol
    ! is in the top layer after deposition, and is not immediately
    ! washed out before radiative calculations are done
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       mss_bcphi(c,snl(c)+1) = mss_bcphi(c,snl(c)+1) + (flx_bc_dep_phi(c)*dtime)
       mss_bcpho(c,snl(c)+1) = mss_bcpho(c,snl(c)+1) + (flx_bc_dep_pho(c)*dtime)
       mss_ocphi(c,snl(c)+1) = mss_ocphi(c,snl(c)+1) + (flx_oc_dep_phi(c)*dtime)
       mss_ocpho(c,snl(c)+1) = mss_ocpho(c,snl(c)+1) + (flx_oc_dep_pho(c)*dtime)
       
       mss_dst1(c,snl(c)+1) = mss_dst1(c,snl(c)+1) + (flx_dst_dep_dry1(c) + flx_dst_dep_wet1(c))*dtime
       mss_dst2(c,snl(c)+1) = mss_dst2(c,snl(c)+1) + (flx_dst_dep_dry2(c) + flx_dst_dep_wet2(c))*dtime
       mss_dst3(c,snl(c)+1) = mss_dst3(c,snl(c)+1) + (flx_dst_dep_dry3(c) + flx_dst_dep_wet3(c))*dtime
       mss_dst4(c,snl(c)+1) = mss_dst4(c,snl(c)+1) + (flx_dst_dep_dry4(c) + flx_dst_dep_wet4(c))*dtime
    end do

  end subroutine SnowWater

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SnowCompaction
!
! !INTERFACE:
  subroutine SnowCompaction(lbc, ubc, num_snowc, filter_snowc)
!
! !DESCRIPTION:
! Determine the change in snow layer thickness due to compaction and
! settling.
! Three metamorphisms of changing snow characteristics are implemented,
! i.e., destructive, overburden, and melt. The treatments of the former
! two are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution
! due to melt metamorphism is simply taken as a ratio of snow ice
! fraction after the melting versus before the melting.
!
! !USES:
    use clmtype
    use clm_time_manager, only : get_step_size
    use clm_varcon      , only : denice, denh2o, tfrz, istice_mec
    use clm_varcon      , only : rpi, isturb, istdlak, istsoil, istcrop
    use clm_varctl      , only : subgridflag
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                ! column bounds
    integer, intent(in) :: num_snowc               ! number of column snow points in column filter
    integer, intent(in) :: filter_snowc(ubc-lbc+1) ! column filter for snow points
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/28/02, Peter Thornton: Migrated to new data structures
! 2/29/08, David Lawrence: Revised snow overburden to be include 0.5 weight of current layer
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
    real(r8), pointer :: frac_sno(:)       !snow covered fraction
    real(r8), pointer :: swe_old(:,:)      !initial swe values
    real(r8), pointer :: int_snow(:)       !integrated snowfall [mm]
    real(r8), pointer :: n_melt(:)         !SCA shape parameter
    integer,  pointer :: snl(:)            !number of snow layers
!
! local pointers to implicit in arguments
!
    integer,  pointer :: imelt(:,:)        !flag for melting (=1), freezing (=2), Not=0
    real(r8), pointer :: frac_iceold(:,:)  !fraction of ice relative to the tot water
    real(r8), pointer :: t_soisno(:,:)     !soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   !liquid water (kg/m2)
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: dz(:,:)           !layer depth (m)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer :: j, l, c, fc                   ! indices
    real(r8):: dtime                         ! land model time step (sec)
    real(r8), parameter :: c2 = 23.e-3_r8    ! [m3/kg]
    real(r8), parameter :: c3 = 2.777e-6_r8  ! [1/s]
    real(r8), parameter :: c4 = 0.04_r8      ! [1/K]
    real(r8), parameter :: c5 = 2.0_r8       !
    real(r8), parameter :: dm = 100.0_r8     ! Upper Limit on Destructive Metamorphism Compaction [kg/m3]
    real(r8), parameter :: eta0 = 9.e+5_r8   ! The Viscosity Coefficient Eta0 [kg-s/m2]
    real(r8) :: burden(lbc:ubc) ! pressure of overlying snow [kg/m2]
    real(r8) :: ddz1   ! Rate of settling of snowpack due to destructive metamorphism.
    real(r8) :: ddz2   ! Rate of compaction of snowpack due to overburden.
    real(r8) :: ddz3   ! Rate of compaction of snowpack due to melt [1/s]
    real(r8) :: dexpf  ! expf=exp(-c4*(273.15-t_soisno)).
    real(r8) :: fi     ! Fraction of ice relative to the total water content at current time step
    real(r8) :: td     ! t_soisno - tfrz [K]
    real(r8) :: pdzdtc ! Nodal rate of change in fractional-thickness due to compaction [fraction/s]
    real(r8) :: void   ! void (1 - vol_ice - vol_liq)
    real(r8) :: wx     ! water mass (ice+liquid) [kg/m2]
    real(r8) :: bi     ! partial density of ice [kg/m3]
    real(r8) :: wsum   ! snowpack total water mass (ice+liquid) [kg/m2]
    real(r8) :: fsno_melt
    integer, pointer :: clandunit(:)       !landunit index for each column
    integer, pointer :: ltype(:)           !landunit type
    real(r8), pointer :: snow_depth(:)         ! snow height (m)

!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (column-level)

    snow_depth      => clm3%g%l%c%cps%snow_depth
    frac_sno    => clm3%g%l%c%cps%frac_sno_eff 
    swe_old     => clm3%g%l%c%cws%swe_old
    int_snow    => clm3%g%l%c%cws%int_snow
    ltype       => clm3%g%l%itype
    clandunit   => clm3%g%l%c%landunit
    n_melt      => clm3%g%l%c%cps%n_melt
    snl         => clm3%g%l%c%cps%snl
    dz          => clm3%g%l%c%cps%dz
    imelt       => clm3%g%l%c%cps%imelt
    frac_iceold => clm3%g%l%c%cps%frac_iceold
    t_soisno    => clm3%g%l%c%ces%t_soisno
    h2osoi_ice  => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq  => clm3%g%l%c%cws%h2osoi_liq
    clandunit   => clm3%g%l%c%landunit
    ltype       => clm3%g%l%itype

    ! Get time step

    dtime = get_step_size()

    ! Begin calculation - note that the following column loops are only invoked if snl(c) < 0

    burden(:) = 0._r8

    do j = -nlevsno+1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then

             wx = h2osoi_ice(c,j) + h2osoi_liq(c,j)
             void = 1._r8 - (h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o) / dz(c,j)
             wx = (h2osoi_ice(c,j) + h2osoi_liq(c,j))
             void = 1._r8 - (h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o)&
                  /(frac_sno(c) * dz(c,j))
             ! If void is negative, then increase dz such that void = 0.
             ! This should be done for any landunit, but for now is done only for glacier_mec 1andunits.
             l = clandunit(c)
             if (ltype(l)==istice_mec .and. void < 0._r8) then
                dz(c,j) = h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o
                void = 0._r8
             endif

             ! Allow compaction only for non-saturated node and higher ice lens node.
             if (void > 0.001_r8 .and. h2osoi_ice(c,j) > .1_r8) then

                bi = h2osoi_ice(c,j) / (frac_sno(c) * dz(c,j))
                fi = h2osoi_ice(c,j) / wx
                td = tfrz-t_soisno(c,j)
                dexpf = exp(-c4*td)

                ! Settling as a result of destructive metamorphism

                ddz1 = -c3*dexpf
                if (bi > dm) ddz1 = ddz1*exp(-46.0e-3_r8*(bi-dm))

                ! Liquid water term

                if (h2osoi_liq(c,j) > 0.01_r8*dz(c,j)*frac_sno(c)) ddz1=ddz1*c5

                ! Compaction due to overburden

                ddz2 = -(burden(c)+wx/2._r8)*exp(-0.08_r8*td - c2*bi)/eta0 

                ! Compaction occurring during melt

                if (imelt(c,j) == 1) then
                    if(subgridflag==1 .and. (ltype(clandunit(c)) == istsoil .or. ltype(clandunit(c)) == istcrop)) then
                      ! first term is delta mass over mass
                      ddz3 = max(0._r8,min(1._r8,(swe_old(c,j) - wx)/wx))

                      ! 2nd term is delta fsno over fsno, allowing for negative values for ddz3
                      wsum = sum(h2osoi_liq(c,snl(c)+1:0)+h2osoi_ice(c,snl(c)+1:0))
                      fsno_melt = 1. - (acos(2.*min(1._r8,wsum/int_snow(c)) - 1._r8)/rpi)**(n_melt(c))

                      ddz3 = ddz3 - max(0._r8,(fsno_melt - frac_sno(c))/frac_sno(c))
                      ddz3 = -1._r8/dtime * ddz3
                   else
                      ddz3 = - 1._r8/dtime * max(0._r8,(frac_iceold(c,j) - fi)/frac_iceold(c,j))
                   endif
                else
                   ddz3 = 0._r8
                end if

                ! Time rate of fractional change in dz (units of s-1)

                pdzdtc = ddz1 + ddz2 + ddz3

                ! The change in dz due to compaction
                ! Limit compaction to be no greater than fully saturated layer thickness

                dz(c,j) = max(dz(c,j) * (1._r8+pdzdtc*dtime),(h2osoi_ice(c,j)/denice+ h2osoi_liq(c,j)/denh2o)/frac_sno(c))
             end if

             ! Pressure of overlying snow

             burden(c) = burden(c) + wx

          end if
       end do
    end do

  end subroutine SnowCompaction

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CombineSnowLayers
!
! !INTERFACE:
  subroutine CombineSnowLayers(lbc, ubc, num_snowc, filter_snowc)
!
! !DESCRIPTION:
! Combine snow layers that are less than a minimum thickness or mass
! If the snow element thickness or mass is less than a prescribed minimum,
! then it is combined with a neighboring element.  The subroutine
! clm\_combo.f90 then executes the combination of mass and energy.
!
! !USES:
    use clmtype
    use clm_varcon, only : istsoil, isturb, istdlak
    use SLakeCon  , only : lsadz
    use clm_varcon, only : istsoil, isturb,istwet,istice, istice_mec
    use clm_varcon, only : istcrop
    use clm_time_manager, only : get_step_size
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)    :: lbc, ubc                    ! column bounds
    integer, intent(inout) :: num_snowc                   ! number of column snow points in column filter
    integer, intent(inout) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/28/02, Peter Thornton: Migrated to new data structures.
! 03/28/08, Mark Flanner: Added aerosol masses and snow grain radius
! 05/20/10, Zack Subin: Adjust minimum thickness for snow layers over lakes to be + lsadz
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: frac_sno(:)       !fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_sno_eff(:)   !fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: int_snow(:)       !integrated snowfall [mm]
    integer, pointer :: clandunit(:)       !landunit index for each column
    integer, pointer :: ltype(:)           !landunit type
!
! local pointers to implicit inout arguments
!
    integer , pointer :: snl(:)            !number of snow layers
    real(r8), pointer :: h2osno(:)         !snow water (mm H2O)
    real(r8), pointer :: snow_depth(:)         !snow height (m)
    real(r8), pointer :: dz(:,:)           !layer depth (m)
    real(r8), pointer :: zi(:,:)           !interface level below a "z" level (m)
    real(r8), pointer :: t_soisno(:,:)     !soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   !liquid water (kg/m2)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: z(:,:)          ! layer thickness (m)
    real(r8), pointer :: mss_bcphi(:,:)  ! hydrophilic BC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bcpho(:,:)  ! hydrophobic BC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocphi(:,:)  ! hydrophilic OC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocpho(:,:)  ! hydrophobic OC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst1(:,:)   ! dust species 1 mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst2(:,:)   ! dust species 2 mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst3(:,:)   ! dust species 3 mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst4(:,:)   ! dust species 4 mass in snow (col,lyr) [kg]
    real(r8), pointer :: snw_rds(:,:)    ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(r8), pointer :: qflx_sl_top_soil(:) ! liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer :: c, fc                 ! column indices
    integer :: i,k                   ! loop indices
    integer :: j,l                   ! node indices
    integer :: msn_old(lbc:ubc)      ! number of top snow layer
    integer :: mssi(lbc:ubc)         ! node index
    integer :: neibor                ! adjacent node selected for combination
    real(r8):: zwice(lbc:ubc)        ! total ice mass in snow
    real(r8):: zwliq (lbc:ubc)       ! total liquid water in snow
    real(r8):: dzmin(5)              ! minimum of top snow layer
    real(r8):: dzminloc(5)           ! minimum of top snow layer (local)
    real(r8) :: dtime                !land model time step (sec)

    data dzmin /0.010_r8, 0.015_r8, 0.025_r8, 0.055_r8, 0.115_r8/
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes (landunit-level)

    ltype      => clm3%g%l%itype

    ! Assign local pointers to derived subtypes (column-level)

    frac_sno       => clm3%g%l%c%cps%frac_sno
    frac_sno_eff   => clm3%g%l%c%cps%frac_sno_eff 
    int_snow       => clm3%g%l%c%cws%int_snow
    clandunit  => clm3%g%l%c%landunit
    snl        => clm3%g%l%c%cps%snl
    snow_depth     => clm3%g%l%c%cps%snow_depth
    h2osno     => clm3%g%l%c%cws%h2osno
    dz         => clm3%g%l%c%cps%dz
    zi         => clm3%g%l%c%cps%zi
    z          => clm3%g%l%c%cps%z
    t_soisno   => clm3%g%l%c%ces%t_soisno
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq
    mss_bcphi  => clm3%g%l%c%cps%mss_bcphi
    mss_bcpho  => clm3%g%l%c%cps%mss_bcpho
    mss_ocphi  => clm3%g%l%c%cps%mss_ocphi
    mss_ocpho  => clm3%g%l%c%cps%mss_ocpho
    mss_dst1   => clm3%g%l%c%cps%mss_dst1
    mss_dst2   => clm3%g%l%c%cps%mss_dst2
    mss_dst3   => clm3%g%l%c%cps%mss_dst3
    mss_dst4   => clm3%g%l%c%cps%mss_dst4
    snw_rds    => clm3%g%l%c%cps%snw_rds
    qflx_sl_top_soil => clm3%g%l%c%cwf%qflx_sl_top_soil

    ! Determine model time step

    dtime = get_step_size()


    ! Check the mass of ice lens of snow, when the total is less than a small value,
    ! combine it with the underlying neighbor.

    dzminloc(:) = dzmin(:) ! dzmin will stay constant between timesteps

    ! Add lsadz to dzmin for lakes
    ! Determine whether called from SLakeHydrology
    ! Note: this assumes that this function is called separately with the lake-snow and non-lake-snow filters.
    if (num_snowc > 0) then
       c = filter_snowc(1)
       l = clandunit(c)
       if (ltype(l) == istdlak) then ! Called from SLakeHydrology
          dzminloc(:) = dzmin(:) + lsadz
       end if
    end if

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       msn_old(c) = snl(c)
       qflx_sl_top_soil(c) = 0._r8
    end do

    ! The following loop is NOT VECTORIZED

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       l = clandunit(c)
       do j = msn_old(c)+1,0
          ! use 0.01 to avoid runaway ice buildup
          if (h2osoi_ice(c,j) <= .01_r8) then
             if (ltype(l) == istsoil .or. ltype(l)==isturb .or. ltype(l) == istcrop) then
                h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
                h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)

                if (j == 0) then
                  qflx_sl_top_soil(c) = (h2osoi_liq(c,j) + h2osoi_ice(c,j))/dtime
                end if

                if (j /= 0) dz(c,j+1) = dz(c,j+1) + dz(c,j)
                
                ! NOTE: Temperature, and similarly snw_rds, of the
                ! underlying snow layer are NOT adjusted in this case. 
                ! Because the layer being eliminated has a small mass, 
                ! this should not make a large difference, but it 
                ! would be more thorough to do so.
                if (j /= 0) then
                   mss_bcphi(c,j+1) = mss_bcphi(c,j+1) + mss_bcphi(c,j)
                   mss_bcpho(c,j+1) = mss_bcpho(c,j+1) + mss_bcpho(c,j)
                   mss_ocphi(c,j+1) = mss_ocphi(c,j+1) + mss_ocphi(c,j)
                   mss_ocpho(c,j+1) = mss_ocpho(c,j+1) + mss_ocpho(c,j)
                   mss_dst1(c,j+1)  = mss_dst1(c,j+1) + mss_dst1(c,j)
                   mss_dst2(c,j+1)  = mss_dst2(c,j+1) + mss_dst2(c,j)
                   mss_dst3(c,j+1)  = mss_dst3(c,j+1) + mss_dst3(c,j)
                   mss_dst4(c,j+1)  = mss_dst4(c,j+1) + mss_dst4(c,j)
                end if

             else if (ltype(l) /= istsoil .and. ltype(l) /= isturb .and. ltype(l) /= istcrop .and. j /= 0) then
                h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
                h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)
                dz(c,j+1) = dz(c,j+1) + dz(c,j)
                
                mss_bcphi(c,j+1) = mss_bcphi(c,j+1) + mss_bcphi(c,j)
                mss_bcpho(c,j+1) = mss_bcpho(c,j+1) + mss_bcpho(c,j)
                mss_ocphi(c,j+1) = mss_ocphi(c,j+1) + mss_ocphi(c,j)
                mss_ocpho(c,j+1) = mss_ocpho(c,j+1) + mss_ocpho(c,j)
                mss_dst1(c,j+1)  = mss_dst1(c,j+1) + mss_dst1(c,j)
                mss_dst2(c,j+1)  = mss_dst2(c,j+1) + mss_dst2(c,j)
                mss_dst3(c,j+1)  = mss_dst3(c,j+1) + mss_dst3(c,j)
                mss_dst4(c,j+1)  = mss_dst4(c,j+1) + mss_dst4(c,j)

             end if

             ! shift all elements above this down one.
             if (j > snl(c)+1 .and. snl(c) < -1) then
                do i = j, snl(c)+2, -1
                   ! If the layer closest to the surface is less than 0.1 mm and the ltype is not
                   ! urban, soil or crop, the h2osoi_liq and h2osoi_ice associated with this layer is sent 
                   ! to qflx_qrgwl later on in the code.  To keep track of this for the snow balance
                   ! error check, we add this to qflx_sl_top_soil here
                   if (ltype(l) /= istsoil .and. ltype(l) /= istcrop .and. ltype(l) /= isturb .and. i == 0) then
                     qflx_sl_top_soil(c) = (h2osoi_liq(c,i) + h2osoi_ice(c,i))/dtime
                   end if

                   t_soisno(c,i)   = t_soisno(c,i-1)
                   h2osoi_liq(c,i) = h2osoi_liq(c,i-1)
                   h2osoi_ice(c,i) = h2osoi_ice(c,i-1)
                   
                   mss_bcphi(c,i)   = mss_bcphi(c,i-1)
                   mss_bcpho(c,i)   = mss_bcpho(c,i-1)
                   mss_ocphi(c,i)   = mss_ocphi(c,i-1)
                   mss_ocpho(c,i)   = mss_ocpho(c,i-1)
                   mss_dst1(c,i)    = mss_dst1(c,i-1)
                   mss_dst2(c,i)    = mss_dst2(c,i-1)
                   mss_dst3(c,i)    = mss_dst3(c,i-1)
                   mss_dst4(c,i)    = mss_dst4(c,i-1)
                   snw_rds(c,i)     = snw_rds(c,i-1)

                   dz(c,i)         = dz(c,i-1)
                end do
             end if
             snl(c) = snl(c) + 1
          end if
       end do
    end do

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       h2osno(c) = 0._r8
       snow_depth(c) = 0._r8
       zwice(c)  = 0._r8
       zwliq(c)  = 0._r8
    end do

    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             h2osno(c) = h2osno(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
             snow_depth(c) = snow_depth(c) + dz(c,j)
             zwice(c)  = zwice(c) + h2osoi_ice(c,j)
             zwliq(c)  = zwliq(c) + h2osoi_liq(c,j)
          end if
       end do
    end do

    ! Check the snow depth - all snow gone
    ! The liquid water assumes ponding on soil surface.

    do fc = 1, num_snowc
       c = filter_snowc(fc)
       l = clandunit(c)
       if (snow_depth(c) > 0._r8) then
          if ((ltype(l) == istdlak .and. snow_depth(c) < 0.01_r8 + lsadz ) .or. &
               ((ltype(l) /= istdlak) .and. ((frac_sno_eff(c)*snow_depth(c) < 0.01_r8)  &
               .or. (h2osno(c)/(frac_sno_eff(c)*snow_depth(c)) < 50._r8)))) then

          snl(c) = 0
          h2osno(c) = zwice(c)

          mss_bcphi(c,:) = 0._r8
          mss_bcpho(c,:) = 0._r8
          mss_ocphi(c,:) = 0._r8
          mss_ocpho(c,:) = 0._r8
          mss_dst1(c,:)  = 0._r8
          mss_dst2(c,:)  = 0._r8
          mss_dst3(c,:)  = 0._r8
          mss_dst4(c,:)  = 0._r8

          if (h2osno(c) <= 0._r8) snow_depth(c) = 0._r8
             ! this is where water is transfered from layer 0 (snow) to layer 1 (soil)
             if (ltype(l) == istsoil .or. ltype(l) == isturb .or. ltype(l) == istcrop) then
                h2osoi_liq(c,0) = 0.0_r8
             h2osoi_liq(c,1) = h2osoi_liq(c,1) + zwliq(c)
          end if
             if (ltype(l) == istwet) then             
                h2osoi_liq(c,0) = 0.0_r8
             endif
             if (ltype(l) == istice .or. ltype(l)==istice_mec) then             
                h2osoi_liq(c,0) = 0.0_r8
             endif
          endif
       end if
          if (h2osno(c) <= 0._r8) then
             snow_depth(c) = 0._r8
             frac_sno(c) = 0._r8
             frac_sno_eff(c) = 0._r8
             int_snow(c) = 0._r8
          endif
    end do

    ! Check the snow depth - snow layers combined
    ! The following loop IS NOT VECTORIZED

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       ! Two or more layers

       if (snl(c) < -1) then

          msn_old(c) = snl(c)
          mssi(c) = 1

          do i = msn_old(c)+1,0
             if ((frac_sno_eff(c)*dz(c,i) < dzminloc(mssi(c))) .or. &
                  ((h2osoi_ice(c,i) + h2osoi_liq(c,i))/(frac_sno_eff(c)*dz(c,i)) < 50._r8)) then
                if (i == snl(c)+1) then
                   ! If top node is removed, combine with bottom neighbor.
                   neibor = i + 1
                else if (i == 0) then
                   ! If the bottom neighbor is not snow, combine with the top neighbor.
                   neibor = i - 1
                else
                   ! If none of the above special cases apply, combine with the thinnest neighbor
                   neibor = i + 1
                   if ((dz(c,i-1)+dz(c,i)) < (dz(c,i+1)+dz(c,i))) neibor = i-1
                end if

                ! Node l and j are combined and stored as node j.
                if (neibor > i) then
                   j = neibor
                   l = i
                else
                   j = i
                   l = neibor
                end if

                ! this should be included in 'Combo' for consistency,
                ! but functionally it is the same to do it here
                mss_bcphi(c,j)=mss_bcphi(c,j)+mss_bcphi(c,l)
                mss_bcpho(c,j)=mss_bcpho(c,j)+mss_bcpho(c,l)
                mss_ocphi(c,j)=mss_ocphi(c,j)+mss_ocphi(c,l)
                mss_ocpho(c,j)=mss_ocpho(c,j)+mss_ocpho(c,l)
                mss_dst1(c,j)=mss_dst1(c,j)+mss_dst1(c,l)
                mss_dst2(c,j)=mss_dst2(c,j)+mss_dst2(c,l)
                mss_dst3(c,j)=mss_dst3(c,j)+mss_dst3(c,l)
                mss_dst4(c,j)=mss_dst4(c,j)+mss_dst4(c,l)
                ! mass-weighted combination of effective grain size:
                snw_rds(c,j) = (snw_rds(c,j)*(h2osoi_liq(c,j)+h2osoi_ice(c,j)) + &
                               snw_rds(c,l)*(h2osoi_liq(c,l)+h2osoi_ice(c,l))) / &
                               (h2osoi_liq(c,j)+h2osoi_ice(c,j)+h2osoi_liq(c,l)+h2osoi_ice(c,l))

                call Combo (dz(c,j), h2osoi_liq(c,j), h2osoi_ice(c,j), &
                   t_soisno(c,j), dz(c,l), h2osoi_liq(c,l), h2osoi_ice(c,l), t_soisno(c,l) )

                ! Now shift all elements above this down one.
                if (j-1 > snl(c)+1) then

                   do k = j-1, snl(c)+2, -1
                      t_soisno(c,k) = t_soisno(c,k-1)
                      h2osoi_ice(c,k) = h2osoi_ice(c,k-1)
                      h2osoi_liq(c,k) = h2osoi_liq(c,k-1)

                      mss_bcphi(c,k) = mss_bcphi(c,k-1)
                      mss_bcpho(c,k) = mss_bcpho(c,k-1)
                      mss_ocphi(c,k) = mss_ocphi(c,k-1)
                      mss_ocpho(c,k) = mss_ocpho(c,k-1)
                      mss_dst1(c,k)  = mss_dst1(c,k-1)
                      mss_dst2(c,k)  = mss_dst2(c,k-1)
                      mss_dst3(c,k)  = mss_dst3(c,k-1)
                      mss_dst4(c,k)  = mss_dst4(c,k-1)
                      snw_rds(c,k)   = snw_rds(c,k-1)

                      dz(c,k) = dz(c,k-1)
                   end do
                end if

                ! Decrease the number of snow layers
                snl(c) = snl(c) + 1
                if (snl(c) >= -1) EXIT

             else

                ! The layer thickness is greater than the prescribed minimum value
                mssi(c) = mssi(c) + 1

             end if
          end do

       end if

    end do

    ! Reset the node depth and the depth of layer interface

    do j = 0, -nlevsno+1, -1
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c) + 1) then
             z(c,j) = zi(c,j) - 0.5_r8*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

  end subroutine CombineSnowLayers

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DivideSnowLayers
!
! !INTERFACE:
  subroutine DivideSnowLayers(lbc, ubc, num_snowc, filter_snowc)
!
! !DESCRIPTION:
! Subdivides snow layers if they exceed their prescribed maximum thickness.
!
! !USES:
    use clmtype
    use clm_varcon,  only : tfrz 
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)    :: lbc, ubc                    ! column bounds
    integer, intent(inout) :: num_snowc                   ! number of column snow points in column filter
    integer, intent(inout) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/28/02, Peter Thornton: Migrated to new data structures.
! 2/29/08, David Lawrence: Snowpack T profile maintained during layer splitting
! 03/28/08, Mark Flanner: Added aerosol masses and snow grain radius
!
! !LOCAL VARIABLES:
!
! local pointers to implicit inout arguments
!
    real(r8), pointer :: frac_sno(:)       !fraction of ground covered by snow (0 to 1)
    integer , pointer :: snl(:)            !number of snow layers
    real(r8), pointer :: dz(:,:)           !layer depth (m)
    real(r8), pointer :: zi(:,:)           !interface level below a "z" level (m)
    real(r8), pointer :: t_soisno(:,:)     !soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   !liquid water (kg/m2)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: z(:,:)          ! layer thickness (m)
    real(r8), pointer :: mss_bcphi(:,:)  ! hydrophilic BC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bcpho(:,:)  ! hydrophobic BC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocphi(:,:)  ! hydrophilic OC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocpho(:,:)  ! hydrophobic OC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst1(:,:)   ! dust species 1 mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst2(:,:)   ! dust species 2 mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst3(:,:)   ! dust species 3 mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst4(:,:)   ! dust species 4 mass in snow (col,lyr) [kg]
    real(r8), pointer :: snw_rds(:,:)    ! effective snow grain radius (col,lyr) [microns, m^-6]
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: j, c, fc               ! indices
    real(r8) :: drr                    ! thickness of the combined [m]
    integer  :: msno                   ! number of snow layer 1 (top) to msno (bottom)
    real(r8) :: dzsno(lbc:ubc,nlevsno) ! Snow layer thickness [m]
    real(r8) :: swice(lbc:ubc,nlevsno) ! Partial volume of ice [m3/m3]
    real(r8) :: swliq(lbc:ubc,nlevsno) ! Partial volume of liquid water [m3/m3]
    real(r8) :: tsno(lbc:ubc ,nlevsno) ! Nodel temperature [K]
    real(r8) :: zwice                  ! temporary
    real(r8) :: zwliq                  ! temporary
    real(r8) :: propor                 ! temporary
    real(r8) :: dtdz                   ! temporary

    ! temporary variables mimicking the structure of other layer division variables
    real(r8) :: mbc_phi(lbc:ubc,nlevsno) ! mass of BC in each snow layer
    real(r8) :: zmbc_phi                 ! temporary
    real(r8) :: mbc_pho(lbc:ubc,nlevsno) ! mass of BC in each snow layer
    real(r8) :: zmbc_pho                 ! temporary
    real(r8) :: moc_phi(lbc:ubc,nlevsno) ! mass of OC in each snow layer
    real(r8) :: zmoc_phi                 ! temporary
    real(r8) :: moc_pho(lbc:ubc,nlevsno) ! mass of OC in each snow layer
    real(r8) :: zmoc_pho                 ! temporary
    real(r8) :: mdst1(lbc:ubc,nlevsno)   ! mass of dust 1 in each snow layer
    real(r8) :: zmdst1                   ! temporary
    real(r8) :: mdst2(lbc:ubc,nlevsno)   ! mass of dust 2 in each snow layer
    real(r8) :: zmdst2                   ! temporary
    real(r8) :: mdst3(lbc:ubc,nlevsno)   ! mass of dust 3 in each snow layer
    real(r8) :: zmdst3                   ! temporary
    real(r8) :: mdst4(lbc:ubc,nlevsno)   ! mass of dust 4 in each snow layer
    real(r8) :: zmdst4                   ! temporary
    real(r8) :: rds(lbc:ubc,nlevsno)

!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtype components (column-level)

    frac_sno   => clm3%g%l%c%cps%frac_sno_eff 
    snl        => clm3%g%l%c%cps%snl
    dz         => clm3%g%l%c%cps%dz
    zi         => clm3%g%l%c%cps%zi
    z          => clm3%g%l%c%cps%z
    t_soisno   => clm3%g%l%c%ces%t_soisno
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq
    mss_bcphi  => clm3%g%l%c%cps%mss_bcphi
    mss_bcpho  => clm3%g%l%c%cps%mss_bcpho
    mss_ocphi  => clm3%g%l%c%cps%mss_ocphi
    mss_ocpho  => clm3%g%l%c%cps%mss_ocpho
    mss_dst1   => clm3%g%l%c%cps%mss_dst1
    mss_dst2   => clm3%g%l%c%cps%mss_dst2
    mss_dst3   => clm3%g%l%c%cps%mss_dst3
    mss_dst4   => clm3%g%l%c%cps%mss_dst4
    snw_rds    => clm3%g%l%c%cps%snw_rds
    

    ! Begin calculation - note that the following column loops are only invoked
    ! for snow-covered columns

    do j = 1,nlevsno
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j <= abs(snl(c))) then
             dzsno(c,j) = frac_sno(c)*dz(c,j+snl(c))
             swice(c,j) = h2osoi_ice(c,j+snl(c))
             swliq(c,j) = h2osoi_liq(c,j+snl(c))
             tsno(c,j)  = t_soisno(c,j+snl(c))

             mbc_phi(c,j) = mss_bcphi(c,j+snl(c))
             mbc_pho(c,j) = mss_bcpho(c,j+snl(c))
             moc_phi(c,j) = mss_ocphi(c,j+snl(c))
             moc_pho(c,j) = mss_ocpho(c,j+snl(c))
             mdst1(c,j)   = mss_dst1(c,j+snl(c))
             mdst2(c,j)   = mss_dst2(c,j+snl(c))
             mdst3(c,j)   = mss_dst3(c,j+snl(c))
             mdst4(c,j)   = mss_dst4(c,j+snl(c))
             rds(c,j)     = snw_rds(c,j+snl(c))

          end if
       end do
    end do

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       msno = abs(snl(c))

       if (msno == 1) then
          ! Specify a new snow layer
          if (dzsno(c,1) > 0.03_r8) then
             msno = 2
             dzsno(c,1) = dzsno(c,1)/2._r8
             swice(c,1) = swice(c,1)/2._r8
             swliq(c,1) = swliq(c,1)/2._r8
             dzsno(c,2) = dzsno(c,1)
             swice(c,2) = swice(c,1)
             swliq(c,2) = swliq(c,1)
             tsno(c,2)  = tsno(c,1)
             
             mbc_phi(c,1) = mbc_phi(c,1)/2._r8
             mbc_phi(c,2) = mbc_phi(c,1)
             mbc_pho(c,1) = mbc_pho(c,1)/2._r8
             mbc_pho(c,2) = mbc_pho(c,1)
             moc_phi(c,1) = moc_phi(c,1)/2._r8
             moc_phi(c,2) = moc_phi(c,1)
             moc_pho(c,1) = moc_pho(c,1)/2._r8
             moc_pho(c,2) = moc_pho(c,1)
             mdst1(c,1) = mdst1(c,1)/2._r8
             mdst1(c,2) = mdst1(c,1)
             mdst2(c,1) = mdst2(c,1)/2._r8
             mdst2(c,2) = mdst2(c,1)
             mdst3(c,1) = mdst3(c,1)/2._r8
             mdst3(c,2) = mdst3(c,1)
             mdst4(c,1) = mdst4(c,1)/2._r8
             mdst4(c,2) = mdst4(c,1)
             rds(c,2) = rds(c,1)

          end if
       end if

       if (msno > 1) then
          if (dzsno(c,1) > 0.02_r8) then
             drr = dzsno(c,1) - 0.02_r8
             propor = drr/dzsno(c,1)
             zwice = propor*swice(c,1)
             zwliq = propor*swliq(c,1)

             zmbc_phi = propor*mbc_phi(c,1)
             zmbc_pho = propor*mbc_pho(c,1)
             zmoc_phi = propor*moc_phi(c,1)
             zmoc_pho = propor*moc_pho(c,1)
             zmdst1 = propor*mdst1(c,1)
             zmdst2 = propor*mdst2(c,1)
             zmdst3 = propor*mdst3(c,1)
             zmdst4 = propor*mdst4(c,1)

             propor = 0.02_r8/dzsno(c,1)
             swice(c,1) = propor*swice(c,1)
             swliq(c,1) = propor*swliq(c,1)

             mbc_phi(c,1) = propor*mbc_phi(c,1)
             mbc_pho(c,1) = propor*mbc_pho(c,1)
             moc_phi(c,1) = propor*moc_phi(c,1)
             moc_pho(c,1) = propor*moc_pho(c,1)
             mdst1(c,1) = propor*mdst1(c,1)
             mdst2(c,1) = propor*mdst2(c,1)
             mdst3(c,1) = propor*mdst3(c,1)
             mdst4(c,1) = propor*mdst4(c,1)

             dzsno(c,1) = 0.02_r8

             mbc_phi(c,2) = mbc_phi(c,2)+zmbc_phi  ! (combo)
             mbc_pho(c,2) = mbc_pho(c,2)+zmbc_pho  ! (combo)
             moc_phi(c,2) = moc_phi(c,2)+zmoc_phi  ! (combo)
             moc_pho(c,2) = moc_pho(c,2)+zmoc_pho  ! (combo)
             mdst1(c,2) = mdst1(c,2)+zmdst1  ! (combo)
             mdst2(c,2) = mdst2(c,2)+zmdst2  ! (combo)
             mdst3(c,2) = mdst3(c,2)+zmdst3  ! (combo)
             mdst4(c,2) = mdst4(c,2)+zmdst4  ! (combo)
             rds(c,2) = rds(c,1) ! (combo)

             call Combo (dzsno(c,2), swliq(c,2), swice(c,2), tsno(c,2), drr, &
                  zwliq, zwice, tsno(c,1))

             ! Subdivide a new layer
             if (msno <= 2 .and. dzsno(c,2) > 0.07_r8) then
                msno = 3
                dtdz = (tsno(c,1) - tsno(c,2))/((dzsno(c,1)+dzsno(c,2))/2._r8) 
                dzsno(c,2) = dzsno(c,2)/2._r8
                swice(c,2) = swice(c,2)/2._r8
                swliq(c,2) = swliq(c,2)/2._r8
                dzsno(c,3) = dzsno(c,2)
                swice(c,3) = swice(c,2)
                swliq(c,3) = swliq(c,2)
                tsno(c,3) = tsno(c,2) - dtdz*dzsno(c,2)/2._r8
                if (tsno(c,3) >= tfrz) then 
                   tsno(c,3)  = tsno(c,2)
                else
                   tsno(c,2) = tsno(c,2) + dtdz*dzsno(c,2)/2._r8 
                endif

                mbc_phi(c,2) = mbc_phi(c,2)/2._r8
                mbc_phi(c,3) = mbc_phi(c,2)
                mbc_pho(c,2) = mbc_pho(c,2)/2._r8
                mbc_pho(c,3) = mbc_pho(c,2)
                moc_phi(c,2) = moc_phi(c,2)/2._r8
                moc_phi(c,3) = moc_phi(c,2)
                moc_pho(c,2) = moc_pho(c,2)/2._r8
                moc_pho(c,3) = moc_pho(c,2)
                mdst1(c,2) = mdst1(c,2)/2._r8
                mdst1(c,3) = mdst1(c,2)
                mdst2(c,2) = mdst2(c,2)/2._r8
                mdst2(c,3) = mdst2(c,2)
                mdst3(c,2) = mdst3(c,2)/2._r8
                mdst3(c,3) = mdst3(c,2)
                mdst4(c,2) = mdst4(c,2)/2._r8
                mdst4(c,3) = mdst4(c,2)
                rds(c,3) = rds(c,2)

             end if
          end if
       end if

       if (msno > 2) then
          if (dzsno(c,2) > 0.05_r8) then
             drr = dzsno(c,2) - 0.05_r8
             propor = drr/dzsno(c,2)
             zwice = propor*swice(c,2)
             zwliq = propor*swliq(c,2)
             
             zmbc_phi = propor*mbc_phi(c,2)
             zmbc_pho = propor*mbc_pho(c,2)
             zmoc_phi = propor*moc_phi(c,2)
             zmoc_pho = propor*moc_pho(c,2)
             zmdst1 = propor*mdst1(c,2)
             zmdst2 = propor*mdst2(c,2)
             zmdst3 = propor*mdst3(c,2)
             zmdst4 = propor*mdst4(c,2)

             propor = 0.05_r8/dzsno(c,2)
             swice(c,2) = propor*swice(c,2)
             swliq(c,2) = propor*swliq(c,2)

             mbc_phi(c,2) = propor*mbc_phi(c,2)
             mbc_pho(c,2) = propor*mbc_pho(c,2)
             moc_phi(c,2) = propor*moc_phi(c,2)
             moc_pho(c,2) = propor*moc_pho(c,2)
             mdst1(c,2) = propor*mdst1(c,2)
             mdst2(c,2) = propor*mdst2(c,2)
             mdst3(c,2) = propor*mdst3(c,2)
             mdst4(c,2) = propor*mdst4(c,2)

             dzsno(c,2) = 0.05_r8

             mbc_phi(c,3) = mbc_phi(c,3)+zmbc_phi  ! (combo)
             mbc_pho(c,3) = mbc_pho(c,3)+zmbc_pho  ! (combo)
             moc_phi(c,3) = moc_phi(c,3)+zmoc_phi  ! (combo)
             moc_pho(c,3) = moc_pho(c,3)+zmoc_pho  ! (combo)
             mdst1(c,3) = mdst1(c,3)+zmdst1  ! (combo)
             mdst2(c,3) = mdst2(c,3)+zmdst2  ! (combo)
             mdst3(c,3) = mdst3(c,3)+zmdst3  ! (combo)
             mdst4(c,3) = mdst4(c,3)+zmdst4  ! (combo)
             rds(c,3) = rds(c,2) ! (combo)

             call Combo (dzsno(c,3), swliq(c,3), swice(c,3), tsno(c,3), drr, &
                  zwliq, zwice, tsno(c,2))

             ! Subdivided a new layer
             if (msno <= 3 .and. dzsno(c,3) > 0.18_r8) then
                msno =  4
                dtdz = (tsno(c,2) - tsno(c,3))/((dzsno(c,2)+dzsno(c,3))/2._r8) 
                dzsno(c,3) = dzsno(c,3)/2._r8
                swice(c,3) = swice(c,3)/2._r8
                swliq(c,3) = swliq(c,3)/2._r8
                dzsno(c,4) = dzsno(c,3)
                swice(c,4) = swice(c,3)
                swliq(c,4) = swliq(c,3)
                tsno(c,4) = tsno(c,3) - dtdz*dzsno(c,3)/2._r8
                if (tsno(c,4) >= tfrz) then 
                   tsno(c,4)  = tsno(c,3)
                else
                   tsno(c,3) = tsno(c,3) + dtdz*dzsno(c,3)/2._r8 
                endif
                
                mbc_phi(c,3) = mbc_phi(c,3)/2._r8
                mbc_phi(c,4) = mbc_phi(c,3)
                mbc_pho(c,3) = mbc_pho(c,3)/2._r8
                mbc_pho(c,4) = mbc_pho(c,3)
                moc_phi(c,3) = moc_phi(c,3)/2._r8
                moc_phi(c,4) = moc_phi(c,3)
                moc_pho(c,3) = moc_pho(c,3)/2._r8
                moc_pho(c,4) = moc_pho(c,3)
                mdst1(c,3) = mdst1(c,3)/2._r8
                mdst1(c,4) = mdst1(c,3)
                mdst2(c,3) = mdst2(c,3)/2._r8
                mdst2(c,4) = mdst2(c,3)
                mdst3(c,3) = mdst3(c,3)/2._r8
                mdst3(c,4) = mdst3(c,3)
                mdst4(c,3) = mdst4(c,3)/2._r8
                mdst4(c,4) = mdst4(c,3)
                rds(c,4) = rds(c,3)

             end if
          end if
       end if

       if (msno > 3) then
          if (dzsno(c,3) > 0.11_r8) then
             drr = dzsno(c,3) - 0.11_r8
             propor = drr/dzsno(c,3)
             zwice = propor*swice(c,3)
             zwliq = propor*swliq(c,3)
             
             zmbc_phi = propor*mbc_phi(c,3)
             zmbc_pho = propor*mbc_pho(c,3)
             zmoc_phi = propor*moc_phi(c,3)
             zmoc_pho = propor*moc_pho(c,3)
             zmdst1 = propor*mdst1(c,3)
             zmdst2 = propor*mdst2(c,3)
             zmdst3 = propor*mdst3(c,3)
             zmdst4 = propor*mdst4(c,3)

             propor = 0.11_r8/dzsno(c,3)
             swice(c,3) = propor*swice(c,3)
             swliq(c,3) = propor*swliq(c,3)

             mbc_phi(c,3) = propor*mbc_phi(c,3)
             mbc_pho(c,3) = propor*mbc_pho(c,3)
             moc_phi(c,3) = propor*moc_phi(c,3)
             moc_pho(c,3) = propor*moc_pho(c,3)
             mdst1(c,3) = propor*mdst1(c,3)
             mdst2(c,3) = propor*mdst2(c,3)
             mdst3(c,3) = propor*mdst3(c,3)
             mdst4(c,3) = propor*mdst4(c,3)

             dzsno(c,3) = 0.11_r8

             mbc_phi(c,4) = mbc_phi(c,4)+zmbc_phi  ! (combo)
             mbc_pho(c,4) = mbc_pho(c,4)+zmbc_pho  ! (combo)
             moc_phi(c,4) = moc_phi(c,4)+zmoc_phi  ! (combo)
             moc_pho(c,4) = moc_pho(c,4)+zmoc_pho  ! (combo)
             mdst1(c,4) = mdst1(c,4)+zmdst1  ! (combo)
             mdst2(c,4) = mdst2(c,4)+zmdst2  ! (combo)
             mdst3(c,4) = mdst3(c,4)+zmdst3  ! (combo)
             mdst4(c,4) = mdst4(c,4)+zmdst4  ! (combo)
             rds(c,4) = rds(c,3) ! (combo)

             call Combo (dzsno(c,4), swliq(c,4), swice(c,4), tsno(c,4), drr, &
                  zwliq, zwice, tsno(c,3))

             ! Subdivided a new layer
             if (msno <= 4 .and. dzsno(c,4) > 0.41_r8) then
                msno = 5
                dtdz = (tsno(c,3) - tsno(c,4))/((dzsno(c,3)+dzsno(c,4))/2._r8) 
                dzsno(c,4) = dzsno(c,4)/2._r8
                swice(c,4) = swice(c,4)/2._r8
                swliq(c,4) = swliq(c,4)/2._r8
                dzsno(c,5) = dzsno(c,4)
                swice(c,5) = swice(c,4)
                swliq(c,5) = swliq(c,4)
                tsno(c,5) = tsno(c,4) - dtdz*dzsno(c,4)/2._r8 
                if (tsno(c,5) >= tfrz) then 
                   tsno(c,5)  = tsno(c,4)
                else
                   tsno(c,4) = tsno(c,4) + dtdz*dzsno(c,4)/2._r8 
                endif

                mbc_phi(c,4) = mbc_phi(c,4)/2._r8
                mbc_phi(c,5) = mbc_phi(c,4)
                mbc_pho(c,4) = mbc_pho(c,4)/2._r8
                mbc_pho(c,5) = mbc_pho(c,4)              
                moc_phi(c,4) = moc_phi(c,4)/2._r8
                moc_phi(c,5) = moc_phi(c,4)
                moc_pho(c,4) = moc_pho(c,4)/2._r8
                moc_pho(c,5) = moc_pho(c,4)
                mdst1(c,4) = mdst1(c,4)/2._r8
                mdst1(c,5) = mdst1(c,4)
                mdst2(c,4) = mdst2(c,4)/2._r8
                mdst2(c,5) = mdst2(c,4)
                mdst3(c,4) = mdst3(c,4)/2._r8
                mdst3(c,5) = mdst3(c,4)
                mdst4(c,4) = mdst4(c,4)/2._r8
                mdst4(c,5) = mdst4(c,4)
                rds(c,5) = rds(c,4)

             end if
          end if
       end if

       if (msno > 4) then
          if (dzsno(c,4) > 0.23_r8) then
             drr = dzsno(c,4) - 0.23_r8
             propor = drr/dzsno(c,4)
             zwice = propor*swice(c,4)
             zwliq = propor*swliq(c,4)
             
             zmbc_phi = propor*mbc_phi(c,4)
             zmbc_pho = propor*mbc_pho(c,4)
             zmoc_phi = propor*moc_phi(c,4)
             zmoc_pho = propor*moc_pho(c,4)
             zmdst1 = propor*mdst1(c,4)
             zmdst2 = propor*mdst2(c,4)
             zmdst3 = propor*mdst3(c,4)
             zmdst4 = propor*mdst4(c,4)

             propor = 0.23_r8/dzsno(c,4)
             swice(c,4) = propor*swice(c,4)
             swliq(c,4) = propor*swliq(c,4)

             mbc_phi(c,4) = propor*mbc_phi(c,4)
             mbc_pho(c,4) = propor*mbc_pho(c,4)
             moc_phi(c,4) = propor*moc_phi(c,4)
             moc_pho(c,4) = propor*moc_pho(c,4)
             mdst1(c,4) = propor*mdst1(c,4)
             mdst2(c,4) = propor*mdst2(c,4)
             mdst3(c,4) = propor*mdst3(c,4)
             mdst4(c,4) = propor*mdst4(c,4)

             dzsno(c,4) = 0.23_r8

             mbc_phi(c,5) = mbc_phi(c,5)+zmbc_phi  ! (combo)
             mbc_pho(c,5) = mbc_pho(c,5)+zmbc_pho  ! (combo)
             moc_phi(c,5) = moc_phi(c,5)+zmoc_phi  ! (combo)
             moc_pho(c,5) = moc_pho(c,5)+zmoc_pho  ! (combo)
             mdst1(c,5) = mdst1(c,5)+zmdst1  ! (combo)
             mdst2(c,5) = mdst2(c,5)+zmdst2  ! (combo)
             mdst3(c,5) = mdst3(c,5)+zmdst3  ! (combo)
             mdst4(c,5) = mdst4(c,5)+zmdst4  ! (combo)
             rds(c,5) = rds(c,4) ! (combo)

             call Combo (dzsno(c,5), swliq(c,5), swice(c,5), tsno(c,5), drr, &
                  zwliq, zwice, tsno(c,4))
          end if
       end if

       snl(c) = -msno

    end do

    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             dz(c,j)         = dzsno(c,j-snl(c))/frac_sno(c)
             h2osoi_ice(c,j) = swice(c,j-snl(c))
             h2osoi_liq(c,j) = swliq(c,j-snl(c))
             t_soisno(c,j)   = tsno(c,j-snl(c))

             mss_bcphi(c,j)   = mbc_phi(c,j-snl(c))
             mss_bcpho(c,j)   = mbc_pho(c,j-snl(c))
             mss_ocphi(c,j)   = moc_phi(c,j-snl(c))
             mss_ocpho(c,j)   = moc_pho(c,j-snl(c))
             mss_dst1(c,j)    = mdst1(c,j-snl(c))
             mss_dst2(c,j)    = mdst2(c,j-snl(c))
             mss_dst3(c,j)    = mdst3(c,j-snl(c))
             mss_dst4(c,j)    = mdst4(c,j-snl(c))
             snw_rds(c,j)     = rds(c,j-snl(c))

          end if
       end do
    end do

    do j = 0, -nlevsno+1, -1
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

  end subroutine DivideSnowLayers

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Combo
!
! !INTERFACE:
  subroutine Combo(dz,  wliq,  wice, t, dz2, wliq2, wice2, t2)
!
! !DESCRIPTION:
! Combines two elements and returns the following combined
! variables: dz, t, wliq, wice.
! The combined temperature is based on the equation:
! the sum of the enthalpies of the two elements =
! that of the combined element.
!
! !USES:
    use clm_varcon,  only : cpice, cpliq, tfrz, hfus
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)    :: dz2   ! nodal thickness of 2 elements being combined [m]
    real(r8), intent(in)    :: wliq2 ! liquid water of element 2 [kg/m2]
    real(r8), intent(in)    :: wice2 ! ice of element 2 [kg/m2]
    real(r8), intent(in)    :: t2    ! nodal temperature of element 2 [K]
    real(r8), intent(inout) :: dz    ! nodal thickness of 1 elements being combined [m]
    real(r8), intent(inout) :: wliq  ! liquid water of element 1
    real(r8), intent(inout) :: wice  ! ice of element 1 [kg/m2]
    real(r8), intent(inout) :: t     ! nodel temperature of elment 1 [K]
!
! !CALLED FROM:
! subroutine CombineSnowLayers in this module
! subroutine DivideSnowLayers in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: dzc   ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
    real(r8) :: wliqc ! Combined liquid water [kg/m2]
    real(r8) :: wicec ! Combined ice [kg/m2]
    real(r8) :: tc    ! Combined node temperature [K]
    real(r8) :: h     ! enthalpy of element 1 [J/m2]
    real(r8) :: h2    ! enthalpy of element 2 [J/m2]
    real(r8) :: hc    ! temporary
!-----------------------------------------------------------------------

    dzc = dz+dz2
    wicec = (wice+wice2)
    wliqc = (wliq+wliq2)
    h = (cpice*wice+cpliq*wliq) * (t-tfrz)+hfus*wliq
    h2= (cpice*wice2+cpliq*wliq2) * (t2-tfrz)+hfus*wliq2

    hc = h + h2
    tc = tfrz + (hc - hfus*wliqc) / (cpice*wicec + cpliq*wliqc)

    dz = dzc
    wice = wicec
    wliq = wliqc
    t = tc

  end subroutine Combo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BuildSnowFilter
!
! !INTERFACE:
  subroutine BuildSnowFilter(lbc, ubc, num_nolakec, filter_nolakec, &
                             num_snowc, filter_snowc, &
                             num_nosnowc, filter_nosnowc)
!
! !DESCRIPTION:
! Constructs snow filter for use in vectorized loops for snow hydrology.
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: lbc, ubc                    ! column bounds
    integer, intent(in)  :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in)  :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(out) :: num_snowc                   ! number of column snow points in column filter
    integer, intent(out) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
    integer, intent(out) :: num_nosnowc                 ! number of column non-snow points in column filter
    integer, intent(out) :: filter_nosnowc(ubc-lbc+1)   ! column filter for non-snow points
!
! !CALLED FROM:
! subroutine Hydrology2 in Hydrology2Mod
! subroutine CombineSnowLayers in this module
!
! !REVISION HISTORY:
! 2003 July 31: Forrest Hoffman
!
! !LOCAL VARIABLES:
! local pointers to implicit in arguments
    integer , pointer :: snl(:)                        ! number of snow layers
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer  :: fc, c
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtype components (column-level)

    snl => clm3%g%l%c%cps%snl

    ! Build snow/no-snow filters for other subroutines

    num_snowc = 0
    num_nosnowc = 0
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       if (snl(c) < 0) then
          num_snowc = num_snowc + 1
          filter_snowc(num_snowc) = c
       else
          num_nosnowc = num_nosnowc + 1
          filter_nosnowc(num_nosnowc) = c
       end if
    end do

  end subroutine BuildSnowFilter

!!!!!!!!!!!!!!!!!!! New subroutine for lakes
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: DivideSnowLayers_Lake
!
! !INTERFACE:
  subroutine DivideSnowLayers_Lake(lbc, ubc, num_snowc, filter_snowc)
!
! !DESCRIPTION:
! Subdivides snow layers if they exceed their prescribed maximum thickness.
!
! !USES:
    use clmtype
    use clm_varcon,  only : tfrz
    use SLakeCon  ,  only : lsadz
    use clm_varctl,  only : iulog
    use abortutils,  only : endrun
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)    :: lbc, ubc                    ! column bounds
    integer, intent(inout) :: num_snowc                   ! number of column snow points in column filter
    integer, intent(inout) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
!
! !CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/28/02, Peter Thornton: Migrated to new data structures.
! 2/29/08, David Lawrence: Snowpack T profile maintained during layer splitting
! 03/28/08, Mark Flanner: Added aerosol masses and snow grain radius
! 05/20/10, Zack Subin: Adjust all thicknesses + lsadz for lakes
!
! !LOCAL VARIABLES:
!
! local pointers to implicit inout arguments
!
    integer , pointer :: snl(:)            !number of snow layers
    real(r8), pointer :: dz(:,:)           !layer depth (m)
    real(r8), pointer :: zi(:,:)           !interface level below a "z" level (m)
    real(r8), pointer :: t_soisno(:,:)     !soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)   !ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)   !liquid water (kg/m2)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: z(:,:)          ! layer thickness (m)
    real(r8), pointer :: mss_bcphi(:,:)  ! hydrophilic BC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bcpho(:,:)  ! hydrophobic BC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocphi(:,:)  ! hydrophilic OC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocpho(:,:)  ! hydrophobic OC mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst1(:,:)   ! dust species 1 mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst2(:,:)   ! dust species 2 mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst3(:,:)   ! dust species 3 mass in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst4(:,:)   ! dust species 4 mass in snow (col,lyr) [kg]
    real(r8), pointer :: snw_rds(:,:)    ! effective snow grain radius (col,lyr) [microns, m^-6]
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: j, c, fc               ! indices
    real(r8) :: drr                    ! thickness of the combined [m]
    integer  :: msno                   ! number of snow layer 1 (top) to msno (bottom)
    real(r8) :: dzsno(lbc:ubc,nlevsno) ! Snow layer thickness [m]
    real(r8) :: swice(lbc:ubc,nlevsno) ! Partial volume of ice [m3/m3]
    real(r8) :: swliq(lbc:ubc,nlevsno) ! Partial volume of liquid water [m3/m3]
    real(r8) :: tsno(lbc:ubc ,nlevsno) ! Nodel temperature [K]
    real(r8) :: zwice                  ! temporary
    real(r8) :: zwliq                  ! temporary
    real(r8) :: propor                 ! temporary
    real(r8) :: dtdz                   ! temporary

    ! temporary variables mimicking the structure of other layer division variables
    real(r8) :: mbc_phi(lbc:ubc,nlevsno) ! mass of BC in each snow layer
    real(r8) :: zmbc_phi                 ! temporary
    real(r8) :: mbc_pho(lbc:ubc,nlevsno) ! mass of BC in each snow layer
    real(r8) :: zmbc_pho                 ! temporary
    real(r8) :: moc_phi(lbc:ubc,nlevsno) ! mass of OC in each snow layer
    real(r8) :: zmoc_phi                 ! temporary
    real(r8) :: moc_pho(lbc:ubc,nlevsno) ! mass of OC in each snow layer
    real(r8) :: zmoc_pho                 ! temporary
    real(r8) :: mdst1(lbc:ubc,nlevsno)   ! mass of dust 1 in each snow layer
    real(r8) :: zmdst1                   ! temporary
    real(r8) :: mdst2(lbc:ubc,nlevsno)   ! mass of dust 2 in each snow layer
    real(r8) :: zmdst2                   ! temporary
    real(r8) :: mdst3(lbc:ubc,nlevsno)   ! mass of dust 3 in each snow layer
    real(r8) :: zmdst3                   ! temporary
    real(r8) :: mdst4(lbc:ubc,nlevsno)   ! mass of dust 4 in each snow layer
    real(r8) :: zmdst4                   ! temporary
    real(r8) :: rds(lbc:ubc,nlevsno)

    ! Variables for consistency check
    real(r8) :: dztot(lbc:ubc), snwicetot(lbc:ubc), snwliqtot(lbc:ubc)
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtype components (column-level)

    snl        => clm3%g%l%c%cps%snl
    dz         => clm3%g%l%c%cps%dz
    zi         => clm3%g%l%c%cps%zi
    z          => clm3%g%l%c%cps%z
    t_soisno   => clm3%g%l%c%ces%t_soisno
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq
    mss_bcphi  => clm3%g%l%c%cps%mss_bcphi
    mss_bcpho  => clm3%g%l%c%cps%mss_bcpho
    mss_ocphi  => clm3%g%l%c%cps%mss_ocphi
    mss_ocpho  => clm3%g%l%c%cps%mss_ocpho
    mss_dst1   => clm3%g%l%c%cps%mss_dst1
    mss_dst2   => clm3%g%l%c%cps%mss_dst2
    mss_dst3   => clm3%g%l%c%cps%mss_dst3
    mss_dst4   => clm3%g%l%c%cps%mss_dst4
    snw_rds    => clm3%g%l%c%cps%snw_rds
    

    ! Initialize for consistency check
    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)

          if (j == -nlevsno+1) then
             dztot(c) = 0._r8
             snwicetot(c) = 0._r8
             snwliqtot(c) = 0._r8
          end if

          if (j >= snl(c)+1) then
             dztot(c) = dztot(c) + dz(c,j)
             snwicetot(c) = snwicetot(c) + h2osoi_ice(c,j)
             snwliqtot(c) = snwliqtot(c) + h2osoi_liq(c,j)
          end if
       end do
    end do


    ! Begin calculation - note that the following column loops are only invoked
    ! for snow-covered columns

    do j = 1,nlevsno
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j <= abs(snl(c))) then
             dzsno(c,j) = dz(c,j+snl(c))
             swice(c,j) = h2osoi_ice(c,j+snl(c))
             swliq(c,j) = h2osoi_liq(c,j+snl(c))
             tsno(c,j)  = t_soisno(c,j+snl(c))

             mbc_phi(c,j) = mss_bcphi(c,j+snl(c))
             mbc_pho(c,j) = mss_bcpho(c,j+snl(c))
             moc_phi(c,j) = mss_ocphi(c,j+snl(c))
             moc_pho(c,j) = mss_ocpho(c,j+snl(c))
             mdst1(c,j)   = mss_dst1(c,j+snl(c))
             mdst2(c,j)   = mss_dst2(c,j+snl(c))
             mdst3(c,j)   = mss_dst3(c,j+snl(c))
             mdst4(c,j)   = mss_dst4(c,j+snl(c))
             rds(c,j)     = snw_rds(c,j+snl(c))

          end if
       end do
    end do

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       msno = abs(snl(c))

       if (msno == 1) then
          ! Specify a new snow layer
          if (dzsno(c,1) > 0.03_r8 + 2._r8 * lsadz) then
             msno = 2
             dzsno(c,1) = dzsno(c,1)/2._r8
             swice(c,1) = swice(c,1)/2._r8
             swliq(c,1) = swliq(c,1)/2._r8
             dzsno(c,2) = dzsno(c,1)
             swice(c,2) = swice(c,1)
             swliq(c,2) = swliq(c,1)
             tsno(c,2)  = tsno(c,1)
             
             mbc_phi(c,1) = mbc_phi(c,1)/2._r8
             mbc_phi(c,2) = mbc_phi(c,1)
             mbc_pho(c,1) = mbc_pho(c,1)/2._r8
             mbc_pho(c,2) = mbc_pho(c,1)
             moc_phi(c,1) = moc_phi(c,1)/2._r8
             moc_phi(c,2) = moc_phi(c,1)
             moc_pho(c,1) = moc_pho(c,1)/2._r8
             moc_pho(c,2) = moc_pho(c,1)
             mdst1(c,1) = mdst1(c,1)/2._r8
             mdst1(c,2) = mdst1(c,1)
             mdst2(c,1) = mdst2(c,1)/2._r8
             mdst2(c,2) = mdst2(c,1)
             mdst3(c,1) = mdst3(c,1)/2._r8
             mdst3(c,2) = mdst3(c,1)
             mdst4(c,1) = mdst4(c,1)/2._r8
             mdst4(c,2) = mdst4(c,1)
             rds(c,2) = rds(c,1)

          end if
       end if

       if (msno > 1) then
          if (dzsno(c,1) > 0.02_r8 + lsadz) then
             drr = dzsno(c,1) - 0.02_r8 - lsadz
             propor = drr/dzsno(c,1)
             zwice = propor*swice(c,1)
             zwliq = propor*swliq(c,1)

             zmbc_phi = propor*mbc_phi(c,1)
             zmbc_pho = propor*mbc_pho(c,1)
             zmoc_phi = propor*moc_phi(c,1)
             zmoc_pho = propor*moc_pho(c,1)
             zmdst1 = propor*mdst1(c,1)
             zmdst2 = propor*mdst2(c,1)
             zmdst3 = propor*mdst3(c,1)
             zmdst4 = propor*mdst4(c,1)

             propor = (0.02_r8+lsadz)/dzsno(c,1)
             swice(c,1) = propor*swice(c,1)
             swliq(c,1) = propor*swliq(c,1)

             mbc_phi(c,1) = propor*mbc_phi(c,1)
             mbc_pho(c,1) = propor*mbc_pho(c,1)
             moc_phi(c,1) = propor*moc_phi(c,1)
             moc_pho(c,1) = propor*moc_pho(c,1)
             mdst1(c,1) = propor*mdst1(c,1)
             mdst2(c,1) = propor*mdst2(c,1)
             mdst3(c,1) = propor*mdst3(c,1)
             mdst4(c,1) = propor*mdst4(c,1)

             dzsno(c,1) = 0.02_r8 + lsadz

             mbc_phi(c,2) = mbc_phi(c,2)+zmbc_phi  ! (combo)
             mbc_pho(c,2) = mbc_pho(c,2)+zmbc_pho  ! (combo)
             moc_phi(c,2) = moc_phi(c,2)+zmoc_phi  ! (combo)
             moc_pho(c,2) = moc_pho(c,2)+zmoc_pho  ! (combo)
             mdst1(c,2) = mdst1(c,2)+zmdst1  ! (combo)
             mdst2(c,2) = mdst2(c,2)+zmdst2  ! (combo)
             mdst3(c,2) = mdst3(c,2)+zmdst3  ! (combo)
             mdst4(c,2) = mdst4(c,2)+zmdst4  ! (combo)
             rds(c,2) = rds(c,1) ! (combo)

             call Combo (dzsno(c,2), swliq(c,2), swice(c,2), tsno(c,2), drr, &
                  zwliq, zwice, tsno(c,1))

             ! Subdivide a new layer
             if (msno <= 2 .and. dzsno(c,2) > 0.07_r8+2._r8*lsadz) then
                msno = 3
                dtdz = (tsno(c,1) - tsno(c,2))/((dzsno(c,1)+dzsno(c,2))/2._r8) 
                dzsno(c,2) = dzsno(c,2)/2._r8
                swice(c,2) = swice(c,2)/2._r8
                swliq(c,2) = swliq(c,2)/2._r8
                dzsno(c,3) = dzsno(c,2)
                swice(c,3) = swice(c,2)
                swliq(c,3) = swliq(c,2)
                tsno(c,3) = tsno(c,2) - dtdz*dzsno(c,2)/2._r8
                if (tsno(c,3) >= tfrz) then 
                   tsno(c,3)  = tsno(c,2)
                else
                   tsno(c,2) = tsno(c,2) + dtdz*dzsno(c,2)/2._r8 
                endif

                mbc_phi(c,2) = mbc_phi(c,2)/2._r8
                mbc_phi(c,3) = mbc_phi(c,2)
                mbc_pho(c,2) = mbc_pho(c,2)/2._r8
                mbc_pho(c,3) = mbc_pho(c,2)
                moc_phi(c,2) = moc_phi(c,2)/2._r8
                moc_phi(c,3) = moc_phi(c,2)
                moc_pho(c,2) = moc_pho(c,2)/2._r8
                moc_pho(c,3) = moc_pho(c,2)
                mdst1(c,2) = mdst1(c,2)/2._r8
                mdst1(c,3) = mdst1(c,2)
                mdst2(c,2) = mdst2(c,2)/2._r8
                mdst2(c,3) = mdst2(c,2)
                mdst3(c,2) = mdst3(c,2)/2._r8
                mdst3(c,3) = mdst3(c,2)
                mdst4(c,2) = mdst4(c,2)/2._r8
                mdst4(c,3) = mdst4(c,2)
                rds(c,3) = rds(c,2)

             end if
          end if
       end if

       if (msno > 2) then
          if (dzsno(c,2) > 0.05_r8 + lsadz) then
             drr = dzsno(c,2) - 0.05_r8 - lsadz
             propor = drr/dzsno(c,2)
             zwice = propor*swice(c,2)
             zwliq = propor*swliq(c,2)
             
             zmbc_phi = propor*mbc_phi(c,2)
             zmbc_pho = propor*mbc_pho(c,2)
             zmoc_phi = propor*moc_phi(c,2)
             zmoc_pho = propor*moc_pho(c,2)
             zmdst1 = propor*mdst1(c,2)
             zmdst2 = propor*mdst2(c,2)
             zmdst3 = propor*mdst3(c,2)
             zmdst4 = propor*mdst4(c,2)

             propor = (0.05_r8+lsadz)/dzsno(c,2)
             swice(c,2) = propor*swice(c,2)
             swliq(c,2) = propor*swliq(c,2)

             mbc_phi(c,2) = propor*mbc_phi(c,2)
             mbc_pho(c,2) = propor*mbc_pho(c,2)
             moc_phi(c,2) = propor*moc_phi(c,2)
             moc_pho(c,2) = propor*moc_pho(c,2)
             mdst1(c,2) = propor*mdst1(c,2)
             mdst2(c,2) = propor*mdst2(c,2)
             mdst3(c,2) = propor*mdst3(c,2)
             mdst4(c,2) = propor*mdst4(c,2)

             dzsno(c,2) = 0.05_r8+lsadz

             mbc_phi(c,3) = mbc_phi(c,3)+zmbc_phi  ! (combo)
             mbc_pho(c,3) = mbc_pho(c,3)+zmbc_pho  ! (combo)
             moc_phi(c,3) = moc_phi(c,3)+zmoc_phi  ! (combo)
             moc_pho(c,3) = moc_pho(c,3)+zmoc_pho  ! (combo)
             mdst1(c,3) = mdst1(c,3)+zmdst1  ! (combo)
             mdst2(c,3) = mdst2(c,3)+zmdst2  ! (combo)
             mdst3(c,3) = mdst3(c,3)+zmdst3  ! (combo)
             mdst4(c,3) = mdst4(c,3)+zmdst4  ! (combo)
             rds(c,3) = rds(c,2) ! (combo)

             call Combo (dzsno(c,3), swliq(c,3), swice(c,3), tsno(c,3), drr, &
                  zwliq, zwice, tsno(c,2))

             ! Subdivided a new layer
             if (msno <= 3 .and. dzsno(c,3) > 0.18_r8+2._r8*lsadz) then
                msno =  4
                dtdz = (tsno(c,2) - tsno(c,3))/((dzsno(c,2)+dzsno(c,3))/2._r8) 
                dzsno(c,3) = dzsno(c,3)/2._r8
                swice(c,3) = swice(c,3)/2._r8
                swliq(c,3) = swliq(c,3)/2._r8
                dzsno(c,4) = dzsno(c,3)
                swice(c,4) = swice(c,3)
                swliq(c,4) = swliq(c,3)
                tsno(c,4) = tsno(c,3) - dtdz*dzsno(c,3)/2._r8
                if (tsno(c,4) >= tfrz) then 
                   tsno(c,4)  = tsno(c,3)
                else
                   tsno(c,3) = tsno(c,3) + dtdz*dzsno(c,3)/2._r8 
                endif
                
                mbc_phi(c,3) = mbc_phi(c,3)/2._r8
                mbc_phi(c,4) = mbc_phi(c,3)
                mbc_pho(c,3) = mbc_pho(c,3)/2._r8
                mbc_pho(c,4) = mbc_pho(c,3)
                moc_phi(c,3) = moc_phi(c,3)/2._r8
                moc_phi(c,4) = moc_phi(c,3)
                moc_pho(c,3) = moc_pho(c,3)/2._r8
                moc_pho(c,4) = moc_pho(c,3)
                mdst1(c,3) = mdst1(c,3)/2._r8
                mdst1(c,4) = mdst1(c,3)
                mdst2(c,3) = mdst2(c,3)/2._r8
                mdst2(c,4) = mdst2(c,3)
                mdst3(c,3) = mdst3(c,3)/2._r8
                mdst3(c,4) = mdst3(c,3)
                mdst4(c,3) = mdst4(c,3)/2._r8
                mdst4(c,4) = mdst4(c,3)
                rds(c,4) = rds(c,3)

             end if
          end if
       end if

       if (msno > 3) then
          if (dzsno(c,3) > 0.11_r8 + lsadz) then
             drr = dzsno(c,3) - 0.11_r8 - lsadz
             propor = drr/dzsno(c,3)
             zwice = propor*swice(c,3)
             zwliq = propor*swliq(c,3)
             
             zmbc_phi = propor*mbc_phi(c,3)
             zmbc_pho = propor*mbc_pho(c,3)
             zmoc_phi = propor*moc_phi(c,3)
             zmoc_pho = propor*moc_pho(c,3)
             zmdst1 = propor*mdst1(c,3)
             zmdst2 = propor*mdst2(c,3)
             zmdst3 = propor*mdst3(c,3)
             zmdst4 = propor*mdst4(c,3)

             propor = (0.11_r8+lsadz)/dzsno(c,3)
             swice(c,3) = propor*swice(c,3)
             swliq(c,3) = propor*swliq(c,3)

             mbc_phi(c,3) = propor*mbc_phi(c,3)
             mbc_pho(c,3) = propor*mbc_pho(c,3)
             moc_phi(c,3) = propor*moc_phi(c,3)
             moc_pho(c,3) = propor*moc_pho(c,3)
             mdst1(c,3) = propor*mdst1(c,3)
             mdst2(c,3) = propor*mdst2(c,3)
             mdst3(c,3) = propor*mdst3(c,3)
             mdst4(c,3) = propor*mdst4(c,3)

             dzsno(c,3) = 0.11_r8 + lsadz

             mbc_phi(c,4) = mbc_phi(c,4)+zmbc_phi  ! (combo)
             mbc_pho(c,4) = mbc_pho(c,4)+zmbc_pho  ! (combo)
             moc_phi(c,4) = moc_phi(c,4)+zmoc_phi  ! (combo)
             moc_pho(c,4) = moc_pho(c,4)+zmoc_pho  ! (combo)
             mdst1(c,4) = mdst1(c,4)+zmdst1  ! (combo)
             mdst2(c,4) = mdst2(c,4)+zmdst2  ! (combo)
             mdst3(c,4) = mdst3(c,4)+zmdst3  ! (combo)
             mdst4(c,4) = mdst4(c,4)+zmdst4  ! (combo)
             rds(c,4) = rds(c,3) ! (combo)

             call Combo (dzsno(c,4), swliq(c,4), swice(c,4), tsno(c,4), drr, &
                  zwliq, zwice, tsno(c,3))

             ! Subdivided a new layer
             if (msno <= 4 .and. dzsno(c,4) > 0.41_r8 + 2._r8*lsadz) then
                msno = 5
                dtdz = (tsno(c,3) - tsno(c,4))/((dzsno(c,3)+dzsno(c,4))/2._r8) 
                dzsno(c,4) = dzsno(c,4)/2._r8
                swice(c,4) = swice(c,4)/2._r8
                swliq(c,4) = swliq(c,4)/2._r8
                dzsno(c,5) = dzsno(c,4)
                swice(c,5) = swice(c,4)
                swliq(c,5) = swliq(c,4)
                tsno(c,5) = tsno(c,4) - dtdz*dzsno(c,4)/2._r8 
                if (tsno(c,5) >= tfrz) then 
                   tsno(c,5)  = tsno(c,4)
                else
                   tsno(c,4) = tsno(c,4) + dtdz*dzsno(c,4)/2._r8 
                endif

                mbc_phi(c,4) = mbc_phi(c,4)/2._r8
                mbc_phi(c,5) = mbc_phi(c,4)
                mbc_pho(c,4) = mbc_pho(c,4)/2._r8
                mbc_pho(c,5) = mbc_pho(c,4)              
                moc_phi(c,4) = moc_phi(c,4)/2._r8
                moc_phi(c,5) = moc_phi(c,4)
                moc_pho(c,4) = moc_pho(c,4)/2._r8
                moc_pho(c,5) = moc_pho(c,4)
                mdst1(c,4) = mdst1(c,4)/2._r8
                mdst1(c,5) = mdst1(c,4)
                mdst2(c,4) = mdst2(c,4)/2._r8
                mdst2(c,5) = mdst2(c,4)
                mdst3(c,4) = mdst3(c,4)/2._r8
                mdst3(c,5) = mdst3(c,4)
                mdst4(c,4) = mdst4(c,4)/2._r8
                mdst4(c,5) = mdst4(c,4)
                rds(c,5) = rds(c,4)

             end if
          end if
       end if

       if (msno > 4) then
          if (dzsno(c,4) > 0.23_r8 + lsadz) then
             drr = dzsno(c,4) - 0.23_r8 - lsadz
             propor = drr/dzsno(c,4)
             zwice = propor*swice(c,4)
             zwliq = propor*swliq(c,4)
             
             zmbc_phi = propor*mbc_phi(c,4)
             zmbc_pho = propor*mbc_pho(c,4)
             zmoc_phi = propor*moc_phi(c,4)
             zmoc_pho = propor*moc_pho(c,4)
             zmdst1 = propor*mdst1(c,4)
             zmdst2 = propor*mdst2(c,4)
             zmdst3 = propor*mdst3(c,4)
             zmdst4 = propor*mdst4(c,4)

             propor = (0.23_r8+lsadz)/dzsno(c,4)
             swice(c,4) = propor*swice(c,4)
             swliq(c,4) = propor*swliq(c,4)

             mbc_phi(c,4) = propor*mbc_phi(c,4)
             mbc_pho(c,4) = propor*mbc_pho(c,4)
             moc_phi(c,4) = propor*moc_phi(c,4)
             moc_pho(c,4) = propor*moc_pho(c,4)
             mdst1(c,4) = propor*mdst1(c,4)
             mdst2(c,4) = propor*mdst2(c,4)
             mdst3(c,4) = propor*mdst3(c,4)
             mdst4(c,4) = propor*mdst4(c,4)

             dzsno(c,4) = 0.23_r8 + lsadz

             mbc_phi(c,5) = mbc_phi(c,5)+zmbc_phi  ! (combo)
             mbc_pho(c,5) = mbc_pho(c,5)+zmbc_pho  ! (combo)
             moc_phi(c,5) = moc_phi(c,5)+zmoc_phi  ! (combo)
             moc_pho(c,5) = moc_pho(c,5)+zmoc_pho  ! (combo)
             mdst1(c,5) = mdst1(c,5)+zmdst1  ! (combo)
             mdst2(c,5) = mdst2(c,5)+zmdst2  ! (combo)
             mdst3(c,5) = mdst3(c,5)+zmdst3  ! (combo)
             mdst4(c,5) = mdst4(c,5)+zmdst4  ! (combo)
             rds(c,5) = rds(c,4) ! (combo)

             call Combo (dzsno(c,5), swliq(c,5), swice(c,5), tsno(c,5), drr, &
                  zwliq, zwice, tsno(c,4))
          end if
       end if

       snl(c) = -msno

    end do

    do j = -nlevsno+1,0
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             dz(c,j)         = dzsno(c,j-snl(c))
             h2osoi_ice(c,j) = swice(c,j-snl(c))
             h2osoi_liq(c,j) = swliq(c,j-snl(c))
             t_soisno(c,j)   = tsno(c,j-snl(c))

             mss_bcphi(c,j)   = mbc_phi(c,j-snl(c))
             mss_bcpho(c,j)   = mbc_pho(c,j-snl(c))
             mss_ocphi(c,j)   = moc_phi(c,j-snl(c))
             mss_ocpho(c,j)   = moc_pho(c,j-snl(c))
             mss_dst1(c,j)    = mdst1(c,j-snl(c))
             mss_dst2(c,j)    = mdst2(c,j-snl(c))
             mss_dst3(c,j)    = mdst3(c,j-snl(c))
             mss_dst4(c,j)    = mdst4(c,j-snl(c))
             snw_rds(c,j)     = rds(c,j-snl(c))

          end if
       end do
    end do

    ! Consistency check
    do j = -nlevsno + 1, 0
       do fc = 1, num_snowc
          c = filter_snowc(fc)

          if (j >= snl(c)+1) then
             dztot(c) = dztot(c) - dz(c,j)
             snwicetot(c) = snwicetot(c) - h2osoi_ice(c,j)
             snwliqtot(c) = snwliqtot(c) - h2osoi_liq(c,j)
          end if

          if (j == 0) then
             if ( abs(dztot(c)) > 1.e-10_r8 .or. abs(snwicetot(c)) > 1.e-7_r8 .or. &
                  abs(snwliqtot(c)) > 1.e-7_r8 ) then
                write(iulog,*)'Inconsistency in SnowDivision_Lake! c, remainders', &
                              'dztot, snwicetot, snwliqtot = ',c,dztot(c),snwicetot(c),snwliqtot(c)
                call endrun()
             end if
          end if
       end do
    end do

    do j = 0, -nlevsno+1, -1
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             z(c,j)    = zi(c,j) - 0.5_r8*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

  end subroutine DivideSnowLayers_Lake
!!!!!!!!!!!!!!!!!!!!! End new subroutine


end module SnowHydrologyMod
