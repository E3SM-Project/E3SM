module Hydrology1Mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  Hydrology1Mod
!
! !DESCRIPTION:
! Calculation of
! (1) water storage of intercepted precipitation
! (2) direct throughfall and canopy drainage of precipitation
! (3) the fraction of foliage covered by water and the fraction
!     of foliage that is dry and transpiring.
! (4) snow layer initialization if the snow accumulation exceeds 10 mm.
!
! !PUBLIC TYPES:
   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: Hydrology1
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
! !IROUTINE: Hydrology1
!
! !INTERFACE:
   subroutine Hydrology1(lbc, ubc, lbp, ubp, num_nolakec, filter_nolakec, &
                         num_nolakep, filter_nolakep)
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
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd   , only : clm_a2l
    use clm_varcon   , only : tfrz, istice, istwet, istsoil, istice_mec, isturb, &
                              icol_roof, icol_sunwall, icol_shadewall
    use clm_varcon   , only : istcrop
    use FracWetMod   , only : FracWet
    use clm_time_manager , only : get_step_size
    use subgridAveMod, only : p2c
    use SNICARMod    , only : snw_rds_min

!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                     ! pft bounds
    integer, intent(in) :: lbc, ubc                     ! column bounds
    integer, intent(in) :: num_nolakec                  ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)    ! column filter for non-lake points
    integer, intent(in) :: num_nolakep                  ! number of pft non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)    ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/15/02, Peter Thornton: Migrated to new data structures. Required
! adding a PFT loop.
! 4/26/05, Peter Thornton: Made the canopy interception factor fpi max=0.25
!   the default behavior
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arrays
!
    real(r8), pointer :: qflx_floodg(:)   ! gridcell flux of flood water from RTM
    real(r8), pointer :: qflx_floodc(:)   ! column flux of flood water from RTM
    integer , pointer :: cgridcell(:)      ! columns's gridcell
    integer , pointer :: clandunit(:)      ! columns's landunit
    integer , pointer :: pgridcell(:)      ! pft's gridcell
    integer , pointer :: plandunit(:)      ! pft's landunit
    integer , pointer :: pcolumn(:)        ! pft's column
    integer , pointer :: npfts(:)          ! number of pfts in column
    integer , pointer :: pfti(:)           ! column's beginning pft index
    integer , pointer :: ltype(:)          ! landunit type
    integer , pointer :: ctype(:)          ! column type
    real(r8), pointer :: forc_rain(:)      ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)      ! snow rate [mm/s]
    real(r8), pointer :: forc_t(:)         ! atmospheric temperature (Kelvin)
    logical , pointer :: do_capsnow(:)     ! true => do snow capping
    real(r8), pointer :: t_grnd(:)         ! ground temperature (Kelvin)
    real(r8), pointer :: dewmx(:)          ! Maximum allowed dew [mm]
    integer , pointer :: frac_veg_nosno(:) ! fraction of veg not covered by snow (0/1 now) [-]
    real(r8), pointer :: elai(:)           ! one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)           ! one-sided stem area index with burying by snow
    real(r8), pointer :: h2ocan_loss(:)    ! canopy water mass balance term (column)
    real(r8), pointer :: irrig_rate(:)     ! current irrigation rate (applied if n_irrig_steps_left > 0) [mm/s]
!
! local pointers to original implicit inout arrays
!
    integer , pointer :: snl(:)            ! number of snow layers
    real(r8), pointer :: snowdp(:)         ! snow height (m)
    real(r8), pointer :: h2osno(:)         ! snow water (mm H2O)
    real(r8), pointer :: h2ocan(:)         ! total canopy water (mm H2O)
    real(r8), pointer :: qflx_irrig(:)          ! irrigation amount (mm/s)
    integer, pointer  :: n_irrig_steps_left(:)  ! number of time steps for which we still need to irrigate today
!
! local pointers to original implicit out arrays
!
    real(r8), pointer :: qflx_prec_intr(:)     ! interception of precipitation [mm/s]
    real(r8), pointer :: qflx_prec_grnd(:)     ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_snwcp_liq(:)     ! excess rainfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_ice(:)     ! excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snow_grnd_pft(:) ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd_col(:) ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_rain_grnd(:)     ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: fwet(:)               ! fraction of canopy that is wet (0 to 1)
    real(r8), pointer :: fdry(:)               ! fraction of foliage that is green and dry [-] (new)
    real(r8), pointer :: zi(:,:)               ! interface level below a "z" level (m)
    real(r8), pointer :: dz(:,:)               ! layer depth (m)
    real(r8), pointer :: z(:,:)                ! layer thickness (m)
    real(r8), pointer :: t_soisno(:,:)         ! soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)       ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)       ! liquid water (kg/m2)
    real(r8), pointer :: frac_iceold(:,:)      ! fraction of ice relative to the tot water
    real(r8), pointer :: snw_rds(:,:)          ! effective snow grain radius (col,lyr) [microns, m^-6]
    real(r8), pointer :: mss_bcpho(:,:)        ! mass of hydrophobic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bcphi(:,:)        ! mass of hydrophilic BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bctot(:,:)        ! total mass of BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bc_col(:)         ! total column mass of BC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_bc_top(:)         ! total top-layer mass of BC (col,lyr) [kg]
    real(r8), pointer :: mss_ocpho(:,:)        ! mass of hydrophobic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_ocphi(:,:)        ! mass of hydrophilic OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_octot(:,:)        ! total mass of OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_oc_col(:)         ! total column mass of OC in snow (col,lyr) [kg]
    real(r8), pointer :: mss_oc_top(:)         ! total top-layer mass of OC (col,lyr) [kg]
    real(r8), pointer :: mss_dst1(:,:)         ! mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst2(:,:)         ! mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst3(:,:)         ! mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst4(:,:)         ! mass of dust species 4 in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dsttot(:,:)       ! total mass of dust in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst_col(:)        ! total column mass of dust in snow (col,lyr) [kg]
    real(r8), pointer :: mss_dst_top(:)        ! total top-layer mass of dust in snow (col,lyr) [kg]
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: f                            ! filter index
    integer  :: pi                           ! pft index
    integer  :: p                            ! pft index
    integer  :: c                            ! column index
    integer  :: l                            ! landunit index
    integer  :: g                            ! gridcell index
    integer  :: newnode                      ! flag when new snow node is set, (1=yes, 0=no)
    real(r8) :: dtime                        ! land model time step (sec)
    real(r8) :: h2ocanmx                     ! maximum allowed water on canopy [mm]
    real(r8) :: fpi                          ! coefficient of interception
    real(r8) :: xrun                         ! excess water that exceeds the leaf capacity [mm/s]
    real(r8) :: dz_snowf                     ! layer thickness rate change due to precipitation [mm/s]
    real(r8) :: bifall                       ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: fracsnow(lbp:ubp)            ! frac of precipitation that is snow
    real(r8) :: fracrain(lbp:ubp)            ! frac of precipitation that is rain
    real(r8) :: qflx_candrip(lbp:ubp)        ! rate of canopy runoff and snow falling off canopy [mm/s]
    real(r8) :: qflx_through_rain(lbp:ubp)   ! direct rain throughfall [mm/s]
    real(r8) :: qflx_through_snow(lbp:ubp)   ! direct snow throughfall [mm/s]
    real(r8) :: qflx_prec_grnd_snow(lbp:ubp) ! snow precipitation incident on ground [mm/s]
    real(r8) :: qflx_prec_grnd_rain(lbp:ubp) ! rain precipitation incident on ground [mm/s]
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    pgridcell          => pft%gridcell
    forc_rain          => clm_a2l%forc_rain
    forc_snow          => clm_a2l%forc_snow

    ! Assign local pointers to derived type members (landunit-level)

    ltype              => lun%itype

    ! Assign local pointers to derived type members (column-level)

    qflx_floodg        => clm_a2l%forc_flood
    qflx_floodc        => cwf%qflx_floodc
    cgridcell          => col%gridcell
    clandunit          => col%landunit
    ctype              => col%itype
    pfti               => col%pfti
    npfts              => col%npfts
    do_capsnow         => cps%do_capsnow
    forc_t             => ces%forc_t
    t_grnd             => ces%t_grnd
    snl                => cps%snl
    snowdp             => cps%snowdp
    h2osno             => cws%h2osno
    zi                 => cps%zi
    dz                 => cps%dz
    z                  => cps%z
    frac_iceold        => cps%frac_iceold
    t_soisno           => ces%t_soisno
    h2osoi_ice         => cws%h2osoi_ice
    h2osoi_liq         => cws%h2osoi_liq
    qflx_snow_grnd_col => pwf_a%qflx_snow_grnd
    h2ocan_loss        => cwf%h2ocan_loss
    snw_rds            => cps%snw_rds
    mss_bcpho          => cps%mss_bcpho
    mss_bcphi          => cps%mss_bcphi
    mss_bctot          => cps%mss_bctot
    mss_bc_col         => cps%mss_bc_col
    mss_bc_top         => cps%mss_bc_top
    mss_ocpho          => cps%mss_ocpho
    mss_ocphi          => cps%mss_ocphi
    mss_octot          => cps%mss_octot
    mss_oc_col         => cps%mss_oc_col
    mss_oc_top         => cps%mss_oc_top
    mss_dst1           => cps%mss_dst1
    mss_dst2           => cps%mss_dst2
    mss_dst3           => cps%mss_dst3
    mss_dst4           => cps%mss_dst4
    mss_dsttot         => cps%mss_dsttot
    mss_dst_col        => cps%mss_dst_col
    mss_dst_top        => cps%mss_dst_top

    ! Assign local pointers to derived type members (pft-level)

    plandunit          => pft%landunit
    pcolumn            => pft%column
    dewmx              => pps%dewmx
    frac_veg_nosno     => pps%frac_veg_nosno
    elai               => pps%elai
    esai               => pps%esai
    h2ocan             => pws%h2ocan
    qflx_prec_intr     => pwf%qflx_prec_intr
    qflx_prec_grnd     => pwf%qflx_prec_grnd
    qflx_snwcp_liq     => pwf%qflx_snwcp_liq
    qflx_snwcp_ice     => pwf%qflx_snwcp_ice
    qflx_snow_grnd_pft => pwf%qflx_snow_grnd
    qflx_rain_grnd     => pwf%qflx_rain_grnd
    fwet               => pps%fwet
    fdry               => pps%fdry
    irrig_rate         => cps%irrig_rate
    n_irrig_steps_left => cps%n_irrig_steps_left
    qflx_irrig         => cwf%qflx_irrig

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

       if (ltype(l)==istsoil .or. ltype(l)==istwet .or. ltype(l)==isturb .or. &
           ltype(l)==istcrop) then
          qflx_candrip(p) = 0._r8      ! rate of canopy runoff
          qflx_through_snow(p) = 0._r8 ! rain precipitation direct through canopy
          qflx_through_rain(p) = 0._r8 ! snow precipitation direct through canopy
          qflx_prec_intr(p) = 0._r8    ! total intercepted precipitation
          fracsnow(p) = 0._r8          ! fraction of input precip that is snow
          fracrain(p) = 0._r8          ! fraction of input precip that is rain

          if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall) then
             if (frac_veg_nosno(p) == 1 .and. (forc_rain(g) + forc_snow(g)) > 0._r8) then

                ! determine fraction of input precipitation that is snow and rain

                fracsnow(p) = forc_snow(g)/(forc_snow(g) + forc_rain(g))
                fracrain(p) = forc_rain(g)/(forc_snow(g) + forc_rain(g))
                
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
                qflx_through_snow(p) = forc_snow(g) * (1._r8-fpi)
                qflx_through_rain(p) = forc_rain(g) * (1._r8-fpi)
                
                ! Intercepted precipitation [mm/s]
                qflx_prec_intr(p) = (forc_snow(g) + forc_rain(g)) * fpi
                
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
       ! PET, 1/18/2005: Added new terms for mass balance correction
       ! due to dynamic pft weight shifting (column-level h2ocan_loss)
       ! Because the fractionation between rain and snow is indeterminate if
       ! rain + snow = 0, I am adding this very small flux only to the rain
       ! components.

       if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall) then
          if (frac_veg_nosno(p) == 0) then
             qflx_prec_grnd_snow(p) = forc_snow(g)
             qflx_prec_grnd_rain(p) = forc_rain(g) + h2ocan_loss(c)
          else
             qflx_prec_grnd_snow(p) = qflx_through_snow(p) + (qflx_candrip(p) * fracsnow(p))
             qflx_prec_grnd_rain(p) = qflx_through_rain(p) + (qflx_candrip(p) * fracrain(p)) + h2ocan_loss(c)
          end if
       ! Urban sunwall and shadewall have no intercepted precipitation
       else
          qflx_prec_grnd_snow(p) = 0.
          qflx_prec_grnd_rain(p) = 0.
       end if

       ! Determine whether we're irrigating here; set qflx_irrig appropriately
       if (n_irrig_steps_left(c) > 0) then
          qflx_irrig(c)         = irrig_rate(c)
          n_irrig_steps_left(c) = n_irrig_steps_left(c) - 1
       else
          qflx_irrig(c) = 0._r8
       end if

       ! Add irrigation water directly onto ground (bypassing canopy interception)
       ! Note that it's still possible that (some of) this irrigation water will runoff (as runoff is computed later)
       qflx_prec_grnd_rain(p) = qflx_prec_grnd_rain(p) + qflx_irrig(c)

       ! Done irrigation

       qflx_prec_grnd(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)

       if (do_capsnow(c)) then
          qflx_snwcp_liq(p) = qflx_prec_grnd_rain(p)
          qflx_snwcp_ice(p) = qflx_prec_grnd_snow(p)

          qflx_snow_grnd_pft(p) = 0._r8
          qflx_rain_grnd(p) = 0._r8
       else
          qflx_snwcp_liq(p) = 0._r8
          qflx_snwcp_ice(p) = 0._r8
          qflx_snow_grnd_pft(p) = qflx_prec_grnd_snow(p)           ! ice onto ground (mm/s)
          qflx_rain_grnd(p)     = qflx_prec_grnd_rain(p)           ! liquid water onto ground (mm/s)
       end if

    end do ! (end pft loop)

    ! Determine the fraction of foliage covered by water and the
    ! fraction of foliage that is dry and transpiring.

    call FracWet(num_nolakep, filter_nolakep)

    ! Update column level state variables for snow.

    call p2c(num_nolakec, filter_nolakec, qflx_snow_grnd_pft, qflx_snow_grnd_col)

!rtm_flood: apply gridcell flood water flux to non-lake columns
!                    no inputs to urban wall columns, as above with atm inputs
!dir$ concurrent
!cdir nodep
    do f = 1, num_nolakec
       c = filter_nolakec(f)
       g = cgridcell(c)
       if (ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall) then      
          qflx_floodc(c) = qflx_floodg(g)
       else
          qflx_floodc(c) = 0._r8
       endif    
    enddo
!rtm_flood

    ! Determine snow height and snow water

    do f = 1, num_nolakec
       c = filter_nolakec(f)
       l = clandunit(c)
       g = cgridcell(c)

       ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
       ! U.S.Department of Agriculture Forest Service, Project F,
       ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

       if (do_capsnow(c)) then
          dz_snowf = 0._r8
       else
          if (forc_t(c) > tfrz + 2._r8) then
             bifall=50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
          else if (forc_t(c) > tfrz - 15._r8) then
             bifall=50._r8 + 1.7_r8*(forc_t(c) - tfrz + 15._r8)**1.5_r8
          else
             bifall=50._r8
          end if
          dz_snowf = qflx_snow_grnd_col(c)/bifall
          snowdp(c) = snowdp(c) + dz_snowf*dtime
          h2osno(c) = h2osno(c) + qflx_snow_grnd_col(c)*dtime  ! snow water equivalent (mm)

       end if

       if (ltype(l)==istwet .and. t_grnd(c)>tfrz) then
          h2osno(c)=0._r8
          snowdp(c)=0._r8
       end if

       ! When the snow accumulation exceeds 10 mm, initialize snow layer
       ! Currently, the water temperature for the precipitation is simply set
       ! as the surface air temperature

       newnode = 0    ! flag for when snow node will be initialized
       if (snl(c) == 0 .and. qflx_snow_grnd_col(c) > 0.0_r8 .and. snowdp(c) >= 0.01_r8) then
          newnode = 1
          snl(c) = -1
          dz(c,0) = snowdp(c)                       ! meter
          z(c,0) = -0.5_r8*dz(c,0)
          zi(c,-1) = -dz(c,0)
          t_soisno(c,0) = min(tfrz, forc_t(c))      ! K
          h2osoi_ice(c,0) = h2osno(c)               ! kg/m2
          h2osoi_liq(c,0) = 0._r8                   ! kg/m2
          frac_iceold(c,0) = 1._r8
       

          ! intitialize SNICAR variables for fresh snow:
          snw_rds(c,0)    = snw_rds_min

          mss_bcpho(c,:)  = 0._r8
          mss_bcphi(c,:)  = 0._r8
          mss_bctot(c,:)  = 0._r8
          mss_bc_col(c)   = 0._r8
          mss_bc_top(c)   = 0._r8

          mss_ocpho(c,:)  = 0._r8
          mss_ocphi(c,:)  = 0._r8
          mss_octot(c,:)  = 0._r8
          mss_oc_col(c)   = 0._r8
          mss_oc_top(c)   = 0._r8

          mss_dst1(c,:)   = 0._r8
          mss_dst2(c,:)   = 0._r8
          mss_dst3(c,:)   = 0._r8
          mss_dst4(c,:)   = 0._r8
          mss_dsttot(c,:) = 0._r8
          mss_dst_col(c)  = 0._r8
          mss_dst_top(c)  = 0._r8
       end if

       ! The change of ice partial density of surface node due to precipitation.
       ! Only ice part of snowfall is added here, the liquid part will be added
       ! later.

       if (snl(c) < 0 .and. newnode == 0) then
          h2osoi_ice(c,snl(c)+1) = h2osoi_ice(c,snl(c)+1)+dtime*qflx_snow_grnd_col(c)
          dz(c,snl(c)+1) = dz(c,snl(c)+1)+dz_snowf*dtime
       end if

    end do

  end subroutine Hydrology1

end module Hydrology1Mod
