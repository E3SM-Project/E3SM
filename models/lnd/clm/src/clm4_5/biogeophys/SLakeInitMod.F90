module SLakeInitMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains time constant initialization code for S Lake scheme.
  !
  ! !USES
  use decompMod         , only : bounds_type
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun
  use spmdMod           , only : masterproc
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initTimeConstSlake ! driver
  public :: initColdSlake
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initTimeConstSlake(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize time invariant clm variables for S Lake code (and h2osoi_vol).
    !
    ! !USES:
    use clmtype          , only : col, cps, cws, lun, namec
    use decompMod        , only : get_proc_bounds, get_proc_global
    use clm_varpar       , only : nlevsoi, nlevlak, nlevgrnd
    use clm_varcon       , only : istdlak, zlak, dzlak, zsoi, dzsoi, zisoi, spval
    use clm_varcon       , only : secspday, denh2o, denice
    use clm_varctl       , only : nsrest, iulog, use_extralakelayers
    use clm_time_manager , only : get_step_size
    use SLakeCon         , only : fcrit, minz0lake, depthcrit, mixfact, pudz
    use SLakeCon         , only : lakepuddling, lake_puddle_thick, deepmixing_depthcrit, deepmixing_mixfact
    use SLakeCon         , only : lake_use_old_fcrit_minz0, lake_melt_icealb, alblakwi
    use CNSharedParamsMod, only : CNParamsShareInst
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: n,i,j,ib,lev,bottom               ! indices
    integer  :: g,l,c,p                           ! indices
    real(r8) :: bd                                ! bulk density of dry soil material [kg/m^3]
    real(r8) :: tkm                               ! mineral conductivity
    real(r8) :: xksat                             ! maximum hydraulic conductivity of soil [mm/s]
    real(r8) :: scalez = 0.025_r8                 ! Soil layer thickness discretization (m)
    real(r8) :: depthratio                        ! ratio of lake depth to standard deep lake depth 
    real(r8) :: clay                              ! temporary
    real(r8) :: sand                              ! temporary
    real(r8) :: om_frac                           ! organic matter fraction
    real(r8) :: om_watsat    = 0.9_r8             ! porosity of organic soil
    real(r8) :: om_hksat     = 0.1_r8             ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8) :: om_tkm       = 0.25_r8            ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(r8) :: om_sucsat    = 10.3_r8            ! saturated suction for organic matter (Letts, 2000)
    real(r8) :: om_csol      = 2.5_r8             ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(r8) :: om_tkd       = 0.05_r8            ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(r8) :: om_b         = 2.7_r8             ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8) :: organic_max                       ! organic matter (kg/m3) where soil is assumed to act like peat
    real(r8) :: csol_bedrock = 2.0e6_r8           ! vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
    real(r8) :: pc           = 0.5_r8             ! percolation threshold
    real(r8) :: pcbeta       = 0.139_r8           ! percolation exponent

    ! Note: These constants should be shared with iniTimeConst.
    real(r8) :: perc_frac                         ! "percolating" fraction of organic soil
    real(r8) :: perc_norm                         ! normalize to 1 when 100% organic soil
    real(r8) :: uncon_hksat                       ! series conductivity of mineral/organic soil
    real(r8) :: uncon_frac                        ! fraction of "unconnected" soil
    real(r8) :: dtime                             ! land model timestep (s)
    real(r8), parameter :: dtime_lsadz = 1800._r8 ! timestep used for nominal lsadz
    real(r8), parameter :: snodzmin = 0.01_r8     ! nominal minimum snow layer thickness
                                                  ! used in SnowHydrology, etc. This should be a global constant.
    !------------------------------------------------------------------------

    if (masterproc) write (iulog,*) 'Attempting to initialize time invariant variables for lakes'

    organic_max = CNParamsShareInst%organic_max

    ! Set SLakeCon constants according to namelist fields
    if (lake_use_old_fcrit_minz0) then
       ! critical dimensionless fetch for Charnock parameter. From Vickers & Mahrt 1997
       ! but converted to use u instead of u* (Form used in Subin et al. 2011)
       fcrit   = 22._r8 

       ! (m) Minimum allowed roughness length for unfrozen lakes.
       ! (Used in Subin et al. 2011)
       minz0lake = 1.e-5_r8        
    else
       ! Vickers & Mahrt 1997
       fcrit   = 100._r8  

       ! (m) Minimum allowed roughness length for unfrozen lakes.
       ! Now set low so it is only to avoid floating point exceptions.
       minz0lake = 1.e-10_r8       
    end if

    if (lakepuddling) then
       ! (m) Minimum total ice thickness required to allow lake puddling. Default is 0.2m.
       ! This option has not been extensively tested.
       ! This option turns on lake_no_melt_icealb, as the decrease in albedo will be based
       ! on whether there is water over nice, not purely a function of ice top temperature.
       pudz = lake_puddle_thick    
    end if

    ! (m) Depth beneath which to increase mixing. See discussion in Subin et al. 2011
    depthcrit = deepmixing_depthcrit 

    ! Mixing increase factor. ! Defaults are 25 m, increase by 10.
    ! Note some other namelists will be used directly in lake physics during model integration.
    mixfact = deepmixing_mixfact 

    ! Set alblakwi
    alblakwi(:) = lake_melt_icealb(:)

    ! Lake layers
    if (.not. use_extralakelayers) then
       dzlak(1) = 0.1_r8
       dzlak(2) = 1._r8
       dzlak(3) = 2._r8
       dzlak(4) = 3._r8
       dzlak(5) = 4._r8
       dzlak(6) = 5._r8
       dzlak(7) = 7._r8
       dzlak(8) = 7._r8
       dzlak(9) = 10.45_r8
       dzlak(10)= 10.45_r8

       zlak(1) =  0.05_r8
       zlak(2) =  0.6_r8
       zlak(3) =  2.1_r8
       zlak(4) =  4.6_r8
       zlak(5) =  8.1_r8
       zlak(6) = 12.6_r8
       zlak(7) = 18.6_r8
       zlak(8) = 25.6_r8
       zlak(9) = 34.325_r8
       zlak(10)= 44.775_r8
    else
       dzlak(1) =0.1_r8
       dzlak(2) =0.25_r8
       dzlak(3) =0.25_r8
       dzlak(4) =0.25_r8
       dzlak(5) =0.25_r8
       dzlak(6) =0.5_r8
       dzlak(7) =0.5_r8
       dzlak(8) =0.5_r8
       dzlak(9) =0.5_r8
       dzlak(10) =0.75_r8
       dzlak(11) =0.75_r8
       dzlak(12) =0.75_r8
       dzlak(13) =0.75_r8
       dzlak(14) =2_r8
       dzlak(15) =2_r8
       dzlak(16) =2.5_r8
       dzlak(17) =2.5_r8
       dzlak(18) =3.5_r8
       dzlak(19) =3.5_r8
       dzlak(20) =3.5_r8
       dzlak(21) =3.5_r8
       dzlak(22) =5.225_r8
       dzlak(23) =5.225_r8
       dzlak(24) =5.225_r8
       dzlak(25) =5.225_r8

       zlak(1) = dzlak(1)/2._r8
       do i=2,nlevlak
          zlak(i) = zlak(i-1) + (dzlak(i-1)+dzlak(i))/2._r8
       end do
    end if

    ! --------------------------------------------------------------------
    ! Initialize soil and lake levels
    ! Initialize soil thermal and hydraulic properties [color is not needed]
    ! --------------------------------------------------------------------

    ! Column level initialization
    do c = bounds%begc, bounds%endc

       ! Set gridcell and landunit indices
       g = col%gridcell(c)
       l = col%landunit(c)

       ! Soil hydraulic and thermal properties
       if (lun%itype(l)==istdlak) then
          do lev = 1,nlevgrnd
             if ( lev <= nlevsoi )then
                clay    = cps%cellclay(c,lev)
                sand    = cps%cellsand(c,lev)
                om_frac =(cps%cellorg(c,lev)/organic_max)**2._r8
             else
                clay    = cps%cellclay(c,nlevsoi)
                sand    = cps%cellsand(c,nlevsoi)
                om_frac = 0.0_r8
             end if
             cps%watsat(c,lev) = 0.489_r8 - 0.00126_r8*sand
             cps%bsw(c,lev)    = 2.91 + 0.159*clay
             cps%sucsat(c,lev) = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )
             bd                = (1._r8-cps%watsat(c,lev))*2.7e3_r8
             cps%watsat(c,lev) = (1._r8 - om_frac)*cps%watsat(c,lev) + om_watsat*om_frac
             tkm               = (1._r8-om_frac)*(8.80_r8*sand+2.92_r8*clay)/(sand+clay)+om_tkm*om_frac ! W/(m K)
             cps%bsw(c,lev)    = (1._r8-om_frac)*(2.91_r8 + 0.159_r8*clay) + om_frac*om_b
             cps%sucsat(c,lev) = (1._r8-om_frac)*cps%sucsat(c,lev) + om_sucsat*om_frac
             xksat             = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s

             ! perc_frac is zero unless perf_frac greater than percolation threshold
             if (om_frac > pc) then
                perc_norm=(1._r8 - pc)**(-pcbeta)
                perc_frac=perc_norm*(om_frac - pc)**pcbeta
             else
                perc_frac=0._r8
             endif

             ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
             uncon_frac=(1._r8-om_frac)+(1._r8-perc_frac)*om_frac

             ! uncon_hksat is series addition of mineral/organic conductivites
             if (om_frac .lt. 1._r8) then
                uncon_hksat=uncon_frac/((1._r8-om_frac)/xksat &
                     +((1._r8-perc_frac)*om_frac)/om_hksat)
             else
                uncon_hksat = 0._r8
             end if

             cps%hksat(c,lev)  = uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat

             cps%tkmg(c,lev)   = tkm ** (1._r8- cps%watsat(c,lev))
             cps%tksatu(c,lev) = cps%tkmg(c,lev)*0.57_r8**cps%watsat(c,lev)
             cps%tkdry(c,lev)  = ((0.135_r8*bd + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd))*(1._r8-om_frac) + &
                  om_tkd*om_frac
             cps%csol(c,lev)   = ((1._r8-om_frac)*(2.128_r8*sand+2.385_r8*clay) / (sand+clay) +   &
                  om_csol*om_frac)*1.e6_r8  ! J/(m3 K)
             if (lev .gt. nlevsoi) then
                cps%csol(c,lev) = csol_bedrock
             endif
             cps%watdry(c,lev) = cps%watsat(c,lev) * (316230._r8/cps%sucsat(c,lev)) ** (-1._r8/cps%bsw(c,lev))
             cps%watopt(c,lev) = cps%watsat(c,lev) * (158490._r8/cps%sucsat(c,lev)) ** (-1._r8/cps%bsw(c,lev))

             !! added by K.Sakaguchi for beta from Lee and Pielke, 1992
             ! water content at field capacity, defined as hk = 0.1 mm/day
             ! used eqn (7.70) in CLM3 technote with k = 0.1 (mm/day) / (# seconds/day)
             cps%watfc(c,lev) = cps%watsat(c,lev) * (0.1_r8 / &
                               (cps%hksat(c,lev)*secspday))**(1._r8/(2._r8*cps%bsw(c,lev)+3._r8))
          end do
       endif

       ! Define levels
       if (lun%itype(l) == istdlak) then

          if (cps%lakedepth(c) == spval) then
             cps%lakedepth(c)         = zlak(nlevlak) + 0.5_r8*dzlak(nlevlak)
             cps%z_lake(c,1:nlevlak)  = zlak(1:nlevlak)
             cps%dz_lake(c,1:nlevlak) = dzlak(1:nlevlak)

          else if (cps%lakedepth(c) > 1._r8 .and. cps%lakedepth(c) < 5000._r8) then

             depthratio                 = cps%lakedepth(c) / (zlak(nlevlak) + 0.5_r8*dzlak(nlevlak)) 
             cps%z_lake(c,1)            = zlak(1)
             cps%dz_lake(c,1)           = dzlak(1)
             cps%dz_lake(c,2:nlevlak-1) = dzlak(2:nlevlak-1)*depthratio
             cps%dz_lake(c,nlevlak)     = dzlak(nlevlak)*depthratio - (cps%dz_lake(c,1) - dzlak(1)*depthratio)
             do lev=2,nlevlak
                cps%z_lake(c,lev) = cps%z_lake(c,lev-1) + (cps%dz_lake(c,lev-1)+cps%dz_lake(c,lev))/2._r8
             end do

          else if (cps%lakedepth(c) > 0._r8 .and. cps%lakedepth(c) <= 1._r8) then

             cps%dz_lake(c,:) = cps%lakedepth(c) / nlevlak;
             cps%z_lake(c,1)  = cps%dz_lake(c,1) / 2._r8;
             do lev=2,nlevlak
                cps%z_lake(c,lev) = cps%z_lake(c,lev-1) + (cps%dz_lake(c,lev-1)+cps%dz_lake(c,lev))/2._r8
             end do

          else

             write(iulog,*)'Bad lake depth: lakedepth: ', cps%lakedepth(c)
             call endrun(decomp_index=c, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))

          end if

          cps%z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
          cps%zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
          cps%dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
       end if

       ! Initialize h2osoi_vol. Moved here from BiogeophysRest to avoid problems with dz not being defined for lakes
       ! there, and to avoid adding dz to the restart file.
       ! Note (Attn EK): if this subroutine is ever called before SLakeRest, then this code should be moved there.

       if (lun%itype(l) == istdlak .and. cws%h2osoi_liq(c,1) /= spval) then
          ! If not spval then this was available on the restart file and it was not an old restart file.
          ! Otherwise it will be set in makearbinit.
          do lev = 1,nlevgrnd
             cws%h2osoi_vol(c,lev) = cws%h2osoi_liq(c,lev)/(cps%dz(c,lev)*denh2o) + &
                                     cws%h2osoi_ice(c,lev)/(cps%dz(c,lev)*denice)
          end do
       end if

    end do

    if (masterproc) write (iulog,*) 'Successfully initialized time invariant variables for lakes'

  end subroutine initTimeConstSlake

  !-----------------------------------------------------------------------
  subroutine initColdSLake(bounds)
    !
    ! !DESCRIPTION:
    ! Initializes "cold start" time varying variables for SLake
    !
    ! !USES:
    use clmtype          , only : col, cps, cws, ces, lun, namec
    use clm_varpar       , only : nlevsno, nlevlak, nlevgrnd, nlevsoi
    use clm_varcon       , only : bdsno, denice, denh2o, spval, tfrz, tkwat
    use clm_varctl       , only : iulog
    use clm_time_manager , only : get_step_size
    use SNICARMod        , only : snw_rds_min
    use SLakeCon         , only : lsadz
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,p                           ! indices
    real(r8) :: dtime                             ! land model timestep (s)
    real(r8), parameter :: dtime_lsadz = 1800._r8 ! timestep used for nominal lsadz
    real(r8), parameter :: snodzmin = 0.01_r8     ! nominal minimum snow layer thickness
                                                  ! used in SnowHydrology, etc. This should be a global constant.
    !-----------------------------------------------------------------------

    ! lun%lakpoi         Input:  [logical (:)    ]  true => landunit is a lake point         
    ! ces%t_lake         Input:  [real(r8) (:,:) ]  lake temperature (Kelvin)  (1:nlevlak)
    ! cps%watsat         Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity) (nlevgrnd)

    ! cps%snl            Output: [integer (:)    ]  number of snow layers                    
    ! cps%dz             Output: [real(r8) (:,:) ]  layer thickness depth (m)             
    ! cws%h2osoi_ice     Output: [real(r8) (:,:) ]  ice lens (kg/m2)                      
    ! cws%h2osoi_liq     Output: [real(r8) (:,:) ]  liquid water (kg/m2)                  
    ! cws%h2osoi_vol     Output: [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    ! ces%t_soisno       Output: [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    ! cws%lake_icefrac   Output: [real(r8) (:,:) ]  mass fraction of lake layer that is frozen
    ! cps%savedtke1      Output: [real(r8) (:)   ]  top level eddy conductivity (W/mK)      
    ! cps%ust_lake       Output: [real(r8) (:)   ]  friction velocity (m/s)                 
    ! cps%z0mg           Output: [real(r8) (:)   ]  roughness length over ground, momentum [m]
    ! cps%snw_rds        Output: [real(r8) (:,:) ]  effective snow grain radius (col,lyr) [microns, m^-6] (SNICAR variable)
    ! cps%snw_rds_top    Output: [real(r8) (:)   ]  snow grain size, top (col) [microns] (SNICAR variable)   
    ! cps%sno_liq_top    Output: [real(r8) (:)   ]  liquid water fraction (mass) in top snow layer (col) [frc] (SNICAR variable)

    if ( masterproc )  then
       write (iulog,*) 'Setting initial data to non-spun up values for lake points for CLM4-LISSS Lake Model.'
    end if

    associate(snl => cps%snl)  ! Output: [integer (:)]  number of snow layers                    
      
    ! Adjust lsadz (added to min. snow layer thickness) if dtime /= 1800 s.

    dtime = get_step_size()
    if (dtime /= dtime_lsadz) then
       lsadz = (lsadz + snodzmin) * (dtime/dtime_lsadz)**0.5_r8  - snodzmin
       lsadz = max(lsadz, 0._r8)
       if (masterproc) then
          write(iulog,*)'CLM timestep appears to differ from 1800s. Adjusting minimum ', &
               'snow layer thickness over lakes to maintain stability. New minimum thickness is ', &
               lsadz + snodzmin, '.' 
       end if
    end if
    
    ! Determine snow levels and interfaces for lake points
    ! lsadz (default 0.03m) is added to min & max thicknesses to not stress lake numerics at 30 min. timestep.
    
    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       if (lun%lakpoi(l)) then

          if (cps%snow_depth(c) < 0.01_r8 + lsadz) then
             snl(c) = 0
             cps%dz(c,-nlevsno+1:0) = 0._r8
             cps%z (c,-nlevsno+1:0) = 0._r8
             cps%zi(c,-nlevsno+0:0) = 0._r8
          else
             if ((cps%snow_depth(c) >= 0.01_r8 + lsadz) .and. (cps%snow_depth(c) <= 0.03_r8 + lsadz)) then
                snl(c) = -1
                cps%dz(c,0)  = cps%snow_depth(c)
             else if ((cps%snow_depth(c) > 0.03_r8 + 2._r8*lsadz) .and. (cps%snow_depth(c) <= 0.04_r8 + 2._r8*lsadz)) then
                snl(c) = -2
                cps%dz(c,-1) = cps%snow_depth(c)/2._r8
                cps%dz(c, 0) = cps%dz(c,-1)
             else if ((cps%snow_depth(c) > 0.04_r8 + 2._r8*lsadz) .and. (cps%snow_depth(c) <= 0.07_r8 + 2._r8*lsadz)) then
                snl(c) = -2
                cps%dz(c,-1) = 0.02_r8 + lsadz
                cps%dz(c, 0) = cps%snow_depth(c) - cps%dz(c,-1)
             else if ((cps%snow_depth(c) > 0.07_r8 + 3._r8*lsadz) .and. (cps%snow_depth(c) <= 0.12_r8 + 3._r8*lsadz)) then
                snl(c) = -3
                cps%dz(c,-2) = 0.02_r8 + lsadz
                cps%dz(c,-1) = (cps%snow_depth(c) - 0.02_r8 - lsadz)/2._r8
                cps%dz(c, 0) = cps%dz(c,-1)
             else if ((cps%snow_depth(c) > 0.12_r8 + 3._r8*lsadz) .and. (cps%snow_depth(c) <= 0.18_r8 + 3._r8*lsadz)) then
                snl(c) = -3
                cps%dz(c,-2) = 0.02_r8 + lsadz
                cps%dz(c,-1) = 0.05_r8 + lsadz
                cps%dz(c, 0) = cps%snow_depth(c) - cps%dz(c,-2) - cps%dz(c,-1)
             else if ((cps%snow_depth(c) > 0.18_r8 + 4._r8*lsadz) .and. (cps%snow_depth(c) <= 0.29_r8 + 4._r8*lsadz)) then
                snl(c) = -4
                cps%dz(c,-3) = 0.02_r8 + lsadz
                cps%dz(c,-2) = 0.05_r8 + lsadz
                cps%dz(c,-1) = (cps%snow_depth(c) - cps%dz(c,-3) - cps%dz(c,-2))/2._r8
                cps%dz(c, 0) = cps%dz(c,-1)
             else if ((cps%snow_depth(c) > 0.29_r8 + 4._r8*lsadz) .and. (cps%snow_depth(c) <= 0.41_r8 + 4._r8*lsadz)) then
                snl(c) = -4
                cps%dz(c,-3) = 0.02_r8 + lsadz
                cps%dz(c,-2) = 0.05_r8 + lsadz
                cps%dz(c,-1) = 0.11_r8 + lsadz
                cps%dz(c, 0) = cps%snow_depth(c) - cps%dz(c,-3) - cps%dz(c,-2) - cps%dz(c,-1)
             else if ((cps%snow_depth(c) > 0.41_r8 + 5._r8*lsadz) .and. (cps%snow_depth(c) <= 0.64_r8 + 5._r8*lsadz)) then
                snl(c) = -5
                cps%dz(c,-4) = 0.02_r8 + lsadz
                cps%dz(c,-3) = 0.05_r8 + lsadz
                cps%dz(c,-2) = 0.11_r8 + lsadz
                cps%dz(c,-1) = (cps%snow_depth(c) - cps%dz(c,-4) - cps%dz(c,-3) - cps%dz(c,-2))/2._r8
                cps%dz(c, 0) = cps%dz(c,-1)
             else if (cps%snow_depth(c) > 0.64_r8 + 5._r8*lsadz) then
                snl(c) = -5
                cps%dz(c,-4) = 0.02_r8 + lsadz
                cps%dz(c,-3) = 0.05_r8 + lsadz
                cps%dz(c,-2) = 0.11_r8 + lsadz
                cps%dz(c,-1) = 0.23_r8 + lsadz
                cps%dz(c, 0)=cps%snow_depth(c)-cps%dz(c,-4)-cps%dz(c,-3)-cps%dz(c,-2)-cps%dz(c,-1)
             endif
          end if
          do j = 0, snl(c)+1, -1
             cps%z(c,j)    = cps%zi(c,j) - 0.5_r8*cps%dz(c,j)
             cps%zi(c,j-1) = cps%zi(c,j) - cps%dz(c,j)
          end do
       end if
    end do

    ! Set snow/soil temperature, note:
    ! t_soisno has valid values over non-lakes and S lakes. 
    ! t_lake   has valid values only over lake
    ! t_veg    has valid values over all land

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%lakpoi(l)) then
          if (snl(c) < 0) then
             do j = snl(c)+1, 0
                ces%t_soisno(c,j) = 250._r8
             end do
          end if

          ! Note that ces%t_lake is set in routine initColdDefault()
          ces%t_soisno(c,1:nlevgrnd) = ces%t_lake(c,nlevlak)

          ! Always initialize with no ice to prevent excessive ice sheets from forming when
          ! starting with old lake model that has unrealistically cold lake temperatures.
          ! Keep lake temperature as is, and the energy deficit below freezing (which is no smaller
          ! than it would have been with prognostic ice, as the temperature would then have been higher
          ! and more heat would have flowed out of the lake) will be converted to ice in the first timestep.

          cws%lake_icefrac(c,1:nlevlak) = 0._r8
          cps%savedtke1(c) = tkwat
          cps%ust_lake(c)  = 0.1_r8
          cps%z0mg(c)      = 0.0004_r8
       end if
    end do

    do j = -nlevsno+1, 0
       do c = bounds%begc, bounds%endc
          l = col%landunit(c)
          if (lun%lakpoi(l)) then
             if (j > snl(c)) then
                cws%h2osoi_ice(c,j) = cps%dz(c,j)*bdsno
                cws%h2osoi_liq(c,j) = 0._r8
             end if
          end if
       end do
    end do

    do j = 1,nlevgrnd
       do c = bounds%begc, bounds%endc
          l = col%landunit(c)
          if (lun%lakpoi(l)) then
             if (j <= nlevsoi) then ! soil
                cws%h2osoi_vol(c,j) = cps%watsat(c,j)
                cws%h2osoi_liq(c,j) = spval
                cws%h2osoi_ice(c,j) = spval
             else ! bedrock
                cws%h2osoi_vol(c,j) = 0._r8
             end if
          end if
          ! soil layers
          if (ces%t_soisno(c,j) <= tfrz) then
             cws%h2osoi_ice(c,j) = cps%dz(c,j)*denice*cws%h2osoi_vol(c,j)
             cws%h2osoi_liq(c,j) = 0._r8
          else
             cws%h2osoi_ice(c,j) = 0._r8
             cws%h2osoi_liq(c,j) = cps%dz(c,j)*denh2o*cws%h2osoi_vol(c,j)
          endif
       end do
    end do

    ! SNICAR fields 
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%lakpoi(l) .and. snl(c) < 0) then
          cps%snw_rds(c,snl(c)+1:0) = snw_rds_min
          cps%snw_rds(c,-nlevsno+1:snl(c)) = 0._r8
          cps%snw_rds_top(c) = snw_rds_min
          cps%sno_liq_top(c) = cws%h2osoi_liq(c,snl(c)+1) / &
                              (cws%h2osoi_liq(c,snl(c)+1)+cws%h2osoi_ice(c,snl(c)+1))
       end if
    end do

  end associate 
  end subroutine initColdSLake

end module SLakeInitMod
