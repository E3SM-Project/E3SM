module clm_cpl_indices
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_cpl_indices
!
! !DESCRIPTION:
!    Module containing the indices for the fields passed between CLM and
!    the driver. Includes the River Transport Model fields (RTM) and the
!    fields needed by the land-ice component (sno).
!
! !USES:

  use shr_sys_mod,    only : shr_sys_abort
  implicit none

  SAVE
  private                              ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:

  public :: clm_cpl_indices_set        ! Set the coupler indices

!
! !PUBLIC DATA MEMBERS:
!
  integer , public :: glc_nec     ! number of elevation classes for glacier_mec landunits 
                                  ! (from coupler) - must equal maxpatch_glcmec from namelist
  integer , parameter, private:: glc_nec_max = 100

  ! lnd -> drv (required)

  integer, public ::index_l2x_Flrl_rofliq     ! lnd->rtm input fluxes
  integer, public ::index_l2x_Flrl_rofice     ! lnd->rtm input fluxes

  integer, public ::index_l2x_Sl_t            ! temperature
  integer, public ::index_l2x_Sl_tref         ! 2m reference temperature
  integer, public ::index_l2x_Sl_qref         ! 2m reference specific humidity
  integer, public ::index_l2x_Sl_avsdr        ! albedo: direct , visible
  integer, public ::index_l2x_Sl_anidr        ! albedo: direct , near-ir
  integer, public ::index_l2x_Sl_avsdf        ! albedo: diffuse, visible
  integer, public ::index_l2x_Sl_anidf        ! albedo: diffuse, near-ir
  integer, public ::index_l2x_Sl_snowh        ! snow height
  integer, public ::index_l2x_Sl_u10          ! 10m wind
  integer, public ::index_l2x_Sl_ddvel        ! dry deposition velocities (optional)
  integer, public ::index_l2x_Sl_fv           ! friction velocity  
  integer, public ::index_l2x_Sl_ram1         ! aerodynamical resistance
  integer, public ::index_l2x_Sl_soilw        ! volumetric soil water
  integer, public ::index_l2x_Fall_taux       ! wind stress, zonal
  integer, public ::index_l2x_Fall_tauy       ! wind stress, meridional
  integer, public ::index_l2x_Fall_lat        ! latent          heat flux
  integer, public ::index_l2x_Fall_sen        ! sensible        heat flux
  integer, public ::index_l2x_Fall_lwup       ! upward longwave heat flux
  integer, public ::index_l2x_Fall_evap       ! evaporation     water flux
  integer, public ::index_l2x_Fall_swnet      ! heat flux       shortwave net       
  integer, public ::index_l2x_Fall_fco2_lnd   ! co2 flux **For testing set to 0
  integer, public ::index_l2x_Fall_flxdst1    ! dust flux size bin 1    
  integer, public ::index_l2x_Fall_flxdst2    ! dust flux size bin 2    
  integer, public ::index_l2x_Fall_flxdst3    ! dust flux size bin 3    
  integer, public ::index_l2x_Fall_flxdst4    ! dust flux size bin 4
  integer, public ::index_l2x_Fall_flxvoc     ! MEGAN fluxes

#ifdef LCH4
  integer, public ::index_x2l_Sa_methane
  integer, public ::index_l2x_Fall_methane
#endif

  integer, public :: nflds_l2x = 0

  ! drv -> lnd (required)

  integer, public ::index_x2l_Sa_z            ! bottom atm level height
  integer, public ::index_x2l_Sa_u            ! bottom atm level zon wind
  integer, public ::index_x2l_Sa_v            ! bottom atm level mer wind
  integer, public ::index_x2l_Sa_ptem         ! bottom atm level pot temp
  integer, public ::index_x2l_Sa_shum         ! bottom atm level spec hum
  integer, public ::index_x2l_Sa_pbot         ! bottom atm level pressure
  integer, public ::index_x2l_Sa_tbot         ! bottom atm level temp
  integer, public ::index_x2l_Faxa_lwdn       ! downward lw heat flux
  integer, public ::index_x2l_Faxa_rainc      ! prec: liquid "convective"
  integer, public ::index_x2l_Faxa_rainl      ! prec: liquid "large scale"
  integer, public ::index_x2l_Faxa_snowc      ! prec: frozen "convective"
  integer, public ::index_x2l_Faxa_snowl      ! prec: frozen "large scale"
  integer, public ::index_x2l_Faxa_swndr      ! sw: nir direct  downward
  integer, public ::index_x2l_Faxa_swvdr      ! sw: vis direct  downward
  integer, public ::index_x2l_Faxa_swndf      ! sw: nir diffuse downward
  integer, public ::index_x2l_Faxa_swvdf      ! sw: vis diffuse downward
  integer, public ::index_x2l_Sa_co2prog      ! bottom atm level prognostic co2
  integer, public ::index_x2l_Sa_co2diag      ! bottom atm level diagnostic co2
  integer, public ::index_x2l_Faxa_bcphidry   ! flux: Black Carbon hydrophilic dry deposition
  integer, public ::index_x2l_Faxa_bcphodry   ! flux: Black Carbon hydrophobic dry deposition
  integer, public ::index_x2l_Faxa_bcphiwet   ! flux: Black Carbon hydrophilic wet deposition
  integer, public ::index_x2l_Faxa_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer, public ::index_x2l_Faxa_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer, public ::index_x2l_Faxa_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer, public ::index_x2l_Faxa_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer, public ::index_x2l_Faxa_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer, public ::index_x2l_Faxa_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer, public ::index_x2l_Faxa_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer, public ::index_x2l_Faxa_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer, public ::index_x2l_Faxa_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer, public ::index_x2l_Faxa_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer, public ::index_x2l_Faxa_dstdry4    ! flux: Size 4 dust -- dry deposition
 
  integer, public ::index_x2l_Flrr_flood      ! rtm->lnd rof (flood) flux

  integer, public ::index_x2l_Slrr_volr      ! rtm->lnd rof volr
  integer, public :: nflds_x2l = 0

  ! sno -> drv (only if land-ice model is NOT a stub model)

  integer, public :: index_s2x_Ss_tsrf(glc_nec_max)   = 0 ! glc MEC temperature
  integer, public :: index_s2x_Ss_topo(glc_nec_max)   = 0 ! glc MEC topo height
  integer, public :: index_s2x_Fgss_qice(glc_nec_max) = 0 ! glc MEC ice flux

  integer, public :: nflds_s2x = 0

  ! drv -> sno (only if land-ice model is NOT a stub model)

  integer, public :: index_x2s_Sg_frac(glc_nec_max)   = 0   ! Fraction of glacier in glc MEC class 1
  integer, public :: index_x2s_Sg_topo(glc_nec_max)   = 0   ! Topo height in glc MEC class 1
  integer, public :: index_x2s_Fsgg_rofi(glc_nec_max) = 0   ! Ice runoff from glc model
  integer, public :: index_x2s_Fsgg_rofl(glc_nec_max) = 0   ! Liquid runoff from glc model
  integer, public :: index_x2s_Fsgg_hflx(glc_nec_max) = 0

  integer, public :: nflds_x2s = 0
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 01/19/2011, Erik Kluzek:         Added protex headers
! 01/20/2011, Erik Kluzek:         Set nflds
! 12/11/2011, Mariana Vertenstein  Updated for new coupler fields spec
!
!EOP
!-----------------------------------------------------------------------

!=======================================================================
contains
!=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_cpl_indices_set
!
! !CALLED FROM: lnd_comp_mct or lnd_comp_esmf
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
! !INTERFACE:
  subroutine clm_cpl_indices_set( )

!
! !DESCRIPTION: 
! Set the coupler indices needed by the land model coupler
! interface.
!
! !USES:
  use seq_flds_mod  , only: seq_flds_x2l_fields, seq_flds_l2x_fields,     &
                            seq_flds_x2s_fields, seq_flds_s2x_fields
  use mct_mod       , only: mct_aVect, mct_aVect_init, mct_avect_indexra, &
                            mct_aVect_clean, mct_avect_nRattr
  use seq_drydep_mod, only: drydep_fields_token, lnd_drydep
  use shr_megan_mod,  only: shr_megan_fields_token, shr_megan_mechcomps_n

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! 01/2011, Erik Kluzek:         Added protex headers
!
! !LOCAL VARIABLES:
    type(mct_aVect)   :: l2x      ! temporary, land to coupler
    type(mct_aVect)   :: x2l      ! temporary, coupler to land
    type(mct_aVect)   :: s2x      ! temporary, glacier to coupler
    type(mct_aVect)   :: x2s      ! temporary, coupler to glacier
    integer           :: num 
    character(len= 2) :: cnum
    character(len=64) :: name
    character(len=32) :: subname = 'clm_cpl_indices_set'  ! subroutine name
!EOP
!
!-----------------------------------------------------------------------

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2l, rList=seq_flds_x2l_fields, lsize=1)
    call mct_aVect_init(l2x, rList=seq_flds_l2x_fields, lsize=1)

    !-------------------------------------------------------------
    ! clm -> drv 
    !-------------------------------------------------------------

    index_l2x_Flrl_rofliq   = mct_avect_indexra(l2x,'Flrl_rofliq')
    index_l2x_Flrl_rofice   = mct_avect_indexra(l2x,'Flrl_rofice')

    index_l2x_Sl_t          = mct_avect_indexra(l2x,'Sl_t')
    index_l2x_Sl_snowh      = mct_avect_indexra(l2x,'Sl_snowh')
    index_l2x_Sl_avsdr      = mct_avect_indexra(l2x,'Sl_avsdr')
    index_l2x_Sl_anidr      = mct_avect_indexra(l2x,'Sl_anidr')
    index_l2x_Sl_avsdf      = mct_avect_indexra(l2x,'Sl_avsdf')
    index_l2x_Sl_anidf      = mct_avect_indexra(l2x,'Sl_anidf')
    index_l2x_Sl_tref       = mct_avect_indexra(l2x,'Sl_tref')
    index_l2x_Sl_qref       = mct_avect_indexra(l2x,'Sl_qref')
    index_l2x_Sl_u10        = mct_avect_indexra(l2x,'Sl_u10')
    index_l2x_Sl_ram1       = mct_avect_indexra(l2x,'Sl_ram1')
    index_l2x_Sl_fv         = mct_avect_indexra(l2x,'Sl_fv')
    index_l2x_Sl_soilw      = mct_avect_indexra(l2x,'Sl_soilw',perrwith='quiet')
    if ( lnd_drydep )then
       index_l2x_Sl_ddvel = mct_avect_indexra(l2x, trim(drydep_fields_token))
    else
       index_l2x_Sl_ddvel = 0
    end if

    index_l2x_Fall_taux     = mct_avect_indexra(l2x,'Fall_taux')
    index_l2x_Fall_tauy     = mct_avect_indexra(l2x,'Fall_tauy')
    index_l2x_Fall_lat      = mct_avect_indexra(l2x,'Fall_lat')
    index_l2x_Fall_sen      = mct_avect_indexra(l2x,'Fall_sen')
    index_l2x_Fall_lwup     = mct_avect_indexra(l2x,'Fall_lwup')
    index_l2x_Fall_evap     = mct_avect_indexra(l2x,'Fall_evap')
    index_l2x_Fall_swnet    = mct_avect_indexra(l2x,'Fall_swnet')
    index_l2x_Fall_flxdst1  = mct_avect_indexra(l2x,'Fall_flxdst1')
    index_l2x_Fall_flxdst2  = mct_avect_indexra(l2x,'Fall_flxdst2')
    index_l2x_Fall_flxdst3  = mct_avect_indexra(l2x,'Fall_flxdst3')
    index_l2x_Fall_flxdst4  = mct_avect_indexra(l2x,'Fall_flxdst4')

    index_l2x_Fall_fco2_lnd = mct_avect_indexra(l2x,'Fall_fco2_lnd',perrwith='quiet')

#ifdef LCH4
    index_l2x_Fall_methane  = mct_avect_indexra(l2x,'Fall_methane',perrWith='quiet')
#endif

    ! MEGAN fluxes
    if (shr_megan_mechcomps_n>0) then
       index_l2x_Fall_flxvoc = mct_avect_indexra(l2x,trim(shr_megan_fields_token))
    else
       index_l2x_Fall_flxvoc = 0
    endif

    nflds_l2x = mct_avect_nRattr(l2x)

    !-------------------------------------------------------------
    ! drv -> clm
    !-------------------------------------------------------------

    index_x2l_Sa_z          = mct_avect_indexra(x2l,'Sa_z')
    index_x2l_Sa_u          = mct_avect_indexra(x2l,'Sa_u')
    index_x2l_Sa_v          = mct_avect_indexra(x2l,'Sa_v')
    index_x2l_Sa_ptem       = mct_avect_indexra(x2l,'Sa_ptem')
    index_x2l_Sa_pbot       = mct_avect_indexra(x2l,'Sa_pbot')
    index_x2l_Sa_tbot       = mct_avect_indexra(x2l,'Sa_tbot')
    index_x2l_Sa_shum       = mct_avect_indexra(x2l,'Sa_shum')
    index_x2l_Sa_co2prog    = mct_avect_indexra(x2l,'Sa_co2prog',perrwith='quiet')
    index_x2l_Sa_co2diag    = mct_avect_indexra(x2l,'Sa_co2diag',perrwith='quiet')

#ifdef LCH4
    index_x2l_Sa_methane    = mct_avect_indexra(x2l,'Sa_methane',perrWith='quiet')
#endif


    index_x2l_Slrr_volr       = mct_avect_indexra(x2l,'Slrr_volr')

    index_x2l_Faxa_lwdn     = mct_avect_indexra(x2l,'Faxa_lwdn')
    index_x2l_Faxa_rainc    = mct_avect_indexra(x2l,'Faxa_rainc')
    index_x2l_Faxa_rainl    = mct_avect_indexra(x2l,'Faxa_rainl')
    index_x2l_Faxa_snowc    = mct_avect_indexra(x2l,'Faxa_snowc')
    index_x2l_Faxa_snowl    = mct_avect_indexra(x2l,'Faxa_snowl')
    index_x2l_Faxa_swndr    = mct_avect_indexra(x2l,'Faxa_swndr')
    index_x2l_Faxa_swvdr    = mct_avect_indexra(x2l,'Faxa_swvdr')
    index_x2l_Faxa_swndf    = mct_avect_indexra(x2l,'Faxa_swndf')
    index_x2l_Faxa_swvdf    = mct_avect_indexra(x2l,'Faxa_swvdf')
    index_x2l_Faxa_bcphidry = mct_avect_indexra(x2l,'Faxa_bcphidry')
    index_x2l_Faxa_bcphodry = mct_avect_indexra(x2l,'Faxa_bcphodry')
    index_x2l_Faxa_bcphiwet = mct_avect_indexra(x2l,'Faxa_bcphiwet')
    index_x2l_Faxa_ocphidry = mct_avect_indexra(x2l,'Faxa_ocphidry')
    index_x2l_Faxa_ocphodry = mct_avect_indexra(x2l,'Faxa_ocphodry')
    index_x2l_Faxa_ocphiwet = mct_avect_indexra(x2l,'Faxa_ocphiwet')
    index_x2l_Faxa_dstdry1  = mct_avect_indexra(x2l,'Faxa_dstdry1')
    index_x2l_Faxa_dstdry2  = mct_avect_indexra(x2l,'Faxa_dstdry2')
    index_x2l_Faxa_dstdry3  = mct_avect_indexra(x2l,'Faxa_dstdry3')
    index_x2l_Faxa_dstdry4  = mct_avect_indexra(x2l,'Faxa_dstdry4')
    index_x2l_Faxa_dstwet1  = mct_avect_indexra(x2l,'Faxa_dstwet1')
    index_x2l_Faxa_dstwet2  = mct_avect_indexra(x2l,'Faxa_dstwet2')
    index_x2l_Faxa_dstwet3  = mct_avect_indexra(x2l,'Faxa_dstwet3')
    index_x2l_Faxa_dstwet4  = mct_avect_indexra(x2l,'Faxa_dstwet4')

    index_x2l_Flrr_flood    = mct_avect_indexra(x2l,'Flrr_flood')

    nflds_x2l = mct_avect_nRattr(x2l)

    call mct_aVect_clean(x2l)
    call mct_aVect_clean(l2x)

    !-------------------------------------------------------------
    ! drv->sno (for cism coupling)
    !-------------------------------------------------------------

    glc_nec = 0

    if (seq_flds_x2s_fields /= ' ') then
       call mct_aVect_init(x2s, rList=seq_flds_x2s_fields, lsize=1)

       do num = 1,glc_nec_max
          write(cnum,'(i2.2)') num
          name = 'Sg_frac' // cnum
          index_x2s_Sg_frac(num)   = mct_avect_indexra(x2s,trim(name),perrwith='quiet') 
          name = 'Sg_topo' // cnum
          index_x2s_Sg_topo(num)   = mct_avect_indexra(x2s,trim(name),perrwith='quiet')
          name = 'Fsgg_rofi' // cnum
          index_x2s_Fsgg_rofi(num) = mct_avect_indexra(x2s,trim(name),perrwith='quiet')
          name = 'Fsgg_rofl' // cnum
          index_x2s_Fsgg_rofl(num) = mct_avect_indexra(x2s,trim(name),perrwith='quiet')
          name = 'Fsgg_hflx' // cnum
          index_x2s_Fsgg_hflx(num) = mct_avect_indexra(x2s,trim(name),perrwith='quiet')
          if ( index_x2s_Sg_frac(num)   == 0 .and. &
               index_x2s_Sg_topo(num)   == 0 .and. &
               index_x2s_Fsgg_rofi(num) == 0 .and. &
               index_x2s_fsgg_rofl(num) == 0 .and. &
               index_x2s_Fsgg_hflx(num) == 0 ) then
             exit
          end if
          glc_nec = num
       end do
       if (glc_nec == glc_nec_max) then
          call shr_sys_abort (subname // 'error: glc_nec_cpl cannot equal glc_nec_max')
       end if

       nflds_x2s = mct_avect_nRattr(x2s)
       call mct_aVect_clean(x2s)
    end if

    !-------------------------------------------------------------
    ! sno -> drv (for cism coupling)
    !-------------------------------------------------------------

    if (seq_flds_s2x_fields /= ' ') then
       call mct_aVect_init(s2x, rList=seq_flds_s2x_fields, lsize=1)

       do num = 1,glc_nec
          write(cnum,'(i2.2)') num

          name = 'Ss_tsrf' // cnum
          index_s2x_Ss_tsrf(num)   = mct_avect_indexra(s2x,trim(name))
          name = 'Ss_topo' // cnum
          index_s2x_Ss_topo(num)   = mct_avect_indexra(s2x,trim(name))
          name = 'Fgss_qice' // cnum
          index_s2x_Fgss_qice(num) = mct_avect_indexra(s2x,trim(name))
       end do

       nflds_s2x = mct_avect_nRattr(s2x)
       call mct_aVect_clean(s2x)
    end if

  end subroutine clm_cpl_indices_set

!=======================================================================

end module clm_cpl_indices
