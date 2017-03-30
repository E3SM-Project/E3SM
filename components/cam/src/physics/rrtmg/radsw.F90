
module radsw
!----------------------------------------------------------------------- 
! 
! Purpose: Solar radiation calculations.
!
!-----------------------------------------------------------------------
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use cam_abortutils,      only: endrun
use cam_history,     only: outfld
use scamMod,         only: single_column,scm_crm_mode,have_asdir, &
                           asdirobs, have_asdif, asdifobs, have_aldir, &
                           aldirobs, have_aldif, aldifobs
use cam_logfile,     only: iulog
use parrrsw,         only: nbndsw, ngptsw
use rrtmg_sw_init,   only: rrtmg_sw_ini
use rrtmg_sw_rad,    only: rrtmg_sw
use perf_mod,        only: t_startf, t_stopf
use radconstants,    only: idx_sw_diag

implicit none

private
save

real(r8) :: fractional_solar_irradiance(1:nbndsw) ! fraction of solar irradiance in each band
real(r8) :: solar_band_irrad(1:nbndsw) ! rrtmg-assumed solar irradiance in each sw band

logical :: pergro_mods = .false.

! Public methods

public ::&
   radsw_init,      &! initialize constants
   rad_rrtmg_sw      ! driver for solar radiation code

!===============================================================================
CONTAINS
!===============================================================================

subroutine rad_rrtmg_sw(lchnk,ncol       ,rrtmg_levs   ,r_state      , &
                    E_pmid   ,E_cld      ,                             &
                    E_aer_tau,E_aer_tau_w,E_aer_tau_w_g,E_aer_tau_w_f, &
                    eccf     ,E_coszrs   ,solin        ,sfac         , &
                    E_asdir  ,E_asdif    ,E_aldir      ,E_aldif      , &
                    qrs      ,qrsc       ,fsnt         ,fsntc        ,fsntoa,fsutoa, &
                    fsntoac  ,fsnirtoa   ,fsnrtoac     ,fsnrtoaq     ,fsns    , &
                    fsnsc    ,fsdsc      ,fsds         ,sols         ,soll    , &
                    solsd    ,solld      ,fns          ,fcns         , &
                    Nday     ,Nnite      ,IdxDay       ,IdxNite      ,E_rngsw , &
                    su       ,sd         ,                             &
                    E_cld_tau, E_cld_tau_w, E_cld_tau_w_g, E_cld_tau_w_f,  &
                    old_convert)


!-----------------------------------------------------------------------
! 
! Purpose: 
! Solar radiation code
! 
! Method: 
! mji/rrtmg
! RRTMG, two-stream, with McICA
! 
! Divides solar spectrum into 14 intervals from 0.2-12.2 micro-meters.
! solar flux fractions specified for each interval. allows for
! seasonally and diurnally varying solar input.  Includes molecular,
! cloud, aerosol, and surface scattering, along with h2o,o3,co2,o2,cloud, 
! and surface absorption. Computes delta-eddington reflections and
! transmissions assuming homogeneously mixed layers. Adds the layers 
! assuming scattering between layers to be isotropic, and distinguishes 
! direct solar beam from scattered radiation.
! 
! Longitude loops are broken into 1 or 2 sections, so that only daylight
! (i.e. coszrs > 0) computations are done.
! 
! Note that an extra layer above the model top layer is added.
! 
! mks units are used.
! 
! Special diagnostic calculation of the clear sky surface and total column
! absorbed flux is also done for cloud forcing diagnostics.
! 
!-----------------------------------------------------------------------

   use cmparray_mod,        only: CmpDayNite, ExpDayNite
   use phys_control,        only: phys_getopts
   use mcica_subcol_gen_sw, only: mcica_subcol_sw
   use physconst,           only: cpair
   use rrtmg_state,         only: rrtmg_state_t
   
   ! Minimum cloud amount (as a fraction of the grid-box area) to 
   ! distinguish from clear sky
   real(r8), parameter :: cldmin = 1.0e-80_r8

   ! Decimal precision of cloud amount (0 -> preserve full resolution;
   ! 10^-n -> preserve n digits of cloud amount)
   real(r8), parameter :: cldeps = 0.0_r8

   ! Input arguments
   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: ncol              ! number of atmospheric columns
   integer, intent(in) :: rrtmg_levs        ! number of levels rad is applied

    type(rrtmg_state_t), intent(in) :: r_state

   integer, intent(in) :: Nday                      ! Number of daylight columns
   integer, intent(in) :: Nnite                     ! Number of night columns
   integer, intent(in), dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
   integer, intent(in), dimension(pcols) :: IdxNite ! Indicies of night coumns

   real(r8), intent(in) :: E_pmid(pcols,pver)  ! Level pressure (Pascals)
   real(r8), intent(in) :: E_cld(pcols,pver)    ! Fractional cloud cover

   real(r8), intent(in) :: E_aer_tau    (pcols, 0:pver, nbndsw)      ! aerosol optical depth
   real(r8), intent(in) :: E_aer_tau_w  (pcols, 0:pver, nbndsw)      ! aerosol OD * ssa
   real(r8), intent(in) :: E_aer_tau_w_g(pcols, 0:pver, nbndsw)      ! aerosol OD * ssa * asm
   real(r8), intent(in) :: E_aer_tau_w_f(pcols, 0:pver, nbndsw)      ! aerosol OD * ssa * fwd

   real(r8), intent(in) :: eccf               ! Eccentricity factor (1./earth-sun dist^2)
   real(r8), intent(in) :: E_coszrs(pcols)    ! Cosine solar zenith angle
   real(r8), intent(in) :: E_asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: E_aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: E_asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
   real(r8), intent(in) :: E_aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad
   real(r8), intent(in) :: sfac(nbndsw)            ! factor to account for solar variability in each band 
   real(r8), intent(in) :: E_rngsw(ngptsw,pcols,pver)            ! rand # for sw!BALLI

   real(r8), optional, intent(in) :: E_cld_tau    (nbndsw, pcols, pver)      ! cloud optical depth
   real(r8), optional, intent(in) :: E_cld_tau_w  (nbndsw, pcols, pver)      ! cloud optical 
   real(r8), optional, intent(in) :: E_cld_tau_w_g(nbndsw, pcols, pver)      ! cloud optical 
   real(r8), optional, intent(in) :: E_cld_tau_w_f(nbndsw, pcols, pver)      ! cloud optical 
   logical, optional, intent(in) :: old_convert

   ! Output arguments

   real(r8), intent(out) :: solin(pcols)     ! Incident solar flux
   real(r8), intent(out) :: qrs (pcols,pver) ! Solar heating rate
   real(r8), intent(out) :: qrsc(pcols,pver) ! Clearsky solar heating rate
   real(r8), intent(out) :: fsns(pcols)      ! Surface absorbed solar flux
   real(r8), intent(out) :: fsnt(pcols)      ! Total column absorbed solar flux
   real(r8), intent(out) :: fsntoa(pcols)    ! Net solar flux at TOA
   real(r8), intent(out) :: fsutoa(pcols)    ! Upward solar flux at TOA
   real(r8), intent(out) :: fsds(pcols)      ! Flux shortwave downwelling surface

   real(r8), intent(out) :: fsnsc(pcols)     ! Clear sky surface absorbed solar flux
   real(r8), intent(out) :: fsdsc(pcols)     ! Clear sky surface downwelling solar flux
   real(r8), intent(out) :: fsntc(pcols)     ! Clear sky total column absorbed solar flx
   real(r8), intent(out) :: fsntoac(pcols)   ! Clear sky net solar flx at TOA
   real(r8), intent(out) :: sols(pcols)      ! Direct solar rad on surface (< 0.7)
   real(r8), intent(out) :: soll(pcols)      ! Direct solar rad on surface (>= 0.7)
   real(r8), intent(out) :: solsd(pcols)     ! Diffuse solar rad on surface (< 0.7)
   real(r8), intent(out) :: solld(pcols)     ! Diffuse solar rad on surface (>= 0.7)
   real(r8), intent(out) :: fsnirtoa(pcols)  ! Near-IR flux absorbed at toa
   real(r8), intent(out) :: fsnrtoac(pcols)  ! Clear sky near-IR flux absorbed at toa
   real(r8), intent(out) :: fsnrtoaq(pcols)  ! Net near-IR flux at toa >= 0.7 microns

   real(r8), intent(out) :: fns(pcols,pverp)   ! net flux at interfaces
   real(r8), intent(out) :: fcns(pcols,pverp)  ! net clear-sky flux at interfaces

   real(r8), pointer, dimension(:,:,:) :: su ! shortwave spectral flux up
   real(r8), pointer, dimension(:,:,:) :: sd ! shortwave spectral flux down

   !---------------------------Local variables-----------------------------

   ! Local and reordered copies of the intent(in) variables
   real(r8) :: rngsw(ngptsw,pcols,pver)            ! rand # for sw!BALLI
   real(r8) :: pmid(pcols,pver)    ! Level pressure (Pascals)

   real(r8) :: cld(pcols,rrtmg_levs-1)    ! Fractional cloud cover
   real(r8) :: cicewp(pcols,rrtmg_levs-1) ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,rrtmg_levs-1) ! in-cloud cloud liquid water path
   real(r8) :: rel(pcols,rrtmg_levs-1)    ! Liquid effective drop size (microns)
   real(r8) :: rei(pcols,rrtmg_levs-1)    ! Ice effective drop size (microns)

   real(r8) :: coszrs(pcols)    ! Cosine solar zenith angle
   real(r8) :: asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
   real(r8) :: aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
   real(r8) :: asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
   real(r8) :: aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad

   real(r8) :: h2ovmr(pcols,rrtmg_levs)   ! h2o volume mixing ratio
   real(r8) :: o3vmr(pcols,rrtmg_levs)    ! o3 volume mixing ratio
   real(r8) :: co2vmr(pcols,rrtmg_levs)   ! co2 volume mixing ratio 
   real(r8) :: ch4vmr(pcols,rrtmg_levs)   ! ch4 volume mixing ratio 
   real(r8) :: o2vmr(pcols,rrtmg_levs)    ! o2  volume mixing ratio 
   real(r8) :: n2ovmr(pcols,rrtmg_levs)   ! n2o volume mixing ratio 

   real(r8) :: tsfc(pcols)          ! surface temperature

   integer :: inflgsw               ! flag for cloud parameterization method
   integer :: iceflgsw              ! flag for ice cloud parameterization method
   integer :: liqflgsw              ! flag for liquid cloud parameterization method
   integer :: icld                  ! Flag for cloud overlap method
                                    ! 0=clear, 1=random, 2=maximum/random, 3=maximum
   integer :: dyofyr                ! Set to day of year for Earth/Sun distance calculation in
                                    ! rrtmg_sw, or pass in adjustment directly into adjes
   real(r8) :: solvar(nbndsw)       ! solar irradiance variability in each band

   integer, parameter :: nsubcsw = ngptsw           ! rrtmg_sw g-point (quadrature point) dimension
   integer :: permuteseed                           ! permute seed for sub-column generator

   real(r8) :: diagnostic_od(pcols, pver)           ! cloud optical depth - diagnostic temp variable

   real(r8) :: tauc_sw(nbndsw, pcols, rrtmg_levs-1)         ! cloud optical depth
   real(r8) :: ssac_sw(nbndsw, pcols, rrtmg_levs-1)         ! cloud single scat. albedo
   real(r8) :: asmc_sw(nbndsw, pcols, rrtmg_levs-1)         ! cloud asymmetry parameter
   real(r8) :: fsfc_sw(nbndsw, pcols, rrtmg_levs-1)         ! cloud forward scattering fraction

   real(r8) :: tau_aer_sw(pcols, rrtmg_levs-1, nbndsw)      ! aer optical depth
   real(r8) :: ssa_aer_sw(pcols, rrtmg_levs-1, nbndsw)      ! aer single scat. albedo
   real(r8) :: asm_aer_sw(pcols, rrtmg_levs-1, nbndsw)      ! aer asymmetry parameter

   real(r8) :: cld_stosw(nsubcsw, pcols, rrtmg_levs-1)      ! stochastic cloud fraction
   real(r8) :: rei_stosw(pcols, rrtmg_levs-1)               ! stochastic ice particle size 
   real(r8) :: rel_stosw(pcols, rrtmg_levs-1)               ! stochastic liquid particle size
   real(r8) :: cicewp_stosw(nsubcsw, pcols, rrtmg_levs-1)   ! stochastic cloud ice water path
   real(r8) :: cliqwp_stosw(nsubcsw, pcols, rrtmg_levs-1)   ! stochastic cloud liquid wter path
   real(r8) :: tauc_stosw(nsubcsw, pcols, rrtmg_levs-1)     ! stochastic cloud optical depth (optional)
   real(r8) :: ssac_stosw(nsubcsw, pcols, rrtmg_levs-1)     ! stochastic cloud single scat. albedo (optional)
   real(r8) :: asmc_stosw(nsubcsw, pcols, rrtmg_levs-1)     ! stochastic cloud asymmetry parameter (optional)
   real(r8) :: fsfc_stosw(nsubcsw, pcols, rrtmg_levs-1)     ! stochastic cloud forward scattering fraction (optional)

   real(r8), parameter :: dps = 1._r8/86400._r8 ! Inverse of seconds per day
 
   real(r8) :: swuflx(pcols,rrtmg_levs+1)       ! Total sky shortwave upward flux (W/m2)
   real(r8) :: swdflx(pcols,rrtmg_levs+1)       ! Total sky shortwave downward flux (W/m2)
   real(r8) :: swhr(pcols,rrtmg_levs)           ! Total sky shortwave radiative heating rate (K/d)
   real(r8) :: swuflxc(pcols,rrtmg_levs+1)      ! Clear sky shortwave upward flux (W/m2)
   real(r8) :: swdflxc(pcols,rrtmg_levs+1)      ! Clear sky shortwave downward flux (W/m2)
   real(r8) :: swhrc(pcols,rrtmg_levs)          ! Clear sky shortwave radiative heating rate (K/d)
   real(r8) :: swuflxs(nbndsw,pcols,rrtmg_levs+1)  ! Shortwave spectral flux up
   real(r8) :: swdflxs(nbndsw,pcols,rrtmg_levs+1)  ! Shortwave spectral flux down

   real(r8) :: dirdnuv(pcols,rrtmg_levs+1)       ! Direct downward shortwave flux, UV/vis
   real(r8) :: difdnuv(pcols,rrtmg_levs+1)       ! Diffuse downward shortwave flux, UV/vis
   real(r8) :: dirdnir(pcols,rrtmg_levs+1)       ! Direct downward shortwave flux, near-IR
   real(r8) :: difdnir(pcols,rrtmg_levs+1)       ! Diffuse downward shortwave flux, near-IR

   ! Added for net near-IR diagnostic
   real(r8) :: ninflx(pcols,rrtmg_levs+1)        ! Net shortwave flux, near-IR
   real(r8) :: ninflxc(pcols,rrtmg_levs+1)       ! Net clear sky shortwave flux, near-IR

   ! Other

   integer :: i, k, ns       ! indices

   ! Cloud radiative property arrays
   real(r8) :: tauxcl(pcols,0:pver) ! water cloud extinction optical depth
   real(r8) :: tauxci(pcols,0:pver) ! ice cloud extinction optical depth
   real(r8) :: wcl(pcols,0:pver) ! liquid cloud single scattering albedo
   real(r8) :: gcl(pcols,0:pver) ! liquid cloud asymmetry parameter
   real(r8) :: fcl(pcols,0:pver) ! liquid cloud forward scattered fraction
   real(r8) :: wci(pcols,0:pver) ! ice cloud single scattering albedo
   real(r8) :: gci(pcols,0:pver) ! ice cloud asymmetry parameter
   real(r8) :: fci(pcols,0:pver) ! ice cloud forward scattered fraction

   ! Aerosol radiative property arrays
   real(r8) :: tauxar(pcols,0:pver) ! aerosol extinction optical depth
   real(r8) :: wa(pcols,0:pver) ! aerosol single scattering albedo
   real(r8) :: ga(pcols,0:pver) ! aerosol assymetry parameter
   real(r8) :: fa(pcols,0:pver) ! aerosol forward scattered fraction

   ! CRM
   real(r8) :: fus(pcols,pverp)   ! Upward flux (added for CRM)
   real(r8) :: fds(pcols,pverp)   ! Downward flux (added for CRM)
   real(r8) :: fusc(pcols,pverp)  ! Upward clear-sky flux (added for CRM)
   real(r8) :: fdsc(pcols,pverp)  ! Downward clear-sky flux (added for CRM)

   integer :: kk

   real(r8) :: pmidmb(pcols,rrtmg_levs)   ! Level pressure (hPa)
   real(r8) :: pintmb(pcols,rrtmg_levs+1) ! Model interface pressure (hPa)
   real(r8) :: tlay(pcols,rrtmg_levs)     ! mid point temperature
   real(r8) :: tlev(pcols,rrtmg_levs+1)   ! interface temperature

   integer :: isub!BALLI
   !-----------------------------------------------------------------------
   ! START OF CALCULATION
   !-----------------------------------------------------------------------

   ! Initialize output fields:

   fsds(1:ncol)     = 0.0_r8

   fsnirtoa(1:ncol) = 0.0_r8
   fsnrtoac(1:ncol) = 0.0_r8
   fsnrtoaq(1:ncol) = 0.0_r8

   fsns(1:ncol)     = 0.0_r8
   fsnsc(1:ncol)    = 0.0_r8
   fsdsc(1:ncol)    = 0.0_r8

   fsnt(1:ncol)     = 0.0_r8
   fsntc(1:ncol)    = 0.0_r8
   fsntoa(1:ncol)   = 0.0_r8
   fsutoa(1:ncol)   = 0.0_r8
   fsntoac(1:ncol)  = 0.0_r8

   solin(1:ncol)    = 0.0_r8

   sols(1:ncol)     = 0.0_r8
   soll(1:ncol)     = 0.0_r8
   solsd(1:ncol)    = 0.0_r8
   solld(1:ncol)    = 0.0_r8

   qrs (1:ncol,1:pver) = 0.0_r8
   qrsc(1:ncol,1:pver) = 0.0_r8
   fns(1:ncol,1:pverp) = 0.0_r8
   fcns(1:ncol,1:pverp) = 0.0_r8
   if (single_column.and.scm_crm_mode) then 
      fus(1:ncol,1:pverp) = 0.0_r8
      fds(1:ncol,1:pverp) = 0.0_r8
      fusc(:ncol,:pverp) = 0.0_r8
      fdsc(:ncol,:pverp) = 0.0_r8
   endif

   if (associated(su)) su(1:ncol,:,:) = 0.0_r8
   if (associated(sd)) sd(1:ncol,:,:) = 0.0_r8

   ! If night everywhere, return:
   if ( Nday == 0 ) then
     return
   endif

   ! Rearrange input arrays
   call CmpDayNite(E_pmid(:,pverp-rrtmg_levs+1:pver), pmid(:,1:rrtmg_levs-1), &
        Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs-1)
   call CmpDayNite(E_cld(:,pverp-rrtmg_levs+1:pver),  cld(:,1:rrtmg_levs-1), &
        Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs-1)

   call CmpDayNite(r_state%pintmb, pintmb, Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs+1)
   call CmpDayNite(r_state%pmidmb, pmidmb, Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs)
   call CmpDayNite(r_state%h2ovmr, h2ovmr, Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs)
   call CmpDayNite(r_state%o3vmr,  o3vmr,  Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs)
   call CmpDayNite(r_state%co2vmr, co2vmr, Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs)

   call CmpDayNite(E_coszrs, coszrs,    Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_asdir,  asdir,     Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_aldir,  aldir,     Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_asdif,  asdif,     Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_aldif,  aldif,     Nday, IdxDay, Nnite, IdxNite, 1, pcols)

   call CmpDayNite(r_state%tlay,   tlay,   Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs)
   call CmpDayNite(r_state%tlev,   tlev,   Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs+1)
   call CmpDayNite(r_state%ch4vmr, ch4vmr, Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs)
   call CmpDayNite(r_state%o2vmr,  o2vmr,  Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs)
   call CmpDayNite(r_state%n2ovmr, n2ovmr, Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, rrtmg_levs)
   !BSINGH - rearrange random numbers
   do isub = 1, ngptsw
      call CmpDayNite(E_rngsw(isub,:,:), rngsw(isub,:,:),Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)      
   enddo


   ! These fields are no longer input by CAM.
   cicewp = 0.0_r8
   cliqwp = 0.0_r8
   rel = 0.0_r8
   rei = 0.0_r8

   ! Aerosol daylight map
   ! Also convert to optical properties of rrtmg interface, even though
   !   these quantities are later multiplied back together inside rrtmg !
   ! Why does rrtmg use the factored quantities?
   ! There are several different ways this factoring could be done.
   ! Other ways might allow for better optimization
   do ns = 1, nbndsw
      do k  = 1, rrtmg_levs-1
         kk=(pverp-rrtmg_levs) + k
         do i  = 1, Nday
            if(E_aer_tau_w(IdxDay(i),kk,ns) > 1.e-80_r8) then
               asm_aer_sw(i,k,ns) = E_aer_tau_w_g(IdxDay(i),kk,ns)/E_aer_tau_w(IdxDay(i),kk,ns)
            else
               asm_aer_sw(i,k,ns) = 0._r8
            endif
            if(E_aer_tau(IdxDay(i),kk,ns) > 0._r8) then
               ssa_aer_sw(i,k,ns) = E_aer_tau_w(IdxDay(i),kk,ns)/E_aer_tau(IdxDay(i),kk,ns)
               tau_aer_sw(i,k,ns) = E_aer_tau(IdxDay(i),kk,ns)
            else
               ssa_aer_sw(i,k,ns) = 1._r8
               tau_aer_sw(i,k,ns) = 0._r8
            endif
         enddo
      enddo
   enddo

   if (scm_crm_mode) then
      ! overwrite albedos for CRM
      if(have_asdir) asdir = asdirobs(1)
      if(have_asdif) asdif = asdifobs(1)
      if(have_aldir) aldir = aldirobs(1)
      if(have_aldif) aldif = aldifobs(1)
   endif

   ! Define solar incident radiation
   do i = 1, Nday
      solin(i)  = sum(sfac(:)*solar_band_irrad(:)) * eccf * coszrs(i)
   end do

   ! Calculate cloud optical properties here if using CAM method, or if using one of the
   ! methods in RRTMG_SW, then pass in cloud physical properties and zero out cloud optical 
   ! properties here

   ! Zero optional cloud optical property input arrays tauc_sw, ssac_sw, asmc_sw, 
   ! if inputting cloud physical properties to RRTMG_SW
   !tauc_sw(:,:,:) = 0.0_r8
   !ssac_sw(:,:,:) = 1.0_r8
   !asmc_sw(:,:,:) = 0.0_r8
   !fsfc_sw(:,:,:) = 0.0_r8
   !
   ! Or, calculate and pass in CAM cloud shortwave optical properties to RRTMG_SW
   !if (present(old_convert)) print *, 'old_convert',old_convert
   !if (present(ancientmethod)) print *, 'ancientmethod',ancientmethod
   if (present(old_convert))then
      if (old_convert)then ! convert without limits
         do i = 1, Nday
         do k = 1, rrtmg_levs-1
         kk=(pverp-rrtmg_levs) + k
         do ns = 1, nbndsw
           if (E_cld_tau_w(ns,IdxDay(i),kk) > 0._r8) then
              fsfc_sw(ns,i,k)=E_cld_tau_w_f(ns,IdxDay(i),kk)/E_cld_tau_w(ns,IdxDay(i),kk)
              asmc_sw(ns,i,k)=E_cld_tau_w_g(ns,IdxDay(i),kk)/E_cld_tau_w(ns,IdxDay(i),kk)
           else
              fsfc_sw(ns,i,k) = 0._r8
              asmc_sw(ns,i,k) = 0._r8
           endif
   
           tauc_sw(ns,i,k)=E_cld_tau(ns,IdxDay(i),kk)
           if (tauc_sw(ns,i,k) > 0._r8) then
              ssac_sw(ns,i,k)=E_cld_tau_w(ns,IdxDay(i),kk)/tauc_sw(ns,i,k)
           else
              tauc_sw(ns,i,k) = 0._r8
              fsfc_sw(ns,i,k) = 0._r8
              asmc_sw(ns,i,k) = 0._r8
              ssac_sw(ns,i,k) = 1._r8
           endif
         enddo
         enddo
         enddo
      else
         ! eventually, when we are done with archaic versions, This set of code will become the default.
         do i = 1, Nday
         do k = 1, rrtmg_levs-1
         kk=(pverp-rrtmg_levs) + k
         do ns = 1, nbndsw
           if (E_cld_tau_w(ns,IdxDay(i),kk) > 0._r8) then
              fsfc_sw(ns,i,k)=E_cld_tau_w_f(ns,IdxDay(i),kk)/max(E_cld_tau_w(ns,IdxDay(i),kk), 1.e-80_r8)
              asmc_sw(ns,i,k)=E_cld_tau_w_g(ns,IdxDay(i),kk)/max(E_cld_tau_w(ns,IdxDay(i),kk), 1.e-80_r8)
           else
              fsfc_sw(ns,i,k) = 0._r8
              asmc_sw(ns,i,k) = 0._r8
           endif
   
           tauc_sw(ns,i,k)=E_cld_tau(ns,IdxDay(i),kk)
           if (tauc_sw(ns,i,k) > 0._r8) then
              ssac_sw(ns,i,k)=max(E_cld_tau_w(ns,IdxDay(i),kk),1.e-80_r8)/max(tauc_sw(ns,i,k),1.e-80_r8)
           else
              tauc_sw(ns,i,k) = 0._r8
              fsfc_sw(ns,i,k) = 0._r8
              asmc_sw(ns,i,k) = 0._r8
              ssac_sw(ns,i,k) = 1._r8
           endif
         enddo
         enddo
         enddo
      endif
   else
      do i = 1, Nday
      do k = 1, rrtmg_levs-1
      kk=(pverp-rrtmg_levs) + k
      do ns = 1, nbndsw
        if (E_cld_tau_w(ns,IdxDay(i),kk) > 0._r8) then
           fsfc_sw(ns,i,k)=E_cld_tau_w_f(ns,IdxDay(i),kk)/max(E_cld_tau_w(ns,IdxDay(i),kk), 1.e-80_r8)
           asmc_sw(ns,i,k)=E_cld_tau_w_g(ns,IdxDay(i),kk)/max(E_cld_tau_w(ns,IdxDay(i),kk), 1.e-80_r8)
        else
           fsfc_sw(ns,i,k) = 0._r8
           asmc_sw(ns,i,k) = 0._r8
        endif

        tauc_sw(ns,i,k)=E_cld_tau(ns,IdxDay(i),kk)
        if (tauc_sw(ns,i,k) > 0._r8) then
           ssac_sw(ns,i,k)=max(E_cld_tau_w(ns,IdxDay(i),kk),1.e-80_r8)/max(tauc_sw(ns,i,k),1.e-80_r8)
        else
           tauc_sw(ns,i,k) = 0._r8
           fsfc_sw(ns,i,k) = 0._r8
           asmc_sw(ns,i,k) = 0._r8
           ssac_sw(ns,i,k) = 1._r8
        endif
      enddo
      enddo
      enddo
   endif

   ! Call mcica sub-column generator for RRTMG_SW

   ! Call sub-column generator for McICA in radiation
   call t_startf('mcica_subcol_sw')

   ! Select cloud overlap approach (1=random, 2=maximum-random, 3=maximum)
   icld = 2
   ! Set permute seed (must be offset between LW and SW by at least 140 to insure 
   ! effective randomization)
   permuteseed = 1


   call mcica_subcol_sw(lchnk, Nday, rrtmg_levs-1, icld, permuteseed, pmid, &
      cld, cicewp, cliqwp, rei, rel, tauc_sw, ssac_sw, asmc_sw, fsfc_sw, &
      cld_stosw, cicewp_stosw, cliqwp_stosw, rei_stosw, rel_stosw, &
      tauc_stosw, ssac_stosw, asmc_stosw, fsfc_stosw, rngsw, pergro_mods) !BSINGH- added rngsw

   call t_stopf('mcica_subcol_sw')

   call t_startf('rrtmg_sw')

   ! Call RRTMG_SW for all layers for daylight columns

   ! Select parameterization of cloud ice and liquid optical depths
   ! Use CAM shortwave cloud optical properties directly
   inflgsw = 0 
   iceflgsw = 0
   liqflgsw = 0
   ! Use E&C param for ice to mimic CAM3 for now
   !   inflgsw = 2 
   !   iceflgsw = 1
   !   liqflgsw = 1
   ! Use merged Fu and E&C params for ice 
   !   inflgsw = 2 
   !   iceflgsw = 3
   !   liqflgsw = 1

   ! Set day of year for Earth/Sun distance calculation in rrtmg_sw, or
   ! set to zero and pass E/S adjustment (eccf) directly into array adjes
   dyofyr = 0

   tsfc(:ncol) = tlev(:ncol,rrtmg_levs+1)

   solvar(1:nbndsw) = sfac(1:nbndsw)

   call rrtmg_sw(lchnk, Nday, rrtmg_levs, icld,         &
                 pmidmb, pintmb, tlay, tlev, tsfc, &
                 h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, &
                 asdir, asdif, aldir, aldif, &
                 coszrs, eccf, dyofyr, solvar, &
                 inflgsw, iceflgsw, liqflgsw, &
                 cld_stosw, tauc_stosw, ssac_stosw, asmc_stosw, fsfc_stosw, &
                 cicewp_stosw, cliqwp_stosw, rei, rel, &
                 tau_aer_sw, ssa_aer_sw, asm_aer_sw, &
                 swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
                 dirdnuv, dirdnir, difdnuv, difdnir, ninflx, ninflxc, swuflxs, swdflxs)

   ! Flux units are in W/m2 on output from rrtmg_sw and contain output for
   ! extra layer above model top with vertical indexing from bottom to top.
   !
   ! Heating units are in J/kg/s on output from rrtmg_sw and contain output 
   ! for extra layer above model top with vertical indexing from bottom to top.  
   !
   ! Reverse vertical indexing to go from top to bottom for CAM output.

   ! Set the net absorted shortwave flux at TOA (top of extra layer)
   fsntoa(1:Nday) = swdflx(1:Nday,rrtmg_levs+1) - swuflx(1:Nday,rrtmg_levs+1)
   fsutoa(1:Nday) = swuflx(1:Nday,rrtmg_levs+1)
   fsntoac(1:Nday) = swdflxc(1:Nday,rrtmg_levs+1) - swuflxc(1:Nday,rrtmg_levs+1)

   ! Set net near-IR flux at top of the model
   fsnirtoa(1:Nday) = ninflx(1:Nday,rrtmg_levs)
   fsnrtoaq(1:Nday) = ninflx(1:Nday,rrtmg_levs)
   fsnrtoac(1:Nday) = ninflxc(1:Nday,rrtmg_levs)

   ! Set the net absorbed shortwave flux at the model top level
   fsnt(1:Nday) = swdflx(1:Nday,rrtmg_levs) - swuflx(1:Nday,rrtmg_levs)
   fsntc(1:Nday) = swdflxc(1:Nday,rrtmg_levs) - swuflxc(1:Nday,rrtmg_levs)

   ! Set the downwelling flux at the surface 
   fsds(1:Nday) = swdflx(1:Nday,1)
   fsdsc(1:Nday) = swdflxc(1:Nday,1)

   ! Set the net shortwave flux at the surface
   fsns(1:Nday) = swdflx(1:Nday,1) - swuflx(1:Nday,1)
   fsnsc(1:Nday) = swdflxc(1:Nday,1) - swuflxc(1:Nday,1)

   ! Set the UV/vis and near-IR direct and dirruse downward shortwave flux at surface
   sols(1:Nday) = dirdnuv(1:Nday,1)
   soll(1:Nday) = dirdnir(1:Nday,1)
   solsd(1:Nday) = difdnuv(1:Nday,1)
   solld(1:Nday) = difdnir(1:Nday,1)


   ! Set the net, up and down fluxes at model interfaces
   fns (1:Nday,pverp-rrtmg_levs+1:pverp) =  swdflx(1:Nday,rrtmg_levs:1:-1) -  swuflx(1:Nday,rrtmg_levs:1:-1)
   fcns(1:Nday,pverp-rrtmg_levs+1:pverp) = swdflxc(1:Nday,rrtmg_levs:1:-1) - swuflxc(1:Nday,rrtmg_levs:1:-1)
   fus (1:Nday,pverp-rrtmg_levs+1:pverp) =  swuflx(1:Nday,rrtmg_levs:1:-1)
   fusc(1:Nday,pverp-rrtmg_levs+1:pverp) = swuflxc(1:Nday,rrtmg_levs:1:-1)
   fds (1:Nday,pverp-rrtmg_levs+1:pverp) =  swdflx(1:Nday,rrtmg_levs:1:-1)
   fdsc(1:Nday,pverp-rrtmg_levs+1:pverp) = swdflxc(1:Nday,rrtmg_levs:1:-1)

   ! Set solar heating, reverse layering
   ! Pass shortwave heating to CAM arrays and convert from K/d to J/kg/s
   qrs (1:Nday,pverp-rrtmg_levs+1:pver) = swhr (1:Nday,rrtmg_levs-1:1:-1)*cpair*dps
   qrsc(1:Nday,pverp-rrtmg_levs+1:pver) = swhrc(1:Nday,rrtmg_levs-1:1:-1)*cpair*dps

   ! Set spectral fluxes, reverse layering
   ! order=(/3,1,2/) maps the first index of swuflxs to the third index of su.
   if (associated(su)) then
      su(1:Nday,pverp-rrtmg_levs+1:pverp,:) = reshape(swuflxs(:,1:Nday,rrtmg_levs:1:-1), &
           (/Nday,rrtmg_levs,nbndsw/), order=(/3,1,2/))
   end if

   if (associated(sd)) then
      sd(1:Nday,pverp-rrtmg_levs+1:pverp,:) = reshape(swdflxs(:,1:Nday,rrtmg_levs:1:-1), &
           (/Nday,rrtmg_levs,nbndsw/), order=(/3,1,2/))
   end if

   call t_stopf('rrtmg_sw')

   ! Rearrange output arrays.
   !
   ! intent(out)

   call ExpDayNite(solin,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(qrs,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call ExpDayNite(qrsc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call ExpDayNite(fns,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call ExpDayNite(fcns,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call ExpDayNite(fsns,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnt,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsntoa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsutoa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsds,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnsc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsdsc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsntc,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsntoac,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(sols,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(soll,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(solsd,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(solld,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnirtoa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnrtoac,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnrtoaq,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)

   if (associated(su)) then
      call ExpDayNite(su,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp, 1, nbndsw)
   end if

   if (associated(sd)) then
      call ExpDayNite(sd,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp, 1, nbndsw)
   end if

   !  these outfld calls don't work for spmd only outfield in scm mode (nonspmd)
   if (single_column .and. scm_crm_mode) then 
      ! Following outputs added for CRM
      call ExpDayNite(fus,Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
      call ExpDayNite(fds,Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
      call ExpDayNite(fusc,Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
      call ExpDayNite(fdsc,Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
      call outfld('FUS     ',fus * 1.e-3_r8 ,pcols,lchnk)
      call outfld('FDS     ',fds * 1.e-3_r8 ,pcols,lchnk)
      call outfld('FUSC    ',fusc,pcols,lchnk)
      call outfld('FDSC    ',fdsc,pcols,lchnk)
   endif

end subroutine rad_rrtmg_sw

!-------------------------------------------------------------------------------

subroutine radsw_init()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize various constants for radiation scheme.
!
!-----------------------------------------------------------------------
    use radconstants,  only: get_solar_band_fraction_irrad, get_ref_solar_band_irrad
    use phys_control,  only: phys_getopts

    call phys_getopts(pergro_mods_out=pergro_mods)
    ! get the reference fractional solar irradiance in each band
    call get_solar_band_fraction_irrad(fractional_solar_irradiance)
    call get_ref_solar_band_irrad( solar_band_irrad )


   ! Initialize rrtmg_sw
   call rrtmg_sw_ini
 
end subroutine radsw_init


!-------------------------------------------------------------------------------

end module radsw
