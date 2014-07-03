
module radsw
!----------------------------------------------------------------------- 
! 
! Purpose: Solar radiation calculations.
!
!-----------------------------------------------------------------------
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use abortutils,      only: endrun
use cam_history,     only: outfld
use scamMod,         only: single_column,scm_crm_mode,have_asdir, &
                           asdirobs, have_asdif, asdifobs, have_aldir, &
                           aldirobs, have_aldif, aldifobs
use cam_logfile,     only: iulog
use radconstants,    only: nswbands, get_sw_spectral_boundaries, &
                           idx_sw_diag, indxsl

implicit none

private
save

! Public methods

public ::&
   radsw_init,   &! initialize constants
   radcswmx       ! driver for solar radiation code

! Private module data

real(r8) :: gravit     ! Acceleration of gravity
real(r8) :: rga        ! 1./gravit
real(r8) :: sslp       ! Standard sea-level pressure

!===============================================================================
CONTAINS
!===============================================================================

subroutine radcswmx(lchnk   ,ncol    ,                         &
                    E_pint    ,E_pmid    ,E_h2ommr  ,E_o3mmr   , &
                    E_o2mmr   ,E_cld     ,E_cicewp  ,E_cliqwp  ,E_rel     , &
                    E_rei     ,eccf      ,E_coszrs  ,solin     , &
                    E_asdir   ,E_asdif   ,E_aldir   ,E_aldif   ,nmxrgn  , &
                    pmxrgn  ,qrs,qrsc,fsnt    ,fsntc  ,fsdtoa,  fsntoa,   &
                    fsutoa ,fsntoac, fsnirtoa,fsnrtoac,fsnrtoaq,fsns    , &
                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                    solsd   ,solld   , fns     ,fcns            , &
                    Nday    ,Nnite   ,IdxDay  ,IdxNite, E_co2mmr, &
                    E_aer_tau, E_aer_tau_w, E_aer_tau_w_g, E_aer_tau_w_f, tauxcl_out, tauxci_out)
!-----------------------------------------------------------------------
! 

! Purpose: 
! Solar radiation code
! 
! Method: 
! Basic method is Delta-Eddington as described in:
! 
! Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
! 
! Five changes to the basic method described above are:
! (1) addition of sulfate aerosols (Kiehl and Briegleb, 1993)
! (2) the distinction between liquid and ice particle clouds 
! (Kiehl et al, 1996);
! (3) provision for calculating TOA fluxes with spectral response to
! match Nimbus-7 visible/near-IR radiometers (Collins, 1998);
! (4) max-random overlap (Collins, 2001)
! (5) The near-IR absorption by H2O was updated in 2003 by Collins, 
!     Lee-Taylor, and Edwards for consistency with the new line data in
!     Hitran 2000 and the H2O continuum version CKD 2.4.  Modifications
!     were optimized by reducing RMS errors in heating rates relative
!     to a series of benchmark calculations for the 5 standard AFGL 
!     atmospheres.  The benchmarks were performed using DISORT2 combined
!     with GENLN3.  The near-IR scattering optical depths for Rayleigh
!     scattering were also adjusted, as well as the correction for
!     stratospheric heating by H2O.
!
! The treatment of maximum-random overlap is described in the
! comment block "INDEX CALCULATIONS FOR MAX OVERLAP".
! 
! Divides solar spectrum into 19 intervals from 0.2-5.0 micro-meters.
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
! cgs units are used.
! 
! Special diagnostic calculation of the clear sky surface and total column
! absorbed flux is also done for cloud forcing diagnostics.
! 
!-----------------------------------------------------------------------
!
! D. Parks (NEC) 09/11/03
! Restructuring of routine to support SX vector architecture.
!
! Possible improvements:
!
! 1. Look at vectorizing index calculations for maximum overlap.
!
! 2. Consider making innermost loop in flux computations the number
!    of spectral intervals.  Given that NS is fixed at 19, the trade-off
!    will be stride one memory accesses of length 19 versus indirect
!    addressing (list vector - gather/scatter) with potential vector
!    lenghts of the number of day light points.  Vectorizing on the number
!    of spectral intervals seems worthwhile for low resolution models (T42),
!    but might be inefficient with higher resolutions.
!
! 3. Move the linearization of daylight points (compression/expansion) out
!    of radcswmx and into d_p_coupling.  This would eliminate the cost of
!    routines CmpDayNite and ExpDayNite.
!
! 4. Look at expliciting computing all streams in upward propagation of
!    radiation. There would be additional floating point operations in
!    exchange for the elimination of indirect addressing.
!
!-----------------------------------------------------------------------

   use rad_solar_var,    only: get_variability
   use cmparray_mod,     only: CmpDayNite, ExpDayNite
   use quicksort,        only: quick_sort
   use phys_control,     only: phys_getopts
   use solar_data,       only: sol_tsi, do_spctrl_scaling, ref_tsi
   use radconstants,     only: frcsol, ph2o, pco2, po2

!-----------------------Constants for new band (640-700 nm)-------------
   real(r8) v_raytau_35
   real(r8) v_raytau_64
   real(r8) v_abo3_35
   real(r8) v_abo3_64
   parameter( &
        v_raytau_35 = 0.155208_r8, &
        v_raytau_64 = 0.0392_r8, &
        v_abo3_35 = 2.4058030e+01_r8, &  
        v_abo3_64 = 2.210e+01_r8 &
        )


!-------------Parameters for accelerating max-random solution-------------
! 
! The solution time scales like prod(j:1->N) (1 + n_j) where 
! N   = number of max-overlap regions (nmxrgn)
! n_j = number of unique cloud amounts in region j
! 
! Therefore the solution cost can be reduced by decreasing n_j.
! cldmin reduces n_j by treating cloud amounts < cldmin as clear sky.
! cldeps reduces n_j by treating cloud amounts identical to log(1/cldeps)
! decimal places as identical
! 
! areamin reduces the cost by dropping configurations that occupy
! a surface area < areamin of the model grid box.  The surface area
! for a configuration C(j,k_j), where j is the region number and k_j is the
! index for a unique cloud amount (in descending order from biggest to
! smallest clouds) in region j, is
! 
! A = prod(j:1->N) [C(j,k_j) - C(j,k_j+1)]
! 
! where C(j,0) = 1.0 and C(j,n_j+1) = 0.0.
! 
! nconfgmax reduces the cost and improves load balancing by setting an upper
! bound on the number of cloud configurations in the solution.  If the number
! of configurations exceeds nconfgmax, the nconfgmax configurations with the
! largest area are retained, and the fluxes are normalized by the total area
! of these nconfgmax configurations.  For the current max/random overlap 
! assumption (see subroutine cldovrlap), 30 levels, and cloud-amount 
! parameterization, the mean and RMS number of configurations are 
! both roughly 5.  nconfgmax has been set to the mean+2*RMS number, or 15.
! 
! Minimum cloud amount (as a fraction of the grid-box area) to 
! distinguish from clear sky
! 
   real(r8) cldmin
   parameter (cldmin = 1.0e-80_r8)
! 
! Minimimum horizontal area (as a fraction of the grid-box area) to retain 
! for a unique cloud configuration in the max-random solution
! 
   real(r8) areamin
   parameter (areamin = 0.01_r8)
! 
! Decimal precision of cloud amount (0 -> preserve full resolution;
! 10^-n -> preserve n digits of cloud amount)
! 
   real(r8) cldeps
   parameter (cldeps = 0.0_r8)
! 
! Maximum number of configurations to include in solution
! 
   integer nconfgmax
   parameter (nconfgmax = 15)
!------------------------------Commons----------------------------------
! 
! Input arguments
! 
   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: ncol              ! number of atmospheric columns

   integer,intent(in) :: Nday                      ! Number of daylight columns
   integer,intent(in) :: Nnite                     ! Number of night columns
   integer,intent(in), dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
   integer,intent(in), dimension(pcols) :: IdxNite ! Indicies of night coumns


   real(r8), intent(in) :: E_pmid(pcols,pver) ! Level pressure
   real(r8), intent(in) :: E_pint(pcols,pverp) ! Interface pressure
   real(r8), intent(in) :: E_h2ommr(pcols,pver) ! Specific humidity (h2o mass mix ratio)
   real(r8), intent(in) :: E_o3mmr(pcols,pver) ! Ozone mass mixing ratio
   real(r8), intent(in) :: E_o2mmr(pcols) ! oxygen mass mixing ratio
! 
   real(r8), intent(in) :: E_cld(pcols,pver)  ! Fractional cloud cover
   real(r8), intent(in) :: E_cicewp(pcols,pver) ! in-cloud cloud ice water path
   real(r8), intent(in) :: E_cliqwp(pcols,pver) ! in-cloud cloud liquid water path
   real(r8), intent(in) :: E_rel(pcols,pver)  ! Liquid effective drop size (microns)
   real(r8), intent(in) :: E_rei(pcols,pver)  ! Ice effective drop size (microns)
! 
   real(r8), intent(in) :: eccf             ! Eccentricity factor (1./earth-sun dist^2)
   real(r8), intent(in) :: E_coszrs(pcols)    ! Cosine solar zenith angle
   real(r8), intent(in) :: E_asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: E_aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: E_asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
   real(r8), intent(in) :: E_aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad
   real(r8), intent(in) :: E_co2mmr(pcols)    ! co2 column mean mmr

! 
! Aerosol radiative property arrays
! 
   real(r8),intent(in) :: E_aer_tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
   real(r8),intent(in) :: E_aer_tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
   real(r8),intent(in) :: E_aer_tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * w * tau
   real(r8),intent(in) :: E_aer_tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * w * tau

! 
! IN/OUT arguments
! 
   real(r8), intent(inout) :: pmxrgn(pcols,pverp) ! Maximum values of pressure for each
!                                                 !    maximally overlapped region. 
!                                                 !    0->pmxrgn(i,1) is range of pressure for
!                                                 !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                                 !    2nd region, etc
   integer, intent(inout) ::  nmxrgn(pcols)    ! Number of maximally overlapped regions
! 
! Output arguments
! 

   real(r8), intent(out) :: solin(pcols)     ! Incident solar flux
   real(r8), intent(out) :: qrs (pcols,pver) ! Solar heating rate
   real(r8), intent(out) :: qrsc(pcols,pver) ! Clearsky solar heating rate
   real(r8), intent(out) :: fsns(pcols)      ! Surface absorbed solar flux
   real(r8), intent(out) :: fsnt(pcols)      ! Total column absorbed solar flux
   real(r8), intent(out) :: fsntoa(pcols)    ! Net solar flux at TOA
   real(r8), intent(out) :: fsutoa(pcols)    ! Upward solar flux at TOA
   real(r8), intent(out) :: fsds(pcols)      ! Flux shortwave downwelling surface
! 
   real(r8), intent(out) :: fsnsc(pcols)     ! Clear sky surface absorbed solar flux
   real(r8), intent(out) :: fsdsc(pcols)     ! Clear sky surface downwelling solar flux
   real(r8), intent(out) :: fsntc(pcols)     ! Clear sky total column absorbed solar flx
   real(r8), intent(out) :: fsdtoa(pcols)    ! Downwelling solar flux at TOA
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

   real(r8), intent(out) :: tauxcl_out(pcols,pver) ! liquid cloud visible sw optical depth
   real(r8), intent(out) :: tauxci_out(pcols,pver) ! ice cloud visible sw optical depth

! 
!---------------------------Local variables-----------------------------
!
! Local and reordered copies of the intent(in) variables
!
   real(r8):: aer_tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
   real(r8):: aer_tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
   real(r8):: aer_tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * w * tau
   real(r8):: aer_tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * w * tau
   real(r8) :: pmid(pcols,pver) ! Level pressure
   real(r8) :: pint(pcols,pverp) ! Interface pressure
   real(r8) :: h2ommr(pcols,pver) ! Specific humidity (h2o mass mix ratio)
   real(r8) :: o3mmr(pcols,pver) ! Ozone mass mixing ratio
! 
   real(r8) :: cld(pcols,pver)  ! Fractional cloud cover
   real(r8) :: cicewp(pcols,pver) ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,pver) ! in-cloud cloud liquid water path
   real(r8) :: rel(pcols,pver)  ! Liquid effective drop size (microns)
   real(r8) :: rei(pcols,pver)  ! Ice effective drop size (microns)
! 
   real(r8) :: coszrs(pcols)    ! Cosine solar zenith angle
   real(r8) :: asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
   real(r8) :: aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
   real(r8) :: asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
   real(r8) :: aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad
   real(r8) :: co2mmr(pcols)    ! co2 column mean mmr
   real(r8) :: o2mmr(pcols)     ! o2 column mean mmr

   real(r8) :: tot_irrad
! 
! Max/random overlap variables
! 
   real(r8) asort(pverp)     ! 1 - cloud amounts to be sorted for max ovrlp.
   real(r8) atmp             ! Temporary storage for sort when nxs = 2
   real(r8) cld0             ! 1 - (cld amt) used to make wstr, cstr, nstr
   real(r8) totwgt(pcols)    ! Total of xwgts = total fractional area of 
!   grid-box covered by cloud configurations
!   included in solution to fluxes

   real(r8) wgtv(nconfgmax)  ! Weights for fluxes
!   1st index is configuration number
   real(r8) wstr(pverp,pverp) ! area weighting factors for streams
!   1st index is for stream #, 
!   2nd index is for region #

   real(r8) xexpt            ! solar direct beam trans. for layer above
   real(r8) xrdnd            ! diffuse reflectivity for layer above
   real(r8) xrupd            ! diffuse reflectivity for layer below
   real(r8) xrups            ! direct-beam reflectivity for layer below
   real(r8) xtdnt            ! total trans for layers above

   real(r8) xwgt             ! product of cloud amounts

   real(r8) yexpt            ! solar direct beam trans. for layer above
   real(r8) yrdnd            ! diffuse reflectivity for layer above
   real(r8) yrupd            ! diffuse reflectivity for layer below
   real(r8) ytdnd            ! dif-beam transmission for layers above
   real(r8) ytupd            ! dif-beam transmission for layers below

   real(r8) zexpt            ! solar direct beam trans. for layer above
   real(r8) zrdnd            ! diffuse reflectivity for layer above
   real(r8) zrupd            ! diffuse reflectivity for layer below
   real(r8) zrups            ! direct-beam reflectivity for layer below
   real(r8) ztdnt            ! total trans for layers above

   logical new_term          ! Flag for configurations to include in fluxes
   logical region_found      ! flag for identifying regions

   integer ccon(nconfgmax,0:pverp,pcols)                                
! flags for presence of clouds
!   1st index is for level # (including 
!    layer above top of model and at surface)
!   2nd index is for configuration #
   integer cstr(0:pverp,pverp)                                
! flags for presence of clouds
!   1st index is for level # (including 
!    layer above top of model and at surface)
!   2nd index is for stream #
   integer icond(nconfgmax,0:pverp,pcols)
! Indices for copying rad. properties from
!     one identical downward cld config.
!     to another in adding method (step 2)
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer iconu(nconfgmax,0:pverp,pcols)
! Indices for copying rad. properties from
!     one identical upward configuration
!     to another in adding method (step 2)
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer iconfig           ! Counter for random-ovrlap configurations
   integer irgn              ! Index for max-overlap regions
   integer is0               ! Lower end of stream index range
   integer is1               ! Upper end of stream index range
   integer isn               ! Stream index
   integer istr(pverp+1)     ! index for stream #s during flux calculation
   integer istrtd(0:nconfgmax+1,0:pverp,pcols)
! indices into icond 
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer istrtu(0:nconfgmax+1,0:pverp,pcols)
! indices into iconu 
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer j                 ! Configuration index
   integer jj                ! Configuration index
   integer k1                ! Level index
   integer k2                ! Level index
   integer ksort(pverp)      ! Level indices of cloud amounts to be sorted
   integer ktmp              ! Temporary storage for sort when nxs = 2
   integer kx1(0:pverp)      ! Level index for top of max-overlap region
   integer kx2(0:pverp)      ! Level index for bottom of max-overlap region
   integer l                 ! Index 
   integer l0                ! Index
   integer mrgn              ! Counter for nrgn
   integer mstr              ! Counter for nstr
   integer n0                ! Number of configurations with ccon(:,k,:)==0
   integer n1                ! Number of configurations with ccon(:,k,:)==1
   integer nconfig(pcols)    ! Number of random-ovrlap configurations
   integer nconfigm          ! Value of config before testing for areamin,
!    nconfgmax
   integer npasses           ! number of passes over the indexing loop
   integer nrgn              ! Number of max overlap regions at current 
!    longitude
   integer nstr(pverp)       ! Number of unique cloud configurations
!   ("streams") in a max-overlapped region
!   1st index is for region #
   integer nuniq             ! # of unique cloud configurations
   integer nuniqd(0:pverp,pcols)   ! # of unique cloud configurations: TOA 
!   to level k
   integer nuniqu(0:pverp,pcols)   ! # of unique cloud configurations: surface
!   to level k 
   integer nxs               ! Number of cloudy layers between k1 and k2 
   integer ptr0(nconfgmax)   ! Indices of configurations with ccon(:,k,:)==0
   integer ptr1(nconfgmax)   ! Indices of configurations with ccon(:,k,:)==1
   integer ptrc(nconfgmax)   ! Pointer for configurations sorted by wgtv
   integer, dimension(1) :: min_idx  ! required for return val of func minloc

! 
! Other
! 
   integer ns                ! Spectral loop index
   integer i                 ! Longitude loop index
   integer k                 ! Level loop index
   integer km1               ! k - 1
   integer kp1               ! k + 1
   integer n                 ! Loop index for daylight
   integer ksz               ! dust size bin index
   integer kaer              ! aerosol group index
! 
! A. Slingo's data for cloud particle radiative properties (from 'A GCM
! Parameterization for the Shortwave Properties of Water Clouds' JAS
! vol. 46 may 1989 pp 1419-1427)
! 

   real(r8), parameter :: abarl(4) = (/ 2.817e-02_r8, 2.682e-02_r8,2.264e-02_r8,1.281e-02_r8/)
   real(r8), parameter :: bbarl(4) = (/ 1.305_r8    , 1.346_r8    ,1.454_r8    ,1.641_r8    /)
   real(r8), parameter :: cbarl(4) = (/-5.62e-08_r8 ,-6.94e-06_r8 ,4.64e-04_r8 ,0.201_r8    /)
   real(r8), parameter :: dbarl(4) = (/ 1.63e-07_r8 , 2.35e-05_r8 ,1.24e-03_r8 ,7.56e-03_r8 /)
   real(r8), parameter :: ebarl(4) = (/ 0.829_r8    , 0.794_r8    ,0.754_r8    ,0.826_r8    /)
   real(r8), parameter :: fbarl(4) = (/ 2.482e-03_r8, 4.226e-03_r8,6.560e-03_r8,4.353e-03_r8/)

   real(r8) :: abarli           ! A coefficient for current spectral band
   real(r8) :: bbarli           ! B coefficient for current spectral band
   real(r8) :: cbarli           ! C coefficient for current spectral band
   real(r8) :: dbarli           ! D coefficient for current spectral band
   real(r8) :: ebarli           ! E coefficient for current spectral band
   real(r8) :: fbarli           ! F coefficient for current spectral band
! 
! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
! greater than 20 micro-meters
! 
! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
! 
   real(r8), parameter :: abari(4) = (/ 3.448e-03_r8, 3.448e-03_r8,3.448e-03_r8,3.448e-03_r8/)
   real(r8), parameter :: bbari(4) = (/ 2.431_r8    , 2.431_r8    ,2.431_r8    ,2.431_r8    /)
   real(r8), parameter :: cbari(4) = (/ 1.00e-05_r8 , 1.10e-04_r8 ,1.861e-02_r8,.46658_r8   /)
   real(r8), parameter :: dbari(4) = (/ 0.0_r8      , 1.405e-05_r8,8.328e-04_r8,2.05e-05_r8 /)
   real(r8), parameter :: ebari(4) = (/ 0.7661_r8   , 0.7730_r8   ,0.794_r8    ,0.9595_r8   /)
   real(r8), parameter :: fbari(4) = (/ 5.851e-04_r8, 5.665e-04_r8,7.267e-04_r8,1.076e-04_r8/)

   real(r8) :: abarii           ! A coefficient for current spectral band
   real(r8) :: bbarii           ! B coefficient for current spectral band
   real(r8) :: cbarii           ! C coefficient for current spectral band
   real(r8) :: dbarii           ! D coefficient for current spectral band
   real(r8) :: ebarii           ! E coefficient for current spectral band
   real(r8) :: fbarii           ! F coefficient for current spectral band

!
! UPDATE TO H2O NEAR-IR: Delta optimized for Hitran 2K and CKD 2.4
!
   real(r8), parameter :: delta = 0.0014257179260883_r8
!
! END UPDATE
!
   real(r8) :: albdir(pcols,nswbands) ! Current spc intrvl srf alb to direct rad
   real(r8) :: albdif(pcols,nswbands) ! Current spc intrvl srf alb to diffuse rad
! 
! Next series depends on spectral interval
! 
   real(r8) :: wgtint           ! Weight for specific spectral interval

! 
! weight for 0.64 - 0.7 microns  appropriate to clear skies over oceans
! 
   real(r8), parameter :: nirwgt(nswbands) = &
              (/  0.0_r8,   0.0_r8,   0.0_r8,      0.0_r8,   0.0_r8, &
                  0.0_r8,   0.0_r8,   0.0_r8, 0.320518_r8,   1.0_r8,  1.0_r8, &
                  1.0_r8,   1.0_r8,   1.0_r8,      1.0_r8,   1.0_r8, &
                  1.0_r8,   1.0_r8,   1.0_r8 /)

!
! UPDATE TO H2O NEAR-IR: Rayleigh scattering optimized for Hitran 2K & CKD 2.4
!
   real(r8), parameter :: raytau(nswbands) = &
               (/ 4.020_r8, 2.180_r8, 1.700_r8, 1.450_r8, 1.250_r8, &
                  1.085_r8, 0.730_r8, v_raytau_35, v_raytau_64, &
                  0.02899756_r8, 0.01356763_r8, 0.00537341_r8, &
                  0.00228515_r8, 0.00105028_r8, 0.00046631_r8, &
                  0.00025734_r8, &
                 .0001_r8, .0001_r8, .0001_r8/)
!
! END UPDATE
!

! 
! Absorption coefficients
! 
!
! UPDATE TO H2O NEAR-IR: abh2o optimized for Hitran 2K and CKD 2.4
!
   real(r8), parameter :: abh2o(nswbands) = &
             (/    .000_r8,     .000_r8,    .000_r8,    .000_r8,    .000_r8, &
                   .000_r8,     .000_r8,    .000_r8,    .000_r8,    &
                   0.00256608_r8,  0.06310504_r8,   0.42287445_r8, 2.45397941_r8, &
                  11.20070807_r8, 47.66091389_r8, 240.19010243_r8, &
                   .000_r8,    .000_r8,    .000_r8/)
!
! END UPDATE
!

   real(r8), parameter :: abo3(nswbands) = &
             (/5.370e+04_r8, 13.080e+04_r8,  9.292e+04_r8, 4.530e+04_r8, 1.616e+04_r8, &
               4.441e+03_r8,  1.775e+02_r8, v_abo3_35, v_abo3_64,      .000_r8, &
               .000_r8,   .000_r8    ,   .000_r8   ,   .000_r8   ,      .000_r8, &
               .000_r8,   .000_r8    ,   .000_r8   ,   .000_r8    /)

   real(r8), parameter :: abco2(nswbands) = &
              (/   .000_r8,     .000_r8,    .000_r8,    .000_r8,    .000_r8, &
                   .000_r8,     .000_r8,    .000_r8,    .000_r8,    .000_r8, &
                   .000_r8,     .000_r8,    .000_r8,    .000_r8,    .000_r8, &
                   .000_r8,     .094_r8,    .196_r8,   1.963_r8/)

   real(r8), parameter :: abo2(nswbands) = &
             (/    .000_r8,     .000_r8,    .000_r8,    .000_r8,    .000_r8, &
                   .000_r8,     .000_r8,    .000_r8,1.11e-05_r8,6.69e-05_r8, &
                   .000_r8,     .000_r8,    .000_r8,    .000_r8,    .000_r8, &  
                   .000_r8,     .000_r8,    .000_r8,    .000_r8/)
! 
! Diagnostic and accumulation arrays; note that fswup, and
! fswdn are not used in the computation,but are retained for future use.
! 
   real(r8) solflx(pcols)    ! Solar flux in current interval
   real(r8) totfld (pcols,0:pver)  ! Spectrally summed flux divergence
   real(r8) totfldc(pcols,0:pver)  ! Spectrally summed flux divergence (clearsky)
   real(r8) fswup(pcols,0:pverp)   ! Spectrally summed up flux
   real(r8) fswdn(pcols,0:pverp)   ! Spectrally summed down flux
! 
! Cloud radiative property arrays
! 
   real(r8) tauxcl(pcols,0:pver) ! water cloud extinction optical depth
   real(r8) tauxci(pcols,0:pver) ! ice cloud extinction optical depth
   real(r8) wcl(pcols,0:pver) ! liquid cloud single scattering albedo
   real(r8) gcl(pcols,0:pver) ! liquid cloud asymmetry parameter
   real(r8) fcl(pcols,0:pver) ! liquid cloud forward scattered fraction
   real(r8) wci(pcols,0:pver) ! ice cloud single scattering albedo
   real(r8) gci(pcols,0:pver) ! ice cloud asymmetry parameter
   real(r8) fci(pcols,0:pver) ! ice cloud forward scattered fraction

! 
! Various arrays and other constants:
! 
   real(r8) pflx(pcols,0:pverp) ! Interface press, including extra layer
   real(r8) zenfac(pcols)    ! Square root of cos solar zenith angle
   real(r8) sqrco2(pcols)    ! Square root of the co2 mass mixg ratio
   real(r8) tmp1             ! Temporary constant array
   real(r8) tmp2             ! Temporary constant array
   real(r8) pdel             ! Pressure difference across layer
   real(r8) path             ! Mass path of layer
   real(r8) ptop             ! Lower interface pressure of extra layer
   real(r8) ptho2            ! Used to compute mass path of o2
   real(r8) ptho3            ! Used to compute mass path of o3
   real(r8) pthco2           ! Used to compute mass path of co2
   real(r8) pthh2o           ! Used to compute mass path of h2o
   real(r8) h2ostr           ! Inverse sq. root h2o mass mixing ratio

   real(r8) wavmin(nswbands) ! Spectral interval minimum wavelength
   real(r8) wavmax(nswbands) ! Spectral interval maximum wavelength
   real(r8) wavmid(nswbands) ! Spectral interval middle wavelength
   real(r8) trayoslp         ! Rayleigh optical depth/standard pressure
   real(r8) tmp1l            ! Temporary constant array
   real(r8) tmp2l            ! Temporary constant array
   real(r8) tmp3l            ! Temporary constant array
   real(r8) tmp1i            ! Temporary constant array
   real(r8) tmp2i            ! Temporary constant array
   real(r8) tmp3i            ! Temporary constant array
   real(r8) rdenom           ! Multiple scattering term
   real(r8) rdirexp          ! layer direct ref times exp transmission
   real(r8) tdnmexp          ! total transmission - exp transmission
   real(r8) psf(nswbands)      ! Frac of solar flux in spect interval
! 
! Layer absorber amounts; note that 0 refers to the extra layer added
! above the top model layer
! 
   real(r8) uh2o(pcols,0:pver) ! Layer absorber amount of h2o
   real(r8) uo3(pcols,0:pver) ! Layer absorber amount of  o3
   real(r8) uco2(pcols,0:pver) ! Layer absorber amount of co2
   real(r8) uo2(pcols,0:pver) ! Layer absorber amount of  o2
! 
! Total column absorber amounts:
! 
   real(r8) uth2o(pcols)     ! Total column  absorber amount of  h2o
   real(r8) uto3(pcols)      ! Total column  absorber amount of  o3
   real(r8) utco2(pcols)     ! Total column  absorber amount of  co2
   real(r8) uto2(pcols)      ! Total column  absorber amount of  o2
! 
! These arrays are defined for pver model layers; 0 refers to the extra
! layer on top:
! 
   real(r8) rdir(nswbands,pcols,0:pver) ! Layer reflectivity to direct rad
   real(r8) rdif(nswbands,pcols,0:pver) ! Layer reflectivity to diffuse rad
   real(r8) tdir(nswbands,pcols,0:pver) ! Layer transmission to direct rad
   real(r8) tdif(nswbands,pcols,0:pver) ! Layer transmission to diffuse rad
   real(r8) explay(nswbands,pcols,0:pver) ! Solar beam exp trans. for layer

   real(r8) rdirc(nswbands,pcols,0:pver) ! Clear Layer reflec. to direct rad
   real(r8) rdifc(nswbands,pcols,0:pver) ! Clear Layer reflec. to diffuse rad
   real(r8) tdirc(nswbands,pcols,0:pver) ! Clear Layer trans. to direct rad
   real(r8) tdifc(nswbands,pcols,0:pver) ! Clear Layer trans. to diffuse rad
   real(r8) explayc(nswbands,pcols,0:pver) ! Solar beam exp trans. clear layer

   real(r8) fus(pcols,pverp)   ! Upward flux (added for CRM)
   real(r8) fds(pcols,pverp)   ! Downward flux (added for CRM)
   real(r8) fusc(pcols,pverp)  ! Upward clear-sky flux (added for CRM)
   real(r8) fdsc(pcols,pverp) ! Downward clear-sky flux (added for CRM)

   real(r8) flxdiv           ! Flux divergence for layer

!
! Temporary arrays for either clear or cloudy values.
!
   real(r8), dimension(nswbands) :: Trdir
   real(r8), dimension(nswbands) :: Trdif
   real(r8), dimension(nswbands) :: Ttdir
   real(r8), dimension(nswbands) :: Ttdif
   real(r8), dimension(nswbands) :: Texplay
!cdir vreg(Trdir)
!cdir vreg(Trdif)
!cdir vreg(Ttdir)
!cdir vreg(Ttdif)
!cdir vreg(Texplay)
! 
! 
! Radiative Properties:
! 
! There are 1 classes of properties:
! (1. All-sky bulk properties
! (2. Clear-sky properties
! 
! The first set of properties are generated during step 2 of the solution.
! 
! These arrays are defined at model interfaces; in 1st index (for level #),
! 0 is the top of the extra layer above the model top, and
! pverp is the earth surface.  2nd index is for cloud configuration
! defined over a whole column.
! 
   real(r8) exptdn(nswbands,0:pverp,nconfgmax,pcols) ! Sol. beam trans from layers above
   real(r8) rdndif(nswbands,0:pverp,nconfgmax,pcols) ! Ref to dif rad for layers above
   real(r8) rupdif(nswbands,0:pverp,nconfgmax,pcols) ! Ref to dif rad for layers below
   real(r8) rupdir(nswbands,0:pverp,nconfgmax,pcols) ! Ref to dir rad for layers below
   real(r8) tdntot(nswbands,0:pverp,nconfgmax,pcols) ! Total trans for layers above
! 
! Bulk properties used during the clear-sky calculation.
! 
   real(r8) exptdnc(pcols,0:pverp) ! clr: Sol. beam trans from layers above
   real(r8) rdndifc(pcols,0:pverp) ! clr: Ref to dif rad for layers above
   real(r8) rupdifc(pcols,0:pverp) ! clr: Ref to dif rad for layers below
   real(r8) rupdirc(pcols,0:pverp) ! clr: Ref to dir rad for layers below
   real(r8) tdntotc(pcols,0:pverp) ! clr: Total trans for layers above

   real(r8) fluxup(nswbands,0:pverp,pcols)  ! Up   flux at model interface
   real(r8) fluxdn(nswbands,0:pverp,pcols)  ! Down flux at model interface
   real(r8) wexptdn(nswbands,pcols)   ! Direct solar beam trans. to surface
!
! Scalars used in vectorization
!
  integer :: kk
!
! Arrays used in vectorization
!
   real(r8) v_wgtv(nconfgmax,pcols)  ! Weights for fluxes

   real(r8) :: rdiff, ro, rn
   rdiff(ro,rn) = abs((ro-rn)/merge(ro,1.0_r8,ro /= 0.0_r8))

!  solar variability factor
   real(r8) :: sfac(nswbands)

   character(len=16) :: microp_scheme  ! microphysics scheme

   !-----------------------------------------------------------------------
   ! START OF CALCULATION
   !-----------------------------------------------------------------------

   call phys_getopts(microp_scheme_out=microp_scheme)

! 
! Initialize output fields:
! 
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
   fsdtoa(1:ncol)   = 0.0_r8
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

   tauxcl_out(1:pcols,1:pver) = 0.0_r8
   tauxci_out(1:pcols,1:pver) = 0.0_r8

! 
! If night everywhere, return:
! 
   if ( Nday == 0 ) then
     return
   endif

!
! Rearrange input arrays
!

   call CmpDayNite(E_pmid, pmid,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_pint, pint,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call CmpDayNite(E_h2ommr, h2ommr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_o3mmr, o3mmr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_cld, cld,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_cicewp, cicewp,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_cliqwp, cliqwp,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_rel, rel, 		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_rei, rei,		Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call CmpDayNite(E_coszrs, coszrs,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_asdir, asdir,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_aldir, aldir,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_asdif, asdif,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_aldif, aldif,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_co2mmr, co2mmr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_o2mmr, o2mmr,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)

   call CmpDayNite(pmxrgn,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call CmpDayNite(nmxrgn,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call CmpDayNite(E_aer_tau,     aer_tau,     Nday, IdxDay, Nnite, IdxNite, 1,pcols, 0,pver, 1,nswbands)
   call CmpDayNite(E_aer_tau_w,   aer_tau_w,   Nday, IdxDay, Nnite, IdxNite, 1,pcols, 0,pver, 1,nswbands)
   call CmpDayNite(E_aer_tau_w_g, aer_tau_w_g, Nday, IdxDay, Nnite, IdxNite, 1,pcols, 0,pver, 1,nswbands)
   call CmpDayNite(E_aer_tau_w_f, aer_tau_w_f, Nday, IdxDay, Nnite, IdxNite, 1,pcols, 0,pver, 1,nswbands)

   if (scm_crm_mode) then
   ! overwrite albedos for CRM
      if(have_asdir) asdir = asdirobs(1)
      if(have_asdif) asdif = asdifobs(1)
      if(have_aldir) aldir = aldirobs(1)
      if(have_aldif) aldif = aldifobs(1)
   endif
! 
! Perform other initializations
! 
   tmp1   = 0.5_r8/(gravit*sslp)
   tmp2   = delta/gravit

   do i = 1, Nday
      sqrco2(i) = sqrt(co2mmr(i))
   end do

   if ( do_spctrl_scaling ) then
      call get_variability(sfac)
      tot_irrad = ref_tsi
   else
      tot_irrad = sol_tsi
   endif

   do k=1,pverp
      do i=1,Nday
         pflx(i,k) = pint(i,k)
      end do
   end do
   
   do i=1,Nday
! 
! Define solar incident radiation and interface pressures:
! 

         solin(i)  = tot_irrad*1.e3_r8*eccf*coszrs(i)

         pflx(i,0) = 0._r8
! 
! Compute optical paths:
! 
         ptop      = pflx(i,1)
         ptho2     = o2mmr(i) * ptop / gravit
         ptho3     = o3mmr(i,1) * ptop / gravit
         pthco2    = sqrco2(i) * (ptop / gravit)
         h2ostr    = sqrt( 1._r8 / h2ommr(i,1) )
         zenfac(i) = sqrt(coszrs(i))
         pthh2o    = ptop**2*tmp1 + (ptop*rga)* &
                    (h2ostr*zenfac(i)*delta)
         uh2o(i,0) = h2ommr(i,1)*pthh2o
         uco2(i,0) = zenfac(i)*pthco2
         uo2 (i,0) = zenfac(i)*ptho2
         uo3 (i,0) = ptho3
! 
! End  do i=1,Nday
! 
   end do

   do k=1,pver

!cdir nodep
      do i=1,Nday

         pdel      = pflx(i,k+1) - pflx(i,k)
         path      = pdel / gravit
         ptho2     = o2mmr(i) * path
         ptho3     = o3mmr(i,k) * path
         pthco2    = sqrco2(i) * path
         h2ostr    = sqrt(1.0_r8/h2ommr(i,k))
         pthh2o    = (pflx(i,k+1)**2 - pflx(i,k)**2)*tmp1 + pdel*h2ostr*zenfac(i)*tmp2
         uh2o(i,k) = h2ommr(i,k)*pthh2o
         uco2(i,k) = zenfac(i)*pthco2
         uo2 (i,k) = zenfac(i)*ptho2
         uo3 (i,k) = ptho3

! 
! End  do i=1,Nday
! 
      end do
! 
! End  k=1,pver
! 
   end do
! 
! Compute column absorber amounts for the clear sky computation:
! 
   do i=1,Nday

      uth2o(i) = 0.0_r8
      uto3(i)  = 0.0_r8
      utco2(i) = 0.0_r8
      uto2(i)  = 0.0_r8

!cdir expand=pver
      do k=1,pver
         uth2o(i) = uth2o(i) + uh2o(i,k)
         uto3(i)  = uto3(i)  + uo3(i,k)
         utco2(i) = utco2(i) + uco2(i,k)
         uto2(i)  = uto2(i)  + uo2(i,k)
! 
! End  k=1,pver
! 
      end do
! 
! End  do i=1,Nday
! 
   end do
! 
! Set cloud properties for top (0) layer; so long as tauxcl is zero,
! there is no cloud above top of model; the other cloud properties
! are arbitrary:
! 
      do i=1,Nday

         tauxcl(i,0)  = 0._r8
         wcl(i,0)     = 0.999999_r8
         gcl(i,0)     = 0.85_r8
         fcl(i,0)     = 0.725_r8
         tauxci(i,0)  = 0._r8
         wci(i,0)     = 0.999999_r8
         gci(i,0)     = 0.85_r8
         fci(i,0)     = 0.725_r8
! 
! End  do i=1,Nday
! 
      end do
! 
! Begin spectral loop
! 
   do ns=1,nswbands
! 
! Set cloud extinction optical depth, single scatter albedo,
! asymmetry parameter, and forward scattered fraction:
! 
      abarli = abarl(indxsl(ns))
      bbarli = bbarl(indxsl(ns))
      cbarli = cbarl(indxsl(ns))
      dbarli = dbarl(indxsl(ns))
      ebarli = ebarl(indxsl(ns))
      fbarli = fbarl(indxsl(ns))
! 
      abarii = abari(indxsl(ns))
      bbarii = bbari(indxsl(ns))
      cbarii = cbari(indxsl(ns))
      dbarii = dbari(indxsl(ns))
      ebarii = ebari(indxsl(ns))
      fbarii = fbari(indxsl(ns))
! 
! adjustfraction within spectral interval to allow for the possibility of
! sub-divisions within a particular interval:
! 
      psf(ns) = 1.0_r8
      if(ph2o(ns)/=0._r8) psf(ns) = psf(ns)*ph2o(ns)
      if(pco2(ns)/=0._r8) psf(ns) = psf(ns)*pco2(ns)
      if(po2 (ns)/=0._r8) psf(ns) = psf(ns)*po2 (ns)

      do k=1,pver

         do i=1,Nday

            ! liquid
            ! note that optical properties for liquid valid only
            ! in range of 4.2 > rel > 16 micron (Slingo 89)
            if ( microp_scheme == 'MG' ) then
               tmp2l = 1._r8 - cbarli - dbarli*min(max(4.2_r8,rel(i,k)),16._r8)
               tmp3l = fbarli*min(max(4.2_r8,rel(i,k)),16._r8)
            else
               tmp2l = 1._r8 - cbarli - dbarli*rel(i,k)
               tmp3l = fbarli*rel(i,k)
            endif

            ! ice
            ! note that optical properties for ice valid only
            ! in range of 13 > rei > 130 micron (Ebert and Curry 92)
            if ( microp_scheme == 'MG' ) then
               tmp2i = 1._r8 - cbarii - dbarii*min(max(13._r8,rei(i,k)),130._r8)
               tmp3i = fbarii*min(max(13._r8,rei(i,k)),130._r8)
            else
               tmp2i = 1._r8 - cbarii - dbarii*rei(i,k)
               tmp3i = fbarii*rei(i,k)
            endif

            if (cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then

               ! liquid
               if ( microp_scheme == 'MG' ) then
                  tmp1l = abarli + bbarli/min(max(4.2_r8,rel(i,k)),16._r8)
               else
                  tmp1l = abarli + bbarli/rel(i,k)
               endif

               ! ice
               if ( microp_scheme == 'MG' ) then
                  tmp1i = abarii + bbarii/max(13._r8,min(rei(i,k),130._r8))
               else
                  tmp1i = abarii + bbarii/rei(i,k)
               endif

               tauxcl(i,k) = cliqwp(i,k)*tmp1l
               tauxci(i,k) = cicewp(i,k)*tmp1i

	       if (ns .eq. idx_sw_diag) then
               	  tauxcl_out(i,k) = cliqwp(i,k)*tmp1l
               	  tauxci_out(i,k) = cicewp(i,k)*tmp1i
	       endif
            else
               tauxcl(i,k) = 0.0_r8
               tauxci(i,k) = 0.0_r8
            endif

            ! Do not let single scatter albedo be 1.  Delta-eddington solution
            ! for non-conservative case has different analytic form from solution
            ! for conservative case, and raddedmx is written for non-conservative case.
            wcl(i,k) = min(tmp2l,.999999_r8)
            gcl(i,k) = ebarli + tmp3l
            fcl(i,k) = gcl(i,k)*gcl(i,k)

            wci(i,k) = min(tmp2i,.999999_r8)
            gci(i,k) = ebarii + tmp3i
            fci(i,k) = gci(i,k)*gci(i,k)

         end do ! End do i=1,Nday
      end do    ! End do k=1,pver


! 
! Set reflectivities for surface based on mid-point wavelength
! 
      call get_sw_spectral_boundaries(wavmin, wavmax, 'micrometer')
      wavmid(ns) = 0.5_r8*(wavmin(ns) + wavmax(ns))
! 
! Wavelength less  than 0.7 micro-meter
! 
      if (wavmid(ns) < 0.7_r8 ) then
         do i=1,Nday
               albdir(i,ns) = asdir(i)
               albdif(i,ns) = asdif(i)
         end do
! 
! Wavelength greater than 0.7 micro-meter
! 
      else
         do i=1,Nday
               albdir(i,ns) = aldir(i)
               albdif(i,ns) = aldif(i)
         end do
      end if
      trayoslp = raytau(ns)/sslp
! 
! Layer input properties now completely specified; compute the
! delta-Eddington solution reflectivities and transmissivities
! for each layer
! 
      call raddedmx(coszrs   ,Nday    , &
              abh2o(ns),abo3(ns) ,abco2(ns),abo2(ns) , &
              uh2o     ,uo3      ,uco2     ,uo2      , &
              trayoslp ,pflx     ,ns       , &
              tauxcl   ,wcl      ,gcl      ,fcl      , &
              tauxci   ,wci      ,gci      ,fci      , &
              aer_tau(:,:,ns) ,aer_tau_w(:,:,ns)   ,aer_tau_w_g(:,:,ns) ,aer_tau_w_f(:,:,ns) , &
              rdir     ,rdif     ,tdir     ,tdif     ,explay  , &
              rdirc    ,rdifc    ,tdirc    ,tdifc    ,explayc )
! 
! End spectral loop
! 
   end do
! 
!----------------------------------------------------------------------
! 
! Solution for max/random cloud overlap.  
! 
! Steps:
! (1. delta-Eddington solution for each layer (called above)
! 
! (2. The adding method is used to
! compute the reflectivity and transmissivity to direct and diffuse
! radiation from the top and bottom of the atmosphere for each
! cloud configuration.  This calculation is based upon the
! max-random overlap assumption.
! 
! (3. to solve for the fluxes, combine the
! bulk properties of the atmosphere above/below the region.
! 
! Index calculations for steps 2-3 are performed outside spectral
! loop to avoid redundant calculations.  Index calculations (with
! application of areamin & nconfgmax conditions) are performed 
! first to identify the minimum subset of terms for the configurations 
! satisfying the areamin & nconfgmax conditions. This minimum set is 
! used to identify the corresponding minimum subset of terms in 
! steps 2 and 3.
! 
   do iconfig = 1, nconfgmax
      ccon(iconfig,0,1:Nday)      = 0
      ccon(iconfig,pverp,1:Nday)  = 0

      icond(iconfig,0,1:Nday)     = iconfig
      iconu(iconfig,pverp,1:Nday) = iconfig
   end do
! 
! Construction of nuniqu/d, istrtu/d, iconu/d using binary tree 
! 
         nuniqd(0,1:Nday) = 1
         nuniqu(pverp,1:Nday) = 1

         istrtd(1,0,1:Nday) = 1
         istrtu(1,pverp,1:Nday) = 1


!CSD$ PARALLEL DO PRIVATE( npasses, kx2, mrgn, region_found, k1, k2, kx1, nxs, ksort, asort ) &
!CSD$ PRIVATE ( ktmp, atmp, cstr, mstr, nstr, cld0, wstr, nrgn, nconfigm, istr, new_term, xwgt ) &
!CSD$ PRIVATE ( j, ptrc, wgtv, km1, nuniq, is0, is1, n0, n1, ptr0, ptr1, kp1, i, irgn ) &
!CSD$ PRIVATE ( k, l, iconfig, l0, isn )
   do i=1,Nday

!----------------------------------------------------------------------
! INDEX CALCULATIONS FOR MAX OVERLAP
! 
! The column is divided into sets of adjacent layers, called regions, 
! in which the clouds are maximally overlapped.  The clouds are
! randomly overlapped between different regions.  The number of
! regions in a column is set by nmxrgn, and the range of pressures
! included in each region is set by pmxrgn.  
! 
! The following calculations determine the number of unique cloud 
! configurations (assuming maximum overlap), called "streams",
! within each region. Each stream consists of a vector of binary
! clouds (either 0 or 100% cloud cover).  Over the depth of the region, 
! each stream requires a separate calculation of radiative properties. These
! properties are generated using the adding method from
! the radiative properties for each layer calculated by raddedmx.
! 
! The upward and downward-propagating streams are treated
! separately.
! 
! We will refer to a particular configuration of binary clouds
! within a single max-overlapped region as a "stream".  We will 
! refer to a particular arrangement of binary clouds over the entire column
! as a "configuration".
! 
! This section of the code generates the following information:
! (1. nrgn    : the true number of max-overlap regions (need not = nmxrgn)
! (2. nstr    : the number of streams in a region (>=1)
! (3. cstr    : flags for presence of clouds at each layer in each stream
! (4. wstr    : the fractional horizontal area of a grid box covered
! by each stream
! (5. kx1,2   : level indices for top/bottom of each region
! 
! The max-overlap calculation proceeds in 3 stages:
! (1. compute layer radiative properties in raddedmx.
! (2. combine these properties between layers 
! (3. combine properties to compute fluxes at each interface.  
! 
! Most of the indexing information calculated here is used in steps 2-3
! after the call to raddedmx.
! 
! Initialize indices for layers to be max-overlapped
! 
! Loop to handle fix in totwgt=0. For original overlap config 
! from npasses = 0.
! 
         npasses = 0
         do
!cdir novector
            do irgn = 0, nmxrgn(i)
               kx2(irgn) = 0
            end do
            mrgn = 0
! 
! Outermost loop over regions (sets of adjacent layers) to be max overlapped
! 
            do irgn = 1, nmxrgn(i)
! 
! Calculate min/max layer indices inside region.  
! 
               region_found = .false.
               if (kx2(irgn-1) < pver) then
                  k1 = kx2(irgn-1)+1
                  kx1(irgn) = k1
                  kx2(irgn) = k1-1
!cdir novector
                  do k2 = pver, k1, -1
                     if (pmid(i,k2) <= pmxrgn(i,irgn)) then
                        kx2(irgn) = k2
                        mrgn = mrgn+1
                        region_found = .true.
                        exit
                     end if
                  end do
               else
                  exit
               endif

               if (region_found) then
! 
! Sort cloud areas and corresponding level indices.  
! 
                  nxs = 0
                  if (cldeps > 0) then 
                     do k = k1,k2
                        if (cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then
                           nxs = nxs+1
                           ksort(nxs) = k
! 
! We need indices for clouds in order of largest to smallest, so
! sort 1-cld in ascending order
! 
                           asort(nxs) = 1.0_r8-(floor(cld(i,k)/cldeps)*cldeps)
                        end if
                     end do
                  else
!cdir novector
                     do k = k1,k2
                        if (cld(i,k) >= cldmin) then
                           nxs = nxs+1
                           ksort(nxs) = k
! 
! We need indices for clouds in order of largest to smallest, so
! sort 1-cld in ascending order
! 
                           asort(nxs) = 1.0_r8-cld(i,k)
                        end if
                     end do
                  endif
! 
! If nxs eq 1, no need to sort. 
! If nxs eq 2, sort by swapping if necessary
! If nxs ge 3, sort using local sort routine
! 
                  if (nxs == 2) then
                     if (asort(2) < asort(1)) then
                        ktmp = ksort(1)
                        ksort(1) = ksort(2)
                        ksort(2) = ktmp

                        atmp = asort(1)
                        asort(1) = asort(2)
                        asort(2) = atmp
                     endif
                  else if (nxs >= 3) then
                     call quick_sort(asort(1:nxs),ksort(1:nxs))
                  endif
! 
! Construct wstr, cstr, nstr for this region
! 
!cdir novector
                  cstr(k1:k2,1:nxs+1) = 0
                  mstr = 1
                  cld0 = 0.0_r8
                  do l = 1, nxs
                     if (asort(l) /= cld0) then
                        wstr(mstr,mrgn) = asort(l) - cld0
                        cld0 = asort(l)
                        mstr = mstr + 1
                     endif
!cdir novector
                     cstr(ksort(l),mstr:nxs+1) = 1
                  end do
                  nstr(mrgn) = mstr
                  wstr(mstr,mrgn) = 1.0_r8 - cld0
! 
! End test of region_found = true
! 
               endif
! 
! End loop over regions irgn for max-overlap
! 
            end do
            nrgn = mrgn
! 
! Finish construction of cstr for additional top layer
! 
!cdir novector
            cstr(0,1:nstr(1)) = 0
! 
! INDEX COMPUTATIONS FOR STEP 2-3
! This section of the code generates the following information:
! (1. totwgt     step 3     total frac. area of configurations satisfying
! areamin & nconfgmax criteria
! (2. wgtv       step 3     frac. area of configurations 
! (3. ccon       step 2     binary flag for clouds in each configuration
! (4. nconfig    steps 2-3  number of configurations
! (5. nuniqu/d   step 2     Number of unique cloud configurations for
! up/downwelling rad. between surface/TOA
! and level k
! (6. istrtu/d   step 2     Indices into iconu/d
! (7. iconu/d    step 2     Cloud configurations which are identical
! for up/downwelling rad. between surface/TOA
! and level k
! 
! Number of configurations (all permutations of streams in each region)
! 
            nconfigm = product(nstr(1: nrgn))
! 
! Construction of totwgt, wgtv, ccon, nconfig
! 
!cdir novector
            istr(1: nrgn) = 1
            nconfig(i) = 0
            totwgt(i) = 0.0_r8
            new_term = .true.
            do iconfig = 1, nconfigm
               xwgt = 1.0_r8
!cdir novector
               do mrgn = 1,  nrgn
                  xwgt = xwgt * wstr(istr(mrgn),mrgn)
               end do
               if (xwgt >= areamin) then
                  nconfig(i) = nconfig(i) + 1
                  if (nconfig(i) <= nconfgmax) then
                     j = nconfig(i)
                     ptrc(nconfig(i)) = nconfig(i)
                  else
                     nconfig(i) = nconfgmax
                     if (new_term) then
                        min_idx = minloc(wgtv)
                        j = min_idx(1)
                     endif
                     if (wgtv(j) < xwgt) then
                        totwgt(i) = totwgt(i) - wgtv(j)
                        new_term = .true.
                     else
                        new_term = .false.
                     endif
                  endif
                  if (new_term) then
                     wgtv(j) = xwgt
                     totwgt(i) = totwgt(i) + xwgt
!cdir novector
                     do mrgn = 1, nrgn
!cdir novector
                        ccon(j,kx1(mrgn):kx2(mrgn),i) = cstr(kx1(mrgn):kx2(mrgn),istr(mrgn))
                     end do
                  endif
               endif

               mrgn =  nrgn
               istr(mrgn) = istr(mrgn) + 1
               do while (istr(mrgn) > nstr(mrgn) .and. mrgn > 1)
                  istr(mrgn) = 1
                  mrgn = mrgn - 1
                  istr(mrgn) = istr(mrgn) + 1
               end do
! 
! End do iconfig = 1, nconfigm
! 
            end do
! 
! If totwgt(i) = 0 implement maximum overlap and make another pass
! if totwgt(i) = 0 on this second pass then terminate.
! 
            if (totwgt(i) > 0._r8) then
               exit
            else
               npasses = npasses + 1
               if (npasses >= 2 ) then
                  write(iulog,*)'RADCSWMX: Maximum overlap of column ','failed'
                  call endrun('RADCSWMX')
               endif
               nmxrgn(i)=1
               pmxrgn(i,1)=1.0e30_r8
            end if
!
! End npasses = 0, do
!
         end do
! 
! Finish construction of ccon
! 

         istrtd(2,0,i) = nconfig(i)+1
         istrtu(2,pverp,i) = nconfig(i)+1

         do k = 1, pverp
            km1 = k-1
            nuniq = 0
            istrtd(1,k,i) = 1
!cdir novector
            do l0 = 1, nuniqd(km1,i)
               is0 = istrtd(l0,km1,i)
               is1 = istrtd(l0+1,km1,i)-1
               n0 = 0
               n1 = 0
!cdir novector
               do isn = is0, is1
                  j = icond(isn,km1,i)
                  if (ccon(j,k,i) == 0) then
                     n0 = n0 + 1
                     ptr0(n0) = j
                  else       ! if (ccon(j,k,i) == 1) then
                     n1 = n1 + 1
                     ptr1(n1) = j
                  endif
               end do
               if (n0 > 0) then
                  nuniq = nuniq + 1
                  istrtd(nuniq+1,k,i) = istrtd(nuniq,k,i)+n0
!cdir novector
                  icond(istrtd(nuniq,k,i):istrtd(nuniq+1,k,i)-1,k,i) =  ptr0(1:n0)
               endif
               if (n1 > 0) then
                  nuniq = nuniq + 1
                  istrtd(nuniq+1,k,i) = istrtd(nuniq,k,i)+n1
!cdir novector
                  icond(istrtd(nuniq,k,i):istrtd(nuniq+1,k,i)-1,k,i) =  ptr1(1:n1)
               endif
            end do
            nuniqd(k,i) = nuniq
         end do
!
!  Find 'transition point' in downward configurations where the number
!  of 'configurations' changes from 1.  This is used to optimize the
!  construction of the upward configurations.
!  Note: k1 == transition point
!

         do k = pverp,0,-1
           if ( nuniqd(k,i) == 1) then
              k1 = k
              exit
           end if
         end do

         do k = pver, k1+1,-1
            kp1 = k+1
            nuniq = 0
            istrtu(1,k,i) = 1
!cdir novector
            do l0 = 1, nuniqu(kp1,i)
               is0 = istrtu(l0,kp1,i)
               is1 = istrtu(l0+1,kp1,i)-1
               n0 = 0
               n1 = 0
!cdir novector
               do isn = is0, is1
                  j = iconu(isn,kp1,i)
                  if (ccon(j,k,i) == 0) then
                     n0 = n0 + 1
                     ptr0(n0) = j
                  else       ! if (ccon(j,k,i) == 1) then
                     n1 = n1 + 1
                     ptr1(n1) = j
                  endif
               end do
               if (n0 > 0) then
                  nuniq = nuniq + 1
                  istrtu(nuniq+1,k,i) = istrtu(nuniq,k,i)+n0
!cdir novector
                  iconu(istrtu(nuniq,k,i):istrtu(nuniq+1,k,i)-1,k,i) =  ptr0(1:n0)
               endif
               if (n1 > 0) then
                  nuniq = nuniq + 1
                  istrtu(nuniq+1,k,i) = istrtu(nuniq,k,i)+n1
!cdir novector
                  iconu(istrtu(nuniq,k,i):istrtu(nuniq+1,k,i)-1,k,i) = ptr1(1:n1)
               endif
            end do
            nuniqu(k,i) = nuniq
         end do
!
!  Copy identical configurations from 'transition point' to surface.
!
         k1 = min(pverp-1,k1)
         nuniq = nuniqu(k1+1,i)
         do k = k1,0,-1
            nuniqu(k,i) = nuniq
!cdir novector
            iconu(1:nuniq,k,i) = iconu(1:nuniq,k1+1,i)
!cdir novector
            istrtu(1:nuniq+1,k,i) = istrtu(1:nuniq+1,k1+1,i)
         end do

!cdir novector
         v_wgtv(1:nconfig(i),i) = wgtv(1:nconfig(i))

! 
!----------------------------------------------------------------------
! End of index calculations
!----------------------------------------------------------------------
! 
! End do i=1,Nday
! 
   end do
!CSD$ END PARALLEL 

!----------------------------------------------------------------------
! Start of flux calculations
!----------------------------------------------------------------------
!
! Initialize spectrally integrated totals:
! 
         totfld (1:Nday,0:pver) = 0.0_r8
         totfldc(1:Nday,0:pver) = 0.0_r8
         fswup  (1:Nday,0:pver) = 0.0_r8
         fswdn  (1:Nday,0:pver) = 0.0_r8

         fswup (1:Nday,pverp)  = 0.0_r8
         fswdn (1:Nday,pverp)  = 0.0_r8
! 
! Start spectral interval
! 
!old   do ns = 1,nswbands
!old     wgtint = nirwgt(ns)

     do i=1,Nday

!----------------------------------------------------------------------
! STEP 2
! 
! 
! Apply adding method to solve for radiative properties
! 
! first initialize the bulk properties at toa
! 

! nswbands, 0:pverp, nconfgmax, pcols

            rdndif(:,0,1:nconfig(i),i) = 0.0_r8
            exptdn(:,0,1:nconfig(i),i) = 1.0_r8
            tdntot(:,0,1:nconfig(i),i) = 1.0_r8
! 
! End do i=1,Nday
! 
     end do
! 
! solve for properties involving downward propagation of radiation.
! the bulk properties are:
! 
! (1. exptdn   sol. beam dwn. trans from layers above
! (2. rdndif   ref to dif rad for layers above
! (3. tdntot   total trans for layers above
! 

!CSD$ PARALLEL DO PRIVATE( km1, is0, is1, j, jj, Ttdif, Trdif, Trdir, Ttdir, Texplay ) &
!CSD$ PRIVATE( xexpt, xrdnd, tdnmexp,  ytdnd, yrdnd, rdenom, rdirexp, zexpt, zrdnd, ztdnt ) &
!CSD$ PRIVATE( i, k, l0, ns, isn )
         do i = 1, Nday
            do k = 1, pverp
               km1 = k - 1
!cdir nodep
               do l0 = 1, nuniqd(km1,i)
                  is0 = istrtd(l0,km1,i)
                  is1 = istrtd(l0+1,km1,i)-1

                  j = icond(is0,km1,i)

! 
! If cloud in layer, use cloudy layer radiative properties (ccon == 1)
! If clear layer, use clear-sky layer radiative properties (ccon /= 1)
! 
                  if ( ccon(j,km1,i) == 1 ) then
                     Ttdif(:) = tdif(:,i,km1)
                     Trdif(:) = rdif(:,i,km1)
                     Trdir(:) = rdir(:,i,km1)
                     Ttdir(:) = tdir(:,i,km1)
                     Texplay(:) = explay(:,i,km1)
                  else
                     Ttdif(:) = tdifc(:,i,km1)
                     Trdif(:) = rdifc(:,i,km1)
                     Trdir(:) = rdirc(:,i,km1)
                     Ttdir(:) = tdirc(:,i,km1)
                     Texplay(:) = explayc(:,i,km1)
                  end if

                  do ns = 1, nswbands
                  xexpt   = exptdn(ns,km1,j,i)
                  xrdnd   = rdndif(ns,km1,j,i)
                  tdnmexp = tdntot(ns,km1,j,i) - xexpt

                  ytdnd = Ttdif(ns)
                  yrdnd = Trdif(ns)

                  rdenom  = 1._r8/(1._r8-yrdnd*xrdnd)
                  rdirexp = Trdir(ns)*xexpt

                  zexpt = xexpt * Texplay(ns)
                  zrdnd = yrdnd + xrdnd*(ytdnd**2)*rdenom
                  ztdnt = xexpt*Ttdir(ns) + ytdnd* &
                          (tdnmexp + xrdnd*rdirexp)*rdenom

                  exptdn(ns,k,j,i) = zexpt
                  rdndif(ns,k,j,i) = zrdnd
                  tdntot(ns,k,j,i) = ztdnt
                  end do ! ns = 1, nswbands
!
! If 2 or more configurations share identical properties at a given level k,
! the properties (at level k) are computed once and copied to
! all the configurations for efficiency.
!
                  do isn = is0+1, is1
                     jj = icond(isn,km1,i)
                     exptdn(:,k,jj,i) = exptdn(:,k,j,i)
                     rdndif(:,k,jj,i) = rdndif(:,k,j,i)
                     tdntot(:,k,jj,i) = tdntot(:,k,j,i)
                  end do

! 
! end do l0 = 1, nuniqd(k,i)
! 
               end do
! 
! end do k = 1, pverp
! 
            end do
! 
! end do i = 1, Nday
! 
         end do
!CSD$ END PARALLEL 
! 
! Solve for properties involving upward propagation of radiation.
! The bulk properties are:
! 
! (1. rupdif   Ref to dif rad for layers below
! (2. rupdir   Ref to dir rad for layers below
! 
! Specify surface boundary conditions (surface albedos)
! 

! nswbands, 0:pverp, nconfgmax, pcols
   rupdir = 0._r8
   rupdif = 0._r8
   do i = 1, Nday
      do ns = 1, nswbands
         rupdir(ns,pverp,1:nconfig(i),i) = albdir(i,ns)
         rupdif(ns,pverp,1:nconfig(i),i) = albdif(i,ns)
      end do
   end do

         do i = 1, Nday
            do k = pver, 0, -1
               do l0 = 1, nuniqu(k,i)
                  is0 = istrtu(l0,k,i)
                  is1 = istrtu(l0+1,k,i)-1

                  j = iconu(is0,k,i)

! 
! If cloud in layer, use cloudy layer radiative properties (ccon == 1)
! If clear layer, use clear-sky layer radiative properties (ccon /= 1)
! 
                  if ( ccon(j,k,i) == 1 ) then
                     Ttdif(:) = tdif(:,i,k)
                     Trdif(:) = rdif(:,i,k)
                     Trdir(:) = rdir(:,i,k)
                     Ttdir(:) = tdir(:,i,k)
                     Texplay(:) = explay(:,i,k)
                  else
                     Ttdif(:) = tdifc(:,i,k)
                     Trdif(:) = rdifc(:,i,k)
                     Trdir(:) = rdirc(:,i,k)
                     Ttdir(:) = tdirc(:,i,k)
                     Texplay(:) = explayc(:,i,k)
                  end if

                  do ns = 1, nswbands
                  xrupd = rupdif(ns,k+1,j,i)
                  xrups = rupdir(ns,k+1,j,i)

! 
! If cloud in layer, use cloudy layer radiative properties (ccon == 1)
! If clear layer, use clear-sky layer radiative properties (ccon /= 1)
! 
                  yexpt = Texplay(ns)
                  yrupd = Trdif(ns)
                  ytupd = Ttdif(ns)

                  rdenom  = 1._r8/( 1._r8 - yrupd*xrupd)
                  tdnmexp = (Ttdir(ns)-yexpt)
                  rdirexp = xrups*yexpt

                  zrupd = yrupd + xrupd*(ytupd**2)*rdenom
                  zrups = Trdir(ns) + ytupd*(rdirexp + xrupd*tdnmexp)*rdenom

                  rupdif(ns,k,j,i) = zrupd
                  rupdir(ns,k,j,i) = zrups
                  end do ! ns = 1, nswbands
!
! If 2 or more configurations share identical properties at a given level k,
! the properties (at level k) are computed once and copied to
! all the configurations for efficiency.
!
                  do isn = is0+1, is1
                     jj = iconu(isn,k,i)
                     rupdif(:,k,jj,i) = rupdif(:,k,j,i)
                     rupdir(:,k,jj,i) = rupdir(:,k,j,i)
                  end do

! 
! end do l0 = 1, nuniqu(k,i)
! 
               end do
! 
! end do k = pver,0,-1
! 
            end do
! 
! end do i = 1, Nday
! 
         end do
! 
!----------------------------------------------------------------------
! 
! STEP 3
! 
! Compute up and down fluxes for each interface k.  This requires
! adding up the contributions from all possible permutations
! of streams in all max-overlap regions, weighted by the
! product of the fractional areas of the streams in each region
! (the random overlap assumption).  The adding principle has been
! used in step 2 to combine the bulk radiative properties 
! above and below the interface.
! 

! 
! Initialize the fluxes
! 
            fluxup = 0.0_r8
            fluxdn = 0.0_r8

            do i = 1, Nday
!cdir novector
            do iconfig = 1, nconfig(i)
               xwgt = v_wgtv(iconfig,i)

!cdir collapse
               do k = 0, pverp
                  do ns = 1, nswbands
                  xexpt = exptdn(ns,k,iconfig,i)
                  xtdnt = tdntot(ns,k,iconfig,i)
                  xrdnd = rdndif(ns,k,iconfig,i)
                  xrupd = rupdif(ns,k,iconfig,i)
                  xrups = rupdir(ns,k,iconfig,i)
! 
! Flux computation
! 
                  rdenom = 1._r8/(1._r8 - xrdnd * xrupd)

                  fluxup(ns,k,i) = fluxup(ns,k,i) + xwgt *  &
                              ((xexpt * xrups + (xtdnt - xexpt) * xrupd) * rdenom)
                  fluxdn(ns,k,i) = fluxdn(ns,k,i) + xwgt *  &
                              (xexpt + (xtdnt - xexpt + xexpt * xrups * xrdnd) * rdenom)
                  end do ! do ns = 1, nswbands
               end do
! 
! End do iconfig = 1, nconfig(i)
! 
            end do
! 
! End do iconfig = 1, Nday
! 
            end do
! 
! Normalize by total area covered by cloud configurations included
! in solution
! 
            do i = 1, Nday
            do k = 0, pverp
            do ns = 1, nswbands
               fluxup(ns,k,i)=fluxup(ns,k,i) / totwgt(i)
               fluxdn(ns,k,i)=fluxdn(ns,k,i) / totwgt(i)
            end do ! do i = 1, nday
            end do ! do k = 0, pverp
            end do ! do i = 1, nday


! 
! Initialize the direct-beam flux at surface
! 
            wexptdn(:,1:Nday) = 0.0_r8

   do ns = 1,nswbands
     wgtint = nirwgt(ns)


            do i=1,Nday
            do iconfig = 1, nconfig(i)
!
! Note: exptdn can be directly indexed by iconfig at k=pverp.
!
               wexptdn(ns,i) =  wexptdn(ns,i) + v_wgtv(iconfig,i) * exptdn(ns,pverp,iconfig,i)
            end do
            end do

            do i=1,Nday
               wexptdn(ns,i) = wexptdn(ns,i) / totwgt(i)
! 
! Monochromatic computation completed; accumulate in totals
! 
            if ( do_spctrl_scaling ) then
               solflx(i)   = solin(i)*frcsol(ns)*psf(ns)*sfac(ns)
            else
               solflx(i)   = solin(i)*frcsol(ns)*psf(ns) 
            endif
            fsnt(i)  = fsnt(i) + solflx(i)*(fluxdn(ns,1,i) - fluxup(ns,1,i))
            fsntoa(i)= fsntoa(i) + solflx(i)*(fluxdn(ns,0,i) - fluxup(ns,0,i))
            fsutoa(i)= fsutoa(i) + solflx(i)*(fluxup(ns,0,i))
            fsns(i)  = fsns(i) + solflx(i)*(fluxdn(ns,pverp,i)-fluxup(ns,pverp,i))
            fsdtoa(i)  = fsdtoa(i) + solflx(i)
            fswup(i,0) = fswup(i,0) + solflx(i)*fluxup(ns,0,i)
            fswdn(i,0) = fswdn(i,0) + solflx(i)*fluxdn(ns,0,i)
! 
! Down spectral fluxes need to be in mks; thus the .001 conversion factors
! 
            if (wavmid(ns) < 0.7_r8) then
               sols(i)  = sols(i) + wexptdn(ns,i)*solflx(i)*0.001_r8
               solsd(i) = solsd(i)+(fluxdn(ns,pverp,i)-wexptdn(ns,i))*solflx(i)*0.001_r8 
            else
               soll(i)  = soll(i) + wexptdn(ns,i)*solflx(i)*0.001_r8
               solld(i) = solld(i)+(fluxdn(ns,pverp,i)-wexptdn(ns,i))*solflx(i)*0.001_r8 
               fsnrtoaq(i) = fsnrtoaq(i) + solflx(i)*(fluxdn(ns,0,i) - fluxup(ns,0,i))
            end if
            fsnirtoa(i) = fsnirtoa(i) + wgtint*solflx(i)*(fluxdn(ns,0,i) - fluxup(ns,0,i))

! 
! End do i=1,Nday
! 
   end do


            do k=0,pver
            do i=1,Nday
! 
! Compute flux divergence in each layer using the interface up and down
! fluxes:
! 
               kp1 = k+1
               flxdiv = (fluxdn(ns,k,i) - fluxdn(ns,kp1,i)) + (fluxup(ns,kp1,i) - fluxup(ns,k,i))
               totfld(i,k)  = totfld(i,k)  + solflx(i)*flxdiv
               fswdn(i,kp1) = fswdn(i,kp1) + solflx(i)*fluxdn(ns,kp1,i)
               fswup(i,kp1) = fswup(i,kp1) + solflx(i)*fluxup(ns,kp1,i)
               fns(i,kp1)   = fswdn(i,kp1) - fswup(i,kp1)
               if (single_column.and.scm_crm_mode) then 
                  fus(i,kp1)=fswup(i,kp1)
                  fds(i,kp1)=fswdn(i,kp1)
               endif
            end do
            end do
! 
! Perform clear-sky calculation
! 

            exptdnc(1:Nday,0) =   1.0_r8
            rdndifc(1:Nday,0) =   0.0_r8
            tdntotc(1:Nday,0) =   1.0_r8
            rupdirc(1:Nday,pverp) = albdir(1:Nday,ns)
            rupdifc(1:Nday,pverp) = albdif(1:Nday,ns)


!cdir expand=pverp
            do k = 1, pverp
            do i=1,Nday
               km1 = k - 1
               xexpt = exptdnc(i,km1)
               xrdnd = rdndifc(i,km1)
               yrdnd = rdifc(ns,i,km1)
               ytdnd = tdifc(ns,i,km1)

               exptdnc(i,k) = xexpt*explayc(ns,i,km1)

               rdenom  = 1._r8/(1._r8 - yrdnd*xrdnd)
               rdirexp = rdirc(ns,i,km1)*xexpt
               tdnmexp = tdntotc(i,km1) - xexpt

               tdntotc(i,k) = xexpt*tdirc(ns,i,km1) + ytdnd*(tdnmexp + xrdnd*rdirexp)* &
                                rdenom
               rdndifc(i,k) = yrdnd + xrdnd*(ytdnd**2)*rdenom
! 
! End do i=1,Nday
! 
            end do
            end do

            do k=pver,0,-1
            do i=1,Nday
               xrupd = rupdifc(i,k+1)
               yexpt = explayc(ns,i,k)
               yrupd = rdifc(ns,i,k)
               ytupd = tdifc(ns,i,k)

               rdenom = 1._r8/( 1._r8 - yrupd*xrupd)

               rupdirc(i,k) = rdirc(ns,i,k) + ytupd*(rupdirc(i,k+1)*yexpt + &
                            xrupd*(tdirc(ns,i,k)-yexpt))*rdenom
               rupdifc(i,k) = yrupd + xrupd*ytupd**2*rdenom
! 
! End do i=1,Nday
! 
            end do
            end do

            do k=0,pverp
            do i=1,Nday
               rdenom    = 1._r8/(1._r8 - rdndifc(i,k)*rupdifc(i,k))
               fluxup(ns,k,i) = (exptdnc(i,k)*rupdirc(i,k) + (tdntotc(i,k)-exptdnc(i,k))*rupdifc(i,k))* &
                           rdenom
               fluxdn(ns,k,i) = exptdnc(i,k) + &
                           (tdntotc(i,k) - exptdnc(i,k) + exptdnc(i,k)*rupdirc(i,k)*rdndifc(i,k))* &
                           rdenom
! 
! End do i=1,Nday
! 
            end do
            end do

            do k=0,pver
            do i=1,Nday
! 
! Compute flux divergence in each layer using the interface up and down
! fluxes:
! 
               kp1 = k+1
               flxdiv        = (fluxdn(ns,k,i) - fluxdn(ns,kp1,i)) + (fluxup(ns,kp1,i) - fluxup(ns,k,i))
               totfldc(i,k)  = totfldc(i,k)  + solflx(i)*flxdiv
            end do
            end do

            do i=1,Nday
            fsntc(i)    = fsntc(i)+solflx(i)*(fluxdn(ns,1,i)-fluxup(ns,1,i))
            fsntoac(i)  = fsntoac(i)+solflx(i)*(fluxdn(ns,0,i)-fluxup(ns,0,i))
            fsnsc(i)    = fsnsc(i)+solflx(i)*(fluxdn(ns,pverp,i)-fluxup(ns,pverp,i))
            fsdsc(i)    = fsdsc(i)+solflx(i)*(fluxdn(ns,pverp,i))
            fsnrtoac(i) = fsnrtoac(i)+wgtint*solflx(i)*(fluxdn(ns,0,i)-fluxup(ns,0,i))
            if (single_column.and.scm_crm_mode) then 
               do k = 1,pverp
                  fusc(i,k)=fusc(i,k) + solflx(i) * fluxup(ns,k,i)
                  fdsc(i,k)=fdsc(i,k) + solflx(i) * fluxdn(ns,k,i)
               enddo
            endif
! 
! End do i=1,Nday
! 
            end do

            do k = 1,pverp
            do i=1,Nday
               fcns(i,k)=fcns(i,k) + solflx(i)*(fluxdn(ns,k,i)-fluxup(ns,k,i))
            enddo
            enddo
! 
! End of clear sky calculation
! 

! 
! End of spectral interval loop
! 
         end do

   do i=1,Nday

! 
! Compute solar heating rate (J/kg/s)
! 
!cdir expand=pver
         do k=1,pver
            qrs(i,k) = -1.E-4_r8*gravit*totfld(i,k)/(pint(i,k) - pint(i,k+1))
            qrsc(i,k) = -1.E-4_r8*gravit*totfldc(i,k)/(pint(i,k) - pint(i,k+1))
         end do
! 
! Set the downwelling flux at the surface 
! 
         fsds(i) = fswdn(i,pverp)
! 
! End do i=1,Nday
! 
   end do
!
! Rearrange output arrays.
!
! intent(inout)
!
   call ExpDayNite(pmxrgn,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pverp)
   call ExpDayNite(nmxrgn,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
!
! intent(out)
!
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
   call ExpDayNite(fsdtoa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsntoac,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(sols,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(soll,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(solsd,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(solld,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnirtoa,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnrtoac,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(fsnrtoaq,	Nday, IdxDay, Nnite, IdxNite, 1, pcols)
   call ExpDayNite(tauxcl_out,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)
   call ExpDayNite(tauxci_out,	Nday, IdxDay, Nnite, IdxNite, 1, pcols, 1, pver)

!  these outfld calls don't work for spmd only outfield in scm mode (nonspmd)
   if (single_column.and.scm_crm_mode) then 
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
!  write(iulog, '(a, x, i3)') 'radcswmx : exiting, chunk identifier', lchnk

   return
end subroutine radcswmx

!-------------------------------------------------------------------------------

subroutine raddedmx(coszrs  ,ndayc   ,abh2o   , &
                    abo3    ,abco2   ,abo2    ,uh2o    ,uo3     , &
                    uco2    ,uo2     ,trayoslp,pflx    ,ns      , &
                    tauxcl  ,wcl     ,gcl     ,fcl     ,tauxci  , &
                    wci     ,gci     ,fci     ,aer_tau ,aer_tau_w, &
                    aer_tau_w_g, aer_tau_w_f      ,rdir    ,rdif    ,tdir    , &
                    tdif    ,explay  ,rdirc   ,rdifc   ,tdirc   , &
                    tdifc   ,explayc )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes layer reflectivities and transmissivities, from the top down
! to the surface using the delta-Eddington solutions for each layer
! 
! Method: 
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
! Modified for maximum/random cloud overlap by Bill Collins and John
!    Truesdale
! 
! Author: Bill Collins
! 
!-----------------------------------------------------------------------

!
! Minimum total transmission below which no layer computation are done:
!
   real(r8) trmin                ! Minimum total transmission allowed
   real(r8) wray                 ! Rayleigh single scatter albedo
   real(r8) gray                 ! Rayleigh asymetry parameter
   real(r8) fray                 ! Rayleigh forward scattered fraction

   parameter (trmin = 1.e-3_r8)
   parameter (wray = 0.999999_r8)
   parameter (gray = 0.0_r8)
   parameter (fray = 0.1_r8)
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: coszrs(pcols)        ! Cosine zenith angle
   real(r8), intent(in) :: trayoslp             ! Tray/sslp
   real(r8), intent(in) :: pflx(pcols,0:pverp)  ! Interface pressure
   real(r8), intent(in) :: abh2o                ! Absorption coefficiant for h2o
   real(r8), intent(in) :: abo3                 ! Absorption coefficiant for o3
   real(r8), intent(in) :: abco2                ! Absorption coefficiant for co2
   real(r8), intent(in) :: abo2                 ! Absorption coefficiant for o2
   real(r8), intent(in) :: uh2o(pcols,0:pver)   ! Layer absorber amount of h2o
   real(r8), intent(in) :: uo3(pcols,0:pver)    ! Layer absorber amount of  o3
   real(r8), intent(in) :: uco2(pcols,0:pver)   ! Layer absorber amount of co2
   real(r8), intent(in) :: uo2(pcols,0:pver)    ! Layer absorber amount of  o2
   real(r8), intent(in) :: tauxcl(pcols,0:pver) ! Cloud extinction optical depth (liquid)
   real(r8), intent(in) :: wcl(pcols,0:pver)    ! Cloud single scattering albedo (liquid)
   real(r8), intent(in) :: gcl(pcols,0:pver)    ! Cloud asymmetry parameter (liquid)
   real(r8), intent(in) :: fcl(pcols,0:pver)    ! Cloud forward scattered fraction (liquid)
   real(r8), intent(in) :: tauxci(pcols,0:pver) ! Cloud extinction optical depth (ice)
   real(r8), intent(in) :: wci(pcols,0:pver)    ! Cloud single scattering albedo (ice)
   real(r8), intent(in) :: gci(pcols,0:pver)    ! Cloud asymmetry parameter (ice)
   real(r8), intent(in) :: fci(pcols,0:pver)    ! Cloud forward scattered fraction (ice)
   real(r8), intent(inout) :: aer_tau(pcols,0:pver) ! Aerosol extinction optical depth
   real(r8), intent(inout) :: aer_tau_w(pcols,0:pver)     ! Aerosol single scattering albedo * tau
   real(r8), intent(inout) :: aer_tau_w_g(pcols,0:pver)     ! Aerosol asymmetry parameter * w * t
   real(r8), intent(inout) :: aer_tau_w_f(pcols,0:pver)     ! Aerosol forward scattered fraction * w * tau

   integer, intent(in) :: ndayc                 ! Number of daylight columns
   integer, intent(in) :: ns                    ! Index of spectral interval
!
! Input/Output arguments
!
! Following variables are defined for each layer; 0 refers to extra
! layer above top of model:
!
   real(r8), intent(inout) :: rdir(nswbands,pcols,0:pver)   ! Layer reflectivity to direct rad
   real(r8), intent(inout) :: rdif(nswbands,pcols,0:pver)   ! Layer reflectivity to diffuse rad
   real(r8), intent(inout) :: tdir(nswbands,pcols,0:pver)   ! Layer transmission to direct rad
   real(r8), intent(inout) :: tdif(nswbands,pcols,0:pver)   ! Layer transmission to diffuse rad
   real(r8), intent(inout) :: explay(nswbands,pcols,0:pver) ! Solar beam exp transm for layer
!
! Corresponding quantities for clear-skies
!
   real(r8), intent(inout) :: rdirc(nswbands,pcols,0:pver)  ! Clear layer reflec. to direct rad
   real(r8), intent(inout) :: rdifc(nswbands,pcols,0:pver)  ! Clear layer reflec. to diffuse rad
   real(r8), intent(inout) :: tdirc(nswbands,pcols,0:pver)  ! Clear layer trans. to direct rad
   real(r8), intent(inout) :: tdifc(nswbands,pcols,0:pver)  ! Clear layer trans. to diffuse rad
   real(r8), intent(inout) :: explayc(nswbands,pcols,0:pver)! Solar beam exp transm clear layer
!
!---------------------------Local variables-----------------------------
!
   integer i                 ! Column indices
   integer k                 ! Level index
   integer nn                ! Index of column loops (max=ndayc)

   real(r8) taugab(pcols)        ! Layer total gas absorption optical depth
   real(r8) tauray(pcols)        ! Layer rayleigh optical depth
   real(r8) taucsc               ! Layer cloud scattering optical depth
   real(r8) tautot               ! Total layer optical depth
   real(r8) wtot                 ! Total layer single scatter albedo
   real(r8) gtot                 ! Total layer asymmetry parameter
   real(r8) ftot                 ! Total layer forward scatter fraction
   real(r8) wtau                 !  rayleigh layer scattering optical depth
   real(r8) wt                   !  layer total single scattering albedo
   real(r8) ts                   !  layer scaled extinction optical depth
   real(r8) ws                   !  layer scaled single scattering albedo
   real(r8) gs                   !  layer scaled asymmetry parameter
!
!---------------------------Statement functions-------------------------
!
! Statement functions and other local variables
!
   real(r8) alpha                ! Term in direct reflect and transmissivity
   real(r8) gamma                ! Term in direct reflect and transmissivity
   real(r8) el                   ! Term in alpha,gamma,n,u
   real(r8) taus                 ! Scaled extinction optical depth
   real(r8) omgs                 ! Scaled single particle scattering albedo
   real(r8) asys                 ! Scaled asymmetry parameter
   real(r8) u                    ! Term in diffuse reflect and
!    transmissivity
   real(r8) n                    ! Term in diffuse reflect and
!    transmissivity
   real(r8) lm                   ! Temporary for el
   real(r8) ne                   ! Temporary for n
   real(r8) w                    ! Dummy argument for statement function
   real(r8) uu                   ! Dummy argument for statement function
   real(r8) g                    ! Dummy argument for statement function
   real(r8) e                    ! Dummy argument for statement function
   real(r8) f                    ! Dummy argument for statement function
   real(r8) t                    ! Dummy argument for statement function
   real(r8) et                   ! Dummy argument for statement function
!
! Intermediate terms for delta-eddington solution
!
   real(r8) alp                  ! Temporary for alpha
   real(r8) gam                  ! Temporary for gamma
   real(r8) ue                   ! Temporary for u
   real(r8) arg                  ! Exponential argument
   real(r8) extins               ! Extinction
   real(r8) amg                  ! Alp - gam
   real(r8) apg                  ! Alp + gam
!
! ssa <=1 limit for aerosol
!
   real(r8) :: w_limited(pcols,0:pver)       ! Aerosol ssa (limited to < 0.999999)
   real(r8) :: aer_g_limit(pcols,0:pver)     ! Aerosol tau_w_g (limited ssa)
   real(r8) :: aer_f_limit(pcols,0:pver)     ! Aerosol tau_w_f (limited ssa)
!
   alpha(w,uu,g,e) = .75_r8*w*uu*((1._r8 + g*(1._r8-w))/(1._r8 - e*e*uu*uu))
   gamma(w,uu,g,e) = .50_r8*w*((3._r8*g*(1._r8-w)*uu*uu + 1._r8)/(1._r8-e*e*uu*uu))
   el(w,g)         = sqrt(3._r8*(1._r8-w)*(1._r8 - w*g))
   taus(w,f,t)     = (1._r8 - w*f)*t
   omgs(w,f)       = (1._r8 - f)*w/(1._r8 - w*f)
   asys(g,f)       = (g - f)/(1._r8 - f)
   u(w,g,e)        = 1.5_r8*(1._r8 - w*g)/e
   n(uu,et)        = ((uu+1._r8)*(uu+1._r8)/et ) - ((uu-1._r8)*(uu-1._r8)*et)
!
!-----------------------------------------------------------------------
!
! Compute layer radiative properties
!
! Compute radiative properties (reflectivity and transmissivity for
!    direct and diffuse radiation incident from above, under clear
!    and cloudy conditions) and transmission of direct radiation
!    (under clear and cloudy conditions) for each layer.
!
   do k=0,pver
      do i=1,ndayc
         if(aer_tau(i,k) > 0._r8) then !where(aer_tau > 0._r8)
            aer_g_limit(i,k) = aer_tau_w_g(i,k) / aer_tau_w(i,k)
            aer_f_limit(i,k) = aer_tau_w_f(i,k) / aer_tau_w(i,k)
            aer_tau_w(i,k) = aer_tau(i,k) * min(aer_tau_w(i,k)/aer_tau(i,k) , 0.999999_r8)
         else
            aer_tau(i,k) = 0._r8
            aer_tau_w(i,k) = 0._r8
            aer_g_limit(i,k) = 0._r8
            aer_f_limit(i,k) = 0._r8
         endif
         aer_tau_w_g(i,k) = aer_tau_w(i,k) * aer_g_limit(i,k)
         aer_tau_w_f(i,k) = aer_tau_w(i,k) * aer_f_limit(i,k)
      enddo
   enddo

   do k=0,pver
      do i=1,ndayc
            tauray(i) = trayoslp*(pflx(i,k+1)-pflx(i,k))
            taugab(i) = abh2o*uh2o(i,k) + abo3*uo3(i,k) + abco2*uco2(i,k) + abo2*uo2(i,k)
            tautot = tauxcl(i,k) + tauxci(i,k) + tauray(i) + taugab(i) + aer_tau(i,k)
            taucsc = tauxcl(i,k)*wcl(i,k) + tauxci(i,k)*wci(i,k) + aer_tau_w(i,k)
            wtau   = wray*tauray(i)
            wt     = wtau + taucsc
            wtot   = wt/tautot
            gtot   = (wtau*gray + gcl(i,k)*wcl(i,k)*tauxcl(i,k) &
                     + gci(i,k)*wci(i,k)*tauxci(i,k) + aer_tau_w_g(i,k))/wt
            ftot   = (wtau*fray + fcl(i,k)*wcl(i,k)*tauxcl(i,k) &
                     + fci(i,k)*wci(i,k)*tauxci(i,k) + aer_tau_w_f(i,k))/wt
            ts   = taus(wtot,ftot,tautot)
            ws   = omgs(wtot,ftot)
            gs   = asys(gtot,ftot)
            lm   = el(ws,gs)
            alp  = alpha(ws,coszrs(i),gs,lm)
            gam  = gamma(ws,coszrs(i),gs,lm)
            ue   = u(ws,gs,lm)
!
!     Limit argument of exponential to 25, in case lm very large:
!
            arg  = min(lm*ts,25._r8)
            extins = exp(-arg)
            ne = n(ue,extins)
            rdif(ns,i,k) = (ue+1._r8)*(ue-1._r8)*(1._r8/extins - extins)/ne
            tdif(ns,i,k)   =   4._r8*ue/ne
!
!     Limit argument of exponential to 25, in case coszrs is very small:
!
            arg       = min(ts/coszrs(i),25._r8)
            explay(ns,i,k) = exp(-arg)
            apg = alp + gam
            amg = alp - gam
            rdir(ns,i,k) = amg*(tdif(ns,i,k)*explay(ns,i,k)-1._r8) + apg*rdif(ns,i,k)
            tdir(ns,i,k) = apg*tdif(ns,i,k) + (amg*rdif(ns,i,k)-(apg-1._r8))*explay(ns,i,k)
!
!     Under rare conditions, reflectivies and transmissivities can be
!     negative; zero out any negative values
!
            rdir(ns,i,k) = max(rdir(ns,i,k),0.0_r8)
            tdir(ns,i,k) = max(tdir(ns,i,k),0.0_r8)
            rdif(ns,i,k) = max(rdif(ns,i,k),0.0_r8)
            tdif(ns,i,k) = max(tdif(ns,i,k),0.0_r8)
!
!     Clear-sky calculation
!
            if (tauxcl(i,k) == 0.0_r8 .and. tauxci(i,k) == 0.0_r8) then

               rdirc(ns,i,k) = rdir(ns,i,k)
               tdirc(ns,i,k) = tdir(ns,i,k)
               rdifc(ns,i,k) = rdif(ns,i,k)
               tdifc(ns,i,k) = tdif(ns,i,k)
               explayc(ns,i,k) = explay(ns,i,k)
            else
               tautot = tauray(i) + taugab(i) + aer_tau(i,k)
               taucsc = aer_tau_w(i,k)
!
! wtau already computed for all-sky
!
               wt     = wtau + taucsc
               wtot   = wt/tautot
               gtot   = (wtau*gray + aer_tau_w_g(i,k))/wt
               ftot   = (wtau*fray + aer_tau_w_f(i,k))/wt
               ts   = taus(wtot,ftot,tautot)
               ws   = omgs(wtot,ftot)
               gs   = asys(gtot,ftot)
               lm   = el(ws,gs)
               alp  = alpha(ws,coszrs(i),gs,lm)
               gam  = gamma(ws,coszrs(i),gs,lm)
               ue   = u(ws,gs,lm)
!
!     Limit argument of exponential to 25, in case lm very large:
!
               arg  = min(lm*ts,25._r8)
               extins = exp(-arg)
               ne = n(ue,extins)
               rdifc(ns,i,k) = (ue+1._r8)*(ue-1._r8)*(1._r8/extins - extins)/ne
               tdifc(ns,i,k)   =   4._r8*ue/ne
!
!     Limit argument of exponential to 25, in case coszrs is very small:
!
               arg       = min(ts/coszrs(i),25._r8)
               explayc(ns,i,k) = exp(-arg)
               apg = alp + gam
               amg = alp - gam
               rdirc(ns,i,k) = amg*(tdifc(ns,i,k)*explayc(ns,i,k)-1._r8)+ &
                               apg*rdifc(ns,i,k)
               tdirc(ns,i,k) = apg*tdifc(ns,i,k) + (amg*rdifc(ns,i,k) - (apg-1._r8))* &
                               explayc(ns,i,k)
!
!     Under rare conditions, reflectivies and transmissivities can be
!     negative; zero out any negative values
!
               rdirc(ns,i,k) = max(rdirc(ns,i,k),0.0_r8)
               tdirc(ns,i,k) = max(tdirc(ns,i,k),0.0_r8)
               rdifc(ns,i,k) = max(rdifc(ns,i,k),0.0_r8)
               tdifc(ns,i,k) = max(tdifc(ns,i,k),0.0_r8)
            end if
         end do
   end do

end subroutine raddedmx

!-------------------------------------------------------------------------------

subroutine radsw_init(gravx)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize various constants for radiation scheme; note that
! the radiation scheme uses cgs units.
! 
! Author: W. Collins (H2O parameterization) and J. Kiehl
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: gravx      ! Acceleration of gravity (MKS)

   real(r8), parameter :: ref_tsi = 1367._r8 ! value supplied by Dan Marsh -- see solvar_woods.F90
!
!-----------------------------------------------------------------------
!
! Set general radiation consts; convert to cgs units where appropriate:
!
   gravit  =  100._r8*gravx
   rga     =  1._r8/gravit
   sslp    =  1.013250e6_r8

end subroutine radsw_init

!-------------------------------------------------------------------------------

end module radsw
