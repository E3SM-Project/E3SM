module micro_p3_utils

  use physical_constants,     only: pi => dd_pi
  
  implicit none
  private
  save

    public :: get_latent_heat, micro_p3_utils_init, &
              avg_diameter, calculate_incloud_mixingratios

    integer, public :: iulog_e3sm
    logical, public :: masterproc_e3sm

    ! Signaling NaN bit pattern that represents a limiter that's turned off.
    !integer(itype), parameter :: limiter_off = int(Z'7FF1111111111111', itype)
    integer, parameter :: limiter_off = 0

    real, public, parameter :: qsmall = 1.e-14
    real, public, parameter :: nsmall = 1.e-16

    real :: latent_heat_vapor, latent_heat_sublim, latent_heat_fusion

    real,public :: rho_1000mb,rho_600mb,rho_h2o,  &
       cpw,cons1,cons2,cons3,cons4,cons5,cons6,cons7,    &
       inv_rho_h2o,inv_dropmass,cp,g,rd,rv,ep_2,inv_cp

    real, public, parameter :: thrd  = 1./3.
    real, public, parameter :: sxth  = 1./6.
    real, public, parameter :: piov3 = pi*thrd
    real, public, parameter :: piov6 = pi*sxth

    ! maximum total ice concentration (sum of all categories)
    real, public, parameter :: max_total_ni = 500.e+3  ! (m)

    ! parameters for Seifert and Beheng (2001) autoconversion/accretion
    real, public, parameter :: kc     = 9.44e+9
    real, public, parameter :: kr     = 5.78e+3

    real, public, parameter :: ar     = 841.99667 
    real, public, parameter :: br     = 0.8
    real, public, parameter :: f1r    = 0.78
    real, public, parameter :: f2r    = 0.32
    real, public, parameter :: ecr    = 1.

    ! limits for rime density [kg m-3]
    real, public, parameter :: rho_rimeMin     =  50.
    real, public, parameter :: rho_rimeMax     = 900.
    real, public, parameter :: inv_rho_rimeMax =   1./rho_rimeMax

    ! Barklie and Gokhale (1959)
    real, public, parameter :: bimm   = 2.
    real, public, parameter :: aimm   = 0.65
    real, public, parameter :: rin    = 0.1e-6
    real, public, parameter :: mi0    = 4.*piov3*900.*1.e-18

    real, public, parameter :: eci    = 0.5
    real, public, parameter :: eri    = 1.
    real, public, parameter :: bcn    = 2.

    ! mean size for soft lambda_r limiter [microns]
    real, public, parameter :: dbrk   = 600.e-6
    ! ratio of rain number produced to ice number loss from melting
    real, public, parameter :: nmltratio = 1.0
    
    ! calibration factors for ice deposition and sublimation
    !   These are adjustable ad hoc factors used to increase or decrease deposition and/or
    !   sublimation rates.  The representation of the ice capacitances are highly simplified
    !   and the appropriate values in the diffusional growth equation are uncertain.
    real, public, parameter :: clbfact_dep = 1.
    real, public, parameter :: clbfact_sub = 1.
    
    logical,public  :: do_Cooper_inP3   ! Use prescribed CCN       

    real,dimension(16), public :: dnu

    real, public, parameter :: mu_r_constant = 0.0 !1.0
    real, public, parameter :: lookup_table_1a_dum1_c =  4.135985029041767d+00 ! 1.0/(0.1*log10(261.7))

    real,public :: T_zerodegc  ! Temperature at zero degree celcius ~K
    real,public :: T_rainfrz  ! Contact and immersion freexing temp, -4C  ~K
    real,public :: T_homogfrz ! Homogeneous freezing temperature, -40C  ~K
    real,public :: T_icenuc   ! Ice nucleation temperature, -5C ~K

    real,public :: pi_e3sm
    ! ice microphysics lookup table array dimensions
    integer, public,parameter :: isize        = 50
    integer, public,parameter :: densize      =  5
    integer, public,parameter :: rimsize      =  4
    integer, public,parameter :: rcollsize    = 30
    integer, public,parameter :: ice_table_size      = 12  ! number of quantities used from lookup table
    integer, public,parameter :: collect_table_size  =  2  ! number of ice-rain collection  quantities used from lookup table
    ! switch for warm-rain parameterization
    ! = 1 Seifert and Beheng 2001
    ! = 2 Beheng 1994
    ! = 3 Khairoutdinov and Kogan 2000
    integer, public,parameter :: iparam = 3

    real, parameter, public :: mincld=0.0001
    real, parameter, public :: rho_h2os = 917.  ! bulk density water solid
    real, parameter, public :: dropmass = 5.2e-7

    ! particle mass-diameter relationship
    ! currently we assume spherical particles for cloud ice/snow
    ! m = cD^d
    ! exponent
    real, parameter :: dsph = 3.

    ! Bounds for mean diameter for different constituents.
    ! real, parameter :: lam_bnd_rain(2) = 1./[500.e-6, 20.e-6]
    ! real, parameter :: lam_bnd_snow(2) = 1./[2000.e-6, 10.e-6]

    ! Minimum average mass of particles.
    real, parameter :: min_mean_mass_liq = 1.e-20
    real, parameter :: min_mean_mass_ice = 1.e-20

    ! in-cloud values
    REAL, PARAMETER :: min_cld_frac   = 1.e-20 !! threshold min value for cloud fraction
    real, parameter :: incloud_limit = 5.1E-3
    real, parameter :: precip_limit  = 1.0E-2

    contains
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine micro_p3_utils_init(cpair,rair,rh2o,rhoh2o,mwh2o,mwdry,gravit,latvap,latice, &
                   cpliq,tmelt,pi,iulog,masterproc)

    real, intent(in) :: cpair
    real, intent(in) :: rair
    real, intent(in) :: rh2o
    real, intent(in) :: rhoh2o
    real, intent(in) :: mwh2o
    real, intent(in) :: mwdry
    real, intent(in) :: gravit
    real, intent(in) :: latvap
    real, intent(in) :: latice
    real, intent(in) :: cpliq
    real, intent(in) :: tmelt
    real, intent(in) :: pi
    integer, intent(in)     :: iulog
    logical, intent(in)     :: masterproc

    ! logfile info
    iulog_e3sm      = iulog
    masterproc_e3sm = masterproc

    ! mathematical/optimization constants
     
    pi_e3sm = pi

    ! Temperature parameters
    T_zerodegc  = tmelt 
    T_homogfrz = tmelt-40.
    T_icenuc   = tmelt-15.
    T_rainfrz  = tmelt-4.

    ! physical constants
    cp     = cpair ! specific heat of dry air (J/K/kg) !1005.
    inv_cp = 1./cp ! inverse of cp
    g      = gravit ! Gravity (m/s^2) !9.816
    rd     = rair ! Dry air gas constant     ~ J/K/kg
    rv     = rh2o ! Water vapor gas constant ~ J/K/kg     !461.51
    ep_2   = mwh2o/mwdry  ! ratio of molecular mass of water to the molecular mass of dry air !0.622
    rho_1000mb = 100000./(rd*T_zerodegc) ! density of air at surface
    rho_600mb = 60000./(rd*253.15)
    rho_h2o   = rhoh2o ! Density of liquid water (STP) !997.
    cpw    = cpliq  ! specific heat of fresh h2o (J/K/kg) !4218.
    inv_rho_h2o = 1./rho_h2o  !inverse of (max.) density of liquid water
    inv_dropmass = 1./dropmass  !inverse of dropmass

    latent_heat_vapor = latvap           ! latent heat of vaporization
    latent_heat_sublim = latvap + latice  ! latent heat of sublimation
    latent_heat_fusion  = latice           ! latent heat of fusion

    ! Bigg (1953)
    !bimm   = 100.
    !aimm   = 0.66

    cons1 = piov6*rho_h2o
    cons2 = 4.*piov3*rho_h2o
    cons3 = 1./(cons2)  ! moved embryonic_rain_size out of cons3 into cloud_water_autoconversion
    cons4 = 1./(dbrk**3*pi*rho_h2o)
    cons5 = piov6*bimm
    cons6 = piov6**2*rho_h2o*bimm
    cons7 = 4.*piov3*rho_h2o*1.e-18

    ! droplet spectral shape parameter for mass spectra, used for Seifert and Beheng (2001)
    ! warm rain autoconversion/accretion option only (iparam = 1)
!    allocate(dnu(16))
    dnu(1)  =  0.000
    dnu(2)  = -0.557
    dnu(3)  = -0.430
    dnu(4)  = -0.307
    dnu(5)  = -0.186
    dnu(6)  = -0.067
    dnu(7)  = -0.050
    dnu(8)  = -0.167
    dnu(9)  = -0.282
    dnu(10) = -0.397
    dnu(11) = -0.512
    dnu(12) = -0.626
    dnu(13) = -0.739
    dnu(14) = -0.853
    dnu(15) = -0.966
    dnu(16) = -0.966

    
    return
    end subroutine micro_p3_utils_init
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine get_latent_heat(kts,kte,v,s,f)

       integer,intent(in) :: kts,kte
       real,dimension(kts:kte),intent(out) :: v,s,f


       v(:) = latent_heat_vapor !latvap           ! latent heat of vaporization
       s(:) = latent_heat_sublim !latvap + latice  ! latent heat of sublimation
       f(:) = latent_heat_fusion  !latice           ! latent heat of fusion
 
! Original P3 definition of latent heats:   
!       do i = its,ite
!          do k = kts,kte
!          latent_heat_vapor(i,k)    = 3.1484e6-2370.*t(i,k)
!          latent_heat_sublim(i,k)    = latent_heat_vapor(i,k)+0.3337e6
!          latent_heat_fusion(i,k)     = latent_heat_sublim(i,k)-latent_heat_vapor(i,k)
!          end do
!       end do
       return
    end subroutine get_latent_heat

!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine calculate_incloud_mixingratios(qc,qr,qi,qm,nc,nr,ni,bm, &
          inv_cld_frac_l,inv_cld_frac_i,inv_cld_frac_r, &
          qc_incld,qr_incld,qi_incld,qm_incld,nc_incld,nr_incld,ni_incld,bm_incld)

       real,intent(in)   :: qc, qr, qi, qm
       real,intent(in)   :: nc, nr, ni, bm
       real,intent(in)   :: inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r
       real,intent(out)  :: qc_incld, qr_incld, qi_incld, qm_incld
       real,intent(out)  :: nc_incld, nr_incld, ni_incld, bm_incld

       if (qc.ge.qsmall) then
          qc_incld = qc*inv_cld_frac_l
          nc_incld = max(nc*inv_cld_frac_l,0.)
       else
          qc_incld = 0.
          nc_incld = 0.
       end if 
       if (qi.ge.qsmall) then
          qi_incld = qi*inv_cld_frac_i
          ni_incld = max(ni*inv_cld_frac_i,0.)
       else
          qi_incld = 0.
          ni_incld = 0.
       end if 
       if (qm.ge.qsmall.and.qi.ge.qsmall) then
          qm_incld = qm*inv_cld_frac_i
          bm_incld = max(bm*inv_cld_frac_l,0.)
       else
          qm_incld = 0.
          bm_incld = 0.
       end if 
       if (qr.ge.qsmall) then
          qr_incld = qr*inv_cld_frac_r
          nr_incld = max(nr*inv_cld_frac_r,0.)
       else
          qr_incld = 0.
          nr_incld = 0.
       end if
       if (qc_incld.gt.incloud_limit .or.qi_incld.gt.incloud_limit &
            .or. qr_incld.gt.precip_limit .or.bm_incld.gt.incloud_limit) then
          qc_incld    = min(qc_incld,incloud_limit)
          qi_incld = min(qi_incld,incloud_limit)
          bm_incld = min(bm_incld,incloud_limit)
          qr_incld    = min(qr_incld,precip_limit)
       end if
    end subroutine calculate_incloud_mixingratios
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!

real elemental function avg_diameter(q, n, rho_air, rho_sub)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  ! Assumes that diameter follows an exponential distribution.
  real, intent(in) :: q         ! mass mixing ratio
  real, intent(in) :: n         ! number concentration (per volume)
  real, intent(in) :: rho_air   ! local density of the air
  real, intent(in) :: rho_sub   ! density of the particle substance

  avg_diameter = (pi_e3sm * rho_sub * n/(q*rho_air))**(-1./3.)

end function avg_diameter


end module micro_p3_utils
