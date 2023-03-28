module micro_p3_utils

  use physics_utils, only: rtype, rtype8, itype, btype

  implicit none
  private
  save

    public :: get_latent_heat, micro_p3_utils_init, &
              avg_diameter, calculate_incloud_mixingratios

    integer, public :: iulog_e3sm
    logical(btype), public :: masterproc_e3sm

    ! Signaling NaN bit pattern that represents a limiter that's turned off.
    integer(itype), parameter :: limiter_off = int(Z'7FF1111111111111', itype)

    real(rtype), public, parameter :: qsmall = 1.e-14_rtype
    real(rtype), public, parameter :: nsmall = 1.e-16_rtype

    real(rtype) :: latent_heat_vapor, latent_heat_sublim, latent_heat_fusion

    real(rtype),public :: rho_1000mb,rho_600mb,ar,br,f1r,f2r,ecr,rho_h2o,kr,kc,aimm,bimm,rin,mi0,nccnst,  &
       eci,eri,bcn,cpw,cons1,cons2,cons3,cons4,cons5,cons6,cons7,         &
       inv_rho_h2o,inv_dropmass,cp,g,rd,rv,ep_2,inv_cp,   &
       thrd,sxth,piov3,piov6,rho_rimeMin,     &
       rho_rimeMax,inv_rho_rimeMax,max_total_ni,dbrk,nmltratio,clbfact_sub,  &
       clbfact_dep

    logical,public  :: do_Cooper_inP3   ! Use prescribed CCN       

    real(rtype),dimension(16), public :: dnu

    real(rtype), public, parameter :: mu_r_constant = 0.0_rtype !1.0_rtype
    real(rtype), public, parameter :: lookup_table_1a_dum1_c =  4.135985029041767d+00 ! 1.0/(0.1*log10(261.7))

    real(rtype),public :: T_zerodegc  ! Temperature at zero degree celcius ~K
    real(rtype),public :: T_rainfrz  ! Contact and immersion freexing temp, -4C  ~K
    real(rtype),public :: T_homogfrz ! Homogeneous freezing temperature, -40C  ~K
    real(rtype),public :: T_icenuc   ! Ice nucleation temperature, -5C ~K

    real(rtype),public :: pi_e3sm
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

    real(rtype), parameter, public :: mincld=0.0001_rtype
    real(rtype), parameter, public :: rho_h2os = 917._rtype  ! bulk density water solid
    real(rtype), parameter, public :: dropmass = 5.2e-7_rtype

    ! particle mass-diameter relationship
    ! currently we assume spherical particles for cloud ice/snow
    ! m = cD^d
    ! exponent
    real(rtype), parameter :: dsph = 3._rtype

    ! Bounds for mean diameter for different constituents.
    ! real(rtype), parameter :: lam_bnd_rain(2) = 1._rtype/[500.e-6_rtype, 20.e-6_rtype]
    ! real(rtype), parameter :: lam_bnd_snow(2) = 1._rtype/[2000.e-6_rtype, 10.e-6_rtype]

    ! Minimum average mass of particles.
    real(rtype), parameter :: min_mean_mass_liq = 1.e-20_rtype
    real(rtype), parameter :: min_mean_mass_ice = 1.e-20_rtype

    ! in-cloud values
    REAL(rtype), PARAMETER :: min_cld_frac   = 1.e-20_rtype !! threshold min value for cloud fraction
    real(rtype), parameter :: incloud_limit = 5.1E-3
    real(rtype), parameter :: precip_limit  = 1.0E-2

    contains
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine micro_p3_utils_init(cpair,rair,rh2o,rhoh2o,mwh2o,mwdry,gravit,latvap,latice, &
                   cpliq,tmelt,pi,iulog,masterproc)

    real(rtype), intent(in) :: cpair
    real(rtype), intent(in) :: rair
    real(rtype), intent(in) :: rh2o
    real(rtype), intent(in) :: rhoh2o
    real(rtype), intent(in) :: mwh2o
    real(rtype), intent(in) :: mwdry
    real(rtype), intent(in) :: gravit
    real(rtype), intent(in) :: latvap
    real(rtype), intent(in) :: latice
    real(rtype), intent(in) :: cpliq
    real(rtype), intent(in) :: tmelt
    real(rtype), intent(in) :: pi
    integer, intent(in)     :: iulog
    logical(btype), intent(in)     :: masterproc


    return
    end subroutine micro_p3_utils_init
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine get_latent_heat(its,ite,kts,kte,v,s,f)

       integer,intent(in) :: its,ite,kts,kte
       real(rtype),dimension(its:ite,kts:kte),intent(out) :: v,s,f

!       integer i,k

    end subroutine get_latent_heat

!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine calculate_incloud_mixingratios(qc,qr,qi,qm,nc,nr,ni,bm, &
          inv_cld_frac_l,inv_cld_frac_i,inv_cld_frac_r, &
          qc_incld,qr_incld,qi_incld,qm_incld,nc_incld,nr_incld,ni_incld,bm_incld)

       real(rtype),intent(in)   :: qc, qr, qi, qm
       real(rtype),intent(in)   :: nc, nr, ni, bm
       real(rtype),intent(in)   :: inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r
       real(rtype),intent(out)  :: qc_incld, qr_incld, qi_incld, qm_incld
       real(rtype),intent(out)  :: nc_incld, nr_incld, ni_incld, bm_incld

    end subroutine calculate_incloud_mixingratios
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!

real(rtype) elemental function avg_diameter(q, n, rho_air, rho_sub)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  ! Assumes that diameter follows an exponential distribution.
  real(rtype), intent(in) :: q         ! mass mixing ratio
  real(rtype), intent(in) :: n         ! number concentration (per volume)
  real(rtype), intent(in) :: rho_air   ! local density of the air
  real(rtype), intent(in) :: rho_sub   ! density of the particle substance


end function avg_diameter


end module micro_p3_utils
