module micro_params

  use grid, only: nzm
  use params, only: crm_rknd

  implicit none

  !  Microphysics stuff:

  ! Densities of hydrometeors

  real(crm_rknd), parameter :: rhor = 1000. ! Density of water, kg/m3
  real(crm_rknd), parameter :: rhos = 100.  ! Density of snow, kg/m3
  real(crm_rknd), parameter :: rhog = 400.  ! Density of graupel, kg/m3
  !real(crm_rknd), parameter :: rhog = 917.  ! hail - Lin 1983

  ! Temperatures limits for various hydrometeors

  real(crm_rknd), parameter :: tbgmin = 253.16    ! Minimum temperature for cloud water., K
  real(crm_rknd), parameter :: tbgmax = 273.16    ! Maximum temperature for cloud ice, K
  real(crm_rknd), parameter :: tprmin = 268.16    ! Minimum temperature for rain, K
  real(crm_rknd), parameter :: tprmax = 283.16    ! Maximum temperature for snow+graupel, K
  real(crm_rknd), parameter :: tgrmin = 223.16    ! Minimum temperature for snow, K
  real(crm_rknd), parameter :: tgrmax = 283.16    ! Maximum temperature for graupel, K

  ! Terminal velocity coefficients

  real(crm_rknd), parameter :: a_rain = 842. ! Coeff.for rain term vel
  real(crm_rknd), parameter :: b_rain = 0.8  ! Fall speed exponent for rain
  real(crm_rknd), parameter :: a_snow = 4.84 ! Coeff.for snow term vel
  real(crm_rknd), parameter :: b_snow = 0.25 ! Fall speed exponent for snow
  !real(crm_rknd), parameter :: a_grau = 40.7! Krueger (1994) ! Coef. for graupel term vel
  real(crm_rknd), parameter :: a_grau = 94.5 ! Lin (1983) (rhog=400)
  !real(crm_rknd), parameter :: a_grau = 127.94! Lin (1983) (rhog=917)
  real(crm_rknd), parameter :: b_grau = 0.5  ! Fall speed exponent for graupel

  ! Autoconversion
#ifdef CLUBB_CRM   /*microphysical tuning for CLUBB*/
  real(crm_rknd), parameter :: qcw0 = 0.6e-3      ! Threshold for water autoconversion, g/g
  real(crm_rknd), parameter :: qci0 = 1.e-4     ! Threshold for ice autoconversion, g/g
  real(crm_rknd), parameter :: alphaelq = 10.e-3  ! autoconversion of cloud water rate coef
  real(crm_rknd), parameter :: betaelq = 6.0e-3   ! autoconversion of cloud ice rate coef
#else
  real(crm_rknd), parameter :: qcw0 = 1.e-3      ! Threshold for water autoconversion, g/g
  real(crm_rknd), parameter :: qci0 = 1.e-4     ! Threshold for ice autoconversion, g/g
  real(crm_rknd), parameter :: alphaelq = 1.e-3  ! autoconversion of cloud water rate coef
  real(crm_rknd), parameter :: betaelq = 1.e-3   ! autoconversion of cloud ice rate coef
#endif /*CLUBB_CRM*/

  ! Accretion

  real(crm_rknd), parameter :: erccoef = 1.0   ! Rain/Cloud water collection efficiency
  real(crm_rknd), parameter :: esccoef = 1.0   ! Snow/Cloud water collection efficiency
  real(crm_rknd), parameter :: esicoef = 0.1   ! Snow/cloud ice collection efficiency
  real(crm_rknd), parameter :: egccoef = 1.0   ! Graupel/Cloud water collection efficiency
  real(crm_rknd), parameter :: egicoef = 0.1   ! Graupel/Cloud ice collection efficiency

  ! Interseption parameters for exponential size spectra

  real(crm_rknd), parameter :: nzeror = 8.e6   ! Intercept coeff. for rain
  real(crm_rknd), parameter :: nzeros = 3.e6   ! Intersept coeff. for snow
  real(crm_rknd), parameter :: nzerog = 4.e6   ! Intersept coeff. for graupel
  !real(crm_rknd), parameter :: nzerog = 4.e4   ! hail - Lin 1993

  real(crm_rknd), parameter :: qp_threshold = 1.e-8 ! minimal rain/snow water content


  ! Misc. microphysics variables

  real*4 gam3       ! Gamma function of 3
  real*4 gams1      ! Gamma function of (3 + b_snow)
  real*4 gams2      ! Gamma function of (5 + b_snow)/2
  real*4 gams3      ! Gamma function of (4 + b_snow)
  real*4 gamg1      ! Gamma function of (3 + b_grau)
  real*4 gamg2      ! Gamma function of (5 + b_grau)/2
  real*4 gamg3      ! Gamma function of (4 + b_grau)
  real*4 gamr1      ! Gamma function of (3 + b_rain)
  real*4 gamr2      ! Gamma function of (5 + b_rain)/2
  real*4 gamr3      ! Gamma function of (4 + b_rain)

  real(crm_rknd), allocatable :: accrsc (:,:)
  real(crm_rknd), allocatable :: accrsi (:,:)
  real(crm_rknd), allocatable :: accrrc (:,:)
  real(crm_rknd), allocatable :: coefice(:,:)
  real(crm_rknd), allocatable :: accrgc (:,:)
  real(crm_rknd), allocatable :: accrgi (:,:)
  real(crm_rknd), allocatable :: evaps1 (:,:)
  real(crm_rknd), allocatable :: evaps2 (:,:)
  real(crm_rknd), allocatable :: evapr1 (:,:)
  real(crm_rknd), allocatable :: evapr2 (:,:)
  real(crm_rknd), allocatable :: evapg1 (:,:)
  real(crm_rknd), allocatable :: evapg2 (:,:)

  real(crm_rknd) a_bg, a_pr, a_gr


contains


  subroutine allocate_micro_params(ncrms)
#if defined(_OPENACC)
    use openacc_utils
#endif
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: zero

    allocate( accrsc(ncrms,nzm) )
    allocate( accrsi(ncrms,nzm) )
    allocate( accrrc(ncrms,nzm) )
    allocate( coefice(ncrms,nzm) )
    allocate( accrgc(ncrms,nzm) )
    allocate( accrgi(ncrms,nzm) )
    allocate( evaps1(ncrms,nzm) )
    allocate( evaps2(ncrms,nzm) )
    allocate( evapr1(ncrms,nzm) )
    allocate( evapr2(ncrms,nzm) )
    allocate( evapg1(ncrms,nzm) )
    allocate( evapg2(ncrms,nzm) )
#if defined(_OPENACC)
    call prefetch( accrsc  )
    call prefetch( accrsi  )
    call prefetch( accrrc  )
    call prefetch( coefice )
    call prefetch( accrgc  )
    call prefetch( accrgi  )
    call prefetch( evaps1  )
    call prefetch( evaps2  )
    call prefetch( evapr1  )
    call prefetch( evapr2  )
    call prefetch( evapg1  )
    call prefetch( evapg2  )
#elif defined(_OPENMP)
    !$omp target enter data map(alloc: accrsc  )
    !$omp target enter data map(alloc: accrsi  )
    !$omp target enter data map(alloc: accrrc  )
    !$omp target enter data map(alloc: coefice )
    !$omp target enter data map(alloc: accrgc  )
    !$omp target enter data map(alloc: accrgi  )
    !$omp target enter data map(alloc: evaps1  )
    !$omp target enter data map(alloc: evaps2  )
    !$omp target enter data map(alloc: evapr1  )
    !$omp target enter data map(alloc: evapr2  )
    !$omp target enter data map(alloc: evapg1  )
    !$omp target enter data map(alloc: evapg2  )
#endif
    zero = 0

    accrsc  = zero
    accrsi  = zero
    accrrc  = zero
    coefice = zero
    accrgc  = zero
    accrgi  = zero
    evaps1  = zero
    evaps2  = zero
    evapr1  = zero
    evapr2  = zero
    evapg1  = zero
    evapg2  = zero
  end subroutine allocate_micro_params

#if defined(_OPENMP)
  subroutine update_device_micro_params()
    implicit none
    !$omp target update to( accrsc  )
    !$omp target update to( accrsi  )
    !$omp target update to( accrrc  )
    !$omp target update to( coefice )
    !$omp target update to( accrgc  )
    !$omp target update to( accrgi  )
    !$omp target update to( evaps1  )
    !$omp target update to( evaps2  )
    !$omp target update to( evapr1  )
    !$omp target update to( evapr2  )
    !$omp target update to( evapg1  )
    !$omp target update to( evapg2  )
  end subroutine update_device_micro_params

  subroutine update_host_micro_params()
    implicit none
    !$omp target update from( accrsc  )
    !$omp target update from( accrsi  )
    !$omp target update from( accrrc  )
    !$omp target update from( coefice )
    !$omp target update from( accrgc  )
    !$omp target update from( accrgi  )
    !$omp target update from( evaps1  )
    !$omp target update from( evaps2  )
    !$omp target update from( evapr1  )
    !$omp target update from( evapr2  )
    !$omp target update from( evapg1  )
    !$omp target update from( evapg2  )
  end subroutine update_host_micro_params
#endif

  subroutine deallocate_micro_params()
    implicit none
#if defined(_OPENMP)
    !$omp target exit data map(delete: accrsc  )
    !$omp target exit data map(delete: accrsi  )
    !$omp target exit data map(delete: accrrc  )
    !$omp target exit data map(delete: coefice )
    !$omp target exit data map(delete: accrgc  )
    !$omp target exit data map(delete: accrgi  )
    !$omp target exit data map(delete: evaps1  )
    !$omp target exit data map(delete: evaps2  )
    !$omp target exit data map(delete: evapr1  )
    !$omp target exit data map(delete: evapr2  )
    !$omp target exit data map(delete: evapg1  )
    !$omp target exit data map(delete: evapg2  )
#endif
    deallocate( accrsc  )
    deallocate( accrsi  )
    deallocate( accrrc  )
    deallocate( coefice )
    deallocate( accrgc  )
    deallocate( accrgi  )
    deallocate( evaps1  )
    deallocate( evaps2  )
    deallocate( evapr1  )
    deallocate( evapr2  )
    deallocate( evapg1  )
    deallocate( evapg2  )
  end subroutine deallocate_micro_params

end module micro_params
