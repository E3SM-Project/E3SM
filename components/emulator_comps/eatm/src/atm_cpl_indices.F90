module atm_cpl_indices

  use seq_flds_mod
  use mct_mod

  implicit none

  SAVE
  public                               ! By default make data private

  !JW list of all fields from EAM except removed BGC
  integer :: index_a2x_Sa_z       = 0  ! bottom atm level height
  integer :: index_a2x_Sa_u       = 0  ! bottom atm level zon wind
  integer :: index_a2x_Sa_v       = 0  ! bottom atm level mer wind
  integer :: index_a2x_Sa_wsresp  = 0  ! first order response of wind to stress
  integer :: index_a2x_Sa_tau_est = 0  ! estimated stress in equilibrium with ubot/vbot
  integer :: index_a2x_Sa_ugust   = 0  ! magnitude of gustiness at surface
  integer :: index_a2x_Sa_topo    = 0  ! topography
  integer :: index_a2x_Sa_tbot    = 0  ! bottom atm level temp
  integer :: index_a2x_Sa_ptem    = 0  ! bottom atm level pot temp
  integer :: index_a2x_Sa_shum    = 0  ! bottom atm level spec hum
  integer :: index_a2x_Sa_dens    = 0  ! bottom atm level air den
  integer :: index_a2x_Sa_pbot    = 0  ! bottom atm level pressure
  integer :: index_a2x_Sa_pslv    = 0  ! sea level atm pressure
  integer :: index_a2x_Sa_uovern  = 0  ! ratio of wind speed/brunt vaisalla frequency
  integer :: index_a2x_Faxa_lwdn  = 0  ! downward lw heat flux
  integer :: index_a2x_Faxa_rainc = 0  ! prec: liquid "convective"
  integer :: index_a2x_Faxa_rainl = 0  ! prec: liquid "large scale"
  integer :: index_a2x_Faxa_snowc = 0  ! prec: frozen "convective"
  integer :: index_a2x_Faxa_snowl = 0  ! prec: frozen "large scale"
  integer :: index_a2x_Faxa_swndr = 0  ! sw: nir direct  downward
  integer :: index_a2x_Faxa_swvdr = 0  ! sw: vis direct  downward
  integer :: index_a2x_Faxa_swndf = 0  ! sw: nir diffuse downward
  integer :: index_a2x_Faxa_swvdf = 0  ! sw: vis diffuse downward
  integer :: index_a2x_Faxa_swnet = 0  ! sw: net

  integer :: index_x2a_Sx_t       = 0  ! surface temperature
  integer :: index_x2a_So_t       = 0  ! sea surface temperature
  integer :: index_x2a_Sf_lfrac   = 0  ! surface land fraction
  integer :: index_x2a_Sf_ifrac   = 0  ! surface ice fraction
  integer :: index_x2a_Sf_ofrac   = 0  ! surface ocn fraction
  integer :: index_x2a_Sx_tref    = 0  ! 2m reference temperature
  integer :: index_x2a_Sx_qref    = 0  ! 2m reference specific humidity
  integer :: index_x2a_Sx_avsdr   = 0  ! albedo, visible, direct
  integer :: index_x2a_Sx_anidr   = 0  ! albedo, near-ir, direct
  integer :: index_x2a_Sx_avsdf   = 0  ! albedo, visible, diffuse
  integer :: index_x2a_Sx_anidf   = 0  ! albedo, near-ir, diffuse
  integer :: index_x2a_Sl_snowh   = 0  ! surface snow depth over land
  integer :: index_x2a_Si_snowh   = 0  ! surface snow depth over ice
  integer :: index_x2a_Sl_fv      = 0  ! friction velocity
  integer :: index_x2a_Sl_ram1    = 0  ! aerodynamical resistance
  integer :: index_x2a_Sl_soilw   = 0  ! volumetric soil water
  integer :: index_x2a_Faxx_taux  = 0  ! wind stress, zonal
  integer :: index_x2a_Faxx_tauy  = 0  ! wind stress, meridional
  integer :: index_x2a_Faxx_lat   = 0  ! latent          heat flux
  integer :: index_x2a_Faxx_sen   = 0  ! sensible        heat flux
  integer :: index_x2a_Faxx_lwup  = 0  ! upward longwave heat flux
  integer :: index_x2a_Faxx_evap  = 0  ! evaporation    water flux
  integer :: index_x2a_So_ustar   = 0  ! surface friction velocity in ocean
  integer :: index_x2a_So_re      = 0  ! square of atm/ocn exch. coeff
  integer :: index_x2a_So_ssq     = 0  ! surface saturation specific humidity in ocean
  integer :: index_x2a_Sx_u10     = 0  ! 10m wind
  integer :: index_x2a_Sx_u10withgusts = 0 ! 10m wind with gusts

contains

  subroutine atm_cpl_indices_set( )

    type(mct_aVect) :: a2x      ! temporary
    type(mct_aVect) :: x2a      ! temporary

    integer, parameter :: tot_mon_in_year = 12
    integer :: imon, ier
    character(len=2) :: monstr ! month string

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2a, rList=seq_flds_x2a_fields, lsize=1)
    call mct_aVect_init(a2x, rList=seq_flds_a2x_fields, lsize=1)

    ! Initialize av indices
    index_x2a_Sx_avsdr      = mct_avect_indexra(x2a,'Sx_avsdr')
    index_x2a_Sx_anidr      = mct_avect_indexra(x2a,'Sx_anidr')
    index_x2a_Sx_avsdf      = mct_avect_indexra(x2a,'Sx_avsdf')
    index_x2a_Sx_anidf      = mct_avect_indexra(x2a,'Sx_anidf')
    index_x2a_Sx_t          = mct_avect_indexra(x2a,'Sx_t')
    index_x2a_So_t          = mct_avect_indexra(x2a,'So_t')
    index_x2a_Sl_snowh      = mct_avect_indexra(x2a,'Sl_snowh')
    index_x2a_Si_snowh      = mct_avect_indexra(x2a,'Si_snowh')

    index_x2a_Sl_fv         = mct_avect_indexra(x2a,'Sl_fv')
    index_x2a_Sl_ram1       = mct_avect_indexra(x2a,'Sl_ram1')
    index_x2a_Sl_soilw      = mct_avect_indexra(x2a,'Sl_soilw',perrWith='quiet')

    index_x2a_Sx_tref       = mct_avect_indexra(x2a,'Sx_tref')
    index_x2a_Sx_qref       = mct_avect_indexra(x2a,'Sx_qref')

    index_x2a_Sf_ifrac      = mct_avect_indexra(x2a,'Sf_ifrac')
    index_x2a_Sf_ofrac      = mct_avect_indexra(x2a,'Sf_ofrac')
    index_x2a_Sf_lfrac      = mct_avect_indexra(x2a,'Sf_lfrac')

    index_x2a_Sx_u10        = mct_avect_indexra(x2a,'Sx_u10')
    index_x2a_Sx_u10withgusts = mct_avect_indexra(x2a,'Sx_u10withgusts')
    index_x2a_Faxx_taux     = mct_avect_indexra(x2a,'Faxx_taux')
    index_x2a_Faxx_tauy     = mct_avect_indexra(x2a,'Faxx_tauy')
    index_x2a_Faxx_lat      = mct_avect_indexra(x2a,'Faxx_lat')
    index_x2a_Faxx_sen      = mct_avect_indexra(x2a,'Faxx_sen')
    index_x2a_Faxx_lwup     = mct_avect_indexra(x2a,'Faxx_lwup')
    index_x2a_Faxx_evap     = mct_avect_indexra(x2a,'Faxx_evap')
    index_x2a_So_ustar      = mct_avect_indexra(x2a,'So_ustar')
    index_x2a_So_re         = mct_avect_indexra(x2a,'So_re')
    index_x2a_So_ssq        = mct_avect_indexra(x2a,'So_ssq')
    index_x2a_Sl_fv         = mct_avect_indexra(x2a,'Sl_fv')
    index_x2a_Sl_ram1       = mct_avect_indexra(x2a,'Sl_ram1')

    index_a2x_Sa_z          = mct_avect_indexra(a2x,'Sa_z')
    index_a2x_Sa_u          = mct_avect_indexra(a2x,'Sa_u')
    index_a2x_Sa_v          = mct_avect_indexra(a2x,'Sa_v')
    index_a2x_Sa_wsresp     = mct_avect_indexra(a2x,'Sa_wsresp', perrWith='quiet')
    index_a2x_Sa_tau_est    = mct_avect_indexra(a2x,'Sa_tau_est', perrWith='quiet')
    index_a2x_Sa_ugust      = mct_avect_indexra(a2x,'Sa_ugust', perrWith='quiet')
    index_a2x_Sa_tbot       = mct_avect_indexra(a2x,'Sa_tbot')
    index_a2x_Sa_ptem       = mct_avect_indexra(a2x,'Sa_ptem')
    index_a2x_Sa_pbot       = mct_avect_indexra(a2x,'Sa_pbot')
    index_a2x_Sa_pslv       = mct_avect_indexra(a2x,'Sa_pslv')
    index_a2x_Sa_shum       = mct_avect_indexra(a2x,'Sa_shum')
    index_a2x_Sa_dens       = mct_avect_indexra(a2x,'Sa_dens')
    index_a2x_Sa_uovern     = mct_avect_indexra(a2x,'Sa_uovern')
    index_a2x_Faxa_swnet    = mct_avect_indexra(a2x,'Faxa_swnet')
    index_a2x_Faxa_lwdn     = mct_avect_indexra(a2x,'Faxa_lwdn')
    index_a2x_Faxa_rainc    = mct_avect_indexra(a2x,'Faxa_rainc')
    index_a2x_Faxa_rainl    = mct_avect_indexra(a2x,'Faxa_rainl')
    index_a2x_Faxa_snowc    = mct_avect_indexra(a2x,'Faxa_snowc')
    index_a2x_Faxa_snowl    = mct_avect_indexra(a2x,'Faxa_snowl')
    index_a2x_Faxa_swndr    = mct_avect_indexra(a2x,'Faxa_swndr')
    index_a2x_Faxa_swvdr    = mct_avect_indexra(a2x,'Faxa_swvdr')
    index_a2x_Faxa_swndf    = mct_avect_indexra(a2x,'Faxa_swndf')
    index_a2x_Faxa_swvdf    = mct_avect_indexra(a2x,'Faxa_swvdf')

    call mct_aVect_clean(x2a)
    call mct_aVect_clean(a2x)

  end subroutine atm_cpl_indices_set

end module atm_cpl_indices
