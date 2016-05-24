module ice_cpl_indices
  
  use seq_flds_mod
  use mct_mod

  implicit none

  SAVE
  public                               ! By default make data private

  ! ice -> drv 

  integer :: index_i2x_Si_ifrac        ! fractional ice coverage wrt ocean
  integer :: index_i2x_Si_snowh        ! snow height (m)
  integer :: index_i2x_Si_t            ! temperature                     
  integer :: index_i2x_Si_tref         ! 2m reference temperature        
  integer :: index_i2x_Si_qref         ! 2m reference specific humidity  
  integer :: index_i2x_Si_avsdr        ! albedo: visible, direct         
  integer :: index_i2x_Si_avsdf        ! albedo: near ir, direct         
  integer :: index_i2x_Si_anidr        ! albedo: visible, diffuse        
  integer :: index_i2x_Si_anidf        ! albedo: near ir, diffuse        
  integer :: index_i2x_Si_u10          ! 10m wind
  integer :: index_i2x_Faii_lwup       ! upward longwave heat flux  
  integer :: index_i2x_Faii_lat        ! latent          heat flux  
  integer :: index_i2x_Faii_sen        ! sensible        heat flux      
  integer :: index_i2x_Faii_evap       ! evaporation    water flux      
  integer :: index_i2x_Faii_taux       ! wind stress, zonal            
  integer :: index_i2x_Faii_tauy       ! wind stress, meridional       
  integer :: index_i2x_Faii_swnet      ! sw: net
  integer :: index_i2x_Fioi_swpen      ! sw: net penetrating ice
  integer :: index_i2x_Fioi_melth      ! heat  flux from melting ice (<0)
  integer :: index_i2x_Fioi_meltw      ! water flux from melting ice
  integer :: index_i2x_Fioi_salt       ! salt  flux from meting  ice
  integer :: index_i2x_Fioi_taux       ! ice/ocn stress, zonal
  integer :: index_i2x_Fioi_tauy       ! ice/ocn stress, zonal

  ! drv -> ice

  integer :: index_x2i_So_t            ! ocn layer temperature
  integer :: index_x2i_So_s            ! ocn salinity
  integer :: index_x2i_So_u            ! ocn u velocity
  integer :: index_x2i_So_v            ! ocn v velocity
  integer :: index_x2i_Sa_z            ! bottom atm level height
  integer :: index_x2i_Sa_u            ! bottom atm level zon wind
  integer :: index_x2i_Sa_v            ! bottom atm level mer wind
  integer :: index_x2i_Sa_tbot         ! bottom atm level temp
  integer :: index_x2i_Sa_pbot         ! bottom atm level pressure
  integer :: index_x2i_Sa_ptem         ! bottom atm level pot temp
  integer :: index_x2i_Sa_shum         ! bottom atm level spec hum
  integer :: index_x2i_Sa_dens         ! bottom atm level air den
  integer :: index_x2i_So_dhdx         ! ocn surface slope, zonal
  integer :: index_x2i_So_dhdy         ! ocn surface slope, meridional
  integer :: index_x2i_Faxa_lwdn       ! downward lw heat flux
  integer :: index_x2i_Faxa_rain       ! prec: liquid 
  integer :: index_x2i_Faxa_snow       ! prec: frozen 
  integer :: index_x2i_Faxa_swndr      ! sw: nir direct  downward
  integer :: index_x2i_Faxa_swvdr      ! sw: vis direct  downward
  integer :: index_x2i_Faxa_swndf      ! sw: nir diffuse downward
  integer :: index_x2i_Faxa_swvdf      ! sw: vis diffuse downward
  integer :: index_x2i_Faxa_swnet      ! sw: net
  integer :: index_x2i_Fioo_q          ! ocn freeze or melt heat  
  integer :: index_x2i_Faxa_bcphidry   ! flux: Black Carbon hydrophilic dry deposition
  integer :: index_x2i_Faxa_bcphodry   ! flux: Black Carbon hydrophobic dry deposition
  integer :: index_x2i_Faxa_bcphiwet   ! flux: Black Carbon hydrophilic wet deposition
  integer :: index_x2i_Faxa_ocphidry   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2i_Faxa_ocphodry   ! flux: Organic Carbon hydrophobic dry deposition
  integer :: index_x2i_Faxa_ocphiwet   ! flux: Organic Carbon hydrophilic dry deposition
  integer :: index_x2i_Faxa_dstwet1    ! flux: Size 1 dust -- wet deposition
  integer :: index_x2i_Faxa_dstwet2    ! flux: Size 2 dust -- wet deposition
  integer :: index_x2i_Faxa_dstwet3    ! flux: Size 3 dust -- wet deposition
  integer :: index_x2i_Faxa_dstwet4    ! flux: Size 4 dust -- wet deposition
  integer :: index_x2i_Faxa_dstdry1    ! flux: Size 1 dust -- dry deposition
  integer :: index_x2i_Faxa_dstdry2    ! flux: Size 2 dust -- dry deposition
  integer :: index_x2i_Faxa_dstdry3    ! flux: Size 3 dust -- dry deposition
  integer :: index_x2i_Faxa_dstdry4    ! flux: Size 4 dust -- dry deposition

contains

  subroutine ice_cpl_indices_set( )

    type(mct_aVect) :: i2x      ! temporary
    type(mct_aVect) :: x2i      ! temporary

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2i, rList=seq_flds_x2i_fields, lsize=1)
    call mct_aVect_init(i2x, rList=seq_flds_i2x_fields, lsize=1)

    index_i2x_Si_t          = mct_avect_indexra(i2x,'Si_t')
    index_i2x_Si_tref       = mct_avect_indexra(i2x,'Si_tref')
    index_i2x_Si_qref       = mct_avect_indexra(i2x,'Si_qref')
    index_i2x_Si_ifrac      = mct_avect_indexra(i2x,'Si_ifrac')
    index_i2x_Si_avsdr      = mct_avect_indexra(i2x,'Si_avsdr')
    index_i2x_Si_anidr      = mct_avect_indexra(i2x,'Si_anidr')
    index_i2x_Si_avsdf      = mct_avect_indexra(i2x,'Si_avsdf')
    index_i2x_Si_anidf      = mct_avect_indexra(i2x,'Si_anidf')
    index_i2x_Si_snowh      = mct_avect_indexra(i2x,'Si_snowh')
    index_i2x_Si_u10        = mct_avect_indexra(i2x,'Si_u10')
    index_i2x_Faii_taux     = mct_avect_indexra(i2x,'Faii_taux')
    index_i2x_Faii_tauy     = mct_avect_indexra(i2x,'Faii_tauy')
    index_i2x_Faii_lat      = mct_avect_indexra(i2x,'Faii_lat')
    index_i2x_Faii_sen      = mct_avect_indexra(i2x,'Faii_sen')
    index_i2x_Faii_lwup     = mct_avect_indexra(i2x,'Faii_lwup')
    index_i2x_Faii_evap     = mct_avect_indexra(i2x,'Faii_evap')
    index_i2x_Faii_swnet    = mct_avect_indexra(i2x,'Faii_swnet')
    index_i2x_Fioi_swpen    = mct_avect_indexra(i2x,'Fioi_swpen')
    index_i2x_Fioi_melth    = mct_avect_indexra(i2x,'Fioi_melth')
    index_i2x_Fioi_meltw    = mct_avect_indexra(i2x,'Fioi_meltw')
    index_i2x_Fioi_salt     = mct_avect_indexra(i2x,'Fioi_salt')
    index_i2x_Fioi_taux     = mct_avect_indexra(i2x,'Fioi_taux')
    index_i2x_Fioi_tauy     = mct_avect_indexra(i2x,'Fioi_tauy')

    index_x2i_So_t          = mct_avect_indexra(x2i,'So_t')
    index_x2i_So_s          = mct_avect_indexra(x2i,'So_s')
    index_x2i_So_u          = mct_avect_indexra(x2i,'So_u')
    index_x2i_So_v          = mct_avect_indexra(x2i,'So_v')
    index_x2i_Sa_z          = mct_avect_indexra(x2i,'Sa_z')
    index_x2i_Sa_u          = mct_avect_indexra(x2i,'Sa_u')
    index_x2i_Sa_v          = mct_avect_indexra(x2i,'Sa_v')
    index_x2i_Sa_tbot       = mct_avect_indexra(x2i,'Sa_tbot')
    index_x2i_Sa_ptem       = mct_avect_indexra(x2i,'Sa_ptem')
    index_x2i_Sa_pbot       = mct_avect_indexra(x2i,'Sa_pbot')
    index_x2i_Sa_shum       = mct_avect_indexra(x2i,'Sa_shum')
    index_x2i_Sa_dens       = mct_avect_indexra(x2i,'Sa_dens')
    index_x2i_So_dhdx       = mct_avect_indexra(x2i,'So_dhdx')
    index_x2i_So_dhdy       = mct_avect_indexra(x2i,'So_dhdy')
    index_x2i_Faxa_lwdn     = mct_avect_indexra(x2i,'Faxa_lwdn')
    index_x2i_Faxa_rain     = mct_avect_indexra(x2i,'Faxa_rain')
    index_x2i_Faxa_snow     = mct_avect_indexra(x2i,'Faxa_snow')
    index_x2i_Faxa_swndr    = mct_avect_indexra(x2i,'Faxa_swndr')
    index_x2i_Faxa_swvdr    = mct_avect_indexra(x2i,'Faxa_swvdr')
    index_x2i_Faxa_swndf    = mct_avect_indexra(x2i,'Faxa_swndf')
    index_x2i_Faxa_swvdf    = mct_avect_indexra(x2i,'Faxa_swvdf')
    index_x2i_Fioo_q        = mct_avect_indexra(x2i,'Fioo_q')
    index_x2i_Faxa_bcphidry = mct_avect_indexra(x2i,'Faxa_bcphidry')
    index_x2i_Faxa_bcphodry = mct_avect_indexra(x2i,'Faxa_bcphodry')
    index_x2i_Faxa_bcphiwet = mct_avect_indexra(x2i,'Faxa_bcphiwet')
    index_x2i_Faxa_ocphidry = mct_avect_indexra(x2i,'Faxa_ocphidry')
    index_x2i_Faxa_ocphodry = mct_avect_indexra(x2i,'Faxa_ocphodry')
    index_x2i_Faxa_ocphiwet = mct_avect_indexra(x2i,'Faxa_ocphiwet')
    index_x2i_Faxa_dstdry1  = mct_avect_indexra(x2i,'Faxa_dstdry1')
    index_x2i_Faxa_dstdry2  = mct_avect_indexra(x2i,'Faxa_dstdry2')
    index_x2i_Faxa_dstdry3  = mct_avect_indexra(x2i,'Faxa_dstdry3')
    index_x2i_Faxa_dstdry4  = mct_avect_indexra(x2i,'Faxa_dstdry4')
    index_x2i_Faxa_dstwet1  = mct_avect_indexra(x2i,'Faxa_dstwet1')
    index_x2i_Faxa_dstwet2  = mct_avect_indexra(x2i,'Faxa_dstwet2')
    index_x2i_Faxa_dstwet3  = mct_avect_indexra(x2i,'Faxa_dstwet3')
    index_x2i_Faxa_dstwet4  = mct_avect_indexra(x2i,'Faxa_dstwet4')

    call mct_aVect_clean(x2i)
    call mct_aVect_clean(i2x)

  end subroutine ice_cpl_indices_set

end module ice_cpl_indices
