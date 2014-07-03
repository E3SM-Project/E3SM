
module mo_synoz

  !--------------------------------------------------------------------
  !	... synthetic stratospheric ozone emission source
  !--------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog

  implicit none

  save 

  real(r8), allocatable :: po3(:,:,:)

  private
  public :: synoz_inti
  public :: po3

contains

  subroutine synoz_inti( )
    !-----------------------------------------------------------------------
    ! 	... initialize synoz emissions
    !	    note: the emissions are in in units of molecules/cm**3/s
    !-----------------------------------------------------------------------

    use dyn_grid,     only : get_dyn_grid_parm
    use ppgrid,       only : pcols, begchunk, endchunk, pver, pverp
    use ref_pres,     only : pref_mid, pref_edge
    use phys_grid,    only : scatter_field_to_chunk
    use chem_mods,    only : adv_mass
    use mo_chem_utls, only : get_spc_ndx
    use spmd_utils,   only : masterproc
    use physconst,    only : pi, &
                             grav => gravit, & ! m/s^2
                             dry_mass => mwdry, & ! kg/kmole
                             seconds_in_day => cday, &
                             rearth  ! m
    use abortutils,   only : endrun
    use dyn_grid,     only : get_dyn_grid_parm_real1d

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------
    real(r8), parameter :: latmin                = -30._r8
    real(r8), parameter :: latmax                = 30._r8
    real(r8), parameter :: prsmin                = 1000._r8
    real(r8), parameter :: prsmax                = 7000._r8
    real(r8), parameter :: ozone_source_per_year = 500._r8     ! global/annual average of ozone source (Tg/yr)
    real(r8), parameter :: tg2kg                 = 1.e9_r8
    real(r8), parameter :: days_per_year         = 365._r8

    integer :: i, j, k, jl, ju
    integer :: jmin, jmax, kmin, kmax
    integer :: synoz_ndx
    integer :: astat
    integer :: plon, plat, plev

    real(r8)    :: diff, diff_min, diff_max
    real(r8)    :: total_mass = 0._r8
    real(r8)    :: seq
    real(r8), allocatable    :: sf(:)
    real(r8), allocatable :: prs(:)
    real(r8), allocatable :: dp(:)
    real(r8), allocatable :: wk(:,:,:)
    real(r8), allocatable :: lat(:),latwts(:)

    allocate( po3(pcols,pver,begchunk:endchunk), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'synoz_inti: failed to allocate po3; error = ',astat
       call endrun
    end if

    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')
    plev = get_dyn_grid_parm('plev')

    allocate(lat(plat),latwts(plat) )

    lat = get_dyn_grid_parm_real1d('latdeg')
    latwts = get_dyn_grid_parm_real1d('w')

    allocate( wk(plon,plev,plat), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'synoz_inti: failed to allocate wk; error = ',astat
       call endrun
    end if

    Masterproc_only :  if( masterproc ) then
       !-----------------------------------------------------------------------
       ! 	... allocate memory -- for global grid
       !-----------------------------------------------------------------------
       allocate(  prs(plev), dp(plev), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'synoz_inti: failed to allocate prs, dp; error = ',astat
          call endrun
       end if

       !-----------------------------------------------------------------------
       ! 	... find indices of the latitudinal box
       !-----------------------------------------------------------------------
       jmin = -99
       jmax = -99
       diff_min = 1.e20_r8
       diff_max = 1.e20_r8
       do j = 1,plat
          diff = abs(lat(j)-latmin)
          if( diff < diff_min ) then
             diff_min = diff
             jmin = j
          end if
          diff = abs(lat(j)-latmax)
          if( diff < diff_max ) then
             diff_max = diff
             jmax = j
          end if
       end do

       !-----------------------------------------------------------------------
       ! 	... make sure we found them
       !-----------------------------------------------------------------------
       if( jmin < 0 .or. jmax < 0 ) then
          write(iulog,*) 'synoz_inti: problem finding the min/max lat in synoz_inti',jmin,jmax
          call endrun
       end if
       !-----------------------------------------------------------------------
       ! 	... define pressure arrays, assuming surface pressure = 1000 hPa
       !-----------------------------------------------------------------------
       prs(:) = pref_mid(:)
       dp(:)  = pref_edge(2:pverp) - pref_edge(1:pver)

       !-----------------------------------------------------------------------
       ! 	... find indices of the pressure box
       !-----------------------------------------------------------------------
       kmin = -99
       kmax = -99
       diff_min = 1.e20_r8
       diff_max = 1.e20_r8
       do k = 1,plev
          diff = abs( prs(k) - prsmin )
          if( diff < diff_min ) then
             diff_min = diff
             kmin = k
          end if
          diff = abs( prs(k) - prsmax )
          if( diff < diff_max ) then
             diff_max = diff
             kmax = k
          end if
       end do

       !-----------------------------------------------------------------------
       ! 	... make sure we found them
       !-----------------------------------------------------------------------
       if( kmin < 0 .or. kmax < 0 ) then
          write(iulog,*) 'synoz_inti: problem finding the min/max prs in synoz_inti',kmin,kmax
          call endrun
       end if

       !-----------------------------------------------------------------------
       ! 	... define geometric factors (in SI)
       !-----------------------------------------------------------------------
       seq = 2._r8*pi*rearth**2/real(plon)
       allocate(sf(plat))
       do j = 1,plat
          sf(j) = seq*latwts(j)
       end do

       !-----------------------------------------------------------------------
       ! 	... find index of synoz
       !-----------------------------------------------------------------------
       synoz_ndx = get_spc_ndx('SYNOZ')
       has_synoz : if( synoz_ndx > 0 ) then
          !-----------------------------------------------------------------------
          ! 	... compute total mass (in kg) over the domain for which
          !           the ozone source will be defined
          !-----------------------------------------------------------------------
          total_mass = 0._r8
          do k = kmin,kmax
             do j = jmin,jmax
                total_mass = total_mass + dp(k)/grav * sf(j)
             end do
          end do
          total_mass = total_mass * plon
#ifdef DEBUG
          write(iulog,*)'synoz_inti: total mass = ',total_mass
#endif
          !-----------------------------------------------------------------------
          ! 	... define the location of the ozone source
          !-----------------------------------------------------------------------
          wk(:,:,:) = 0._r8
          do k = kmin,kmax
             do j = jmin,jmax
                wk(1:plon,k,j) = 1._r8
             end do
          end do
          !-----------------------------------------------------------------------
          ! 	... define the ozone source as vmr/s (what is needed for the chemistry solver)
          !           note : a change in chemdr is made to avoid the division by invariants
          !-----------------------------------------------------------------------
          wk(:,:,:) = wk(:,:,:) * (ozone_source_per_year*tg2kg/total_mass) &
                                / (seconds_in_day*days_per_year) &
                                * (dry_mass/adv_mass(synoz_ndx))
          write(iulog,*) 'synoz_inti: max wk = ',maxval( wk(:,:,:) )
          deallocate( prs, dp )
       end if has_synoz

       deallocate(sf)

    endif Masterproc_only

    !-----------------------------------------------------------------------
    ! 	... scatter to mpi tasks
    !-----------------------------------------------------------------------
    call scatter_field_to_chunk(1, plev, 1, plon, wk, po3)

    !-----------------------------------------------------------------------
    ! 	... deallocate memory
    !-----------------------------------------------------------------------
    deallocate( wk )
    deallocate( lat, latwts )

  end subroutine synoz_inti

end module mo_synoz
