module mo_airplane
  !--------------------------------------------------------------------
  !	... Airplane insitu emission sources
  !--------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use pio,          only : pio_inq_dimid, pio_inq_dimlen, pio_get_var, &
       file_desc_t, var_desc_t, pio_inq_vardimid, pio_inq_varndims, pio_nowrite, &
       pio_inq_varid, pio_closefile
       
  use cam_pio_utils,only : cam_pio_openfile
  use cam_logfile,  only : iulog
  implicit none

  private

  save 

  real(r8), allocatable :: &
       pno(:,:,:), &
       pco(:,:,:), &
       air_altitude(:)
  public :: airpl_set, airpl_src



  logical :: has_airpl_src = .false.

contains

  subroutine airpl_set( lchnk, ncol, no_ndx, co_ndx, xno_ndx, cldtop, zint_abs, extfrc)
    use ppgrid,       only : pver
    use cam_history,  only : outfld

    implicit none
    integer, intent(in) :: lchnk, ncol, no_ndx, co_ndx, xno_ndx
    real(r8), intent(in) :: cldtop(:), zint_abs(:,:)
    real(r8), intent(inout) :: extfrc(:,:,:)
    

! Local Variables
    real(r8), dimension(ncol,pver) :: no_air, co_air
    real(r8) :: ztab_top, ztab_bot, zdel, zdeli, frac, zlev_top, zlev_bot
    integer :: nlev
    integer :: cldind, kk, i, k

    no_air(:,:) = 0._r8
    co_air(:,:) = 0._r8

    if(has_airpl_src) then
       !---------------------------------------------------------------------
       !     ... Add the airplane emissions; must do vertical interpolation
       !---------------------------------------------------------------------
       ztab_top = maxval( air_altitude )
       ztab_bot = minval( air_altitude )
       nlev     = size(air_altitude) - 1

       !---------------------------------------------------------------------
       !     ... add the airplane emissions; must do vertical interpolation
       !         Note: the interpolation code is conserving and assumes the
       !               aircraft emission vertical grid is uniform with a
       !               one kilometer spacing
       !---------------------------------------------------------------------
       level_loop : do k = 1,pver
          long_loop : do i = 1,ncol
             zlev_top = zint_abs(i,k)                 ! altitude at top of model level (km)
             zlev_bot = zint_abs(i,k+1)               ! altitude at bottom of model level (km)
             zdel = (zlev_top - zlev_bot) * 1.e5_r8  ! model level thickness (cm)
             zdeli = 1._r8/zdel
             if( zlev_bot <= ztab_top .and. zlev_top >= ztab_bot ) then
                do kk = 1,nlev
                   if( zlev_bot <= air_altitude(kk+1) .and. zlev_top >= air_altitude(kk) ) then
                      frac = (min( zlev_top, air_altitude(kk+1) ) - max( zlev_bot, air_altitude(kk) )) &
                           /(air_altitude(kk+1) - air_altitude(kk)) ! *del_alti(kk)
                      if( no_ndx > 0 ) then
                         extfrc(i,k,no_ndx) = extfrc(i,k,no_ndx) + frac * pno(i,kk,lchnk) * zdeli
                         no_air(i,k) = frac * pno(i,kk,lchnk) * zdeli
                      end if
                      if( xno_ndx > 0 ) then
                         extfrc(i,k,xno_ndx) = extfrc(i,k,xno_ndx) + frac * pno(i,kk,lchnk) * zdeli
                      end if
                      if( co_ndx > 0 ) then
                         extfrc(i,k,co_ndx) = extfrc(i,k,co_ndx) + frac * pco(i,kk,lchnk) * zdeli
                         co_air(i,k) = frac * pco(i,kk,lchnk) * zdeli
                      end if
                   end if
                end do
             end if
             if( k == pver ) then
                do kk = 1,nlev
                   if( zlev_bot > air_altitude(kk) ) then
                      frac = (min( zlev_bot, air_altitude(kk+1) ) - air_altitude(kk)) &
                           /(air_altitude(kk+1) - air_altitude(kk)) ! *del_alti(kk)
                      if( no_ndx > 0 ) then
                         extfrc(i,k,no_ndx) = extfrc(i,k,no_ndx) + frac * pno(i,kk,lchnk) * zdeli
                         no_air(i,k) = frac * pno(i,kk,lchnk) * zdeli
                      end if
                      if( xno_ndx > 0 ) then
                         extfrc(i,k,xno_ndx) = extfrc(i,k,xno_ndx) + frac * pno(i,kk,lchnk) * zdeli
                      end if
                      if( co_ndx > 0 ) then
                         extfrc(i,k,co_ndx) = extfrc(i,k,co_ndx) + frac * pco(i,kk,lchnk) * zdeli
                         co_air(i,k) = frac * pco(i,kk,lchnk) * zdeli
                      end if
                   else
                      exit
                   end if
                end do
             end if
          end do long_loop
       end do level_loop
      
    end if
    call outfld( 'NO_Aircraft',  no_air(:ncol,:), ncol, lchnk )
    call outfld( 'CO_Aircraft',  co_air(:ncol,:), ncol, lchnk )

  end subroutine airpl_set


  subroutine airpl_src( airpl_emis_file )
    !-----------------------------------------------------------------------
    ! 	... Initialize airplane emissions
    !	    Note: The emissions are read in in units of molecules/cm**2/s
    !	          on a vertically resolved grid.
    !		  Conversion to units of molec/cm**3/s is done in SETEXT
    !-----------------------------------------------------------------------
    use spmd_utils,    only : masterproc
    use interpolate_data, only : lininterp_init, lininterp, lininterp_finish, &
         interp_type
    use chem_mods,     only : adv_mass
    use ioFileMod,     only : getfil
    use mo_chem_utls,  only : get_spc_ndx, get_extfrc_ndx
    use phys_grid,     only : get_ncols_p, get_rlat_all_p, get_rlon_all_p, ngcols_p
    use ppgrid,        only : begchunk, endchunk, pcols
    use mo_constants,  only : pi, d2r, rearth
    use phys_gmean,    only : gmean
    implicit none

    !-----------------------------------------------------------------------
    ! 	... Dummy args
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: airpl_emis_file

    !-----------------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------------
    real(r8), parameter :: msq2cmsq = 1.e4_r8, zero=0._r8, twopi=2._r8*pi
    integer  :: ios, k, j
    integer  :: nlat, nlon, nlev, ndims
    integer  :: ierr
    type(file_desc_t) :: piofile
    type(var_desc_t) :: vid
    integer  :: dimid_lat, dimid_lon, dimid_lev
    integer  :: dimid(3)
    real(r8), allocatable :: lat(:), lon(:)
    real(r8), allocatable :: pno_in(:,:,:), pco_in(:,:,:)
    real(r8) :: total(2), tmp
    real(r8) :: factor
    character(len=256) :: locfn
    integer :: co_ndx, no_ndx
    type(interp_type) :: lon_wgts, lat_wgts    
    real(r8) :: to_lats(pcols), to_lons(pcols)
    integer :: ncols, c

    co_ndx = get_extfrc_ndx('CO')
    no_ndx = get_extfrc_ndx('NO')

    if ( co_ndx < 0 .and. no_ndx < 0 ) then
       if( masterproc ) then
          write(iulog,*) 'airpl_src: NO and CO do not have external source --> no aircraft sources will be applied'      
       endif
       return
    endif

    if ( len_trim(airpl_emis_file) == 0 ) then
       return
    endif

    has_airpl_src = .true.

    co_ndx = get_spc_ndx('CO')
    no_ndx = get_spc_ndx('NO')


    !-----------------------------------------------------------------------
    !	... Open NetCDF file
    !-----------------------------------------------------------------------
    call getfil (airpl_emis_file, locfn, 0)
    call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)

    !-----------------------------------------------------------------------
    !       ... Get grid dimensions from file
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( piofile, 'lat', dimid_lat )
    ierr = pio_inq_dimlen( piofile, dimid_lat, nlat )
    allocate( lat(nlat), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: lat allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( piofile, 'lat', vid )
    ierr = pio_get_var( piofile, vid, lat )
    lat(:nlat) = lat(:nlat) * d2r

    ierr = pio_inq_dimid( piofile, 'lon', dimid_lon )
    ierr = pio_inq_dimlen( piofile, dimid_lon, nlon )
    allocate( lon(nlon), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: lon allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( piofile, 'lon', vid )
    ierr = pio_get_var( piofile, vid, lon )
    lon(:nlon) = lon(:nlon) * d2r

    ierr = pio_inq_dimid( piofile, 'altitude', dimid_lev )
    ierr = pio_inq_dimlen( piofile, dimid_lev, nlev )
    allocate( air_altitude(nlev+1), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: air_altitude allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( piofile, 'altitude', vid )
    ierr = pio_get_var( piofile, vid, air_altitude(1:nlev) )
    air_altitude(nlev+1) = air_altitude(nlev) + (air_altitude(nlev) - air_altitude(nlev-1))

    !-----------------------------------------------------------------------
    !       ... Set up regridding
    !-----------------------------------------------------------------------

    allocate( pno_in(nlon,nlat,nlev), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: pno_in allocation error = ',ierr
       call endrun
    end if
    allocate( pco_in(nlon,nlat,nlev), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: pco_in allocation error = ',ierr
       call endrun
    end if
    allocate(pno(pcols,nlev,begchunk:endchunk), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: pno allocation error = ',ierr
       call endrun
    end if
    allocate( pco(pcols,nlev,begchunk:endchunk), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: pco allocation error = ',ierr
       call endrun
    end if

    !-----------------------------------------------------------------------
    !	... Read emissions
    !-----------------------------------------------------------------------
    ierr = pio_inq_varid( piofile, 'nox', vid )
    ierr = pio_inq_varndims( piofile, vid, ndims )
    if( ndims /= 3 ) then
       write(iulog,*) 'airpl_src: variable nox has ndims = ',ndims,', expecting 3'
       call endrun
    end if
    ierr = pio_inq_vardimid( piofile, vid, dimid )
    if( dimid(1) /= dimid_lon  .or. dimid(2) /= dimid_lat .or.  dimid(3) /= dimid_lev ) then
       write(iulog,*) 'airpl_src: Dimensions in wrong order for variable nox'
       write(iulog,*) '...      Expecting (lon, lat, lev)'
       call endrun
    end if
    ierr = pio_get_var( piofile, vid, &
         (/ 1, 1, 1/), &                    ! start
         (/ nlon, nlat, nlev /), &   ! count
         pno_in )

    ierr = pio_inq_varid( piofile, 'co', vid )
    ierr = pio_inq_varndims( piofile, vid, ndims )

    if( ndims /= 3 ) then
       write(iulog,*) 'READ_SFLX: variable co has ndims = ',ndims,', expecting 3'
       call endrun
    end if
    ierr = pio_inq_vardimid( piofile, vid, dimid )
    if( dimid(1) /= dimid_lon .or. dimid(2) /= dimid_lat .or. dimid(3) /= dimid_lev ) then
       write(iulog,*) 'airpl_src: Dimensions in wrong order for variable co'
       write(iulog,*) '...      Expecting (lon, lat, lev)'
       call endrun
    end if
    ierr = pio_get_var( piofile, vid, &
         (/ 1, 1, 1/), &                    ! start
         (/ nlon, nlat, nlev /), &   ! count
         pco_in )
    call pio_closefile( piofile )

    !-----------------------------------------------------------------------
    !	... Regrid emissions
    !-----------------------------------------------------------------------
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)
       call lininterp_init(lon, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(lat, nlat, to_lats, ncols, 1, lat_wgts)

       do k = 1,nlev
          call lininterp(pno_in(:,:,k), nlon, nlat, pno(:,k,c), ncols, lon_wgts, lat_wgts) 
          call lininterp(pco_in(:,:,k), nlon, nlat, pco(:,k,c), ncols, lon_wgts, lat_wgts) 
       enddo
       call lininterp_finish(lon_wgts)
       call lininterp_finish(lat_wgts)
    enddo

    deallocate( pno_in, pco_in, lon, lat, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'airpl_src: Failed to deallocate pno_in,pco_in; ierr = ',ierr
       call endrun
    end if
    !-----------------------------------------------------------------------
    !       ... Get global emission from this source
    !-----------------------------------------------------------------------
    total = zero
    do k=1,nlev
       call gmean(pno(:,k,:), tmp)
       total(1)= total(1)+tmp
       call gmean(pco(:,k,:), tmp)
       total(2)= total(2)+tmp
    end do

    if(masterproc) then

       factor = 86400._r8 * 365._r8 &   ! sec / year
            / 6.022e23_r8 &           ! molec / mole
            * 1.e-12_r8   &            ! Tg / g
            * msq2cmsq    &             ! meters**2 to cm**2
            * 4._r8*pi*rearth*rearth      ! global mean to global total

       write(iulog,*) 'airpl_src: nlev = ',nlev
       !-----------------------------------------------------------------------
       !       ... Convert totals from molec cm^-2 s^-1 to Tg y^-1
       !-----------------------------------------------------------------------
       if (no_ndx .gt. 0) then
          total(1) = total(1) * adv_mass(no_ndx) * factor
          write(iulog,'('' airpl_src Aircraft emissions: '',a6,'' = '',f10.3,1X,a6)') 'NO',total(1),'TgN/y'
       endif
       if (co_ndx .gt. 0) then
          total(2) = total(2) * adv_mass(co_ndx) * factor
          write(iulog,'('' airpl_src Aircraft emissions: '',a6,'' = '',f10.3,1X,a6)') 'CO',total(2),'Tg/y'
       endif
    end if

  end subroutine airpl_src

end module mo_airplane
