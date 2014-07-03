      module mo_snoe
!----------------------------------------------------------------------
! An empirical model of nitric oxide (NO) in the lower thermosphere
! (100 - 150 km altitude), based on measurements from the Student
! Nitric Oxide Explorer (SNOE). Model uses empirical orthogonal functions
! (EOFs) derived from the SNOE dataset to describe spatial variability
! in NO. Model NO is the sum of a mean distribution and EOFs multiplied
! by coefficients based on geophysical parameters. These geophysical
! parameters are day of year, Kp magnetic index and F10.7 solar uv index.
!
! Model is utilized by calling subroutine snoe_3d(), which returns
! a 3-D distribution of NO on geographic coordinates (lon/lat are
! supplied by user). Altitude is fixed to SNOE grid (every 3.33 km).
!
! 
!  DRM 4/03
!----------------------------------------------------------------------

      use shr_kind_mod, only : r8 => shr_kind_r8
      use ppgrid,       only : pcols
      use physconst,    only : pi
      use mo_constants, only : r2d, d2r
      use abortutils,   only : endrun
      use cam_logfile,  only : iulog

      implicit none

      private

      public :: snoe_inti, snoe_timestep_init, set_no_ubc
      public :: ndx_no

      save

      integer,  parameter    :: nmodes  = 3
      real(r8), parameter    :: delz    = 104.6_r8             ! delta z (km)
      real(r8), parameter    :: re      = 6471._r8             ! radial distance from ED at 100 km altitude (km)
      real(r8), parameter    :: delx    = -399.1_r8            ! delta x (km)
      real(r8), parameter    :: dely    = -286.1_r8            ! delta y (km)
      real(r8), parameter    :: twopi   = 2._r8*pi, pid2=0.5_r8*pi
      real(r8), parameter    :: theta_n = d2r*(90._r8 - 78.98_r8)    ! co-latitude of CD North Pole (radians)
      real(r8), parameter    :: phi_n   = d2r*289.1_r8               ! longitude of CD North Pole (radians)

      integer :: nmlat
      integer :: nlev
      integer :: ndx_no

      
!----------------------------------------------------------------------
!	... snoe mean and eof data
!----------------------------------------------------------------------
      real(r8), allocatable :: mlat(:,:)                       ! magnetic latitudes corresponding to geo latitudes (radians)
      real(r8), allocatable :: mlat2d(:)                       ! snoe latitudes (degrees)
      real(r8), allocatable :: lev(:)                          ! snoe levels (km)
      real(r8), allocatable :: no_mean(:,:)                    ! mean no
      real(r8), allocatable :: eofs(:,:,:)                     ! empirical orthogonal ftns
      real(r8), allocatable :: snoe_no(:,:,:)                  ! snoe no interpolated to model lons (molecules/cm^3)

      real(r8) :: cthetan
      real(r8) :: sthetan

      contains

      subroutine snoe_inti(snoe_ubc_file)
!----------------------------------------------------------------------
!	... initialize snoe
!----------------------------------------------------------------------

      use ppgrid,       only : begchunk, endchunk
      use constituents, only : cnst_get_ind, cnst_fixed_ubc

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      character(len=*), intent(in) :: snoe_ubc_file

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer  :: astat

!----------------------------------------------------------------------
!	... do we have no with a fixed ubc ?
!----------------------------------------------------------------------
      call cnst_get_ind( 'NO', ndx_no, abort=.false. )
      if( ndx_no > 0 ) then
         if( .not. cnst_fixed_ubc(ndx_no) ) then
            return
         end if
      else
         return
      end if



      cthetan = cos( theta_n )	
      sthetan = sin( theta_n )
!----------------------------------------------------------------------
!	... read snoe netcdf file
!----------------------------------------------------------------------
      call snoe_rdeof(snoe_ubc_file)
      allocate( snoe_no(pcols,nlev,begchunk:endchunk), &
                mlat(pcols,begchunk:endchunk),stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'snoe_inti: failed to allocate snoe_no,mlat; error = ',astat
	 call endrun
      end if

!----------------------------------------------------------------------
!	... get lon/lat transformed to magnetic coords
!----------------------------------------------------------------------
      call geo2mag

      end subroutine snoe_inti

      subroutine snoe_rdeof(snoe_ubc_file)
!----------------------------------------------------------------------
! 	... read in eofs from netcdf file
!----------------------------------------------------------------------

      use ioFileMod,     only : getfil
      use cam_pio_utils, only : cam_pio_openfile
      use pio,           only : file_desc_t, pio_get_var, pio_closefile, &
           pio_nowrite, pio_inq_varid, pio_inq_dimid, pio_inq_dimlen
      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      character(len=*), intent(in) :: snoe_ubc_file

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer :: istat, ierr
      integer :: dim_id, var_id 
      type(file_desc_t) :: ncid
      character(len=256) :: locfn

!----------------------------------------------------------------------

#ifdef SNOE_DIAGS
      write(iulog,*) 'snoe_rdeof: entered routine'
#endif

!----------------------------------------------------------------------
!	... open the netcdf file
!----------------------------------------------------------------------
      call getfil(snoe_ubc_file, locfn, 0)
      call cam_pio_openfile( ncid, trim(locfn), PIO_NOWRITE)

!----------------------------------------------------------------------
!	... read the snoe dimensions
!----------------------------------------------------------------------
      ierr = pio_inq_dimid( ncid, 'lat', dim_id )
      ierr = pio_inq_dimlen( ncid, dim_id, nmlat )
      ierr = pio_inq_dimid( ncid, 'lev', dim_id )
      ierr = pio_inq_dimlen( ncid, dim_id, nlev )

!----------------------------------------------------------------------
!	... allocate snoe variables
!----------------------------------------------------------------------
      allocate( mlat2d(nmlat), lev(nlev), &
                no_mean(nmlat,nlev), eofs(nmlat,nlev,nmodes), stat=istat )
      if( istat /= 0 ) then
         write(iulog,*) 'snoe_rdeof: failed to allocate mlat2d ... eofs; error = ',istat
	 call endrun
      end if

!----------------------------------------------------------------------
!	... read the snoe variables
!----------------------------------------------------------------------
      ierr = pio_inq_varid( ncid, 'NO', var_id )
      ierr = pio_get_var( ncid, var_id, no_mean )

      ierr = pio_inq_varid( ncid, 'EOF', var_id )
      ierr = pio_get_var( ncid, var_id, (/1,1,1/), (/nmlat, nlev, nmodes/), eofs )

      ierr = pio_inq_varid( ncid, 'lat', var_id )
      ierr = pio_get_var( ncid, var_id, mlat2d )
      mlat2d(:) = d2r*mlat2d(:)

      ierr = pio_inq_varid( ncid, 'z', var_id )
      ierr = pio_get_var( ncid, var_id, lev )

!----------------------------------------------------------------------
!	... close the netcdf file
!----------------------------------------------------------------------
      call pio_closefile( ncid)

      end subroutine snoe_rdeof

      subroutine geo2mag
!----------------------------------------------------------------------
!	... converts geographic latitude and longitude pair
!           to eccentric dipole (ED) geomagnetic coordinates
!----------------------------------------------------------------------

      use ppgrid,    only : begchunk, endchunk
      use phys_grid, only : get_ncols_p, get_rlon_all_p, get_rlat_all_p

      implicit none

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer  :: c, i, j
      integer  :: ncol
      real(r8) :: cdlat
      real(r8) :: wrk
      real(r8) :: smaglat, cmaglat, cmaglon, smaglon
      real(r8) :: edlon, sedlon, cedlon, edlat, cedlat
      real(r8) :: singlat, cosglat
      real(r8) :: rlons(pcols)
      real(r8) :: rlats(pcols)

!----------------------------------------------------------------------
!	... start with centered dipole (CD) coordinates
!----------------------------------------------------------------------
chunk_loop : &
      do c = begchunk,endchunk
         ncol = get_ncols_p( c )
         call get_rlon_all_p( c, ncol, rlons )
         call get_rlat_all_p( c, ncol, rlats )
column_loop : &
         do i = 1,ncol
            singlat = sin( rlats(i) )
            cosglat = cos( rlats(i) )
            wrk     = rlons(i) - phi_n
            cdlat = acos( singlat * cthetan + cosglat * sthetan * cos(wrk) )
            smaglat = sin( cdlat )
            cmaglat = cos( cdlat )
            cmaglon = -(singlat - cthetan * cmaglat) / (sthetan * smaglat)
            smaglon = cosglat * sin(wrk) / smaglat
!----------------------------------------------------------------------
!	... convert CD coords to ED coords- equation (39 a,b) of Fraser-Smith
!----------------------------------------------------------------------
            sedlon = re * smaglat * smaglon - dely
            cedlon = re * smaglat * cmaglon - delx
            edlon  = atan2( sedlon, cedlon )	   	! ED magnetic long. (degrees)
            cedlat = (re * cmaglat - delz) * cos(edlon)
            edlat  = atan2( cedlon, cedlat )		! ED magnetic latitude (degrees)
!----------------------------------------------------------------------
!	... convert co-latitudes into latitudes 
!----------------------------------------------------------------------
            if( edlat < 0._r8 ) then
               edlat = edlat + pi
            end if
            mlat(i,c) = pid2 - edlat
         end do column_loop
      end do chunk_loop

      end subroutine geo2mag 

      subroutine snoe_timestep_init( kp, f107 )
!----------------------------------------------------------------------
!	... driver routine that provides access to empirical model
!           data returned is three dimensional (on geographic coordinates)
!----------------------------------------------------------------------

      use time_manager, only : is_first_step, is_first_restart_step, &
                               is_end_curr_day, get_curr_calday
      use ppgrid,       only : pcols, begchunk, endchunk
      use phys_grid,    only : get_ncols_p
      use spmd_utils,   only : masterproc

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      real(r8), intent(in)  :: kp                   ! solar activity index
      real(r8), intent(in)  :: f107                 ! solar activity index

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer                  :: astat, c
      integer                  :: ncol
      real(r8)                 :: doy                ! day of year
      real(r8), allocatable    :: zm(:,:)            ! zonal mean nitric oxide distribution (molecules/cm^3)

      if( is_first_step() .or. is_first_restart_step() .or. is_end_curr_day() ) then
         allocate( zm(nmlat,nlev),stat=astat )
         if( astat /= 0 ) then
            write(iulog,*) 'snoe_3d: failed to allocate zm; error = ',astat
	    call endrun
         end if
         doy = aint( get_curr_calday() )
#ifdef SNOE_DIAGS
         if( masterproc ) then
	    write(iulog,*) ' '
	    write(iulog,*) 'set_snoe_no: doy,kp,f107 = ',doy,kp,f107
	    write(iulog,*) ' '
         end if
#endif
!----------------------------------------------------------------------
!	... obtain SNOE zonal mean data in geomagnetic coordinates
!----------------------------------------------------------------------
         call snoe_zm( doy, kp, f107, zm )

!----------------------------------------------------------------------
!	... map mean to model longitudes
!----------------------------------------------------------------------
         do c = begchunk,endchunk
            ncol = get_ncols_p( c )
            call snoe_zmto3d( c, ncol, zm, nmlat, nlev )
         end do
         deallocate( zm )
      end if

      end subroutine snoe_timestep_init

      subroutine set_no_ubc( lchunk, ncol, zint, mmr, rho )
!----------------------------------------------------------------------
!	... set no mixing ratio at the model top interface level
!----------------------------------------------------------------------

      use constituents, only : pcnst, cnst_mw

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      integer,  intent(in)    :: lchunk
      integer,  intent(in)    :: ncol
      real(r8), intent(in)    :: zint(pcols)          ! geopot height (km)
      real(r8), intent(in)    :: rho(pcols)           ! total atm density (kg/m^3)
      real(r8), intent(inout) :: mmr(pcols,pcnst)     ! concentration at top interface (kg/kg)

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      real(r8), parameter :: amu_fac  = 1.65979e-27_r8 ! kg/amu
      real(r8), parameter :: cm3_2_m3  = 1.e6_r8
      integer    :: astat, i, k, ks, ks1
      real(r8)   :: zinterp, zfac, mfac

!----------------------------------------------------------------------
!	... map to model levels
!----------------------------------------------------------------------
      if( ndx_no > 0 ) then
column_loop : &
         do i = 1,ncol
	    zinterp = zint(i)
	    if( zinterp >= lev(1) ) then
	       mmr(i,ndx_no) = snoe_no(i,1,lchunk)
	    else if( zint(i) < lev(nlev) ) then
	       mmr(i,ndx_no) = 0._r8
	    else
	       do ks = 2,nlev
	          if( zinterp >= lev(ks) ) then
		     ks1 = ks - 1
		     zfac = (zinterp - lev(ks))/(lev(ks1) - lev(ks))
	             mmr(i,ndx_no) = snoe_no(i,ks,lchunk) &
                                     + zfac*(snoe_no(i,ks1,lchunk) - snoe_no(i,ks,lchunk))
		     exit
	          end if
	       end do
	    end if
         end do column_loop
         mfac = amu_fac * cm3_2_m3 * cnst_mw(ndx_no)
         mmr(:ncol,ndx_no) = mmr(:ncol,ndx_no) * mfac /rho(:ncol)
      end if

      end subroutine set_no_ubc 

      subroutine snoe_zm( doy, kp, f107, zm )
!----------------------------------------------------------------------
!	... calculates zonal mean nitric oxide distribution on a given day 
!           and solar conditions (represented by the f10.7 and kp indices)
!----------------------------------------------------------------------

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      real(r8), intent(in)  :: doy
      real(r8), intent(in)  :: kp 
      real(r8), intent(in)  :: f107
      real(r8), intent(out) :: zm(:,:)

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer  :: k
      real(r8) :: theta0         ! day number in degrees
      real(r8) :: dec            ! solar declination angle
      real(r8) :: m1, m2, m3     ! coefficients for first 3 eofs

#ifdef SNOE_DIAGS
      write(iulog,*) 'snoe_zm: doy, kp, f107 = ',doy,kp,f107
#endif

!----------------------------------------------------------------------
!	... calculate coefficients (m1 to m3) for eofs based on
!           geophysical parametes eof1 - kp 
!----------------------------------------------------------------------

      m1 =  kp * 0.785760_r8 - 1.94262_r8

!----------------------------------------------------------------------
!	... eof2 - declination
!----------------------------------------------------------------------
      theta0 = twopi * (doy - 1._r8) / 365._r8

      dec = .006918_r8 &
          - .399912_r8 * cos(theta0)   + .070257_r8 * sin(theta0) &
          - .006758_r8 * cos(2._r8*theta0) + .000907_r8 * sin(2._r8*theta0) & 
          - .002697_r8 * cos(3._r8*theta0) + .001480_r8 * sin(3._r8*theta0)

      dec = dec * r2d

#ifdef SNOE_DIAGS
      write(iulog,*) 'snoe_zm: dec = ',dec
#endif
      m2 = -.319782_r8 + dec*(.0973109_r8 + dec*(.000489814_r8 - dec*.000103608_r8))
      
!----------------------------------------------------------------------
!	... eof3 - f107 
!----------------------------------------------------------------------
      m3 =  log10(f107) * 6.44069_r8 - 13.9832_r8

#ifdef SNOE_DIAGS
      write(iulog,*) 'snoe_zm: m1, m2, m3 = ',m1,m2,m3
#endif

!----------------------------------------------------------------------
!	... zonal mean distrib. is sum of mean and eofs
!----------------------------------------------------------------------
      do k = 1,nlev
         zm(:,k) = no_mean(:,k) - m1 * eofs(:,k,1) + m2 * eofs(:,k,2) - m3 * eofs(:,k,3) 
      end do

      end subroutine snoe_zm 

      subroutine snoe_zmto3d( lchunk, ncol, zm, nmlat, nlev )
!----------------------------------------------------------------------
!	... interpolate zonal mean on magnetic grid 
!           to 3d field in geographic coords
!----------------------------------------------------------------------

      implicit none

!----------------------------------------------------------------------
!	... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)  :: lchunk
      integer, intent(in)  :: ncol
      integer, intent(in)  :: nmlat
      integer, intent(in)  :: nlev
      real(r8), intent(in) :: zm(nmlat,nlev)             ! zonal mean snoe no concentration

!----------------------------------------------------------------------
!	... local variables
!----------------------------------------------------------------------
      integer :: k

      do k = 1,nlev
         call interpol( nmlat, mlat2d, zm(1,k), ncol, mlat(1,lchunk), snoe_no(1,k,lchunk) )
      end do

      end subroutine snoe_zmto3d

      subroutine interpol( nin, xin, yin, nout, xout, yout )
!-----------------------------------------------------------------------
!	... linear interpolation
!           does not extrapolate, but repeats edge values
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: nin, nout
      real(r8), intent(in)    :: xin(nin)
      real(r8), intent(in)    :: yin(nin)
      real(r8), intent(in)    :: xout(nout)
      real(r8), intent(out)   :: yout(nout)

!-----------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------
      integer :: i, j

      do j = 1,nout
        if( xout(j) < xin(1) ) then
          yout(j) = yin(1)
        else
          if( xout(j) > xin(nin) ) then
            yout(j) = yin(nin)
          else
            do i = 1, nin-1
              if ((xout(j) >= xin(i)) .and. (xout(j) < xin(i+1)) ) then
                yout(j) =  yin(i) &
		           + (yin(i+1) - yin(i)) * (xout(j) - xin(i)) / (xin(i+1) - xin(i))
              end if
            end do
          end if
        end if
      end do

      end subroutine interpol

      end module mo_snoe
