!===============================================================================
!===============================================================================
module soil_erod_mod
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc
  use abortutils,   only: endrun

  implicit none
  private

  public :: soil_erod_init
  public :: soil_erodibility
  public :: soil_erod_fact

  real(r8), allocatable ::  soil_erodibility(:,:)  ! soil erodibility factor
  real(r8) :: soil_erod_fact                       ! tuning parameter for dust emissions

contains

  !=============================================================================
  !=============================================================================
  subroutine soil_erod_init( dust_emis_fact, soil_erod_file )
    use interpolate_data, only: lininterp_init, lininterp, lininterp_finish, interp_type
    use ppgrid,           only: begchunk, endchunk, pcols
    use mo_constants,     only: pi, d2r
    use pio,              only: file_desc_t,pio_inq_dimid,pio_inq_dimlen,pio_get_var,pio_inq_varid, PIO_NOWRITE
    use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use cam_pio_utils,    only: cam_pio_openfile
    use ioFileMod,        only: getfil

    real(r8),         intent(in) :: dust_emis_fact
    character(len=*), intent(in) :: soil_erod_file

    real(r8), allocatable ::  soil_erodibility_in(:,:)  ! temporary input array
    real(r8), allocatable :: dst_lons(:)
    real(r8), allocatable :: dst_lats(:)
    character(len=cl)     :: infile
    integer :: did, vid, nlat, nlon
    type(file_desc_t) :: ncid

    type(interp_type) :: lon_wgts, lat_wgts
    real(r8) :: to_lats(pcols), to_lons(pcols)
    integer :: c, ncols, ierr
    real(r8), parameter :: zero=0._r8, twopi=2._r8*pi

    soil_erod_fact = dust_emis_fact

    ! Summary to log file
    if (masterproc) then
       write(iulog,*) 'soil_erod_mod: soil erodibility dataset: ', trim(soil_erod_file)
       write(iulog,*) 'soil_erod_mod: soil_erod_fact = ', soil_erod_fact
    end if

    ! for soil erodibility in mobilization, apply inside CAM instead of lsm.
    ! read in soil erodibility factors, similar to Zender's boundary conditions

    ! Get file name.  
    call getfil(soil_erod_file, infile, 0)
    call cam_pio_openfile (ncid, trim(infile), PIO_NOWRITE)

    ! Get input data resolution.
    ierr = pio_inq_dimid( ncid, 'lon', did )
    ierr = pio_inq_dimlen( ncid, did, nlon )

    ierr = pio_inq_dimid( ncid, 'lat', did )
    ierr = pio_inq_dimlen( ncid, did, nlat )

    allocate(dst_lons(nlon))
    allocate(dst_lats(nlat))
    allocate(soil_erodibility_in(nlon,nlat))

    ierr = pio_inq_varid( ncid, 'lon', vid )
    ierr = pio_get_var( ncid, vid, dst_lons  )

    ierr = pio_inq_varid( ncid, 'lat', vid )
    ierr = pio_get_var( ncid, vid, dst_lats  )

    ierr = pio_inq_varid( ncid, 'mbl_bsn_fct_geo', vid )
    ierr = pio_get_var( ncid, vid, soil_erodibility_in )

    !-----------------------------------------------------------------------
    !     	... convert to radians and setup regridding
    !-----------------------------------------------------------------------
    dst_lats(:) = d2r * dst_lats(:)
    dst_lons(:) = d2r * dst_lons(:)

    allocate( soil_erodibility(pcols,begchunk:endchunk), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to allocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to allocate soil_erodibility_in')
    end if

    !-----------------------------------------------------------------------
    !     	... regrid ..
    !-----------------------------------------------------------------------
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)

       call lininterp_init(dst_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(dst_lats, nlat, to_lats, ncols, 1, lat_wgts)

       call lininterp(soil_erodibility_in(:,:), nlon,nlat , soil_erodibility(:,c), ncols, lon_wgts,lat_wgts)

       call lininterp_finish(lat_wgts)
       call lininterp_finish(lon_wgts)
    end do
    deallocate( soil_erodibility_in, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to deallocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to deallocate soil_erodibility_in')
    end if

    deallocate( dst_lats )
    deallocate( dst_lons )

  end  subroutine soil_erod_init

end module soil_erod_mod
