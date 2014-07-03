!===============================================================================
! Dust for Bulk Aerosol Model
!===============================================================================
module dust_model 
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,   only: masterproc
  use abortutils,   only: endrun

  implicit none
  private

  public :: dust_active
  public :: dust_names
  public :: dust_nbin
  public :: dust_indices
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init

  public :: dust_depvel

  logical :: dust_active = .false.

  integer, parameter :: dust_nbin = 4
  integer, parameter :: dust_nnum = 0

  character(len=6), parameter :: dust_names(dust_nbin) &
       = (/'DST01 ', 'DST02 ', 'DST03 ', 'DST04 '/)

  real(r8), parameter :: dust_dmt_grd(dust_nbin+1) &
       = (/ 0.1e-6_r8, 1.0e-6_r8, 2.5e-6_r8, 5.0e-6_r8, 10.0e-6_r8 /)

  integer  :: dust_indices(dust_nbin)
  real(r8) :: dust_dmt_vwr(dust_nbin)
  real(r8) :: dust_stk_crc(dust_nbin)

  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset
contains

  !=============================================================================
  ! reads dust namelist options
  !=============================================================================
  subroutine dust_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(dust_emis_fact, 1,                   mpir8,   0, mpicom)
    call mpibcast(soil_erod_file, len(soil_erod_file), mpichar, 0, mpicom)
#endif

  end subroutine dust_readnl

  !=============================================================================
  !=============================================================================
  subroutine dust_init()
    use soil_erod_mod, only: soil_erod_init
    use constituents,  only: cnst_get_ind
    use dust_common,   only: dust_set_params

    integer :: n

    do n = 1, dust_nbin
       call cnst_get_ind(dust_names(n), dust_indices(n),abort=.false.)
    end do
    dust_active = any(dust_indices(:) > 0)
    if (.not.dust_active) return

    call  soil_erod_init( dust_emis_fact, soil_erod_file )

    call dust_set_params( dust_nbin, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )

  end subroutine dust_init

  !==============================================================================
  !==============================================================================
  subroutine dust_emis( ncol, lchnk, dust_flux_in, cflx, soil_erod )
    use soil_erod_mod, only : soil_erod_fact
    use soil_erod_mod, only : soil_erodibility

   ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: dust_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

   ! local vars
    integer :: i, m, idst
    real(r8), parameter :: soil_erod_threshold = 0.001_r8
    real(r8), parameter :: dust_emis_sclfctr(dust_nbin) &
         = (/ 0.038_r8/0.032456_r8, 0.11_r8/0.174216_r8, 0.17_r8/0.4085517_r8, 0.67_r8/0.384811_r8 /)

    ! set dust emissions

    col_loop: do i =1,ncol

       soil_erod(i) = soil_erodibility( i, lchnk )

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

       ! adjust emissions based on soil erosion 
       do m = 1,dust_nbin

          idst = dust_indices(m)
          cflx(i,idst) = -dust_flux_in(i,m) &
               * dust_emis_sclfctr(m)*soil_erod(i)/soil_erod_fact*1.15_r8

       enddo

    end do col_loop

  end subroutine dust_emis

  !===============================================================================
  !===============================================================================
  subroutine dust_depvel( temp, pmid, ram1, fv, ncol,  vlc_dry,vlc_trb,vlc_grv )
    use aerosol_depvel, only: aerosol_depvel_compute
    use mo_constants,   only: dust_density
    use ppgrid,         only: pver

    real(r8), intent(in) :: temp(:,:)  ! temperature
    real(r8), intent(in) :: pmid(:,:)  ! mid point pressure
    real(r8), intent(in) :: ram1(:)    ! aerodynamical resistance (s/m)
    real(r8), intent(in) :: fv(:)      ! friction velocity (m/s)
    integer,  intent(in) :: ncol

    real(r8), intent(out) :: vlc_trb(:,:)    !Turbulent deposn velocity (m/s)
    real(r8), intent(out) :: vlc_grv(:,:,:)  !grav deposn velocity (m/s)
    real(r8), intent(out) :: vlc_dry(:,:,:)  !dry deposn velocity (m/s)

    real(r8) :: diam(ncol,pver,dust_nbin)
    integer :: m

    do m=1,dust_nbin
       diam(:,:,m) = dust_dmt_vwr(m)
    enddo
    call aerosol_depvel_compute( ncol, pver, dust_nbin, temp, pmid, ram1, fv, diam, dust_stk_crc, dust_density, &
                                 vlc_dry,vlc_trb,vlc_grv)
  endsubroutine dust_depvel

end module dust_model
