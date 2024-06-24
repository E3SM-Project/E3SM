!===============================================================================
! Dust for Modal Aerosol Model
!===============================================================================
module dust_model 
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,   only: masterproc
  use cam_abortutils,   only: endrun
  use shr_dust_mod,     only: dust_emis_scheme

  implicit none
  private

  public :: dust_names
  public :: dust_nbin
  public :: dust_nnum
  public :: dust_indices
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init
  public :: dust_active

#if ( defined MOSAIC_SPECIES )
  integer, parameter :: dust_nbin = 6
  integer, parameter :: dust_nnum = 2
#else
  integer, parameter :: dust_nbin = 2
  integer, parameter :: dust_nnum = 2
#endif

#if ( ( defined MODAL_AERO_3MODE || defined MODAL_AERO_4MODE || defined MODAL_AERO_4MODE_MOM || defined MODAL_AERO_5MODE ) && ( defined MOSAIC_SPECIES ) )
  character(len=6), parameter :: dust_names(dust_nbin+dust_nnum) = (/ 'dst_a1', 'dst_a3', &
                                                                      'ca_a1 ', 'ca_a3 ', &
                                                                      'co3_a1', 'co3_a3', &
                                                                      'num_a1', 'num_a3' /)
  real(r8),         parameter :: dust_dmt_grd(dust_nnum+1) = (/0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = (/ 0.011_r8, 0.989_r8, &
                                                                   0.011_r8, 0.989_r8, &
                                                                   0.011_r8, 0.989_r8 /)
#elif ( defined MODAL_AERO_3MODE || defined MODAL_AERO_4MODE || defined MODAL_AERO_4MODE_MOM )
  character(len=6), parameter :: dust_names(dust_nbin+dust_nnum) = (/ 'dst_a1', 'dst_a3', 'num_a1', 'num_a3' /)
  real(r8),         parameter :: dust_dmt_grd(dust_nbin+1) = (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
! Zender03: fractions of bin (0.1-1) and bin (1-10) in size 0.1-10
!  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = (/ 0.032_r8,0.968_r8 /)
! Kok11: fractions of bin (0.1-1) and bin (1-10) in size 0.1-10
  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = (/ 0.011_r8,0.989_r8 /)
#elif ( defined MODAL_AERO_5MODE )
  character(len=6), parameter :: dust_names(dust_nbin+dust_nnum) = (/ 'dst_a1', 'dst_a3', 'num_a1', 'num_a3' /)
  real(r8),         parameter :: dust_dmt_grd(dust_nbin+1) = (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = (/ 0.011_r8,0.989_r8 /)
#elif ( defined MODAL_AERO_7MODE || defined MODAL_AERO_9MODE )
  character(len=6), parameter :: dust_names(dust_nbin+dust_nnum) = (/ 'dst_a5', 'dst_a7', 'num_a5', 'num_a7' /)
  real(r8),         parameter :: dust_dmt_grd(dust_nbin+1) = (/ 0.1e-6_r8, 2.0e-6_r8, 10.0e-6_r8/)
  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = (/ 0.13_r8, 0.87_r8 /)
#endif

  integer  :: dust_indices(dust_nbin+dust_nnum)
#if ( defined MOSAIC_SPECIES )
  real(r8) :: dust_dmt_vwr(dust_nnum)
  real(r8) :: dust_stk_crc(dust_nnum)
#else
  real(r8) :: dust_dmt_vwr(dust_nbin)
  real(r8) :: dust_stk_crc(dust_nbin)
#endif

  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

  logical :: dust_active = .false.

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
    use cam_logfile,   only: iulog

    integer :: n

    do n = 1, dust_nbin
       call cnst_get_ind(dust_names(n), dust_indices(n),abrtf=.false.)
    end do
    do n = 1, dust_nnum
       call cnst_get_ind(dust_names(dust_nbin+n), dust_indices(dust_nbin+n),abrtf=.false.)
    enddo 
    dust_active = any(dust_indices(:) > 0)
    if (.not.dust_active) return
   
    call  soil_erod_init( dust_emis_fact, soil_erod_file )

#if ( defined MOSAIC_SPECIES )
    call dust_set_params( dust_nnum, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )
#else
    call dust_set_params( dust_nbin, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )
#endif

    if (masterproc) write(iulog,*) "modal_aero, dust_init: dust_emis_scheme = ",dust_emis_scheme

  end subroutine dust_init

  !===============================================================================
  !===============================================================================
  subroutine dust_emis( ncol, lchnk, dust_flux_in, cflx, soil_erod )
    use soil_erod_mod, only : soil_erod_fact
    use soil_erod_mod, only : soil_erodibility
    use mo_constants,  only : dust_density
    use physconst,     only : pi

  ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: dust_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

  ! local vars
    integer :: i, m, idst, inum
    real(r8) :: x_mton
    real(r8),parameter :: soil_erod_threshold = 0.1_r8
#if ( defined MOSAIC_SPECIES )
    real(r8),parameter :: frc_caco3 = 0.05_r8                           ! fraction of dust emitted as caco3
    real(r8),parameter :: frc_ca    = frc_caco3 * 0.4004308_r8          ! fraction of dust emitted as ca
    real(r8),parameter :: frc_co3   = frc_caco3 - frc_ca                ! fraction of dust emitted as co3
    real(r8),parameter :: frc_oin   = 1.0_r8 - frc_caco3                ! fraction of dust emitted as oin (=dst)
#endif

    ! set dust emissions

    col_loop: do i =1,ncol

       soil_erod(i) = soil_erodibility( i, lchnk )

       if (dust_emis_scheme == 2) soil_erod(i) = 1._r8

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

#if ( defined MOSAIC_SPECIES )
       idst = dust_indices(1)
       cflx(i,idst) = sum( -dust_flux_in(i,:) ) * 0.73_r8/0.87_r8 * frc_oin * dust_emis_sclfctr(1) * soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(2)
       cflx(i,idst) = sum( -dust_flux_in(i,:) ) * 0.73_r8/0.87_r8 * frc_oin * dust_emis_sclfctr(2) * soil_erod(i)/soil_erod_fact*1.15_r8

       idst = dust_indices(3)
       cflx(i,idst) = sum( -dust_flux_in(i,:) ) * 0.73_r8/0.87_r8 * frc_ca  * dust_emis_sclfctr(3) * soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(4)
       cflx(i,idst) = sum( -dust_flux_in(i,:) ) * 0.73_r8/0.87_r8 * frc_ca  * dust_emis_sclfctr(4) * soil_erod(i)/soil_erod_fact*1.15_r8

       idst = dust_indices(5)
       cflx(i,idst) = sum( -dust_flux_in(i,:) ) * 0.73_r8/0.87_r8 * frc_co3 * dust_emis_sclfctr(5) * soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(6)
       cflx(i,idst) = sum( -dust_flux_in(i,:) ) * 0.73_r8/0.87_r8 * frc_co3 * dust_emis_sclfctr(6) * soil_erod(i)/soil_erod_fact*1.15_r8

       inum = dust_indices(dust_nbin+1)
       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(1)**3._r8))
       cflx(i,inum) = ( cflx(i,dust_indices(1)) + cflx(i,dust_indices(3)) + cflx(i,dust_indices(5)) ) * x_mton

       inum = dust_indices(dust_nbin+2)
       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(2)**3._r8))
       cflx(i,inum) = ( cflx(i,dust_indices(2)) + cflx(i,dust_indices(4)) + cflx(i,dust_indices(6)) ) * x_mton
#else
       ! rebin and adjust dust emissons..
       do m = 1,dust_nbin

          idst = dust_indices(m)

       ! Correct the dust input flux calculated by CLM, which uses size distribution in Zender03
       ! to calculate fraction of bin (0.1-10um) in range (0.1-20um) = 0.87
       ! based on Kok11, that fraction is 0.73
!          cflx(i,idst) = sum( -dust_flux_in(i,:) ) &
          cflx(i,idst) = sum( -dust_flux_in(i,:) ) * 0.73_r8/0.87_r8 &
               * dust_emis_sclfctr(m)*soil_erod(i)/soil_erod_fact*1.15_r8

          x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(m)**3._r8))                

          inum = dust_indices(m+dust_nbin)

          cflx(i,inum) = cflx(i,idst)*x_mton

       enddo
#endif

    end do col_loop

  end subroutine dust_emis

end module dust_model
