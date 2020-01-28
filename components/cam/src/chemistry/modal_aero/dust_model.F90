!===============================================================================
! Dust for Modal Aerosol Model
!===============================================================================
module dust_model 
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,   only: masterproc
  use cam_abortutils,   only: endrun

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
!LXu@06/2019+++
#if  ( defined MODAL_AERO_4MODE_MOM_BIOP)
  public :: phos_dust_active
  public :: phos_dust_indices
#endif
!LXu@06/2019---
  integer, parameter :: dust_nbin = 2
  integer, parameter :: dust_nnum = 2

!LXu@05/2019+++
#if  ( defined MODAL_AERO_3MODE || defined MODAL_AERO_4MODE || defined MODAL_AERO_4MODE_MOM || defined MODAL_AERO_4MODE_MOM_PFIRE || defined MODAL_AERO_4MODE_MOM_BIOP)
  character(len=6), parameter :: dust_names(dust_nbin+dust_nnum) = (/ 'dst_a1', 'dst_a3', 'num_a1', 'num_a3' /)
  real(r8),         parameter :: dust_dmt_grd(dust_nbin+1) = (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = (/ 0.032_r8,0.968_r8 /)

#if  ( defined MODAL_AERO_4MODE_MOM_BIOP)
!LXu@06/2019+++
  character(len=10), parameter :: phos_dust_names(dust_nbin+2) = (/ 'sp_a1', 'sp_a3', 'isp_a1', 'isp_a3'/)
  real(r8),         parameter :: sp_emis_factor(dust_nbin+2) = (/ 0.1_r8, 0.1_r8, 0.9_r8, 0.9_r8/)
  integer                     :: phos_dust_indices(dust_nbin+2)
#endif
!LXu@06/2019---
#elif ( defined MODAL_AERO_7MODE || defined MODAL_AERO_9MODE )
  character(len=6), parameter :: dust_names(dust_nbin+dust_nnum) = (/ 'dst_a5', 'dst_a7', 'num_a5', 'num_a7' /)
  real(r8),         parameter :: dust_dmt_grd(dust_nbin+1) = (/ 0.1e-6_r8, 2.0e-6_r8, 10.0e-6_r8/)
  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = (/ 0.13_r8, 0.87_r8 /)
#endif

  integer  :: dust_indices(dust_nbin+dust_nnum)
  real(r8) :: dust_dmt_vwr(dust_nbin)
  real(r8) :: dust_stk_crc(dust_nbin)

  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

  logical :: dust_active = .false.

!LXu@06/2019+++
  logical :: phos_dust_active = .false.
!LXu@06/2019---

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

!LXu@04/2019 turn off dust emissions
    use cam_logfile,    only: iulog

    integer :: n

    do n = 1, dust_nbin
       call cnst_get_ind(dust_names(n), dust_indices(n),abort=.false.)
!LXu@06/2019+++
#if  ( defined MODAL_AERO_4MODE_MOM_BIOP)
       call cnst_get_ind(phos_dust_names(n), phos_dust_indices(n),abort=.false.)
#endif
!LXu@06/2019---
    end do
!LXu@06/2019+++
#if  ( defined MODAL_AERO_4MODE_MOM_BIOP)
    phos_dust_indices(dust_nbin+1) = phos_dust_indices(1) + 1
    phos_dust_indices(dust_nbin+2) = phos_dust_indices(2) + 1
#endif
!LXu@06/2019---
    do n = 1, dust_nnum
       call cnst_get_ind(dust_names(dust_nbin+n), dust_indices(dust_nbin+n),abort=.false.)
    enddo 

    dust_active = any(dust_indices(:) > 0)
!LXu@06/2019+++
#if  ( defined MODAL_AERO_4MODE_MOM_BIOP)
    phos_dust_active = any(phos_dust_indices(:) > 0)
    if (masterproc) then
       write(iulog,*) 'dust_active = ', dust_active, dust_indices, 'phos_dust_active = ', phos_dust_active, phos_dust_indices
    end if
    phos_dust_active = .false.
!LXu@06/2019---

!LXu@04/2019 turn off dust emissions
    if (masterproc) then
       write(iulog,*) 'dust_active = ', dust_active, dust_indices, 'phos_dust_active = ', phos_dust_active, phos_dust_indices
    end if
#endif
!    dust_active = .false.
!    if (masterproc) then
!       write(iulog,*) 'dust_active = ', dust_active
!    end if


    if (.not.dust_active) return

    call  soil_erod_init( dust_emis_fact, soil_erod_file )

    call dust_set_params( dust_nbin, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )

  end subroutine dust_init

  !===============================================================================
  !===============================================================================
  subroutine dust_emis( ncol, lchnk, dust_flux_in, cflx, soil_erod )
    use soil_erod_mod, only : soil_erod_fact
    use soil_erod_mod, only : soil_erodibility
    use mo_constants,  only : dust_density
    use physconst,     only : pi
!LXu@06/2019+++
#if  ( defined MODAL_AERO_4MODE_MOM_BIOP)
    use soil_erod_mod, only : frac_sol_p_dust,    frac_insol_p_dust
    use ppgrid,        only : pcols
#endif
!LXu@06/2019---

  ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: dust_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

  ! local vars
    integer :: i, m, idst, inum
    real(r8) :: x_mton
    real(r8),parameter :: soil_erod_threshold = 0.1_r8
!LXu@06/2019+++
#if  ( defined MODAL_AERO_4MODE_MOM_BIOP)
    integer  :: i_sp, i_isp
    real(r8) :: frac_sol_p(pcols)
    real(r8) :: frac_insol_p(pcols)
#endif
!LXu@06/2019---

!LXu@06/2019+++
#if  ( defined MODAL_AERO_4MODE_MOM_BIOP)
      frac_sol_p(:) = 0.0_r8
      frac_insol_p(:) = 0.0_r8
#endif
!LXu@06/2019---
    ! set dust emissions
!LXu@06/2019+++
#if  ( defined MODAL_AERO_4MODE_MOM_BIOP)
    col_loop: do i =1,ncol

       soil_erod(i) = soil_erodibility( i, lchnk )

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

       frac_sol_p(i)   = frac_sol_p_dust(i,lchnk)
       frac_insol_p(i) = frac_insol_p_dust(i,lchnk)
      
        ! rebin and adjust dust emissons..
       do m = 1,dust_nbin
       
	  idst = dust_indices(m)

	  cflx(i,idst) = sum( -dust_flux_in(i,:) ) &
              * dust_emis_sclfctr(m)*soil_erod(i)/soil_erod_fact*1.15_r8

	  x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(m)**3._r8))                

	  inum = dust_indices(m+dust_nbin)

          ! set bioavailable(biop) and unbioavailable (unbiop) phosphorus emission from dust
	  i_sp = phos_dust_indices(m)

	  cflx(i,i_sp) = cflx(i,idst) * frac_sol_p(i) 

	  i_isp = phos_dust_indices(m+dust_nbin)

	  cflx(i,i_isp) = cflx(i,idst) * frac_insol_p(i) 

	  cflx(i,inum) = ( cflx(i,idst) + cflx(i,i_sp) + cflx(i,i_isp) )*x_mton
 
       enddo	  

      end do col_loop
#else
    col_loop: do i =1,ncol

       soil_erod(i) = soil_erodibility( i, lchnk )

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

       ! rebin and adjust dust emissons..
       do m = 1,dust_nbin

          idst = dust_indices(m)

          cflx(i,idst) = sum( -dust_flux_in(i,:) ) &
               * dust_emis_sclfctr(m)*soil_erod(i)/soil_erod_fact*1.15_r8

          x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(m)**3._r8))                

          inum = dust_indices(m+dust_nbin)

          cflx(i,inum) = cflx(i,idst)*x_mton
       enddo	  

     end do col_loop
#endif
!LXu@06/2019--


  end subroutine dust_emis

end module dust_model
