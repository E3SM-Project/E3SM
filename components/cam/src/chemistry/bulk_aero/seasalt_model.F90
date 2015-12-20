!===============================================================================
! Seasalt for Bulk Aerosol Model
!===============================================================================
module seasalt_model
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,       only: pcols, pver

  implicit none
  private

  public :: seasalt_nbin
  public :: seasalt_nnum
  public :: seasalt_names
  public :: seasalt_indices
  public :: seasalt_init
  public :: seasalt_emis
  public :: seasalt_active

  public :: seasalt_depvel

  logical :: seasalt_active = .false.

  integer, parameter :: seasalt_nbin = 4
  integer, parameter :: seasalt_nnum = 0

  character(len=6), parameter :: seasalt_names(seasalt_nbin) &
       = (/'SSLT01', 'SSLT02', 'SSLT03', 'SSLT04'/)

  integer :: seasalt_indices(seasalt_nbin)

 contains

   !=============================================================================
   !=============================================================================
   subroutine seasalt_init
     use cam_history,   only: addfld, add_default, phys_decomp, fieldname_len
     use constituents,  only: cnst_get_ind

     character(len=fieldname_len) :: dummy
     integer :: m

     do m = 1, seasalt_nbin
        call cnst_get_ind(seasalt_names(m), seasalt_indices(m),abort=.false.)
     enddo
     seasalt_active = any(seasalt_indices(:) > 0)

     if (.not.seasalt_active) return

     dummy = 'RH'
     call addfld (dummy,'frac',pver, 'A','RH in dry dep calc',phys_decomp)
     do m = 1,seasalt_nbin
        dummy = trim(seasalt_names(m)) // 'DI'
        call addfld (dummy,'m/s ',pver, 'A',trim(seasalt_names(m))//' deposition diameter',phys_decomp)
     enddo

   end subroutine seasalt_init

  !=============================================================================
  !=============================================================================
  subroutine seasalt_emis( u10cubed,  srf_temp, ocnfrc, ncol, cflx )

    ! dummy arguments
    real(r8), intent(in) :: u10cubed(:)
    real(r8), intent(in) :: srf_temp(:)
    real(r8), intent(in) :: ocnfrc(:)
    integer,  intent(in) :: ncol
    real(r8), intent(inout) :: cflx(:,:)

    ! local vars
    integer :: ix,m
    real(r8), parameter :: sslt_source(seasalt_nbin) = (/ 4.77e-15_r8, 5.19e-14_r8, 1.22e-13_r8, 6.91e-14_r8 /)

    do m = 1, seasalt_nbin
       ix = seasalt_indices(m)
       cflx(:ncol,ix) = sslt_source(m) * u10cubed(:ncol) * ocnfrc(:ncol)
    enddo
    
  end subroutine seasalt_emis

  !=============================================================================
  !=============================================================================
  subroutine seasalt_depvel( temp, pmid, q, ram1, fv, ncol, lchnk, vlc_dry,vlc_trb,vlc_grv )
    use aerosol_depvel, only: aerosol_depvel_compute
    use wv_saturation,  only: qsat
    use cam_history,    only: outfld
    use mo_constants,   only: dns_aer_sst=>seasalt_density

    integer,  intent(in) :: ncol, lchnk
    real(r8), intent(in) :: temp(:,:)  ! temperature
    real(r8), intent(in) :: pmid(:,:)  ! mid point pressure
    real(r8), intent(in) :: q(:,:)     ! water vapor
    real(r8), intent(in) :: ram1(:)    ! aerodynamical resistance (s/m)
    real(r8), intent(in) :: fv(:)      ! friction velocity (m/s)

    real(r8), intent(out) :: vlc_trb(:,:)    !Turbulent deposn velocity (m/s)
    real(r8), intent(out) :: vlc_grv(:,:,:)  !grav deposn velocity (m/s)
    real(r8), intent(out) :: vlc_dry(:,:,:)  !dry deposn velocity (m/s)

    real(r8) :: r
    real(r8) :: wetdia(pcols,pver,seasalt_nbin)
    real(r8) :: RH(pcols,pver),es(pcols,pver),qs(pcols,pver)  ! for wet radius calculation
    real(r8),parameter:: c1=0.7674_r8, c2=3.0790_r8, c3=2.57e-11_r8,c4=-1.424_r8  ! wet radius calculation constants
    integer :: m, i,k

    ! set stokes correction to 1.0 for now not a bad assumption for our size range)
    real(r8), parameter :: sslt_stk_crc(seasalt_nbin) = (/ 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8 /)
    real(r8), parameter :: sslt_smt_vwr(seasalt_nbin) = (/0.52e-6_r8,2.38e-6_r8,4.86e-6_r8,15.14e-6_r8/) 

    !-----------------------------------------------------------------------

    call qsat(temp(:ncol,:),pmid(:ncol,:),es(:ncol,:),qs(:ncol,:))
    RH(:ncol,:)=q(:ncol,:)/qs(:ncol,:)
    RH(:ncol,:)=max(0.01_r8,min(0.99_r8,RH(:ncol,:)))
    ! set stokes correction to 1.0 for now not a bad assumption for our size range)
    do m=1,seasalt_nbin
       r=sslt_smt_vwr(m)/2.0_r8
       do k=1,pver
          do i=1,ncol
             wetdia(i,k,m)=((r**3+c1*r**c2/(c3*r**c4-log(RH(i,k))))**(1._r8/3._r8))*2.0_r8
          enddo
       enddo
       call outfld( trim(seasalt_names(m))//'DI',wetdia(:,:,m), pcols, lchnk)
    enddo
    call outfld( 'RH',RH(:,:), pcols, lchnk)

    call aerosol_depvel_compute( ncol, pver, seasalt_nbin, temp, pmid, ram1, fv, wetdia, sslt_stk_crc, dns_aer_sst, &
                                 vlc_dry,vlc_trb,vlc_grv)

  endsubroutine seasalt_depvel

end module seasalt_model
