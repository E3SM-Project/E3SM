
!-------------------------------------------------------------------
! rebins the 4 sea salt bins into 2 bins for the radiation
!
!  N.B. This code looks for the constituents of SSLTA and SSLTC
!       in the physics buffer first, and uses those if found.
!       Consequently, it is not possible to have prognostic sea
!       salt be radiatively active if the prescribed sea salt is
!       also present.  The current (cam3_5_52) chemistry configurations
!       don't allow both prescribed and prognostic to be present
!       simultaneously, but a more flexible chemistry package that
!       allows this would break this code.
!
! Created by: Francis Vitt
! Date: 9 May 2008
!-------------------------------------------------------------------
module sslt_rebin

  use shr_kind_mod,   only: r8 => shr_kind_r8

  implicit none

  integer :: indices(4)
  integer :: sslta_idx, ssltc_idx

  logical :: has_sslt = .false.
  character(len=1) :: source
  character(len=1), parameter :: DATA = 'D'
  character(len=1), parameter :: PROG = 'P'

  private
  public :: sslt_rebin_init, sslt_rebin_adv, sslt_rebin_register
contains


!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sslt_rebin_register
    use ppgrid,       only : pver,pcols
    
    use physics_buffer, only : pbuf_add_field, dtype_r8

    ! add SSLTA and SSLTC to physics buffer
    call pbuf_add_field('SSLTA','physpkg',dtype_r8,(/pcols,pver/),sslta_idx)
    call pbuf_add_field('SSLTC','physpkg',dtype_r8,(/pcols,pver/),ssltc_idx)

  endsubroutine sslt_rebin_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sslt_rebin_init()

    use constituents, only : cnst_get_ind
    
    use physics_buffer, only : pbuf_get_index, pbuf_set_field, physics_buffer_desc
    use ppgrid,       only : pver
    use cam_history,  only : addfld

    implicit none

    integer :: errcode


    indices(1) = pbuf_get_index('sslt1',errcode)
    indices(2) = pbuf_get_index('sslt2',errcode)
    indices(3) = pbuf_get_index('sslt3',errcode)
    indices(4) = pbuf_get_index('sslt4',errcode)

    has_sslt = all( indices(:) > 0 )
    if ( has_sslt ) source = DATA

    if ( .not. has_sslt ) then
       call cnst_get_ind ('SSLT01', indices(1), abort=.false.)
       call cnst_get_ind ('SSLT02', indices(2), abort=.false.)
       call cnst_get_ind ('SSLT03', indices(3), abort=.false.)
       call cnst_get_ind ('SSLT04', indices(4), abort=.false.)
       has_sslt = all( indices(:) > 0 )
       if ( has_sslt ) source = PROG
    endif

    if ( has_sslt ) then
       call addfld('SSLTA', (/ 'lev' /), 'A','kg/kg', 'sea salt' )
       call addfld('SSLTC', (/ 'lev' /), 'A','kg/kg', 'sea salt' )
    endif

  end subroutine sslt_rebin_init
  
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sslt_rebin_adv(pbuf,  phys_state)

    use physics_types,only : physics_state
    
    use ppgrid,       only : pver, pcols
    use cam_history,  only : outfld
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    implicit none

    
    type(physics_state), target, intent(in) :: phys_state
    type(physics_buffer_desc), pointer :: pbuf(:)

!++ changed wgt_sscm declaration for roundoff validation with earlier code
!    real(r8), parameter :: wgt_sscm = 6.0_r8 / 7.0_r8 ! Fraction of total seasalt mass in coarse mode 
    real(r8), parameter :: wgt_sscm = 6.0_r8 / 7.0_r8 ! Fraction of total seasalt mass in coarse mode 

    real(r8), dimension(:,:), pointer :: sslt1, sslt2, sslt3, sslt4
    real(r8), dimension(:,:), pointer :: sslta, ssltc
    integer :: lchnk, ncol
    real(r8) :: sslt_sum(pcols,pver)

    lchnk = phys_state%lchnk
    ncol = phys_state%ncol

    if (.not. has_sslt) return

    select case( source )
    case (PROG)
       sslt1 => phys_state%q(:,:,indices(1))
       sslt2 => phys_state%q(:,:,indices(2))
       sslt3 => phys_state%q(:,:,indices(3))
       sslt4 => phys_state%q(:,:,indices(4))
    case (DATA)
       call pbuf_get_field(pbuf, indices(1), sslt1)
       call pbuf_get_field(pbuf, indices(2), sslt2)
       call pbuf_get_field(pbuf, indices(3), sslt3)
       call pbuf_get_field(pbuf, indices(4), sslt4)
    end select

    call pbuf_get_field(pbuf, sslta_idx, sslta )
    call pbuf_get_field(pbuf, ssltc_idx, ssltc )

    sslt_sum(:ncol,:) = sslt1(:ncol,:) + sslt2(:ncol,:) + sslt3(:ncol,:) + sslt4(:ncol,:)
    sslta(:ncol,:) = (1._r8-wgt_sscm)*sslt_sum(:ncol,:) ! fraction of seasalt mass in accumulation mode
    ssltc(:ncol,:) = wgt_sscm*sslt_sum(:ncol,:) ! fraction of seasalt mass in coagulation mode

    call outfld( 'SSLTA', sslta(:ncol,:), ncol, lchnk )
    call outfld( 'SSLTC', ssltc(:ncol,:), ncol, lchnk )

  end subroutine sslt_rebin_adv

end module sslt_rebin
