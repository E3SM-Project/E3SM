!----------------------------------------------------------------------------------
! Bulk aerosol implementation
!----------------------------------------------------------------------------------
module sox_cldaero_mod

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use ppgrid,       only : pcols, pver
  use mo_chem_utls, only : get_spc_ndx
  use cldaero_mod,  only : cldaero_conc_t, cldaero_allocate, cldaero_deallocate

  implicit none
  private

  public :: sox_cldaero_init
  public :: sox_cldaero_create_obj
  public :: sox_cldaero_update
  public :: sox_cldaero_destroy_obj

  integer :: id_so2, id_so4, id_h2o2

  real(r8), parameter :: small_value = 1.e-20_r8

contains
  
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

  subroutine sox_cldaero_init

     id_so2 = get_spc_ndx( 'SO2' )
     id_so4 = get_spc_ndx( 'SO4' )
     id_h2o2 = get_spc_ndx( 'H2O2' )

     if ( id_so2<1 ) then 
        call endrun('sox_cldaero_init: SO2 is not included in chemistry -- should not invoke sox_cldaero_mod...')
     endif
   
  end subroutine sox_cldaero_init

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  function sox_cldaero_create_obj(cldfrc, qcw, lwc, cfact, ncol, loffset) result( conc_obj )

    real(r8), intent(in) :: cldfrc(:,:)
    real(r8), intent(in) :: qcw(:,:,:)
    real(r8), intent(in) :: lwc(:,:)
    real(r8), intent(in) :: cfact(:,:)
    integer,  intent(in) :: ncol
    integer,  intent(in) :: loffset

    type(cldaero_conc_t), pointer :: conc_obj

    conc_obj => cldaero_allocate()

    conc_obj%xlwc(:ncol,:) = lwc(:ncol,:)*cfact(:ncol,:) ! cloud water L(water)/L(air)

  end function sox_cldaero_create_obj

!----------------------------------------------------------------------------------
! Update the mixing ratios
!----------------------------------------------------------------------------------
  subroutine sox_cldaero_update( &
       ncol, lchnk, loffset, dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, xlwc, &
       delso4_hprxn, xh2so4, xso4, xso4_init, nh3g, hno3g, xnh3, xhno3, xnh4c,  xno3c, xmsa, xso2, xh2o2, qcw, qin )
    
    ! args 
    
    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk ! chunk id
    integer,  intent(in) :: loffset

    real(r8), intent(in) :: dtime ! time step (sec)

    real(r8), intent(in) :: mbar(:,:) ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: pdel(:,:) 
    real(r8), intent(in) :: press(:,:)
    real(r8), intent(in) :: tfld(:,:)

    real(r8), intent(in) :: cldnum(:,:)
    real(r8), intent(in) :: cldfrc(:,:)
    real(r8), intent(in) :: cfact(:,:)
    real(r8), intent(in) :: xlwc(:,:)

    real(r8), intent(in) :: delso4_hprxn(:,:)
    real(r8), intent(in) :: xh2so4(:,:)
    real(r8), intent(in) :: xso4(:,:)
    real(r8), intent(in) :: xso4_init(:,:)
    real(r8), intent(in) :: nh3g(:,:)
    real(r8), intent(in) :: hno3g(:,:)
    real(r8), intent(in) :: xnh3(:,:)
    real(r8), intent(in) :: xhno3(:,:)
    real(r8), intent(in) :: xnh4c(:,:)
    real(r8), intent(in) :: xmsa(:,:)
    real(r8), intent(in) :: xso2(:,:)
    real(r8), intent(in) :: xh2o2(:,:)
    real(r8), intent(in) :: xno3c(:,:)

    real(r8), intent(inout) :: qcw(:,:,:) ! cloud-borne aerosol (vmr)
    real(r8), intent(inout) :: qin(:,:,:) ! xported species ( vmr )
    
    ! local vars ...
    
    integer :: k
    
    !==============================================================
    !       ... Update the mixing ratios
    !==============================================================
    do k = 1,pver

       if (id_so2>0) then
          qin(:,k,id_so2) =  MAX( xso2(:,k), small_value )
       endif
       if (id_h2o2>0) then
          qin(:,k,id_h2o2)=  MAX( xh2o2(:,k), small_value ) 
       endif

       qin(:,k,id_so4) =  MAX( xso4(:,k), small_value )

    end do

  end subroutine sox_cldaero_update

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  subroutine sox_cldaero_destroy_obj( conc_obj )
    type(cldaero_conc_t), pointer :: conc_obj

    call cldaero_deallocate( conc_obj )

  end subroutine sox_cldaero_destroy_obj

end module sox_cldaero_mod
