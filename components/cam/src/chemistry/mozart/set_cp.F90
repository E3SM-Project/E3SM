
module set_cp

  use shr_kind_mod,      only : r8 => shr_kind_r8
  use physconst,         only : r_universal

  implicit none

  private
  public :: calc_cp

  save

  real(r8), parameter :: ur = .5_r8 * r_universal

contains

  subroutine calc_cp( ncol, vmr, cpairz )
    !-----------------------------------------------------------------------      
    !        ... force ion/electron balance
    !-----------------------------------------------------------------------      

    use ppgrid,       only : pver
    use physconst,    only : cpair
    use chem_mods,    only : adv_mass
    use mo_chem_utls, only : get_spc_ndx, get_inv_ndx

    implicit none

    !-----------------------------------------------------------------------      
    !        ... dummy arguments
    !-----------------------------------------------------------------------      
    integer, intent(in)           :: ncol
    real(r8), intent(in)          :: vmr(:,:,:)         ! species vmrentrations (mol/mol)
    real(r8), intent(inout)       :: cpairz(:,:)

    !-----------------------------------------------------------------------      
    !        ... local variables
    !-----------------------------------------------------------------------      
    integer  :: k, n
    real(r8) :: ro_mw, ro2_mw, rn2_mw                    ! inverse molecular weights
    real(r8) :: hvmr(ncol,pver)                            ! h  vmrentration (mol/mol)
    real(r8) :: n2vmr(ncol,pver)                           ! n2 vmrentration (mol/mol)
    real(r8) :: o3vmr(ncol,pver)                           ! o3 vmrentration (mol/mol)
    real(r8) :: o2vmr(ncol,pver)                           ! o2 vmrentration (mol/mol)
    real(r8) :: ovmr(ncol,pver)                            ! o  vmrentration (mol/mol)

    logical, parameter :: fixed_cp = .true.

    if( fixed_cp ) then
       !-----------------------------------------------------------------------      
       !        ... use same cp as rest of CAM 
       !-----------------------------------------------------------------------      
       cpairz(:ncol,:pver) = cpair 
    else
       !-----------------------------------------------------------------------      
       !        ... caculate cp based on ratio of molecules to atoms 
       !-----------------------------------------------------------------------      
       n = get_spc_ndx( 'O' )
       if( n > 0 ) then
          ro_mw = 1._r8/adv_mass(n)
          do k = 1,pver
             ovmr(:,k) = vmr(:ncol,k,n)
          end do
       else
          ro_mw = 1._r8
          do k = 1,pver
             ovmr(:,k) = 0._r8
          end do
       end if
       n = get_spc_ndx( 'O2' )
       if( n > 0 ) then
          ro2_mw = 1._r8/adv_mass(n)
          do k = 1,pver
             o2vmr(:,k) = vmr(:ncol,k,n)
          end do
       else
          ro2_mw = 1._r8
          do k = 1,pver
             o2vmr(:,k) = 0._r8
          end do
       end if
       n = get_spc_ndx( 'O3' )
       if( n > 0 ) then
          do k = 1,pver
             o3vmr(:,k) = vmr(:ncol,k,n)
          end do
       else
          do k = 1,pver
             o3vmr(:,k) = 0._r8
          end do
       end if
       n = get_spc_ndx( 'H' )
       if( n > 0 ) then
          do k = 1,pver
             hvmr(:,k) = vmr(:ncol,k,n)
          end do
       else
          do k = 1,pver
             hvmr(:,k) = 0._r8
          end do
       end if
       !-----------------------------------------------------------------------      
       !        ... calculate n2 concentration
       !-----------------------------------------------------------------------      
       do k = 1,pver
          n2vmr(:,k) = 1._r8 - (ovmr(:,k) + o2vmr(:,k) + hvmr(:,k))
       end do
       n = get_spc_ndx( 'N' )
       if( n > 0 ) then
          rn2_mw = .5_r8/adv_mass(n)
       else
          rn2_mw = 1._r8
       end if

       !-----------------------------------------------------------------------      
       !        ... calculate cp
       !-----------------------------------------------------------------------      
       do k = 1,pver
          cpairz(:ncol,k) = &
               ur *(7._r8*(o2vmr(:ncol,k)*ro2_mw + n2vmr(:ncol,k)*rn2_mw) &
               + 5._r8*(ovmr(:ncol,k) + o3vmr(:ncol,k))*ro_mw)
       end do
    endif

  end subroutine calc_cp

end module set_cp
