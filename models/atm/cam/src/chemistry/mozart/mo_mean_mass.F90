
module mo_mean_mass

  implicit none

  private
  public :: set_mean_mass, init_mean_mass

  integer :: id_o2, id_o, id_h, id_n

contains

  subroutine init_mean_mass
    use mo_chem_utls, only : get_spc_ndx

    implicit none

    id_o2 = get_spc_ndx('O2')
    id_o  = get_spc_ndx('O')
    id_h  = get_spc_ndx('H')
    id_n  = get_spc_ndx('N')

  endsubroutine init_mean_mass

  subroutine set_mean_mass( ncol, mmr, mbar )
    !-----------------------------------------------------------------
    !        ... Set the invariant densities (molecules/cm**3)
    !-----------------------------------------------------------------

    use shr_kind_mod, only : r8 => shr_kind_r8
    use ppgrid,       only : pver, pcols
    use chem_mods,    only : adv_mass, gas_pcnst
    use physconst,    only : mwdry                   ! molecular weight of dry air
    use abortutils,   only : endrun
    use phys_control, only : waccmx_is               !WACCM-X runtime switch

    implicit none

    !-----------------------------------------------------------------
    !        ... Dummy arguments
    !-----------------------------------------------------------------
    integer, intent(in)   ::      ncol
    real(r8), intent(in)  ::      mmr(:,:,:)           ! species concentrations (kg/kg)
    real(r8), intent(out) ::      mbar(:,:)            ! mean mass (g/mole)

    !-----------------------------------------------------------------
    !        ... Local variables
    !-----------------------------------------------------------------
    integer  :: k
    real(r8) :: xn2(ncol)                                  ! n2 mmr
    real(r8) :: fn2(ncol)                                  ! n2 vmr
    real(r8) :: fo(ncol)                                   ! o  vmr
    real(r8) :: fo2(ncol)                                  ! o2 vmr
    real(r8) :: fh(ncol)                                   ! h vmr
    real(r8) :: ftot(ncol)                                 ! total vmr
    real(r8) :: mean_mass(ncol)                            ! wrk variable

    logical  :: fixed_mbar                                 ! Fixed mean mass flag

    !-------------------------------------------
    !  Mean mass not fixed for WACCM-X
    !-------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
      fixed_mbar = .false.
    else
      fixed_mbar = .true.
    endif

    if( fixed_mbar ) then
       !-----------------------------------------------------------------
       !	... use CAM meam molecular weight 
       !-----------------------------------------------------------------
       mbar(:ncol,:pver) = mwdry  
    else
       if ( id_o2 > 0 .and. id_o > 0 .and. id_h > 0 .and. id_n > 0 ) then
          !-----------------------------------------------------------------
          !	... set the mean mass
          !-----------------------------------------------------------------
          do k = 1,pver
             xn2(:)    = 1._r8 - (mmr(:ncol,k,id_o2) + mmr(:ncol,k,id_o) + mmr(:ncol,k,id_h))
             fn2(:)    = .5_r8 * xn2(:) / adv_mass(id_n)
             fo2(:)    = mmr(:ncol,k,id_o2) / adv_mass(id_o2)
             fo(:)     = mmr(:ncol,k,id_o) / adv_mass(id_o)
             fh(:)     = mmr(:ncol,k,id_h) / adv_mass(id_h)
             mbar(:ncol,k) = 1._r8 / (fn2(:) + fo2(:) + fo(:) + fh(:))
          end do
       else
          call endrun('set_mean_mass: not able to compute mean mass')
       endif
    endif

  end subroutine set_mean_mass

end module mo_mean_mass
