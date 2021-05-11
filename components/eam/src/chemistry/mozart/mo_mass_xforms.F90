module mo_mass_xforms

  use ppgrid,       only : pcols, pver
  use shr_kind_mod, only : r8 => shr_kind_r8


  private
  public :: mmr2vmr, mmr2vmri, vmr2mmr, vmr2mmri, h2o_to_vmr, h2o_to_mmr, init_mass_xforms
  save

  real(r8) :: adv_mass_h2o = 18._r8

contains

  subroutine init_mass_xforms
    use mo_chem_utls, only : get_spc_ndx
    use chem_mods,    only : adv_mass

    implicit none

    integer  :: id_h2o

    id_h2o = get_spc_ndx('H2O')

    if ( id_h2o > 0 ) then
       adv_mass_h2o = adv_mass(id_h2o)
    else
       adv_mass_h2o = 18._r8
    endif

  endsubroutine init_mass_xforms

  subroutine mmr2vmr( mmr, vmr, mbar, ncol )
    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    use chem_mods, only : adv_mass, gas_pcnst

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: mmr(pcols,pver,gas_pcnst)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m

    do m = 1,gas_pcnst
       if( adv_mass(m) /= 0._r8 ) then
          do k = 1,pver
             vmr(:ncol,k,m) = mbar(:ncol,k) * mmr(:ncol,k,m) / adv_mass(m)
          end do
       end if
    end do

  end subroutine mmr2vmr

  subroutine mmr2vmri( mmr, vmr, mbar, mi, ncol )
    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    real(r8), intent(in)    :: mi
    real(r8), intent(in)    :: mbar(:,:)
    real(r8), intent(in)    :: mmr(:,:)
    real(r8), intent(inout) :: vmr(:,:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer  :: k
    real(r8) :: rmi

    rmi = 1._r8/mi
    do k = 1,pver
       vmr(:ncol,k) = mbar(:ncol,k) * mmr(:ncol,k) * rmi
    end do

  end subroutine mmr2vmri

  subroutine vmr2mmr( vmr, mmr, mbar, ncol )
    !-----------------------------------------------------------------
    !	... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    use m_spc_id
    use chem_mods, only : adv_mass, gas_pcnst

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: vmr(ncol,pver,gas_pcnst)
    real(r8), intent(inout) :: mmr(pcols,pver,gas_pcnst)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m

    !-----------------------------------------------------------------
    !	... The non-group species
    !-----------------------------------------------------------------
    do m = 1,gas_pcnst
       if( adv_mass(m) /= 0._r8 ) then
          do k = 1,pver
             mmr(:ncol,k,m) = adv_mass(m) * vmr(:ncol,k,m) / mbar(:ncol,k)
          end do
       end if
    end do

  end subroutine vmr2mmr

  subroutine vmr2mmri( vmr, mmr, mbar, mi, ncol )
    !-----------------------------------------------------------------
    !	... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------
    !	... dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    real(r8), intent(in)    :: mi
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: vmr(ncol,pver)
    real(r8), intent(inout) :: mmr(pcols,pver)

    !-----------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------
    integer :: k, m

    !-----------------------------------------------------------------
    !	... mass to volume mixing for individual species
    !-----------------------------------------------------------------
    do k = 1,pver
       mmr(:ncol,k) = mi * vmr(:ncol,k) / mbar(:ncol,k)
    end do

  end subroutine vmr2mmri

  subroutine h2o_to_vmr( h2o_mmr, h2o_vmr, mbar, ncol )
    !-----------------------------------------------------------------------
    !     ... Transform water vapor from mass to volumetric mixing ratio
    !-----------------------------------------------------------------------

    use chem_mods, only : adv_mass

    implicit none

    !-----------------------------------------------------------------------
    !	... Dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) ::    ncol
    real(r8), dimension(pcols,pver), intent(in) :: &
         h2o_mmr                ! specific humidity ( mmr )
    real(r8), dimension(ncol,pver), intent(in)  :: &
         mbar                   ! atmos mean mass
    real(r8), dimension(ncol,pver), intent(out) :: &
         h2o_vmr                ! water vapor vmr

    !-----------------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------------
    integer ::   k

    do k = 1,pver
       h2o_vmr(:ncol,k) = mbar(:ncol,k) * h2o_mmr(:ncol,k) / adv_mass_h2o
    end do

  end subroutine h2o_to_vmr

  subroutine h2o_to_mmr( h2o_vmr, h2o_mmr, mbar, ncol )
    !-----------------------------------------------------------------------
    !     ... Transform water vapor from volumetric to mass mixing ratio
    !-----------------------------------------------------------------------

    use chem_mods, only : adv_mass

    implicit none

    !-----------------------------------------------------------------------
    !	... Dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) ::    ncol
    real(r8), dimension(ncol,pver), intent(in)  :: &
         mbar                   ! atmos mean mass
    real(r8), dimension(ncol,pver), intent(in)  :: &
         h2o_vmr               ! water vapor vmr
    real(r8), dimension(pcols,pver), intent(out) :: &
         h2o_mmr                ! specific humidity ( mmr )

    !-----------------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------------
    integer ::   k

    do k = 1,pver
       h2o_mmr(:ncol,k) = h2o_vmr(:ncol,k) * adv_mass_h2o / mbar(:ncol,k)
    end do

  end subroutine h2o_to_mmr

end module mo_mass_xforms
