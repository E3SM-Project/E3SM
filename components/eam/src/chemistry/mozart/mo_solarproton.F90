
module mo_solarproton

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use physconst,      only: pi

  implicit none

  save

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine spe_init

    use spedata, only : spedata_init

    implicit none

    !-----------------------------------------------------------------------
    !      ... read in SPE ionization rates
    !-----------------------------------------------------------------------

    call spedata_init()

  end subroutine spe_init

  !-----------------------------------------------------------------------
  !
  !     ... calculates NO production on calday (output in molec/cm3/s)
  !
  !-----------------------------------------------------------------------
  subroutine spe_prod( noxprod, hoxprod, pmid, zmid, lchnk, ncol)

    use mo_apex, only : alatm              ! magnetic latitude grid (radians)
    use ppgrid,  only : pcols, pver
    use spedata, only : get_ionpairs_profile
    use spehox,  only : hox_prod_factor

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) ::  &
         ncol, &                           ! column count
         lchnk                             ! chunk index
    real(r8), intent(in) :: &
         pmid(pcols,pver)                 ! midpoint pressure (Pa)
    real(r8), intent(in) :: &
         zmid(ncol,pver)                  ! midpoint altitude (km)
    real(r8), intent(out) :: &
         noxprod(ncol,pver)                ! NO production
    real(r8), intent(out) :: &
         hoxprod(ncol,pver)               ! HOx production

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------

    integer  :: i
    real(r8) :: dlat_aur
    logical  :: do_spe(ncol)

    real(r8) :: ion_pairs(pver)
    real(r8), parameter :: noxprod_factor = 1._r8 
    real(r8) :: hoxprod_factor(pver)

    !-----------------------------------------------------------------------
    ! 	... intialize NO production
    !-----------------------------------------------------------------------

    noxprod(:ncol,:pver) = 0._r8
    hoxprod(:ncol,:pver) = 0._r8

    !-----------------------------------------------------------------------
    ! 	... check magnetic latitudes, and return if all below 60 deg
    !-----------------------------------------------------------------------
    do i = 1,ncol
       dlat_aur = alatm(i,lchnk)
       do_spe(i) = abs( dlat_aur ) > pi/3._r8
    enddo

    if( all( .not. do_spe(:) ) ) then
       return
    end if

    do i = 1,ncol
       if( do_spe(i) ) then
          call get_ionpairs_profile( pmid(i,:pver), ion_pairs(:pver) )
          noxprod(i,:pver) = noxprod_factor * ion_pairs(:pver)
          hoxprod_factor(:pver) = hox_prod_factor( ion_pairs(:pver), zmid(i,:pver) )
          hoxprod(i,:pver) = hoxprod_factor(:pver)* ion_pairs(:pver)
       end if
    end do

  end subroutine spe_prod

end module mo_solarproton

