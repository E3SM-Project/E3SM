! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates the total production of gases due to nucleation,
!! growth, and evaporation <gasprod> [g/x_units/y_units/z_units/s].
!! It also calculates the latent heating rate from a condensing gas
!! <rlprod> [deg_K/s]
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine gasexchange(carma, cstate, iz, rc)

	! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

	implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(in)                  :: iz      !! z index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: igroup  !! group index
  integer                              :: iepart
  integer                              :: igas    !! gasindex
  integer                              :: i
  integer                              :: i2
  integer                              :: ig2
  integer                              :: ienuc2
  integer                              :: ielem   !! element index
  real(kind=f)                         :: rlh
  real(kind=f)                         :: gasgain
  real(kind=f)                         :: gprod_nuc(NGROUP,NGAS)
  real(kind=f)                         :: gprod_grow(NGROUP,NGAS)


  ! Initialize local variables for keeping track of gas changes due
  ! to nucleation and growth in each particle group.
  gprod_nuc(:,:) = 0._f
  gprod_grow(:,:) = 0._f

  ! First calculate gas loss and latent heat gain rates due to nucleation.
  do igroup = 1,NGROUP

    igas = inucgas(igroup)      ! condensing gas
    ielem = ienconc(igroup)     ! element of particle number concentration

    if( igas .ne. 0 .and. nnuc2elem(ielem) .gt. 0 )then

      do ienuc2 = 1,NELEM

        ig2 = igelem( ienuc2 )    ! target particle group

        if( if_nuc(ielem,ienuc2) ) then

          do i = 1,NBIN
      
            ! If there is no place for the nucleating particle bin to fit in the
            ! nucleated particle, then just skip it. 
            !
            ! This could be an error if significant nucleation really happens from
            ! these bins, but also more flexibility in setting up particle grids.
            gprod_nuc(igroup,igas) = gprod_nuc(igroup,igas) - &
              rhompe(i,ielem) * rmass(i,igroup)

            i2 = inuc2bin(i,igroup,ig2)            ! target bin
            if (i2 /= 0) then
              gprod_nuc(igroup,igas) = gprod_nuc(igroup,igas) - &
                pc(iz,i,ielem) * rnuclg(i,igroup,ig2) * diffmass(i2,ig2,i,igroup)
            end if
          enddo    

          ! Latent heating rate from condensing gas: <rlh> is latent heat of evaporation 
          ! ( + fusion, for ice deposition ) [erg/g]
!          if(( inucproc(ielem,ienuc2) .eq. I_DROPACT ) .or. &
!             ( inucproc(ielem,ienuc2) .eq. I_HOMNUC  )) then
!            rlh = rlhe(iz,igas)
!          elseif(( inucproc(ielem,ienuc2) .eq. I_AERFREEZE ) .or. &
!                 ( inucproc(ielem,ienuc2) .eq. I_HETNUC ))then
!            rlh = rlhe(iz,igas) + rlhm(iz,igas)
!          endif

!          rlprod = rlprod - rlh * gprod_nuc(igroup,igas) / ( CP * rhoa(iz) )
        endif
      enddo     ! ienuc2 = 1,NELEM
    endif      ! (igas = inucgas(ielem) .ne. 0

    ! Next calculate gas lost/gained due to and heat gained/lost from 
    ! growth/evaporation.
    igas = igrowgas(ielem)     ! condensing gas

    if( igas .ne. 0 )then

      do i = 1,NBIN-1

        ! Calculate <gasgain>, mass concentration of gas gained due to evaporation
        ! from each droplet in bin <i+1>.  First check for total evaporation.
        if( totevap(i+1,igroup) )then
          gasgain = ( 1._f - cmf(i+1,igroup) )*rmass(i+1,igroup)
        else
          gasgain = diffmass(i+1,igroup,i,igroup)
        endif

        gprod_grow(igroup,igas) = gprod_grow(igroup,igas) &
          + evaplg(i+1,igroup) * pc(iz,i+1,ielem) * &
            gasgain &
          - growlg(i,igroup) * pc(iz,i,ielem) * &
            diffmass(i+1,igroup,i,igroup)
      enddo    

      ! Add evaporation out of smallest bin (always total evaporation).
      gprod_grow(igroup,igas) = gprod_grow(igroup,igas) + &
        evaplg(1,igroup) * pc(iz,1,ielem) * &
        ( 1._f - cmf(1,igroup) ) * rmass(1,igroup)

      ! Latent heating rate from condensing gas: <rlh> is latent heat of evaporation 
      ! ( + fusion, for ice deposition ) [erg/g]
!      if( is_grp_ice(igroup) )then
!        rlh = rlhe(iz,igas) + rlhm(iz,igas)
!      else
!        rlh = rlhe(iz,igas)
!      endif

!      rlprod = rlprod - rlh * gprod_grow(igroup,igas) / &
!               ( CP * rhoa(iz) )
    endif      ! (igas = igrowgas(ielem)) .ne. 0
  enddo        ! igroup=1,NGROUP

  ! Sum up gas production from nucleation and growth terms.
  do igas = 1,NGAS
    do igroup = 1,NGROUP
      gasprod(igas) = gasprod(igas) + &
         gprod_nuc(igroup,igas) + gprod_grow(igroup,igas)
    enddo
  enddo

  ! Return to caller with <gasprod> evaluated.
  return
end
