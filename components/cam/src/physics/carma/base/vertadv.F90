! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates vertrical advection rates using
!!  Piecewise Polynomial Method [Colela and Woodard, J. Comp. Phys.,
!!  54, 174-201, 1984]
!!
!! @author Eric Jensen
!! @version Dec-1996
subroutine vertadv(carma, cstate, vtrans, cvert, itbnd, ibbnd, cvert_tbnd, cvert_bbnd, vertadvu, vertadvd, rc)

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
  real(kind=f), intent(inout)          :: vtrans(NZP1)    !! vertical velocity
  real(kind=f), intent(in)             :: cvert(NZ)       !! quantity being transported
  integer, intent(in)                  :: itbnd           !! top boundary condition
  integer, intent(in)                  :: ibbnd           !! bottom boundary condition
  real(kind=f), intent(in)             :: cvert_tbnd      !! quantity at top boundary
  real(kind=f), intent(in)             :: cvert_bbnd      !! quantity at bottom boundary
  real(kind=f), intent(out)            :: vertadvu(NZP1)  !! upward vertical transport rate into level k from level k-1 [cm/s]
  real(kind=f), intent(out)            :: vertadvd(NZP1)  !! downward vertical transport rate into level k from level k-1 [cm/s]
  integer, intent(inout)               :: rc              !! return code, negative indicates failure

  ! Local declarations
  integer                        :: k
  integer                        :: nzm1
  integer                        :: nzm2
  integer                        :: itwo
  real(kind=f)                   :: dela(NZ)
  real(kind=f)                   :: delma(NZ)
  real(kind=f)                   :: aju(NZ)
  real(kind=f)                   :: ar(NZ)
  real(kind=f)                   :: al(NZ)
  real(kind=f)                   :: a6(NZ)
  real(kind=f)                   :: dpc, dpc1, dpcm1
  real(kind=f)                   :: ratt1, ratt2, ratt3, rat1, rat2, rat3, rat4, den1
  real(kind=f)                   :: com2, x, xpos
  real(kind=f)                   :: cvert0, cvertnzp1


  ! Initialize fluxes to zero
  vertadvu(:) = 0._f
  vertadvd(:) = 0._f
  
  ! If doing explicit sedimentation then do a simple sorting of positive and negative
  ! velocities into up and down components.
  if (do_explised) then  
    where (vtrans < 0._f)
      vertadvd = -vtrans
    elsewhere
      vertadvu = vtrans
    end where
  else
  
  
    if( ibbnd .eq. I_FLUX_SPEC ) vtrans(1) = 0._f
    if( itbnd .eq. I_FLUX_SPEC ) vtrans(NZP1) = 0._f
            
    ! Set some constants
    nzm1 = max( 1, NZ-1 )
    nzm2 = max( 1, NZ-2 )
    itwo = min( 2, NZ   )
  
    ! First, use cubic fits to estimate concentration values at layer
    ! boundaries
    do k = 2,NZ-1
      dpc = cvert(k) / dz(k)
      dpc1 = cvert(k+1) / dz(k+1)
      dpcm1 = cvert(k-1) / dz(k-1)
      ratt1 = dz(k) / ( dz(k-1) + dz(k) + dz(k+1) )
      ratt2 = ( 2._f*dz(k-1) + dz(k) ) / ( dz(k+1) + dz(k) )
      ratt3 = ( 2._f*dz(k+1) + dz(k) ) / ( dz(k-1) + dz(k) )
      dela(k) = ratt1 * ( ratt2*(dpc1-dpc) + ratt3*(dpc-dpcm1) )
  
      if( (dpc1-dpc)*(dpc-dpcm1) .gt. 0._f .and. dela(k) .ne. 0._f ) then
        delma(k) = min(abs(dela(k)), 2._f*abs(dpc-dpc1), 2._f*abs(dpc-dpcm1)) * abs(dela(k))/dela(k)
      else
        delma(k) = 0._f
      endif
    enddo     ! k = 2,NZ-2
  
    do k = 2,NZ-2
      dpc = cvert(k) / dz(k)
      dpc1 = cvert(k+1) / dz(k+1)
      dpcm1 = cvert(k-1) / dz(k-1)
      rat1 = dz(k) / ( dz(k) + dz(k+1) )
      rat2 = 2._f * dz(k+1) * dz(k) / ( dz(k) + dz(k+1) )
      rat3 = ( dz(k-1) + dz(k) ) / ( 2._f*dz(k) + dz(k+1) )
      rat4 = ( dz(k+2) + dz(k+1) ) / ( 2._f*dz(k+1) + dz(k) )
      den1 = dz(k-1) + dz(k) + dz(k+1) + dz(k+2)
  
      ! <aju(k)> is the estimate for concentration (dn/dz) at layer
      ! boundary <k>+1/2.
      aju(k) = dpc + rat1*(dpc1-dpc) + 1._f/den1 * ( rat2*(rat3-rat4)*(dpc1-dpc) - &
        dz(k)*rat3*delma(k+1) + dz(k+1)*rat4*delma(k) )
  
    enddo     ! k = 2,NZ-2
  
    ! Now construct polynomial functions in each layer
    do k = 3,NZ-2
      al(k) = aju(k-1)
      ar(k) = aju(k)
    enddo
  
    ! Use linear functions in first two and last two layers
    ar(itwo) = aju(itwo)
    al(itwo) = cvert(1)/dz(1) + (zl(itwo)-zc(1)) / &
      (zc(itwo)-zc(1)) * (cvert(itwo)/dz(itwo)-cvert(1)/dz(1))
    ar(1) = al(itwo)
    al(1) = cvert(1)/dz(1) - (zc(1)-zl(1)) / &
      (zc(itwo)-zc(1)) * (cvert(itwo)/dz(itwo)-cvert(1)/dz(1))
  
    al(nzm1) = aju(nzm2)
    ar(nzm1) = cvert(nzm1)/dz(nzm1) + (zl(NZ)-zc(nzm1)) &
      / (zc(NZ)-zc(nzm1)) * (cvert(NZ)/dz(NZ)-cvert(nzm1)/dz(nzm1))
    al(NZ) = ar(nzm1)
    ar(NZ) = cvert(nzm1)/dz(nzm1) + (zl(NZ+1)-zc(nzm1)) &
      / (zc(NZ)-zc(nzm1)) * (cvert(NZ)/dz(NZ)-cvert(nzm1)/dz(nzm1))
  
    ! Ensure that boundary values are not negative
    al(1) = max( al(1), 0._f )
    ar(NZ) = max( ar(NZ), 0._f )
  
    ! Next, ensure that polynomial functions do not deviate beyond the
    ! range [<al(k)>,<ar(k)>]
    do k = 1,NZ
      dpc = cvert(k) / dz(k)
      if( (ar(k)-dpc)*(dpc-al(k)) .le. 0._f ) then
        al(k) = dpc
        ar(k) = dpc
      endif
  
      if( (ar(k)-al(k))*( dpc - 0.5_f*(al(k)+ar(k)) ) .gt. 1._f/6._f*(ar(k)-al(k))**2 ) &
        al(k) = 3._f*dpc - 2._f*ar(k)
  
      if( (ar(k)-al(k))*( dpc - 0.5_f*(al(k)+ar(k)) ) .lt. -1._f/6._f*(ar(k)-al(k))**2 ) &
        ar(k) = 3._f*dpc - 2._f*al(k)
    enddo
  
    ! Calculate fluxes across each layer boundary
    do k = 1,NZ
      dpc = cvert(k) / dz(k)
      dela(k) = ar(k) - al(k)
      a6(k) = 6._f * ( dpc - 0.5_f*(ar(k)+al(k)) )
    enddo
  
    do k = 1,NZ-1
      com2  = ( dz(k) + dz(k+1) ) / 2._f
      x = vtrans(k+1)*dtime/dz(k)
      xpos = abs(x)
  
      ! Upward transport rate
      if( vtrans(k+1) .gt. 0._f )then
  
        if( x .lt. 1._f .and. cvert(k) .ne. 0._f )then
           vertadvu(k+1) = ( vtrans(k+1) * com2 ) * ( ( ar(k) - 0.5_f*dela(k)*x + &
             (x/2._f - (x**2)/3._f)*a6(k) ) / cvert(k) )
  
        ! If Courant # > 1, use upwind advection
        else
          vertadvu(k+1) = vtrans(k+1)
        endif
  
      ! Downward transport rate
      elseif( vtrans(k+1) .lt. 0._f )then
  
        if( x .gt. -1._f .and. cvert(k+1) .ne. 0._f )then
          vertadvd(k+1) = ( -vtrans(k+1) * com2 ) * &
                  ( ( al(k+1) + 0.5_f*dela(k+1)*xpos + &
                  ( xpos/2._f - (xpos**2)/3._f)*a6(k+1) ) / cvert(k+1) )
        else
          vertadvd(k+1) = -vtrans(k+1)
        endif
      endif
  
    enddo    ! k = 1,NZ-1
  
    ! Lower boundary transport rates:  If I_FIXED_CONC boundary
    ! condtion is selected, then use concentration assumed just beyond
    ! the lowest layer edge to calculate the transport rate across
    ! the bottom boundary of the model.
    if( ibbnd .eq. I_FIXED_CONC ) then
  
      com2  = ( dz(1) + dz(itwo) ) / 2._f
      x = vtrans(1)*dtime/dz(1)
      xpos = abs(x)
      cvert0 = cvert_bbnd
      if( vtrans(1) .gt. 0._f )then
  
        if( x .lt. 1._f .and. cvert0 .ne. 0._f )then
          vertadvu(1) = vtrans(1)/cvert0*com2 &
                     * ( ar(1) - 0.5_f*dela(1)*x + &
                     (x/2._f - (x**2)/3._f)*a6(1) )
        else
          vertadvu(1) = vtrans(1)
        endif
  
      elseif( vtrans(1) .lt. 0._f )then
  
        if( x .gt. -1._f .and. cvert(1) .ne. 0._f )then
          vertadvd(1) = -vtrans(1)/ &
                     cvert(1)*com2 &
                     * ( al(1) + 0.5_f*dela(1)*xpos + &
                     (xpos/2._f - (xpos**2)/3._f)*a6(1) )
        else
          vertadvd(1) = -vtrans(1)
        endif
      endif
    endif
  
    ! Upper boundary transport rates
    if( itbnd .eq. I_FIXED_CONC ) then
  
      com2  = ( dz(NZ) + dz(nzm1) ) / 2._f
      x = vtrans(NZ+1)*dtime/dz(NZ)
      xpos = abs(x)
      cvertnzp1 = cvert_tbnd
  
      if( vtrans(NZ+1) .gt. 0._f )then
  
        if( x .lt. 1._f .and. cvert(NZ) .ne. 0._f )then
          vertadvu(NZ+1) = vtrans(NZ+1)/cvert(NZ)*com2 &
                   * ( ar(NZ) - 0.5_f*dela(NZ)*x + &
                   (x/2._f - (x**2)/3._f)*a6(NZ) )
        else
          vertadvu(NZ+1) = vtrans(NZ+1)
        endif
  
      elseif( vtrans(NZ+1) .lt. 0._f )then
  
        if( x .gt. -1._f .and. cvertnzp1 .ne. 0._f )then
          vertadvd(NZ+1) = -vtrans(NZ+1)/ &
                   cvertnzp1*com2 &
                   * ( al(NZ) + 0.5_f*dela(NZ)*xpos + &
                   (xpos/2._f - (xpos**2)/3._f)*a6(NZ) )
        else
          vertadvd(NZ+1) = -vtrans(NZ+1)
        endif
      endif
    endif
  endif
      
  ! Return to caller with vertical transport rates.
  return
end
