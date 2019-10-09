! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine defines radius-dependent but time-independent parameters
!! used to calculate condensational growth of particles.  Growth rates
!! are calculated at bin boundaries: the parameters calculated here 
!! ( <gro>, <gro1>, <gro2>, and <akelvin> ) 
!! are defined at lower bin boundaries through the growth rate expression
!! (for one particle) used in growevapl.f:
!!>
!!   dm = gro*pvap*( S + 1 - Ak*As - gro1*gro2*qrad )
!!   --   -------------------------------------------
!!   dt               1 + gro*gro1*pvap
!!
!!  where 
!!
!!  S    = supersaturation
!!  Ak   = exp(akelvin/r)
!!  As   = exp(-sol_ions * solute_mass/solwtmol * gwtmol/condensate_mass)
!!  pvap = saturation vapor pressure [dyne cm**-2]
!!  qrad = radiative energy absorbed
!!<
!! This routine requires that vertical profiles of temperature <T>,
!! and pressure <p> are defined.
!!
!! This routine also requires that particle Reynolds' numbers are
!! defined (setupvfall.f must be called before this).
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine setupgkern(carma, cstate, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc       !! return code, negative indicates failure
  
  ! Local declarations
  integer                        :: igas     !! gas index
  integer                        :: ielem    !! element index
  integer                        :: k        !! z index
  integer                        :: igroup   !! group index
  integer                        :: i
  real(kind=f)                   :: gstick
  real(kind=f)                   :: cor
  real(kind=f)                   :: phish
  real(kind=f)                   :: esh1
  real(kind=f)                   :: a1
  real(kind=f)                   :: br
  real(kind=f)                   :: rknudn
  real(kind=f)                   :: rknudnt
  real(kind=f)                   :: rlam
  real(kind=f)                   :: rlamt
  real(kind=f)                   :: rhoa_cgs(NZ, NGAS)
  real(kind=f)                   :: freep(NZ, NGAS)
  real(kind=f)                   :: freept(NZ, NGAS)
  real(kind=f)                   :: rlh
  real(kind=f)                   :: diffus1
  real(kind=f)                   :: thcond1
  real(kind=f)                   :: reyn_shape
  real(kind=f)                   :: schn
  real(kind=f)                   :: prnum
  real(kind=f)                   :: x1
  real(kind=f)                   :: x2
  real(kind=f)                   :: fv
  real(kind=f)                   :: surf_tens  ! surface tension of H2SO4 particle
  real(kind=f)                   :: rho_H2SO4  ! wet density of H2SO4 particle
  

  ! Calculate gas properties for all of the gases. Better to do them all once, than to
  ! repeat this for multiple groups.
  do igas = 1, NGAS
 
    ! Radius-independent parameters for condensing gas
    !
    ! This is <rhoa> in cgs units.
    !
    rhoa_cgs(:, igas) = rhoa(:) / (xmet(:)*ymet(:)*zmet(:))

    if (igas .eq. igash2o) then
        
      ! Condensing gas is water vapor
      !
      ! <surfctwa> is surface tension of water-air interface (valid from 0 to 40 C)
      ! from Pruppacher and Klett (eq. 5-12).
      surfctwa(:) = 76.10_f - 0.155_f*( t(:) - 273.16_f )

      ! <surfctiw> is surface tension of water-ice interface
      ! from Pruppacher and Klett (eq. 5-48).!
      surfctiw(:) = 28.5_f + 0.25_f*( t(:) - 273.16_f )

      ! <surfctiw> is surface tension of water-ice interface
      ! from Hale and Plummer [J. Chem. Phys., 61, 1974].
      surfctia(:) = 141._f - 0.15_f * t(:)

      ! <akelvin> is argument of exponential in kelvin curvature term.
      akelvin(:,igas) = 2._f*gwtmol(igas)*surfctwa(:) &
                        / ( t(:)*RHO_W*RGAS )

      akelvini(:,igas) = 2._f*gwtmol(igas)*surfctia(:) & 
                        / ( t(:)*RHO_W*RGAS )
                        
    ! condensing gas is H2SO4                    
    else if (igas .eq. igash2so4) then
    
      ! Calculate Kelvin curvature factor for H2SO4 interactively with temperature:
      do k = 1, NZ  
        surf_tens = sulfate_surf_tens(carma, wtpct(k), t(k), rc)
        rho_H2SO4 = sulfate_density(carma, wtpct(k), t(k), rc)
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens / (t(k) * rho_H2SO4 * RGAS)
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvini(k, igash2o)
      end do   
    else 

      ! Condensing gas is not yet configured.
      if (do_print) write(LUNOPRT,*) 'setupgkern::ERROR - invalid igas'
      rc = RC_ERROR
      return
    endif

    ! Molecular free path of condensing gas 
    freep(:,igas)  = 3._f*diffus(:,igas) &
             * sqrt( ( PI*gwtmol(igas) ) / ( 8._f*RGAS*t(:) ) )

    ! Thermal free path of condensing gas
    freept(:,igas) = freep(:,igas)*thcond(:) / &
               ( diffus(:,igas) * rhoa_cgs(:, igas) &
             * ( CP - RGAS/( 2._f*WTMOL_AIR ) ) )
  end do
  

  ! Loop over aerosol groups only (no radius, gas, or spatial dependence).
  do igroup = 1, NGROUP

    ! Use gstickl or gsticki, depending on whether group is ice or not
    if( is_grp_ice(igroup) ) then
      gstick = gsticki
    else
      gstick = gstickl
    endif

    ! Non-spherical corrections (need a reference for these)
    if( ishape(igroup) .eq. I_SPHERE )then

      !   Spheres
      cor = 1._f
      phish = 1._f
    else 

      if( ishape(igroup) .eq. I_HEXAGON )then
        
        ! Hexagons
        phish = 6._f/PI*tan(PI/6._f)*( eshape(igroup) + 0.5_f ) &
               * ( PI / ( 9._f*eshape(igroup)*tan(PI/6._f) ) )**(2._f/3._f)

      else if( ishape(igroup) .eq. I_CYLINDER )then

        ! Spheroids
        phish = ( eshape(igroup) + 0.5_f ) &
               * ( 2._f / ( 3._f*eshape(igroup) ) )**(2._f/3._f)
      endif

      if( eshape(igroup) .lt. 1._f )then

        ! Oblate spheroids
        esh1 = 1._f / eshape(igroup)
        a1 = sqrt(esh1**2 - 1._f)
        cor = a1 / asin( a1 / esh1 ) / esh1**(2._f/3._f)
      else 

        ! Prolate spheroids
        a1 = sqrt( eshape(igroup)**2 - 1._f )
            cor = a1 / log( eshape(igroup) + a1 ) &
                / eshape(igroup)**(ONE/3._f)
      endif
    endif

    ! Evaluate growth terms only for particle elements that grow.
    ! particle number concentration element
    ielem = ienconc(igroup)
    
    ! condensing gas is <igas>
    igas = igrowgas(ielem)
    
    ! If the group doesn't grow, but is involved in aerosol
    ! freezing, then the gas properties still need to be calculated.
    if( igas .eq. 0 ) igas = inucgas(igroup)

    if( igas .ne. 0 )then

      do k = 1, NZ

        ! Latent heat of condensing gas 
        if( is_grp_ice(igroup) )then
          rlh = rlhe(k,igas) + rlhm(k,igas)
        else
          rlh = rlhe(k,igas)
        endif

        ! Radius-dependent parameters 
        do i = 1, NBIN

          br = rlow_wet(k,i,igroup)     ! particle bin Boundary Radius

          ! These are Knudsen numbers
          rknudn  = freep(k,igas) / br
          rknudnt = freept(k,igas) / br

          ! These are "lambdas" used in correction for gas kinetic effects.
          rlam  = ( 1.33_f*rknudn  + 0.71_f ) / ( rknudn  + 1._f ) &
                + ( 4._f*( 1._f - gstick ) ) / ( 3._f*gstick )

          rlamt = ( 1.33_f*rknudnt + 0.71_f ) / ( rknudnt + 1._f ) &
                + ( 4._f*( 1._f - tstick ) ) / ( 3._f*tstick )

          ! Diffusion coefficient and thermal conductivity modified for
          ! free molecular limit and for particle shape.
          diffus1 = diffus(k,igas)*cor / ( 1._f + rlam*rknudn*cor/phish )
          thcond1 = thcond(k)*cor / ( 1._f + rlamt*rknudnt*cor/phish )

          ! Save the modified thermal conductivity off so it can be used in pheat.
          thcondnc(k,i,igroup) = thcond1
          
          ! Reynolds' number based on particle shape <reyn_shape>
          if( ishape(igroup) .eq. I_SPHERE )then
            reyn_shape = re(k,i,igroup)

          else if( eshape(igroup) .lt. 1._f )then
            reyn_shape = re(k,i,igroup) * ( 1._f + 2._f*eshape(igroup) )

          else
            reyn_shape = re(k,i,igroup) * PI*( 1._f+2._f*eshape(igroup) ) &
                       / ( 2._f*( 1._f + eshape(igroup) ) )
          endif

          ! Particle Schmidt number
          schn = rmu(k) / ( rhoa_cgs(k,igas) * diffus1 )

          ! Prandtl number
          prnum = rmu(k)*CP/thcond1

          ! Ventilation factors <fv> and <ft> from Pruppacher and Klett
          x1 = schn **(ONE/3._f) * sqrt( reyn_shape )
          x2 = prnum**(ONE/3._f) * sqrt( reyn_shape )

          if( is_grp_ice(igroup) )then

            ! Ice crystals
            if( x1 .le. 1._f )then
              fv = 1._f   + 0.14_f*x1**2
            else
              fv = 0.86_f + 0.28_f*x1
            endif

            if( x2 .le. 1._f )then
              ft(k,i,igroup) = 1._f   + 0.14_f*x2**2
            else
              ft(k,i,igroup) = 0.86_f + 0.28_f*x2
            endif
          else
          
            ! Liquid water drops
            if( x1 .le. 1.4_f  )then
              fv = 1._f   + 0.108_f*x1**2
            else
              fv = 0.78_f + 0.308_f*x1
            endif

            if( x2 .le. 1.4_f )then
              ft(k,i,igroup) = 1._f   + 0.108_f*x2**2
            else
              ft(k,i,igroup) = 0.78_f + 0.308_f*x2
            endif
          endif

          ! Growth kernel for particle without radiation or heat conduction at
          ! radius lower boundary [g cm^3 / erg / s]
          gro(k,i,igroup) = 4._f*PI*br &
                        * diffus1*fv*gwtmol(igas) &
                        / ( BK*t(k)*AVG )
  
          ! Coefficient for conduction term in growth kernel [s/g]
          gro1(k,i,igroup) = gwtmol(igas)*rlh**2 &
                / ( RGAS*t(k)**2*ft(k,i,igroup)*thcond1 ) &
                / ( 4._f*PI*br )
  
          ! Coefficient for radiation term in growth kernel [g/erg]
          ! (note: no radial dependence).
          if( i .eq. 1 )then
            gro2(k,igroup) = 1._f / rlh
          endif
 
        enddo   ! i=1,NBIN
      enddo    ! k=1,NZ
    endif     ! igas ne 0
  enddo       ! igroup=1,NGROUP

  ! Return to caller with time-independent particle growth 
  ! parameters initialized.
  return
end
