! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle fall velocities, vf(k) [cm s^-1]
!! and reynolds' numbers based on fall velocities, re(j,i,k) [dimensionless].
!! indices correspond to vertical level <k>, bin index <i>, and aerosol
!! group <j>.
!!
!! Method: first use Stokes flow (with Fuchs' size corrections, 
!! valid only for Stokes flow) to estimate fall velocity, then calculate
!! Reynolds' number (Re) (for spheres, Stokes drag coefficient is 24/Re).
!! Then for Re > 1, correct drag coefficient (Cd) for turbulent boundary
!! layer through standard trick to solving the drag problem: 
!! fit y = log( Re ) as a function of x = log( Cd Re^2 ).  
!! We use the data for rigid spheres taken from Figure 10-6 of
!! Pruppacher and Klett (1978):
!!
!!  Re     Cd
!! -----  ------
!!    1    24
!!   10     4.3
!!  100     1.1
!! 1000     0.45
!!
!! Note that we ignore the "drag crisis" at Re > 200,000
!! (as discussed on p. 341 and shown in Fig 10-36 of P&K 1978), where
!! Cd drops dramatically to 0.2 for smooth, rigid spheres, and instead 
!! assume Cd = 0.45 for Re > 1,000
!!
!! Note that we also ignore hydrodynamic deformation of liquid droplets
!! as well as any breakup due to Rayleigh-Taylor instability.  
!!
!! This routine requires that vertical profiles of temperature <t>,
!! air density <rhoa>, and viscosity <rmu> are defined (i.e., initatm.f
!! must be called before this).  The vertical profile with ix = iy = 1
!! is used.
!!
!! We assume spherical particles -- call setupvf_std_shape() to use legacy
!! code from old Toon model for non-spherical effects -- use (better
!! yet, fix) at own risk.
!!
!! Added support for the particle radius being dependent on the relative
!! humidity according to the parameterizations of Gerber [1995] and
!! Fitzgerald [1975]. The fall velocity is then based upon the wet radius
!! rather than the dry radius. For particles that are not subject to
!! swelling, the wet and dry radii are the same. 
!!
!! @author  Chuck Bardeen, Pete Colarco from Andy Ackerman
!! @version Mar-2010 from Nov-2000 
subroutine setupvf_std(carma, cstate, j, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma    !! the carma object
  type(carmastate_type), intent(inout) :: cstate   !! the carma state object
  integer, intent(in)                  :: j        !! group index
  integer, intent(inout)               :: rc       !! return code, negative indicates failure

  ! Local declarations
  integer                 :: i, k
  real(kind=f)            :: x, y, cdrag
  real(kind=f)            :: rhoa_cgs, vg, rmfp, rkn, expon
                                   
  ! Define formats
  1 format(/,'Non-spherical particles specified for group ',i3, &
      ' (ishape=',i3,') but spheres assumed in I_FALLRTN_STD.', &
      ' Suggest using non-spherical code in I_FALLRTN_STD_SHAPE.')

  !  Warning message for non-spherical particles!
  if( ishape(j) .ne. 1 )then
    if (do_print) write(LUNOPRT,1) j, ishape(j)
  endif
  
  ! Loop over all atltitudes.
  do k = 1, NZ

    ! This is <rhoa> in cartesian coordinates (good old cgs units)
    rhoa_cgs = rhoa(k) / (xmet(k)*ymet(k)*zmet(k))

    ! <vg> is mean thermal velocity of air molecules [cm/s]
    vg = sqrt(8._f / PI * R_AIR * t(k))

    ! <rmfp> is mean free path of air molecules [cm]
    rmfp = 2._f * rmu(k) / (rhoa_cgs * vg)

    ! Loop over particle size bins.
    do i = 1,NBIN
    
      ! <rkn> is knudsen number
      rkn = rmfp / (r_wet(k,i,j) * rrat(i,j))

      ! <bpm> is the slip correction factor, the correction term for
      ! non-continuum effects.  Also used to calculate coagulation kernels
      ! and diffusion coefficients.
      expon = -.87_f / rkn
      expon = max(-POWMAX, expon)
      bpm(k,i,j) = 1._f + (1.246_f*rkn + 0.42_f*rkn*exp(expon))

      ! Stokes fall velocity and Reynolds' number
      vf(k,i,j) = (ONE * 2._f / 9._f) * rhop_wet(k,i,j) * r_wet(k,i,j)**2 * GRAV * bpm(k,i,j) / rmu(k) / rprat(i,j)
      re(k,i,j) = 2. * rhoa_cgs * r_wet(k,i,j) * rprat(i,j) * vf(k,i,j) / rmu(k)

      if (re(k,i,j) .ge. 1._f) then

        ! Correct drag coefficient for turbulence 
        x = log(re(k,i,j) / bpm(k,i,j))
        y = x*(0.83_f - 0.013_f*x)

        re(k,i,j) = exp(y) * bpm(k,i,j)

        if (re(k,i,j) .le. 1.e3_f) then

          ! drag coefficient from quadratic fit y(x) when Re < 1,000
          vf(k,i,j) = re(k,i,j) * rmu(k) / (2._f * r_wet(k,i,j) * rprat(i,j) * rhoa_cgs)
        else
        
          ! drag coefficient = 0.45 independent of Reynolds number when Re > 1,000
          cdrag = 0.45_f 
          vf(k,i,j) = bpm(k,i,j) * &
                      sqrt( 8._f * rhop_wet(k,i,j) * r_wet(k,i,j) * GRAV / &
                      (3._f * cdrag * rhoa_cgs * rprat(i,j)**2.) )
        endif
      endif
    enddo      ! <i=1,NBIN>
  enddo      ! <k=1,NZ>

  ! Return to caller with particle fall velocities evaluated.
  return
end
