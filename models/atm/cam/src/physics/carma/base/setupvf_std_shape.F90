! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine evaluates particle fall velocities, vf(k) [cm s^-1]
!!  and reynolds' numbers based on fall velocities, re(j,i,k) [dimensionless].
!!  indices correspond to vertical level <k>, bin index <i>, and aerosol
!!  group <j>.
!!
!!  Non-spherical particles are treated through shape factors <ishape>
!!  and <eshape>.
!!  
!!  General method is to first use Stokes' flow to estimate fall
!!!  velocity, then calculate reynolds' number, then use "y function" 
!!  (defined in Pruppacher and Klett) to reevaluate reynolds' number,
!!  from which the fall velocity is finally obtained.
!!
!!  This routine requires that vertical profiles of temperature <t>,
!!  air density <rhoa>, and viscosity <rmu> are defined (i.e., initatm.f
!!  must be called before this).
!!
!! @author  Chuck Bardeen, Pete Colarco from Andy Ackerman
!! @version Mar-2010 from Oct-1995


subroutine setupvf_std_shape(carma, cstate, j, rc)

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
  integer, intent(in)                  :: j       !! group index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                 :: i, k, ilast
  real(kind=f)            :: x, y
  real(kind=f)            :: rhoa_cgs, vg, rmfp, rkn, expon
  real(kind=f)            :: f1, f2, f3, ex, exx, exy, xcc, xa, bxx, r_shape, rfix, b0, bb1, bb2, bb3, z

  !  Define formats
  1 format('setupvfall::ERROR - ishape != 1, no fall velocity algorithm')


  ! First evaluate factors that depend upon particle shape (used in correction
  ! factor <bpm> below).  
  if (ishape(j) .eq. I_SPHERE) then

    ! Spheres
    f1 = 1.0_f
    f2 = 1.0_f

  else if (ishape(j) .eq. I_HEXAGON) then

    ! Hexagons: taken from Turco et al (Planet. Space Sci. Rev. 30, 1147-1181, 1982)
    ! with diffuse reflection of air molecules assumed 
    f2 = (PI / 9._f / tan(PI / 6._f))**(ONE/3._f) * eshape(j)**(ONE/6._f)

  else if (ishape(j) .eq. I_CYLINDER)then

    ! Spheroids: also from Turco et al. [1982]
    f2 = (2._f / 3._f)**(ONE/3._f) * eshape(j)**(ONE/6._f)
  endif

  ! (following statement yields <f3> = 1.0 for <eshape> = I_SPHERE)
  f3 = 1.39_f / sqrt((1.14_f + 0.25_f / eshape(j)) * (0.89_f + eshape(j) / 2._f))
  f2 = f2 * f3

  if (eshape(j) .gt. 1._f) then

    ! For Stokes regime there is no separate data for hexagonal plates or columns,
    ! so we use prolate spheroids.  This is from Fuchs' book.
    exx = eshape(j)**2 - 1._f
    exy = sqrt(exx)
    xcc = 1.333_f * exx / ((2._f * eshape(j)**2 - 1._f) * log(eshape(j) + exy) / exy-eshape(j))
    xa  = 2.666_f * exx / ((2._f * eshape(j)**2 - 3._f) * log(eshape(j) + exy) / exy+eshape(j))
!      f1  = eshape(j)**(-ONE/3._f) * (xcc + 2._f*xa) / 3._f
    f1  = eshape(j)**(-2._f/3._f) * (xcc + 2._f*xa) / 3._f

  elseif (eshape(j) .lt. 1._f) then

    ! Use oblate spheroids for disks (eshape < 1.).  Also from Fuchs' book.
    bxx = 1._f / eshape(j)
    exx = bxx**2 - 1._f
    exy = sqrt(exx)
    xcc = 1.333_f * exx / (bxx * (bxx**2 - 2._f) * atan(exy) / exy + bxx)
    xa  = 2.666_f * exx / (bxx * (3._f * bxx**2 - 2._f) * atan(exy) / exy - bxx)
    f1  = bxx**(ONE/3._f) * (xcc + 2._f * xa) / 3._f
  endif


  ! Loop over column with ixy = 1
  do k = 1,NZ

    ! This is <rhoa> in cartesian coordinates (good old cgs units)
    rhoa_cgs = rhoa(k) / (xmet(k)*ymet(k)*zmet(k))

    ! <vg> is mean thermal velocity of air molecules [cm/s]
    vg = sqrt(8._f / PI * R_AIR * t(k))

    ! <rmfp> is mean free path of air molecules [cm]
    rmfp = 2._f * rmu(k) / (rhoa_cgs * vg)

    ! Loop over particle size bins.
    do i = 1,NBIN

      ! <r_shape> is radius of particle used to calculate <re>.
      if (ishape(j) .eq. I_SPHERE) then
        r_shape = r_wet(k,i,j)
      else if (ishape(j) .eq. I_HEXAGON) then
        r_shape = r_wet(k,i,j) * 0.8456_f * eshape(j)**(-ONE/3._f)
      else if(ishape(j) .eq. I_CYLINDER) then
!        r_shape = r_wet(k,i,j) * eshape(j)**(-ONE/3._f)

        ! Shouldn't this have a factor related to being a cylinder vs a
        ! sphere in addition to the aspect ratio factor?
        r_shape = r_wet(k,i,j) * 0.8736_f * eshape(j)**(-ONE/3._f)
      endif

      ! <rkn> is knudsen number
      rkn = rmfp / r_wet(k,i,j)

      ! <bpm> is the slip correction factor, the correction term for
      ! non-continuum effects.  Also used to calculate coagulation kernels
      ! and diffusion coefficients.
      expon = -.87_f / rkn
      expon = max(-POWMAX, expon)
      bpm(k,i,j) = 1._f + f1*f2*(1.246_f*rkn + 0.42_f*rkn*exp(expon))

      ! These are first guesses for fall velocity and Reynolds' number, 
      ! valid for Reynolds' number < 0.01
      !
      ! This is "regime 1" in Pruppacher and Klett (chap. 10, pg 416).
      vf(k,i,j) = (2._f / 9._f) * rhop_wet(k,i,j) *(r_wet(k,i,j)**2) * GRAV * bpm(k,i,j) / (f1 * rmu(k))
      re(k,i,j) = 2._f * rhoa_cgs * r_shape * vf(k,i,j) / rmu(k)


      ! <rfix> is used in drag coefficient.
      rfix = vol(i,j) * rhop_wet(k,i,j) * GRAV * rhoa_cgs / rmu(k)**2

      if ((re(k,i,j) .ge. 0.01_f) .and. (re(k,i,j) .le. 300._f)) then

        ! This is "regime 2" in Pruppacher and Klett (chap. 10, pg 417).
        ! 
        ! NOTE: This sphere case is not the same solution used when
        ! interpolating other shape factors. This seems potentially inconsistent.
        if (ishape(j) .eq. I_SPHERE) then

          x = log(24._f * re(k,i,j) / bpm(k,i,j))
          y = -0.3318657e1_f + x * 0.992696_f - x**2 * 0.153193e-2_f - &
              x**3 * 0.987059e-3_f - x**4 * 0.578878e-3_f + &
              x**5 * 0.855176E-04_f - x**6 * 0.327815E-05_f 
              
          if (y .lt. -675._f) y = -675._f
          if (y .ge.  741._f) y =  741._f
          
          re(k,i,j) = exp(y) * bpm(k,i,j)

        else if (eshape(j) .le. 1._f) then

          ! P&K pg. 427
          if (ishape(j) .eq. I_HEXAGON) then
            x = log10(16._f * rfix / (3._f * sqrt(3._f)))
          else if (ishape(j) .eq. I_CYLINDER) then
            x = log10(8._f * rfix / PI)
          endif

          if (eshape(j) .le. 0.2_f) then
            
            ! P&K, page 424-427
            b0 = -1.33_f
            bb1 = 1.0217_f
            bb2 = -0.049018_f
            bb3 = 0.0_f
          else if (eshape(j) .le. 0.5_f) then

            ! NOTE: This interpolation/extrapolation method is
            ! not discussed in P&K; although, the solution for
            ! eshape = 0.5 is shown. Does this really work?
            ex = (eshape(j) - 0.2_f) / 0.3_f
            b0 = -1.33_f      + ex * (-1.3247_f   + 1.33_f)
            bb1 = 1.0217_f    + ex * (1.0396_f    - 1.0217_f)
            bb2 = -0.049018_f + ex * (-0.047556_f + 0.049018_f)
            bb3 =               ex * (-0.002327_f)
          else
          
            ! Extrapolating to cylinder cases on 436.
            ex = (eshape(j) - 0.5_f) / 0.5_f
            b0 = -1.3247_f    + ex * (-1.310_f    + 1.3247_f)
            bb1 = 1.0396_f    + ex * (0.98968_f   - 1.0396_f)
            bb2 = -0.047556_f + ex * (-0.042379_f + 0.047556_f)
            bb3 = -0.002327_f + ex * (              0.002327_f)
          endif

          y = b0 + x * bb1 + x**2 * bb2 + x**3 * bb3
          re(k,i,j) = 10._f**y * bpm(k,i,j)

        else if (eshape(j) .gt. 1._f) then
          ! Why is this so different from the oblate case?
          ! This seems wrong.
!            x = log10(2._f * rfix / eshape(j))
          if (ishape(j) .eq. I_CYLINDER) then
            x = log10(8._f * rfix / PI)
          endif
          
          ! P&K pg 430
          if( eshape(j) .le. 2._f )then
            ex = eshape(j) - 1._f
            b0 = -1.310_f     + ex * (-1.11812_f  + 1.310_f)
            bb1 = 0.98968_f   + ex * (0.97084_f   - 0.98968_f)
            bb2 = -0.042379_f + ex * (-0.058810_f + 0.042379_f)
            bb3 =               ex * (0.002159_f)
          else if (eshape(j) .le. 10._f) then
            ex = (eshape(j) - 2._f) / 8.0_f
            b0 = -1.11812_f   + ex * (-0.90629_f  + 1.11812_f)
            bb1 = 0.97084_f   + ex * (0.90412_f   - 0.97084_f)
            bb2 = -0.058810_f + ex * (-0.059312_f + 0.058810_f)
            bb3 = 0.002159_f  + ex * (0.0029941_f - 0.002159_f)
          else
          
            ! This is interpolating to a solution for an infinite
            ! cylinder, so it may not be the greatest estimate.
            ex = 10._f / eshape(j)
            b0 = -0.79888_f   + ex * (-0.90629_f  + 0.79888_f)
            bb1 = 0.80817_f   + ex * (0.90412_f   - 0.80817_f)
            bb2 = -0.030528_f + ex * (-0.059312_f + 0.030528_f)
            bb3 =               ex * (0.0029941_f)
          endif
          
          y = b0 + x * bb1 + x**2 * bb2 + x**3 * bb3
          re(k,i,j) = 10._f**y * bpm(k,i,j)

        endif

        !  Adjust <vf> for non-sphericicity.
        vf(k,i,j) = re(k,i,j) * rmu(k) / (2._f * r_shape * rhoa_cgs)

      endif

      if (re(k,i,j) .gt. 300._f) then

        ! This is "regime 3" in Pruppacher and Klett (chap. 10, pg 418).

!          if ((do_print) .and. (ishape(j) .ne. I_SPHERE)) write(LUNOPRT,1)
!          if ((do_print) .and. (ishape(j) .ne. I_SPHERE)) write(LUNOPRT,*) "setupvfall:", j, i, k, re(k,i,j)
!          rc = RC_ERROR
!          return
        
        z  = ((1.e6_f * rhoa_cgs**2) / (GRAV * rhop_wet(k,i,j) * rmu(k)**4))**(ONE/6._f)
        b0 = (24._f * vf(k,i,j) * rmu(k)) / 100._f
        x  = log(z * b0)
        y  = -5.00015_f + x * (5.23778_f   - x * (2.04914_f - x * (0.475294_f - &
                x * (0.0542819_f - x * 0.00238449_f))))
        
        if (y .lt. -675._f) y = -675.0_f
        if (y .ge.  741._f) y =  741.0_f

        re(k,i,j) = z * exp(y) * bpm(k,i,j)
        vf(k,i,j) = re(k,i,j) * rmu(k) / ( 2._f * r_wet(k,i,j) * rhoa_cgs)

        ! Values should not decrease with diameter, but instead should
        ! reach a limiting velocity that is independent of size (see
        ! Figure 10-25 of Pruppacher and Klett, 1997)
        ilast = max(1,i-1)
        if ((vf(k,i,j) .lt.  vf(k,ilast,j)) .or. (re(k,i,j) .gt. 4000._f)) then
          vf(k,i,j) = vf(k,ilast,j) 
        endif
      endif
    enddo    ! <i=1,NBIN>
  enddo      ! <k=1,NZ>

  ! Return to caller with particle fall velocities evaluated.
  return
end
