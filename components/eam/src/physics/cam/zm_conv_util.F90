module zm_conv_util
   !----------------------------------------------------------------------------
   ! Purpose: utility methods for ZM deep convection scheme
   !----------------------------------------------------------------------------
   use shr_kind_mod,     only: r8=>shr_kind_r8
   use cam_abortutils,   only: endrun
   use zm_conv_types,    only: zm_const_t

   public :: entropy       ! calculate entropy
   public :: ientropy      ! invert entropy equation to get temperature and saturated vapor mixing ratio
   public :: qsat_hPa      ! wrapper for qsat_water that translates between Pa and hPa

!===================================================================================================
contains
!===================================================================================================

real(r8) function entropy(TK, p, qtot, zm_const)
   !----------------------------------------------------------------------------
   ! Purpose: function to calculate entropy following:
   ! 
   !    Raymond, D. J., and A. M. Blyth, 1992: Extension of the Stochastic Mixing
   !       Model to Cumulonimbus Clouds. J. Atmos. Sci., 49, 1968–1983
   !----------------------------------------------------------------------------
   ! Arguments
   real(r8),         intent(in) :: TK        ! temperature              [K]
   real(r8),         intent(in) :: p         ! pressure                 [mb]
   real(r8),         intent(in) :: qtot      ! total water mixing ratio [kg/kg]
   type(zm_const_t), intent(in) :: zm_const  ! derived type to hold ZM constants
   !----------------------------------------------------------------------------
   ! Local variables
   real(r8) :: qv    ! water vapor mixing ratio
   real(r8) :: qst   ! saturated vapor mixing ratio
   real(r8) :: e     ! water vapor pressure
   real(r8) :: est   ! saturated vapor pressure
   real(r8) :: L     ! latent heat of vaporization
   real(r8), parameter :: pref = 1000._r8
   !----------------------------------------------------------------------------

   ! Calculate latent heat of vaporization - note T is converted to centigrade
   L = zm_const%latvap - (zm_const%cpliq - zm_const%cpwv)*(TK-zm_const%tfreez)

   ! Use saturation mixing ratio to partition qtot into vapor part only
   call qsat_hPa(TK, p, est, qst)
   qv = min(qtot,qst)
   e = qv*p / (zm_const%epsilo+qv)

   ! calculate entropy per unit mass of dry air - Eq. 1
   entropy = (  zm_const%cpair &
              + qtot*zm_const%cpliq)*log(TK/zm_const%tfreez) &
              - zm_const%rdair*log( (p-e)/pref &
             ) + L*qv/TK - qv*zm_const%rh2o*log(qv/qst)

end function entropy

!===================================================================================================

subroutine ientropy(rcall, s, p, qt, T, qst, Tfg, zm_const)
   !----------------------------------------------------------------------------
   ! Purpose: invert the entropy equation to return temperature and saturated
   ! vapor mixing ratio following Richard Brent's method::
   ! 
   !    Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives.
   !       Englewood Cliffs, NJ: Prentice-Hall, 1973.
   !----------------------------------------------------------------------------
   ! Arguments
   integer,          intent(in)  :: rcall   ! call index
   real(r8),         intent(in)  :: s       ! entropy                           [J/kg]
   real(r8),         intent(in)  :: p       ! pressure                          [mb]
   real(r8),         intent(in)  :: qt      ! total water mixing ratio          [kg/kg]
   real(r8),         intent(in)  :: Tfg     ! input temperature for first guess [K]
   real(r8),         intent(out) :: qst     ! saturation vapor mixing ratio     [kg/kg]
   real(r8),         intent(out) :: T       ! temperature                       [k]
   type(zm_const_t), intent(in)  :: zm_const  ! derived type to hold ZM constants
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i           ! loop iterator
   logical  :: converged   ! flag for convergence
   real(r8) :: est         ! saturation vapor pressure
   real(r8) :: this_lat    ! local latitude
   real(r8) :: this_lon    ! local logitude
   real(r8) :: tolerance
   real(r8) :: a, b, c, d, ebr, fa, fb, fc, pbr, qbr, rbr, sbr, xm
   integer,  parameter :: LOOPMAX   = 100      ! Max number of iteration loops
   real(r8), parameter :: tol_coeff = 0.001_r8 ! tolerance coeficient
   real(r8), parameter :: tol_eps   = 3.e-8_r8 ! small value for tolerance calculation
   !----------------------------------------------------------------------------
   ! initialize variables
   converged = .false.
   T = Tfg            ! first guess based on input temperature
   a = Tfg-10         ! low bracket
   b = Tfg+10         ! high bracket

   fa = entropy(a, p, qt, zm_const) - s
   fb = entropy(b, p, qt, zm_const) - s

   c = b
   fc = fb
   !----------------------------------------------------------------------------
   ! 
   converge: do i=0, LOOPMAX

      if ((fb > 0.0_r8 .and. fc > 0.0_r8) .or. &
          (fb < 0.0_r8 .and. fc < 0.0_r8)) then
         c   = a
         d   = b-a
         fc  = fa
         ebr = d
      end if

      if (abs(fc) < abs(fb)) then
         a  = b
         b  = c
         c  = a
         fa = fb
         fb = fc
         fc = fa
      end if

      tolerance = 2.0_r8*tol_eps*abs(b) + 0.5_r8*tol_coeff
      xm = 0.5_r8*(c-b)
      
      converged = (abs(xm) <= tolerance .or. fb == 0.0_r8)
      if (converged) exit converge

      if (abs(ebr) >= tolerance .and. abs(fa) > abs(fb)) then
         sbr=fb/fa
         if (a == c) then
            pbr = 2.0_r8*xm*sbr
            qbr = 1.0_r8-sbr
         else
            qbr = fa/fc
            rbr = fb/fc
            pbr = sbr*(2.0_r8*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1.0_r8))
            qbr = (qbr-1.0_r8)*(rbr-1.0_r8)*(sbr-1.0_r8)
         end if
         if (pbr > 0.0_r8) qbr=-qbr
         pbr=abs(pbr)
         if (2.0_r8*pbr  <  min(3.0_r8*xm*qbr-abs(tolerance*qbr),abs(ebr*qbr))) then
            ebr = d
            d = pbr/qbr
         else
            d = xm
            ebr = d
         end if
      else
         d = xm
         ebr = d
      end if
      a = b
      fa = fb
      b = b + merge( d, sign(tolerance,xm), abs(d)>tolerance )

      fb = entropy(b, p, qt, zm_const) - s

   end do converge

   T = b
   call qsat_hPa(T, p, est, qst)

   if (.not. converged) then
      write(iulog,*) '*** ZM_CONV: IENTROPY: Failed and about to exit, info follows ****'
      write(iulog,100) 'ZM_CONV: IENTROPY Details:', &
                       ' call#: ',rcall, &
                       ' P(mb): ',p, &
                       ' Tfg(K): ', Tfg, &
                       ' qt(g/kg): ', 1000._r8*qt, &
                       ' qst(g/kg): ', 1000._r8*qst,&
                       ' s(J/kg): ',s
     call endrun('**** ZM_CONV IENTROPY: Tmix did not converge ****')
  end if

100 format (A,I1,I4,I4,7(A,F6.2))

end subroutine ientropy

!===================================================================================================

elemental subroutine qsat_hPa(t, p, es, qm)
   !----------------------------------------------------------------------------
   ! Purpose: wrapper for qsat_water that translates between Pa and hPa
   ! qsat_water uses Pa internally, so pass in Pa and set es back to hPa after
   use wv_saturation, only: qsat_water
   !----------------------------------------------------------------------------
   ! Arguments
   real(r8), intent(in)  :: t   ! Temperature                  [K]
   real(r8), intent(in)  :: p   ! Pressure                     [hPa]
   real(r8), intent(out) :: es  ! Saturation vapor pressure    [hPa]
   real(r8), intent(out) :: qm  ! Saturation mass mixing ratio [kg/kg] (vapor mass over dry mass)
   !----------------------------------------------------------------------------

   call qsat_water(t, p*100._r8, es, qm)

   es = es*0.01_r8

end subroutine qsat_hPa

!===================================================================================================

end module zm_conv_util
