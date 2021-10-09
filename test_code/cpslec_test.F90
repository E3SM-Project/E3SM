subroutine cpslec_test (ncol, pmid, phis, ps, T, Q, psl, gravit, rair, diag_psl_opt)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Hybrid coord version:  Compute sea level pressure for a latitude line
! 
! Method: 
! CCM2 hybrid coord version using ECMWF formulation
! Algorithm: See section 3.3 in ykfpos46t1r1 "FULL-POS IN 
! THE CYCLE 46T1R1 OF ARPEGE/IFS". 2D dynamical variables 
! which need extrapolations 
!
! Author: Shixuan Zhang 
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver
  use cam_abortutils,  only: endrun
  use shr_log_mod, only:errMsg => shr_log_errMsg

  implicit none

!-----------------------------Arguments---------------------------------
  integer , intent(in) :: ncol             ! longitude dimension

  real(r8), intent(in) :: pmid(pcols,pver) ! Atmospheric pressure (pascals)
  real(r8), intent(in) :: phis(pcols)      ! Surface geopotential (m**2/sec**2)
  real(r8), intent(in) :: ps(pcols)        ! Surface pressure (pascals)
  real(r8), intent(in) :: T(pcols,pver)    ! Vertical slice of temperature (top to bot)
  real(r8), intent(in) :: Q(pcols,pver)    ! Vertical slice of temperature (top to bot)
  real(r8), intent(in) :: gravit           ! Gravitational acceleration
  real(r8), intent(in) :: rair             ! gas constant for dry air

  real(r8), intent(out) :: psl(pcols)      ! Sea level pressures (pascals)
  integer, intent(in) :: diag_psl_opt      ! option for different method for psl 

!-----------------------------------------------------------------------

!-----------------------------Parameters--------------------------------
  real(r8), parameter :: clapse   = 6.5e-3_r8 ! Temperature lapse rate (K/m)
  integer, parameter  :: nlextrap = 1  ! the reference level for the
                                       ! extropolation 
!-----------------------------------------------------------------------

!-----------------------------Local Variables---------------------------
  integer i              ! Loop index
  real(r8) alpha         ! Temperature lapse rate in terms of pressure ratio (unitless)
  real(r8) Tstar         ! Computed surface temperature
  real(r8) TT0           ! Computed temperature at sea-level
  real(r8) alph          ! Power to raise P/Ps to get rate of increase of T with pressure
  real(r8) beta          ! alpha*phis/(R*T) term used in approximation of PSL
  real(r8) Tref          ! Temperature (K) at reference level (pver - nlextrap)
  real(r8) Pref          ! Pressure (hPa) at reference level (pver - nlextrap)
  real(r8) Qref          ! Humidity(kg/kg) at bottom model level  
  real(r8) Tv            ! Virtual temperaure 
  character(len=2000) err_str

!-----------------------------------------------------------------------
!
  alpha = rair*clapse/gravit

  do i=1,ncol

     !Apply topography adjustment on temperature 
     !The reference level nlextrap indicate that the extropolation will happen 
     !below the pver-nlextrap layer 
     Tref  = T(i,pver-nlextrap)
     Qref  = Q(i,pver)
     Pref  = pmid(i,pver-nlextrap)

     !calculate the surface temperature 
     Tstar = Tref*(1._r8+alpha*(ps(i)/pref - 1._r8))          ! pg 8 eq (1)
     !calculate the temperature at mean sea level 
     TT0   = Tstar + clapse*phis(i)/gravit                      ! pg 11 

     !calculate the pressure at mean seal level
     if ( abs(phis(i)/gravit) < 1.e-3_r8 )then

        psl(i)=ps(i)

     else
 
        if ( Tstar<=290.5_r8 .and. TT0>290.5_r8 ) then           ! pg 11 eq (7)
           alph = rair / phis(i) * ( 290.5_r8 - Tstar )  
        else if (Tstar > 290.5_r8 .and. TT0 > 290.5_r8) then     ! pg 11 
           alph=0._r8
           Tstar= 0.5_r8 * (290.5_r8 + Tstar)  
        else  
           alph=alpha  
           if (Tstar < 255._r8) then  
              Tstar= 0.5_r8 * (255._r8 + Tstar)                  ! pg 11
           endif
        endif

       ! Verify that version numbers are valid.
       if ( diag_psl_opt == 1) then 

         !Computes sea level pressure using the hypsometric equation.
         Tv = Tstar * (1.0_r8 + 0.608_r8 * Qref)  
         psl(i)= ps(i)*exp( phis(i) / (rair * Tv) )

       else if ( diag_psl_opt ==2) then  

         !Computes sea level pressure using the ECMWF's algorithm
         !Note the calculate is the same as the default with cpslec ()
         !Except that the threshold and the layer for reference are changed 
         beta = phis(i)/(rair*Tstar)
         psl(i)=ps(i)*exp( beta*(1._r8-alph*beta/2._r8+((alph*beta)**2)/3._r8)) ! pg 11 eq (8)
       
       else

         write(err_str,*) 'cpslec_test:  diag_psl_opt is invalid', errmsg(__FILE__, __LINE__)
         call endrun(err_str)

       end if         
       
     end if
  enddo

  return
end subroutine cpslec_test

