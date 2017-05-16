#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

! To remove the effect of moisture comment the following
!#define _COMPUTE_MOISTURE_
module physics_mod
  
  ! =======================
  use kinds,              only : real_kind
  ! =======================
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, Rd_on_Rv, Cp, Cpd_on_Cpv, cpwater_vapor
  ! =======================
  use physical_constants, only : rearth,p0
  ! =======================
  use dimensions_mod, only : np, nlev
  ! =======================
  implicit none
  
  private
  
  public :: Saturation_Vapor_Pressure
  public :: Specific_Humidity
  public :: Saturation_Specific_Humidity
  public :: Relative_Humidity
  public :: Vapor_Pressure
  public :: Mixing_Ratio
  public :: Prim_Condense
  public :: getsurfpress
  public :: Temp2PotTemp
  public :: Virtual_Temperature
  public :: Virtual_Specific_Heat
  public :: kappastar  

 interface Virtual_Temperature
    module procedure Virtual_Temperature1d
    module procedure Virtual_Temperature3d
 end interface


contains
  
  !===========================
  !
  ! For help or information:
  ! 
  ! Amik St-Cyr
  ! 
  ! e-mail: amik@ucar.edu
  !
  !===========================
 
  !================================
  ! For reference see Emanuel 1994 
  !================================
  
  function Virtual_Temperature1d(Tin,rin) result(Tv)
    
    real (kind=real_kind),intent(in) :: Tin
    real (kind=real_kind),intent(in) :: rin
    real (kind=real_kind)            :: Tv

!    Tv = Tin*(1_real_kind + rin/Rd_on_Rv)/(1_real_kind + rin)

    Tv = Tin*(1_real_kind + (Rwater_vapor/Rgas - 1.0_real_kind)*rin)


  end function Virtual_Temperature1d

  function Virtual_Temperature3d(T,Q) result(T_v)
    real (kind=real_kind),intent(in) :: T(np,np,nlev)
    real (kind=real_kind),intent(in) :: Q(np,np,nlev)
    real (kind=real_kind) :: T_v(np,np,nlev)
    integer :: i, j, k

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
    do k=1,nlev
       do j=1,np
          do i=1,np
             T_v(i,j,k) = Virtual_Temperature1d(T(i,j,k), Q(i,j,k))
          end do
       end do
    end do
  end function Virtual_Temperature3d

  function kappastar(Q) result(ks)
    real(kind=real_kind), intent(in) :: Q(np,np,nlev)
    real(kind=real_kind) :: ks(np,np,nlev)
    integer i,j,k

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j)
#endif
    do k=1,nlev
       do j=1,np
          do i=1,np
             ks(i,j,k) =  Rgas/Virtual_Specific_Heat(Q(i,j,k))
          end do
       end do
    end do
  end function kappastar


  function Virtual_Specific_Heat(rin) result(Cp_star)
    
    real (kind=real_kind),intent(in) :: rin
    real (kind=real_kind)            :: Cp_star

!    Cp_star = Cp*(1_real_kind + rin/Cpd_on_Cpv)/(1_real_kind + rin)
 
    Cp_star = Cp*(1.0_real_kind + (Cpwater_vapor/Cp - 1.0_real_kind)*rin)
   
  end function Virtual_Specific_Heat


  !=================================================
  ! Approx. Solution to the Clausius-Clapeyron eqn.
  !=================================================
  function Saturation_Vapor_Pressure(Tin) result(estar)

    real (kind=real_kind),intent(in) :: Tin
    real (kind=real_kind)            :: estar

    ! 4.4.13 p. 116 Emanuel
#ifdef USE_SAT_TABLE
    logical, save :: firstcall=.true.
    integer, parameter :: len=(10000)
    real(kind=real_kind), parameter :: Tmin=0.01, Tmax=300.0
    real(kind=real_kind) :: Table(len), t, dt
    integer :: i

    dt = (Tmax-Tmin)/real(len-1)
    if(firstcall) then
       do i=1,len
          t = Tmin + dt*real(i-1)
          Table(i)  = exp(53.67957_real_kind - 6743.769_real_kind/T &
	- 4.8451_real_kind * log(T))
       enddo
    else
       i = (tin-tmin)/dt + 1
       estar = table(i)
    endif

#else
    
    estar = exp(53.67957_real_kind - 6743.769_real_kind/abs(Tin) &
	- 4.8451_real_kind * log(ABS(Tin)))
#endif
    ! convert to code units
    if (p0 <  2000 ) then
       ! code is using mb, do nothing
    else
       ! code is using Pa.  convert from mb to Pa
       estar = estar*100
    endif

  end function Saturation_Vapor_Pressure

  function Specific_Humidity(r) result(q)
    
    real (kind=real_kind),intent(in) :: r
    real (kind=real_kind)            :: q
    
    q = r/( r + 1_real_kind )

  end  function Specific_Humidity
  
  function Saturation_Specific_Humidity(p,T) result(qstar)
 
    real (kind=real_kind),intent(in) :: p,T
    real (kind=real_kind)            :: qstar,estar

    estar = Saturation_Vapor_Pressure(T)
    qstar = Rd_on_Rv * estar / (p - estar * (1._real_kind - Rd_on_Rv))
  end function Saturation_Specific_Humidity
  !==================================
  ! Fraction of humid air in dry air 
  !==================================
  function Relative_Humidity(e,T) result(h)
    
    real (kind=real_kind),intent(in) :: e,T
    real (kind=real_kind)            :: estar
    real (kind=real_kind)            :: h
    
    estar = Saturation_Vapor_Pressure(T)

    h=e/estar

  end function Relative_Humidity

  !=================================
  ! Partial pressure of water vapor
  !=================================
  function Vapor_Pressure(r,p) result(e_out)

    real (kind=real_kind),intent(in) :: r,p
    real (kind=real_kind)            :: e_out

    ! 4.1.2 p. 108 Emanuel

    ! With p -> dry pressure
    e_out = abs(r*p)/Rd_on_Rv

  end function Vapor_Pressure

  !==============================================
  ! Mass of water vapor per unit mass of dry air
  !==============================================
  function Mixing_Ratio(ein,p) result(r)
    
    real (kind=real_kind),intent(in) :: ein,p
        
    real (kind=real_kind)            :: r
    
    ! 4.1.2 p. 108 Emanuel
    ! p is supposed DRY
    
    r = Rd_on_Rv*abs(ein/p)

  end function Mixing_Ratio

  !===============================================
  ! This function creates rain.
  ! If e > e_saturation then
  ! e = e_saturation.
  ! The latent heat is not computed here for now.
  !===============================================
  subroutine Prim_Condense(r,Tin,pin)
    
    real (kind=real_kind),intent(inout)  :: r
    real (kind=real_kind),intent(in)     :: Tin,pin

    real (kind=real_kind)                :: estar,e_vapor
    real*8                               :: st,et

    ! returns pressure in mb or Pa?  
    estar   = Saturation_Vapor_Pressure(Tin)     

    ! returns pressure in same units as pin. 
    e_vapor = Vapor_Pressure(r,pin)

    if(e_vapor/estar>1_real_kind)then
       e_vapor = estar
       r       = Mixing_Ratio(e_vapor,pin)
    endif

  end subroutine Prim_Condense

  function Temp2PotTemp(pr3d,t3d) result(pt3d)
    real (kind=real_kind),intent(in) :: pr3d(np,np,nlev),t3d(np,np,nlev)
    real (kind=real_kind)            :: pt3d(np,np,nlev)
    integer:: i,j,k    
    real (kind=real_kind):: pp

    !
    ! dry
    !    
    do k=1,nlev
       do j=1,np
          do i=1,np
             pp = (pr3d(i,j,k) + pr3d(i,j,k+1))*0.5D0
             pt3d(i,j,k)=  t3d(i,j,k)*(p0/pp)**kappa
          enddo
       enddo
    enddo
  end function Temp2PotTemp

  function Exner_function(pr3d) result(exner)
    real (kind=real_kind),intent(in) :: pr3d(np,np,nlev)
    real (kind=real_kind)            :: exner(np,np,nlev)
    integer:: i,j,k    
    real (kind=real_kind):: pp
    !
    ! dry
    !    
    do k=1,nlev
       do j=1,np
          do i=1,np
             pp = (pr3d(i,j,k) + pr3d(i,j,k+1))*0.5D0
             exner(i,j,k)=  (pp/p0)**kappa
          enddo
       enddo
    enddo
  end function Exner_function

  function getsurfpress(lnps) result (press)
    real (kind=real_kind) :: press(np,np)
    real (kind=real_kind) :: lnps(np,np)

    press(:,:) = 0.0
    
  end function getsurfpress
  
     
  end module physics_mod
