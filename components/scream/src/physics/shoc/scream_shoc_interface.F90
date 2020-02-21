#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module scream_shoc_interface_mod

  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool,C_NULL_CHAR, c_float

  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
  integer,parameter,public :: rtype = c_double ! 8 byte real, compatible with c type double
#else
# define c_real c_float
  integer,parameter,public :: rtype = c_float ! 4 byte real, compatible with c type float
#endif

  public :: shoc_init_f90
  public :: shoc_main_f90
  public :: shoc_finalize_f90

  real   :: test

  integer(kind=c_int) :: pcols = 32
  integer(kind=c_int) :: pver  = 72
  integer(kind=c_int) :: qsize = 9

  real(kind=c_real) :: cpair  =    1004.64000000000
  real(kind=c_real) :: rair   =    287.042311365049
  real(kind=c_real) :: rh2o   =    461.504639820160
  real(kind=c_real) :: rhoh2o =    1000.00000000000
  real(kind=c_real) :: mwh2o  =    18.0160000000000
  real(kind=c_real) :: mwdry  =    28.9660000000000
  real(kind=c_real) :: gravit =    9.80616000000000
  real(kind=c_real) :: latvap =    2501000.00000000
  real(kind=c_real) :: latice =    333700.000000000
  real(kind=c_real) :: cpliq  =    4188.00000000000
  real(kind=c_real) :: tmelt  =    273.150000000000
  real(kind=c_real) :: pi     =    3.14159265358979
  real(kind=c_real) :: karman =    0.40000000000000
  real(kind=c_real) :: zvir   =    0.60779307282415

contains

  !====================================================================!
  subroutine shoc_init_f90 (q) bind(c)

    use shoc,                   only: shoc_init,r8
 
    real(kind=c_real), intent(inout) :: q(pcols,pver,9) ! State array  kg/kg
    
    real(kind=c_real) :: pref_mid(pcols,pver)           ! pressure at midlevel hPa
    integer(kind=c_int) :: its, ite, kts, kte

    kts     = 1
    kte     = pver

    do k = kte,kts,-1 
      pref_mid(:,k)    = 1e3_rtype - (1e3_rtype-0.1)/real(pver)!state%pmid(:,:)
    end do

    call shoc_init(& 
          integer(pver),&
          real(gravit,kind=r8),&
          real(rair,kind=r8),  &
          real(rh2o,kind=r8),  &
          real(cpair,kind=r8), &
	  real(zvir,kind=r8),  &
          real(latvap,kind=r8),&
	  real(latice,kind=r8),&
	  real(karman,kind=r8),&
	  pref_mid,            &
	  integer(kte),&
	  integer(kts))   

    q(:,:,:) = 0.0_rtype
    q(:,:,1) = 1.0e-5_rtype!state%q(:,:,1)
    q(:,:,2) = 1.0e-6_rtype!state%q(:,:,ixcldliq)
    q(:,:,3) = 1.0e-7_rtype!state%q(:,:,ixcldice)
    q(:,:,4) = 1.0e6_rtype!state%q(:,:,ixnumliq)
    q(:,:,5) = 1.0e5_rtype!state%q(:,:,ixnumice)
    q(:,:,6) = 1.0e-5_rtype!state%q(:,:,ixrain)
    q(:,:,7) = 1.0e5_rtype!state%q(:,:,ixnumrain)
    q(:,:,8) = 1.0e-8_rtype!state%q(:,:,ixcldrim) !Aaron, changed ixqirim to ixcldrim to match Kai's code
    q(:,:,9) = 1.0e4_rtype!state%q(:,:,ixrimvol)
     

    test = 0.0
    print '(a15,f16.8,e16.8,i8,i8)', 'SHOC init = ', test, sum(q(1,:,1)), pcols, pver

  end subroutine shoc_init_f90
  !====================================================================!
  subroutine shoc_main_f90 (dtime,q,FQ,qdp) bind(c)

    use shoc,           only: shoc_main

!    real, intent(in) :: q(pcols,pver,9) ! Tracer mass concentrations from SCREAM      kg/kg
    real(kind=c_real), intent(in)    :: dtime ! Timestep 
    real(kind=c_real), intent(inout) :: q(pcols,pver,qsize) ! Tracer mass concentrations from SCREAM kg/kg
    real(kind=c_real), intent(inout) :: FQ(pcols,4,pver) ! Tracer mass concentrations from SCREAM kg/kg
    real(kind=c_real), intent(in)    :: qdp(pcols,2,4,pver) ! Tracer mass concentrations from SCREAM kg/kg

    real(kind=c_real) :: cldliq(pcols,pver)     !cloud liquid water mixing ratio        kg/kg
    real(kind=c_real) :: numliq(pcols,pver)     !cloud liquid water drop concentraiton  #/kg
    real(kind=c_real) :: rain(pcols,pver)       !rain water mixing ratio                kg/kg
    real(kind=c_real) :: numrain(pcols,pver)    !rain water number concentration        #/kg
    real(kind=c_real) :: qv(pcols,pver)         !water vapor mixing ratio               kg/kg
    real(kind=c_real) :: ice(pcols,pver)        !total ice water mixing ratio           kg/kg
    real(kind=c_real) :: qirim(pcols,pver)      !rime ice mixing ratio                  kg/kg
    real(kind=c_real) :: numice(pcols,pver)     !total ice crystal number concentration #/kg
    real(kind=c_real) :: rimvol(pcols,pver)     !rime volume mixing ratio               m3/kg
    real(kind=c_real) :: pres(pcols,pver)       !pressure at midlevel                   hPa

    real(kind=c_real) :: inv_cp 

    integer(kind=c_int) :: k, ncol
    integer :: i

    integer(kind=c_int) :: its, ite, kts, kte

    real(kind=c_real) :: qtest

    inv_cp = 1.0_rtype/cpair

    qtest = sum(q)
    ncol = pcols


    its     = 1
    ite     = ncol
    kts     = 1
    kte     = pver

    do i = its,ite
      do k = kts,kte
        qv(i,k)      = q(i,k,1) !1.0e-4_rtype!state%q(:,:,1)
        cldliq(i,k)  = q(i,k,2) !1.0e-6_rtype!state%q(:,:,ixcldliq)
        ice(i,k)     = q(i,k,3) !1.0e-7_rtype!state%q(:,:,ixcldice)
        numliq(i,k)  = q(i,k,4) !1.0e6_rtype!state%q(:,:,ixnumliq)
        numice(i,k)  = q(i,k,5) !1.0e5_rtype!state%q(:,:,ixnumice)
        rain(i,k)    = q(i,k,6) !1.0e-5_rtype!state%q(:,:,ixrain)
        numrain(i,k) = q(i,k,7) !1.0e5_rtype!state%q(:,:,ixnumrain)
        qirim(i,k)   = q(i,k,8) !1.0e-8_rtype!state%q(:,:,ixcldrim) !Aaron, changed ixqirim to ixcldrim to match Kai's code
        rimvol(i,k)  = q(i,k,9) !1.0e4_rtype!state%q(:,:,ixrimvol)
      end do
    end do 

    do k = kte,kts,-1 
      pres(:,k)    = 1e3_rtype - (1e3_rtype-0.1)/real(pver)!state%pmid(:,:)
    end do

    do i = its,ite
      do k = kts,kte
        q(i,k,1) = qv(i,k)*1.01 
        q(i,k,2) = cldliq(i,k)  
        q(i,k,3) = ice(i,k)     
        q(i,k,4) = numliq(i,k)  
        q(i,k,5) = numice(i,k)  
        q(i,k,6) = rain(i,k)    
        q(i,k,7) = numrain(i,k) 
        q(i,k,8) = qirim(i,k)   
        q(i,k,9) = rimvol(i,k)
      end do
    end do  

    test = test + dtime
    FQ(1,1,1) = 9e9
    print '(a15,f16.8,5e16.8)', 'SHOC run = ', test, qtest, sum(q), sum(qv), sum(FQ(:,:,:)), sum(qdp)

  end subroutine shoc_main_f90
  !====================================================================!
  subroutine shoc_finalize_f90 () bind(c)

    test = 42.
    print '(a15,f16.8)', 'SHOC final = ', test

  end subroutine shoc_finalize_f90
  !====================================================================!

end module scream_shoc_interface_mod
