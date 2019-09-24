#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module scream_shoc_interface_mod

  use iso_c_binding

  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
  integer,parameter,public :: rtype = c_double ! 8 byte real, compatible with c type double
#else
# define c_real c_float
  integer,parameter,public :: rtype = c_float ! 4 byte real, compatible with c type float
#endif

  integer(kind=c_int) :: pcols = 32
  integer(kind=c_int) :: pver  = 72
  integer(kind=c_int) :: pveri = 73
  integer(kind=c_int) :: qsize = 9

  real   :: test
contains

  subroutine hi(pcols, pver, pveri, nqtracers, dtime, arr) bind(c)
    integer, value, intent(in) :: pcols, pver, pveri, nqtracers
    real(kind=c_real), value, intent(in) :: dtime
    real(kind=c_real), dimension(pcols,pver), intent(inout) :: arr
    integer :: i, j

    print *, 'hi', pcols, pver, pveri, nqtracers, dtime
    do j = 1, pver
       do i = 1, pcols
          arr(i,j) = 10*i + j
       end do
    end do
  end subroutine hi

  subroutine shoc_c_init() bind(c)
    use shoc, only: shoc_init
    implicit none

    real(kind=c_real) :: gravit =    9.80616000000000
    real(kind=c_real) :: rair   =    287.042311365049
    real(kind=c_real) :: rh2o   =    461.504639820160
    real(kind=c_real) :: cpair  =    1004.64000000000
    real(kind=c_real) :: latvap =    2501000.00000000

    call shoc_init(gravit, rair, rh2o, cpair, latvap)

    test = 10.0_rtype
  end subroutine shoc_c_init

  subroutine shoc_c_main( &
     q,dtime &          ! Input
     ) bind(c)

    use shoc, only: shoc_main
    implicit none

    real(kind=c_real), intent(inout) :: q(pcols,pver,qsize) ! Tracer mass concentrations from SCREAM kg/kg
    real(kind=c_real), intent(in) :: dtime
    real(kind=c_real), dimension(pcols) :: &
         host_dx, host_dy, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc
    real(kind=c_real), dimension(pcols,pver) :: &
         zt_grid, pres, pdel, thv, cldliq, w_field
    real(kind=c_real), dimension(pcols,pveri) :: &
         zi_grid
    real(kind=c_real), dimension(pcols,qsize) :: &
         wtracer_sfc

    real(kind=c_real), dimension(pcols,pver) :: &
         tke, thetal, qw, u_wind, v_wind, wthv_sec, tk, tkh
    real(kind=c_real), dimension(pcols,pver,qsize) :: &
         qtracers

    real(kind=c_real), dimension(pcols,pver) :: &
         shoc_cldfrac, shoc_ql, shoc_mix, w_sec, wqls_sec, brunt, isotropy
    real(kind=c_real), dimension(pcols,pveri) :: &
         thl_sec, qw_sec, qwthl_sec, wthl_sec, wqw_sec, wtke_sec, uw_sec, vw_sec, w3

     host_dx(:)         = 0.0_rtype 
     host_dy(:)         = 0.0_rtype
     thv(:,:)           = 0.0_rtype
     cldliq(:,:)        = 0.0_rtype      
     zt_grid(:,:)       = 0.0_rtype 
     zi_grid(:,:)       = 0.0_rtype 
     pres(:,:)          = 0.0_rtype 
     pdel(:,:)          = 0.0_rtype          
     wthl_sfc (:)       = 0.0_rtype
     wqw_sfc(:)         = 0.0_rtype
     uw_sfc(:)          = 0.0_rtype
     vw_sfc(:)          = 0.0_rtype 
     wtracer_sfc(:,:)   = 0.0_rtype
     w_field(:,:)       = 0.0_rtype
     tke(:,:)           = 0.0_rtype
     thetal(:,:)        = 0.0_rtype
     qw(:,:)            = 0.0_rtype            
     u_wind(:,:)        = 0.0_rtype
     v_wind(:,:)        = 0.0_rtype
     qtracers(:,:,:)    = 0.0_rtype
     wthv_sec(:,:)      = 0.0_rtype
     tkh(:,:)           = 0.0_rtype
     tk(:,:)            = 0.0_rtype                    
     shoc_cldfrac(:,:)  = 0.0_rtype
     shoc_ql(:,:)       = 0.0_rtype          
     shoc_mix(:,:)      = 0.0_rtype 
     isotropy(:,:)      = 0.0_rtype                 
     w_sec(:,:)         = 0.0_rtype
     thl_sec(:,:)       = 0.0_rtype
     qw_sec(:,:)        = 0.0_rtype 
     qwthl_sec(:,:)     = 0.0_rtype
     wthl_sec(:,:)      = 0.0_rtype 
     wqw_sec(:,:)       = 0.0_rtype
     wtke_sec(:,:)      = 0.0_rtype
     uw_sec(:,:)        = 0.0_rtype
     vw_sec(:,:)        = 0.0_rtype
     w3(:,:)            = 0.0_rtype              
     wqls_sec(:,:)      = 0.0_rtype
     brunt(:,:)         = 0.0_rtype                  

    call shoc_main( &
         pcols, pver, pveri, dtime,&          ! Input
         host_dx, host_dy,thv, cldliq,&       ! Input
         zt_grid, zi_grid, pres, pdel,&       ! Input
         wthl_sfc, wqw_sfc, uw_sfc, vw_sfc,&  ! Input
         wtracer_sfc, qsize, w_field,& ! Input     
         tke, thetal, qw, &                   ! Input/Output
         u_wind, v_wind, qtracers,&           ! Input/Output
         wthv_sec,tkh,tk,&                    ! Input/Output
         shoc_cldfrac, shoc_ql,&              ! Output
         shoc_mix, isotropy,&                 ! Output (diagnostic)
         w_sec, thl_sec, qw_sec, qwthl_sec,&  ! Output (diagnostic)
         wthl_sec, wqw_sec, wtke_sec,&        ! Output (diagnostic)
         uw_sec, vw_sec, w3,&                 ! Output (diagnostic)    
         wqls_sec, brunt)
 
    test = test + dtime
    print '(a15,f16.8,e16.8)', 'SHOC Run = ', test, dtime
  end subroutine shoc_c_main

  subroutine shoc_c_main_ntimes(dtime,ntimes) bind(c)

    use shoc, only: shoc_main
    implicit none

    integer(kind=c_int), value, intent(in) :: ntimes
    real(kind=c_real), value, intent(in) :: dtime
    real(kind=c_real), dimension(pcols) :: &
         host_dx, host_dy, wthl_sfc, wqw_sfc
    real(kind=c_real), dimension(pcols) :: uw_sfc, vw_sfc
    real(kind=c_real), dimension(pcols,pver) :: &
         zt_grid, pres, pdel, thv, cldliq, w_field
    real(kind=c_real), dimension(pcols,pveri) :: &
         zi_grid
    real(kind=c_real), dimension(pcols,qsize) :: &
         wtracer_sfc

    real(kind=c_real), dimension(pcols,pver) :: &
         tke, thetal, qw, u_wind, v_wind, wthv_sec, tk, tkh
    real(kind=c_real), dimension(pcols,pver,qsize) :: &
         qtracers

    real(kind=c_real), dimension(pcols,pver) :: &
         shoc_cldfrac, shoc_ql, shoc_mix, w_sec, wqls_sec, brunt, isotropy
    real(kind=c_real), dimension(pcols,pveri) :: &
         thl_sec, qw_sec, qwthl_sec, wthl_sec, wqw_sec, wtke_sec, uw_sec, vw_sec, w3

    integer :: i, c
    real(8), parameter :: ustar2 = 0.28d0**2
    real(8) :: speed

    do i = 1, ntimes
       do c = 1, pcols
          speed = sqrt(u_wind(c,1)**2 + v_wind(c,1)**2)
          uw_sfc(c) = -ustar2*(u_wind(c,1)/speed)
          vw_sfc(c) = -ustar2*(v_wind(c,1)/speed)
       end do

       call shoc_main( &
            pcols, pver, pveri, dtime,&          ! Input
            host_dx, host_dy,thv, cldliq,&       ! Input
            zt_grid, zi_grid, pres, pdel,&       ! Input
            wthl_sfc, wqw_sfc, uw_sfc, vw_sfc,&  ! Input
            wtracer_sfc, qsize, w_field,& ! Input     
            tke, thetal, qw, &                   ! Input/Output
            u_wind, v_wind, qtracers,&           ! Input/Output
            wthv_sec,tkh,tk,&                    ! Input/Output
            shoc_cldfrac, shoc_ql,&              ! Output
            shoc_mix, isotropy,&                 ! Output (diagnostic)
            w_sec, thl_sec, qw_sec, qwthl_sec,&  ! Output (diagnostic)
            wthl_sec, wqw_sec, wtke_sec,&        ! Output (diagnostic)
            uw_sec, vw_sec, w3,&                 ! Output (diagnostic)    
            wqls_sec, brunt)
    end do
  end subroutine shoc_c_main_ntimes

  subroutine shoc_c_finalize() bind(c)
    use shoc, only: shoc_finalize
    implicit none

    call shoc_finalize()

  end subroutine shoc_c_finalize
end module scream_shoc_interface_mod
