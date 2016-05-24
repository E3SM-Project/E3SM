!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Simple 2d periodic domain flow model
!c  with stretched vertical grid and 1D BL model
!c  also, with surface flux spatial variability allowed
!c
!c  Author: W. Grabowski (grabow@ncar.ucar.edu) 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc
PROGRAM moist_convec 

  use crcp_mod,       only :  crcp, crcp_input_parameters

  implicit none
  integer, parameter :: nx=101, nz=61, nxz=nx*nz

#ifdef _BGL
#include <mpif.h>
      integer :: info
#endif


  integer :: istep=0

  real, parameter :: dt=15.0,dx=2.0e3,dz=400.0
  real, parameter :: t_strt=0.,t_end=1800.
  real, parameter :: sst=300.48, ssqv= 25.43e-3

  integer, PARAMETER :: npin = 23


  real :: zin(npin), press(npin),temp(npin),vap(npin)
  real :: uu(npin), vv(npin)
  real :: dux(npin),duy(npin),dth(npin),dqv(npin),dqc(npin),dqr(npin) 

  !ccc overwrite with st calculated for deep atmosphere:                  
  ! taken from global model...      

  real, parameter :: st = 1.52E-05

  DATA zin / npin * 0. / 
  DATA press / 1008.00, 991.25, 945.50, 893.79, 836.06, 772.82,     &
       705.22, 635.05, 564.48, 495.73, 430.71, 370.78, 316.72, 268.82,   &
       226.98, 190.82, 159.87, 133.55, 111.29, 92.56, 52.31, 22.08, 9.32 &
       /                                                                 
  DATA temp / 25.26, 24.13, 21.04, 18.66, 16.50, 13.41, 9.06, 3.73, &
       - 1.51, - 6.97, - 14.09, - 22.44, - 30.57, - 39.60, - 48.69,      &
       - 57.40, - 65.21, - 72.58, - 76.71, - 74.98, - 74.98, - 74.98,    &
       - 74.98 /                                                         
  DATA vap / 0.178E+02, 0.172E+02, 0.156E+02, 0.134E+02, 0.111E+02, &
       0.888E+01, 0.631E+01, 0.487E+01, 0.396E+01, 0.200E+01, 0.984E+00, &
       0.806E+00, 0.370E+00, 0.135E+00, 0.599E-01, 0.258E-01, 0.123E-01, &
       0.582E-02, 0.367E-02, 0.589E-02, 0.104E-02, 0.247E-02, 0.585E-02 /

  DATA uu / npin * 4. / 
  DATA vv / npin * 0. / 



!     call opngks
!     call gsclip(0)
#ifdef _BGL
      call mpi_init(info)
#endif
#ifdef TESTMODE
      call tbeg('main')
#endif
! Read the crcp Namelist (crcp.nl in the run directory)

      call crcp_input_parameters()

      call crcp(press,temp,zin,vap,uu,vv,st,dt,dx,dz,npin,dux,duy, &
         dth, dqv, dqc, dqr, sst, ssqv, t_strt, t_end, istep, nx, nz)  

!cc    finished...

!     call clsgks
#ifdef TESTMODE
      call tend('main')
#endif
      call tprt()
#ifdef _BGL
      call mpi_finalize(info)
#endif
      stop
    end PROGRAM moist_convec 
