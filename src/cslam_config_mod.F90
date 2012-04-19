#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
 
!MODULE CSLAM_CONFIG_MOD------------------------------------------------CE-for CSLAM!
! AUTHOR: CHRISTOPH ERATH, 11.June 2011                                             !
! This module sets parameters to run CSLAM with HOMME                               ! 
!                                                                                   !
! set in configure.ac NC                                                            !
! nc  : nc*nc is the number of finite volume cells on each element (we have         !
!       6*ne*ne*nc*nc cells on the sphere)                                          !
! nhe ... halo/depth of the extended element (CFL number), now only nhe=1 is tested !
!         this number sets where we have to calculate the reconstruction in the halo!
! nhe ... halo/depth for the tracer values, only cubic reconstruction is supported  !
!         now, therefore nhc=nhe+3                                                  !
! ngpc... number of Gausspoints for the CSLAM integral approximation                !
! _CSLAM_ON_GAUSS: is not supported by CSLAM now                                    !
!                  we can define a fvm mesh build on the Gaussian points NP         ! 
!                  then nc is defined by NP (nc=np-1)                               !
!                                                                                   !
!-----------------------------------------------------------------------------------!  
module cslam_config_mod
  use kinds, only : real_kind
  
  implicit none
  private
  !nhc has to be at least 4 (because of the cubic reconstruction)
  integer, parameter, public :: nhe=1     !number/depth of the extended element (CFL number)
  integer, parameter, public :: nhc=nhe+3 !number/depth of halos (ghost cells) for cslam tracers
  integer, parameter, public :: ngpc=4    !number of Gausspoints for the CSLAM integral approximation

#ifdef _CSLAM_ON_GAUSS
  integer, parameter, public :: nc=NP-1
#else
  integer, parameter, public :: nc=NC       !number (one side) of fvm cells on an element
#endif

end module cslam_config_mod
