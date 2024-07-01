! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_manufactured_solution
!
!> \brief Computes tendency terms for a manufactured solution
!> \author Steven Brus, Carolyn Begeman
!> \date   June 2023
!> \details
!>  This module contains the routines for computing the thickness and
!>  normal velocity tendencies for a manufactured solution case. See
!>  Bishnu et al. 2022 (https://doi.org/10.22541/essoar.167100170.03833124/v1)
!
!-----------------------------------------------------------------------

module ocn_manufactured_solution 

   use mpas_constants
   use mpas_timer
   use mpas_timekeeping
   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_manufactured_solution_tend_thick, &
             ocn_manufactured_solution_tend_vel, &
             ocn_manufactured_solution_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND) :: kx, ky
   real (kind=RKIND) :: ang_freq
   real (kind=RKIND) :: eta0
   real (kind=RKIND) :: H0

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_manufactured_solution_tend_thick 
!
!> \brief  Computes manufactured solution thickness tendency
!> \author Steven Brus, Carolyn Begeman
!> \date   June 2023 
!> \details
!>  This routine computes the thickness tendency for the manufactured solution
!
!-----------------------------------------------------------------------

   subroutine ocn_manufactured_solution_tend_thick(tend, err)!{{{

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: thickness tendency 

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: iCell, kmin, kmax, k
      real (kind=RKIND) :: phase, time

      ! End preamble
      !-------------
      ! Begin code

      if (.not. config_use_manufactured_solution) return

      time = daysSinceStartOfSim*86400.0_RKIND

      do iCell = 1,nCellsOwned

         kmin = minLevelCell(iCell)
         kmax = maxLevelCell(iCell)

         phase = kx*xCell(iCell) + ky*yCell(iCell) - ang_freq*time
         do k = kmin, kmax
            tend(k,iCell) = tend(k,iCell) + eta0*(-H0*(kx + ky)*sin(phase) &
                                                  - ang_freq*cos(phase) &
                                                  + eta0*(kx + ky)*cos(2.0_RKIND*phase))
         enddo

      enddo

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_manufactured_solution_tend_thick!}}}

!***********************************************************************
!
!  routine ocn_manufactured_solution_tend_vel 
!
!> \brief  Computes manufactured solution velocity tendency 
!> \author Steven Brus, Carolyn Begeman
!> \date   June 2023 
!> \details
!>  This routine computes the velocity tendency for the manufactured solution
!
!-----------------------------------------------------------------------
   subroutine ocn_manufactured_solution_tend_vel(tend, err)!{{{

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency 

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: iEdge, kmin, kmax, k
      real (kind=RKIND) :: phase, u, v, time

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      if (.not. config_use_manufactured_solution) return

      time = daysSinceStartOfSim*86400.0_RKIND

      do iEdge = 1, nEdgesOwned

         kmin = minLevelEdgeBot(iEdge)
         kmax = maxLevelEdgeTop(iEdge)

         phase = kx*xEdge(iEdge) + ky*yEdge(iEdge) - ang_freq*time

         do k = kmin, kmax

            u = eta0*((-fEdge(iEdge) + gravity*kx)*cos(phase) &
                       + ang_freq*sin(phase) &
                       - 0.5_RKIND*eta0*(kx + ky)*sin(2.0_RKIND*(phase)))
            v = eta0*((fEdge(iEdge) + gravity*ky)*cos(phase) &
                      + ang_freq*sin(phase) &
                      - 0.5_RKIND*eta0*(kx + ky)*sin(2.0_RKIND*(phase))) 

            tend(k,iEdge) = tend(k,iEdge) + u*cos(angleEdge(iEdge)) + v*sin(angleEdge(iEdge))
         enddo

      enddo

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_manufactured_solution_tend_vel!}}}

!***********************************************************************
!
!  routine ocn_manufactured_solution_init
!
!> \brief  Initialize the manufactured solution tendencies 
!> \author Steven Brus, Carolyn Begeman
!> \date   June 2023 
!> \details
!>  This routine initializes constants related to the manufactured
!>  solution tendencies
!
!-----------------------------------------------------------------------
    subroutine ocn_manufactured_solution_init(domain, err)!{{{

      type (domain_type), intent(inout) :: domain      
      integer, intent(out) :: err !< Output: Error flag

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: verticalMeshPool
      real (kind=RKIND), dimension(:,:), pointer :: restingThickness

      if (.not. config_use_manufactured_solution) return

      kx = 2.0_RKIND*pi / config_manufactured_solution_wavelength_x
      ky = 2.0_RKIND*pi / config_manufactured_solution_wavelength_y

      eta0 = config_manufactured_solution_amplitude


      ! This test case assumes that the restingThickness is horizontally uniform
      ! and that only one vertical level is used so only one set of indices is
      ! used here.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)

         call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

         H0 = restingThickness(1,1)

         block => block % next
      enddo

      ang_freq = sqrt(H0*gravity * (kx**2+ky**2))

      err = 0

    end subroutine ocn_manufactured_solution_init!}}}

!***********************************************************************

end module ocn_manufactured_solution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
