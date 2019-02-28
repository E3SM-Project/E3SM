! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_sediment_transport
!
!> \brief MPAS ocean analysis mode member: sediment_transport
!> \author Zhendong Cao and Phillip J. Wolfram
!> \date   2019/03/07
!> \details
!>  MPAS ocean analysis mode member: sediment_transport
!-----------------------------------------------------------------------

module ocn_sediment_transport

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines

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

   public :: ocn_init_sediment_transport, &
             ocn_compute_sediment_transport, &
             ocn_restart_sediment_transport, &
             ocn_finalize_sediment_transport

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_sediment_transport
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Zhendong Cao and Phillip J. Wolfram
!> \date    2019/03/07
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_sediment_transport(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_init_sediment_transport!}}}

!***********************************************************************
!
!  routine ocn_compute_sediment_transport
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Zhendong Cao and Phillip J. Wolfram
!> \date    2019/03/07
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_sediment_transport(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: sedimentTransportAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      real (kind=RKIND), pointer :: salpha, SD50, rho0, rhoS, Vh, Dv,Ws
      real (kind=RKIND), pointer :: tau_ce, Erate, poro, tau_cd, Cd
  !!    real (kind=RKIND), dimension(:,:) :: tau_bot
      real (kind=RKIND), dimension(:,:), pointer :: velX, velY, velZ
      real (kind=RKIND), dimension(:), pointer :: posX, posY
      real (kind=RKIND), dimension(:), pointer :: sedFluxVAx, sedFluxVAy, sedFluxBx, sedFluxBy
      real (kind=RKIND), dimension(:), pointer :: ero_flux, dep_flux
      real (kind=RKIND), dimension(:), pointer :: bedld_x, bedld_y
      logical, pointer :: on_a_sphere, use_lat_lon_coords
      logical, pointer :: bedload, suspended
      integer, pointer :: nCells, nVertLevels, nCellsSolve
      integer k, iCell, i
      real ratio1, rho_R, ND50, Umag
      real cff1, cff2, cff3
      real phicw, w_asym   !! parameters related to waves, assigned ZEROs for just now
      real tau_tide, tau_wave, tau_mean, theta_mean, theta_wave, theta_ce
      real phi_x1, phi_x2, phi_x, phi_y
      real, parameter :: g = 9.81_RKIND
      real, parameter :: eps = 1.0E-14_RKIND

      err = 0

      dminfo = domain % dminfo
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_alpha',salpha)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_grain_size',SD50)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_drag_coefficient',Cd)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_grain_porosity',poro)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_erate',Erate)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_tau_ce',tau_ce)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_tau_cd',tau_cd)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_water_density',rho0)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_grain_density',rhoS)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_kinematic_viscosity',Vh)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_vertical_diffusion_coefficient',Dv)

      block => domain % blocklist
      do while (associated(block))
         ! get dimensions
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'sedimentTransportAM', sedimentTransportAMPool)

         call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
         call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_use_lat_lon_coords', use_lat_lon_coords)
         call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_bedload', bedload)
         call mpas_pool_get_config(ocnConfigs, 'config_AM_sedimentTransport_suspended', suspended)

         if (.not. on_a_sphere) then
            use_lat_lon_coords = .false.
         end if

         if (use_lat_lon_coords) then
            call mpas_pool_get_array(meshPool, 'lonCell', posX)
            call mpas_pool_get_array(meshPool, 'latCell', posY)
            call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velX)
            call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velY)
         else
            call mpas_pool_get_array(meshPool, 'xCell', posX)
            call mpas_pool_get_array(meshPool, 'yCell', posY)
            call mpas_pool_get_array(diagnosticsPool, 'velocityX', velX)
            call mpas_pool_get_array(diagnosticsPool, 'velocityY', velY)
         end if 

         call mpas_pool_get_array(sedimentTransportAMPool, 'sedimentFallVelocity', Ws)
         call mpas_pool_get_array(sedimentTransportAMPool, 'sedimentErosionFlux', ero_flux)
         call mpas_pool_get_array(sedimentTransportAMPool, 'sedimentDepositionFlux', dep_flux)
         call mpas_pool_get_array(sedimentTransportAMPool, 'sedimentFluxVAX', sedFluxVAx)
         call mpas_pool_get_array(sedimentTransportAMPool, 'sedimentFluxVAY', sedFluxVAy)
         call mpas_pool_get_array(sedimentTransportAMPool, 'sedimentFluxBX', sedFluxBx)
         call mpas_pool_get_array(sedimentTransportAMPool, 'sedimentFluxBY', sedFluxBy)
         call mpas_pool_get_array(sedimentTransportAMPool, 'sedimentBedloadX', bedld_x)
         call mpas_pool_get_array(sedimentTransportAMPool, 'sedimentBedloadY', bedld_y)

     ! compute sedimentFallVelocity
         ! sieve diameter (SD50) to nominal diameter (ND50), Raudkivi [1990]
         ND50 = SD50*1.1_RKIND
         ! sediment/water density comparsion  -> specific graivity
         rho_R = rhoS/rho0-1
         ! sedimentFallVelocity formula GP technique, Goldstein & Coco [2013]
         Ws = (37.8_RKIND*ND50*rho_R*(1+100.0_RKIND*ND50))/   &
              (0.383_RKIND+1E4*rho_R*Vh+1E2*ND50*rho_R**2)
        ! compute ratio1 in sedFlux = ratio1*U^3, Grawe et al. [2014]
        ratio1 = salpha*Dv/Ws**2
        do iCell = 1,nCellsSolve
            sedFluxVAx(iCell) = ratio1*(sum(velX(:,iCell))/float(nVertLevels))**3.0_RKIND
            sedFluxVAy(iCell) = ratio1*(sum(velY(:,iCell))/float(nVertLevels))**3.0_RKIND
            sedFluxBx(iCell) = ratio1*velX(size(velX,1),iCell)**3.0_RKIND
            sedFluxBy(iCell) = ratio1*velY(size(velY,1),iCell)**3.0_RKIND
        end do

    ! compute bedload and suspended transport
        cff1 = rho_R*g*SD50
        cff2 = sqrt(rho_R*g*SD50)*SD50*rhoS
        !wave-related parameters
        w_asym = 0.0_RKIND
        phicw = 0.0_RKIND

        do iCell = 1,nCellsSolve
            Umag = sqrt(velX(size(velX,1),iCell)**2 + velY(size(velY,1),iCell)**2)
            tau_tide = rho0*Cd*Umag**2
            tau_wave = 0.0_RKIND
            tau_mean = tau_tide*(1.0_RKIND+1.2_RKIND*(tau_wave/(tau_wave+tau_tide))**1.5_RKIND)

        ! compute suspended transport <--- to be continued
            if (suspended) then
               if (tau_mean .gt. tau_ce ) then
                 ero_flux(iCell) = (1.0_RKIND-poro)*Erate*(tau_mean/tau_ce -1.0_RKIND)
               else
                 ero_flux(iCell) = 0.0_RKIND
               endif
            end if

        ! compute bedload transport, using formulae by Soulsby and Damgaard [2005]
            if (bedload) then
               theta_mean = tau_mean/cff1
               theta_wave = tau_wave/cff1
               theta_ce = tau_ce/cff1

               cff3 = 0.5*(1.0_RKIND+SIGN(1.0_RKIND,theta_mean/theta_ce-1.0_RKIND))

               phi_x1 = 12.0_RKIND*sqrt(theta_mean)*MAX(theta_mean-theta_ce,0.0_RKIND)
               phi_x2 = 12.0_RKIND*(0.9534_RKIND+0.1907_RKIND*COS(2.0_RKIND*phicw))* &
                        sqrt(theta_wave)*theta_mean +                                       &
                        12.0_RKIND*(0.229_RKIND*w_asym*theta_wave**1.5_RKIND*COS(phicw))

               if (ABS(phi_x2) .gt. phi_x1) then
                 phi_x = phi_x2
               else
                 phi_x = phi_x1
               end if

               bedld_x(iCell) = phi_x*cff2*cff3

               phi_y = 12.0_RKIND*0.1907_RKIND*theta_wave**2* &
                       (theta_mean*SIN(2.0_RKIND*phicw)+1.2_RKIND*w_asym*theta_wave*SIN(phicw)) &
                       /(theta_wave**1.5_RKIND + 1.5_RKIND*theta_mean**1.5_RKIND+eps)
               bedld_y = phi_y*cff2*cff3

            end if
        end do

         block => block % next
      end do

   end subroutine ocn_compute_sediment_transport!}}}

!***********************************************************************
!
!  routine ocn_restart_sediment_transport
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Zhendong Cao and Phillip J. Wolfram
!> \date    2019/03/07
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_sediment_transport(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_sediment_transport!}}}

!***********************************************************************
!
!  routine ocn_finalize_sediment_transport
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Zhendong Cao and Phillip J. Wolfram
!> \date    2019/03/07
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_sediment_transport(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_finalize_sediment_transport!}}}

end module ocn_sediment_transport

! vim: foldmethod=marker
