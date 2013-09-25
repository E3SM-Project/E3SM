! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routine for the split explicit
!>  time integration scheme
!
!-----------------------------------------------------------------------


module ocn_time_integration_split

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_tendency
   use ocn_diagnostics

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average

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

   public :: ocn_time_integrator_split

   type (timer_node), pointer :: timer_main, timer_prep, timer_bcl_vel, timer_btr_vel, timer_diagnostic_update, timer_implicit_vmix, &
                                 timer_halo_diagnostic, timer_halo_normalBarotropicVelocity, timer_halo_ssh, timer_halo_f, timer_halo_thickness, & 
                                 timer_halo_tracers, timer_halo_normalBaroclinicVelocity

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine integrates a single time step (dt) using a
!>  split explicit time integrator.
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integrator_split(domain, dt)!{{{
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Advance model state forward in time by the specified time step using 
    !   Split_Explicit timestepping scheme
    !
    ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
    !                 plus mesh meta-data
    ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
    !                  model state advanced forward in time by dt seconds
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      type (dm_info) :: dminfo
      integer :: iCell, i,k,j, iEdge, cell1, cell2, split_explicit_step, split, &
                 eoe, oldBtrSubcycleTime, newBtrSubcycleTime, uPerpTime, BtrCorIter, &
                 n_bcl_iter(config_n_ts_iter), stage1_tend_time, startIndex, endIndex
      type (block_type), pointer :: block
      real (kind=RKIND) :: normalThicknessFluxSum, thicknessSum, flux, sshEdge, hEdge1, &
                 CoriolisTerm, uCorr, temp, temp_h, coef, barotropicThicknessFlux_coeff, sshCell1, sshCell2
      integer :: num_tracers, ucorr_coef, err
      real (kind=RKIND), dimension(:,:), pointer :: &
                 u, h, layerThicknessEdge, vertViscTopOfEdge, vertDiffTopOfCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer, dimension(:), pointer :: & 
                 maxLevelCell, maxLevelEdgeTop
      real (kind=RKIND), dimension(:), allocatable:: uTemp
      real (kind=RKIND), dimension(:,:), allocatable:: tracersTemp

      call mpas_timer_start("se timestep", .false., timer_main)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !  Prep variables before first iteration
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpas_timer_start("se prep", .false., timer_prep)
      block => domain % blocklist
      do while (associated(block))

         ! Initialize * variables that are used to compute baroclinic tendencies below.
         do iEdge=1,block % mesh % nEdges
            do k=1,block % mesh % nVertLevels !maxLevelEdgeTop % array(iEdge)

               ! The baroclinic velocity needs be recomputed at the beginning of a 
               ! timestep because the implicit vertical mixing is conducted on the
               ! total u.  We keep normalBarotropicVelocity from the previous timestep.
               ! Note that normalBaroclinicVelocity may now include a barotropic component, because the 
               ! weights layerThickness have changed.  That is OK, because the barotropicForcing variable
               ! subtracts out the barotropic component from the baroclinic.
                 block % state % time_levs(1) % state % normalBaroclinicVelocity % array(k,iEdge) &
               = block % state % time_levs(1) % state % normalVelocity    % array(k,iEdge) &
               - block % state % time_levs(1) % state % normalBarotropicVelocity % array(  iEdge)

                 block % state % time_levs(2) % state % normalVelocity % array(k,iEdge) &
               = block % state % time_levs(1) % state % normalVelocity % array(k,iEdge)

                 block % state % time_levs(2) % state % normalBaroclinicVelocity % array(k,iEdge) &
               = block % state % time_levs(1) % state % normalBaroclinicVelocity % array(k,iEdge)

                 block % diagnostics % layerThicknessEdge % array(k,iEdge) &
               = block % diagnostics % layerThicknessEdge % array(k,iEdge)

            end do 
         end do 

           block % state % time_levs(2) % state % ssh % array(:) &
         = block % state % time_levs(1) % state % ssh % array(:)

         do iCell=1,block % mesh % nCells  
            do k=1,block % mesh % maxLevelCell % array(iCell)

                 block % state % time_levs(2) % state % layerThickness % array(k,iCell) &
               = block % state % time_levs(1) % state % layerThickness % array(k,iCell)

                 block % state % time_levs(2) % state % tracers % array(:,k,iCell) & 
               = block % state % time_levs(1) % state % tracers % array(:,k,iCell) 

            end do
         end do

         if (config_use_freq_filtered_thickness) then

              block % state % time_levs(2) % state % highFreqThickness % array(:,:) &
            = block % state % time_levs(1) % state % highFreqThickness % array(:,:)

              block % state % time_levs(2) % state % lowFreqDivergence % array(:,:) &
            = block % state % time_levs(1) % state % lowFreqDivergence % array(:,:)

         endif

         block => block % next
      end do

      call mpas_timer_stop("se prep", timer_prep)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN large iteration loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      n_bcl_iter = config_n_bcl_iter_mid
      n_bcl_iter(1) = config_n_bcl_iter_beg
      n_bcl_iter(config_n_ts_iter) = config_n_bcl_iter_end

      do split_explicit_step = 1, config_n_ts_iter
         stage1_tend_time = min(split_explicit_step,2)

         ! ---  update halos for diagnostic variables
         call mpas_timer_start("se halo diag", .false., timer_halo_diagnostic)
         call mpas_dmpar_exch_halo_field(domain % blocklist % diagnostics % normalizedRelativeVorticityEdge)
         if (config_mom_del4 > 0.0) then
           call mpas_dmpar_exch_halo_field(domain % blocklist % diagnostics % divergence)
           call mpas_dmpar_exch_halo_field(domain % blocklist % diagnostics % relativeVorticity)
         end if
         call mpas_timer_stop("se halo diag", timer_halo_diagnostic)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 1: Baroclinic velocity (3D) prediction, explicit with long timestep
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (config_use_freq_filtered_thickness) then
            call mpas_timer_start("se freq-filtered-thick computations")
            block => domain % blocklist
            do while (associated(block))
               call ocn_tend_freq_filtered_thickness(block % tend, &
                  block % state % time_levs(stage1_tend_time) % state, block % diagnostics, block % mesh)
               block => block % next
            end do
            call mpas_timer_stop("se freq-filtered-thick computations")

            call mpas_timer_start("se freq-filtered-thick halo update")
            call mpas_dmpar_exch_halo_field(domain % blocklist % tend % highFreqThickness)
            call mpas_dmpar_exch_halo_field(domain % blocklist % tend % lowFreqDivergence)
            call mpas_timer_stop("se freq-filtered-thick halo update")

            block => domain % blocklist
            do while (associated(block))
               do iCell=1,block % mesh % nCells
                  do k=1,block % mesh % maxLevelCell % array(iCell)
                     ! this is h^{hf}_{n+1}
                        block % state % time_levs(2) % state % highFreqThickness % array(k,iCell) &
                      = block % state % time_levs(1) % state % highFreqThickness % array(k,iCell) &
                      + dt* block % tend % highFreqThickness % array(k,iCell) 
                  end do
               end do
               block => block % next
            end do

         endif


         ! compute velocity tendencies, T(u*,w*,p*)
         call mpas_timer_start("se bcl vel", .false., timer_bcl_vel)

         block => domain % blocklist
         do while (associated(block))

           ! compute vertTransportVelocityTop.  Use u (rather than uTransport) for momentum advection.
           ! Use the most recent time level available.
           call ocn_vert_transport_velocity_top(block % mesh, block % verticalMesh, &
              block % state % time_levs(1) % state % layerThickness % array, &
              block % diagnostics % layerThicknessEdge % array, &
              block % state % time_levs(stage1_tend_time) % state % normalVelocity % array, &
              block % state % time_levs(1) % state % ssh % array, &
              block % state % time_levs(2) % state % highFreqThickness % array, dt, &
              block % diagnostics % vertTransportVelocityTop % array, err)

            call ocn_tend_vel(block % tend, block % state % time_levs(stage1_tend_time) % state, block % diagnostics, block % mesh)

            block => block % next
         end do

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do j=1,n_bcl_iter(split_explicit_step)

            ! Use this G coefficient to avoid an if statement within the iEdge loop.
            if (trim(config_time_integrator) == 'unsplit_explicit') then
               split = 0
            elseif (trim(config_time_integrator) == 'split_explicit') then
               split = 1
            endif

            block => domain % blocklist
            do while (associated(block))
               allocate(uTemp(block % mesh % nVertLevels))

               ! Put f*normalBaroclinicVelocity^{perp} in uNew as a work variable
               call ocn_fuperp(block % state % time_levs(2) % state , block % mesh)

               do iEdge=1,block % mesh % nEdges
                  cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
                  cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

                  uTemp = 0.0  ! could put this after with uTemp(maxleveledgetop+1:nvertlevels)=0
                  do k=1,block % mesh % maxLevelEdgeTop % array(iEdge)

                     ! normalBaroclinicVelocityNew = normalBaroclinicVelocityOld + dt*(-f*normalBaroclinicVelocityPerp + T(u*,w*,p*) + g*grad(SSH*) )
                     ! Here uNew is a work variable containing -fEdge(iEdge)*normalBaroclinicVelocityPerp(k,iEdge)
                      uTemp(k) = block % state % time_levs(1) % state % normalBaroclinicVelocity % array(k,iEdge) &
                         + dt * (block % tend % normalVelocity % array (k,iEdge) &
                         + block % state % time_levs(2) % state % normalVelocity % array (k,iEdge) &  ! this is f*normalBaroclinicVelocity^{perp}
                         + split * gravity * (  block % state % time_levs(2) % state % ssh % array(cell2) &
                         - block % state % time_levs(2) % state % ssh % array(cell1) ) &
                          /block % mesh % dcEdge % array(iEdge) )
                  enddo

                  ! thicknessSum is initialized outside the loop because on land boundaries 
                  ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a 
                  ! nonzero value to avoid a NaN.
                  normalThicknessFluxSum = block % diagnostics % layerThicknessEdge % array(1,iEdge) * uTemp(1)
                  thicknessSum  = block % diagnostics % layerThicknessEdge % array(1,iEdge)

                  do k=2,block % mesh % maxLevelEdgeTop % array(iEdge)
                     normalThicknessFluxSum = normalThicknessFluxSum + block % diagnostics % layerThicknessEdge % array(k,iEdge) * uTemp(k)
                     thicknessSum  =  thicknessSum + block % diagnostics % layerThicknessEdge % array(k,iEdge)
                  enddo
                  block % diagnostics % barotropicForcing % array(iEdge) = split*normalThicknessFluxSum/thicknessSum/dt


                  do k=1,block % mesh % maxLevelEdgeTop % array(iEdge)
                     ! These two steps are together here:
                     !{\bf u}'_{k,n+1} = {\bf u}'_{k,n} - \Delta t {\overline {\bf G}}
                     !{\bf u}'_{k,n+1/2} = \frac{1}{2}\left({\bf u}^{'}_{k,n} +{\bf u}'_{k,n+1}\right) 
                     ! so that normalBaroclinicVelocityNew is at time n+1/2
                       block % state % time_levs(2) % state % normalBaroclinicVelocity % array(k,iEdge) &
                     = 0.5*( &
                       block % state % time_levs(1) % state % normalBaroclinicVelocity % array(k,iEdge) &
                     + uTemp(k) - dt * block % diagnostics % barotropicForcing % array(iEdge))

                  enddo
 
               enddo ! iEdge

               deallocate(uTemp)

               block => block % next
            end do

            call mpas_timer_start("se halo normalBaroclinicVelocity", .false., timer_halo_normalBaroclinicVelocity)
            call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % normalBaroclinicVelocity)
            call mpas_timer_stop("se halo normalBaroclinicVelocity", timer_halo_normalBaroclinicVelocity)

         end do  ! do j=1,config_n_bcl_iter

         call mpas_timer_stop("se bcl vel", timer_bcl_vel)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! END baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 2: Barotropic velocity (2D) prediction, explicitly subcycled
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("se btr vel", .false., timer_btr_vel)

         oldBtrSubcycleTime = 1
         newBtrSubcycleTime = 2

         if (trim(config_time_integrator) == 'unsplit_explicit') then

            block => domain % blocklist
            do while (associated(block))

               ! For Split_Explicit unsplit, simply set normalBarotropicVelocityNew=0, normalBarotropicVelocitySubcycle=0, and uNew=normalBaroclinicVelocityNew
               block % state % time_levs(2) % state % normalBarotropicVelocity % array(:) = 0.0

               block % state % time_levs(2) % state % normalVelocity % array(:,:)  = block % state % time_levs(2) % state % normalBaroclinicVelocity % array(:,:) 

               do iEdge=1,block % mesh % nEdges
                  do k=1,block % mesh % nVertLevels

                     ! uTranport = normalBaroclinicVelocity + uBolus 
                     ! This is u used in advective terms for layerThickness and tracers 
                     ! in tendency calls in stage 3.
                     block % diagnostics % uTransport % array(k,iEdge) &
                           = block % mesh % edgeMask % array(k,iEdge) &
                           *(  block % state % time_levs(2) % state % normalBaroclinicVelocity       % array(k,iEdge) &
                           + block % diagnostics % uBolusGM   % array(k,iEdge) )

                  enddo
               end do  ! iEdge
   
               block => block % next
            end do  ! block

         elseif (trim(config_time_integrator) == 'split_explicit') then

            ! Initialize variables for barotropic subcycling
            block => domain % blocklist
            do while (associated(block))

               if (config_filter_btr_mode) then
                  block % diagnostics % barotropicForcing % array(:) = 0.0
               endif

               do iCell=1,block % mesh % nCells
                  ! sshSubcycleOld = sshOld  
                    block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(iCell) &
                  = block % state % time_levs(1) % state % ssh % array(iCell)  
               end do

               do iEdge=1,block % mesh % nEdges

                  ! normalBarotropicVelocitySubcycleOld = normalBarotropicVelocityOld 
                    block % state % time_levs(oldBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge) &
                  = block % state % time_levs(1) % state % normalBarotropicVelocity % array(iEdge) 

                  ! normalBarotropicVelocityNew = BtrOld  This is the first for the summation
                    block % state % time_levs(2) % state % normalBarotropicVelocity % array(iEdge) &
                  = block % state % time_levs(1) % state % normalBarotropicVelocity % array(iEdge) 

                  ! barotropicThicknessFlux = 0  
                  block % diagnostics % barotropicThicknessFlux % array(iEdge) = 0.0
               end do

               block => block % next
            end do  ! block

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! BEGIN Barotropic subcycle loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do j=1,config_n_btr_subcycles*config_btr_subcycle_loop_factor

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: VELOCITY PREDICTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if (config_btr_gam1_velWt1>1.0e-12) then  ! only do this part if it is needed in next SSH solve
                  uPerpTime = oldBtrSubcycleTime

                  block => domain % blocklist
                  do while (associated(block))

                     do iEdge=1,block % mesh % nEdges

                        cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
                        cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

                        ! Compute the barotropic Coriolis term, -f*uPerp
                        CoriolisTerm = 0.0
                        do i = 1,block % mesh % nEdgesOnEdge % array(iEdge)
                           eoe = block % mesh % edgesOnEdge % array(i,iEdge)
                           CoriolisTerm = CoriolisTerm &
                             + block % mesh % weightsOnEdge % array(i,iEdge) &
                             * block % state % time_levs(uPerpTime) % state % normalBarotropicVelocitySubcycle % array(eoe) &
                             * block % mesh % fEdge % array(eoe)
                        end do
      
                        ! normalBarotropicVelocityNew = normalBarotropicVelocityOld + dt/J*(-f*normalBarotropicVelocityoldPerp - g*grad(SSH) + G)
                        block % state % time_levs(newBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge) &
                          = (block % state % time_levs(oldBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge) &
                          + dt / config_n_btr_subcycles * (CoriolisTerm - gravity &
                          * (block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell2) &
                           - block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell1) ) &
                          / block % mesh % dcEdge % array(iEdge) &
                          + block % diagnostics % barotropicForcing % array(iEdge))) * block % mesh % edgeMask % array(1, iEdge)
                     end do

                     block => block % next
                  end do  ! block

                !   boundary update on normalBarotropicVelocityNew
                call mpas_timer_start("se halo normalBarotropicVelocity", .false., timer_halo_normalBarotropicVelocity)
                call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(newBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle)
                call mpas_timer_stop("se halo normalBarotropicVelocity", timer_halo_normalBarotropicVelocity)
              endif ! config_btr_gam1_velWt1>1.0e-12

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic subcycle: SSH PREDICTOR STEP 
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              block => domain % blocklist
              do while (associated(block))
      
                block % tend % ssh % array(:) = 0.0
      
                if (config_btr_solve_SSH2) then
                   ! If config_btr_solve_SSH2=.true., then do NOT accumulate barotropicThicknessFlux in this SSH predictor 
                   ! section, because it will be accumulated in the SSH corrector section.
                   barotropicThicknessFlux_coeff = 0.0
                else
                   ! otherwise, DO accumulate barotropicThicknessFlux in this SSH predictor section
                   barotropicThicknessFlux_coeff = 1.0
                endif
      
                ! config_btr_gam1_velWt1 sets the forward weighting of velocity in the SSH computation
                ! config_btr_gam1_velWt1=  1     flux = normalBarotropicVelocityNew*H
                ! config_btr_gam1_velWt1=0.5     flux = 1/2*(normalBarotropicVelocityNew+normalBarotropicVelocityOld)*H
                ! config_btr_gam1_velWt1=  0     flux = normalBarotropicVelocityOld*H

                do iCell = 1, block % mesh % nCells
                  do i = 1, block % mesh % nEdgesOnCell % array(iCell)
                    iEdge = block % mesh % edgesOnCell % array(i, iCell)

                    cell1 = block % mesh % cellsOnEdge % array(1, iEdge)
                    cell2 = block % mesh % cellsOnEdge % array(2, iEdge)

                    sshEdge = 0.5 * (block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell1) &
                              + block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell2) )

                   ! method 0: orig, works only without pbc:      
                   !thicknessSum = sshEdge + block % mesh % refBottomDepthTopOfCell % array (block % mesh % maxLevelEdgeTop % array(iEdge)+1)
 
                   ! method 1, matches method 0 without pbcs, works with pbcs.
                   thicknessSum = sshEdge + min(block % mesh % bottomDepth % array(cell1), &
                                        block % mesh % bottomDepth % array(cell2))

                   ! method 2: may be better than method 1.
                   ! Take average  of full thickness at two neighboring cells.
                   !thicknessSum = sshEdge + 0.5 *(  block % mesh % bottomDepth % array(cell1) &
                   !                       + block % mesh % bottomDepth % array(cell2) )


                    flux = ((1.0-config_btr_gam1_velWt1) * block % state % time_levs(oldBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge) &
                           + config_btr_gam1_velWt1 * block % state % time_levs(newBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge)) &
                           * thicknessSum 

                    block % tend % ssh % array(iCell) = block % tend % ssh % array(iCell) + block % mesh % edgeSignOncell % array(i, iCell) * flux &
                           * block % mesh % dvEdge % array(iEdge)

                  end do
                end do

                do iEdge=1,block % mesh % nEdges
                   cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
                   cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

                   sshEdge = 0.5 * (block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell1) &
                             + block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell2) )

                   ! method 0: orig, works only without pbc:      
                   !thicknessSum = sshEdge + block % mesh % refBottomDepthTopOfCell % array (block % mesh % maxLevelEdgeTop % array(iEdge)+1)
 
                   ! method 1, matches method 0 without pbcs, works with pbcs.
                   thicknessSum = sshEdge + min(block % mesh % bottomDepth % array(cell1), &
                                        block % mesh % bottomDepth % array(cell2))

                   ! method 2: may be better than method 1.
                   ! take average  of full thickness at two neighboring cells
                   !thicknessSum = sshEdge + 0.5 *(  block % mesh % bottomDepth % array(cell1) &
                   !                       + block % mesh % bottomDepth % array(cell2) )

                   flux = ((1.0-config_btr_gam1_velWt1) * block % state % time_levs(oldBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge) &
                          + config_btr_gam1_velWt1 * block % state % time_levs(newBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge)) &
                          * thicknessSum 

                   block % diagnostics % barotropicThicknessFlux % array(iEdge) = block % diagnostics % barotropicThicknessFlux % array(iEdge) &
                     + barotropicThicknessFlux_coeff*flux
                end do
      
                ! SSHnew = SSHold + dt/J*(-div(Flux))
                do iCell=1,block % mesh % nCells 
      
                   block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(iCell) & 
                       = block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(iCell) & 
                       + dt/config_n_btr_subcycles * block % tend % ssh % array(iCell) / block % mesh % areaCell % array (iCell)
      
                end do
      
                block => block % next
              end do  ! block
      
              !   boundary update on SSHnew
              call mpas_timer_start("se halo ssh", .false., timer_halo_ssh)
              call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle)
              call mpas_timer_stop("se halo ssh", timer_halo_ssh)
      
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic subcycle: VELOCITY CORRECTOR STEP
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do BtrCorIter=1,config_n_btr_cor_iter
                uPerpTime = newBtrSubcycleTime
      
                block => domain % blocklist
                do while (associated(block))
                   allocate(utemp(block % mesh % nEdges+1))
                   uTemp(:) = block % state % time_levs(newBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(:)
                   do iEdge=1,block % mesh % nEdges 
                     cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
                     cell2 = block % mesh % cellsOnEdge % array(2,iEdge)
      
                     ! Compute the barotropic Coriolis term, -f*uPerp
                     CoriolisTerm = 0.0
                     do i = 1,block % mesh % nEdgesOnEdge % array(iEdge)
                         eoe = block % mesh % edgesOnEdge % array(i,iEdge)
                       CoriolisTerm = CoriolisTerm + block % mesh % weightsOnEdge % array(i,iEdge) &
                             !* block % state % time_levs(uPerpTime) % state % normalBarotropicVelocitySubcycle % array(eoe) &
                             * uTemp(eoe) &
                             * block % mesh % fEdge  % array(eoe) 
                     end do
      
                     ! In this final solve for velocity, SSH is a linear
                     ! combination of SSHold and SSHnew.
                     sshCell1 = (1-config_btr_gam2_SSHWt1)*block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell1) &
                               +   config_btr_gam2_SSHWt1 *block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(cell1)
                     sshCell2 = (1-config_btr_gam2_SSHWt1)*block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell2) &
                               +   config_btr_gam2_SSHWt1 *block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(cell2)
    
                     ! normalBarotropicVelocityNew = normalBarotropicVelocityOld + dt/J*(-f*normalBarotropicVelocityoldPerp - g*grad(SSH) + G)
                     block % state % time_levs(newBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge) & 
                         = (block % state % time_levs(oldBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge) & 
                         + dt/config_n_btr_subcycles *(CoriolisTerm - gravity *(sshCell2 - sshCell1) /block % mesh % dcEdge % array(iEdge) &
                         + block % diagnostics % barotropicForcing % array(iEdge))) * block % mesh % edgeMask % array(1,iEdge)
                   end do
                   deallocate(uTemp)
      
                   block => block % next
                end do  ! block
      
                !   boundary update on normalBarotropicVelocityNew
                call mpas_timer_start("se halo normalBarotropicVelocity", .false., timer_halo_normalBarotropicVelocity)
                call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(newBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle)
                call mpas_timer_stop("se halo normalBarotropicVelocity", timer_halo_normalBarotropicVelocity)
              end do !do BtrCorIter=1,config_n_btr_cor_iter
      
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic subcycle: SSH CORRECTOR STEP
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (config_btr_solve_SSH2) then
      
                block => domain % blocklist
                do while (associated(block))
                   block % tend % ssh % array(:) = 0.0
      
                  ! config_btr_gam3_velWt2 sets the forward weighting of velocity in the SSH computation
                  ! config_btr_gam3_velWt2=  1     flux = normalBarotropicVelocityNew*H
                  ! config_btr_gam3_velWt2=0.5     flux = 1/2*(normalBarotropicVelocityNew+normalBarotropicVelocityOld)*H
                  ! config_btr_gam3_velWt2=  0     flux = normalBarotropicVelocityOld*H

                  do iCell = 1, block % mesh % nCells
                    do i = 1, block % mesh % nEdgesOnCell % array(iCell)
                      iEdge = block % mesh % edgesOnCell % array(i, iCell)

                      cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
                      cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

                      ! SSH is a linear combination of SSHold and SSHnew.
                      sshCell1 = (1-config_btr_gam2_SSHWt1)*block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell1) &
                                +   config_btr_gam2_SSHWt1 *block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(cell1)
                      sshCell2 = (1-config_btr_gam2_SSHWt1)*block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell2) &
                                +   config_btr_gam2_SSHWt1 *block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(cell2)
 
                      sshEdge = 0.5 * (sshCell1 + sshCell2)

                     ! method 0: orig, works only without pbc:      
                     !thicknessSum = sshEdge + block % mesh % refBottomDepthTopOfCell % array (block % mesh % maxLevelEdgeTop % array(iEdge)+1)
 
                     ! method 1, matches method 0 without pbcs, works with pbcs.
                     thicknessSum = sshEdge + min(block % mesh % bottomDepth % array(cell1), &
                                          block % mesh % bottomDepth % array(cell2))

                     ! method 2: may be better than method 1.
                     ! take average  of full thickness at two neighboring cells
                     !thicknessSum = sshEdge + 0.5 *(  block % mesh % bottomDepth % array(cell1) &
                     !                       + block % mesh % bottomDepth % array(cell2) )
      
       
                      flux = ((1.0-config_btr_gam3_velWt2) * block % state % time_levs(oldBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge) &
                             + config_btr_gam3_velWt2 * block % state % time_levs(newBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge)) &
                             * thicknessSum

                      block % tend % ssh % array(iCell) = block % tend % ssh % array(iCell) + block % mesh % edgeSignOnCell % array(i, iCell) * flux &
                             * block % mesh % dvEdge % array(iEdge)

                    end do
                  end do

                  do iEdge=1,block % mesh % nEdges
                     cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
                     cell2 = block % mesh % cellsOnEdge % array(2,iEdge)
      
                     ! SSH is a linear combination of SSHold and SSHnew.
                     sshCell1 = (1-config_btr_gam2_SSHWt1)*block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell1) &
                               +   config_btr_gam2_SSHWt1 *block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(cell1)
                     sshCell2 = (1-config_btr_gam2_SSHWt1)*block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(cell2) &
                               +   config_btr_gam2_SSHWt1 *block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(cell2)
                     sshEdge = 0.5 * (sshCell1 + sshCell2)

                     ! method 0: orig, works only without pbc:      
                     !thicknessSum = sshEdge + block % mesh % refBottomDepthTopOfCell % array (block % mesh % maxLevelEdgeTop % array(iEdge)+1)
 
                     ! method 1, matches method 0 without pbcs, works with pbcs.
                     thicknessSum = sshEdge + min(block % mesh % bottomDepth % array(cell1), &
                                          block % mesh % bottomDepth % array(cell2))

                     ! method 2, better, I think.
                     ! take average  of full thickness at two neighboring cells
                     !thicknessSum = sshEdge + 0.5 *(  block % mesh % bottomDepth % array(cell1) &
                     !                       + block % mesh % bottomDepth % array(cell2) )
      
                     flux = ((1.0-config_btr_gam3_velWt2) * block % state % time_levs(oldBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge) &
                            + config_btr_gam3_velWt2 * block % state % time_levs(newBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge)) &
                            * thicknessSum
      
                     block % diagnostics % barotropicThicknessFlux % array(iEdge) = block % diagnostics % barotropicThicknessFlux % array(iEdge) + flux
                  end do
      
                  ! SSHnew = SSHold + dt/J*(-div(Flux))
                  do iCell=1,block % mesh % nCells 
                    block % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle % array(iCell) & 
                          = block % state % time_levs(oldBtrSubcycleTime) % state % sshSubcycle % array(iCell) & 
                          + dt/config_n_btr_subcycles * block % tend % ssh % array(iCell) / block % mesh % areaCell % array (iCell)
                  end do
      
                  block => block % next
                end do  ! block
      
                !   boundary update on SSHnew
                call mpas_timer_start("se halo ssh", .false., timer_halo_ssh)
                call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(newBtrSubcycleTime) % state % sshSubcycle)
                call mpas_timer_stop("se halo ssh", timer_halo_ssh)
               endif ! config_btr_solve_SSH2
      
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: Accumulate running sums, advance timestep pointers
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
               block => domain % blocklist
               do while (associated(block))
      
                  ! normalBarotropicVelocityNew = normalBarotropicVelocityNew + normalBarotropicVelocitySubcycleNEW
                  ! This accumulates the sum.
                  ! If the Barotropic Coriolis iteration is limited to one, this could 
                  ! be merged with the above code.
                  do iEdge=1,block % mesh % nEdges 
      
                       block % state % time_levs(2) % state % normalBarotropicVelocity % array(iEdge) &
                     = block % state % time_levs(2) % state % normalBarotropicVelocity % array(iEdge) & 
                     + block % state % time_levs(newBtrSubcycleTime) % state % normalBarotropicVelocitySubcycle % array(iEdge)  
      
                  end do  ! iEdge
                  block => block % next
               end do  ! block
      
               ! advance time pointers
               oldBtrSubcycleTime = mod(oldBtrSubcycleTime,2)+1
               newBtrSubcycleTime = mod(newBtrSubcycleTime,2)+1
      
            end do ! j=1,config_n_btr_subcycles
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! END Barotropic subcycle loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Normalize Barotropic subcycle sums: ssh, normalBarotropicVelocity, and F
            block => domain % blocklist
            do while (associated(block))
      
               do iEdge=1,block % mesh % nEdges
                  block % diagnostics % barotropicThicknessFlux % array(iEdge) = block % diagnostics % barotropicThicknessFlux % array(iEdge) &
                      / (config_n_btr_subcycles*config_btr_subcycle_loop_factor)
      
                  block % state % time_levs(2) % state % normalBarotropicVelocity % array(iEdge) = block % state % time_levs(2) % state % normalBarotropicVelocity % array(iEdge) & 
                     / (config_n_btr_subcycles*config_btr_subcycle_loop_factor + 1)
               end do
      
               block => block % next
            end do  ! block
      
      
            ! boundary update on F
            call mpas_timer_start("se halo F", .false., timer_halo_f)
            call mpas_dmpar_exch_halo_field(domain % blocklist % diagnostics % barotropicThicknessFlux)
            call mpas_timer_stop("se halo F", timer_halo_f)


            ! Check that you can compute SSH using the total sum or the individual increments
            ! over the barotropic subcycles.
            ! efficiency: This next block of code is really a check for debugging, and can 
            ! be removed later.
            block => domain % blocklist
            do while (associated(block))

               allocate(uTemp(block % mesh % nVertLevels))

               ! Correction velocity    uCorr = (Flux - Sum(h u*))/H
               ! or, for the full latex version:
               !{\bf u}^{corr} = \left( {\overline {\bf F}} 
               !  - \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}  {\bf u}_k^{avg} \right)
               ! \left/ \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}   \right. 

               if (config_vel_correction) then
                  ucorr_coef = 1
               else
                  ucorr_coef = 0
               endif

               do iEdge=1,block % mesh % nEdges

                  ! velocity for uCorrection is normalBarotropicVelocity + normalBaroclinicVelocity + uBolus
                  uTemp(:) &
                     = block % state % time_levs(2) % state % normalBarotropicVelocity % array(  iEdge) &
                     + block % state % time_levs(2) % state % normalBaroclinicVelocity % array(:,iEdge) &
                     + block % diagnostics % uBolusGM % array(:,iEdge)

                  ! thicknessSum is initialized outside the loop because on land boundaries 
                  ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a 
                  ! nonzero value to avoid a NaN.
                  normalThicknessFluxSum = block % diagnostics % layerThicknessEdge % array(1,iEdge) * uTemp(1)
                  thicknessSum  = block % diagnostics % layerThicknessEdge % array(1,iEdge)

                  do k=2,block % mesh % maxLevelEdgeTop % array(iEdge)
                     normalThicknessFluxSum = normalThicknessFluxSum + block % diagnostics % layerThicknessEdge % array(k,iEdge) * uTemp(k)
                     thicknessSum  =  thicknessSum + block % diagnostics % layerThicknessEdge % array(k,iEdge)
                  enddo

                  uCorr =   ucorr_coef*(( block % diagnostics % barotropicThicknessFlux % array(iEdge) - normalThicknessFluxSum)/thicknessSum)

                  do k=1,block % mesh % nVertLevels

                     ! uTranport = normalBarotropicVelocity + normalBaroclinicVelocity + uBolus + uCorrection
                     ! This is u used in advective terms for layerThickness and tracers 
                     ! in tendency calls in stage 3.
                     block % diagnostics % uTransport % array(k,iEdge) &
                           = block % mesh % edgeMask % array(k,iEdge) &
                           *(  block % state % time_levs(2) % state % normalBarotropicVelocity % array(  iEdge) &
                           + block % state % time_levs(2) % state % normalBaroclinicVelocity % array(k,iEdge) &
                           + block % diagnostics % uBolusGM   % array(k,iEdge) &
                           + uCorr )

                  enddo

               end do ! iEdge

               deallocate(uTemp)

               block => block % next
            end do  ! block

         endif ! split_explicit  

         call mpas_timer_stop("se btr vel", timer_btr_vel)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 3: Tracer, density, pressure, vertical velocity prediction
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! Thickness tendency computations and thickness halo updates are completed before tracer 
         ! tendency computations to allow monotonic advection.
         block => domain % blocklist
         do while (associated(block))

            ! compute vertTransportVelocityTop.  Use uTransport for advection of layerThickness and tracers.
            ! Use time level 1 values of layerThickness and layerThicknessEdge because 
            ! layerThickness has not yet been computed for time level 2.
            call ocn_vert_transport_velocity_top(block % mesh, block % verticalMesh, &
               block % state % time_levs(1) % state % layerThickness % array, &
               block % diagnostics % layerThicknessEdge % array, &
               block % diagnostics % uTransport % array, &
               block % state % time_levs(1) % state % ssh % array, &
               block % state % time_levs(2) % state % highFreqThickness % array, dt, &
               block % diagnostics % vertTransportVelocityTop % array, err)

            call ocn_tend_thick(block % tend, block % state % time_levs(2) % state, block % diagnostics, block % mesh)

            block => block % next
         end do

         ! update halo for thickness tendencies
         call mpas_timer_start("se halo thickness", .false., timer_halo_thickness)
         call mpas_dmpar_exch_halo_field(domain % blocklist % tend % layerThickness)
         call mpas_timer_stop("se halo thickness", timer_halo_thickness)

         block => domain % blocklist
         do while (associated(block))
            call ocn_tend_tracer(block % tend, block % state % time_levs(2) % state, block % diagnostics, block % mesh, dt)

            block => block % next
         end do

         ! update halo for tracer tendencies
         call mpas_timer_start("se halo tracers", .false., timer_halo_tracers)
         call mpas_dmpar_exch_halo_field(domain % blocklist % tend % tracers)
         call mpas_timer_stop("se halo tracers", timer_halo_tracers)

         block => domain % blocklist
         do while (associated(block))

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !  If iterating, reset variables for next iteration
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (split_explicit_step < config_n_ts_iter) then

               ! Get indices for dynamic tracers (Includes T&S).
               startIndex = block % state % time_levs(1) % state % dynamics_start
               endIndex = block % state % time_levs(1) % state % dynamics_end

               ! Only need T & S for earlier iterations,
               ! then all the tracers needed the last time through.
               do iCell=1,block % mesh % nCells
                  ! sshNew is a pointer, defined above.
                  do k=1,block % mesh % maxLevelCell % array(iCell)

                     ! this is h_{n+1}
                     temp_h &
                        = block % state % time_levs(1) % state % layerThickness % array(k,iCell) &
                        + dt* block % tend % layerThickness % array(k,iCell) 

                     ! this is h_{n+1/2}
                       block % state % time_levs(2) % state % layerThickness % array(k,iCell) &
                     = 0.5*(  &
                       block % state % time_levs(1) % state % layerThickness % array(k,iCell) &
                       + temp_h)

                     do i=startIndex, endIndex
                        ! This is Phi at n+1
                        temp = (  &
                           block % state % time_levs(1) % state % tracers % array(i,k,iCell) &
                         * block % state % time_levs(1) % state % layerThickness % array(k,iCell) &
                         + dt * block % tend % tracers % array(i,k,iCell)) &
                              / temp_h
  
                        ! This is Phi at n+1/2
                          block % state % time_levs(2) % state % tracers % array(i,k,iCell) &
                        = 0.5*( &
                          block % state % time_levs(1) % state % tracers % array(i,k,iCell) &
                          + temp )
                     end do
                  end do
               end do ! iCell

               if (config_use_freq_filtered_thickness) then
                  do iCell=1,block % mesh % nCells
                     do k=1,block % mesh % maxLevelCell % array(iCell)

                        ! h^{hf}_{n+1} was computed in Stage 1

                        ! this is h^{hf}_{n+1/2}
                        block % state % time_levs(2) % state % highFreqThickness % array(k,iCell) &
                           = 0.5*(block % state % time_levs(1) % state % highFreqThickness % array(k,iCell) &
                                + block % state % time_levs(2) % state % highFreqThickness % array(k,iCell))

                        ! this is D^{lf}_{n+1}
                        temp = block % state % time_levs(1) % state % lowFreqDivergence % array(k,iCell) &
                         + dt* block % tend % lowFreqDivergence % array(k,iCell) 

                        ! this is D^{lf}_{n+1/2}
                        block % state % time_levs(2) % state % lowFreqDivergence % array(k,iCell) &
                           = 0.5*(block % state % time_levs(1) % state % lowFreqDivergence % array(k,iCell) + temp)

                     end do
                  end do
               end if

               do iEdge=1,block % mesh % nEdges

                  do k=1,block % mesh % nVertLevels

                     ! u = normalBarotropicVelocity + normalBaroclinicVelocity 
                     ! here normalBaroclinicVelocity is at time n+1/2
                     ! This is u used in next iteration or step
                       block % state % time_levs(2) % state % normalVelocity    % array(k,iEdge) &
                     = block % mesh % edgeMask % array(k,iEdge) &
                     *(  block % state % time_levs(2) % state % normalBarotropicVelocity % array(  iEdge) &
                       + block % state % time_levs(2) % state % normalBaroclinicVelocity % array(k,iEdge) )

                  enddo

               end do ! iEdge

               ! Efficiency note: We really only need this to compute layerThicknessEdge, density, pressure, and SSH 
               ! in this diagnostics solve.
               call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % forcing, block % mesh, block % diagnostics, block % scratch)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !  If large iteration complete, compute all variables at time n+1
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif (split_explicit_step == config_n_ts_iter) then

               do iCell=1,block % mesh % nCells
                  do k=1,block % mesh % maxLevelCell % array(iCell)

                     ! this is h_{n+1}
                        block % state % time_levs(2) % state % layerThickness % array(k,iCell) &
                      = block % state % time_levs(1) % state % layerThickness % array(k,iCell) &
                      + dt* block % tend % layerThickness % array(k,iCell) 

                     ! This is Phi at n+1
                     do i=1,block % state % time_levs(1) % state % num_tracers
                           block % state % time_levs(2) % state % tracers % array(i,k,iCell)  &
                        = (block % state % time_levs(1) % state % tracers % array(i,k,iCell) &
                         * block % state % time_levs(1) % state % layerThickness % array(k,iCell) &
                         + dt * block % tend % tracers % array(i,k,iCell)) &
                         / block % state % time_levs(2) % state % layerThickness % array(k,iCell)

                     enddo
                  end do
               end do

               if (config_use_freq_filtered_thickness) then
                  do iCell=1,block % mesh % nCells
                     do k=1,block % mesh % maxLevelCell % array(iCell)

                        ! h^{hf}_{n+1} was computed in Stage 1

                        ! this is D^{lf}_{n+1}
                           block % state % time_levs(2) % state % lowFreqDivergence % array(k,iCell) &
                         = block % state % time_levs(1) % state % lowFreqDivergence % array(k,iCell) &
                         + dt* block % tend % lowFreqDivergence % array(k,iCell) 

                     end do
                  end do
               end if

               ! Recompute final u to go on to next step.
               ! u_{n+1} = normalBarotropicVelocity_{n+1} + normalBaroclinicVelocity_{n+1} 
               ! Right now normalBaroclinicVelocityNew is at time n+1/2, so back compute to get normalBaroclinicVelocity at time n+1
               !   using normalBaroclinicVelocity_{n+1/2} = 1/2*(normalBaroclinicVelocity_n + u_Bcl_{n+1})
               ! so the following lines are
               ! u_{n+1} = normalBarotropicVelocity_{n+1} + 2*normalBaroclinicVelocity_{n+1/2} - normalBaroclinicVelocity_n
               ! note that normalBaroclinicVelocity is recomputed at the beginning of the next timestep due to Imp Vert mixing,
               ! so normalBaroclinicVelocity does not have to be recomputed here.
      
               do iEdge=1,block % mesh % nEdges
                  do k=1,block % mesh % maxLevelEdgeTop % array(iEdge)
                       block % state % time_levs(2) % state % normalVelocity    % array(k,iEdge) &
                     = block % state % time_levs(2) % state % normalBarotropicVelocity % array(  iEdge) &
                    +2*block % state % time_levs(2) % state % normalBaroclinicVelocity % array(k,iEdge) &
                     - block % state % time_levs(1) % state % normalBaroclinicVelocity % array(k,iEdge)
                  end do
               end do ! iEdges

            endif ! split_explicit_step

            block => block % next
         end do



      end do  ! split_explicit_step = 1, config_n_ts_iter
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END large iteration loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_start("se implicit vert mix")
      block => domain % blocklist
      do while(associated(block))

        ! Call ocean diagnostic solve in preparation for vertical mixing.  Note 
        ! it is called again after vertical mixing, because u and tracers change.
        ! For Richardson vertical mixing, only density, layerThicknessEdge, and kineticEnergyCell need to 
        ! be computed.  For kpp, more variables may be needed.  Either way, this
        ! could be made more efficient by only computing what is needed for the
        ! implicit vmix routine that follows.
        call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % forcing, block % mesh, block % diagnostics, block % scratch)

        call ocn_vmix_implicit(dt, block % mesh, block % diagnostics, block % state % time_levs(2) % state, err)

        block => block % next
      end do

      ! Update halo on u and tracers, which were just updated for implicit vertical mixing.  If not done, 
      ! this leads to lack of volume conservation.  It is required because halo updates in stage 3 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to 
      ! communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start("se implicit vert mix halos")
      call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % normalVelocity)
      call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % tracers)
      call mpas_timer_stop("se implicit vert mix halos")

      call mpas_timer_stop("se implicit vert mix")

      block => domain % blocklist
      do while (associated(block))

         if (config_prescribe_velocity) then
            block % state % time_levs(2) % state % normalVelocity % array(:,:) = block % state % time_levs(1) % state % normalVelocity % array(:,:)
         end if

         if (config_prescribe_thickness) then
            block % state % time_levs(2) % state % layerThickness % array(:,:) = block % state % time_levs(1) % state % layerThickness % array(:,:)
         end if

         call ocn_diagnostic_solve(dt, block % state % time_levs(2) % state, block % forcing, block % mesh, block % diagnostics, block % scratch)

         ! Compute velocity transport, used in advection terms of layerThickness and tracer tendency
         block % diagnostics % uTransport % array(:,:) &
               = block % state % time_levs(2) % state % normalVelocity % array(:,:) &
               + block % diagnostics % uBolusGM % array(:,:)

         call mpas_reconstruct(block % mesh, block % state % time_levs(2) % state % normalVelocity % array,          &
                          block % diagnostics % normalVelocityX % array,            &
                          block % diagnostics % normalVelocityY % array,            &
                          block % diagnostics % normalVelocityZ % array,            &
                          block % diagnostics % normalVelocityZonal % array,        &
                          block % diagnostics % normalVelocityMeridional % array    &
                         )

         call ocn_time_average_accumulate(block % average, block % state % time_levs(2) % state, block % diagnostics)

         block => block % next
      end do

      call mpas_timer_stop("se timestep", timer_main)

   end subroutine ocn_time_integrator_split!}}}

end module ocn_time_integration_split

! vim: foldmethod=marker
