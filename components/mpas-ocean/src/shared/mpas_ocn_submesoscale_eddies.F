! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

module ocn_submesoscale_eddies

   use mpas_pool_routines
   use mpas_derived_types
   use mpas_constants
   use mpas_threading
   use mpas_timer

   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_mesh

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

   public :: ocn_submesoscale_compute_velocity, &
             ocn_submesoscale_init,             &
             ocn_submesoscale_finalize,         &
             ocn_MLE_add_to_transport_vel

   real(kind=RKIND) :: tau, LfMin, Ce, dsMax
   real(kind=RKIND),dimension(:),allocatable :: time_scale

  !***********************************************************************

  contains

  !***********************************************************************
  !
  !  routine ocn_submesoscale_compute_velocity
  !
  !> \brief   Computes velocity from submesoscale eddies
  !> \details
  !>  This routine implements the Fox-Kemper et al 2011 submesoscale 
  !>  parameterization (https://doi.org/10.1016/j.ocemod.2010.09.002).
  !>  The transport velocity from the submesoscales is given by a 
  !>  stream function which is given as
  !>
  !>  Psi = C_e*Delta_S/L_f*H^2 int(bGrad_H)/sqrt(f^2+tau^-2)*mu(z)
  !>
  !>  here C_e is a specified efficiency, Delta_s is taken as 
  !>  min(dcEdge,dsMax), H is the mixed layer depth, L_f is the subgrid
  !>  frontal width that is scaled to the coarse grid, f is the 
  !>  coriolis parameter, and tau is a time scale specified to prevent 
  !>  the singularity near the equator.  the horizontal buoyancy gradient
  !>  is integrated over the mixed layer depth and mu(z) is a non-dimensional
  !>  shape function that is zero below the mixed layer depth
  !>
  !***********************************************************************

  subroutine ocn_submesoscale_compute_velocity()

     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------

     integer :: k, nEdges, nCells, iCell, iEdge, cell1, cell2, indMLDedge

     real(kind=RKIND),dimension(:),allocatable :: &
        streamFunction, &
        zEdge

     real(kind=RKIND) :: mu, zMLD, bvfML, hML, gradBuoyML, ds, Lf, bvfAv, &
        hAv, bldEdge, ustarEdge

     nEdges = nEdgesHalo(2)

     allocate(streamFunction(nVertLevels+1))
     allocate(zEdge(nVertLevels+1))

     !compute the mixed layer average buoyancy gradient and bvf
     !This is thickness weighted to allow for non uniform grid spacing, 
     !and centered on layer interfaces. The mixed layer goes from the 
     !ocean top layer to mid-depth of layer(indMLDedge)
     !NOTE, gradBuoyEddy and BFV are defined at the top cell interface 
     !and not valid for minLevelCell, so properties for the top layer 
     !rely on the values at (minLevelCell+1).
     !$omp parallel
     !$omp do schedule(runtime) &
     !$omp private(cell1,cell2,bvfML,zEdge,streamFunction,gradBuoyML,hML,bvfAv, &
     !$omp         zMLD,Lf,ds,mu,k)
     do iEdge=1,nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        zEdge(:) = 0.0_RKIND
        streamFunction(:) = 0.0_RKIND
        bvfML = 0.0_RKIND
        gradBuoyML = 0.0_RKIND
        !the BFV is defined at cell interfaces in the vertical so there is an offset in index
        !this is the top half layer, assuming homogeneity of top layer
        bvfAv = sqrt(0.5_RKIND*(max(BruntVaisalaFreqTop(minLevelCell(cell1)+1,cell1),1.0E-20_RKIND) + &
                       max(BruntVaisalaFreqTop(minLevelCell(cell2)+1,cell2),1.0E-20_RKIND)))
        hML = 0.5_RKIND*layerThickEdgeMean(minLevelEdgeBot(iEdge),iEdge)
        bvfML = bvfML + hML*bvfAv
        if (minLevelEdgeTop(iEdge) .ge. 1) then
           gradBuoyML = hML*gradBuoyEddy(minLevelEdgeTop(iEdge)+1,iEdge)
           bvfML = hML*bvfAv
        else
           cycle
        end if

        indMLDedge = min(indMLD(cell1),indMLD(cell2))
        do k = minLevelEdgeBot(iEdge)+1,indMLDedge
           hAv = 0.5_RKIND*(layerThickEdgeMean(k,iEdge) + layerThickEdgeMean(k-1,iEdge))
           bvfAv = sqrt(0.5_RKIND*(max(BruntVaisalaFreqTop(k,cell1),1.0E-20_RKIND) + &
                              max(BruntVaisalaFreqTop(k,cell2),1.0E-20_RKIND)))
           bvfML = bvfML + hAv*bvfAv
           hML = hML + hAv
           gradBuoyML = gradBuoyML + hAv*gradBuoyEddy(k,iEdge)
        end do
        bvfML = bvfML / (1.0E-20_RKIND + hML)
        gradBuoyML = gradBuoyML / (1.0E-20_RKIND + hML)

        !compute depths and shape function

        do k = minLevelEdgeTop(iEdge)+1,maxLevelEdgeTop(iEdge)+1
           zEdge(k) = zEdge(k-1) - layerThickEdgeMean(k-1,iEdge)
         end do

        zMLD = 0.5_RKIND*(dThreshMLD(cell1)+dThreshMLD(cell2))
        Lf = max(Lfmin, abs(gradBuoyML)*zMLD / time_scale(iEdge)**2.0, &
                 bvfML*zMLD / time_scale(iEdge))

        ds = min(dcEdge(iEdge),dsMax)

        do k = minLevelEdgeTop(iEdge)+1,maxLevelEdgeTop(iEdge)
           mu = max(0.0_RKIND,(1.0_RKIND - (2.0_RKIND*zEdge(k) / zMLD + 1.0_RKIND)**2.0)* &
                    (1.0_RKIND + 5.0_RKIND/21.0_RKIND*(2.0_RKIND*zEdge(k)/zMLD + 1.0_RKIND)**2.0))
           streamFunction(k) = Ce*ds/Lf*zMLD**2.0*gradBuoyML/time_scale(iEdge)*mu
        end do

        ! integrate in vertical to get the velocity
        do k = minLevelEdgeTop(iEdge),maxLevelEdgeTop(iEdge)
           normalMLEvelocity(k,iEdge) = -(streamFunction(k) - streamFunction(k+1)) / layerThickEdgeMean(k,iEdge)
        end do

    end do!iEdge loop
    !$omp end do
    !$omp end parallel

    deallocate(streamFunction)
    deallocate(zEdge)

  end subroutine ocn_submesoscale_compute_velocity

  !***********************************************************************
  !
  !  routine ocn_add_MLE_to_transport_vel
  !
  !> \brief   Submesoscale parameterization add to transport Vel
  !> \details
  !>  adds the MLE induced velocity to the transport velocity 
  !>  
  !***********************************************************************

  subroutine ocn_MLE_add_to_transport_vel(normalTransportVelocity, nEdges)

     real(kind=RKIND), dimension(:,:), intent(inout) :: &
         normalTransportVelocity

     integer, intent(in) :: nEdges

     integer :: iEdge, k

     if (.not. config_submesoscale_enable) return

     !$omp parallel
     !$omp do schedule(runtime) &
     !$omp private(k)
     do iEdge = 1, nEdges
        do k = 1, nVertLevels
           normalTransportVelocity(k,iEdge) = normalTransportVelocity(k,iEdge) + &
               normalMLEvelocity(k,iEdge)
        end do
     end do
     !$omp end do
     !$omp end parallel

  end subroutine ocn_MLE_add_to_transport_vel

  !***********************************************************************
  !
  !  routine ocn_submesoscale_init
  !
  !> \brief   Submesoscale parameterization init
  !> \details
  !>  Initializes parameters related to the submesoscale eddy parameterization 
  !>  
  !***********************************************************************

  subroutine ocn_submesoscale_init(err)

     integer, intent(out) :: err !< Output: error flag

     integer :: iEdge

     err = 0
     Ce = config_submesoscale_ce
     Lfmin = config_submesoscale_lfmin
     dsMax = config_submesoscale_ds_max
     tau = config_submesoscale_tau

     allocate(time_scale(nEdgesAll))

     !$omp parallel
     !$omp do schedule(runtime)
     do iEdge=1,nEdgesAll
        time_scale(iEdge) = sqrt(fEdge(iEdge)**2 + tau**(-2.0))
     end do
     !$omp end do
     !$omp end parallel

   end subroutine ocn_submesoscale_init

   subroutine ocn_submesoscale_finalize()

     deallocate(time_scale)

   end subroutine ocn_submesoscale_finalize

end module ocn_submesoscale_eddies
