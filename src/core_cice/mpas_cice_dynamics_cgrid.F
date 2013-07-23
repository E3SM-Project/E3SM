module cice_dynamics_cgrid

  use mpas_grid_types

  implicit none

  private
  public :: 



contains

  !-------------------------------------------------------------
  ! Initialization
  !-------------------------------------------------------------



  !-------------------------------------------------------------
  ! time step
  !-------------------------------------------------------------
  
  subroutine strain_tensor_cgrid(mesh, &
                                 strain11,       strain22,       &
                                 strain12,                       &
                                 velocity,                       &
                                 solveStress)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), dimension(:), intent(in) :: &
         velocity

    integer, dimension(:), intent(in) :: &
         solveStress

    integer :: &
         iCell, &
         iEdgeOnCell, &
         iEdge

    real(kind=RKIND) :: &
         uVelocityEdge, &
         vVelocityEdge
    
    ! loop over cells
    do iCell = 1, mesh % nCells

       strain11(iCell) = 0.0_RKIND
       strain22(iCell) = 0.0_RKIND
       strain12(iCell) = 0.0_RKIND

       if (solveStress(iCell) == 1) then

          do iEdgeOnCell = 1, mesh % nEdgesOnCell % array(iCell)

             iEdge = mesh % edgesOnCell % array(iEdgeOnCell,iCell)

             uVelocityEdge = velocity(iEdge) * normalVectorPolygon(1,iEdgeOnCell,iCell)
             vVelocityEdge = velocity(iEdge) * normalVectorPolygon(2,iEdgeOnCell,iCell)

             strain11(iCell) = strain11(iCell) + uVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) * mesh % dvEdge % array(iEdge)
             strain22(iCell) = strain22(iCell) + vVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) * mesh % dvEdge % array(iEdge)
             strain12(iCell) = strain12(iCell) + 0.5_RKIND * ( &
                  uVelocityEdge * normalVectorPolygon(2,iEdgeOnCell,iCell) + &
                  vVelocityEdge * normalVectorPolygon(1,iEdgeOnCell,iCell) ) * mesh % dvEdge % array(iEdge)

          enddo ! iEdgeOnCell

          strain11(iCell) = strain11(iCell) / mesh % areaCell % array(iCell)
          strain22(iCell) = strain22(iCell) / mesh % areaCell % array(iCell)
          strain12(iCell) = strain12(iCell) / mesh % areaCell % array(iCell)

       endif ! solveStress

    enddo ! iCell

  end subroutine strain_tensor_cgrid

  !-------------------------------------------------------------

  subroutine stress_divergence(mesh, &
                               stressDivergence,                        &
                               stress11,             stress22,          &
                               stress12,                                &
                               solveVelocity)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergence

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &
         stress22, &
         stress12

    integer :: &
         iEdge

    do iEdge = 1, mesh % nEdges

       stressDivergence(iEdge) = 0.0_RKIND

       if (solveVelocity(iEdge) == 1) then




       endif ! solveVelocity

    enddo ! iEdge

  end subroutine stress_divergence

  !-------------------------------------------------------------

end module cice_dynamics_cgrid
