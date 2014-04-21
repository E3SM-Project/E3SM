module cice_dyn_hex_variational

  use mpas_grid_types

  implicit none

  private
  public :: &
       stress_tensor_hex_variational

contains

  !-------------------------------------------------------------

  subroutine stress_tensor_hex_variational(mesh, &
                                           stress11, stress22, stress12,     &
                                           strain11, strain22, strain12,     &
                                           icePressure, replacementPressure, &
                                           solveStress,                      &
                                           dtElastic,   evpDamping)

    use cice_dynamics_shared, only: evp_constitutive_relation, &
                                    linear_constitutive_relation

    use cice_diagnostics, only: iCellTest

    type(MPAS_pool_type), pointer, intent(in) :: mesh

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), dimension(:), intent(out) :: &
         replacementPressure

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), dimension(:), intent(in) :: &
         icePressure

    integer, dimension(:), intent(in) :: &
         solveStress

    real(kind=RKIND), intent(in) :: &
         dtElastic

    logical, intent(in) :: &
         evpDamping

    real(kind=RKIND) :: &
         replacementPressureCell

    integer :: &
         iCell, &
         iVertexOnCell

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    do iCell = 1, nCells

       if (solveStress(iCell) == 1) then

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             call evp_constitutive_relation(stress11(iVertexOnCell,iCell), stress22(iVertexOnCell,iCell), stress12(iVertexOnCell,iCell), &
                  strain11(iVertexOnCell,iCell), strain22(iVertexOnCell,iCell), strain12(iVertexOnCell,iCell), &
                  icePressure(iCell), replacementPressureCell, &
                  areaCell(iCell),              &
                  dtElastic, evpDamping)

             if (iVertexOnCell == 1) then

                replacementPressure(iCell) = replacementPressureCell

             endif

             !call linear_constitutive_relation(stress11(iVertexOnCell,iCell), stress22(iVertexOnCell,iCell), stress12(iVertexOnCell,iCell), &
             !                                  strain11(iVertexOnCell,iCell), strain22(iVertexOnCell,iCell), strain12(iVertexOnCell,iCell))

          enddo ! iVertexOnCell

       else

          stress11(:,iCell) = 0.0_RKIND
          stress22(:,iCell) = 0.0_RKIND
          stress12(:,iCell) = 0.0_RKIND

       endif ! solveStress

    enddo ! iCell

  end subroutine stress_tensor_hex_variational

  !-------------------------------------------------------------

end module cice_dyn_hex_variational
