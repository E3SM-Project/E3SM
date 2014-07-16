!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver_variational
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_velocity_solver_variational

  use mpas_grid_types

  implicit none

  private
  save

  public :: &
       cice_stress_tensor_variational

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_stress_tensor_variational
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_stress_tensor_variational(&
       mesh, &
       stress11, &
       stress22, &
       stress12, &
       strain11, &
       strain22, &
       strain12, &
       icePressure, &
       replacementPressure, &
       solveStress, &
       dtElastic, &
       evpDamping)!{{{
    
    use cice_velocity_solver_shared, only: &
         cice_evp_constitutive_relation, &
         cice_linear_constitutive_relation

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input: 

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         stress11, & !< Input/Output: 
         stress22, & !< Input/Output: 
         stress12    !< Input/Output: 

    real(kind=RKIND), dimension(:), intent(out) :: &
         replacementPressure !< Output: 

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         strain11, & !< Input: 
         strain22, & !< Input: 
         strain12    !< Input: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         icePressure !< Input: 

    integer, dimension(:), intent(in) :: &
         solveStress !< Input: 

    real(kind=RKIND), intent(in) :: &
         dtElastic !< Input: 

    logical, intent(in) :: &
         evpDamping !< Input: 

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
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    do iCell = 1, nCells

       if (solveStress(iCell) == 1) then

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             call cice_evp_constitutive_relation(&
                  stress11(iVertexOnCell,iCell), &
                  stress22(iVertexOnCell,iCell), &
                  stress12(iVertexOnCell,iCell), &
                  strain11(iVertexOnCell,iCell), &
                  strain22(iVertexOnCell,iCell), &
                  strain12(iVertexOnCell,iCell), &
                  icePressure(iCell), &
                  replacementPressureCell, &
                  areaCell(iCell), &
                  dtElastic, evpDamping)

             if (iVertexOnCell == 1) then

                replacementPressure(iCell) = replacementPressureCell

             endif

             !call cice_linear_constitutive_relation(&
             !     stress11(iVertexOnCell,iCell), &
             !     stress22(iVertexOnCell,iCell), &
             !     stress12(iVertexOnCell,iCell), &
             !     strain11(iVertexOnCell,iCell), &
             !     strain22(iVertexOnCell,iCell), &
             !     strain12(iVertexOnCell,iCell))

          enddo ! iVertexOnCell

       else

          stress11(:,iCell) = 0.0_RKIND
          stress22(:,iCell) = 0.0_RKIND
          stress12(:,iCell) = 0.0_RKIND

       endif ! solveStress

    enddo ! iCell

  end subroutine cice_stress_tensor_variational!}}}

!-----------------------------------------------------------------------

end module cice_velocity_solver_variational
