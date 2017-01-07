! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_test
!
!> \brief Driver for testing MPAS ocean core
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   October 2013
!> \details
!>  This module contains routines to test various components of
!>  the MPAS ocean core.
!
!-----------------------------------------------------------------------

module ocn_test

   use mpas_constants
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_field_routines
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_tensor_operations

   use ocn_constants

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

   public :: ocn_test_suite, &
             ocn_init_gm_test_functions

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: hmixOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_test_suite
!
!> \brief   Call all internal start-up tests
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    October 2013
!> \details 
!>  Call all routines to test various MPAS-Ocean components.
!
!-----------------------------------------------------------------------

   subroutine ocn_test_suite(domain, err)!{{{

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

      integer :: err1

      err=0

      call ocn_prep_test_tensor(domain,err1)
      err = ior(err1,err)

   end subroutine ocn_test_suite!}}}

!***********************************************************************
!
!  routine ocn_prep_test_tensor
!
!> \brief   set up scratch variables to test strain rate and tensor divergence operators
!> \author  Mark Petersen
!> \date    May 2013
!> \details 
!>  This routine sets up scratch variables to test strain rate and tensor divergence operators.
!
!-----------------------------------------------------------------------

   subroutine ocn_prep_test_tensor(domain,err)!{{{

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

      type (mpas_pool_type), pointer :: scratchPool, meshPool
      logical, pointer :: config_test_tensors
      character (len=StrKIND), pointer :: config_tensor_test_function

      type (field2DInteger), pointer :: edgeSignOnCellField
      type (field2DReal), pointer :: edgeTangentVectorsField
      type (field2DReal), pointer :: normalVelocityTestField
      type (field2DReal), pointer :: tangentialVelocityTestField
      type (field3DReal), pointer :: strainRateR3CellField
      type (field3DReal), pointer :: strainRateR3CellSolutionField
      type (field3DReal), pointer :: strainRateR3EdgeField
      type (field3DReal), pointer :: strainRateLonLatRCellField
      type (field3DReal), pointer :: strainRateLonLatRCellSolutionField
      type (field3DReal), pointer :: strainRateLonLatREdgeField
      type (field3DReal), pointer :: divTensorR3CellField
      type (field3DReal), pointer :: divTensorR3CellSolutionField
      type (field3DReal), pointer :: divTensorLonLatRCellField
      type (field3DReal), pointer :: divTensorLonLatRCellSolutionField
      type (field3DReal), pointer :: outerProductEdgeField

      call mpas_pool_get_config(domain % configs, 'config_test_tensors', config_test_tensors)

      if (.not.config_test_tensors) return

      call mpas_pool_get_config(domain % configs, 'config_tensor_test_function', config_tensor_test_function)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)

      call mpas_pool_get_field(meshPool, 'edgeSignOnCell', edgeSignOnCellField)
      call mpas_pool_get_field(meshPool, 'edgeTangentVectors', edgeTangentVectorsField)

      call mpas_pool_get_field(scratchPool, 'normalVelocityTest', normalVelocityTestField)
      call mpas_pool_get_field(scratchPool, 'tangentialVelocityTest', tangentialVelocityTestField)
      call mpas_pool_get_field(scratchPool, 'strainRateR3Cell', strainRateR3CellField)
      call mpas_pool_get_field(scratchPool, 'strainRateR3CellSolution', strainRateR3CellSolutionField)
      call mpas_pool_get_field(scratchPool, 'strainRateR3Edge', strainRateR3EdgeField)
      call mpas_pool_get_field(scratchPool, 'strainRateLonLatRCell', strainRateLonLatRCellField)
      call mpas_pool_get_field(scratchPool, 'strainRateLonLatRCellSolution', strainRateLonLatRCellSolutionField)
      call mpas_pool_get_field(scratchPool, 'strainRateLonLatREdge', strainRateLonLatREdgeField)
      call mpas_pool_get_field(scratchPool, 'divTensorR3Cell', divTensorR3CellField)
      call mpas_pool_get_field(scratchPool, 'divTensorR3CellSolution', divTensorR3CellSolutionField)
      call mpas_pool_get_field(scratchPool, 'divTensorLonLatRCell', divTensorLonLatRCellField)
      call mpas_pool_get_field(scratchPool, 'divTensorLonLatRCellSolution', divTensorLonLatRCellSolutionField)
      call mpas_pool_get_field(scratchPool, 'outerProductEdge', outerProductEdgeField)

      call mpas_allocate_scratch_field(normalVelocityTestField, .false.)
      call mpas_allocate_scratch_field(tangentialVelocityTestField, .false.)
      call mpas_allocate_scratch_field(strainRateR3CellField, .false.)
      call mpas_allocate_scratch_field(strainRateR3CellSolutionField, .false.)
      call mpas_allocate_scratch_field(strainRateR3EdgeField, .false.)
      call mpas_allocate_scratch_field(strainRateLonLatRCellField, .false.)
      call mpas_allocate_scratch_field(strainRateLonLatRCellSolutionField, .false.)
      call mpas_allocate_scratch_field(strainRateLonLatREdgeField, .false.)
      call mpas_allocate_scratch_field(divTensorR3CellField, .false.)
      call mpas_allocate_scratch_field(divTensorR3CellSolutionField, .false.)
      call mpas_allocate_scratch_field(divTensorLonLatRCellField, .false.)
      call mpas_allocate_scratch_field(divTensorLonLatRCellSolutionField, .false.)
      call mpas_allocate_scratch_field(outerProductEdgeField, .false.)

      call mpas_test_tensor(domain, config_tensor_test_function, &
         edgeSignOnCellField, &
         edgeTangentVectorsField, &
         normalVelocityTestField, &
         tangentialVelocityTestField, &
         strainRateR3CellField, &
         strainRateR3CellSolutionField, &
         strainRateR3EdgeField, &
         strainRateLonLatRCellField, &
         strainRateLonLatRCellSolutionField, &
         strainRateLonLatREdgeField, &
         divTensorR3CellField, &
         divTensorR3CellSolutionField, &
         divTensorLonLatRCellField, &
         divTensorLonLatRCellSolutionField, &
         outerProductEdgeField )

      call mpas_deallocate_scratch_field(normalVelocityTestField, .false.)
      call mpas_deallocate_scratch_field(tangentialVelocityTestField, .false.)
      call mpas_deallocate_scratch_field(strainRateR3CellField, .false.)
      call mpas_deallocate_scratch_field(strainRateR3CellSolutionField, .false.)
      call mpas_deallocate_scratch_field(strainRateR3EdgeField, .false.)
      call mpas_deallocate_scratch_field(strainRateLonLatRCellField, .false.)
      call mpas_deallocate_scratch_field(strainRateLonLatRCellSolutionField, .false.)
      call mpas_deallocate_scratch_field(strainRateLonLatREdgeField, .false.)
      call mpas_deallocate_scratch_field(divTensorR3CellField, .false.)
      call mpas_deallocate_scratch_field(divTensorR3CellSolutionField, .false.)
      call mpas_deallocate_scratch_field(divTensorLonLatRCellField, .false.)
      call mpas_deallocate_scratch_field(divTensorLonLatRCellSolutionField, .false.)
      call mpas_deallocate_scratch_field(outerProductEdgeField, .false.)

      err = 0

   end subroutine ocn_prep_test_tensor!}}}

!***********************************************************************
!
!  routine ocn_init_gm_test_functions
!
!> \brief   Initialize Gent-McWilliams test functions
!> \author  Mark Petersen
!> \date    May 2014
!> \details 
!>  For the initial temperature distribution
!>  T = T_1 + T_2*y/y_{max} + T_3*z/z_{max}
!>  and linear EOS with T coefficient alpha, this subroutine computes
!>  the instantaneous analytic solution for:
!>    - the Bolus stream function
!>    - horizontal Bolus velocity
!
!-----------------------------------------------------------------------

   subroutine ocn_init_gm_test_functions(diagnosticsPool, meshPool, scratchPool)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool            !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: diagnosticsPool  !< Input: diagnostic fields derived from State
      type (mpas_pool_type), intent(in) :: scratchPool !< Input: scratch variables

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: k, iCell
      integer, pointer  :: nCells
      integer, dimension(:), pointer   :: maxLevelCell

      real(kind=RKIND) :: zTop, config_gm_analytic_temperature2, config_gm_analytic_temperature3, config_gm_analytic_ymax, &
         config_gm_analytic_bottom_depth, L, R, c1, c2, zMax, zBot
      real (kind=RKIND), pointer :: config_gravWaveSpeed_trunc, config_density0, config_standardGM_tracer_kappa, config_eos_linear_alpha

      real(kind=RKIND), dimension(:), pointer   :: bottomDepth, refBottomDepthTopOfCell, yCell, yEdge
      real(kind=RKIND), dimension(:,:), pointer :: yRelativeSlopeSolution, yGMStreamFuncSolution, yGMBolusVelocitySolution, zMid

      type (field2DReal), pointer :: yRelativeSlopeSolutionField, yGMStreamFuncSolutionField, yGMBolusVelocitySolutionField

      call mpas_pool_get_config(ocnConfigs, 'config_density0',config_density0)
      call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_alpha', config_eos_linear_alpha)
      call mpas_pool_get_config(ocnConfigs, 'config_gravWaveSpeed_trunc',config_gravWaveSpeed_trunc)
      call mpas_pool_get_config(ocnConfigs, 'config_standardGM_tracer_kappa',config_standardGM_tracer_kappa)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell',refBottomDepthTopOfCell)
      call mpas_pool_get_array(meshPool, 'yCell',yCell)
      call mpas_pool_get_array(meshPool, 'yEdge',yEdge)
      call mpas_pool_get_array(meshPool, 'bottomDepth',bottomDepth)

      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

      ! move solution variables to diagnostics pool to include them in output.
      call mpas_pool_get_field(scratchPool, 'yRelativeSlopeSolution',yRelativeSlopeSolutionField)
      call mpas_pool_get_field(scratchPool, 'yGMStreamFuncSolution',yGMStreamFuncSolutionField)
      call mpas_pool_get_field(scratchPool, 'yGMBolusVelocitySolution',yGMBolusVelocitySolutionField)

      yRelativeSlopeSolution   => yRelativeSlopeSolutionField % array
      yGMStreamFuncSolution    => yGMStreamFuncSolutionField % array
      yGMBolusVelocitySolution => yGMBolusVelocitySolutionField % array 

      ! These are flags that must match your initial conditions settings.  See gm_analytic initial condition in mode_init.
      config_gm_analytic_temperature2 = 10;
      config_gm_analytic_temperature3 = -10;
      config_gm_analytic_ymax = 500000;
      config_gm_analytic_bottom_depth = 1000;

      ! zMax is associated with linear temperature profile in z
      zMax = -config_gm_analytic_bottom_depth;  
      ! zBot is location we apply boundary conditions on the ODE for stream function.  
      zBot = zMax;  

      L = config_gravWaveSpeed_trunc * sqrt(config_density0*zMax/gravity/config_eos_linear_alpha/config_gm_analytic_temperature3);
      R = - config_standardGM_tracer_kappa * config_gm_analytic_temperature2 * zMax / config_gm_analytic_temperature3 / config_gm_analytic_ymax;
      c1 = R*(1-exp(-zBot/L))/(exp(zBot/L) - exp(-zBot/L));
      c2 = R-c1;

      do iCell = 1, nCells

         do k = 1, maxLevelCell(iCell)
            ! placed at mid-depth of cell center:
            yGMBolusVelocitySolution(k,iCell) = 1/L*(c1*exp(zMid(k,iCell)/L) - c2*exp(-zMid(k,iCell)/L) );
         end do

         do k = 1, maxLevelCell(iCell)
            ! placed at top interface, cell center.
            zTop = - refBottomDepthTopOfCell(k)
            yGMStreamFuncSolution(k,iCell) = c1*exp(zTop/L) + c2*exp(-zTop/L) - R;

         end do

         k = maxLevelCell(iCell)+1
            ! placed at top interface, cell center.
            zTop = zBot
            yGMStreamFuncSolution(k,iCell) = c1*exp(zTop/L) + c2*exp(-zTop/L) - R;

      end do

   end subroutine ocn_init_gm_test_functions!}}}

end module ocn_test

! vim: foldmethod=marker
