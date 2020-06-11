! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection
!
!> \brief MPAS ocean tracer advection driver
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This module contains initialization and driver routines for computing tracer
!>  advection tendencies. It primarily calls submodule routines based on user
!>  choices for advection options.
!
!-------------------------------------------------------------------------------

module ocn_tracer_advection

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_sort
   use mpas_hash
   use mpas_timer

   use ocn_tracer_advection_std
   use ocn_tracer_advection_mono

   use ocn_constants

   implicit none
   private
   save

   ! public module method interfaces
   public :: ocn_tracer_advection_init,         &
             ocn_tracer_advection_tend

   ! privat module variables
   logical :: tracerAdvOn !< flag to turn on tracer advection
   logical :: monotonicOn !< flag to choose a monotone advection scheme
   logical :: budgetDiagsOn !< flag to compute active tracer budgets

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_tracer_advection_tend
!
!> \brief MPAS ocean tracer advection tendency
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This routine is the driver routine for computing tracer advection
!>  tendencies. It simply calls submodule tendency routines based on choice of
!>  algorithm.
!
!-------------------------------------------------------------------------------

   subroutine ocn_tracer_advection_tend(tracers, normalThicknessFlux, w,       &
                                        layerThickness, dt, meshPool,          &
                                        scratchPool, diagnosticsPool, tend,    &
                                        tracerGroupName)!{{{

      !*** Input/Output parameters

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend            !< [in,out] Tracer tendency to which advection added

      !*** Input parameters

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers             !< [in] Current tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalThicknessFlux !< [in] Thickness weighted horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         w                   !< [in] Vertical velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness      !< [in] Thickness field
      real (kind=RKIND), intent(in) :: &
         dt                  !< [in] Time step
      type (mpas_pool_type), intent(in) :: &
         meshPool            !< [in] Mesh information
      type (mpas_pool_type), intent(in) :: &
         scratchPool         !< [in] Scratch fields
      type (mpas_pool_type), intent(in) :: &
         diagnosticsPool     !< [in] Diagnostic fields
      character (len=*), intent(in) :: &
         tracerGroupName     !< [in] Tracer name to check if active tracers

      !*** Local variables

      logical :: &
         computeBudgets      ! flag to compute active tracer budget

      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         advCoefs, advCoefs3rd ! advection coefficients

      integer, dimension(:), pointer, contiguous :: &
         maxLevelCell,    &! index of max level at each cell
         maxLevelEdgeTop, &! max level at edge with both cells active
         nAdvCellsForEdge  ! number of advective cells for each edge
      integer, dimension(:,:), pointer, contiguous :: &!
         highOrderAdvectionMask, &! mask for higher order contributions
         edgeSignOnCell,  &! sign at cell edge for fluxes
         advCellsForEdge   ! index of advective cells for each edge

      ! end of preamble
      !----------------
      ! begin code

      ! immediate return if tracer advection not selected
      if(.not. tracerAdvOn) return

      call mpas_timer_start("tracer adv")

      ! extract pool variables
      call mpas_pool_get_array(meshPool, 'advCoefs', advCoefs)
      call mpas_pool_get_array(meshPool, 'advCoefs3rd', advCoefs3rd)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask', &
                                          highOrderAdvectionMask)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(meshPool, 'advCellsForEdge', advCellsForEdge)

      ! determine whether active tracer budgets should be computed
      computeBudgets = (budgetDiagsOn .and. tracerGroupName == 'activeTracers')

      ! call specific advection routine based on choice of monotonicity
      if (monotonicOn) then
         call ocn_tracer_advection_mono_tend(tend, tracers, layerThickness,    &
                                             normalThicknessFlux, w, dt,       &
                                             advCoefs, advCoefs3rd,            &
                                             nAdvCellsForEdge, advCellsForEdge,&
                                             maxLevelCell, maxLevelEdgeTop,    &
                                             highOrderAdvectionMask,           &
                                             edgeSignOnCell, meshPool,         &
                                             diagnosticsPool, computeBudgets)

      else
         call ocn_tracer_advection_std_tend(tracers, advCoefs, advCoefs3rd, &
            nAdvCellsForEdge, advCellsForEdge, normalThicknessFlux, w, layerThickness, &
            layerThickness, dt, meshPool, tend, maxLevelCell, maxLevelEdgeTop, &
            highOrderAdvectionMask, edgeSignOnCell)
      endif

      call mpas_timer_stop("tracer adv")

   end subroutine ocn_tracer_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_tracer_advection_init
!
!> \brief MPAS ocean tracer advection tendency
!> \author Mark Petersen, David Lee, Doug Jacobsen, Phil Jones
!> \date   October 2017, updated May 2019
!> \details
!>  This routine is the driver routine for initializing various tracer
!>  advection choices and variables.
!
!-------------------------------------------------------------------------------

   subroutine ocn_tracer_advection_init(err)!{{{

      !*** output parameters

      integer, intent(out) :: err !< [out] Error flag

      !*** local variables

      integer :: errTmp ! temporary error flag

      ! pointer variables needed to extract options from config
      logical, pointer ::        &
         disableTracerAdvection, &! flag to disable advection
         checkMonotonicity,      &! flag to check monotonicity
         useMonotone,            &! flag to choose monotone scheme
         computeBudgetDiags       ! flag to compute active tracer budgets

      integer, pointer :: &
         horzAdvOrder,    &! choice of horizontal advection order
         vertAdvOrder,    &! choice of vertical advection order
         numHalos          ! width of halo region

      real (kind=RKIND), pointer :: &
         coef3rdOrder  ! coefficient for blending high-order horz adv terms

      ! end preamble
      !-------------
      ! begin code

      err = 0 ! initialize error code to success

      ! extract options from config
      call mpas_pool_get_config(ocnConfigs,'config_num_halos', numHalos)
      call mpas_pool_get_config(ocnConfigs,'config_disable_tr_adv', &
                                            disableTracerAdvection)
      call mpas_pool_get_config(ocnConfigs,'config_check_tracer_monotonicity', &
                                            checkMonotonicity)
      call mpas_pool_get_config(ocnConfigs,'config_horiz_tracer_adv_order', &
                                            horzAdvOrder)
      call mpas_pool_get_config(ocnConfigs,'config_vert_tracer_adv_order', &
                                            vertAdvOrder)
      call mpas_pool_get_config(ocnConfigs,'config_coef_3rd_order', &
                                            coef3rdOrder)
      call mpas_pool_get_config(ocnConfigs,'config_monotonic', &
                                            useMonotone)
      call mpas_pool_get_config(ocnConfigs, &
                                'config_compute_active_tracer_budgets', &
                                 computeBudgetDiags)

      ! set some basic flags for options
      tracerAdvOn = .not. disableTracerAdvection
      monotonicOn = useMonotone
      budgetDiagsOn = computeBudgetDiags

      ! set all other options from submodule initialization routines
      call ocn_tracer_advection_std_init(horzAdvOrder, vertAdvOrder,      &
                                         coef3rdOrder, checkMonotonicity, &
                                         errTmp)
      call ocn_tracer_advection_mono_init(numHalos, horzAdvOrder,    &
                                          vertAdvOrder, coef3rdOrder, &
                                          checkMonotonicity, err)

      ! if an error is returned from init routines, write an error
      ! message and return a non-zero error code
      if (err /= 0 .or. errTmp /= 0) then
         err = 1
         call mpas_log_write(                                 &
            'Error encountered during tracer advection init', &
            MPAS_LOG_ERR, masterOnly=.true.)
      endif

   end subroutine ocn_tracer_advection_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_tracer_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
