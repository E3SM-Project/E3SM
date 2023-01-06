!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_gotm
!
!> \brief MPAS ocean vertical mixing interface to GOTM
!> \author Qing Li
!> \date   11 May 2020
!> \details
!>  This module contains the routines for calls into GOTM
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_gotm

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_log
   use mpas_timekeeping

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

   use mpas_io_units

   use turbulence, only: init_turbulence, do_turbulence, clean_turbulence
   use turbulence, only: gotm_vertVisc=>num, gotm_vertDiffT=>nuh, &
                         gotm_vertDiffS=>nus, gotm_TKE=>tke,      &
                         gotm_eps=>eps, gotm_L=>L,                &
                         gotm_kb=>kb, gotm_epsb=>epsb
   use mtridiagonal, only: init_tridiagonal

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

   public :: ocn_vmix_coefs_gotm_build, &
             ocn_vmix_gotm_init,        &
             ocn_vmix_gotm_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: gotmOn
   real (kind=RKIND) :: gotm_dt
   real (kind=RKIND) :: gotmBottomDragCoef, gotmBottomRoughnessLength, &
                        gotmSurfaceRoughnessLength


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs_gotm_build
!
!> \brief   Computes mixing coefficients using GOTM
!> \author  Qing Li
!> \date    11 May 2020
!> \details
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers by calling GOTM routines.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs_gotm_build(statePool, forcingPool, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in), optional ::     &
         timeLevelIn !< Input: time level for state pool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool         !< Input/Output: state information

      type (mpas_pool_type), intent(inout) :: &
         forcingPool   !< Input/Output: forcing information

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

      real (kind=RKIND), dimension(:,:), pointer :: &
        normalVelocity, layerThickness, &
        gotmVertViscTopOfCellCur, gotmVertDiffTopOfCellCur, &
        gotmTKETopOfCellCur, gotmDissTopOfCellCur, gotmLengthTopOfCellCur, &
        gotmKbTopOfCellCur, gotmEpsbTopOfCellCur, &
        gotmVertViscTopOfCellNew, gotmVertDiffTopOfCellNew, &
        gotmTKETopOfCellNew, gotmDissTopOfCellNew, gotmLengthTopOfCellNew, &
        gotmKbTopOfCellNew, gotmEpsbTopOfCellNew

      integer :: i, iCell, iEdge, k, ki, timeLevel, nCells

      real (kind=RKIND) :: invAreaCell, factor, delU2, bottomFrictionVelocity, waterDepth
      real (kind=RKIND), dimension(:), allocatable :: &
        shearSquaredTop, NSquaredTop, hLayer

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing mixing-related fields
      ! note that the user can choose multiple options and the
      !   mixing fields have to be added/merged together
      !
      !-----------------------------------------------------------------

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      !
      ! only build up viscosity/diffusivity if GOTM is turned on
      !
      if ( .not. gotmOn ) return

      !
      ! set parameters
      !
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

      !
      ! set pointers for GOTM viscosity/diffusivity, TKE, rate of dissipation and turbulent length scale
      ! (current and new state)
      !
      call mpas_pool_get_array(statePool, 'gotmVertViscTopOfCell', gotmVertViscTopOfCellCur, 1)
      call mpas_pool_get_array(statePool, 'gotmVertViscTopOfCell', gotmVertViscTopOfCellNew, 2)
      call mpas_pool_get_array(statePool, 'gotmVertDiffTopOfCell', gotmVertDiffTopOfCellCur, 1)
      call mpas_pool_get_array(statePool, 'gotmVertDiffTopOfCell', gotmVertDiffTopOfCellNew, 2)
      call mpas_pool_get_array(statePool, 'gotmTKETopOfCell', gotmTKETopOfCellCur, 1)
      call mpas_pool_get_array(statePool, 'gotmTKETopOfCell', gotmTKETopOfCellNew, 2)
      call mpas_pool_get_array(statePool, 'gotmDissTopOfCell', gotmDissTopOfCellCur, 1)
      call mpas_pool_get_array(statePool, 'gotmDissTopOfCell', gotmDissTopOfCellNew, 2)
      call mpas_pool_get_array(statePool, 'gotmLengthTopOfCell', gotmLengthTopOfCellCur, 1)
      call mpas_pool_get_array(statePool, 'gotmLengthTopOfCell', gotmLengthTopOfCellNew, 2)
      call mpas_pool_get_array(statePool, 'gotmKbTopOfCell', gotmKbTopOfCellCur, 1)
      call mpas_pool_get_array(statePool, 'gotmKbTopOfCell', gotmKbTopOfCellNew, 2)
      call mpas_pool_get_array(statePool, 'gotmEpsbTopOfCell', gotmEpsbTopOfCellCur, 1)
      call mpas_pool_get_array(statePool, 'gotmEpsbTopOfCell', gotmEpsbTopOfCellNew, 2)

      nCells = nCellsAll

      !$omp parallel
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         vertViscTopOfCell(:, iCell) = 0.0_RKIND
         vertDiffTopOfCell(:, iCell) = 0.0_RKIND
      end do
      !$omp end do
      !$omp end parallel

      nCells = nCellsHalo( 2 )

      ! allocate arrays for squared shear, squared buoyancy frequency and layer thickness
      ! to be used in GOTM, the indices are in the GOTM order (bottom to top)
      allocate(shearSquaredTop(0:nVertLevels))
      allocate(NSquaredTop(0:nVertLevels))
      allocate(hLayer(0:nVertLevels))

      call mpas_timer_start('gotm cell loop', .false.)

      ! MDT, Mar 2021: Removed OpenMP directives because the subroutine
      !                do_turbulence is not thread safe
      do iCell = 1, nCells

         ! initialize arrays
         shearSquaredTop = 0.0_RKIND
         NSquaredTop = 0.0_RKIND
         hLayer = 0.0_RKIND

         ! compute shear frequency squared
         invAreaCell = 1.0_RKIND / areaCell(iCell)
         do k = 1, maxLevelCell(iCell) - minLevelCell(iCell)
            ki = maxLevelCell(iCell) - k + 1
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               factor = 0.5_RKIND * dcEdge(iEdge) * dvEdge(iEdge) * invAreaCell
               delU2 = (normalVelocity(ki-1,iEdge) - normalVelocity(ki,iEdge))**2
               shearSquaredTop(k) = shearSquaredTop(k) + factor * delU2
            end do
            shearSquaredTop(k) = shearSquaredTop(k) / (zMid(ki-1,iCell) - zMid(ki,iCell))**2
         end do
         shearSquaredTop(0) = shearSquaredTop(1)
         shearSquaredTop(maxLevelCell(iCell)-minLevelCell(iCell)+1) = shearSquaredTop(maxLevelCell(iCell)-minLevelCell(iCell))

         ! fill arrays of the squared buoyancy frequency and layer thickness
         ! and compute the water depth
         waterDepth = 0.0_RKIND
         do k = 1, maxLevelCell(iCell) - minLevelCell(iCell) + 1
            ki = maxLevelCell(iCell) - k + 1
            NSquaredTop(k) = BruntVaisalaFreqTop(ki,iCell)
            hLayer(k) = layerThickness(ki,iCell)
            waterDepth = waterDepth + layerThickness(ki,iCell)
         end do

         ! compute the bottom friction velocity
         bottomFrictionVelocity = sqrt(2.0_RKIND * gotmBottomDragCoef *  &
                                       kineticEnergyCell(maxLevelCell(iCell),iCell))

         ! fill GOTM arrays with current values at iCell
         do k = 0, maxLevelCell(iCell) - minLevelCell(iCell) + 1
            ki = maxLevelCell(iCell) - k + 1
            gotm_vertVisc(k) = gotmVertViscTopOfCellCur(ki, iCell)
            gotm_vertDiffT(k) = gotmVertDiffTopOfCellCur(ki, iCell)
            gotm_vertDiffS(k) = gotmVertDiffTopOfCellCur(ki, iCell)
            gotm_TKE(k) = gotmTKETopOfCellCur(ki, iCell)
            gotm_eps(k) = gotmDissTopOfCellCur(ki, iCell)
            gotm_L(k) = gotmLengthTopOfCellCur(ki, iCell)
            gotm_kb(k) = gotmKbTopOfCellCur(ki, iCell)
            gotm_epsb(k) = gotmEpsbTopOfCellCur(ki, iCell)
         end do

         !
         ! interface with do_turbulence():
         !
         !   subroutine do_turbulence(nlev,dt,depth,u_taus,u_taub,z0s,z0b,h,NN,SS,xP)
         !
         !   nlev          : number of vertical layers
         !   dt            : time step (s)
         !   depth         : water depth (m)
         !   u_taus        : surface friction velocity (m/s)
         !   u_taub        : bottom friction velocity (m/s)
         !   z0s           : surface roughness length (m)
         !   z0b           : bottom roughness length (m)
         !   h             : layer thickness (m)
         !   NN            : bouyancy frequency squared (1/s^2)
         !   SS            : shear-frequency squared (1/s^2)
         !   xP (optional) : extra TKE production due to, e.g., seagrass (m^2/s^3)
         !
         ! MPAS-O (1:nlev+1) <-> GOTM (nlev:0)
         !
         call do_turbulence(maxLevelCell(iCell) - minLevelCell(iCell) + 1,                &
                            gotm_dt,                                                      &
                            waterDepth,                                                   &
                            surfaceFrictionVelocity(iCell),                               &
                            bottomFrictionVelocity,                                       &
                            gotmSurfaceRoughnessLength,                                   &
                            gotmBottomRoughnessLength,                                    &
                            hLayer(0:maxLevelCell(iCell)-minLevelCell(iCell)+1),          &
                            NSquaredTop(0:maxLevelCell(iCell)-minLevelCell(iCell)+1),     &
                            shearSquaredTop(0:maxLevelCell(iCell)-minLevelCell(iCell)+1))

         ! update GOTM arrays for next time step
         ! MPAS-O (1:nlev+1) <-> GOTM (nlev:0)
         do k = minLevelCell(iCell), maxLevelCell(iCell)+1
            ki = maxLevelCell(iCell) - k + 1
            gotmVertViscTopOfCellNew(k, iCell) = gotm_vertVisc(ki)
            gotmVertDiffTopOfCellNew(k, iCell) = gotm_vertDiffT(ki)
            gotmTKETopOfCellNew(k, iCell) = gotm_TKE(ki)
            gotmDissTopOfCellNew(k, iCell) = gotm_eps(ki)
            gotmLengthTopOfCellNew(k, iCell) = gotm_L(ki)
            gotmKbTopOfCellNew(k, iCell) = gotm_kb(ki)
            gotmEpsbTopOfCellNew(k, iCell) = gotm_epsb(ki)
            ! update viscosity/diffusivity
            vertViscTopOfCell(k, iCell) = gotm_vertVisc(ki)
            vertDiffTopOfCell(k, iCell) = gotm_vertDiffT(ki)
         end do

      end do

      call mpas_timer_stop('gotm cell loop')

      deallocate(shearSquaredTop)
      deallocate(NSquaredTop)
      deallocate(hLayer)

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_gotm_build!}}}

!***********************************************************************
!
!  routine ocn_vmix_gotm_init
!
!> \brief   Initializes ocean vertical mixing quantities using GOTM
!> \author  Qing Li
!> \date    11 May 2020
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical mixing in the ocean. Parameters are set by calling into GOTM
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_gotm_init(domain,err)!{{{

   !--------------------------------------------------------------------

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool

      integer, pointer :: nVertLevels, nVertLevelsP1
      integer, dimension(:), pointer :: nCellsArray

      real (kind=RKIND), dimension(:,:), pointer :: &
        gotmVertViscTopOfCellCur, gotmVertDiffTopOfCellCur, &
        gotmTKETopOfCellCur, gotmDissTopOfCellCur, gotmLengthTopOfCellCur, &
        gotmKbTopOfCellCur, gotmEpsbTopOfCellCur

      logical :: testFile
      integer :: namlstUnit !< unit for GOTM namelist file
      integer :: k, ki, iCell, nCells
      type (MPAS_timeInterval_type) :: timeStep

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      gotmOn = config_use_gotm

      !
      ! only initialize if GOTM is turned on
      !
      if (.not. gotmOn) return

      !
      ! When GOTM is turned on, all other vertical mixing schemes should be off
      ! Test to make sure this is the case.
      !
      ! test here, err=1 if a problem

      ! use MPAS-O dt here. Sub-stepping is possible, though it may not help too much
      ! without updating the shear and stratification in the substeps
      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err)
      call mpas_get_timeInterval(timeStep, dt=gotm_dt)
      !
      ! pull nVertLevels out of the mesh structure
      !
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

      ! initialize GOTM parameters
      gotmSurfaceRoughnessLength = config_gotm_constant_surface_roughness_length
      gotmBottomRoughnessLength = config_gotm_constant_bottom_roughness_length
      gotmBottomDragCoef = config_gotm_constant_bottom_drag_coeff

      !
      ! initialize GOTM turbulence module
      ! check if GOTM namelist file exist or not
      testFile = .false.
      inquire(file=trim(config_gotm_namelist_file), exist=testFile)
      if (testFile) then
         ! read in GOTM namelist, required by init_turbulence
         call mpas_new_unit(namlstUnit)
         call init_turbulence(namlstUnit, trim(config_gotm_namelist_file), nVertLevels)
         call mpas_release_unit(namlstUnit)
      else
         call mpas_log_write(' fatal error: gotm namelist file does not exist: config_gotm_namelist_file = '//config_gotm_namelist_file, MPAS_LOG_ERR)
         err = 1
         return
      end if

      ! initialize GOTM tridiagnoal solver
      call init_tridiagonal(nVertLevels)

      !
      ! no need to initialize GOTM arrays for restart runs
      ! as these arrays are read from the restart file
      !
      if (config_do_restart) return

      ! initialize GOTM arrays
      block => domain%blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
         call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

         ! set pointers for GOTM viscosity/diffusivity, TKE, rate of dissipation and turbulent length scale
         ! (current and new state)
         call mpas_pool_get_array(statePool, 'gotmVertViscTopOfCell', gotmVertViscTopOfCellCur, 1)
         call mpas_pool_get_array(statePool, 'gotmVertDiffTopOfCell', gotmVertDiffTopOfCellCur, 1)
         call mpas_pool_get_array(statePool, 'gotmTKETopOfCell', gotmTKETopOfCellCur, 1)
         call mpas_pool_get_array(statePool, 'gotmDissTopOfCell', gotmDissTopOfCellCur, 1)
         call mpas_pool_get_array(statePool, 'gotmLengthTopOfCell', gotmLengthTopOfCellCur, 1)
         call mpas_pool_get_array(statePool, 'gotmKbTopOfCell', gotmKbTopOfCellCur, 1)
         call mpas_pool_get_array(statePool, 'gotmEpsbTopOfCell', gotmEpsbTopOfCellCur, 1)

         nCells = nCellsArray( size(nCellsArray) )
         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp private(k, ki)
         do iCell = 1, nCells
            ! MPAS-O (1:nlev+1) <-> GOTM (nlev:0)
            do k = 1, nVertLevelsP1
               ki = nVertLevelsP1 - k
               gotmVertViscTopOfCellCur(k, iCell) = gotm_vertVisc(ki)
               gotmVertDiffTopOfCellCur(k, iCell) = gotm_vertDiffT(ki)
               gotmTKETopOfCellCur(k, iCell) = gotm_TKE(ki)
               gotmDissTopOfCellCur(k, iCell) = gotm_eps(ki)
               gotmLengthTopOfCellCur(k, iCell) = gotm_L(ki)
               gotmKbTopOfCellCur(k, iCell) = gotm_kb(ki)
               gotmEpsbTopOfCellCur(k, iCell) = gotm_epsb(ki)
            end do
         end do
         !$omp end do
         !$omp end parallel

         block => block%next
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_gotm_init!}}}

!***********************************************************************
!
!  routine ocn_vmix_gotm_finalize
!
!> \brief   Finalizes ocean vertical mixing quantities in GOTM
!> \author  Qing Li
!> \date    11 May 2020
!> \details
!>  This routine deallocates all memory allocated in GOTM
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_gotm_finalize(domain,err)!{{{

   !--------------------------------------------------------------------

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      !
      ! only finalize if GOTM is turned on
      !
      if (.not. gotmOn) return

      ! deallocate all memory allocated in GOTM
      call clean_turbulence()

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_gotm_finalize!}}}


!***********************************************************************

end module ocn_vmix_gotm

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
