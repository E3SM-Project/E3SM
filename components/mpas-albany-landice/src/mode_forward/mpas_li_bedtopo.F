! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_bedtopo
!
!> \MPAS land-ice bedtopo driver
!> \author Matt Hoffman and Holly Han (modified)
!> \date   20 June 2019, January 2022 (modified)
!> \details
!>  This module contains the routines for
!>  bed topography for solid earth changes
!
!-----------------------------------------------------------------------

module li_bedtopo
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log
   use li_mask
   use li_setup
   use netcdf
!#ifdef USE_SEALEVELMODEL
   use mpi 
!#endif   
   implicit none
   private

   !--------------------------------------------------------------------
   ! Public parameters
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! Public member functions
   !--------------------------------------------------------------------
   public :: li_bedtopo_init, &
             li_bedtopo_finalize, &
             li_bedtopo_block_init, &
             li_bedtopo_solve

   !--------------------------------------------------------------------
   ! Private module variables
   !--------------------------------------------------------------------
   ! sea-level model timestep
   integer, save :: slmTimeStep

   ! Interpolation weights variables
   integer, dimension(:), allocatable :: toRowValues, toColValues
   integer, dimension(:), allocatable :: fromRowValues, fromColValues
   real, dimension(:), allocatable :: toSValues, fromSValues
   integer:: nMpas, nGrid
   
   ! MPI variables
   integer :: nCellsGlobal
   integer, dimension(:), allocatable :: nCellsDisplacement
   integer, dimension(:), allocatable :: indexToCellIDGathered
   integer, dimension(:), allocatable :: nCellsPerProc
   integer, pointer :: nCellsAll
   integer, pointer ::  nCellsOwned
   integer :: iCell, ilm, curProc
   
   real (kind=RKIND), dimension(:), allocatable :: globalArrayIce, gatheredArrayIce
   real (kind=RKIND), dimension(:), allocatable :: globalArrayBedTopo, gatheredArrayBedTopo
   real (kind=RKIND), dimension(:), allocatable :: globalArraySLC, gatheredArraySLC

!***********************************************************************

contains

!***********************************************************************
!
!  routine li_bedtopo_init
!
!> \brief   Initializes bedtopo solver
!> \author  Matt Hoffman and Holly Han (modified)
!> \date    20 June 2019 (original), December 2021 (modified)
!> \details
!>  This routine initializes the bedtopo solver.
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_init(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      ! No init is needed.
      err = 0

#ifdef USE_SEALEVELMODEL

      ! initialize the 1D sea-level model
      call slmodel_init(domain, err)

#endif

   !--------------------------------------------------------------------

   end subroutine li_bedtopo_init



!***********************************************************************
!
!  routine li_bedtopo_block_init
!
!> \brief   Initializes blocks for bedtopo solver
!> \author  Matt Hoffman
!> \date    20 June 2019
!> \details
!>  This routine initializes each block of the bedtopo solver.
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_block_init(block, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: &
         block          !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! No block init needed.
      err = 0

   !--------------------------------------------------------------------
   end subroutine li_bedtopo_block_init



!***********************************************************************
!
!  subroutine li_bedtopo_solve
!
!> \brief   Updates bed topography
!> \author  Matt Hoffman and Holly Han (modified)
!> \date    20 June 2019 (original), December 2021 (modified)
!> \details
!>  This routine updates the bed topography.  Currently the only option
!>  is a data field passed in as input.
!
!-----------------------------------------------------------------------
subroutine li_bedtopo_solve(domain, err)

      use mpas_timekeeping 
      use li_mask
      use li_advection

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block
	   character (len=StrKIND), pointer :: config_slm_coupling_interval
      character (len=StrKIND), pointer :: config_uplift_method
      type (mpas_pool_type), pointer :: meshPool          !< mesh information
      type (mpas_pool_type), pointer :: geometryPool      !< geometry information
      type (mpas_pool_type), pointer :: velocityPool      !< velocity information

      real (kind=RKIND), dimension(:), pointer :: bedTopography, upliftRate
      real (kind=RKIND), pointer :: deltat
      integer :: err_tmp

      err = 0
      err_tmp = 0

      ! Set needed variables and pointers
      call mpas_pool_get_config(liConfigs, 'config_uplift_method', config_uplift_method)
      if (trim(config_uplift_method)=='none') then
         ! do nothing
      elseif (trim(config_uplift_method)=='data') then

         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)

            call mpas_pool_get_array(meshPool, 'deltat', deltat)
            call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
            call mpas_pool_get_array(geometryPool, 'upliftRate', upliftRate)

            bedTopography(:) = bedTopography(:) + upliftRate(:) * deltat

            call li_update_geometry(geometryPool)
            call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)

            block => block % next
         end do

      elseif (trim(config_uplift_method)=='sealevelmodel') then

			call mpas_log_write("Sea-level model is used to calculate changes in bedTopography")

      else
         call mpas_log_write("Unknown option selected for 'config_uplift_method'", MPAS_LOG_ERR)
      endif

#ifdef USE_SEALEVELMODEL

      if (trim(config_uplift_method)=='sealevelmodel') then
			call mpas_pool_get_config(liConfigs, 'config_slm_coupling_interval', config_slm_coupling_interval)

	      if (mpas_is_alarm_ringing(domain % clock, 'slmCouplingInterval', ierr=err_tmp)) then
	         err = ior(err, err_tmp)

				slmTimeStep = slmTimeStep + 1

	         call mpas_log_write("<<< SL alarm is ringing. slmTimeStep: $i", intArgs=(/slmTimeStep/))
            call slmodel_solve(slmTimeStep, domain)

            call mpas_reset_clock_alarm(domain % clock, 'slmCouplingInterval', ierr=err_tmp)
            err = ior(err, err_tmp)

	      else
	        ! do nothing for now, but could calculate uplift rate here later instead.
	      endif
		else
         call mpas_log_write("'sealevelmodel' should be selected for option 'config_uplift_method'", MPAS_LOG_ERR)
		endif
#endif

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_bedtopo_solve.")
      endif

   !--------------------------------------------------------------------
   end subroutine li_bedtopo_solve



!***********************************************************************
!
!  routine li_bedtopo_finalize
!
!> \brief   finalizes bedtopo solver
!> \author  Matt Hoffman
!> \date    20 June 2019
!> \details
!>  This routine finalizes the bedtopo solver.
!
!-----------------------------------------------------------------------

   subroutine li_bedtopo_finalize(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      err = 0

   !--------------------------------------------------------------------
   end subroutine li_bedtopo_finalize


   ! private subroutines

!***********************************************************************
!
!  routine slmodel_init
!
!> \brief   Initializes the sea-level model
!> \author  Holly Kyeore Han
!> \date    January 2022
!> \details
!>  This wrapper routine initializes the sea-level solver
!> (Han et al., 2022, GMD, https://doi.org/10.5281/zenodo.5775235)
!
!-----------------------------------------------------------------------

   subroutine slmodel_init(domain, err)
      use sl_model_mod
      use user_specs_mod, only: nglv
      use spharmt
      use planets_mod
      use io_mod !<<HH delete if not

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object
      type (block_type), pointer :: block         !< Input/output: Block object
      type (mpas_pool_type), pointer :: geometryPool !< mesh information

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography
		real (kind=RKIND), dimension(:), allocatable :: mask_mesh
      real (kind=RKIND), dimension(nglv,2*nglv) :: ism_iceload, ism_bedrock, ism_mask
		real (kind=RKIND), dimension(nglv*2*nglv) :: mask_slGrid1D
      real (kind=RKIND), dimension(nglv*2*nglv) :: thickness_slGrid1D
      real (kind=RKIND), dimension(nglv*2*nglv) :: bedrock_slGrid1D

      integer :: unit_num
      integer :: slmTimeStep, itersl_sh, dtime_sh !get these values from namelist
      real    :: starttime_sh

      ! initialize coupling time step number. initial time is 0
      slmTimeStep = 0  !<<<HH: hardcode it until making the namelist file
      itersl_sh = 1      ! <<<HH* note: hardcode itersl_sh == 1.
      starttime_sh = 0 ! <<<HH get itersl_sh,, starttime_sh (from MALI simulation, needs manipulation)
      dtime_sh = 5
      ! * note: when reading in namelist values, SLM should return dtime from its own namelist file and MALI stores the values

      ! initialize interpolation 
      call interpolate_init(domain, err)

      ! Allocate globalArray and gatheredArray only on process 0
      call MPI_COMM_RANK( domain % dminfo % comm, curProc, err)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)

      if (curProc.eq.0) then
          allocate(globalArrayIce(nCellsGlobal), gatheredArrayIce(nCellsGlobal))
          allocate(globalArrayBedTopo(nCellsGlobal), gatheredArrayBedTopo(nCellsGlobal))
			 allocate(mask_mesh(nCellsGlobal))
      endif

      ! Gather only the nCellsOwned from thickness and bedtopo (does not include Halos)
      call MPI_GATHERV(thickness, nCellsOwned, MPI_DOUBLE, gatheredArrayIce, nCellsPerProc, &
                       nCellsDisplacement, MPI_DOUBLE, 0, domain % dminfo % comm, err)
      call MPI_GATHERV(bedTopography, nCellsOwned, MPI_DOUBLE, gatheredArrayBedTopo, nCellsPerProc, &
                       nCellsDisplacement, MPI_DOUBLE, 0, domain % dminfo % comm, err)

      if (curProc.eq.0) then

         ! Rearrange data into CellID order
         do iCell = 1,nCellsGlobal
            globalArrayIce(indexToCellIDGathered(iCell)) = gatheredArrayIce(iCell)
            globalArrayBedTopo(indexToCellIDGathered(iCell)) = gatheredArrayBedTopo(iCell)
				mask_mesh(indexToCellIDGathered(iCell)) = 1
         enddo

         ! interpolate thickness, bedTopograpy, mesh mask to the Gaussian grid
         call interpolate(toColValues, toRowValues, toSvalues, globalArrayIce, thickness_slGrid1D)
         call interpolate(toColValues, toRowValues, toSvalues, globalArrayBedTopo, bedrock_slGrid1D)
			call interpolate(toColValues, toRowValues, toSvalues, mask_mesh, mask_slGrid1D)

         ! reformat the interpolated data
         ism_iceload = reshape(thickness_slGrid1D, [nglv,2*nglv])
         ism_bedrock = reshape(bedrock_slGrid1D, [nglv, 2*nglv])
         ism_mask = reshape(mask_slGrid1D, [nglv, 2*nglv])

         ! define unit_num for writing output messages
         unit_num = domain % logInfo % outputLog % unitNum

         call sl_set_unit_num(unit_num)
         call sl_solver_checkpoint(itersl_sh, dtime_sh) !Uncomment this after fixing the spharmt error
         call sl_timewindow(slmTimeStep)
         call sl_solver_init(itersl_sh, starttime_sh, ism_iceload, ism_bedrock, ism_mask)
         call deallocate_slarrays ! Deallocate arrays that are set from sl_timewindow()

         deallocate(globalArrayIce, gatheredArrayIce)
         deallocate(globalArrayBedTopo, gatheredArrayBedTopo)
			deallocate(mask_mesh)
      endif

   !--------------------------------------------------------------------
   end subroutine slmodel_init



!***********************************************************************
!
!  routine slmodel_solve
!
!> \brief   Solves gravitationally consistent sea-level change
!> \author  Holly Kyeore Han
!> \date    January 2022
!> \details
!>  This wrapper routine calls the sea-level solver that takes in
!>  ice thickness and provides sea-level change (i.e., changes in the
!>  heights of the sea surface and the solid Earth surface associated
!>  with ice sheet changes. The sea-levle model is taken and modified
!>  from Han et al. (2021, GMD, https://doi.org/10.5281/zenodo.5775235)
!
!-----------------------------------------------------------------------

   subroutine slmodel_solve(slmTimeStep, domain)
      use li_advection
      use sl_model_mod
      use user_specs_mod, only: nglv
      use io_mod !<<HH delete if not needed

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object
	   integer, intent(in) :: slmTimeStep

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: meshPool          !< mesh information
      type (mpas_pool_type), pointer :: geometryPool      !< geometry information
      type (mpas_pool_type), pointer :: velocityPool      !< velocity information

      real (kind=RKIND), dimension(:), pointer :: bedTopography, thickness
		real (kind=RKIND), dimension(:), pointer :: upliftDiff

      real (kind=RKIND), dimension(:), allocatable :: slChange_maliMesh, mask_mesh
      real (kind=RKIND), dimension(nglv,2*nglv) :: ism_iceload, ism_mask
	   real (kind=RKIND), dimension(nglv,2*nglv) :: slm_slchange
      real (kind=RKIND), dimension(nglv*2*nglv) :: thickness_slGrid1D, slChange_slGrid1D
		real (kind=RKIND), dimension(nglv*2*nglv) :: mask_slGrid1D

		integer :: err_tmp

      integer ::  itersl_sh, dtime_sh
      real    :: starttime_sh

		err_tmp = 0
      itersl_sh = 1      ! * note: hardcode itersl_sh == 1.
      starttime_sh = 0 ! get itersl_sh,, starttime_sh (from MALI simulation, needs manipulation)
      dtime_sh = 5

      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'upliftDiff', upliftDiff)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'velocity', velocityPool)

      ! Allocate globalArray and gatheredArray only on process 0
      call MPI_COMM_RANK(domain % dminfo % comm, curProc, err)

      if (curProc.eq.0) then
         allocate(globalArrayIce(nCellsGlobal), gatheredArrayIce(nCellsGlobal))
         allocate(globalArrayBedTopo(nCellsGlobal), gatheredArrayBedTopo(nCellsGlobal))
         allocate(globalArraySLC(nCellsGlobal), gatheredArraySLC(nCellsGlobal))
			allocate(mask_mesh(nCellsGlobal))
      endif

      ! Gather only the nCellsOwned from ice thickness (does not include Halos)
      call MPI_GATHERV(thickness, nCellsOwned, MPI_DOUBLE, gatheredArrayIce, nCellsPerProc, &
                       nCellsDisplacement, MPI_DOUBLE, 0, domain % dminfo % comm, err)

      if (curProc.eq.0) then

         ! Rearrange thickness into CellID order
         do iCell = 1,nCellsGlobal
            globalArrayIce(indexToCellIDGathered(iCell)) = gatheredArrayIce(iCell)
				mask_mesh(indexToCellIDGathered(iCell)) = 1
         enddo

         ! interpolate thickness to Gaussian grid
         call interpolate(toColValues, toRowValues, toSvalues, globalArrayIce, thickness_slGrid1D)
			call interpolate(toColValues, toRowValues, toSvalues, mask_mesh, mask_slGrid1D)
         call write_txt(thickness_slGrid1D, 'mali_iceload','')

         ! reshape the interpolated data
         ism_iceload(:,:) = reshape(thickness_slGrid1D, [nglv,2*nglv])
			ism_mask = reshape(mask_slGrid1D, [nglv, 2*nglv])

		   ! call the sea-level model
		   call sl_solver_checkpoint(itersl_sh, dtime_sh) !Uncomment this after fixing the spharmt boundcheck error
		   call sl_timewindow(slmTimeStep)
		   call sl_solver(itersl_sh, slmTimeStep, dtime_sh, starttime_sh, ism_iceload, ism_mask, slm_slchange)
		   call deallocate_slarrays

         ! interpolate slchange from Gaussian grid to MALI mesh and reshape
         slChange_slGrid1D = reshape(slm_slchange, [nglv*2*nglv])

         ! interpolate sea-level change from GL grid to MALI mesh.
         call interpolate(fromColValues, fromRowValues, fromSValues, slChange_slGrid1D, globalArraySLC)
          
         ! Rearrange back to index order
         do iCell = 1,nCellsGlobal
             gatheredArraySLC(iCell) = globalArraySLC(indexToCellIDGathered(iCell))
         enddo
            
      endif

      ! scatter bedTopography and sea-level change to processors
      call MPI_SCATTERV(gatheredArrayBedTopo, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
                  bedTopography, nCellsAll, MPI_DOUBLE, 0, domain % dminfo % comm, err) 
						
      call MPI_SCATTERV(gatheredArraySLC, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
                  upliftDiff, nCellsAll, MPI_DOUBLE, 0, domain % dminfo % comm, err) 

      ! update bedTopography
      bedTopography(:) = bedTopography(:) - upliftDiff(:)

      call li_update_geometry(geometryPool)
      call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)

      ! deallocate memory
      if (curProc.eq.0) then
         deallocate(globalArrayIce, gatheredArrayIce)
         deallocate(globalArrayBedTopo, gatheredArrayBedTopo)
         deallocate(globalArraySLC, gatheredArraySLC)
			deallocate(mask_mesh)
      endif

      ! Perform Halo exchange update
      call mpas_dmpar_field_halo_exch(domain,'bedTopography')
      call mpas_dmpar_field_halo_exch(domain,'upliftDiff')

   !--------------------------------------------------------------------
   end subroutine slmodel_solve



!***********************************************************************
!
!  routine interpolate
!
!> \brief   Perform interpolation between MALI mesh and SLM grid
!> \author  Holly Han
!> \date    December 2021
!> \details
!>  This routine contains the sparse matrix multiplication
!>  algorithm to interpolate between MPAS and Gaussian Grid.
!>  Note: This routine is an exact copy of the inerpolation 
!>  routine written by Kristin Barton in the code
!>  mpas_ocn_vel_self_attraction_loading.F in MPAS-Ocean 
!
!-----------------------------------------------------------------------

   subroutine interpolate(colValues, rowValues, sValues, dataIn, dataOut)
      use io_mod
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      integer, dimension(:) :: rowValues, colValues
      real (kind=RKIND), dimension(:) :: sValues, dataIn

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:) :: dataOut

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real (kind=RKIND) :: rhs = 0
      integer :: n_S, n, nRow, nCol

      n_S = size(sValues)
      n = 1
      do while (n .LE. n_S)
          nRow = rowValues(n)
          do while ( rowValues(n) .EQ. nRow )
              nCol = colValues(n)
              rhs = rhs + dataIn(nCol) * sValues(n)
              n = n + 1
          end do
          dataOut(nRow) = rhs
          rhs = 0
      end do

   !--------------------------------------------------------------------
   end subroutine interpolate



!***********************************************************************
!
!  routine interpolate_init
!
!> \brief   Sets up interpolation between MALI and SLM native grids
!> \author  Holly Han
!> \date    December 2021
!> \details
!>  This routine reads in map (weight) files needed to interpolate
!>  values of ice thicknesss, bedTopography, sea-level change
!>  between native grid of MALI (unstructured) and SLM (Gaussian).
!>  It also gathers and scatters data from and to multiple processors.
!>  Note: A big portion of the routine is copied from routine
!>  'ocn_vel_self_attraction_loading_init' written by Kristin Barton
!>  in the code smpas_ocn_vel_self_attraction_loading.F in MPAS-Ocean
!
!-----------------------------------------------------------------------

   subroutine interpolate_init(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

   
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain    !< Input/output: Domain
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------


      character (len=StrKIND), pointer :: config_MALI_to_SLM_weights_file
      character (len=StrKIND), pointer :: config_SLM_to_MALI_weights_file
      type (mpas_pool_type), pointer :: meshPool          !< mesh information

      ! NetCDF and weights file variables
      integer :: toNcId, toNsDimId, toRowId, toColId, toSId
      integer :: fromNcId, fromNsDimId, fromRowId, fromColId, fromSId
      integer:: nMpasDimId, nGridDimId, toNsLen, fromNsLen
      character (len = NF90_MAX_NAME) :: toNsName, fromNsName, nMpasName, nGridName
      integer, pointer :: n_s
      character(len=StrKIND) :: mpasToGridFile, gridToMpasFile

      ! MPI variables
      integer :: curProc
     ! integer, pointer ::  nCells
      integer, dimension(:), pointer :: indexToCellID
      integer :: iProc, l, ilm, nProcs

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_MALI_to_SLM_weights_file', config_MALI_to_SLM_weights_file)
      call mpas_pool_get_config(liConfigs, 'config_SLM_to_MALI_weights_file', config_SLM_to_MALI_weights_file)
      
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCellsAll)   
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsOwned)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

      ! Begin MPI portion
      call MPI_COMM_RANK( domain % dminfo % comm, curProc, err)
      call MPI_COMM_SIZE( domain % dminfo % comm, nProcs, err)

      
      ! perform the initialization on the head processor
      if (curProc.eq.0) then
          allocate(nCellsPerProc(nProcs))
         allocate(nCellsDisplacement(nProcs))
      endif
      
      ! Gather nCellsOwned
      call MPI_GATHER( nCellsOwned, 1, MPI_INTEGER, nCellsPerProc, 1, MPI_INTEGER, &
                       0, domain % dminfo % comm, err)
      
      ! Set Displacement variable for GATHERV command
      if (curProc.eq.0) then
         nCellsGlobal = sum(nCellsPerProc)
         allocate(indexToCellIDGathered(nCellsGlobal))
         nCellsDisplacement(1) = 0
         if (nProcs > 1) then
            do iProc=2,nProcs
               nCellsDisplacement(iProc) = nCellsDisplacement(iProc-1) + nCellsPerProc(iProc-1)
            enddo
         endif
      endif

      ! Gather indexToCellID
      call MPI_GATHERV( indexToCellID, nCellsOwned, MPI_INTEGER, indexToCellIDGathered, &
              nCellsPerProc, nCellsDisplacement, MPI_INTEGER, 0, domain % dminfo % comm, err)

      !initialize interpoation
       if (curProc.eq.0) then
         mpasToGridFile = trim(config_MALI_to_SLM_weights_file)
         gridToMpasFile = trim(config_MALI_to_SLM_weights_file)

         ! Open netcdf weights files
         call check( nf90_open(path = mpasToGridFile, mode = nf90_nowrite, ncid = toNcId) ,err)
         call check( nf90_open(path = gridToMpasFile, mode = nf90_nowrite, ncid = fromNcId) ,err)

         ! Get dimension ID
         call check( nf90_inq_dimid(fromNcId, "n_s", fromNsDimId) ,err)
         call check( nf90_inq_dimid(toNcId, "n_a", nMpasDimId) ,err)
         call check( nf90_inq_dimid(toNcId, "n_s", toNsDimId) ,err)
         call check( nf90_inq_dimid(fromNcId, "n_a", nGridDimId) ,err)

         ! Get Variable IDs
         call check( nf90_inq_varid(toNcId, "row", toRowId) ,err)
         call check( nf90_inq_varid(toNcId, "col", toColId) ,err)
         call check( nf90_inq_varid(toNcId, "S", toSId) ,err)
         call check( nf90_inq_varid(fromNcId, "row", fromRowId) ,err)
         call check( nf90_inq_varid(fromNcId, "col", fromColId) ,err)
         call check( nf90_inq_varid(fromNcId, "S", fromSId) ,err)

         ! Get Dimension Length
         call check( nf90_inquire_dimension(toNcId, toNsDimId, toNsName, toNsLen) ,err)
         call check( nf90_inquire_dimension(fromNcId, fromNsDimId, fromNsName, fromNsLen) ,err)
         call check( nf90_inquire_dimension(toNcId, nMpasDimId, nMpasName, nMpas) ,err)
         call check( nf90_inquire_dimension(fromNcId, nGridDimId, nGridName, nGrid) ,err)

         ! Allocate matrices to read data into
         allocate ( toRowValues (toNsLen) )
         allocate ( toColValues (toNsLen) )
         allocate ( toSValues (toNsLen) )
         allocate ( fromRowValues (fromNsLen) )
         allocate ( fromColValues (fromNsLen) )
         allocate ( fromSValues (fromNsLen) )

         ! Retrieve data
         call check( nf90_get_var(toNcId, toColId, toColValues(:) ) ,err)
         call check( nf90_get_var(toNcId, toRowId, toRowValues(:) ) ,err)
         call check( nf90_get_var(toNcId, toSId, toSValues(:) ) ,err)
         call check( nf90_get_var(fromNcId, fromColId, fromColValues(:) ) ,err)
         call check( nf90_get_var(fromNcId, fromRowId, fromRowValues(:) ) ,err)
         call check( nf90_get_var(fromNcId, fromSId, fromSValues(:) ) ,err)

      endif

   !--------------------------------------------------------------------
   end subroutine interpolate_init



!***********************************************************************
!
!  routine check
!
!> \brief   Check status of netcdf operations
!> \author  Holly Han
!> \date    December 2021
!> \details
!>  This routine checks to status of the netcdf file
!
!-----------------------------------------------------------------------

   subroutine check(status, err)
      integer, intent ( in) :: status
      integer, intent(inout) :: err

      if(status /= nf90_noerr) then
         err = 1
      end if
   !--------------------------------------------------------------------   
   end subroutine check


!***********************************************************************

end module li_bedtopo

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
