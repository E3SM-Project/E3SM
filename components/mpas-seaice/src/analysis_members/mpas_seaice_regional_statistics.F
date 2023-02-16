!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_regional_statistics
!
!> \brief MPAS sea ice analysis mode member: regional_statistics
!> \author Adrian K. Turner
!> \date   6th September 2015
!> \details
!>  MPAS sea ice analysis mode member: regional_statistics
!>  calculates regional statistics
!>
!-----------------------------------------------------------------------

module seaice_regional_statistics

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_log, only: mpas_log_write

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

   public :: seaice_bootstrap_regional_statistics, &
             seaice_init_regional_statistics, &
             seaice_precompute_regional_statistics, &
             seaice_compute_regional_statistics, &
             seaice_restart_regional_statistics, &
             seaice_finalize_regional_statistics

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   ! reduction function interface
   abstract interface

      subroutine reduction_interface(&
           domain, &
           outputFieldName, &
           inputFieldName)

        use mpas_derived_types

        type(domain_type), intent(inout) :: &
             domain

        character(len=*), intent(in) :: &
             outputFieldName, &
             inputFieldName

      end subroutine reduction_interface

   end interface

   ! reduction type
   type :: reduction_type

      ! reduction name
      character(len=strKIND) :: reductionName

      ! reduction subroutines
      procedure (reduction_interface), pointer, nopass :: reduction_ptr_1D => null ()
      procedure (reduction_interface), pointer, nopass :: reduction_ptr_2D => null ()

      ! linked list pointer
      type(reduction_type), pointer :: next => null()

   end type reduction_type

   ! reduction types linked list head
   type(reduction_type), pointer :: reductionHead

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_regional_statistics
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    10th November 2015
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_regional_statistics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine seaice_bootstrap_regional_statistics!}}}

!***********************************************************************
!
!  routine seaice_init_regional_statistics
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_regional_statistics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      ! check if regions have been read in
      if (region_masks_are_empty(domain)) then

         ! initialize hemispheric regions
         call init_hemisphere_regions(domain)

      endif

      ! init vertex masks
      call init_vertex_masks(domain)

      ! initialize the runtime regional statistics system
      call init_runtime_regional_statistics(domain)

   end subroutine seaice_init_regional_statistics!}}}

!***********************************************************************
!
!  routine region_masks_are_empty
!
!> \brief   Determine of the region mask is all zeros
!> \author  Adrian K. Turner
!> \date    7th April 2016
!> \details This function returns a logical determining whether the
!>   region mask is full of zeros, essentially returning whether an
!>   external region file has been read in at initialization.
!
!-----------------------------------------------------------------------

   function region_masks_are_empty(domain) result(masksAreEmpty)

     type(domain_type), intent(inout) :: &
          domain

     logical :: &
          masksAreEmpty

     type(block_type), pointer :: &
          block

     type(MPAS_pool_type), pointer :: &
          regionsPool

     integer, dimension(:,:), pointer :: &
          regionCellMasks

     integer, pointer :: &
          nCellsSolve

     integer :: &
          nCellsMasked, &
          nCellsMaskedGlobal

     nCellsMasked = 0

     block => domain % blocklist
     do while (associated(block))

        call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

        call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

        call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

        nCellsMasked = nCellsMasked + sum(regionCellMasks(:,1:nCellsSolve))

        block => block % next
     enddo

     call MPAS_dmpar_sum_int(domain % dminfo, nCellsMasked, nCellsMaskedGlobal)

     if (nCellsMaskedGlobal == 0) then
        masksAreEmpty = .true.
     else
        masksAreEmpty = .false.
     endif

   end function region_masks_are_empty

!***********************************************************************
!
!  routine init_hemisphere_regions
!
!> \brief   Calculate hemispheric regions
!> \author  Adrian K. Turner
!> \date    7th April 2016
!> \details If no region file is read in at initialization we calculate
!>   three default regions: global, northern hemisphere and southern
!>   hemisphere.
!
!-----------------------------------------------------------------------

   subroutine init_hemisphere_regions(domain)

     use seaice_constants, only: &
          seaiceDegreesToRadians

     type(domain_type), intent(inout) :: &
          domain

     type(block_type), pointer :: &
          block

     type(MPAS_pool_type), pointer :: &
          regionsPool, &
          meshPool

     integer, dimension(:,:), pointer :: &
          regionCellMasks

     character(len=strKIND), dimension(:), pointer :: &
          regionNames

     real(kind=RKIND), dimension(:), pointer :: &
          latCell

     integer, pointer :: &
          nCells

     integer :: &
          iCell

     real(kind=RKIND), parameter :: &
          northernHemisphereLatitudeLimit =  40.0_RKIND, &
          southernHemisphereLatitudeLimit = -40.0_RKIND

     block => domain % blocklist
     do while (associated(block))

        call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)
        call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

        call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)
        call MPAS_pool_get_array(regionsPool, "regionNames", regionNames)

        call MPAS_pool_get_array(meshPool, "latCell", latCell)

        call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

        ! set region names
        regionNames(1) = "Global"
        regionNames(2) = "Northern Hemisphere"
        regionNames(3) = "Southern Hemisphere"

        do iCell = 1, nCells

           regionCellMasks(:,iCell) = 0

           ! global
           if (latCell(iCell) >= northernHemisphereLatitudeLimit * seaiceDegreesToRadians .or. &
               latCell(iCell) <= southernHemisphereLatitudeLimit * seaiceDegreesToRadians) then
              regionCellMasks(1,iCell) = 1
           endif

           ! northern hemisphere
           if (latCell(iCell) >= northernHemisphereLatitudeLimit * seaiceDegreesToRadians) then
              regionCellMasks(2,iCell) = 1
           endif

           ! southern hemisphere
           if (latCell(iCell) <= southernHemisphereLatitudeLimit * seaiceDegreesToRadians) then
              regionCellMasks(3,iCell) = 1
           endif

        enddo ! iCell

        block => block % next
     enddo

   end subroutine init_hemisphere_regions

!***********************************************************************
!
!  routine init_vertex_masks
!
!> \brief   Determine of the region mask is all zeros
!> \author  Adrian K. Turner
!> \date    13th May 2016
!> \details Calculate the vertex masks from the cell masks
!
!-----------------------------------------------------------------------

   subroutine init_vertex_masks(domain)

     type(domain_type), intent(inout) :: &
          domain

     type(block_type), pointer :: &
          block

     type(MPAS_pool_type), pointer :: &
          regionsPool, &
          meshPool

     integer, dimension(:,:), pointer :: &
          regionCellMasks, &
          regionVertexMasks, &
          cellsOnVertex

     integer, pointer :: &
          nVerticesSolve, &
          nRegions, &
          vertexDegree

     integer :: &
          iVertex, &
          iRegion, &
          iCellOnVertex

     block => domain % blocklist
     do while (associated(block))

        call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)
        call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

        call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)
        call MPAS_pool_get_array(regionsPool, "regionVertexMasks", regionVertexMasks)

        call MPAS_pool_get_array(meshPool, "cellsOnVertex", cellsOnVertex)

        call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
        call MPAS_pool_get_dimension(block % dimensions, "nRegions", nRegions)
        call MPAS_pool_get_dimension(block % dimensions, "vertexDegree", vertexDegree)

        do iVertex = 1, nVerticesSolve
           do iRegion = 1, nRegions

              regionVertexMasks(iRegion, iVertex) = 0

              do iCellOnVertex = 1, vertexDegree

                 if (regionCellMasks(iRegion,cellsOnVertex(iCellOnVertex,iVertex)) == 1) then
                    regionVertexMasks(iRegion, iVertex) = 1
                 endif

              enddo ! iCellOnVertex

           enddo ! iRegion
        enddo ! iVertex

        block => block % next
     enddo

   end subroutine init_vertex_masks

!***********************************************************************
!
!  routine init_runtime_regional_statistics
!
!> \brief   Initialize runtime configurable regional statistics
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine iterates through the input regional statistics
!>  streams and creates an output field based on the input field
!>  parameters. This field is also added to the output stream.
!
!-----------------------------------------------------------------------

  subroutine init_runtime_regional_statistics(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(reduction_type), pointer :: &
         reduction

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool

    type(MPAS_pool_field_info_type) :: &
         inputFieldInfo

    character(len=strKIND) :: &
         streamID, &
         inputFieldName, &
         outputFieldName

    type(field1DReal), pointer :: &
         inputField1D, &
         outputField1D

    type(field2DReal), pointer :: &
         inputField2D, &
         outputField2D

    integer :: &
         ierrField, &
         ierrStream

    integer, pointer :: &
         nRegions

    logical :: &
         isActive

    ! initialize the reduction type linked list
    call init_reduction_linked_list()

    ! loop over blocks
    block => domain % blocklist
    blockLoop: do while (associated(block))

       ! get the region number
       call MPAS_pool_get_dimension(block % dimensions, "nRegions", nRegions)

       ! create a new pool for derived fields
       call MPAS_pool_create_pool(derivedFieldPool)

       ! loop over reduction types
       reduction => reductionHead
       reductionLoop: do while (associated(reduction))

          ! get the stream id
          streamID = "regionalStatistics_"//trim(reduction % reductionName)
          !call mpas_log_write("init_runtime_regional_statistics: streams: "//trim(streamID))

          ! test stream exists
          if (MPAS_stream_mgr_stream_exists(domain % streamManager, trim(streamID))) then

             ! loop over entries in the input stream
             call MPAS_stream_mgr_begin_iteration(domain % streamManager, trim(streamID), ierr=ierrField)
             streamLoop: do while (MPAS_stream_mgr_get_next_field(&
                  domain % streamManager, trim(streamID), inputFieldName, isActive))

                !call mpas_log_write("init_runtime_regional_statistics: fields: "//trim(inputFieldName))

                ! check the input field is of the correct type - must be real(nCells), single time level
                call MPAS_pool_get_field_info(block % allFields, trim(inputFieldName), inputFieldInfo)
                if (inputFieldInfo % fieldType   /= MPAS_POOL_REAL .or. &
                    inputFieldInfo % nDims       >  2 .or. &
                    inputFieldInfo % nTimeLevels /= 1) then
                   call mpas_log_write(&
                        "init_runtime_regional_statistics: input field wrong type: "//trim(inputFieldName), &
                        MPAS_LOG_CRIT)
                endif

                !call mpas_log_write("init_runtime_regional_statistics: field active: $l, nDims: $i ", intArgs=(/inputFieldInfo % nDims/), logicArgs=(/inputFieldInfo % isActive/))

                ! check if field is active
                if (inputFieldInfo % isActive) then

                   ! set the output fieldname
                   outputFieldName = "regionalStatistics_"//trim(inputFieldName)//"_"//trim(reduction % reductionName)

                   !call mpas_log_write("init_runtime_regional_statistics: outputFieldName: "//trim(outputFieldName))

                   ! select array rank
                   if (inputFieldInfo % nDims == 1) then

                      ! get the input field
                      call MPAS_pool_get_field(block % allFields, trim(inputFieldName), inputField1D)
                      if (trim(inputField1D % dimNames(1)) /= "nCells") then
                         call mpas_log_write(&
                              "init_runtime_regional_statistics: input field wrong dimensions: "//&
                              trim(inputFieldName), MPAS_LOG_CRIT)
                      endif

                      ! create a new field
                      call create_new_output_field_1D(outputField1D, inputField1D, block, nRegions, outputFieldName)

                      ! add the new field to the output pool
                      call MPAS_pool_add_field(derivedFieldPool,  trim(outputFieldName), outputField1D)
                      call MPAS_pool_add_field(block % allFields, trim(outputFieldName), outputField1D)

                   else if (inputFieldInfo % nDims == 2) then

                      ! get the input field
                      call MPAS_pool_get_field(block % allFields, trim(inputFieldName), inputField2D)
                      if (trim(inputField2D % dimNames(2)) /= "nCells") then
                         call mpas_log_write(&
                              "init_runtime_regional_statistics: input field wrong dimensions: "//&
                              trim(inputFieldName), MPAS_LOG_CRIT)
                      endif

                      ! create a new field
                      call create_new_output_field_2D(outputField2D, inputField2D, block, nRegions, outputFieldName)

                      ! add the new field to the output pool
                      call MPAS_pool_add_field(derivedFieldPool,  trim(outputFieldName), outputField2D)
                      call MPAS_pool_add_field(block % allFields, trim(outputFieldName), outputField2D)

                   endif

                   ! add the field to the output stream
                   call MPAS_stream_mgr_add_field(domain % streamManager, "regionalStatisticsOutput", outputFieldName)

                endif ! isActive

             enddo streamLoop

          endif ! stream exists

          reduction => reduction % next
       enddo reductionLoop

       ! add the new pool to the block pool list
       call MPAS_pool_add_subpool(block % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)

       block => block % next
    enddo blockLoop

  end subroutine init_runtime_regional_statistics

!***********************************************************************
!
!  routine init_reduction_linked_list
!
!> \brief   Initialize the reduction linked list
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine adds the defined reduction types to the
!>  reduction linked list.
!
!-----------------------------------------------------------------------

  subroutine init_reduction_linked_list()

    call add_reduction(reductionHead, "ice_aggregate",      ice_aggregate_1D,      ice_aggregate_2D)
    call add_reduction(reductionHead, "cell_aggregate",     cell_aggregate_1D,     cell_aggregate_2D)
    call add_reduction(reductionHead, "ice_areal_average",  ice_areal_average_1D,  ice_areal_average_2D)
    call add_reduction(reductionHead, "cell_areal_average", cell_areal_average_1D, cell_areal_average_2D)
    call add_reduction(reductionHead, "min",                min_1D,                min_2D)
    call add_reduction(reductionHead, "max",                max_1D,                max_2D)

  end subroutine init_reduction_linked_list

!***********************************************************************
!
!  routine add_reduction
!
!> \brief   Add a reduction type to the reduction linked list
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine adds an individual reduction type to the
!>  reduction linked list. Required information is the reduction name,
!>  which is used to generate the output field name, and function
!>  pointers to the subroutines that will perform the reduction.
!
!-----------------------------------------------------------------------

  subroutine add_reduction(&
       reductionHead, &
       reductionName, &
       reductionSubroutine1D, &
       reductionSubroutine2D)

    type(reduction_type), pointer :: &
         reductionHead

    character(len=*), intent(in) :: &
         reductionName

    interface
      subroutine reductionSubroutine1D(&
           domain, &
           outputFieldName, &
           inputFieldName)

        use mpas_derived_types

        type(domain_type), intent(inout) :: &
             domain

        character(len=*), intent(in) :: &
             outputFieldName, &
             inputFieldName

      end subroutine reductionSubroutine1D
   end interface

   interface
      subroutine reductionSubroutine2D(&
           domain, &
           outputFieldName, &
           inputFieldName)

        use mpas_derived_types

        type(domain_type), intent(inout) :: &
             domain

        character(len=*), intent(in) :: &
             outputFieldName, &
             inputFieldName

      end subroutine reductionSubroutine2D
   end interface

    type(reduction_type), pointer :: &
         reductionNew

    ! loop through the linked list
    if (.not. associated(reductionHead)) then
       allocate(reductionHead)
       nullify(reductionHead % next)
       reductionNew => reductionHead
    else
       reductionNew => reductionHead
       do while (associated(reductionNew % next))
          if (trim(reductionNew % reductionName) == trim(reductionName)) then
             call mpas_log_write('MPAS-seaice: add_reduction', MPAS_LOG_CRIT)
          endif
          reductionNew => reductionNew % next
       enddo
       allocate(reductionNew % next)
       reductionNew => reductionNew % next
       nullify(reductionNew % next)
    endif

    ! now actually add the reduction
    ! reduction name
    reductionNew % reductionName = trim(reductionName)

    ! reduction subroutine pointers
    reductionNew % reduction_ptr_1D => reductionSubroutine1D
    reductionNew % reduction_ptr_2D => reductionSubroutine2D

  end subroutine add_reduction

!***********************************************************************
!
!  routine create_new_output_field_1D
!
!> \brief   Creates a new 1D output field
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine creates a regional statistics output field
!>  based on the input field. This subroutine operates on 1D (nCells)
!>  input fields only
!
!-----------------------------------------------------------------------

  subroutine create_new_output_field_1D(&
       outputField, &
       inputField1D, &
       block, &
       nRegions, &
       fieldName)

    type(field1DReal), pointer :: &
         outputField, &
         inputField1D

    type(block_type), pointer :: &
         block

    integer, intent(in) :: &
         nRegions

    character(len=*), intent(in) :: &
         fieldName

    ! allocate the actual field
    allocate(outputField)

    ! set block back pointer
    outputField % block => block

    ! allocate the array
    allocate(outputField % array(nRegions))

    ! initialize the array
    outputField % array = 0.0_RKIND

    ! fieldname
    outputField % fieldName       = trim(fieldName)
    outputField % outputFieldName = trim(fieldName)

    ! I/O layer
    outputField % dimNames(1)      = "nRegions"
    outputField % dimSizes(1)      = nRegions
    outputField % defaultValue     = 0.0_RKIND
    outputField % isDecomposed     = .false.
    outputField % hasTimeDimension = .true.
    outputField % isActive         = .true.
    outputField % isVarArray       = .false.
    outputField % isPersistent     = .true.

    allocate(outputField % attLists(1))

  end subroutine create_new_output_field_1D

!***********************************************************************
!
!  routine create_new_output_field_2D
!
!> \brief   Creates a new 2D output field
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine creates a regional statistics output field
!>  based on the input field. This subroutine operates on 2D (:,nCells)
!>  input fields only
!
!-----------------------------------------------------------------------

  subroutine create_new_output_field_2D(&
       outputField, &
       inputField, &
       block, &
       nRegions, &
       fieldName)

    type(field2DReal), pointer :: &
         outputField, &
         inputField

    type(block_type), pointer :: &
         block

    integer, intent(in) :: &
         nRegions

    character(len=*), intent(in) :: &
         fieldName

    ! allocate the actual field
    allocate(outputField)

    ! set block back pointer
    outputField % block => block

    ! allocate the array
    allocate(outputField % array(size(inputField % array,1),nRegions))

    ! initialize the array
    outputField % array = 0.0_RKIND

    ! fieldname
    outputField % fieldName       = trim(fieldName)
    outputField % outputFieldName = trim(fieldName)

    ! I/O layer
    outputField % dimNames(1)      = inputField % dimNames(1)
    outputField % dimSizes(1)      = inputField % dimSizes(1)
    outputField % dimNames(2)      = "nRegions"
    outputField % dimSizes(2)      = nRegions
    outputField % defaultValue     = 0.0_RKIND
    outputField % isDecomposed     = .false.
    outputField % hasTimeDimension = .true.
    outputField % isActive         = .true.
    outputField % isVarArray       = .false.
    outputField % isPersistent     = .true.

    allocate(outputField % attLists(1))

  end subroutine create_new_output_field_2D

!***********************************************************************
!
!  routine seaice_precompute_regional_statistics
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_regional_statistics(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      err = 0

   end subroutine seaice_precompute_regional_statistics!}}}

!***********************************************************************
!
!  routine seaice_compute_regional_statistics
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_regional_statistics(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

      err = 0

      ! fixed regional stats
      call fixed_regional_statistics(domain)

      ! runtime regional statistics
      call runtime_regional_statistics(domain)

   end subroutine seaice_compute_regional_statistics!}}}

!***********************************************************************
!
!  routine fixed_regional_statistics
!
!> \brief   Calculate the non-runtime configurable regional statistics
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routines calculates the standard fixed regional
!>  statistics during the time stepping.
!
!-----------------------------------------------------------------------

   subroutine fixed_regional_statistics(domain)

      use seaice_mesh, only: &
           seaice_interpolate_vertex_to_cell

      use ice_constants_colpkg, only: &
           awtvdr, &
           awtidr, &
           awtvdf, &
           awtidf, &
           rhoi, &
           rhos

      type(domain_type), intent(inout) :: &
          domain

      type (dm_info) :: &
           dminfo

      type (block_type), pointer :: &
           block

      type (mpas_pool_type), pointer :: &
           meshPool, &
           regionsPool, &
           tracersAggregatePool, &
           shortwavePool, &
           velocitySolverPool, &
           boundaryPool, &
           regionalStatisticsAMPool

      integer, pointer :: &
           nRegions, &
           nCellsSolve, &
           nVerticesSolve

      integer :: &
           iRegion, &
           iCell, &
           iVertex

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           iceAreaCell, &
           iceVolumeCell, &
           snowVolumeCell, &
           solarZenithAngleCosine, &
           albedoVisibleDirectCell, &
           albedoIRDirectCell, &
           albedoVisibleDiffuseCell, &
           albedoIRDiffuseCell, &
           icePressure, &
           uVelocity, &
           vVelocity, &
           uVelocityCell, &
           vVelocityCell

      integer, dimension(:,:), pointer :: &
           regionCellMasks, &
           regionVertexMasks

      integer, dimension(:), pointer :: &
           dynamicallyLockedCellsMask

      real(kind=RKIND), dimension(:), pointer :: &
           totalIceArea, &
           totalIceExtent, &
           totalIceVolume, &
           totalSnowVolume, &
           totalKineticEnergy, &
           rmsIceSpeed, &
           averageAlbedo, &
           maximumIceVolume, &
           maximumIceVolumeLocked, &
           maximumIceVolumeNotlocked, &
           maximumIcePressure, &
           maximumIceSpeed

      real(kind=RKIND), dimension(:), allocatable :: &
           globalSumsIn, &
           globalMaxsIn, &
           globalSumsOut, &
           globalMaxsOut

      integer, parameter :: &
           nSums = 8, &
           nMaxs = 5

      integer :: &
           iSum, &
           iMax

      real(kind=RKIND), pointer :: &
           iceExtentLimit

      real(kind=RKIND), parameter :: &
           m2_to_km2 = 1.0e-6_RKIND, &
           m3_to_km3 = 1.0e-9_RKIND, &
           Nm_to_kNm = 1.0e-3_RKIND

      dminfo = domain % dminfo

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

      allocate(globalSumsIn(nRegions*nSums))
      allocate(globalMaxsIn(nRegions*nMaxs))
      globalSumsIn = 0.0_RKIND
      globalMaxsIn = -1.0e30_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_config(block % configs, "config_AM_regionalStatistics_ice_extent_limit", iceExtentLimit)

         call MPAS_pool_get_subpool(block % structs, 'mesh', meshPool)
         call MPAS_pool_get_subpool(block % structs, 'regions', regionsPool)
         call MPAS_pool_get_subpool(block % structs, 'tracers_aggregate', tracersAggregatePool)
         call MPAS_pool_get_subpool(block % structs, 'shortwave', shortwavePool)
         call MPAS_pool_get_subpool(block % structs, 'velocity_solver', velocitySolverPool)
         call MPAS_pool_get_subpool(block % structs, 'boundary', boundaryPool)
         call MPAS_pool_get_subpool(block % structs, 'regionalStatisticsAM', regionalStatisticsAMPool)

         call MPAS_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)

         call MPAS_pool_get_array(meshPool, 'areaCell', areaCell)

         call MPAS_pool_get_array(tracersAggregatePool, 'iceAreaCell', iceAreaCell)
         call MPAS_pool_get_array(tracersAggregatePool, 'iceVolumeCell', iceVolumeCell)
         call MPAS_pool_get_array(tracersAggregatePool, 'snowVolumeCell', snowVolumeCell)

         call MPAS_pool_get_array(shortwavePool, "solarZenithAngleCosine", solarZenithAngleCosine)
         call MPAS_pool_get_array(shortwavePool, "albedoVisibleDirectCell", albedoVisibleDirectCell)
         call MPAS_pool_get_array(shortwavePool, "albedoIRDirectCell", albedoIRDirectCell)
         call MPAS_pool_get_array(shortwavePool, "albedoVisibleDiffuseCell", albedoVisibleDiffuseCell)
         call MPAS_pool_get_array(shortwavePool, "albedoIRDiffuseCell", albedoIRDiffuseCell)

         call MPAS_pool_get_array(velocitySolverPool, "icePressure", icePressure)
         call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
         call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
         call MPAS_pool_get_array(velocitySolverPool, "uVelocityCell", uVelocityCell)
         call MPAS_pool_get_array(velocitySolverPool, "vVelocityCell", vVelocityCell)
         call MPAS_pool_get_array(velocitySolverPool, "dynamicallyLockedCellsMask", dynamicallyLockedCellsMask)

         call MPAS_pool_get_array(regionsPool, 'regionCellMasks', regionCellMasks)
         call MPAS_pool_get_array(regionsPool, 'regionVertexMasks', regionVertexMasks)

         ! cell centre velocity for velocity statistics
         call seaice_interpolate_vertex_to_cell(meshPool, boundaryPool, uVelocityCell, uVelocity)
         call seaice_interpolate_vertex_to_cell(meshPool, boundaryPool, vVelocityCell, vVelocity)

         ! quantities on cells
         do iCell = 1, nCellsSolve
            do iRegion = 1, nRegions

               if (regionCellMasks(iRegion,iCell) == 1) then

                  ! total ice area
                  iSum = (iRegion-1) * nSums + 1
                  globalSumsIn(iSum) = globalSumsIn(iSum) + iceAreaCell(iCell) * areaCell(iCell)

                  ! total ice extent
                  if (iceAreaCell(iCell) > iceExtentLimit) then

                     iSum = (iRegion-1) * nSums + 2
                     globalSumsIn(iSum) = globalSumsIn(iSum) + areaCell(iCell)

                  endif

                  ! total ice volume
                  iSum = (iRegion-1) * nSums + 3
                  globalSumsIn(iSum) = globalSumsIn(iSum) + iceVolumeCell(iCell) * areaCell(iCell)

                  ! total snow volume
                  iSum = (iRegion-1) * nSums + 4
                  globalSumsIn(iSum) = globalSumsIn(iSum) + snowVolumeCell(iCell) * areaCell(iCell)

                  ! kinetic energy
                  iSum = (iRegion-1) * nSums + 5
                  globalSumsIn(iSum) = globalSumsIn(iSum) + 0.5_RKIND * areaCell(iCell) * &
                       (snowVolumeCell(iCell) * rhos + iceVolumeCell(iCell) * rhoi) * &
                       (uVelocityCell(iCell)**2 + vVelocityCell(iCell)**2)

                  ! total mass (for RMS ice speed)
                  iSum = (iRegion-1) * nSums + 6
                  globalSumsIn(iSum) = globalSumsIn(iSum) + areaCell(iCell) * &
                       (snowVolumeCell(iCell) * rhos + iceVolumeCell(iCell) * rhoi)

                  ! average albedo
                  if (solarZenithAngleCosine(iCell) > 0.0_RKIND) then

                     iSum = (iRegion-1) * nSums + 7
                     globalSumsIn(iSum) = globalSumsIn(iSum) + &
                          areaCell(iCell) * &
                          (awtvdr * albedoVisibleDirectCell(iCell) + &
                           awtidr * albedoIRDirectCell(iCell) + &
                           awtvdf * albedoVisibleDiffuseCell(iCell) + &
                           awtidf * albedoIRDiffuseCell(iCell))

                     iSum = (iRegion-1) * nSums + 8
                     globalSumsIn(iSum) = globalSumsIn(iSum) + &
                          areaCell(iCell)

                  endif

                  ! maximum ice volume
                  iMax = (iRegion-1) * nMaxs + 1
                  globalMaxsIn(iMax) = max(globalMaxsIn(iMax), iceVolumeCell(iCell))

                  ! maximum locked ice volume
                  iMax = (iRegion-1) * nMaxs + 2
                  if (dynamicallyLockedCellsMask(iCell) == 1) then
                     globalMaxsIn(iMax) = max(globalMaxsIn(iMax), iceVolumeCell(iCell))
                  endif

                  ! maximum un-locked ice volume
                  iMax = (iRegion-1) * nMaxs + 3
                  if (dynamicallyLockedCellsMask(iCell) == 0) then
                     globalMaxsIn(iMax) = max(globalMaxsIn(iMax), iceVolumeCell(iCell))
                  endif

                  ! maximum ice pressure
                  iMax = (iRegion-1) * nMaxs + 4
                  globalMaxsIn(iMax) = max(globalMaxsIn(iMax), icePressure(iCell))

               endif ! regionCellMasks == 1

            enddo ! iRegion
         enddo ! iCell

         ! quantities on vertices
         do iVertex = 1, nVerticesSolve
            do iRegion = 1, nRegions

               if (regionVertexMasks(iRegion,iVertex) == 1) then

                  ! maximum ice speed
                  iMax = (iRegion-1) * nMaxs + 5
                  globalMaxsIn(iMax) = max(globalMaxsIn(iMax), sqrt(uVelocity(iVertex)**2 + vVelocity(iVertex)**2))

               endif ! regionVertexMasks == 1

            enddo ! iRegion
         enddo ! iVertex

         block => block % next
      enddo

      ! MPI calls
      allocate(globalSumsOut(nRegions*nSums))
      allocate(globalMaxsOut(nRegions*nMaxs))

      call mpas_dmpar_sum_real_array(dminfo, nRegions*nSums, globalSumsIn, globalSumsOut)
      call mpas_dmpar_max_real_array(dminfo, nRegions*nMaxs, globalMaxsIn, globalMaxsOut)

      ! set Registry variables
      call MPAS_pool_get_subpool(domain % blocklist % structs, 'regionalStatisticsAM', regionalStatisticsAMPool)

      call MPAS_pool_get_array(regionalStatisticsAMPool, 'totalIceArea', totalIceArea)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'totalIceExtent', totalIceExtent)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'totalIceVolume', totalIceVolume)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'totalSnowVolume', totalSnowVolume)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'totalKineticEnergy', totalKineticEnergy)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'rmsIceSpeed', rmsIceSpeed)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'averageAlbedo', averageAlbedo)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'maximumIceVolume', maximumIceVolume)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'maximumIceVolumeLocked', maximumIceVolumeLocked)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'maximumIceVolumeNotLocked', maximumIceVolumeNotLocked)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'maximumIcePressure', maximumIcePressure)
      call MPAS_pool_get_array(regionalStatisticsAMPool, 'maximumIceSpeed', maximumIceSpeed)

      do iRegion = 1, nRegions

         totalIceArea(iRegion)              = globalSumsOut((iRegion-1) * nSums + 1)
         totalIceExtent(iRegion)            = globalSumsOut((iRegion-1) * nSums + 2)
         totalIceVolume(iRegion)            = globalSumsOut((iRegion-1) * nSums + 3)
         totalSnowVolume(iRegion)           = globalSumsOut((iRegion-1) * nSums + 4)
         totalKineticEnergy(iRegion)        = globalSumsOut((iRegion-1) * nSums + 5)
         rmsIceSpeed(iRegion)               = globalSumsOut((iRegion-1) * nSums + 5) / &
                                              max(globalSumsOut((iRegion-1) * nSums + 6), 1e-11_RKIND)
         averageAlbedo(iRegion)             = globalSumsOut((iRegion-1) * nSums + 7) / &
                                              max(globalSumsOut((iRegion-1) * nSums + 8), 1e-11_RKIND)
         maximumIceVolume(iRegion)          = globalMaxsOut((iRegion-1) * nMaxs + 1)
         maximumIceVolumeLocked(iRegion)    = globalMaxsOut((iRegion-1) * nMaxs + 2)
         maximumIceVolumeNotLocked(iRegion) = globalMaxsOut((iRegion-1) * nMaxs + 3)
         maximumIcePressure(iRegion)        = globalMaxsOut((iRegion-1) * nMaxs + 4)
         maximumIceSpeed(iRegion)           = globalMaxsOut((iRegion-1) * nMaxs + 5)

      enddo ! iRegion

      ! unit conversions
      totalIceArea    = totalIceArea    * m2_to_km2
      totalIceExtent  = totalIceExtent  * m2_to_km2
      totalIceVolume  = totalIceVolume  * m3_to_km3
      totalSnowVolume = totalSnowVolume * m3_to_km3

      maximumIcePressure = maximumIcePressure * Nm_to_kNm ! N/m -> kN/m

      deallocate(globalSumsOut)
      deallocate(globalMaxsOut)

      deallocate(globalSumsIn)
      deallocate(globalMaxsIn)

   end subroutine fixed_regional_statistics

!***********************************************************************
!
!  routine runtime_regional_statistics
!
!> \brief   Calculate the runtime configurable regional statistics
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine calculates the runtime-configurable regional
!>  statistics during time stepping. It loops through the reduction
!>  types defined in the reduction linked list and for each type
!>  calculates the appropriate statistic for the user defined input
!>  fields.
!
!-----------------------------------------------------------------------

   subroutine runtime_regional_statistics(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(reduction_type), pointer :: &
         reduction

    character(len=strKIND) :: &
         streamID, &
         inputFieldName, &
         outputFieldName

    type(MPAS_pool_field_info_type) :: &
         inputFieldInfo

    integer :: &
         ierrField

    logical :: &
         isActive

    ! loop over reduction types
    reduction => reductionHead
    reductionLoop: do while (associated(reduction))

       ! get the stream id
       streamID = "regionalStatistics_"//trim(reduction % reductionName)

       !call mpas_log_write("runtime_regional_statistics: streamID: "//trim(streamID)//", $l", logicArgs=(/MPAS_stream_mgr_stream_exists(domain % streamManager, trim(streamID))/))

       ! test stream exists
       if (MPAS_stream_mgr_stream_exists(domain % streamManager, trim(streamID))) then

          ! loop over entries in input stream
          call MPAS_stream_mgr_begin_iteration(domain % streamManager, trim(streamID), ierr=ierrField)
          streamLoop: do while (MPAS_stream_mgr_get_next_field(domain % streamManager, trim(streamID), inputFieldName, isActive))

             !call mpas_log_write("runtime_regional_statistics: inputFieldName: "//trim(inputFieldName))

             ! check the input field is active
             call MPAS_pool_get_field_info(domain % blocklist % allFields, trim(inputFieldName), inputFieldInfo)
             if (inputFieldInfo % isActive) then

                ! set the output fieldname
                outputFieldName = "regionalStatistics_"//trim(inputFieldName)//"_"//trim(reduction % reductionName)

                !call mpas_log_write("runtime_regional_statistics: outputFieldName: "//trim(outputFieldName))

                if (inputFieldInfo % nDims == 1) then

                   ! aggregate the input field array
                   call reduction % reduction_ptr_1D(domain, inputFieldName, outputFieldName)

                else if (inputFieldInfo % nDims == 2) then

                   ! aggregate the input field array
                   call reduction % reduction_ptr_2D(domain, inputFieldName, outputFieldName)

                endif

             endif ! field is active

          enddo streamLoop

       endif ! stream exists

       reduction => reduction % next
    enddo reductionLoop

  end subroutine runtime_regional_statistics

!***********************************************************************
!
!  routine ice_aggregate_1D
!
!> \brief   Perform aggregation for a 1D ice defined variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine performs aggregation by region for a 1D
!>  (nCells) field that is defined only over the ice covered part of
!>  the cell.
!
!-----------------------------------------------------------------------

  subroutine ice_aggregate_1D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         meshPool, &
         tracersAggregatePool, &
         regionsPool

    real(kind=RKIND), dimension(:), pointer :: &
         outputFieldArray, &
         inputFieldArray, &
         iceAreaCell, &
         areaCell

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:), allocatable :: &
         outputFieldArrayTmp

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = 0.0_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get the output array
       call MPAS_pool_get_subpool(block % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
       call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                outputFieldArrayTmp(iRegion) = outputFieldArrayTmp(iRegion) + &
                     inputFieldArray(iCell) * &
                     iceAreaCell(iCell) * &
                     areaCell(iCell)

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, outputFieldArrayTmp, outputFieldArray)

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

  end subroutine ice_aggregate_1D

!***********************************************************************
!
!  routine ice_aggregate_2D
!
!> \brief   Perform aggregation for a 2D ice defined variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine performs aggregation by region for a 2D
!>  (:,nCells) field that is defined only over the ice covered part of
!>  the cell.
!
!-----------------------------------------------------------------------

  subroutine ice_aggregate_2D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         meshPool, &
         tracersAggregatePool, &
         regionsPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         outputFieldArray, &
         inputFieldArray

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         areaCell

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:,:), allocatable :: &
         outputFieldArrayTmp

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iDim, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(size(outputFieldArray,1),nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = 0.0_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                do iDim = 1, size(outputFieldArray,1)

                   outputFieldArrayTmp(iDim,iRegion) = outputFieldArrayTmp(iDim,iRegion) + &
                        inputFieldArray(iDim,iCell) * &
                        iceAreaCell(iCell) * &
                        areaCell(iCell)

                enddo ! iDim

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    do iDim = 1, size(outputFieldArray,1)
       call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, outputFieldArrayTmp(iDim,:), outputFieldArray(iDim,:))
    enddo ! iDim

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

  end subroutine ice_aggregate_2D

!***********************************************************************
!
!  routine cell_aggregate_1D
!
!> \brief   Perform aggregation for a 1D cell defined variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine performs aggregation by region for a 1D
!>  (nCells) field that is defined over all of the cell.
!
!-----------------------------------------------------------------------

  subroutine cell_aggregate_1D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         meshPool, &
         regionsPool

    real(kind=RKIND), dimension(:), pointer :: &
         outputFieldArray, &
         inputFieldArray, &
         areaCell

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:), allocatable :: &
         outputFieldArrayTmp

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = 0.0_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get the output array
       call MPAS_pool_get_subpool(block % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
       call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                outputFieldArrayTmp(iRegion) = outputFieldArrayTmp(iRegion) + &
                     inputFieldArray(iCell) * &
                     areaCell(iCell)

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, outputFieldArrayTmp, outputFieldArray)

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

  end subroutine cell_aggregate_1D

!***********************************************************************
!
!  routine cell_aggregate_2D
!
!> \brief   Perform aggregation for a 2D cell defined variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine performs aggregation by region for a 2D
!>  (:,nCells) field that is defined over all of the cell.
!
!-----------------------------------------------------------------------

  subroutine cell_aggregate_2D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         meshPool, &
         regionsPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         outputFieldArray, &
         inputFieldArray

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:,:), allocatable :: &
         outputFieldArrayTmp

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iDim, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(size(outputFieldArray,1),nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = 0.0_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                do iDim = 1, size(outputFieldArray,1)

                   outputFieldArrayTmp(iDim,iRegion) = outputFieldArrayTmp(iDim,iRegion) + &
                        inputFieldArray(iDim,iCell) * &
                        areaCell(iCell)

                enddo ! iDim

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    do iDim = 1, size(outputFieldArray,1)
       call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, outputFieldArrayTmp(iDim,:), outputFieldArray(iDim,:))
    enddo ! iDim

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

  end subroutine cell_aggregate_2D

!***********************************************************************
!
!  routine ice_areal_average_1D
!
!> \brief   Perform areal averaging for a 1D ice defined variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine performs areal averaging by region for a 1D
!>  (nCells) field that is defined only over the ice covered part of
!>  the cell.
!
!-----------------------------------------------------------------------

  subroutine ice_areal_average_1D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         meshPool, &
         tracersAggregatePool, &
         regionsPool

    real(kind=RKIND), dimension(:), pointer :: &
         outputFieldArray, &
         inputFieldArray, &
         iceAreaCell, &
         areaCell

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:), allocatable :: &
         outputFieldArrayTmp

    real(kind=RKIND), dimension(:), allocatable :: &
         denominatorTmp, &
         denominator

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = 0.0_RKIND

    ! allocate temporary denominator
    allocate(denominatorTmp(nRegions))

    ! initialize temporary denominator
    denominatorTmp = 0.0_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get the output array
       call MPAS_pool_get_subpool(block % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
       call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                outputFieldArrayTmp(iRegion) = outputFieldArrayTmp(iRegion) + &
                     inputFieldArray(iCell) * &
                     iceAreaCell(iCell) * &
                     areaCell(iCell)

                denominatorTmp(iRegion) = denominatorTmp(iRegion) + &
                     iceAreaCell(iCell) * &
                     areaCell(iCell)

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, outputFieldArrayTmp, outputFieldArray)

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

    ! allocate denominator across processors
    allocate(denominator(nRegions))

    ! aggregate denominator across processors
    call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, denominatorTmp, denominator)

    ! deallocate temporary denominator array
    deallocate(denominatorTmp)

    ! renormalize average
    do iRegion = 1, nRegions
       if (denominator(iRegion) > 0.0_RKIND) then
          outputFieldArray(iRegion) = outputFieldArray(iRegion) / denominator(iRegion)
       endif
    enddo ! iRegion

  end subroutine ice_areal_average_1D

!***********************************************************************
!
!  routine ice_areal_average_2D
!
!> \brief   Perform areal averaging for a 2D ice defined variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine performs areal averaging by region for a 2D
!>  (:,nCells) field that is defined only over the ice covered part of
!>  the cell.
!
!-----------------------------------------------------------------------

  subroutine ice_areal_average_2D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         meshPool, &
         tracersAggregatePool, &
         regionsPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         outputFieldArray, &
         inputFieldArray

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         areaCell

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:,:), allocatable :: &
         outputFieldArrayTmp

    real(kind=RKIND), dimension(:), allocatable :: &
         denominatorTmp, &
         denominator

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iDim, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(size(outputFieldArray,1),nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = 0.0_RKIND

    ! allocate temporary denominator
    allocate(denominatorTmp(nRegions))

    ! initialize temporary denominator
    denominatorTmp = 0.0_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                do iDim = 1, size(outputFieldArray,1)

                   outputFieldArrayTmp(iDim,iRegion) = outputFieldArrayTmp(iDim,iRegion) + &
                        inputFieldArray(iDim,iCell) * &
                        iceAreaCell(iCell) * &
                        areaCell(iCell)

                enddo ! iDim

                denominatorTmp(iRegion) = denominatorTmp(iRegion) + &
                     iceAreaCell(iCell) * &
                     areaCell(iCell)

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    do iDim = 1, size(outputFieldArray,1)
       call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, outputFieldArrayTmp(iDim,:), outputFieldArray(iDim,:))
    enddo ! iDim

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

    ! allocate denominator across processors
    allocate(denominator(nRegions))

    ! aggregate denominator across processors
    call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, denominatorTmp, denominator)

    ! deallocate temporary denominator array
    deallocate(denominatorTmp)

    ! renormalize average
    do iRegion = 1, nRegions
       if (denominator(iRegion) > 0.0_RKIND) then
          outputFieldArray(:,iRegion) = outputFieldArray(:,iRegion) / denominator(iRegion)
       endif
    enddo ! iRegion

  end subroutine ice_areal_average_2D

!***********************************************************************
!
!  routine cell_areal_average_1D
!
!> \brief   Perform areal averaging for a 1D cell defined variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine performs areal averaging by region for a 1D
!>  (nCells) field that is defined over all of the cell.
!
!-----------------------------------------------------------------------

  subroutine cell_areal_average_1D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         meshPool, &
         regionsPool

    real(kind=RKIND), dimension(:), pointer :: &
         outputFieldArray, &
         inputFieldArray, &
         areaCell

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:), allocatable :: &
         outputFieldArrayTmp

    real(kind=RKIND), dimension(:), allocatable :: &
         denominatorTmp, &
         denominator

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = 0.0_RKIND

    ! allocate temporary denominator
    allocate(denominatorTmp(nRegions))

    ! initialize temporary denominator
    denominatorTmp = 0.0_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get the output array
       call MPAS_pool_get_subpool(block % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
       call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                outputFieldArrayTmp(iRegion) = outputFieldArrayTmp(iRegion) + &
                     inputFieldArray(iCell) * &
                     areaCell(iCell)

                denominatorTmp(iRegion) = denominatorTmp(iRegion) + &
                     areaCell(iCell)

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, outputFieldArrayTmp, outputFieldArray)

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

    ! allocate denominator across processors
    allocate(denominator(nRegions))

    ! aggregate denominator across processors
    call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, denominatorTmp, denominator)

    ! deallocate temporary denominator array
    deallocate(denominatorTmp)

    ! renormalize average
    do iRegion = 1, nRegions
       if (denominator(iRegion) > 0.0_RKIND) then
          outputFieldArray(iRegion) = outputFieldArray(iRegion) / denominator(iRegion)
       endif
    enddo ! iRegion

  end subroutine cell_areal_average_1D

!***********************************************************************
!
!  routine cell_areal_average_2D
!
!> \brief   Perform areal averaging for a 2D cell defined variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine performs areal averaging by region for a 2D
!>  (:,nCells) field that is defined over all of the cell.
!
!-----------------------------------------------------------------------

  subroutine cell_areal_average_2D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         meshPool, &
         regionsPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         outputFieldArray, &
         inputFieldArray

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:,:), allocatable :: &
         outputFieldArrayTmp

    real(kind=RKIND), dimension(:), allocatable :: &
         denominatorTmp, &
         denominator

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iDim, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(size(outputFieldArray,1),nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = 0.0_RKIND

    ! allocate temporary denominator
    allocate(denominatorTmp(nRegions))

    ! initialize temporary denominator
    denominatorTmp = 0.0_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                do iDim = 1, size(outputFieldArray,1)

                   outputFieldArrayTmp(iDim,iRegion) = outputFieldArrayTmp(iDim,iRegion) + &
                        inputFieldArray(iDim,iCell) * &
                        areaCell(iCell)

                enddo ! iDim

                denominatorTmp(iRegion) = denominatorTmp(iRegion) + &
                     areaCell(iCell)

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    do iDim = 1, size(outputFieldArray,1)
       call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, outputFieldArrayTmp(iDim,:), outputFieldArray(iDim,:))
    enddo ! iDim

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

    ! allocate denominator across processors
    allocate(denominator(nRegions))

    ! aggregate denominator across processors
    call MPAS_dmpar_sum_real_array(domain % dminfo, nRegions, denominatorTmp, denominator)

    ! deallocate temporary denominator array
    deallocate(denominatorTmp)

    ! renormalize average
    do iRegion = 1, nRegions
       if (denominator(iRegion) > 0.0_RKIND) then
          outputFieldArray(:,iRegion) = outputFieldArray(:,iRegion) / denominator(iRegion)
       endif
    enddo ! iRegion

  end subroutine cell_areal_average_2D

!***********************************************************************
!
!  routine min_1D
!
!> \brief   Determine minimum for a 1D variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine finds the minimum by region for a 1D (nCells)
!>  field.
!
!-----------------------------------------------------------------------

  subroutine min_1D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         regionsPool

    real(kind=RKIND), dimension(:), pointer :: &
         outputFieldArray, &
         inputFieldArray

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:), allocatable :: &
         outputFieldArrayTmp

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = 1.0e34_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get the output array
       call MPAS_pool_get_subpool(block % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
       call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                outputFieldArrayTmp(iRegion) = min(outputFieldArrayTmp(iRegion), inputFieldArray(iCell))

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    call MPAS_dmpar_min_real_array(domain % dminfo, nRegions, outputFieldArrayTmp, outputFieldArray)

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

  end subroutine min_1D

!***********************************************************************
!
!  routine min_2D
!
!> \brief   Determine minimum for a 2D variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine finds the minimum by region for a 2D (:,nCells)
!>  field.
!
!-----------------------------------------------------------------------

  subroutine min_2D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         regionsPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         outputFieldArray, &
         inputFieldArray

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:,:), allocatable :: &
         outputFieldArrayTmp

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iDim, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(size(outputFieldArray,1),nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = 1.0e34_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                do iDim = 1, size(outputFieldArray,1)

                   outputFieldArrayTmp(iDim,iRegion) = min(outputFieldArrayTmp(iDim,iRegion), inputFieldArray(iDim,iCell))

                enddo ! iDim

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    do iDim = 1, size(outputFieldArray,1)
       call MPAS_dmpar_min_real_array(domain % dminfo, nRegions, outputFieldArrayTmp(iDim,:), outputFieldArray(iDim,:))
    enddo ! iDim

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

  end subroutine min_2D

!***********************************************************************
!
!  routine max_1D
!
!> \brief   Determine maximum for a 1D variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine finds the maximum by region for a 1D (nCells)
!>  field.
!
!-----------------------------------------------------------------------

  subroutine max_1D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         regionsPool

    real(kind=RKIND), dimension(:), pointer :: &
         outputFieldArray, &
         inputFieldArray

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:), allocatable :: &
         outputFieldArrayTmp

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = -1.0e34_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get the output array
       call MPAS_pool_get_subpool(block % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
       call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                outputFieldArrayTmp(iRegion) = max(outputFieldArrayTmp(iRegion), inputFieldArray(iCell))

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    call MPAS_dmpar_max_real_array(domain % dminfo, nRegions, outputFieldArrayTmp, outputFieldArray)

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

  end subroutine max_1D

!***********************************************************************
!
!  routine max_2D
!
!> \brief   Determine maximum for a 2D variable
!> \author  Adrian K. Turner
!> \date    28th March 2016
!> \details This routine finds the maximum by region for a 2D (:,nCells)
!>  field.
!
!-----------------------------------------------------------------------

  subroutine max_2D(&
       domain, &
       inputFieldName, &
       outputFieldName)

    type(domain_type), intent(inout) :: &
         domain

    character(len=*), intent(in) :: &
         inputFieldName, &
         outputFieldName

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         derivedFieldPool, &
         regionsPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         outputFieldArray, &
         inputFieldArray

    integer, dimension(:,:), pointer :: &
         regionCellMasks

    real(kind=RKIND), dimension(:,:), allocatable :: &
         outputFieldArrayTmp

    integer, pointer :: &
         nRegions, &
         nCellsSolve

    integer :: &
         iCell, &
         iDim, &
         iRegion

    ! get the number of regions
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nRegions", nRegions)

    ! get the output array
    call MPAS_pool_get_subpool(domain % blocklist % structs, "regionalStatisticsAM_derivedFieldsPool", derivedFieldPool)
    call MPAS_pool_get_array(derivedFieldPool, trim(outputFieldName), outputFieldArray)

    ! allocate temporary array
    allocate(outputFieldArrayTmp(size(outputFieldArray,1),nRegions))

    ! initialize aggregate sum
    outputFieldArrayTmp = -1.0e34_RKIND

    block => domain % blocklist
    do while (associated(block))

       ! get the input array
       call MPAS_pool_get_array(block % allFields, trim(inputFieldName), inputFieldArray)

       ! get other needed pools
       call MPAS_pool_get_subpool(block % structs, "regions", regionsPool)

       ! get other needed arrays
       call MPAS_pool_get_array(regionsPool, "regionCellMasks", regionCellMasks)

       ! get the needed dimensions
       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! aggregate the input field array
       do iCell = 1, nCellsSolve
          do iRegion = 1, nRegions

             if (regionCellMasks(iRegion,iCell) == 1) then

                do iDim = 1, size(outputFieldArray,1)

                   outputFieldArrayTmp(iDim,iRegion) = max(outputFieldArrayTmp(iDim,iRegion), inputFieldArray(iDim,iCell))

                enddo ! iDim

             endif ! region mask

          enddo ! iRegion
       enddo ! iCell

       block => block % next
    enddo

    ! aggregate across processors
    do iDim = 1, size(outputFieldArray,1)
       call MPAS_dmpar_max_real_array(domain % dminfo, nRegions, outputFieldArrayTmp(iDim,:), outputFieldArray(iDim,:))
    enddo ! iDim

    ! deallocate temporary array
    deallocate(outputFieldArrayTmp)

  end subroutine max_2D

!***********************************************************************
!
!  routine seaice_restart_regional_statistics
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_regional_statistics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine seaice_restart_regional_statistics!}}}

!***********************************************************************
!
!  routine seaice_finalize_regional_statistics
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_regional_statistics(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

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

   end subroutine seaice_finalize_regional_statistics!}}}

!-----------------------------------------------------------------------

end module seaice_regional_statistics

! vim: foldmethod=marker
