! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_cice_advection_incremental_remap_tracers
!
!> \brief  tracer setup for incremental remapping transport scheme
!> \author Adrian Turner and William Lipscomb
!> \date   March 2015
!> \details
!>  This module contains routines for setting up the tracer data structures
!>  required by the CICE incremental remapping scheme.
!
!-----------------------------------------------------------------------
module cice_advection_incremental_remap_tracers

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar

   implicit none
   private
   save

   public :: cice_add_tracers_to_linked_list, &
             cice_set_tracer_array_pointers,  &
             cice_update_tracer_halo

   ! the type that makes up each element in the linked list
   type, public :: tracer_type

      character(len=strKIND) :: tracerName
      character(len=strKIND) :: parentName

      integer :: ndims     ! number of dimensions in the tracer array
                           ! (nCategories, nCells) for 2D
                           ! (nLayer1, nCategories, nCells) for 3D
                           ! (nLayer1, nLayer2, nCategories, nCells) for 3D

      ! tracer arrays
      real(kind=RKIND), dimension(:,:),     pointer :: array2D => null()
      real(kind=RKIND), dimension(:,:,:),   pointer :: array3D => null()

      ! mass*tracer products
      ! = mass for mass field
      ! = mass*tracer1 for tracers with one parent
      ! = mass*tracer1*tracer2 for tracers with two parents
      ! = mass*tracer1*tracer2*tracer3 for tracers with three parents
      real(kind=RKIND), dimension(:,:),     pointer :: massTracerProduct2D => null()
      real(kind=RKIND), dimension(:,:,:),   pointer :: massTracerProduct3D => null()

      !TODO - Replace with one field, 2 time levels?
      ! global sums (over all cells)
      real(kind=RKIND), dimension(:), pointer :: globalSumInit2D => null()
      real(kind=RKIND), dimension(:), pointer :: globalSumFinal2D => null()
      real(kind=RKIND), dimension(:,:), pointer :: globalSumInit3D => null()
      real(kind=RKIND), dimension(:,:), pointer :: globalSumFinal3D => null()

      ! max/min tracer values in a local neighborhood
      real(kind=RKIND), dimension(:,:), pointer :: localMin2D => null()
      real(kind=RKIND), dimension(:,:), pointer :: localMax2D => null()
      real(kind=RKIND), dimension(:,:,:), pointer :: localMin3D => null()
      real(kind=RKIND), dimension(:,:,:), pointer :: localMax3D => null()

      ! mask: = 1 if hasChild = T and tracer value is physically meaningful, else = 0
      !       set at runtime in IR subroutine
      integer, dimension(:,:),     pointer :: arrayMask2D => null()
      integer, dimension(:,:,:),   pointer :: arrayMask3D => null()

      ! fluxes across cell edges
      !                   (nCategories, nEdges) for 2D
      !          (nLayers, nCategories, nEdges) for 3D
      real(kind=RKIND), dimension(:,:),     pointer :: edgeFlux2D => null()
      real(kind=RKIND), dimension(:,:,:),   pointer :: edgeFlux3D => null()

      ! tracer values at each quadrature point of departure triangles
      !                  (nCategories, nQuadPoints, nTriPerEdge, nEdges) for 2D
      !          (nLayers,nCategories, nQuadPoints, nTriPerEdge, nEdges) for 3D
      real(kind=RKIND), dimension(:,:,:,:),     pointer :: triangleValue2D => null()
      real(kind=RKIND), dimension(:,:,:,:,:),   pointer :: triangleValue3D => null()

      ! coordinates of barycenter associated with this tracer
      ! The term 'barycenter' refers to a center of mass or related quantity, as distinct from the geometric center.
      ! A mass-type field (0 parents) is defined to have its barycenter at the center of mass.
      ! A tracer with 1 parent has its barycenter at the center of mass*thisTracer.
      ! A tracer with 2 parents has its barycenter at the center of mass*parentTracer*thisTracer.
      ! And so on.
      ! The barycenter for this tracer is the location where the child tracer value (array2D or array3D) is located.
      ! Only mass-type fields (nParents = 0) are located at the geometric cell center.

      real(kind=RKIND), dimension(:,:),     pointer :: xBarycenter2D => null()
      real(kind=RKIND), dimension(:,:,:),   pointer :: xBarycenter3D => null()

      real(kind=RKIND), dimension(:,:),     pointer :: yBarycenter2D => null()
      real(kind=RKIND), dimension(:,:,:),   pointer :: yBarycenter3D => null()

      ! quantities needed for linear reconstruction (value at cell center plus x and y gradients)
      ! Note: The center value is the value at the geometric cell center and generally is
      !  difference from the value at the barycenter.
      ! xGrad and yGrad are the gradient components defined at the cell center

      real(kind=RKIND), dimension(:,:),     pointer :: center2D => null()
      real(kind=RKIND), dimension(:,:,:),   pointer :: center3D => null()

      real(kind=RKIND), dimension(:,:),     pointer :: xGrad2D => null()
      real(kind=RKIND), dimension(:,:,:),   pointer :: xGrad3D => null()

      real(kind=RKIND), dimension(:,:),     pointer :: yGrad2D => null()
      real(kind=RKIND), dimension(:,:,:),   pointer :: yGrad3D => null()

      ! pointer to parent tracer
      type(tracer_type), pointer :: parent => null()

      logical :: hasChild     ! true if this tracer has a child tracer

      integer :: nParents     ! = 1, 2 or 3 depending on place in tracer hierarchy
                              ! = 0 for mass-like field

      logical :: isActive     ! true if tracer is active for this run

      ! initial pointer to the next element in the linked list (not yet ordered)
      type(tracer_type), pointer :: nextInitial => null()

      ! pointer to the next element in the ordered linked list
      type(tracer_type), pointer :: next => null()

      ! some diagnostic arrays defined locally (for a single cell or a triangle within the cell)
      real(kind=RKIND), dimension(:),   pointer :: massTracerProductTriangle2D => null()
      real(kind=RKIND), dimension(:,:), pointer :: massTracerProductTriangle3D => null()

      real(kind=RKIND), dimension(:),   pointer :: massTracerProductCell2D => null()
      real(kind=RKIND), dimension(:,:), pointer :: massTracerProductCell3D => null()

      real(kind=RKIND), dimension(:,:),   pointer :: valueQP2D => null()
      real(kind=RKIND), dimension(:,:,:), pointer :: valueQP3D => null()

   end type tracer_type

   type(tracer_type), pointer, public :: tracersHead ! the linked list, but this always points to the first element

   logical, parameter :: verboseTracers = .false.
!   logical, parameter :: verboseTracers = .true.

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine cice_add_tracers_to_linked_list
!
!> \brief  construct a linked list of tracers
!> \author Adrian Turner and William Lipscomb
!> \date   March 2015
!> \details
!>  This routine constructs a linked list of tracers, with iceAreaCategory
!>  at the head of the list.
!-----------------------------------------------------------------------

  subroutine cice_add_tracers_to_linked_list(domain)

    type(domain_type), intent(in) :: &
         domain   !< Input:

    logical, pointer :: &
         pkgColumnTracerIceAgeActive, &
         pkgColumnTracerFirstYearIceActive, &
         pkgColumnTracerLevelIceActive, &
         pkgColumnTracerPondsActive, &
         pkgColumnTracerLidThicknessActive, &
         pkgColumnTracerAerosolsActive, &
         pkgColumnBiogeochemistryActive, &
         pkgTracerBrineActive, &
         pkgTracerMobileFractionActive, &
         pkgTracerSkeletalAlgaeActive, &
         pkgTracerSkeletalNitrateActive, &
         pkgTracerSkeletalCarbonActive, &
         pkgTracerSkeletalAmmoniumActive, &
         pkgTracerSkeletalSilicateActive, &
         pkgTracerSkeletalDMSActive, &
         pkgTracerSkeletalNonreactiveActive, &
         pkgTracerSkeletalHumicsActive, &
         pkgTracerSkeletalDONActive, &
         pkgTracerSkeletalIronActive, &
         pkgTracerVerticalAlgaeActive, &
         pkgTracerVerticalNitrateActive, &
         pkgTracerVerticalCarbonActive, &
         pkgTracerVerticalAmmoniumActive, &
         pkgTracerVerticalSilicateActive, &
         pkgTracerVerticalDMSActive, &
         pkgTracerVerticalNonreactiveActive, &
         pkgTracerVerticalHumicsActive, &
         pkgTracerVerticalDONActive, &
         pkgTracerVerticalIronActive, &
         pkgTracerZAerosolsActive, &
         pkgTracerZSalinityActive

    logical, pointer :: &
         config_use_level_meltponds

    type (mpas_pool_type), pointer :: configPool

    type(tracer_type), pointer :: thisTracer  ! diagnostic only

    nullify(tracersHead)
    nullify(thisTracer)

    ! add all the tracers to the linked list
    ! Note: These can be added in any order. They are reordered later to the order required for IR.
    ! Note: The fractional ice area (iceAreaCategory) is actually a mass-like field rather than a tracer.
    !       It is included here because it is the parent (or grandparent or great-grandparent) of all tracers
    !        and is part of the tracer pool.
    ! Note: The ice and snow volume are not tracers (since volume = area*thickness, where thickness is the tracer),
    !       but they are converted to thickness before being transported.

    ! First the tracers that are always included:

    call add_tracer_to_tracer_linked_list(tracersHead, 'iceAreaCategory',      'None')
    call add_tracer_to_tracer_linked_list(tracersHead, 'iceVolumeCategory',    'iceAreaCategory')
    call add_tracer_to_tracer_linked_list(tracersHead, 'snowVolumeCategory',   'iceAreaCategory')
    call add_tracer_to_tracer_linked_list(tracersHead, 'surfaceTemperature',   'iceAreaCategory')
    call add_tracer_to_tracer_linked_list(tracersHead, 'iceEnthalpy',          'iceVolumeCategory')
    call add_tracer_to_tracer_linked_list(tracersHead, 'iceSalinity',          'iceVolumeCategory')
    call add_tracer_to_tracer_linked_list(tracersHead, 'snowEnthalpy',         'snowVolumeCategory')

    ! Next the tracers that might be included, depending on column package options:

    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgColumnTracerIceAgeActive', &
                                                               pkgColumnTracerIceAgeActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgColumnTracerFirstYearIceActive', &
                                                               pkgColumnTracerFirstYearIceActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgColumnTracerLevelIceActive', &
                                                               pkgColumnTracerLevelIceActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgColumnTracerPondsActive', &
                                                               pkgColumnTracerPondsActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgColumnTracerLidThicknessActive', &
                                                               pkgColumnTracerLidThicknessActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgColumnTracerAerosolsActive', &
                                                               pkgColumnTracerAerosolsActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgColumnBiogeochemistryActive', &
                                                               pkgColumnBiogeochemistryActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerBrineActive', &
                                                               pkgTracerBrineActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerMobileFractionActive', &
                                                               pkgTracerMobileFractionActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerSkeletalAlgaeActive', &
                                                               pkgTracerSkeletalAlgaeActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerSkeletalNitrateActive', &
                                                               pkgTracerSkeletalNitrateActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerSkeletalCarbonActive', &
                                                               pkgTracerSkeletalCarbonActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerSkeletalAmmoniumActive', &
                                                               pkgTracerSkeletalAmmoniumActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerSkeletalSilicateActive', &
                                                               pkgTracerSkeletalSilicateActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerSkeletalDMSActive', &
                                                               pkgTracerSkeletalDMSActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerSkeletalNonreactiveActive', &
                                                               pkgTracerSkeletalNonreactiveActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerSkeletalHumicsActive', &
                                                               pkgTracerSkeletalHumicsActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerSkeletalDONActive', &
                                                               pkgTracerSkeletalDONActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerSkeletalIronActive', &
                                                               pkgTracerSkeletalIronActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerVerticalAlgaeActive', &
                                                               pkgTracerVerticalAlgaeActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerVerticalNitrateActive', &
                                                               pkgTracerVerticalNitrateActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerVerticalCarbonActive', &
                                                               pkgTracerVerticalCarbonActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerVerticalAmmoniumActive', &
                                                               pkgTracerVerticalAmmoniumActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerVerticalSilicateActive', &
                                                               pkgTracerVerticalSilicateActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerVerticalDMSActive', &
                                                               pkgTracerVerticalDMSActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerVerticalNonreactiveActive', &
                                                               pkgTracerVerticalNonreactiveActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerVerticalHumicsActive', &
                                                               pkgTracerVerticalHumicsActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerVerticalDONActive', &
                                                               pkgTracerVerticalDONActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerVerticalIronActive', &
                                                               pkgTracerVerticalIronActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerZAerosolsActive', &
                                                               pkgTracerZAerosolsActive)
    call MPAS_pool_get_package(domain % blocklist % packages, 'pkgTracerZSalinityActive', &
                                                               pkgTracerZSalinityActive)

    configPool => domain % blocklist % configs
    call MPAS_pool_get_config(configPool, "config_use_level_meltponds", config_use_level_meltponds)

    if (pkgColumnTracerIceAgeActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'iceAge',               'iceVolumeCategory')
    endif

    if (pkgColumnTracerFirstYearIceActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'firstYearIceArea',     'iceAreaCategory')
    endif

    if (pkgColumnTracerLevelIceActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'levelIceArea',         'iceAreaCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'levelIceVolume',       'iceVolumeCategory')
    endif

    if (pkgColumnTracerPondsActive) then
       if (config_use_level_meltponds) then
          call add_tracer_to_tracer_linked_list(tracersHead, 'pondArea',             'levelIceArea')
       else
          call add_tracer_to_tracer_linked_list(tracersHead, 'pondArea',             'iceAreaCategory')
       endif
       call add_tracer_to_tracer_linked_list(tracersHead, 'pondDepth',            'pondArea')
    endif

    if (pkgColumnTracerLidThicknessActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'pondLidThickness',     'pondArea')
    endif

    if (pkgColumnTracerAerosolsActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'snowScatteringAerosol','snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'snowBodyAerosol',      'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'iceScatteringAerosol', 'iceVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'iceBodyAerosol',       'iceVolumeCategory')
    endif

    ! brine height tracer
    if (pkgTracerBrineActive) &
       call add_tracer_to_tracer_linked_list(tracersHead, 'brineFraction',               'iceVolumeCategory')

    ! skeletal layer biogeochemistry
    if (pkgTracerSkeletalAlgaeActive) &
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalAlgaeConc',           'iceAreaCategory')

    ! nitrate
    if (pkgTracerSkeletalNitrateActive) &
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalNitrateConc',         'iceAreaCategory')

    ! carbon (DOC and DIC)
    if (pkgTracerSkeletalCarbonActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalDOCConc',             'iceAreaCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalDICConc',             'iceAreaCategory')
    endif

    ! DON (proteins...)
    if (pkgTracerSkeletalDONActive) &
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalDONConc',             'iceAreaCategory')

    ! ammonium
    if (pkgTracerSkeletalAmmoniumActive) &
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalAmmoniumConc',        'iceAreaCategory')

    ! silicate
    if (pkgTracerSkeletalSilicateActive) &
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalSilicateConc',        'iceAreaCategory')

    ! DMS and products
    if (pkgTracerSkeletalDMSActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalDMSConc',             'iceAreaCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalDMSPpConc',           'iceAreaCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalDMSPdConc',           'iceAreaCategory')
    endif

    ! nonreactive mobile tracer
    if (pkgTracerSkeletalNonreactiveActive) &
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalNonreactiveConc',     'iceAreaCategory')

    ! humic material
    if (pkgTracerSkeletalHumicsActive) &
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalHumicsConc',          'iceAreaCategory')

    ! iron
    if (pkgTracerSkeletalIronActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalParticulateIronConc', 'iceAreaCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'skeletalDissolvedIronConc',   'iceAreaCategory')
    endif

    ! vertical biogeochemistry
    ! fraction of tracer in the mobile phase
    if (pkgTracerMobileFractionActive) &
       call add_tracer_to_tracer_linked_list(tracersHead, 'mobileFraction',              'brineFraction')

    ! algal nitrogen
    if (pkgTracerVerticalAlgaeActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalAlgaeSnow',           'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalAlgaeIce',            'brineFraction')
    endif

    ! nitrate
    if (pkgTracerVerticalNitrateActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalNitrateSnow',         'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalNitrateIce',          'brineFraction')
    endif

    ! carbon
    if (pkgTracerVerticalCarbonActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDOCSnow',             'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDICSnow',             'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDOCIce',              'brineFraction')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDICIce',              'brineFraction')
    endif

    ! DON (proteins...)
    if (pkgTracerVerticalDONActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDONSnow',             'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDONIce',              'brineFraction')
    endif

    ! ammonium
    if (pkgTracerVerticalAmmoniumActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalAmmoniumSnow',        'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalAmmoniumIce',         'brineFraction')
    endif

    ! silicate
    if (pkgTracerVerticalSilicateActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalSilicateSnow',        'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalSilicateIce',         'brineFraction')
    endif

    ! DMS and related products
    if (pkgTracerVerticalDMSActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDMSSnow',             'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDMSPpSnow',           'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDMSPdSnow',           'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDMSIce',              'brineFraction')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDMSPpIce',            'brineFraction')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDMSPdIce',            'brineFraction')
    endif

    ! nonreactive mobile tracer
    if (pkgTracerVerticalNonreactiveActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalNonreactiveSnow',     'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalNonreactiveIce',      'brineFraction')
    endif

    ! humic material
    if (pkgTracerVerticalHumicsActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalHumicsSnow',          'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalHumicsIce',           'brineFraction')
    endif

    ! iron
    if (pkgTracerVerticalIronActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalParticulateIronSnow', 'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDissolvedIronSnow',   'snowVolumeCategory')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalParticulateIronIce',  'brineFraction')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalDissolvedIronIce',    'brineFraction')
    endif

    ! vertical aerosols
    if (pkgTracerZAerosolsActive) then
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalAerosolsIce','brineFraction')
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalAerosolsSnow','snowVolumeCategory')
    endif

    ! vertical salinity
    if (pkgTracerZSalinityActive) &
       call add_tracer_to_tracer_linked_list(tracersHead, 'verticalSalinity',            'brineFraction')

    if (verboseTracers) then
       thisTracer => tracersHead
       write(stderrUnit,*) ' '
       write(stderrUnit,*) 'Tracers in initial list:'
       write(stderrUnit,*) ' '
       do while (associated(thisTracer))
          write(stderrUnit,*) 'Tracer: ', trim(thisTracer % tracerName)
          if (associated(thisTracer % nextInitial)) then
             write(stderrUnit,*) 'next: ', trim(thisTracer % nextInitial % tracerName)
          endif
          write(stderrUnit,*) ' '
          thisTracer => thisTracer % nextInitial
       enddo
    endif

    ! set the connectivity in the list
    call create_tracers_connectivity(tracersHead)

    ! calculate the number of parents for each tracer
    call set_parent_numbers(tracersHead)

    ! remove inactive tracers from the list
    call remove_inactive_tracers(domain, tracersHead)

    ! set the ordered list
    ! The resulting list will have the Type 1 tracers first, followed by Type 2 and then Type 3
    ! IR requires this order so as to compute parent properties before they are needed by child tracers
    call set_tracer_order(tracersHead)

    !WHL - debug
    if (verboseTracers) then
       thisTracer => tracersHead
       write(stderrUnit,*) ' '
       write(stderrUnit,*) 'Tracers in ordered list:'
       write(stderrUnit,*) ' '
       do while (associated(thisTracer))
          write(stderrUnit,*) 'Tracer: ', trim(thisTracer % tracerName)
          write(stderrUnit,*) 'Parent: ', trim(thisTracer % parentName)
          write(stderrUnit,*) 'hasChild: ', thisTracer % hasChild
          write(stderrUnit,*) 'nParents: ', thisTracer % nParents
          if (associated(thisTracer % next)) then
             write(stderrUnit,*) 'next: ', trim(thisTracer % next % tracerName)
          endif
          write(stderrUnit,*) ' '
          thisTracer => thisTracer % next
       enddo
       if (associated(tracersHead)) write(stderrUnit,*) 'tracersHead:', trim(tracersHead % tracerName)
    endif

  end subroutine cice_add_tracers_to_linked_list

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine add_tracer_to_tracer_linked_list
!
!> \brief  add a tracer to the linked list
!> \author Adrian Turner and William Lipscomb
!> \date   March 2015
!> \details
!>  This routine adds a tracer to the linked list.
!-----------------------------------------------------------------------

  subroutine add_tracer_to_tracer_linked_list(tracersHead, &
                                              tracerName,  &
                                              parentName)

    type(tracer_type), pointer :: &
         tracersHead ! pointer to first element of linked list

    character(len=*), intent(in) :: &
         tracerName, &
         parentName

    type(tracer_type), pointer :: &
         newTracer

    ! see if the linked list has any elements yet
    if (.not. associated(tracersHead)) then

       ! no elements yet
       allocate(tracersHead)
       nullify(tracersHead % nextInitial)
       newTracer => tracersHead

    else

       ! already have elements in linked list

       ! do this so we do not alter the original pointer; it is always the head of the linked list
       newTracer => tracersHead

       ! search for the end of the list
       do while (associated(newTracer % nextInitial))
          newTracer => newTracer % nextInitial
       enddo

       ! found the end so allocate
       allocate(newTracer % nextInitial)
       newTracer => newTracer % nextInitial
       nullify(newTracer % nextInitial)

    endif

    ! newTracer is the new linked list element ready to have stuff done to it!
    newTracer % tracerName = trim(tracerName)
    newTracer % parentName = trim(parentName)

    ! initialize hasChild and nParents; these are modified below
    newTracer % hasChild = .false.
    newTracer % nParents = 0

  end subroutine add_tracer_to_tracer_linked_list

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine create_tracers_connectivity
!
!> \brief  set pointer from child to parent tracers
!> \author Adrian Turner and William Lipscomb
!> \date   March 2015
!> \details
!>  This routine sets pointers from each child tracer to its parent tracer.
!-----------------------------------------------------------------------

  subroutine create_tracers_connectivity(tracersHead)

    type(tracer_type), pointer, intent(in) :: &
         tracersHead ! do not alter this pointer - it should always point to the first element

    type(tracer_type), pointer :: &
         thisTracer,  &
         parentTracer

    ! loop over the tracers in the linked list
    thisTracer => tracersHead
    do while (associated(thisTracer))

       ! loop over potential parent tracers
       parentTracer => tracersHead
       do while (associated(parentTracer))

          ! see if the potential parent tracer is the parent of the child tracer
          if (trim(parentTracer % tracerName) == trim(thisTracer % parentName)) then

             ! set the child parent pointer to the parent tracer
             thisTracer % parent => parentTracer

             ! set the parent tracer to have children
             parentTracer % hasChild = .true.

          endif

          parentTracer => parentTracer % nextInitial

       enddo

       thisTracer => thisTracer % nextInitial

    enddo  ! associated(thisTracer)

  end subroutine create_tracers_connectivity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine set_parent_numbers
!
!> \brief  compute the number of parents for each tracer
!> \author Adrian Turner and William Lipscomb
!> \date   March 2015
!> \details
!>  This routine determines the number of parents of each tracer.
!-----------------------------------------------------------------------

  subroutine set_parent_numbers(tracersHead)

    type(tracer_type), pointer, intent(in) :: &
         tracersHead ! do not alter this pointer - it should always point to the first element

    type(tracer_type), pointer :: &
         thisTracer, &
         thisTracerParent

    ! loop over the tracers in the linked list
    thisTracer => tracersHead
    do while (associated(thisTracer))

       thisTracer % nParents = 0

       thisTracerParent => thisTracer % parent

       do while (associated(thisTracerParent))

          thisTracer % nParents = thisTracer % nParents + 1

          thisTracerParent => thisTracerParent % parent

       enddo

       thisTracer => thisTracer % nextInitial

    enddo

  end subroutine set_parent_numbers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine set_tracers_active
!
!> \brief  determine which tracers are active
!> \author Adrian Turner and William Lipscomb
!> \date   March 2015
!> \details
!>  This routine determines whether each tracer is active.
!-----------------------------------------------------------------------

  subroutine set_tracers_active(domain, tracersHead)

    type(domain_type), intent(in) :: domain

    type(tracer_type), pointer, intent(in) :: &
         tracersHead ! do not alter this pointer - it should always point to the first element

    type(tracer_type), pointer :: &
         tracer

    type (mpas_pool_type), pointer :: tracerPool

    type(MPAS_pool_field_info_type) :: fieldInfo

    ! loop over the tracers in the linked list
    tracer => tracersHead
    do while (associated(tracer))

       call MPAS_pool_get_subpool(domain % blocklist % structs, 'tracers', tracerPool)

       call MPAS_pool_get_field_info(tracerPool, trim(tracer % tracerName), fieldInfo)

       tracer % isActive = fieldInfo % isActive

       tracer => tracer % nextInitial

    enddo

  end subroutine set_tracers_active

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine remove_inactive_tracers
!
!> \brief  remove inactive tracers from the linked list
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  This routine removes inactive tracers from the linked list.
!-----------------------------------------------------------------------

  subroutine remove_inactive_tracers(domain, tracersHead)

    type(domain_type), intent(in) :: domain

    type(tracer_type), pointer, intent(inout) :: &
         tracersHead ! this pointer should always point to the first element
                     ! if the first element is inactive, it will be moved to the first active element

    type(tracer_type), pointer :: &
         tracer, &
         tracerPrev,  &
         childTracer

    type (mpas_pool_type), pointer :: tracerPool

    type(MPAS_pool_field_info_type) :: fieldInfo

    logical :: first ! initially true; set to false after the first active tracer is found

    call MPAS_pool_get_subpool(domain % blocklist % structs, 'tracers', tracerPool)

    ! loop over the tracers in the linked list, and identify the active tracers
    tracer => tracersHead
    do while (associated(tracer))

       call MPAS_pool_get_field_info(tracerPool, trim(tracer % tracerName), fieldInfo)

       tracer % isActive = fieldInfo % isActive

       tracer => tracer % nextInitial

    enddo

    ! loop over the tracers again, removing inactive tracers from the list

    first = .true.
    tracer => tracersHead
    tracerPrev => tracersHead

    do while (associated(tracer))

       if (tracer % isActive) then

          if (verboseTracers) then
             write(stderrUnit,*) 'Active tracer: ', trim(tracer % tracerName)
          endif

          ! move pointer to the next tracer
          tracerPrev => tracer
          tracer => tracer % nextInitial
          first = .false.

       else

          ! point the previous tracer to the next tracer (if not null) and
          ! remove this tracer from the list

          if (verboseTracers) then
             write(stderrUnit,*) 'Inactive tracer:', trim(tracer % tracerName)
          endif

          if (associated(tracer % nextInitial)) then
             tracerPrev % nextInitial => tracer % nextInitial

             if (first) then   ! this is the first element; move tracersHead before nullifying the tracer
                tracersHead => tracer % nextInitial
             endif

             ! make sure the tracer to be nullified has no active children
             if (tracer % hasChild) then

                ! loop over potential child tracers
                childTracer => tracersHead
                do while (associated(childTracer))

                   if (childTracer % nParents > 1) then  ! potential child tracer
                      if (trim(childTracer % parent % tracerName) == trim(tracer % tracerName)) then

                         ! this is a child of the tracer to be nullified; make sure it is not active
                         if (childTracer % isActive) then
                            write(stderrUnit,*) 'ERROR, IR tracers: Cannot have an inactive parent with an active child'
                            write(stderrUnit,*) 'Inactive parent:', trim(tracer % tracerName)
                            write(stderrUnit,*) 'Active child:', trim(childTracer % tracerName)
                            call mpas_dmpar_global_abort(&
                                 'MPAS-seaice: IR tracers: Cannot have an inactive parent with an active child')
                         endif
                      endif
                   endif     ! nParents > 1

                   childTracer => childTracer % nextInitial

                enddo

             endif  ! tracer has child

             nullify(tracer)

             tracer => tracerPrev % nextInitial

          else  ! tracer % nextInitial is not associated; this is the end of the list

             nullify(tracerPrev % nextInitial)
             nullify(tracer)

             if (first) then ! the list is empty; something is wrong
                nullify(tracersHead)
                write(stderrUnit,*) 'ERROR: No active tracers for IR'
                call mpas_dmpar_global_abort('MPAS-seaice: IR tracers: No active tracers')
             endif

          endif  ! associated(tracer % nextInitial)

       endif     ! tracer is active

    enddo        ! tracer loop

  end subroutine remove_inactive_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine set_tracer_order
!
!> \brief  order the tracer linked list
!> \author Adrian Turner and William Lipscomb
!> \date   March 2015
!> \details
!>  This routine orders the tracer linked list such that parent tracers
!>  always appear before their children.
!-----------------------------------------------------------------------

  subroutine set_tracer_order(tracersHead)

    type(tracer_type), pointer :: &
         tracersHead ! pointer to the first element of the linked list
                     ! altered here to point to the first element of the newly ordered list

    type(tracer_type), pointer :: &
         tracer, &
         tracerPrev, &
         tracerFirst

    integer :: &
         maxLevels, &
         tracerLevel

    ! determine the number of levels (1, 2 or 3, from max value of nParents)

    maxLevels = 0

    tracer => tracersHead
    do while (associated(tracer))

       maxLevels = max(maxLevels,tracer % nParents)

       tracer => tracer % nextInitial

    enddo

    ! loop over the levels making sure the order is by levels
    nullify(tracerPrev)

    do tracerLevel = 0, maxLevels

       ! loop over tracers
       tracer => tracersHead
       do while (associated(tracer))

          ! see if the tracer is at the current level
          if (tracer % nParents == tracerLevel) then

             ! see if the previous tracer pointer is associated
             if (associated(tracerPrev)) then
                ! the previous tracer is associated so set the next pointer
                tracerPrev % next => tracer   ! make this the next tracer in the linked list
                tracerPrev => tracer          ! make tracerPrev point to this tracer for next time
             else
                ! set the first pointer in the new order
                tracerFirst => tracer
                tracerPrev => tracer
             endif

          endif

          tracer => tracer % nextInitial

       enddo

    enddo

    ! make the tracer object point to the new start of the list
    tracersHead => tracerFirst

  end subroutine set_tracer_order

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine cice_set_tracer_array_pointers
!
!> \brief  set tracer pointers to the appropriate field objects
!> \author Adrian Turner and William Lipscomb
!> \date   March 2015
!> \details
!>  This routine sets tracer pointers to field objects retrieved from
!>  the appropriate pools.
!-----------------------------------------------------------------------

  subroutine cice_set_tracer_array_pointers(tracersHead, block, timeLevelIn)

    type(tracer_type), pointer, intent(in) :: &
         tracersHead ! do not alter this pointer - it should always point to the first element

    type(block_type), intent(inout) :: block

    integer, intent(in), optional :: timeLevelIn

    type(tracer_type), pointer :: &
         thisTracer

    type(MPAS_pool_type), pointer :: &
         tracerPool,              &! ice area plus tracer fields
         tracerMaskPool,          &! integer mask for each tracer
         tracerEdgeFluxPool,      &! edge flux for each tracer
         tracerProductPool,       &! mass-tracer products for each tracer
         tracerTrianglePool,      &! tracer values at quadrature points of departure triangles
         tracerBarycenterPool,    &! coordinates of barycenter associated with each tracer
         tracerReconstructionPool,&! quantities needed for linear reconstuction of each tracer
         tracerConservationPool,  &! conserved quantities for each tracer
         tracerMonotonicityPool    ! local max/min tracer values for checking monotonicity

    type(MPAS_pool_field_info_type) :: thisFieldInfo

    type(field1DReal), pointer :: tracerField1DReal
    type(field2DReal), pointer :: tracerField2DReal
    type(field3DReal), pointer :: tracerField3DReal
    type(field4DReal), pointer :: tracerField4DReal
    type(field5DReal), pointer :: tracerField5DReal

    type(field2DInteger), pointer :: tracerField2DInteger
    type(field3DInteger), pointer :: tracerField3DInteger

    integer :: timeLevel

    if (present(timeLevelIn)) then
       timeLevel = timeLevelIn
    else
       timeLevel = 1
    endif

    call MPAS_pool_get_subpool(block % structs, 'tracers', tracerPool)
    call MPAS_pool_get_subpool(block % structs, 'tracer_masks', tracerMaskPool)
    call MPAS_pool_get_subpool(block % structs, 'tracer_edge_fluxes', tracerEdgeFluxPool)
    call MPAS_pool_get_subpool(block % structs, 'tracer_products', tracerProductPool)
    call MPAS_pool_get_subpool(block % structs, 'tracer_triangles', tracerTrianglePool)
    call MPAS_pool_get_subpool(block % structs, 'tracer_barycenter', tracerBarycenterPool)
    call MPAS_pool_get_subpool(block % structs, 'tracer_reconstruction', tracerReconstructionPool)
    call MPAS_pool_get_subpool(block % structs, 'tracer_conservation', tracerConservationPool)
    call MPAS_pool_get_subpool(block % structs, 'tracer_monotonicity', tracerMonotonicityPool)

    ! loop over the tracers in the linked list

    thisTracer => tracersHead

    do while (associated(thisTracer))

       call MPAS_pool_get_field_info(tracerPool, trim(thisTracer % tracerName), thisFieldInfo)

       !NOTE: As of Aug. 2015, all tracers are formally 3D. But some (e.g., iceAreaCategory) are really 2D,
       !       with a first (layer) index of dimension 1.
       !      For purposes of IR, it is easier to define a 2D pointer to iceAreaCategory(1,:,:) than
       !       a 3D pointer to iceAreaCategory(:,:,:).  This avoids some artificial logic in the main remapping
       !       subroutines to handle parent and child tracers that are both 3D but with different first dimensions.
       !TODO - Remove this extra logic when fields like iceAreaCategory are formally 2D.

       if (thisFieldInfo%nDims == 2) then

          thisTracer % nDims = 2

       elseif (thisFieldInfo%nDims == 3) then   ! check whether the layer index has size > 1

          call MPAS_pool_get_field(tracerPool, trim(thisTracer % tracerName), tracerField3DReal, timeLevel)

          if (size(tracerField3DReal % array, 1) == 1) then   ! treat this as a 2D tracer
             thisTracer % nDims = 2
          else   ! size > 1; treat this as a 3D tracer
             thisTracer % nDims = 3
          endif

       endif

       if (thisFieldInfo%nDims == 2 .and. thisTracer % nDims == 2) then    ! 2D formally and physically

          call MPAS_pool_get_field(tracerPool, trim(thisTracer % tracerName), tracerField2DReal, timeLevel)
          thisTracer % array2D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerMaskPool, trim(thisTracer % tracerName) // 'Mask', tracerField2DInteger)
          thisTracer % arrayMask2D => tracerField2DInteger % array
          call MPAS_pool_get_field(tracerEdgeFluxPool, trim(thisTracer % tracerName) // 'EdgeFlux', tracerField2DReal)
          thisTracer % edgeFlux2D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerProductPool, trim(thisTracer % tracerName) // 'Product', tracerField2DReal)
          thisTracer % massTracerProduct2D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerTrianglePool, trim(thisTracer % tracerName) // 'Triangle', tracerField4DReal)
          thisTracer % triangleValue2D => tracerField4DReal % array
          call MPAS_pool_get_field(tracerBarycenterPool, trim(thisTracer % tracerName) // 'Barycenterx', tracerField2DReal)
          thisTracer % xBarycenter2D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerBarycenterPool, trim(thisTracer % tracerName) // 'Barycentery', tracerField2DReal)
          thisTracer % yBarycenter2D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerReconstructionPool, trim(thisTracer % tracerName) // 'Center', tracerField2DReal)
          thisTracer % center2D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerReconstructionPool, trim(thisTracer % tracerName) // 'Gradx', tracerField2DReal)
          thisTracer % xGrad2D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerReconstructionPool, trim(thisTracer % tracerName) // 'Grady', tracerField2DReal)
          thisTracer % yGrad2D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerConservationPool, trim(thisTracer % tracerName) // 'Cons', tracerField1DReal, 1)
          thisTracer % globalSumInit2D => tracerField1DReal % array
          call MPAS_pool_get_field(tracerConservationPool, trim(thisTracer % tracerName) // 'Cons', tracerField1DReal, 2)
          thisTracer % globalSumFinal2D => tracerField1DReal % array
          call MPAS_pool_get_field(tracerMonotonicityPool, trim(thisTracer % tracerName) // 'LocalMin', tracerField2DReal)
          thisTracer % localMin2D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerMonotonicityPool, trim(thisTracer % tracerName) // 'LocalMax', tracerField2DReal)
          thisTracer % localMax2D => tracerField2DReal % array

       elseif (thisFieldInfo%nDims == 3 .and. thisTracer % nDims == 2) then    ! 3D formally but 2D physically

          call MPAS_pool_get_field(tracerPool, trim(thisTracer % tracerName), tracerField3DReal, timeLevel)
          thisTracer % array2D => tracerField3DReal % array(1,:,:)
          call MPAS_pool_get_field(tracerMaskPool, trim(thisTracer % tracerName) // 'Mask', tracerField3DInteger)
          thisTracer % arrayMask2D => tracerField3DInteger % array(1,:,:)
          call MPAS_pool_get_field(tracerEdgeFluxPool, trim(thisTracer % tracerName) // 'EdgeFlux', tracerField3DReal)
          thisTracer % edgeFlux2D => tracerField3DReal % array(1,:,:)
          call MPAS_pool_get_field(tracerProductPool, trim(thisTracer % tracerName) // 'Product', tracerField3DReal)
          thisTracer % massTracerProduct2D => tracerField3DReal % array(1,:,:)
          call MPAS_pool_get_field(tracerTrianglePool, trim(thisTracer % tracerName) // 'Triangle', tracerField5DReal)
          thisTracer % triangleValue2D => tracerField5DReal % array(1,:,:,:,:)
          call MPAS_pool_get_field(tracerBarycenterPool, trim(thisTracer % tracerName) // 'Barycenterx', tracerField3DReal)
          thisTracer % xBarycenter2D => tracerField3DReal % array(1,:,:)
          call MPAS_pool_get_field(tracerBarycenterPool, trim(thisTracer % tracerName) // 'Barycentery', tracerField3DReal)
          thisTracer % yBarycenter2D => tracerField3DReal % array(1,:,:)
          call MPAS_pool_get_field(tracerReconstructionPool, trim(thisTracer % tracerName) // 'Center', tracerField3DReal)
          thisTracer % center2D => tracerField3DReal % array(1,:,:)
          call MPAS_pool_get_field(tracerReconstructionPool, trim(thisTracer % tracerName) // 'Gradx', tracerField3DReal)
          thisTracer % xGrad2D => tracerField3DReal % array(1,:,:)
          call MPAS_pool_get_field(tracerReconstructionPool, trim(thisTracer % tracerName) // 'Grady', tracerField3DReal)
          thisTracer % yGrad2D => tracerField3DReal % array(1,:,:)
          call MPAS_pool_get_field(tracerConservationPool, trim(thisTracer % tracerName) // 'Cons', tracerField2DReal, 1)
          thisTracer % globalSumInit2D => tracerField2DReal % array(1,:)
          call MPAS_pool_get_field(tracerConservationPool, trim(thisTracer % tracerName) // 'Cons', tracerField2DReal, 2)
          thisTracer % globalSumFinal2D => tracerField2DReal % array(1,:)
          call MPAS_pool_get_field(tracerMonotonicityPool, trim(thisTracer % tracerName) // 'LocalMin', tracerField3DReal)
          thisTracer % localMin2D => tracerField3DReal % array(1,:,:)
          call MPAS_pool_get_field(tracerMonotonicityPool, trim(thisTracer % tracerName) // 'LocalMax', tracerField3DReal)
          thisTracer % localMax2D => tracerField3DReal % array(1,:,:)

       elseif (thisFieldInfo%nDims == 3 .and. thisTracer % nDims == 3) then    ! 3D formally and physically

          call MPAS_pool_get_field(tracerPool, trim(thisTracer % tracerName), tracerField3DReal, timeLevel)
          thisTracer % array3D => tracerField3DReal % array
          call MPAS_pool_get_field(tracerMaskPool, trim(thisTracer % tracerName) // 'Mask', tracerField3DInteger)
          thisTracer % arrayMask3D => tracerField3DInteger % array
          call MPAS_pool_get_field(tracerEdgeFluxPool, trim(thisTracer % tracerName) // 'EdgeFlux', tracerField3DReal)
          thisTracer % edgeFlux3D => tracerField3DReal % array
          call MPAS_pool_get_field(tracerProductPool, trim(thisTracer % tracerName) // 'Product', tracerField3DReal)
          thisTracer % massTracerProduct3D => tracerField3DReal % array
          call MPAS_pool_get_field(tracerTrianglePool, trim(thisTracer % tracerName) // 'Triangle', tracerField5DReal)
          thisTracer % triangleValue3D => tracerField5DReal % array
          call MPAS_pool_get_field(tracerBarycenterPool, trim(thisTracer % tracerName) // 'Barycenterx', tracerField3DReal)
          thisTracer % xBarycenter3D => tracerField3DReal % array
          call MPAS_pool_get_field(tracerBarycenterPool, trim(thisTracer % tracerName) // 'Barycentery', tracerField3DReal)
          thisTracer % yBarycenter3D => tracerField3DReal % array
          call MPAS_pool_get_field(tracerReconstructionPool, trim(thisTracer % tracerName) // 'Center', tracerField3DReal)
          thisTracer % center3D => tracerField3DReal % array
          call MPAS_pool_get_field(tracerReconstructionPool, trim(thisTracer % tracerName) // 'Gradx', tracerField3DReal)
          thisTracer % xGrad3D => tracerField3DReal % array
          call MPAS_pool_get_field(tracerReconstructionPool, trim(thisTracer % tracerName) // 'Grady', tracerField3DReal)
          thisTracer % yGrad3D => tracerField3DReal % array
          call MPAS_pool_get_field(tracerConservationPool, trim(thisTracer % tracerName) // 'Cons', tracerField2DReal, 1)
          thisTracer % globalSumInit3D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerConservationPool, trim(thisTracer % tracerName) // 'Cons', tracerField2DReal, 2)
          thisTracer % globalSumFinal3D => tracerField2DReal % array
          call MPAS_pool_get_field(tracerMonotonicityPool, trim(thisTracer % tracerName) // 'LocalMin', tracerField3DReal)
          thisTracer % localMin3D => tracerField3DReal % array
          call MPAS_pool_get_field(tracerMonotonicityPool, trim(thisTracer % tracerName) // 'LocalMax', tracerField3DReal)
          thisTracer % localMax3D => tracerField3DReal % array

       endif  ! nDims

       thisTracer => thisTracer % next

    enddo   ! associated(tracer)

  end subroutine cice_set_tracer_array_pointers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine cice_update_tracer_halo
!
!> \brief  halo update for tracers
!> \author William Lipscomb
!> \date   March 2015
!> \details
!>  This routine does a halo update for each tracer in the linked list.
!-----------------------------------------------------------------------

  subroutine cice_update_tracer_halo(tracersHead, domain, timeLevelIn)

    use mpas_dmpar

    type(tracer_type), pointer, intent(in) :: &
         tracersHead ! do not alter this pointer - it should always point to the first element

    type(domain_type), intent(inout) :: domain

    integer, intent(in), optional :: timeLevelIn

    type(tracer_type), pointer :: &
         thisTracer

    integer :: timeLevel

    if (present(timeLevelIn)) then
       timeLevel = timeLevelIn
    else
       timeLevel = 1
    endif

    thisTracer => tracersHead   ! point to first element of linked list
    do while(associated(thisTracer))

       if (verboseTracers) write(stderrUnit,*) 'Halo update:', trim(thisTracer % tracerName)

       call MPAS_dmpar_field_halo_exch(domain, trim(thisTracer % tracerName), timeLevel=timeLevel)

       thisTracer => thisTracer % next

    enddo  ! associated(thisTracer)

  end subroutine cice_update_tracer_halo

  !-----------------------------------------------

end module cice_advection_incremental_remap_tracers
