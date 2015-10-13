!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_diagnostics
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date August 20th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_diagnostics

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping

  implicit none

  private
  save 

  public :: &
       cice_run_diagnostics

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_run_diagnostics
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date August 20th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_run_diagnostics(block, clock)

    type(block_type), intent(inout) :: &
         block !< Input/Output: 

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input: 

    ! min max values for various variables
    call cice_diagnostic_min_max(block, clock)

    ! global sums
    call cice_diagnostic_global(block, clock)

  end subroutine cice_run_diagnostics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_diagnostic_min_max
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date August 20th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_diagnostic_min_max(block, clock)

    use cice_constants, only: &
         iceAreaMinimum

    type(block_type), intent(inout) :: &
         block !< Input/Output: 

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input: 

    type(MPAS_pool_type), pointer :: &
         mesh, &
         global_diagnostics, &
         tracers

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCategoryMin, &
         iceAreaCategoryMax, &
         iceVolumeCategoryMin, &
         iceVolumeCategoryMax, &
         snowVolumeCategoryMin, &
         snowVolumeCategoryMax, &
         surfaceTemperatureMin, &
         surfaceTemperatureMax, &
         iceThicknessCategoryMin, &
         iceThicknessCategoryMax, &
         snowThicknessCategoryMin, &
         snowThicknessCategoryMax

    integer, dimension(:), pointer :: &
         snowThicknessCategoryLocMax

    real(kind=RKIND), dimension(:,:,:), pointer :: &
        iceAreaCategory, &
        iceVolumeCategory, &
        snowVolumeCategory, &
        surfaceTemperature

    integer, pointer :: &
        nCells, &
        nCategories

    integer :: &
         iCell, &
         iCategory

    real(kind=RKIND) :: &
         iceThickness, &
         snowThickness

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "global_diagnostics", global_diagnostics)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_array(global_diagnostics, "iceAreaCategoryMin", iceAreaCategoryMin)
    call MPAS_pool_get_array(global_diagnostics, "iceAreaCategoryMax", iceAreaCategoryMax)
    call MPAS_pool_get_array(global_diagnostics, "iceVolumeCategoryMin", iceVolumeCategoryMin)
    call MPAS_pool_get_array(global_diagnostics, "iceVolumeCategoryMax", iceVolumeCategoryMax)
    call MPAS_pool_get_array(global_diagnostics, "snowVolumeCategoryMin", snowVolumeCategoryMin)
    call MPAS_pool_get_array(global_diagnostics, "snowVolumeCategoryMax", snowVolumeCategoryMax)
    call MPAS_pool_get_array(global_diagnostics, "surfaceTemperatureMin", surfaceTemperatureMin)
    call MPAS_pool_get_array(global_diagnostics, "surfaceTemperatureMax", surfaceTemperatureMax)
    call MPAS_pool_get_array(global_diagnostics, "iceThicknessCategoryMin", iceThicknessCategoryMin)
    call MPAS_pool_get_array(global_diagnostics, "iceThicknessCategoryMax", iceThicknessCategoryMax)
    call MPAS_pool_get_array(global_diagnostics, "snowThicknessCategoryMin", snowThicknessCategoryMin)
    call MPAS_pool_get_array(global_diagnostics, "snowThicknessCategoryMax", snowThicknessCategoryMax)
    call MPAS_pool_get_array(global_diagnostics, "snowThicknessCategoryLocMax", snowThicknessCategoryLocMax)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    iceAreaCategoryMin       =  1.0e34_RKIND
    iceAreaCategoryMax       = -1.0e34_RKIND
    iceVolumeCategoryMin     =  1.0e34_RKIND
    iceVolumeCategoryMax     = -1.0e34_RKIND
    snowVolumeCategoryMin    =  1.0e34_RKIND
    snowVolumeCategoryMax    = -1.0e34_RKIND
    surfaceTemperatureMin    =  1.0e34_RKIND
    surfaceTemperatureMax    = -1.0e34_RKIND
    iceThicknessCategoryMin  =  1.0e34_RKIND
    iceThicknessCategoryMax  = -1.0e34_RKIND
    snowThicknessCategoryMin =  1.0e34_RKIND
    snowThicknessCategoryMax = -1.0e34_RKIND

    snowThicknessCategoryLocMax = -1

    do iCell = 1, nCells
       do iCategory = 1, nCategories

          if (iceAreaCategory(1,iCategory,iCell) > iceAreaMinimum) then

             iceAreaCategoryMin(iCategory)    = min(iceAreaCategoryMin(iCategory),    iceAreaCategory(1,iCategory,iCell))
             iceAreaCategoryMax(iCategory)    = max(iceAreaCategoryMax(iCategory),    iceAreaCategory(1,iCategory,iCell))

             iceVolumeCategoryMin(iCategory)  = min(iceVolumeCategoryMin(iCategory),  iceVolumeCategory(1,iCategory,iCell))
             iceVolumeCategoryMax(iCategory)  = max(iceVolumeCategoryMax(iCategory),  iceVolumeCategory(1,iCategory,iCell))

             snowVolumeCategoryMin(iCategory) = min(snowVolumeCategoryMin(iCategory), snowVolumeCategory(1,iCategory,iCell))
             snowVolumeCategoryMax(iCategory) = max(snowVolumeCategoryMax(iCategory), snowVolumeCategory(1,iCategory,iCell))

             surfaceTemperatureMin(iCategory) = min(surfaceTemperatureMin(iCategory), surfaceTemperature(1,iCategory,iCell))
             surfaceTemperatureMax(iCategory) = max(surfaceTemperatureMax(iCategory), surfaceTemperature(1,iCategory,iCell))

             iceThickness  = iceVolumeCategory(1,iCategory,iCell)  / iceAreaCategory(1,iCategory,iCell)
             snowThickness = snowVolumeCategory(1,iCategory,iCell) / iceAreaCategory(1,iCategory,iCell)

             iceThicknessCategoryMin(iCategory)  = min(iceThicknessCategoryMin(iCategory),  iceThickness)
             iceThicknessCategoryMax(iCategory)  = max(iceThicknessCategoryMax(iCategory),  iceThickness)

             if (snowThickness > snowThicknessCategoryMax(iCategory)) snowThicknessCategoryLocMax(iCategory) = iCell

             snowThicknessCategoryMin(iCategory) = min(snowThicknessCategoryMin(iCategory), snowThickness)
             snowThicknessCategoryMax(iCategory) = max(snowThicknessCategoryMax(iCategory), snowThickness)

          endif

       enddo ! iCategory
    enddo ! iCell

  end subroutine cice_diagnostic_min_max

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_diagnostic_global
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date August 20th 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_diagnostic_global(block, clock)

    use cice_constants, only: &
         iceAreaMinimum

    type(block_type), intent(inout) :: &
         block !< Input/Output: 

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    type(MPAS_pool_type), pointer :: &
         mesh, &
         global_diagnostics, &
         tracers

    real(kind=RKIND), dimension(:,:,:), pointer :: &
        iceAreaCategory, &
        iceVolumeCategory, &
        snowVolumeCategory, &
        surfaceTemperature

    real(kind=RKIND), dimension(:), pointer :: &
         globalIceArea, &
         globalIceExtent, &
         globalIceVolume, &
         globalIceThickness, &
         globalSnowVolume, &
         globalSnowThickness

    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         areaCell

    integer, pointer :: &
        nCells, &
        nCategories

    integer :: &
         iCell, &
         iCategory

    integer, parameter :: &
         iGlobal = 1, &
         iNorthern = 2, &
         iSouthern = 3

    real(kind=RKIND), parameter :: &
         iceExtentThreshold = 0.15_RKIND

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "global_diagnostics", global_diagnostics)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_array(global_diagnostics, "globalIceArea", globalIceArea)
    call MPAS_pool_get_array(global_diagnostics, "globalIceExtent", globalIceExtent)
    call MPAS_pool_get_array(global_diagnostics, "globalIceVolume", globalIceVolume)
    call MPAS_pool_get_array(global_diagnostics, "globalIceThickness", globalIceThickness)
    call MPAS_pool_get_array(global_diagnostics, "globalSnowVolume", globalSnowVolume)
    call MPAS_pool_get_array(global_diagnostics, "globalSnowThickness", globalSnowThickness)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
    
    call MPAS_pool_get_array(mesh, "latCell", latCell)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    globalIceArea(:) = 0.0_RKIND
    globalIceExtent(:) = 0.0_RKIND
    globalIceVolume(:) = 0.0_RKIND
    globalIceThickness(:) = 0.0_RKIND
    globalSnowVolume(:) = 0.0_RKIND
    globalSnowThickness(:) = 0.0_RKIND

    do iCell = 1, nCells
       do iCategory = 1, nCategories

          if (iceAreaCategory(1,iCategory,iCell) > iceAreaMinimum) then

             ! global
             ! global ice area
             globalIceArea(iGlobal) = globalIceArea(iGlobal) + &
                  areaCell(iCell) * iceAreaCategory(1,iCategory,iCell)

             ! global ice extent
             if (iceAreaCategory(1,iCategory,iCell) > iceExtentThreshold) &
                  globalIceExtent(iGlobal) = globalIceExtent(iGlobal) + areaCell(iCell)

             ! global ice volume
             globalIceVolume(iGlobal) = globalIceVolume(iGlobal) + &
                  areaCell(iCell) * iceVolumeCategory(1,iCategory,iCell)

             ! global snow volume
             globalSnowVolume(iGlobal) = globalSnowVolume(iGlobal) + &
                  areaCell(iCell) * snowVolumeCategory(1,iCategory,iCell)

             ! northern hemisphere
             if (latCell(iCell) > 0.0_RKIND) then

                ! northern hemisphere ice area
                globalIceArea(iNorthern) = globalIceArea(iNorthern) + &
                     areaCell(iCell) * iceAreaCategory(1,iCategory,iCell)

                ! northern hemisphere ice extent
                if (iceAreaCategory(1,iCategory,iCell) > iceExtentThreshold) &
                     globalIceExtent(iNorthern) = globalIceExtent(iNorthern) + areaCell(iCell)

                ! northern hemisphere ice volume
                globalIceVolume(iNorthern) = globalIceVolume(iNorthern) + &
                     areaCell(iCell) * iceVolumeCategory(1,iCategory,iCell)
                
                ! northern hemisphere snow volume
                globalSnowVolume(iNorthern) = globalSnowVolume(iNorthern) + &
                     areaCell(iCell) * snowVolumeCategory(1,iCategory,iCell)

             endif ! northern hemisphere

             ! southern hemisphere 
             if (latCell(iCell) < 0.0_RKIND) then

                ! southern hemisphere ice area
                globalIceArea(iSouthern) = globalIceArea(iSouthern) + &
                     areaCell(iCell) * iceAreaCategory(1,iCategory,iCell)

                ! southern hemisphere ice extent
                if (iceAreaCategory(1,iCategory,iCell) > iceExtentThreshold) &
                     globalIceExtent(iSouthern) = globalIceExtent(iSouthern) + areaCell(iCell)

                ! southern hemisphere ice volume
                globalIceVolume(iSouthern) = globalIceVolume(iSouthern) + &
                     areaCell(iCell) * iceVolumeCategory(1,iCategory,iCell)

                ! southern hemisphere snow volume
                globalSnowVolume(iSouthern) = globalSnowVolume(iSouthern) + &
                     areaCell(iCell) * snowVolumeCategory(1,iCategory,iCell)

             endif ! southern hemisphere

          endif

       enddo ! iCategory
    enddo ! iCell

    ! ice thickness
    where (globalIceArea > 0.0_RKIND) &
         globalIceThickness = globalIceVolume / globalIceArea

    ! snow thickness
    where (globalIceArea > 0.0_RKIND) &
         globalSnowThickness = globalSnowVolume / globalIceArea

  end subroutine cice_diagnostic_global

!-----------------------------------------------------------------------

end module cice_diagnostics
