!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date
!> \details
!>
!
!-----------------------------------------------------------------------

module cice_history

  use mpas_derived_types
  use mpas_pool_routines

  implicit none

  private
  save

  public :: &
       cice_init_history_thermodynamics, &
!nj-bgc+
       cice_init_history_biogeochemistry
!nj-bgc-

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_init_history_thermodynamics(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         atmosFluxes, &
         meltGrowthRates, &
         shortwave, &
         ponds

    real(kind=RKIND), dimension(:), pointer :: &
         surfaceHeatFlux, &
         surfaceConductiveFlux, &
         congelation, &
         frazilFormation, &
         snowiceFormation, &
         surfaceIceMelt, &
         snowMelt, &
         basalIceMelt, &
         lateralIceMelt, &
         pondFreshWaterFlux, &
         bareIceAlbedoCell, &
         snowAlbedoCell, &
         pondAlbedoCell, &
         snowThicknessChange

    real(kind=RKIND), dimension(:,:), pointer :: &
         surfaceHeatFluxCategory, &
         surfaceConductiveFluxCategory, &
         latentHeatFluxCategory

         
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxes)
       call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", meltGrowthRates)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)

       call MPAS_pool_get_array(atmosFluxes, "surfaceHeatFlux", surfaceHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, "surfaceConductiveFlux", surfaceConductiveFlux)
       call MPAS_pool_get_array(atmosFluxes, "surfaceHeatFluxCategory", surfaceHeatFluxCategory)
       call MPAS_pool_get_array(atmosFluxes, "surfaceConductiveFluxCategory", surfaceConductiveFluxCategory)
       call MPAS_pool_get_array(atmosFluxes, "latentHeatFluxCategory", latentHeatFluxCategory)

       call MPAS_pool_get_array(shortwave, "bareIceAlbedoCell", bareIceAlbedoCell)
       call MPAS_pool_get_array(shortwave, "snowAlbedoCell", snowAlbedoCell)
       call MPAS_pool_get_array(shortwave, "pondAlbedoCell", pondAlbedoCell)

       call MPAS_pool_get_array(ponds, "pondFreshWaterFlux", pondFreshWaterFlux)

       call MPAS_pool_get_array(meltGrowthRates, "congelation", congelation)
       call MPAS_pool_get_array(meltGrowthRates, "frazilFormation", frazilFormation)
       call MPAS_pool_get_array(meltGrowthRates, "snowiceFormation", snowiceFormation)
       call MPAS_pool_get_array(meltGrowthRates, "snowThicknessChange", snowThicknessChange)
       call MPAS_pool_get_array(meltGrowthRates, "surfaceIceMelt", surfaceIceMelt)
       call MPAS_pool_get_array(meltGrowthRates, "snowMelt", snowMelt)
       call MPAS_pool_get_array(meltGrowthRates, "basalIceMelt", basalIceMelt)
       call MPAS_pool_get_array(meltGrowthRates, "lateralIceMelt", lateralIceMelt)

       surfaceHeatFlux = 0.0_RKIND
       surfaceConductiveFlux = 0.0_RKIND
       surfaceHeatFluxCategory = 0.0_RKIND
       surfaceConductiveFluxCategory = 0.0_RKIND
       latentHeatFluxCategory = 0.0_RKIND
       congelation = 0.0_RKIND
       frazilFormation = 0.0_RKIND
       snowiceFormation = 0.0_RKIND
       snowThicknessChange = 0.0_RKIND
       surfaceIceMelt = 0.0_RKIND
       snowMelt = 0.0_RKIND
       basalIceMelt = 0.0_RKIND
       lateralIceMelt = 0.0_RKIND
       pondFreshWaterFlux = 0.0_RKIND
       bareIceAlbedoCell = 0.0_RKIND
       snowAlbedoCell = 0.0_RKIND
       pondAlbedoCell = 0.0_RKIND
       
       block => block % next
    enddo

  end subroutine cice_init_history_thermodynamics

!nj-bgc+
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!> \brief
!> \author Nicole Jeffery LANL
!> \date   21 September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_init_history_biogeochemistry(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         biogeochemistry

    real(kind=RKIND), dimension(:), pointer :: &
         primaryProduction, &
         netSpecificAlgalGrowthRate, &
         netBrineHeight, &
         zSalinityFlux, &
         zSalinityGDFlux
         
    real(kind=RKIND), dimension(:,:), pointer :: &
         oceanBioFluxes, &
         atmosIceBioFluxes, &
         snowIceBioFluxes, &
         totalVerticalBiologyIce, &
         totalVerticalBiologySnow

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         bioTracerShortwave

         
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

       call MPAS_pool_get_array(biogeochemistry, "primaryProduction", primaryProduction)
       call MPAS_pool_get_array(biogeochemistry, "netSpecificAlgalGrowthRate", netSpecificAlgalGrowthRate)
       call MPAS_pool_get_array(biogeochemistry, "netBrineHeight", netBrineHeight)
       call MPAS_pool_get_array(biogeochemistry, "zSalinityFlux", zSalinityFlux)
       call MPAS_pool_get_array(biogeochemistry, "zSalinityGDFlux", zSalinityGDFlux)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "atmosIceBioFluxes", atmosIceBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "snowIceBioFluxes", snowIceBioFluxes)
       call MPAS_pool_get_array(biogeochemistry, "totalVerticalBiologyIce", totalVerticalBiologyIce)
       call MPAS_pool_get_array(biogeochemistry, "totalVerticalBiologySnow", totalVerticalBiologySnow)
       call MPAS_pool_get_array(biogeochemistry, "bioTracerShortwave", bioTracerShortwave)

       primaryProduction = 0.0_RKIND
       netSpecificAlgalGrowthRate = 0.0_RKIND
       netBrineHeight = 0.0_RKIND
       zSalinityFlux = 0.0_RKIND
       zSalinityGDFlux = 0.0_RKIND
       oceanBioFluxes = 0.0_RKIND
       atmosIceBioFluxes = 0.0_RKIND
       snowIceBioFluxes = 0.0_RKIND
       totalVerticalBiologyIce = 0.0_RKIND
       totalVerticalBiologySnow = 0.0_RKIND
       bioTracerShortwave = 0.0_RKIND
       
       block => block % next
    enddo

  end subroutine cice_init_history_biogeochemistry
!nj-bgc-
!-----------------------------------------------------------------------

end module cice_history

