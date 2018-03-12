! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_unit_conversion
!
!> \brief MPAS sea ice analysis mode member: unit_conversion
!> \author Adrian K. Turner
!> \date   19th February 2016
!> \details
!>
!-----------------------------------------------------------------------

module cice_unit_conversion

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

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

   public :: cice_bootstrap_unit_conversion, &
             cice_init_unit_conversion, &
             cice_precompute_unit_conversion, &
             cice_compute_unit_conversion, &
             cice_restart_unit_conversion, &
             cice_finalize_unit_conversion

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine cice_bootstrap_unit_conversion
!
!> \brief   Bootstrap MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    19th February 2016
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_bootstrap_unit_conversion(domain, instance, err)!{{{

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

   end subroutine cice_bootstrap_unit_conversion!}}}

!***********************************************************************
!
!  routine cice_init_unit_conversion
!
!> \brief   Initialize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    19th February 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_init_unit_conversion(domain, instance, err)!{{{

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

   end subroutine cice_init_unit_conversion!}}}

!***********************************************************************
!
!  routine cice_precompute_unit_conversion
!
!> \brief   Precompute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    19th February 2016
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_precompute_unit_conversion(domain, instance, timeLevel, err)!{{{

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

   end subroutine cice_precompute_unit_conversion!}}}

!***********************************************************************
!
!  routine cice_compute_unit_conversion
!
!> \brief   Compute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    19th February 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_compute_unit_conversion(domain, instance, timeLevel, err)!{{{

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

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           meltGrowthRatesPool

      real(kind=RKIND), dimension(:), pointer :: &
           surfaceIceMelt, &
           basalIceMelt, &
           lateralIceMelt, &
           snowMelt, &
           congelation, &
           snowiceFormation, &
           snowThicknessChange, &
           frazilFormation

      real(kind=RKIND), dimension(:,:), pointer :: &
           surfaceIceMeltCategory, &
           basalIceMeltCategory, &
           snowMeltCategory, &
           congelationCategory, &
           snowiceFormationCategory, &
           snowThicknessChangeCategory

      integer :: &
           iCell

      integer, pointer :: &
           nCellsSolve

      real(kind=RKIND), pointer :: &
           config_dt

      real(kind=RKIND) :: &
           meltGrowthRateConversionFactor

      block => domain % blocklist
      do while(associated(block))

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", meltGrowthRatesPool)

         call MPAS_pool_get_array(meltGrowthRatesPool, "surfaceIceMelt", surfaceIceMelt)
         call MPAS_pool_get_array(meltGrowthRatesPool, "surfaceIceMeltCategory", surfaceIceMeltCategory)
         call MPAS_pool_get_array(meltGrowthRatesPool, "basalIceMelt", basalIceMelt)
         call MPAS_pool_get_array(meltGrowthRatesPool, "basalIceMeltCategory", basalIceMeltCategory)
         call MPAS_pool_get_array(meltGrowthRatesPool, "lateralIceMelt", lateralIceMelt)
         call MPAS_pool_get_array(meltGrowthRatesPool, "snowMelt", snowMelt)
         call MPAS_pool_get_array(meltGrowthRatesPool, "snowMeltCategory", snowMeltCategory)
         call MPAS_pool_get_array(meltGrowthRatesPool, "congelation", congelation)
         call MPAS_pool_get_array(meltGrowthRatesPool, "congelationCategory", congelationCategory)
         call MPAS_pool_get_array(meltGrowthRatesPool, "snowiceFormation", snowiceFormation)
         call MPAS_pool_get_array(meltGrowthRatesPool, "snowiceFormationCategory", snowiceFormationCategory)
         call MPAS_pool_get_array(meltGrowthRatesPool, "snowThicknessChange", snowThicknessChange)
         call MPAS_pool_get_array(meltGrowthRatesPool, "snowThicknessChangeCategory", snowThicknessChangeCategory)
         call MPAS_pool_get_array(meltGrowthRatesPool, "frazilFormation", frazilFormation)

         call MPAS_pool_get_config(block % configs, "config_dt", config_dt)

         ! convert growth/melt rates from m/timestep to m/s
         meltGrowthRateConversionFactor = 1.0_RKIND / config_dt

         do iCell = 1, nCellsSolve

            surfaceIceMelt(iCell)                = surfaceIceMelt(iCell)                * meltGrowthRateConversionFactor
            surfaceIceMeltCategory(:,iCell)      = surfaceIceMeltCategory(:,iCell)      * meltGrowthRateConversionFactor
            basalIceMelt(iCell)                  = basalIceMelt(iCell)                  * meltGrowthRateConversionFactor
            basalIceMeltCategory(:,iCell)        = basalIceMeltCategory(:,iCell)        * meltGrowthRateConversionFactor
            lateralIceMelt(iCell)                = lateralIceMelt(iCell)                * meltGrowthRateConversionFactor
            snowMelt(iCell)                      = snowMelt(iCell)                      * meltGrowthRateConversionFactor
            snowMeltCategory(:,iCell)            = snowMeltCategory(:,iCell)            * meltGrowthRateConversionFactor
            congelation(iCell)                   = congelation(iCell)                   * meltGrowthRateConversionFactor
            congelationCategory(:,iCell)         = congelationCategory(:,iCell)         * meltGrowthRateConversionFactor
            snowiceFormation(iCell)              = snowiceFormation(iCell)              * meltGrowthRateConversionFactor
            snowiceFormationCategory(:,iCell)    = snowiceFormationCategory(:,iCell)    * meltGrowthRateConversionFactor
            snowThicknessChange(iCell)           = snowThicknessChange(iCell)           * meltGrowthRateConversionFactor
            snowThicknessChangeCategory(:,iCell) = snowThicknessChangeCategory(:,iCell) * meltGrowthRateConversionFactor
            frazilFormation(iCell)               = frazilFormation(iCell)               * meltGrowthRateConversionFactor

         enddo ! iCell

         block => block % next
      enddo

   end subroutine cice_compute_unit_conversion!}}}

!***********************************************************************
!
!  routine cice_restart_unit_conversion
!
!> \brief   Save restart for MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    19th February 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_restart_unit_conversion(domain, instance, err)!{{{

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

   end subroutine cice_restart_unit_conversion!}}}

!***********************************************************************
!
!  routine cice_finalize_unit_conversion
!
!> \brief   Finalize MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    19th February 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_finalize_unit_conversion(domain, instance, err)!{{{

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

   end subroutine cice_finalize_unit_conversion!}}}

!-----------------------------------------------------------------------

end module cice_unit_conversion

! vim: foldmethod=marker
