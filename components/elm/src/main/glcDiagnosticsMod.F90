module glcDiagnosticsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Computes and outputs a number of glacier-related diagnostic quantities
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type 
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: glc_diagnostics_type
     private

     real(r8), pointer :: gris_mask_grc        (:) ! Greenland ice sheet mask 
     real(r8), pointer :: gris_area_grc        (:) ! Greenland ice-covered area per gridcell (km^2)
     real(r8), pointer :: aais_mask_grc        (:) ! Antarctic ice sheet mask 
     real(r8), pointer :: aais_area_grc        (:) ! Antarctic ice-covered area per gridcell (km^2)

   contains
     
     procedure, public  :: Init

     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: calc_timeconst_diagnostics ! calculate time-constant diagnostics (which can be computed once, at initialization)

  end type glc_diagnostics_type

  !------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    ! !ARGUMENTS:
    class(glc_diagnostics_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------
    
    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)

    call this%calc_timeconst_diagnostics(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    ! !ARGUMENTS:
    class(glc_diagnostics_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! LOCAL VARAIBLES:
    integer  :: begg, endg
    
    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------
    
    begg = bounds%begg; endg=bounds%endg

    allocate(this%gris_mask_grc (begg:endg)) ; this%gris_mask_grc (:) = nan
    allocate(this%gris_area_grc (begg:endg)) ; this%gris_area_grc (:) = nan
    allocate(this%aais_mask_grc (begg:endg)) ; this%aais_mask_grc (:) = nan
    allocate(this%aais_area_grc (begg:endg)) ; this%aais_area_grc (:) = nan

  end subroutine InitAllocate


  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use clm_varctl     , only : create_glacier_mec_landunit
    use histFileMod    , only: hist_addfld1d
    use elm_varcon     , only : spval
    !
    ! !ARGUMENTS:
    class(glc_diagnostics_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begg, endg    

    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------
    
    begg = bounds%begg
    endg = bounds%endg

    if (create_glacier_mec_landunit) then

       this%gris_mask_grc(begg:endg) = spval
       call hist_addfld1d (fname='gris_mask',  units='unitless',  &
            avgflag='A', long_name='Greenland mask', &
            ptr_gcell= this%gris_mask_grc)

       this%gris_area_grc(begg:endg) = spval
       call hist_addfld1d (fname='gris_area',  units='km^2',  &
            avgflag='A', long_name='Greenland ice area', &
            ptr_gcell= this%gris_area_grc)

       this%aais_mask_grc(begg:endg) = spval
       call hist_addfld1d (fname='aais_mask',  units='unitless',  &
            avgflag='A', long_name='Antarctic mask', &
            ptr_gcell= this%aais_mask_grc)

       this%aais_area_grc(begg:endg) = spval
       call hist_addfld1d (fname='aais_area',  units='km^2',  &
            avgflag='A', long_name='Antarctic ice area', &
            ptr_gcell= this%aais_area_grc)

    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine calc_timeconst_diagnostics(this, bounds)
    !
    ! !DESCRIPTION:
    ! Calculate time-constant glacier-related diagnostic fields (this only needs to be
    ! called once, at initialization)
    !
    ! !USES:
    use domainMod, only : ldomain
    !
    ! !ARGUMENTS:
    class(glc_diagnostics_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: gdc  ! gridcell index
    
    character(len=*), parameter :: subname = 'calc_timeconst_diagnostics'
    !-----------------------------------------------------------------------
    
    do gdc = bounds%begg, bounds%endg

       ! Make ice sheet masks

       this%gris_mask_grc(gdc) = 0._r8
       this%gris_area_grc(gdc) = 0._r8
       this%aais_mask_grc(gdc) = 0._r8
       this%aais_area_grc(gdc) = 0._r8

       ! Greenland mask
       if ( (ldomain%latc(gdc) >  58. .and. ldomain%latc(gdc) <= 67.  .and.   &
            ldomain%lonc(gdc) > 302. .and. ldomain%lonc(gdc) < 330.)         &
            .or.                                 &
            (ldomain%latc(gdc) >  67. .and. ldomain%latc(gdc) <= 70. .and.    &
            ldomain%lonc(gdc) > 300. .and. ldomain%lonc(gdc) < 345.)         &
            .or.                                 &
            (ldomain%latc(gdc) >  70. .and. ldomain%latc(gdc) <= 75. .and.    &
            ldomain%lonc(gdc) > 295. .and. ldomain%lonc(gdc) < 350.)         &
            .or.                                 &
            (ldomain%latc(gdc) >  75. .and. ldomain%latc(gdc) <= 79. .and.    &
            ldomain%lonc(gdc) > 285. .and. ldomain%lonc(gdc) < 350.)         &
            .or.                                 &
            (ldomain%latc(gdc) >  79. .and. ldomain%latc(gdc) <  85. .and.    &
            ldomain%lonc(gdc) > 290. .and. ldomain%lonc(gdc) < 355.) ) then

          this%gris_mask_grc(gdc) = 1.0_r8

       elseif (ldomain%latc(gdc) < -60.) then

          this%aais_mask_grc(gdc) = 1.0_r8

       endif  ! Greenland or Antarctic grid cell

    end do

  end subroutine calc_timeconst_diagnostics



end module glcDiagnosticsMod
