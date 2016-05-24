module subgridMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: subgridMod
!
! !DESCRIPTION:
! sub-grid data and mapping types and modules
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun

  implicit none
  private	
  save

! !PUBLIC MEMBER FUNCTIONS:
  public subgrid_get_gcellinfo        ! Returns g,l,c,p properties from wtxy


! !REVISION HISTORY:
! 2006.07.04 T Craig, rename initSubgridMod
!
!
! !PRIVATE MEMBER FUNCTIONS: None
!
! !PRIVATE DATA MEMBERS: None
!EOP
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subgrid_get_gcellinfo
!
! !INTERFACE:
  subroutine subgrid_get_gcellinfo (nw, &
                             nlunits, ncols, npfts, &
                             nveg, wtveg, &
                             ncrop, wtcrop, &
                             nurban, wturban, &
                             nlake, wtlake, &
                             nwetland, wtwetland, &
                             nglacier, wtglacier, &
                             nglacier_mec, wtglacier_mec,  &
                             glcmask)
!
! !DESCRIPTION:
! Obtain gridcell properties
!
! !USES
  use clm_varpar  , only : numpft, maxpatch_pft, numcft, &
                           npatch_lake, npatch_glacier, npatch_wet, npatch_urban
  use clm_varpar  , only : npatch_glacier_mec
  use clm_varctl  , only : allocate_all_vegpfts, create_crop_landunit
  use clm_varctl  , only : create_glacier_mec_landunit, glc_topomax
  use clm_varsur  , only : wtxy
  use clm_varsur  , only : topoxy

! !ARGUMENTS
    implicit none
    integer , intent(in)  :: nw                   ! wtxy cell index
    integer , optional, intent(out) :: nlunits    ! number of landunits
    integer , optional, intent(out) :: ncols      ! number of columns 
    integer , optional, intent(out) :: npfts      ! number of pfts 
    integer , optional, intent(out) :: nveg       ! number of vegetated pfts in naturally vegetated landunit
    real(r8), optional, intent(out) :: wtveg      ! weight (relative to gridcell) of naturally vegetated landunit
    integer , optional, intent(out) :: ncrop      ! number of crop pfts in crop landunit
    real(r8), optional, intent(out) :: wtcrop     ! weight (relative to gridcell) of crop landunit
    integer , optional, intent(out) :: nurban     ! number of urban pfts (columns) in urban landunit
    real(r8), optional, intent(out) :: wturban    ! weight (relative to gridcell) of urban pfts (columns) in urban la
    integer , optional, intent(out) :: nlake      ! number of lake pfts (columns) in lake landunit
    real(r8), optional, intent(out) :: wtlake     ! weight (relative to gridcell) of lake landunitof lake pfts (columns) in lake landunit
    integer , optional, intent(out) :: nwetland   ! number of wetland pfts (columns) in wetland landunit
    real(r8), optional, intent(out) :: wtwetland  ! weight (relative to gridcell) of wetland landunitof wetland pfts (columns) in wetland landunit
    integer , optional, intent(out) :: nglacier   ! number of glacier pfts (columns) in glacier landunit
    real(r8), optional, intent(out) :: wtglacier  ! weight (relative to gridcell) of glacier landunitof glacier pfts (columns) in glacier landunit
    integer , optional, intent(out) :: nglacier_mec  ! number of glacier_mec pfts (columns) in glacier_mec landunit
    real(r8), optional, intent(out) :: wtglacier_mec ! weight (relative to gridcell) of glacier_mec landunitof glacier pfts (columns) in glacier_mec landunit
    integer , optional, intent(in)  :: glcmask  ! = 1 if glc requires surface mass balance in this gridcell
!
! !CALLED FROM:
! subroutines decomp_init, initGridCells
!
! !REVISION HISTORY:
! 2002.09.11  Mariana Vertenstein  Creation.
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: m                ! loop index
    integer  :: n                ! elevation class index
    integer  :: ipfts            ! number of pfts in gridcell
    integer  :: icols            ! number of columns in gridcell
    integer  :: ilunits          ! number of landunits in gridcell
    integer  :: npfts_per_lunit  ! number of pfts in landunit
    real(r8) :: wtlunit          ! weight (relative to gridcell) of landunit
!------------------------------------------------------------------------------

    ! Initialize pfts, columns and landunits counters for gridcell

    ipfts   = 0
    icols   = 0
    ilunits = 0

    ! Set naturally vegetated landunit

    npfts_per_lunit = 0
    wtlunit = 0._r8
    ! If crop should be on separate land units
    if (allocate_all_vegpfts .and. create_crop_landunit) then
       do m = 1, maxpatch_pft-numcft
          if (wtxy(nw,m) > 0.0_r8) then
             npfts_per_lunit = npfts_per_lunit + 1 ! sum natural pfts
             wtlunit = wtlunit + wtxy(nw,m)        ! and their wts
          end if
       end do
       do m = maxpatch_pft-numcft+1, maxpatch_pft
          if (wtxy(nw,m) > 0.0_r8) then
             npfts_per_lunit = npfts_per_lunit + 1 ! sum crops, too, but not
          end if                                   ! their wts for now
       end do
    ! Assume that the vegetated landunit has one column
    else
       do m = 1, maxpatch_pft            
          if (wtxy(nw,m) > 0.0_r8) then
             npfts_per_lunit = npfts_per_lunit + 1
             wtlunit = wtlunit + wtxy(nw,m)
          end if
       end do
    end if
    if (npfts_per_lunit > 0) then ! true even when only crops are present
       if (allocate_all_vegpfts) npfts_per_lunit = numpft+1
       if (allocate_all_vegpfts .and. create_crop_landunit) npfts_per_lunit = numpft+1-numcft
       ilunits = ilunits + 1
       icols = icols + 1  
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nveg )) nveg  = npfts_per_lunit
    if (present(wtveg)) wtveg = wtlunit

    ! Set urban landunit

    npfts_per_lunit = 0
    wtlunit = 0._r8
    do m = npatch_urban, npatch_lake-1
       if (wtxy(nw,m) > 0.0_r8) then
          npfts_per_lunit = npfts_per_lunit + 1
          wtlunit = wtlunit + wtxy(nw,m)
       end if
    end do
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nurban )) nurban  = npfts_per_lunit
    if (present(wturban)) wturban = wtlunit

    ! Set lake landunit

    npfts_per_lunit = 0
    wtlunit = 0._r8
    if (wtxy(nw,npatch_lake) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
       wtlunit = wtlunit + wtxy(nw,npatch_lake)
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nlake )) nlake  = npfts_per_lunit
    if (present(wtlake)) wtlake = wtlunit

    ! Set wetland landunit

    npfts_per_lunit = 0
    wtlunit = 0._r8
    if (wtxy(nw,npatch_wet) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
       wtlunit = wtlunit + wtxy(nw,npatch_wet)
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nwetland )) nwetland  = npfts_per_lunit
    if (present(wtwetland)) wtwetland = wtlunit

    ! Set glacier landunit

    npfts_per_lunit = 0
    wtlunit = 0._r8
    if (wtxy(nw,npatch_glacier) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
       wtlunit = wtlunit + wtxy(nw,npatch_glacier)
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nglacier )) nglacier  = npfts_per_lunit
    if (present(wtglacier)) wtglacier = wtlunit

    ! Set glacier_mec landunit
    ! If glcmask = 1, we create a column for each elevation class even if wtxy = 0.

    if (create_glacier_mec_landunit) then
       npfts_per_lunit = 0
       wtlunit = 0._r8
       do m = npatch_glacier+1, npatch_glacier_mec
          if (wtxy(nw,m) > 0._r8) then
             npfts_per_lunit = npfts_per_lunit + 1
             wtlunit = wtlunit + wtxy(nw,m)
             topoxy(nw,m) = max (topoxy(nw,m), 0._r8)
          elseif (present(glcmask)) then
             if (glcmask == 1) then      ! create a virtual column 
                npfts_per_lunit = npfts_per_lunit + 1
                n = m - npatch_glacier    ! elevation class index
                if (m < npatch_glacier_mec) then   ! classes 1 to maxpatch_glcmec-1 
                   topoxy(nw,m) = 0.5_r8 * (glc_topomax(n-1) + glc_topomax(n))
                else                               ! class maxpatch_glcmec
                   topoxy(nw,m) = 2.0_r8*glc_topomax(n-1) - glc_topomax(n-2)     ! somewhat arbitrary
                endif 
             endif  ! glcmask = 1 
          endif  ! wtxy > 0
       enddo   ! npatch_glacier_mec
       if (npfts_per_lunit > 0) then
          ilunits = ilunits + 1
          icols   = icols + npfts_per_lunit
       end if
       ipfts = ipfts + npfts_per_lunit
       if (present(nglacier_mec )) nglacier_mec  = npfts_per_lunit
       if (present(wtglacier_mec)) wtglacier_mec = wtlunit

    endif    ! create_glacier_mec_landunit

    ! Set crop landunit if appropriate

    npfts_per_lunit = 0
    wtlunit = 0._r8
    if (allocate_all_vegpfts .and. create_crop_landunit) then
       do m = 1, maxpatch_pft-numcft
          if (wtxy(nw,m) > 0.0_r8) then
             npfts_per_lunit = npfts_per_lunit + 1 ! sum natural pfts again
          end if                                   ! not their wts this time
       end do
       do m = maxpatch_pft-numcft+1, maxpatch_pft
          if (wtxy(nw,m) > 0.0_r8) then
             npfts_per_lunit = npfts_per_lunit + 1 ! sum crops
             wtlunit = wtlunit + wtxy(nw,m)        ! and their wts
          end if
       end do
    end if
    if (npfts_per_lunit > 0) then ! true even if only natural veg is present
       if (allocate_all_vegpfts .and. create_crop_landunit) npfts_per_lunit = numcft
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(ncrop )) ncrop  = npfts_per_lunit
    if (present(wtcrop)) wtcrop = wtlunit

    ! Determine return arguments

    if (present(nlunits)) nlunits = ilunits
    if (present(ncols))   ncols   = icols
    if (present(npfts))   npfts   = ipfts

  end subroutine subgrid_get_gcellinfo

!-----------------------------------------------------------------------

end module subgridMod
