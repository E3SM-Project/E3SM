module subgridMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! sub-grid data and mapping types and modules
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun
  use elm_varctl  , only : iulog

  implicit none
  private   
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public subgrid_get_gcellinfo        ! Obtain gridcell properties
  public subgrid_get_topounitinfo     ! Obtain topounit properties
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine subgrid_get_gcellinfo (gi, &
       ntopounits, nlunits, ncols, npfts, ncohorts, &
       nveg, &
       ncrop, &
       nurban_tbd, &
       nurban_hd, &
       nurban_md, &
       nlake, &
       nwetland, &
       nglacier, &
       nglacier_mec,  &
       glcmask)
    !
    ! !DESCRIPTION:
    ! Obtain gridcell properties
    !
    ! !USES
    use elm_varpar  , only : natpft_size, cft_size, maxpatch_urb, maxpatch_glcmec
    use elm_varctl  , only : create_crop_landunit
    use clm_varsur  , only : wt_lunit, urban_valid, wt_glc_mec
    use landunit_varcon  , only : istsoil, istcrop, istice, istice_mec, istdlak, istwet, &
                             isturb_tbd, isturb_hd, isturb_md
    use topounit_varcon  , only : max_topounits
    use FatesInterfaceTypesMod, only : fates_maxElementsPerSite

    !
    ! !ARGUMENTS
    implicit none
    integer , intent(in)  :: gi                   ! grid cell index
    integer , optional, intent(out) :: ntopounits ! number of topographic units
    integer , optional, intent(out) :: nlunits    ! number of landunits
    integer , optional, intent(out) :: ncols      ! number of columns 
    integer , optional, intent(out) :: npfts      ! number of pfts 
    integer , optional, intent(out) :: ncohorts   ! number of cohorts 
    integer , optional, intent(out) :: nveg       ! number of vegetated pfts in naturally vegetated landunit
    integer , optional, intent(out) :: ncrop      ! number of crop pfts in crop landunit
    integer , optional, intent(out) :: nurban_tbd ! number of urban pfts (columns) in urban TBD landunit
    integer , optional, intent(out) :: nurban_hd ! number of urban pfts (columns) in urban HD landunit
    integer , optional, intent(out) :: nurban_md ! number of urban pfts (columns) in urban MD landunit
    integer , optional, intent(out) :: nlake      ! number of lake pfts (columns) in lake landunit
    integer , optional, intent(out) :: nwetland   ! number of wetland pfts (columns) in wetland landunit
    integer , optional, intent(out) :: nglacier   ! number of glacier pfts (columns) in glacier landunit
    integer , optional, intent(out) :: nglacier_mec  ! number of glacier_mec pfts (columns) in glacier_mec landunit
    integer , optional, intent(in)  :: glcmask  ! = 1 if glc requires surface mass balance in this gridcell
    !
    ! !LOCAL VARIABLES:
    integer  :: t, m             ! loop index
    integer  :: n                ! elevation class index
    integer  :: itopounits       ! number of topographic units in gridcell
    integer  :: ilunits          ! number of landunits in gridcell
    integer  :: icols            ! number of columns in gridcell
    integer  :: ipfts            ! number of pfts in gridcell
    integer  :: icohorts         ! number of cohorts in gridcell
    integer  :: npfts_per_lunit  ! number of pfts in landunit
    integer  :: ntopounits_per_gcell  ! number of topounits in this gridcell
    !------------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Initialize topounits, landunits, pfts, columns and cohorts counters for gridcell
    ! -------------------------------------------------------------------------

    itopounits = 0
    ilunits    = 0
    icols      = 0
    ipfts      = 0
    icohorts   = 0
    
    ! -------------------------------------------------------------------------
    ! Initialize other optional counters for gridcell
    ! -------------------------------------------------------------------------
    if (present(nveg)) nveg = 0
    if (present(nurban_tbd)) nurban_tbd = 0
    if (present(nurban_hd)) nurban_hd = 0
    if (present(nurban_md)) nurban_md = 0
    if (present(nlake)) nlake = 0
    if (present(nwetland)) nwetland = 0
    if (present(nglacier)) nglacier = 0
    if (present(nglacier_mec)) nglacier_mec = 0
    if (present(ncrop)) ncrop = 0
    
    ! -------------------------------------------------------------------------
    ! Loop through topographic units
    ! -------------------------------------------------------------------------
    
    ! as a start, everything that used to be allocated on a gridcell will now be allocated
    ! for each topographic unit on each gridcell. That will be wasteful in many cases, but 
    ! is a simple way to make the transition to topounits. We should return to this to see
    ! if there is a smarter way to allocate space that does not preclude land area transitions
    ! of interest.
    ntopounits_per_gcell = max_topounits  ! this will be replaced later with a constant > 1, or with a function call
                                          ! that sets the number of topounits uniquely for each gridcell.
    do t = 1, ntopounits_per_gcell
       
       itopounits = itopounits + 1

       ! -------------------------------------------------------------------------
       ! Set naturally vegetated landunit
       ! -------------------------------------------------------------------------

       ! To support dynamic landunits, we have a naturally vegetated landunit in every grid
       ! cell, because it might need to come into existence even if its weight is 0 at the
       ! start of the run. And to support transient pfts or dynamic vegetation, we always
       ! allocate space for ALL PFTs on this landunit.

       npfts_per_lunit = natpft_size

       ! Assume that the vegetated landunit has one column
       ilunits = ilunits + 1
       icols = icols + 1  

       ipfts = ipfts + npfts_per_lunit

       ! -------------------------------------------------------------------------
       ! Number of cohorts is set here
       ! ED cohorts (via FATES) populate all natural vegetation columns.
       ! Current implementations mostly assume that only one column contains
       ! natural vegetation, which is synonymous with the soil column. 
       ! Because of the single natural vegetation paradigm currently in place,
       ! fates cohort allocations are based off a multiplier on grid-cells.
       ! -------------------------------------------------------------------------

       icohorts = fates_maxElementsPerSite

       if (present(nveg)) nveg = nveg + npfts_per_lunit

       ! -------------------------------------------------------------------------
       ! Set urban landunits
       ! -------------------------------------------------------------------------

       ! To support dynamic landunits, we have all urban landunits in every grid cell that
       ! has valid urban parameters, because they might need to come into existence even if
       ! their weight is 0 at the start of the run. And for simplicity, we always allocate
       ! space for ALL columns on the urban landunits.

       ! Set urban tall building district landunit

       npfts_per_lunit = 0
       if (urban_valid(gi)) then
          npfts_per_lunit = maxpatch_urb
          ilunits = ilunits + 1
          icols   = icols + npfts_per_lunit
          ipfts = ipfts + npfts_per_lunit
       end if
       if (present(nurban_tbd)) nurban_tbd = nurban_tbd + npfts_per_lunit

       ! Set urban high density landunit

       npfts_per_lunit = 0
       if (urban_valid(gi)) then
          npfts_per_lunit = maxpatch_urb
          ilunits = ilunits + 1
          icols   = icols + npfts_per_lunit
          ipfts = ipfts + npfts_per_lunit
       end if
       if (present(nurban_hd)) nurban_hd = nurban_hd + npfts_per_lunit

       ! Set urban medium density landunit

       npfts_per_lunit = 0
       if (urban_valid(gi)) then
          npfts_per_lunit = maxpatch_urb
          ilunits = ilunits + 1
          icols   = icols + npfts_per_lunit
          ipfts = ipfts + npfts_per_lunit
       end if
       if (present(nurban_md)) nurban_md = nurban_md + npfts_per_lunit

       ! -------------------------------------------------------------------------
       ! Set lake landunit
       ! -------------------------------------------------------------------------

       ! We currently do NOT allow the lake landunit to expand via dynamic landunits, so we
       ! only need to allocate space for it where its weight is currently non-zero.

       npfts_per_lunit = 0
       if (wt_lunit(gi, istdlak) > 0.0_r8) then
          npfts_per_lunit = npfts_per_lunit + 1
       end if
       if (npfts_per_lunit > 0) then
          ilunits = ilunits + 1
          icols   = icols + npfts_per_lunit
       end if
       ipfts = ipfts + npfts_per_lunit
       if (present(nlake)) nlake = nlake + npfts_per_lunit

       ! -------------------------------------------------------------------------
       ! Set wetland landunit
       ! -------------------------------------------------------------------------

       ! We currently do NOT allow the wetland landunit to expand via dynamic landunits, so
       ! we only need to allocate space for it where its weight is currently non-zero.

       npfts_per_lunit = 0
       if (wt_lunit(gi, istwet) > 0.0_r8) then
          npfts_per_lunit = npfts_per_lunit + 1
       end if
       if (npfts_per_lunit > 0) then
          ilunits = ilunits + 1
          icols   = icols + npfts_per_lunit
       end if
       ipfts = ipfts + npfts_per_lunit
       if (present(nwetland)) nwetland = nwetland + npfts_per_lunit

       ! -------------------------------------------------------------------------
       ! Set glacier landunit
       ! -------------------------------------------------------------------------

       ! We currently do NOT allow the glacier landunit to expand via dynamic landunits, so
       ! we only need to allocate space for it where its weight is currently non-zero. (If we
       ! have dynamic glacier area, we will be using glacier_mec landunits rather than
       ! glacier landunits.)

       npfts_per_lunit = 0
       if (wt_lunit(gi, istice) > 0.0_r8) then
          npfts_per_lunit = npfts_per_lunit + 1
       end if
       if (npfts_per_lunit > 0) then
          ilunits = ilunits + 1
          icols   = icols + npfts_per_lunit
       end if
       ipfts = ipfts + npfts_per_lunit
       if (present(nglacier)) nglacier = nglacier + npfts_per_lunit

       ! -------------------------------------------------------------------------
       ! Set glacier_mec landunit
       ! -------------------------------------------------------------------------

       ! If glcmask = 1, we create a column for each elevation class even if the weight on
       ! the grid cell is 0. This is needed for coupling to CISM. In addition, this is
       ! currently sufficient to ensure that we have glaciers everywhere they might be
       ! needed with dynamic landunits, since CISM won't be able to create glaciers outside
       ! of the area specified by glcmask. 

       npfts_per_lunit = 0
       do m = 1, maxpatch_glcmec
          ! If the landunit has non-zero weight on the grid cell, and this column has
          ! non-zero weight on the landunit...
          if (wt_lunit(gi, istice_mec) > 0.0_r8 .and. wt_glc_mec(gi, m) > 0.0_r8) then
             npfts_per_lunit = npfts_per_lunit + 1

          elseif (present(glcmask)) then
             if (glcmask == 1) then      ! create a virtual column 
                npfts_per_lunit = npfts_per_lunit + 1
             endif  ! glcmask = 1 
          endif  ! wt > 0
       enddo   ! maxpatch_glcmec
       if (npfts_per_lunit > 0) then
          ilunits = ilunits + 1
          icols   = icols + npfts_per_lunit
       end if
       ipfts = ipfts + npfts_per_lunit
       if (present(nglacier_mec)) nglacier_mec = nglacier_mec + npfts_per_lunit

       ! -------------------------------------------------------------------------
       ! Set crop landunit if appropriate
       ! -------------------------------------------------------------------------

       npfts_per_lunit = 0
       if (create_crop_landunit) then
          ! To support dynamic landunits, we have a crop landunit in every grid cell (if
          ! create_crop_landunit is true), because it might need to come into existence even
          ! if its weight is 0 at the start of the run.
          npfts_per_lunit = cft_size
          ilunits = ilunits + 1
          icols   = icols + npfts_per_lunit
          ipfts = ipfts + npfts_per_lunit
       end if
       if (present(ncrop)) ncrop = ncrop + npfts_per_lunit
    
    enddo ! end of ntopounits loop

    ! -------------------------------------------------------------------------
    ! Determine return arguments
    ! -------------------------------------------------------------------------

    if (present(ntopounits)) ntopounits = itopounits
    if (present(nlunits))       nlunits = ilunits
    if (present(ncols))           ncols = icols
    if (present(npfts))           npfts = ipfts
    if (present(ncohorts))     ncohorts = icohorts

  end subroutine subgrid_get_gcellinfo

  subroutine subgrid_get_topounitinfo (ti, gi, nlunits, ncols, npfts, ncohorts, &
       nveg, &
       ncrop, &
       nurban_tbd, &
       nurban_hd, &
       nurban_md, &
       nlake, &
       nwetland, &
       nglacier, &
       nglacier_mec,  &
       glcmask)
    !
    ! !DESCRIPTION:
    ! Obtain topounit properties
    !
    ! !USES
    use elm_varpar  , only : natpft_size, cft_size, maxpatch_urb, maxpatch_glcmec
    use elm_varctl  , only : create_crop_landunit
    use clm_varsur  , only : wt_lunit, urban_valid, wt_glc_mec
    use landunit_varcon  , only : istsoil, istcrop, istice, istice_mec, istdlak, istwet, &
                             isturb_tbd, isturb_hd, isturb_md
    use FatesInterfaceTypesMod, only : fates_maxElementsPerSite

    !
    ! !ARGUMENTS
    implicit none
    integer , intent(in)  :: ti                   ! topounit index
    integer , intent(in)  :: gi  ! gridcell index for this topounit. Needed in development while
                                 ! each topounit is being given the same properties as the original gridcell.
    integer , optional, intent(out) :: nlunits    ! number of landunits
    integer , optional, intent(out) :: ncols      ! number of columns 
    integer , optional, intent(out) :: npfts      ! number of pfts 
    integer , optional, intent(out) :: ncohorts   ! number of cohorts 
    integer , optional, intent(out) :: nveg       ! number of vegetated pfts in naturally vegetated landunit
    integer , optional, intent(out) :: ncrop      ! number of crop pfts in crop landunit
    integer , optional, intent(out) :: nurban_tbd ! number of urban pfts (columns) in urban TBD landunit
    integer , optional, intent(out) :: nurban_hd ! number of urban pfts (columns) in urban HD landunit
    integer , optional, intent(out) :: nurban_md ! number of urban pfts (columns) in urban MD landunit
    integer , optional, intent(out) :: nlake      ! number of lake pfts (columns) in lake landunit
    integer , optional, intent(out) :: nwetland   ! number of wetland pfts (columns) in wetland landunit
    integer , optional, intent(out) :: nglacier   ! number of glacier pfts (columns) in glacier landunit
    integer , optional, intent(out) :: nglacier_mec  ! number of glacier_mec pfts (columns) in glacier_mec landunit
    integer , optional, intent(in)  :: glcmask  ! = 1 if glc requires surface mass balance in this gridcell
    !
    ! !LOCAL VARIABLES:
    integer  :: t, m             ! loop index
    integer  :: n                ! elevation class index
    integer  :: ilunits          ! number of landunits in topounit
    integer  :: icols            ! number of columns in topounit
    integer  :: ipfts            ! number of pfts in topounit
    integer  :: icohorts         ! number of cohorts in topounit
    integer  :: npfts_per_lunit  ! number of pfts in landunit
    !------------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Initialize landunits, pfts, columns and cohorts counters for topounit
    ! -------------------------------------------------------------------------

    ilunits    = 0
    icols      = 0
    ipfts      = 0
    icohorts   = 0
    
    ! -------------------------------------------------------------------------
    ! Initialize other optional counters for gridcell
    ! -------------------------------------------------------------------------
    if (present(nveg)) nveg = 0
    if (present(nurban_tbd)) nurban_tbd = 0
    if (present(nurban_hd)) nurban_hd = 0
    if (present(nurban_md)) nurban_md = 0
    if (present(nlake)) nlake = 0
    if (present(nwetland)) nwetland = 0
    if (present(nglacier)) nglacier = 0
    if (present(nglacier_mec)) nglacier_mec = 0
    if (present(ncrop)) ncrop = 0
    
    ! -------------------------------------------------------------------------
    ! Set naturally vegetated landunit
    ! -------------------------------------------------------------------------

    ! To support dynamic landunits, we have a naturally vegetated landunit in every grid
    ! cell, because it might need to come into existence even if its weight is 0 at the
    ! start of the run. And to support transient pfts or dynamic vegetation, we always
    ! allocate space for ALL PFTs on this landunit.

    npfts_per_lunit = natpft_size

    ! Assume that the vegetated landunit has one column
    ilunits = ilunits + 1
    icols = icols + 1  

    ipfts = ipfts + npfts_per_lunit

    ! -------------------------------------------------------------------------
    ! Number of cohorts is set here
    ! ED cohorts (via FATES) populate all natural vegetation columns.
    ! Current implementations mostly assume that only one column contains
    ! natural vegetation, which is synonymous with the soil column. 
    ! Because of the single natural vegetation paradigm currently in place,
    ! fates cohort allocations are based off a multiplier on grid-cells.
    ! -------------------------------------------------------------------------

    icohorts = fates_maxElementsPerSite

    if (present(nveg)) nveg = nveg + npfts_per_lunit

    ! -------------------------------------------------------------------------
    ! Set urban landunits
    ! -------------------------------------------------------------------------

    ! To support dynamic landunits, we have all urban landunits in every grid cell that
    ! has valid urban parameters, because they might need to come into existence even if
    ! their weight is 0 at the start of the run. And for simplicity, we always allocate
    ! space for ALL columns on the urban landunits.

    ! Set urban tall building district landunit

    npfts_per_lunit = 0
    if (urban_valid(gi)) then
       npfts_per_lunit = maxpatch_urb
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
       ipfts = ipfts + npfts_per_lunit
    end if
    if (present(nurban_tbd)) nurban_tbd = nurban_tbd + npfts_per_lunit

    ! Set urban high density landunit

    npfts_per_lunit = 0
    if (urban_valid(gi)) then
       npfts_per_lunit = maxpatch_urb
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
       ipfts = ipfts + npfts_per_lunit
    end if
    if (present(nurban_hd)) nurban_hd = nurban_hd + npfts_per_lunit

    ! Set urban medium density landunit

    npfts_per_lunit = 0
    if (urban_valid(gi)) then
       npfts_per_lunit = maxpatch_urb
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
       ipfts = ipfts + npfts_per_lunit
    end if
    if (present(nurban_md)) nurban_md = nurban_md + npfts_per_lunit

    ! -------------------------------------------------------------------------
    ! Set lake landunit
    ! -------------------------------------------------------------------------

    ! We currently do NOT allow the lake landunit to expand via dynamic landunits, so we
    ! only need to allocate space for it where its weight is currently non-zero.

    npfts_per_lunit = 0
    if (wt_lunit(gi, istdlak) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nlake)) nlake = nlake + npfts_per_lunit

    ! -------------------------------------------------------------------------
    ! Set wetland landunit
    ! -------------------------------------------------------------------------

    ! We currently do NOT allow the wetland landunit to expand via dynamic landunits, so
    ! we only need to allocate space for it where its weight is currently non-zero.

    npfts_per_lunit = 0
    if (wt_lunit(gi, istwet) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nwetland)) nwetland = nwetland + npfts_per_lunit

    ! -------------------------------------------------------------------------
    ! Set glacier landunit
    ! -------------------------------------------------------------------------

    ! We currently do NOT allow the glacier landunit to expand via dynamic landunits, so
    ! we only need to allocate space for it where its weight is currently non-zero. (If we
    ! have dynamic glacier area, we will be using glacier_mec landunits rather than
    ! glacier landunits.)

    npfts_per_lunit = 0
    if (wt_lunit(gi, istice) > 0.0_r8) then
       npfts_per_lunit = npfts_per_lunit + 1
    end if
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nglacier)) nglacier = nglacier + npfts_per_lunit

    ! -------------------------------------------------------------------------
    ! Set glacier_mec landunit
    ! -------------------------------------------------------------------------

    ! If glcmask = 1, we create a column for each elevation class even if the weight on
    ! the grid cell is 0. This is needed for coupling to CISM. In addition, this is
    ! currently sufficient to ensure that we have glaciers everywhere they might be
    ! needed with dynamic landunits, since CISM won't be able to create glaciers outside
    ! of the area specified by glcmask. 

    npfts_per_lunit = 0
    do m = 1, maxpatch_glcmec
       ! If the landunit has non-zero weight on the grid cell, and this column has
       ! non-zero weight on the landunit...
       if (wt_lunit(gi, istice_mec) > 0.0_r8 .and. wt_glc_mec(gi, m) > 0.0_r8) then
          npfts_per_lunit = npfts_per_lunit + 1

       elseif (present(glcmask)) then
          if (glcmask == 1) then      ! create a virtual column 
             npfts_per_lunit = npfts_per_lunit + 1
          endif  ! glcmask = 1 
       endif  ! wt > 0
    enddo   ! maxpatch_glcmec
    if (npfts_per_lunit > 0) then
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if (present(nglacier_mec)) nglacier_mec = nglacier_mec + npfts_per_lunit

    ! -------------------------------------------------------------------------
    ! Set crop landunit if appropriate
    ! -------------------------------------------------------------------------

    npfts_per_lunit = 0
    if (create_crop_landunit) then
       ! To support dynamic landunits, we have a crop landunit in every grid cell (if
       ! create_crop_landunit is true), because it might need to come into existence even
       ! if its weight is 0 at the start of the run.
       npfts_per_lunit = cft_size
       ilunits = ilunits + 1
       icols   = icols + npfts_per_lunit
       ipfts = ipfts + npfts_per_lunit
    end if
    if (present(ncrop)) ncrop = ncrop + npfts_per_lunit

    ! -------------------------------------------------------------------------
    ! Determine return arguments
    ! -------------------------------------------------------------------------

    if (present(nlunits))       nlunits = ilunits
    if (present(ncols))           ncols = icols
    if (present(npfts))           npfts = ipfts
    if (present(ncohorts))     ncohorts = icohorts

  end subroutine subgrid_get_topounitinfo
  
end module subgridMod
