module initGridCellsMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initializes sub-grid mapping for each land grid cell. This module handles the high-
  ! level logic that determines how the subgrid structure is set up in a CLM run. It
  ! makes use of lower-level routines in initSubgridMod, which contains stuff that is
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc,iam
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use clm_varcon     , only : namep, namec, namel, nameg
  use decompMod      , only : bounds_type, ldecomp
  use GridcellType   , only : grc                
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : pft                
  use initSubgridMod , only : clm_ptrs_compdown, clm_ptrs_check
  use initSubgridMod , only : add_landunit, add_column, add_patch
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public initGridcells ! initialize sub-grid gridcell mapping 
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private set_cohort_decomp
  private set_landunit_veg_compete
  private set_landunit_wet_ice_lake
  private set_landunit_crop_noncompete
  private set_landunit_urban
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine initGridcells
    !
    ! !DESCRIPTION: 
    ! Initialize sub-grid mapping and allocates space for derived type hierarchy.
    ! For each land gridcell determine landunit, column and pft properties.
    !
    ! !USES
    use domainMod         , only : ldomain
    use decompMod         , only : get_proc_bounds, get_clump_bounds, get_proc_clumps
    use subgridWeightsMod , only : compute_higher_order_weights
    use landunit_varcon   , only : istsoil, istice, istwet, istdlak, istice_mec
    use landunit_varcon   , only : isturb_tbd, isturb_hd, isturb_md, istcrop
    use clm_varctl        , only : create_glacier_mec_landunit, use_ed
    use shr_const_mod     , only : SHR_CONST_PI
    !
    ! !LOCAL VARIABLES:
    integer :: nc,li,ci,pi,gdc      ! indices
    integer :: nclumps              ! number of clumps on this processor
    type(bounds_type) :: bounds_proc
    type(bounds_type) :: bounds_clump
    !------------------------------------------------------------------------

    ! Notes about how this routine is arranged, and its implications for the arrangement
    ! of 1-d vectors in memory: 
    ! 
    ! (1) There is an outer loop over clumps; this results in all of a clump's points (at
    !     the gridcell, landunit, column & pft level) being contiguous. This is important
    !     for the use of begg:endg, etc., and also for performance.
    !
    ! (2) Next, there is a section for each landunit, with the loop over grid cells
    !     happening separately for each landunit. This means that, within a given clump,
    !     points with the same landunit are grouped together (this is true at the
    !     landunit, column and pft levels). Thus, different landunits for a given grid
    !     cell are separated in memory. This improves performance in the many parts of
    !     the code that operate over a single landunit, or two similar landunits. 
    !
    ! Example: landunit-level array: For a processor with 2 clumps, each of which has 2
    ! grid cells, each of which has 3 landunits, the layout of a landunit-level array
    ! looks like the following:
    !
    ! Array index:   1   2   3   4   5   6   7   8   9  10  11  12
    ! ------------------------------------------------------------
    ! Clump index:   1   1   1   1   1   1   2   2   2   2   2   2
    ! Gridcell:      1   2   1   2   1   2   3   4   3   4   3   4
    ! Landunit type: 1   1   2   2   3   3   1   1   2   2   3   3
    !
    ! Example: pft-level array: For a processor with 1 clump, which has 2 grid cells, each
    ! of which has 2 landunits, each of which has 3 pfts, the layout of a pft-level array
    ! looks like the following:
    !
    ! Array index:   1   2   3   4   5   6   7   8   9  10  11  12
    ! ------------------------------------------------------------
    ! Gridcell:      1   1   1   2   2   2   1   1   1   2   2   2
    ! Landunit type: 1   1   1   1   1   1   2   2   2   2   2   2
    ! PATCH type:      1   2   3   1   2   3   1   2   3   1   2   3
    !
    ! So note that clump index is most slowly varying, followed by landunit type,
    ! followed by gridcell, followed by column and pft type.
    ! 
    ! Cohort layout
    ! Array index:   1   2   3   4   5   6   7   8   9  10  11  12
    ! ------------------------------------------------------------
    ! Gridcell:      1   1   2   2   3   3   1   1   2   2   3   3
    ! Cohort:        1   2   1   2   1   2   1   2   1   2   1   2

    nclumps = get_proc_clumps()

    ! FIX(SPM,032414) add private vars for cohort and perhaps patch dimension
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump, li, ci, pi, gdc)
    do nc = 1, nclumps

       call get_clump_bounds(nc, bounds_clump)

       ! For each land gridcell on global grid determine landunit, column and pft properties
       
       li = bounds_clump%begl-1
       ci = bounds_clump%begc-1
       pi = bounds_clump%begp-1

       ! Determine naturally vegetated landunit
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_veg_compete(               &
               ltype=istsoil, gi=gdc, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       ! Determine crop landunit
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_crop_noncompete(           &
               ltype=istcrop, gi=gdc, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       ! Determine urban tall building district landunit
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_urban( &
               ltype=isturb_tbd, gi=gdc, li=li, ci=ci, pi=pi, &
               setdata=.true.)

       end do

       ! Determine urban high density landunit
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_urban( &
               ltype=isturb_hd, gi=gdc, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       ! Determine urban medium density landunit
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_urban( &
               ltype=isturb_md, gi=gdc, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       ! Determine lake, wetland and glacier landunits 
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_wet_ice_lake(              &
               ltype=istdlak, gi=gdc, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_wet_ice_lake(              &
               ltype=istwet, gi=gdc, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_wet_ice_lake(              &
               ltype=istice, gi=gdc, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       if (create_glacier_mec_landunit) then
          do gdc = bounds_clump%begg,bounds_clump%endg
             call set_landunit_wet_ice_lake(              &
                  ltype=istice_mec, gi=gdc, li=li, ci=ci, pi=pi, &
                  setdata=.true., &
                  glcmask = ldomain%glcmask(gdc))
          end do
       endif

       if ( use_ed ) then
          ! cohort decomp
          call set_cohort_decomp( bounds_clump=bounds_clump )
       end if

       ! Ensure that we have set the expected number of pfts, cols and landunits for this clump
       SHR_ASSERT(li == bounds_clump%endl, errMsg(__FILE__, __LINE__))
       SHR_ASSERT(ci == bounds_clump%endc, errMsg(__FILE__, __LINE__))
       SHR_ASSERT(pi == bounds_clump%endp, errMsg(__FILE__, __LINE__))

       ! Set some other gridcell-level variables

       do gdc = bounds_clump%begg,bounds_clump%endg
          grc%gindex(gdc) = ldecomp%gdc2glo(gdc)
          grc%area(gdc)   = ldomain%area(gdc)
          grc%latdeg(gdc) = ldomain%latc(gdc) 
          grc%londeg(gdc) = ldomain%lonc(gdc) 
          grc%lat(gdc)    = grc%latdeg(gdc) * SHR_CONST_PI/180._r8  
          grc%lon(gdc)    = grc%londeg(gdc) * SHR_CONST_PI/180._r8
       enddo

       ! Fill in subgrid datatypes

       call clm_ptrs_compdown(bounds_clump)

       ! By putting this check within the loop over clumps, we ensure that (for example)
       ! if a clump is responsible for landunit L, then that same clump is also
       ! responsible for all columns and pfts in L.
       call clm_ptrs_check(bounds_clump)

       ! Set pft%wtlunit, pft%wtgcell and col%wtgcell
       call compute_higher_order_weights(bounds_clump)

    end do
    !$OMP END PARALLEL DO

  end subroutine initGridcells

  !------------------------------------------------------------------------
  subroutine set_cohort_decomp ( bounds_clump )
    !
    ! !DESCRIPTION: 
    ! Set gridcell decomposition for cohorts
    !
    use EDtypesMod      , only : cohorts_per_gcell
    use EDVecCohortType , only : coh
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds_clump  
    !
    ! !LOCAL VARIABLES:
    integer cohi, gi
    !------------------------------------------------------------------------

    gi = bounds_clump%begg

    do cohi = bounds_clump%begCohort, bounds_clump%endCohort

       coh%gridcell(cohi) = gi
       if ( mod(cohi,cohorts_per_gcell ) == 0 ) gi = gi + 1

     end do

  end subroutine set_cohort_decomp

  !------------------------------------------------------------------------
  subroutine set_landunit_veg_compete (ltype, gi, li, ci, pi, setdata)
    !
    ! !DESCRIPTION: 
    ! Initialize vegetated landunit with competition
    !
    ! !USES
    use clm_varsur, only : wt_lunit, wt_nat_patch
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varpar, only : numpft, maxpatch_pft, numcft, natpft_lb, natpft_ub
    !
    ! !ARGUMENTS:
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! patch index
    logical , intent(in)    :: setdata           ! set info or just compute
    !
    ! !LOCAL VARIABLES:
    integer  :: m                                ! index
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: pitype                           ! patch itype
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    !------------------------------------------------------------------------

    ! Set decomposition properties

    call subgrid_get_gcellinfo(gi, nveg=npfts)
    wtlunit2gcell = wt_lunit(gi, ltype)

    if (npfts > 0) then
       call add_landunit(li=li, gi=gi, ltype=ltype, wtgcell=wtlunit2gcell)
       
       ! Assume one column on the landunit
       call add_column(ci=ci, li=li, ctype=1, wtlunit=1.0_r8)

       do m = natpft_lb,natpft_ub
          call add_patch(pi=pi, ci=ci, ptype=m, wtcol=wt_nat_patch(gi,m))
       end do
    end if

  end subroutine set_landunit_veg_compete
  
  !------------------------------------------------------------------------
  subroutine set_landunit_wet_ice_lake (ltype, gi, li, ci, pi, setdata, glcmask)
    !
    ! !DESCRIPTION: 
    ! Initialize wet_ice_lake landunits that are non-urban (lake, wetland, glacier, glacier_mec)
    !
    ! !USES
    use clm_varpar      , only : maxpatch_glcmec
    use clm_varsur      , only : wt_lunit, wt_glc_mec
    use landunit_varcon , only : istwet, istdlak, istice, istice_mec
    use column_varcon   , only : icemec_class_to_col_itype
    use subgridMod      , only : subgrid_get_gcellinfo
    use pftvarcon       , only : noveg

    !
    ! !ARGUMENTS:
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! patch index
    logical , intent(in)    :: setdata           ! set info or just compute
    integer , intent(in), optional :: glcmask    ! = 1 where glc requires sfc mass balance
    !
    ! !LOCAL VARIABLES:
    integer  :: m                                ! index
    integer  :: c                                ! column loop index
    integer  :: ier                              ! error status 
    integer  :: npfts                            ! number of pfts in landunit
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    real(r8) :: wtcol2lunit                      ! col weight in landunit
    !------------------------------------------------------------------------

    ! Set decomposition properties

    if (ltype == istwet) then
       call subgrid_get_gcellinfo(gi, nwetland=npfts)
    else if (ltype == istdlak) then
       call subgrid_get_gcellinfo(gi, nlake=npfts)
    else if (ltype == istice) then 
       call subgrid_get_gcellinfo(gi, nglacier=npfts)
    else if (ltype == istice_mec) then
       call subgrid_get_gcellinfo(gi, nglacier_mec=npfts, glcmask = glcmask)
    else
       write(iulog,*)' set_landunit_wet_ice_lake: ltype of ',ltype,' not valid'
       write(iulog,*)' only istwet, istdlak, istice and istice_mec ltypes are valid'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    wtlunit2gcell = wt_lunit(gi, ltype)

    if (npfts > 0) then

       if (npfts /=1 .and. ltype /= istice_mec) then
          write(iulog,*)' set_landunit_wet_ice_lake: compete landunit must'// &
               ' have one pft '
          write(iulog,*)' current value of npfts=',npfts
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       if (ltype==istice_mec) then   ! multiple columns per landunit

          call add_landunit(li=li, gi=gi, ltype=ltype, wtgcell=wtlunit2gcell)

          ! Determine column and properties
          ! (Each column has its own pft)
          ! 
          ! For grid cells with glcmask = 1, make sure all the elevations classes
          !  are populated, even if some have zero fractional area.  This ensures that the 
          !  ice sheet component, glc, will receive a surface mass balance in each elevation 
          !  class wherever the SMB is needed.
          ! Columns with zero weight are referred to as "virtual" columns.
 
          do m = 1, maxpatch_glcmec

             wtcol2lunit = wt_glc_mec(gi,m)

             if (wtcol2lunit > 0._r8 .or. glcmask == 1) then
                call add_column(ci=ci, li=li, ctype=icemec_class_to_col_itype(m), wtlunit=wtcol2lunit)
                call add_patch(pi=pi, ci=ci, ptype=noveg, wtcol=1.0_r8)
             endif
          enddo

       else

          ! Currently assume that each landunit only has only one column 
          ! and that each column has its own pft
       
          call add_landunit(li=li, gi=gi, ltype=ltype, wtgcell=wtlunit2gcell)
          call add_column(ci=ci, li=li, ctype=ltype, wtlunit=1.0_r8)
          call add_patch(pi=pi, ci=ci, ptype=noveg, wtcol=1.0_r8)

       end if   ! ltype = istice_mec
    endif       ! npfts > 0       

  end subroutine set_landunit_wet_ice_lake

  !------------------------------------------------------------------------

  subroutine set_landunit_crop_noncompete (ltype, gi, li, ci, pi, setdata)
    !
    ! !DESCRIPTION: 
    ! Initialize crop landunit without competition
    !
    ! Note about the ltype input argument: This provides the value for this landunit index
    ! (i.e., the crop landunit index). This may differ from the landunit's 'itype' value,
    ! since itype is istsoil if we are running with create_crop_landunit but crop_prog = false.
    !
    ! !USES
    use clm_varsur      , only : wt_lunit, wt_cft
    use landunit_varcon , only : istcrop, istsoil
    use subgridMod      , only : subgrid_get_gcellinfo
    use clm_varctl      , only : create_crop_landunit
    use clm_varpar      , only : maxpatch_pft, numcft, crop_prog, cft_lb, cft_ub
    !
    ! !ARGUMENTS:
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! patch index
    logical , intent(in)    :: setdata           ! set info or just compute
    !
    ! !LOCAL VARIABLES:
    integer  :: my_ltype                         ! landunit type for crops
    integer  :: m                                ! index
    integer  :: npfts                            ! number of pfts in landunit
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    !------------------------------------------------------------------------

    ! Set decomposition properties

    call subgrid_get_gcellinfo(gi, ncrop=npfts)
    wtlunit2gcell = wt_lunit(gi, ltype)

    if (npfts > 0) then

       ! Note that we cannot simply use the 'ltype' argument to set itype here,
       ! because ltype will always indicate istcrop
       if ( crop_prog )then
          my_ltype = istcrop
       else
          my_ltype = istsoil
       end if

       call add_landunit(li=li, gi=gi, ltype=my_ltype, wtgcell=wtlunit2gcell)
       
       ! Set column and pft properties for this landunit 
       ! (each column has its own pft)

       if (create_crop_landunit) then
          do m = cft_lb, cft_ub
             call add_column(ci=ci, li=li, ctype=((istcrop*100) + m), wtlunit=wt_cft(gi,m))
             call add_patch(pi=pi, ci=ci, ptype=m, wtcol=1.0_r8)
          end do
       end if

    end if
       
  end subroutine set_landunit_crop_noncompete

  !------------------------------------------------------------------------------

  subroutine set_landunit_urban (ltype, gi, li, ci, pi, setdata)
    !
    ! !DESCRIPTION: 
    ! Initialize urban landunits
    !
    ! !USES
    use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon   , only : icol_road_perv, icol_road_imperv
    use landunit_varcon , only : isturb_tbd, isturb_hd, isturb_md, isturb_MIN
    use clm_varpar      , only : maxpatch_urb
    use clm_varsur      , only : wt_lunit
    use subgridMod      , only : subgrid_get_gcellinfo
    use UrbanParamsType , only : urbinp
    use decompMod       , only : ldecomp
    use pftvarcon       , only : noveg
    !
    ! !ARGUMENTS:
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! patch index
    logical , intent(in)    :: setdata           ! set info or just compute
    !
    ! !LOCAL VARIABLES:
    integer  :: c             ! column loop index
    integer  :: m             ! index
    integer  :: n             ! urban density type index
    integer  :: ctype         ! column type
    integer  :: npfts         ! number of pfts in landunit
    real(r8) :: wtlunit2gcell ! weight relative to gridcell of landunit
    real(r8) :: wtcol2lunit   ! weight of column with respect to landunit
    real(r8) :: wtlunit_roof  ! weight of roof with respect to landunit
    real(r8) :: wtroad_perv   ! weight of pervious road column with respect to total road
    integer  :: ier           ! error status 
    !------------------------------------------------------------------------

    ! Set decomposition properties, and set variables specific to urban density type

    select case (ltype)
    case (isturb_tbd)
       call subgrid_get_gcellinfo(gi, nurban_tbd=npfts)
    case (isturb_hd)
       call subgrid_get_gcellinfo(gi, nurban_hd=npfts)
    case (isturb_md)
       call subgrid_get_gcellinfo(gi, nurban_md=npfts)
    case default
       write(iulog,*)' set_landunit_urban: unknown ltype: ', ltype
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    wtlunit2gcell = wt_lunit(gi, ltype)

    n = ltype - isturb_MIN + 1
    wtlunit_roof = urbinp%wtlunit_roof(gi,n)
    wtroad_perv  = urbinp%wtroad_perv(gi,n)

    if (npfts > 0) then

       call add_landunit(li=li, gi=gi, ltype=ltype, wtgcell=wtlunit2gcell)

       ! Loop through columns for this landunit and set the column and pft properties
       ! For the urban landunits it is assumed that each column has its own pft
       
       do m = 1, maxpatch_urb
          
          if (m == 1) then
             ctype = icol_roof
             wtcol2lunit = wtlunit_roof
          else if (m == 2) then
             ctype = icol_sunwall
             wtcol2lunit = (1. - wtlunit_roof)/3
          else if (m == 3) then
             ctype = icol_shadewall
             wtcol2lunit = (1. - wtlunit_roof)/3
          else if (m == 4) then
             ctype = icol_road_imperv
             wtcol2lunit = ((1. - wtlunit_roof)/3) * (1.-wtroad_perv)
          else if (m == 5) then
             ctype = icol_road_perv
             wtcol2lunit = ((1. - wtlunit_roof)/3) * (wtroad_perv)
          end if

          call add_column(ci=ci, li=li, ctype=ctype, wtlunit=wtcol2lunit)

          call add_patch(pi=pi, ci=ci, ptype=noveg, wtcol=1.0_r8)

       end do   ! end of loop through urban columns-pfts
    end if

  end subroutine set_landunit_urban

end module initGridCellsMod
