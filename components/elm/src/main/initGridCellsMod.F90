module initGridCellsMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initializes sub-grid mapping for each land grid cell. This module handles the high-
  ! level logic that determines how the subgrid structure is set up in an ELM run. 
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc,iam
  use abortutils     , only : endrun
  use elm_varctl     , only : iulog
  use elm_varctl     , only : use_fates, use_fates_sp, use_polygonal_tundra
  use elm_varcon     , only : namep, namec, namel, nameg
  use decompMod      , only : bounds_type, ldecomp
  use GridcellType   , only : grc_pp
  use TopounitType   , only : top_pp  
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp                
  use VegetationType , only : veg_pp                
  use initSubgridMod , only : elm_ptrs_compdown, elm_ptrs_check
  use initSubgridMod , only : add_topounit, add_landunit, add_polygon_landunit, add_column, add_patch
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public initGridcells ! initialize sub-grid gridcell mapping 
  public initGhostGridcells
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private set_landunit_veg_compete
  private set_landunit_wet_ice_lake
  private set_landunit_crop_noncompete
  private set_landunit_urban
  private set_topounit
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine initGridcells
    !
    ! !DESCRIPTION: 
    ! Initialize sub-grid mapping and allocates space for derived type hierarchy.
    ! For each land gridcell determine topounit, landunit, column and pft properties.
    !
    ! !USES
    use domainMod         , only : ldomain
    use decompMod         , only : get_proc_bounds, get_clump_bounds, get_proc_clumps
    use subgridWeightsMod , only : compute_higher_order_weights
    use topounit_varcon   , only : max_topounits, has_topounit 
    !!use elm_varsur       , only : wt_tunit, elv_tunit, slp_tunit,asp_tunit,num_tunit_per_grd
    use landunit_varcon   , only : istsoil, istice, istwet, istdlak, istice_mec
    use landunit_varcon   , only : isturb_tbd, isturb_hd, isturb_md, istcrop
    use elm_varctl        , only : create_glacier_mec_landunit
    use shr_const_mod     , only : SHR_CONST_PI
    !
    ! !LOCAL VARIABLES:
    integer :: nc,ti,li,ci,pi,gdc,topounit, topo_ind      ! indices
    integer :: nclumps                           ! number of clumps on this processor
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
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump, ti, li, ci, pi, gdc, topounit)
    do nc = 1, nclumps

       call get_clump_bounds(nc, bounds_clump)
       
       ! Initialize indexing counters for each subgrid level
       ti = bounds_clump%begt-1
       li = bounds_clump%begl-1
       ci = bounds_clump%begc-1
       pi = bounds_clump%begp-1
       
       ! For each gridcell in clump, create the correct number of topounits       
       do gdc = bounds_clump%begg, bounds_clump%endg
          call set_topounit(gdc, ti, ldomain%num_tunits_per_grd(gdc) )
       end do

       ! With all topounits defined, next place landunits

       ! Determine naturally vegetated landunit
       do topounit = bounds_clump%begt,bounds_clump%endt
          topo_ind = top_pp%topo_grc_ind(topounit)
          call set_landunit_veg_compete(               &
               ltype=istsoil, gi=top_pp%gridcell(topounit), ti=topounit,topo_ind=topo_ind, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       ! Determine crop landunit
       do topounit = bounds_clump%begt,bounds_clump%endt
          topo_ind = top_pp%topo_grc_ind(topounit)
          call set_landunit_crop_noncompete(           &
               ltype=istcrop, gi=top_pp%gridcell(topounit), ti=topounit,topo_ind=topo_ind, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       ! Determine urban tall building district landunit
       do topounit = bounds_clump%begt,bounds_clump%endt
          topo_ind = top_pp%topo_grc_ind(topounit)
          call set_landunit_urban( &
               ltype=isturb_tbd, gi=top_pp%gridcell(topounit), ti=topounit,topo_ind=topo_ind, li=li, ci=ci, pi=pi, &
               setdata=.true.)

       end do

       ! Determine urban high density landunit
       do topounit = bounds_clump%begt,bounds_clump%endt
          topo_ind = top_pp%topo_grc_ind(topounit)
          call set_landunit_urban( &
               ltype=isturb_hd, gi=top_pp%gridcell(topounit), ti=topounit,topo_ind=topo_ind, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       ! Determine urban medium density landunit
       do topounit = bounds_clump%begt,bounds_clump%endt
          topo_ind = top_pp%topo_grc_ind(topounit)
          call set_landunit_urban( &
               ltype=isturb_md, gi=top_pp%gridcell(topounit), ti=topounit,topo_ind=topo_ind, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       ! Determine lake, wetland and glacier landunits 
       do topounit = bounds_clump%begt,bounds_clump%endt
          topo_ind = top_pp%topo_grc_ind(topounit)
          call set_landunit_wet_ice_lake(              &
               ltype=istdlak, gi=top_pp%gridcell(topounit), ti=topounit,topo_ind=topo_ind, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       do topounit = bounds_clump%begt,bounds_clump%endt
          topo_ind = top_pp%topo_grc_ind(topounit)
          call set_landunit_wet_ice_lake(              &
               ltype=istwet, gi=top_pp%gridcell(topounit), ti=topounit,topo_ind=topo_ind, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       do topounit = bounds_clump%begt,bounds_clump%endt
          topo_ind = top_pp%topo_grc_ind(topounit)
          call set_landunit_wet_ice_lake(              &
               ltype=istice, gi=top_pp%gridcell(topounit), ti=topounit,topo_ind=topo_ind, li=li, ci=ci, pi=pi, &
               setdata=.true.)
       end do

       if (create_glacier_mec_landunit) then
          do topounit = bounds_clump%begt,bounds_clump%endt
             topo_ind = top_pp%topo_grc_ind(topounit)
             gdc = top_pp%gridcell(topounit)
             call set_landunit_wet_ice_lake(              &
                  ltype=istice_mec, gi=gdc, ti=topounit,topo_ind=topo_ind, li=li, ci=ci, pi=pi, &
                  setdata=.true., &
                  glcmask = ldomain%glcmask(gdc))
          end do
       endif

       ! Ensure that we have set the expected number of pfts, cols and landunits for this clump
       SHR_ASSERT(ti == bounds_clump%endt, errMsg(__FILE__, __LINE__))
       SHR_ASSERT(li == bounds_clump%endl, errMsg(__FILE__, __LINE__))
       SHR_ASSERT(ci == bounds_clump%endc, errMsg(__FILE__, __LINE__))
       SHR_ASSERT(pi == bounds_clump%endp, errMsg(__FILE__, __LINE__))

       ! Set some other gridcell-level variables

       do gdc = bounds_clump%begg,bounds_clump%endg
          grc_pp%gindex(gdc) = ldecomp%gdc2glo(gdc)
          grc_pp%area(gdc)   = ldomain%area(gdc)
          grc_pp%latdeg(gdc) = ldomain%latc(gdc) 
          grc_pp%londeg(gdc) = ldomain%lonc(gdc) 
          grc_pp%lat(gdc)    = grc_pp%latdeg(gdc) * SHR_CONST_PI/180._r8  
          grc_pp%lon(gdc)    = grc_pp%londeg(gdc) * SHR_CONST_PI/180._r8

          grc_pp%stdev_elev(gdc)     = ldomain%stdev_elev(gdc)
          grc_pp%sky_view(gdc)       = ldomain%sky_view(gdc)
          grc_pp%terrain_config(gdc) = ldomain%terrain_config(gdc)
          grc_pp%sinsl_cosas(gdc)    = ldomain%sinsl_cosas(gdc)
          grc_pp%sinsl_sinas(gdc)    = ldomain%sinsl_sinas(gdc)
          
       enddo

       ! Fill in subgrid datatypes

       call elm_ptrs_compdown(bounds_clump)

       ! By putting this check within the loop over clumps, we ensure that (for example)
       ! if a clump is responsible for landunit L, then that same clump is also
       ! responsible for all columns and pfts in L.
       !$OMP CRITICAL
       call elm_ptrs_check(bounds_clump)
       !$OMP END CRITICAL

       ! Set veg_pp%wtlunit, veg_pp%wtgcell and col_pp%wtgcell
       call compute_higher_order_weights(bounds_clump)

    end do
   !$OMP END PARALLEL DO 
  end subroutine initGridcells

  !----------------------------------------------------------------------
  subroutine set_topounit(gdc, ti, num_tunits_per_grd ) 
    !
    ! !DESCRIPTION:
    ! Initialize each topounit for a gridcell.
    !
    ! !USES
    use elm_varsur , only : wt_tunit, elv_tunit, slp_tunit, asp_tunit, num_tunit_per_grd 
    use elm_varctl , only : use_IM2_hillslope_hydrology
    use topounit_varcon   , only : max_topounits, has_topounit 
    ! !ARGUMENTS
    integer, intent(in) :: gdc
    integer, intent(inout) :: ti
    integer, intent(in) :: num_tunits_per_grd
    ! !LOCAL VARIABLES
    integer :: topounit, ntopos,topo_ind, num_topo_tmp,tmp_tpu
    real(r8) :: wttopounit2gridcell, elv, slp                  ! topounit weight on gridcell, elevation and slope
    integer :: asp                                             ! aspect
    integer :: t1, t2, begt, endt, dn_index, min_index         ! local topounit indexing
    real(r8):: t1_elev, t2_elev, min_elev, dn_elev             ! for finding downhill neighbor
    logical :: is_tpu_active                                   ! Check if topounit is active
     
    tmp_tpu = num_tunits_per_grd       ! Actual number of topounits per grid
    if(max_topounits > 1) then
       ntopos = tmp_tpu                                
    else 
       ntopos = max_topounits
    endif

    begt = ti+1
    endt = begt + ntopos - 1
    
    do topounit = 1, ntopos                    ! use actual/valid # of topounits per grid intead of max_topounits
       if (max_topounits == 1) then
           wttopounit2gridcell = 1.0           ! The weight of topounit is 1 if only 1 topounit per grid
           is_tpu_active = .true.              ! Make topounit active if only one topounit is in a grid
       else
           wttopounit2gridcell = wt_tunit(gdc,topounit) !grc_pp%tfrc_area(gdc,topounit) 
           !if (topounit <= num_topo_tmp) then
           if (wttopounit2gridcell > 0.0) then
               is_tpu_active = .true.
           else
               is_tpu_active = .false.
           endif                    
       endif
       elv = elv_tunit(gdc,topounit) !grc_pp%televation(gdc,topounit) 
       slp = slp_tunit(gdc,topounit) !grc_pp%tslope(gdc,topounit) 
       asp = asp_tunit(gdc,topounit) !grc_pp%taspect(gdc,topounit) 
       topo_ind = topounit
       call add_topounit(ti=ti, gi=gdc, wtgcell=wttopounit2gridcell, elv=elv, slp=slp, asp=asp,topo_ind=topo_ind,is_tpu_active = is_tpu_active)
    end do

    ! Loop through topounits again to find its nearest downhill topounit on this gridcell
    ! part of the IM2 hillslope hydrology implementation
    if (ntopos > 1 .and. use_IM2_hillslope_hydrology) then
      ! find the minimum elevation over all topounits on the gridcell
      min_elev = top_pp%elevation(begt)
      min_index = begt
      do t1 = begt+1, endt
         if (top_pp%elevation(t1) < min_elev) then
            min_elev = top_pp%elevation(t1)
            min_index = t1
         endif
      end do
      ! find the closest downhill neighbor for each topounit
      dn_index = -1  ! value of -1 indicates no downhill neighbor
      do t1 = begt, endt
         t1_elev = top_pp%elevation(t1)
         dn_elev = min_elev
         do t2 = begt, endt
            t2_elev = top_pp%elevation(t2)
            if ((t2_elev < t1_elev) .and. (t2_elev >= dn_elev)) then
               dn_elev = t2_elev
               dn_index = t2
            endif
         end do
         top_pp%downhill_ti(t1)=dn_index
      end do
   endif
 
  end subroutine set_topounit  
 
 
  !------------------------------------------------------------------------
  subroutine set_landunit_veg_compete (ltype, gi, ti,topo_ind, li, ci, pi, setdata)
    !
    ! !DESCRIPTION: 
    ! Initialize vegetated landunit with competition
    !
    ! !USES
    use elm_varsur, only : wt_lunit, wt_nat_patch, wt_polygon
    use subgridMod, only : subgrid_get_topounitinfo
    use elm_varpar, only : numpft, maxpatch_pft, numcft, natpft_lb, natpft_ub, natpft_size
    use landunit_varcon, only: istlowcenpoly, istflatcenpoly, isthighcenpoly, max_non_poly_lunit
    !
    ! !ARGUMENTS:    
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(in)    :: ti                ! topounit index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! patch index
    integer , intent(inout) :: topo_ind                ! topounit index within each grid
    logical , intent(in)    :: setdata           ! set info or just compute
    !
    ! !LOCAL VARIABLES:
    integer  :: m,tgi,z                          ! index
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: pitype                           ! patch itype
    real(r8) :: wtlunit2topounit                 ! landunit weight on topounit
    real(r8) :: p_wt                             ! patch weight (0-1)
    real(r8) :: wtpoly2lndunit                   ! weight of polygon type wrt nat. veg. landunit
    !------------------------------------------------------------------------

    ! Set decomposition properties

    ! Initial topounit implementation: use the pft information provided at the gridcell
    ! level to assign PFTs on veg landunit for each topounit. Also, use the existing landunit weights on the 
    ! gridcell as the new landunit weights on each topounit.
    ! Later, this information will come from new surface datasat.
    call subgrid_get_topounitinfo(ti, gi,tgi=topo_ind, nveg=npfts)
    wtlunit2topounit = wt_lunit(gi,topo_ind, ltype)

    ! for polygonal tundra, we need to adjust the weight of the landunit as
    ! this wtlunit2topounit now corresponds to 4 landunits:
    ! a) natural vegetation, no polygonal tundra
    ! b) natural vegetation, high centered polygons
    ! c) natural vegetation, flat centered polygons
    ! d) natural vegetation, low centered polygons

    ! For FATES: the total number of patches may not match what is in the surface
    ! file, and therefor the weighting can't be used. The weightings in
    ! wt_nat_patch may be meaningful (like with fixed biogeography), but they
    ! they need a mapping table to connect to the allocated patches (in fates)
    ! so the wt_nat_patch array is not applicable to these area weights
    ! A subsequent call, via the clmfates interface will update these weights
    ! by using said mapping table
    
    if (npfts > 0) then
      ! do standard veg landunit first
       call add_landunit(li=li, ti=ti, ltype=ltype, wttopounit=wtlunit2topounit)
       
       ! Assume one column on the landunit
       call add_column(ci=ci, li=li, ctype=1, wtlunit=1.0_r8, is_soil=.true.)
       do m = natpft_lb,natpft_ub
          if(use_fates .and. .not.use_fates_sp)then
             p_wt = 1.0_r8/real(natpft_size,r8)
          else
             p_wt = wt_nat_patch(gi,topo_ind,m)
          end if
          call add_patch(pi=pi, ci=ci, ptype=m, wtcol=p_wt, is_on_soil_col=.true.)
       end do

       ! add polygonal landunits and columns if feature turned on
       ! continue to assume one column per landunit.
       if (use_polygonal_tundra) then
         ! loop over polygon types:
         do z = istlowcenpoly,isthighcenpoly
           ! get new weight for wttopounit:
           wtpoly2lndunit = wt_lunit(gi, topo_ind, z) 
           call add_polygon_landunit(li=li, ti=ti, ltype=ltype, wttopounit=wtpoly2lndunit, polytype = z - max_non_poly_lunit)
           call add_column(ci=ci, li=li, ctype=1, wtlunit=1.0_r8, is_soil=.true.)
           ! add patch:
           do m = natpft_lb,natpft_ub
             if(use_fates .and. .not.use_fates_sp) then
               p_wt = 1.0_r8/real(natpft_size,r8)
             else
               p_wt = wt_nat_patch(gi,topo_ind,m)
             end if
             call add_patch(pi=pi, ci=ci, ptype=m, wtcol=p_wt, is_on_soil_col=.true.)
           end do
         end do
       end if 

    end if

  end subroutine set_landunit_veg_compete
  
  !------------------------------------------------------------------------
  subroutine set_landunit_wet_ice_lake (ltype, gi, ti,topo_ind, li, ci, pi, setdata, glcmask)
    !
    ! !DESCRIPTION: 
    ! Initialize wet_ice_lake landunits that are non-urban (lake, wetland, glacier, glacier_mec)
    !
    ! !USES
    use elm_varpar      , only : maxpatch_glcmec
    use elm_varsur      , only : wt_lunit, wt_glc_mec
    use landunit_varcon , only : istwet, istdlak, istice, istice_mec
    use column_varcon   , only : icemec_class_to_col_itype
    use subgridMod      , only : subgrid_get_topounitinfo
    use pftvarcon       , only : noveg

    !
    ! !ARGUMENTS:    
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(in)    :: ti                ! topounit index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! patch index
    integer , intent(inout) :: topo_ind                ! topounit index within each grid
    logical , intent(in)    :: setdata           ! set info or just compute
    integer , intent(in), optional :: glcmask    ! = 1 where glc requires sfc mass balance
    !
    ! !LOCAL VARIABLES:
    integer  :: m,tgi                                ! index
    integer  :: c                                ! column loop index
    integer  :: ier                              ! error status 
    integer  :: npfts                            ! number of pfts in landunit
    real(r8) :: wtlunit2topounit                 ! landunit weight in topounit
    real(r8) :: wtcol2lunit                      ! col weight in landunit
    logical  :: is_lake_col
    !------------------------------------------------------------------------

    ! Set decomposition properties
    ! Initial topounit implementation: use the pft information provided at the gridcell
    ! level to assign PFTs on landunit for each topounit. Also, use the existing landunit weights on the 
    ! gridcell as the new landunit weights on each topounit.
    ! Later, this information will come from new surface datasat.

    is_lake_col = .false.
    if (ltype == istwet) then
       call subgrid_get_topounitinfo(ti, gi,tgi=topo_ind, nwetland=npfts)
    else if (ltype == istdlak) then
       is_lake_col = .true.
       call subgrid_get_topounitinfo(ti, gi,tgi=topo_ind, nlake=npfts)
    else if (ltype == istice) then 
       call subgrid_get_topounitinfo(ti, gi,tgi=topo_ind, nglacier=npfts)
    else if (ltype == istice_mec) then
       call subgrid_get_topounitinfo(ti, gi,tgi=topo_ind, nglacier_mec=npfts, glcmask = glcmask)
    else
       write(iulog,*)' set_landunit_wet_ice_lake: ltype of ',ltype,' not valid'
       write(iulog,*)' only istwet, istdlak, istice and istice_mec ltypes are valid'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    wtlunit2topounit = wt_lunit(gi,topo_ind, ltype)

    if (npfts > 0) then

       if (npfts /=1 .and. ltype /= istice_mec) then
          write(iulog,*)' set_landunit_wet_ice_lake: landunit must'// &
               ' have one pft '
          write(iulog,*)' current value of npfts=',npfts
          write(iulog,*)' landunit type = ',ltype
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       if (ltype==istice_mec) then   ! multiple columns per landunit

          call add_landunit(li=li, ti=ti, ltype=ltype, wttopounit=wtlunit2topounit)

          ! Determine column and properties
          ! (Each column has its own pft)
          ! 
          ! For grid cells with glcmask = 1, make sure all the elevations classes
          !  are populated, even if some have zero fractional area.  This ensures that the 
          !  ice sheet component, glc, will receive a surface mass balance in each elevation 
          !  class wherever the SMB is needed.
          ! Columns with zero weight are referred to as "virtual" columns.
 
          do m = 1, maxpatch_glcmec

             wtcol2lunit = wt_glc_mec(gi,topo_ind, m)

             if (wtcol2lunit > 0._r8 .or. glcmask == 1) then
                call add_column(ci=ci, li=li, ctype=icemec_class_to_col_itype(m), wtlunit=wtcol2lunit)
                call add_patch(pi=pi, ci=ci, ptype=noveg, wtcol=1.0_r8)
             endif
          enddo

       else

          ! Currently assume that each landunit only has only one column 
          ! and that each column has its own pft
       
          call add_landunit(li=li, ti=ti, ltype=ltype, wttopounit=wtlunit2topounit)
          call add_column(ci=ci, li=li, ctype=ltype, wtlunit=1.0_r8, is_lake=is_lake_col)
          call add_patch(pi=pi, ci=ci, ptype=noveg, wtcol=1.0_r8)

       end if   ! ltype = istice_mec
    endif       ! npfts > 0       

  end subroutine set_landunit_wet_ice_lake

  !------------------------------------------------------------------------

  subroutine set_landunit_crop_noncompete (ltype, gi, ti,topo_ind, li, ci, pi, setdata)
    !
    ! !DESCRIPTION: 
    ! Initialize crop landunit without competition
    !
    ! Note about the ltype input argument: This provides the value for this landunit index
    ! (i.e., the crop landunit index). This may differ from the landunit's 'itype' value,
    ! since itype is istsoil if we are running with create_crop_landunit but crop_prog = false.
    !
    ! !USES
    use elm_varsur      , only : wt_lunit, wt_cft
    use landunit_varcon , only : istcrop, istsoil
    use subgridMod      , only : subgrid_get_topounitinfo
    use elm_varctl      , only : create_crop_landunit
    use elm_varpar      , only : maxpatch_pft, numcft, crop_prog, cft_lb, cft_ub
    !
    ! !ARGUMENTS:
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(in)    :: ti                ! topounit index
    integer , intent(inout) :: topo_ind                ! topounit index within each grid
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! patch index
    logical , intent(in)    :: setdata           ! set info or just compute
    !
    ! !LOCAL VARIABLES:
    integer  :: my_ltype                         ! landunit type for crops
    integer  :: m,tgi,z                                ! index
    integer  :: npfts                            ! number of pfts in landunit
    real(r8) :: wtlunit2topounit                 ! landunit weight in topounit
    !------------------------------------------------------------------------

    ! Set decomposition properties

    ! Initial topounit implementation: use the pft information provided at the gridcell
    ! level to assign PFTs on landunit for each topounit. Also, use the existing landunit weights on the 
    ! gridcell as the new landunit weights on each topounit.
    ! Later, this information will come from new surface datasat.
    call subgrid_get_topounitinfo(ti, gi,tgi=topo_ind, ncrop=npfts)
    wtlunit2topounit = wt_lunit(gi,topo_ind, ltype)

    if (npfts > 0) then

       ! Note that we cannot simply use the 'ltype' argument to set itype here,
       ! because ltype will always indicate istcrop
       if (create_crop_landunit) then
          my_ltype = istcrop
       else
          my_ltype = istsoil
       end if

       call add_landunit(li=li, ti=ti, ltype=my_ltype, wttopounit=wtlunit2topounit)

       ! Set column and pft properties for this landunit 
       ! (each column has its own pft)

       if (create_crop_landunit) then
          do m = cft_lb, cft_ub
             call add_column(ci=ci, li=li, ctype=((istcrop*100) + m), wtlunit=wt_cft(gi,topo_ind,m), is_crop=.true.)
             call add_patch(pi=pi, ci=ci, ptype=m, wtcol=1.0_r8, is_on_crop_col=.true.)
          end do
       end if

    end if
       
  end subroutine set_landunit_crop_noncompete

  !------------------------------------------------------------------------------

  subroutine set_landunit_urban (ltype, gi, ti,topo_ind, li, ci, pi, setdata)
    !
    ! !DESCRIPTION: 
    ! Initialize urban landunits
    !
    ! !USES
    use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon   , only : icol_road_perv, icol_road_imperv
    use landunit_varcon , only : isturb_tbd, isturb_hd, isturb_md, isturb_MIN
    use elm_varpar      , only : maxpatch_urb
    use elm_varsur      , only : wt_lunit
    use subgridMod      , only : subgrid_get_topounitinfo
    use UrbanParamsType , only : urbinp
    use decompMod       , only : ldecomp
    use pftvarcon       , only : noveg
    !
    ! !ARGUMENTS:
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(in)    :: ti                ! topounit index
    integer , intent(inout) :: topo_ind                ! topounit index within each grid
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! patch index
    logical , intent(in)    :: setdata           ! set info or just compute
    !
    ! !LOCAL VARIABLES:
    integer  :: c                ! column loop index
    integer  :: m,tgi                ! index
    integer  :: n                ! urban density type index
    integer  :: ctype            ! column type
    integer  :: npfts            ! number of pfts in landunit
    real(r8) :: wtlunit2topounit ! weight relative to topounit of landunit
    real(r8) :: wtcol2lunit      ! weight of column with respect to landunit
    real(r8) :: wtlunit_roof     ! weight of roof with respect to landunit
    real(r8) :: wtroad_perv      ! weight of pervious road column with respect to total road
    integer  :: ier              ! error status 
    !------------------------------------------------------------------------

    ! Set decomposition properties, and set variables specific to urban density type
    ! Initial topounit implementation: use the pft information provided at the gridcell
    ! level to assign PFTs on landunit for each topounit. Also, use the existing landunit weights on the 
    ! gridcell as the new landunit weights on each topounit.
    ! Later, this information will come from new surface datasat.

    select case (ltype)
    case (isturb_tbd)
       call subgrid_get_topounitinfo(ti, gi,tgi=topo_ind, nurban_tbd=npfts)
    case (isturb_hd)
       call subgrid_get_topounitinfo(ti, gi,tgi=topo_ind, nurban_hd=npfts)
    case (isturb_md)
       call subgrid_get_topounitinfo(ti, gi,tgi=topo_ind, nurban_md=npfts)
    case default
       write(iulog,*)' set_landunit_urban: unknown ltype: ', ltype
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    wtlunit2topounit = wt_lunit(gi,topo_ind, ltype)

    n = ltype - isturb_MIN + 1
    wtlunit_roof = urbinp%wtlunit_roof(gi,topo_ind,n)
    wtroad_perv  = urbinp%wtroad_perv(gi,topo_ind,n)

    if (npfts > 0) then

       call add_landunit(li=li, ti=ti, ltype=ltype, wttopounit=wtlunit2topounit)

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

  !------------------------------------------------------------------------
  subroutine initGhostGridcells()
    !
    ! !DESCRIPTION:
    !   Initialize ghost/halo subgrid categroies.
    !

#ifdef MOAB_LATERAL
    call initGhostColumnsMOAB()
#endif
  end subroutine initGhostGridcells

#if MOAB_LATERAL
  !------------------------------------------------------------------------
  subroutine initGhostColumnsMOAB()
    !
    ! DESCRIPTION
    !
    use decompMod        , only : get_proc_bounds
    use domainLateralMod , only : ldomain_lateral, GridLevelIntegerDataHaloExchange
    !
    implicit none
    !
    ! LOCAL VARIABLES
    type(bounds_type)  :: bounds_proc ! temporary
    integer, pointer   :: num_nat_veg_columns(:)
    integer            :: g, c
    integer, parameter :: icol_nat_veg = 1

    call get_proc_bounds(bounds_proc)

    allocate(num_nat_veg_columns(bounds_proc%begg:bounds_proc%endg_all))
    num_nat_veg_columns(:) = 0

    ! for owned grid cells, save the number of naturally vegetated soil columns
    do c = bounds_proc%begc, bounds_proc%endc
       g = col_pp%gridcell(c)

       if (col_pp%itype(c) == icol_nat_veg) then
          if (num_nat_veg_columns(g) /= 0) then
             call endrun(msg='ERROR: initGhostColumnsMOAB only supports the one natural soil column per grid cell')
          else
             num_nat_veg_columns(g) = num_nat_veg_columns(g) + 1
          end if
       end if

    end do

    ! do the halo exchange
    call GridLevelIntegerDataHaloExchange(ldomain_lateral, bounds_proc%begg, bounds_proc%endg, bounds_proc%endg_all, num_nat_veg_columns)

    ! for ghost columns, set mapping to ghost grid cell
    c = bounds_proc%endc
    do g = bounds_proc%endg + 1, bounds_proc%endg_all
       if (num_nat_veg_columns(g) > 0) then
          c = c + 1
          col_pp%gridcell(c) = g
          col_pp%itype(c)    = icol_nat_veg
       endif
    end do

  end subroutine initGhostColumnsMOAB
#endif
end module initGridCellsMod
