module initGridCellsMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initializes sub-grid mapping for each land grid cell
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use spmdMod        , only : masterproc,iam,mpicom
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use decompMod      , only : bounds_type, ldecomp
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public initGridcells ! initialize sub-grid gridcell mapping 

  ! The following need to be public for setting up unit tests. They should not generally
  ! be called directly by other modules in the production code
  public add_landunit  ! add an entry in the landunit-level arrays
  public add_column    ! add an entry in the column-level arrays
  public add_patch     ! add an entry in the patch-level arrays
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private clm_ptrs_compdown
  private clm_ptrs_check
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
    use clmtype 
    use domainMod         , only : ldomain
    use decompMod         , only : get_proc_bounds, get_clump_bounds, get_proc_clumps
    use subgridWeightsMod , only : compute_higher_order_weights
    use clm_varcon        , only : istsoil, istice, istwet, istdlak, istice_mec, &
                                   isturb_tbd, isturb_hd, isturb_md, istcrop
    use clm_varctl        , only : create_glacier_mec_landunit
    use shr_const_mod     , only : SHR_CONST_PI
    !
    ! !ARGUMENTS:
    implicit none
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
    ! PFT type:      1   2   3   1   2   3   1   2   3   1   2   3
    !
    ! So note that clump index is most slowly varying, followed by landunit type,
    ! followed by gridcell, followed by column and pft type.


    nclumps = get_proc_clumps()

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
               ltype=istsoil, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
       end do

       ! Determine crop landunit
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_crop_noncompete(           &
               ltype=istcrop, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
       end do

       ! Determine urban tall building district landunit
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_urban( &
               ltype=isturb_tbd, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)

       end do

       ! Determine urban high density landunit
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_urban( &
               ltype=isturb_hd, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
       end do

       ! Determine urban medium density landunit
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_urban( &
               ltype=isturb_md, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
       end do

       ! Determine lake, wetland and glacier landunits 
       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_wet_ice_lake(              &
               ltype=istdlak, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
       end do

       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_wet_ice_lake(              &
               ltype=istwet, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
       end do

       do gdc = bounds_clump%begg,bounds_clump%endg
          call set_landunit_wet_ice_lake(              &
               ltype=istice, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true.)
       end do

       if (create_glacier_mec_landunit) then
          do gdc = bounds_clump%begg,bounds_clump%endg
             call set_landunit_wet_ice_lake(              &
                  ltype=istice_mec, gi=gdc, li=li, ci=ci, pi=pi, setdata=.true., &
                  glcmask = ldomain%glcmask(gdc))
          end do
       endif

       ! Ensure that we have set the expected number of pfts, cols and landunits for this clump
       SHR_ASSERT(li == bounds_clump%endl, errMsg(__FILE__, __LINE__))
       SHR_ASSERT(ci == bounds_clump%endc, errMsg(__FILE__, __LINE__))
       SHR_ASSERT(pi == bounds_clump%endp, errMsg(__FILE__, __LINE__))

       ! Set some other gridcell-level variables

       do gdc = bounds_clump%begg,bounds_clump%endg

          ! Make ice sheet masks

          grc%gris_mask(gdc) = 0._r8
          grc%gris_area(gdc) = 0._r8
          grc%aais_mask(gdc) = 0._r8
          grc%aais_area(gdc) = 0._r8

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

             grc%gris_mask(gdc) = 1.0_r8

          elseif (ldomain%latc(gdc) < -60.) then

             grc%aais_mask(gdc) = 1.0_r8

          endif  ! Greenland or Antarctic grid cell

          grc%gindex(gdc) = ldecomp%gdc2glo(gdc)
          grc%latdeg(gdc) = ldomain%latc(gdc) 
          grc%londeg(gdc) = ldomain%lonc(gdc) 
          grc%lat(gdc)    = grc%latdeg(gdc) * SHR_CONST_PI/180._r8  
          grc%lon(gdc)    = grc%londeg(gdc) * SHR_CONST_PI/180._r8
          grc%area(gdc)   = ldomain%area(gdc)

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

  !------------------------------------------------------------------------------
  subroutine clm_ptrs_compdown(bounds)
    !
    ! !DESCRIPTION:
    ! Assumes the part of the subgrid pointing up has been set.  Fills 
    ! in the data pointing down.  Up is p_c, p_l, p_g, c_l, c_g, and l_g.
    !
    ! This algorithm assumes all indices besides grid cell are monotonically
    ! increasing.  (Note that grid cell index is NOT monotonically increasing,
    ! hence we cannot set initial & final indices at the grid cell level - 
    ! grc%luni, grc%lunf, etc.)
    !
    ! Algorithm works as follows.  The p, c, and l loops march through
    ! the full arrays (nump, numc, and numl) checking the "up" indexes.
    ! As soon as the "up" index of the current (p,c,l) cell changes relative 
    ! to the previous (p,c,l) cell, the *i array will be set to point down 
    ! to that cell.  The *f array follows the same logic, so it's always the
    ! last "up" index from the previous cell when an "up" index changes.
    !
    ! For example, a case where p_c(1:4) = 1 and p_c(5:12) = 2.  This 
    ! subroutine will set c_pi(1) = 1, c_pf(1) = 4, c_pi(2) = 5, c_pf(2) = 12.
    !
    ! !USES
    use clmtype
    use clm_varcon, only : ispval
    !
    ! !ARGUMENTS
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: l,c,p               ! loop counters
    integer :: curg,curl,curc,curp ! tracks g,l,c,p indexes in arrays
    integer :: ltype               ! landunit type
    !------------------------------------------------------------------------------

    !--- Set the current c,l (curc, curl) to zero for initialization,
    !---   these indices track the current "up" index.
    !--- Take advantage of locality of l/c/p cells
    !--- Loop p through full local begp:endp length
    !--- Separately check the p_c, p_l, and p_g indexes for a change in
    !---   the "up" index.
    !--- If there is a change, verify that the current c,l,g is within the 
    !---   valid range, and set c_pi, l_pi, or g_pi to that current c,l,g
    !--- Constantly update the c_pf, l_pf, and g_pf array.  When the
    !---   g, l, c index changes, the *_pf array will be set correctly
    !--- Do the same for cols setting c_li, c_gi, c_lf, c_gf and
    !---   lunits setting l_gi, l_gf.

    curc = 0
    curl = 0
    do p = bounds%begp,bounds%endp
       if (pft%column(p) /= curc) then
          curc = pft%column(p)
          if (curc < bounds%begc .or. curc > bounds%endc) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: pcolumn ',p,curc,bounds%begc,bounds%endc
             call endrun(decomp_index=p, clmlevel=namep, msg=errMsg(__FILE__, __LINE__))
          endif
          col%pfti(curc) = p
       endif
       col%pftf(curc) = p
       col%npfts(curc) = col%pftf(curc) - col%pfti(curc) + 1
       if (pft%landunit(p) /= curl) then
          curl = pft%landunit(p)
          if (curl < bounds%begl .or. curl > bounds%endl) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: plandunit ',p,curl,bounds%begl,bounds%endl
             call endrun(decomp_index=p, clmlevel=namep, msg=errMsg(__FILE__, __LINE__))
          endif
          lun%pfti(curl) = p
       endif
       lun%pftf(curl) = p
       lun%npfts(curl) = lun%pftf(curl) - lun%pfti(curl) + 1
    enddo

    curl = 0
    do c = bounds%begc,bounds%endc
       if (col%landunit(c) /= curl) then
          curl = col%landunit(c)
          if (curl < bounds%begl .or. curl > bounds%endl) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: clandunit ',c,curl,bounds%begl,bounds%endl
             call endrun(decomp_index=c, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
          endif
          lun%coli(curl) = c
       endif
       lun%colf(curl) = c
       lun%ncolumns(curl) = lun%colf(curl) - lun%coli(curl) + 1
    enddo

    ! Determine landunit_indices: indices into landunit-level arrays for each grid cell.
    ! Note that landunits not present in a given grid cell are set to ispval.
    grc%landunit_indices(:,bounds%begg:bounds%endg) = ispval
    do l = bounds%begl,bounds%endl
       ltype = lun%itype(l)
       curg = lun%gridcell(l)
       if (curg < bounds%begg .or. curg > bounds%endg) then
          write(iulog,*) 'clm_ptrs_compdown ERROR: landunit_indices ', l,curg,bounds%begg,bounds%endg
          call endrun(decomp_index=l, clmlevel=namel, msg=errMsg(__FILE__, __LINE__))
       end if

       if (grc%landunit_indices(ltype, curg) == ispval) then
          grc%landunit_indices(ltype, curg) = l
       else
          write(iulog,*) 'clm_ptrs_compdown ERROR: This landunit type has already been set for this gridcell'
          write(iulog,*) 'l, ltype, curg = ', l, ltype, curg
          call endrun(decomp_index=l, clmlevel=namel, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine clm_ptrs_compdown

  !------------------------------------------------------------------------------
  subroutine clm_ptrs_check(bounds)
    !
    ! !DESCRIPTION:
    ! Checks and writes out a summary of subgrid data
    !
    ! !USES
    use clmtype
    use clm_varcon, only : ispval, max_lunit
    !
    ! !ARGUMENTS
    implicit none
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g,l,c,p       ! loop counters
    integer :: l_prev        ! l value of previous point
    integer :: ltype         ! landunit type
    logical :: error         ! error flag
    !------------------------------------------------------------------------------

    associate( &
         begg => bounds%begg, &
         endg => bounds%endg, &
         begl => bounds%begl, &
         endl => bounds%endl, &
         begc => bounds%begc, &
         endc => bounds%endc, &
         begp => bounds%begp, &
         endp => bounds%endp  &
         )
    
    if (masterproc) write(iulog,*) ' '
    if (masterproc) write(iulog,*) '---clm_ptrs_check:'

    !--- check index ranges ---
    error = .false.
    do g = begg, endg
       do ltype = 1, max_lunit
          l = grc%landunit_indices(ltype, g)
          if (l /= ispval) then
             if (l < begl .or. l > endl) error = .true.
          end if
       end do
    end do
    if (error) then
       write(iulog,*) '   clm_ptrs_check: g index ranges - ERROR'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    if (masterproc) write(iulog,*) '   clm_ptrs_check: g index ranges - OK'

    error = .false.
    if (minval(lun%gridcell(begl:endl)) < begg .or. maxval(lun%gridcell(begl:endl)) > endg) error=.true.
    if (minval(lun%coli(begl:endl)) < begc .or. maxval(lun%coli(begl:endl)) > endc) error=.true.
    if (minval(lun%colf(begl:endl)) < begc .or. maxval(lun%colf(begl:endl)) > endc) error=.true.
    if (minval(lun%pfti(begl:endl)) < begp .or. maxval(lun%pfti(begl:endl)) > endp) error=.true.
    if (minval(lun%pftf(begl:endl)) < begp .or. maxval(lun%pftf(begl:endl)) > endp) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: l index ranges - ERROR'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: l index ranges - OK'

    error = .false.
    if (minval(col%gridcell(begc:endc)) < begg .or. maxval(col%gridcell(begc:endc)) > endg) error=.true.
    if (minval(col%landunit(begc:endc)) < begl .or. maxval(col%landunit(begc:endc)) > endl) error=.true.
    if (minval(col%pfti(begc:endc)) < begp .or. maxval(col%pfti(begc:endc)) > endp) error=.true.
    if (minval(col%pftf(begc:endc)) < begp .or. maxval(col%pftf(begc:endc)) > endp) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: c index ranges - ERROR'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: c index ranges - OK'

    error = .false.
    if (minval(pft%gridcell(begp:endp)) < begg .or. maxval(pft%gridcell(begp:endp)) > endg) error=.true.
    if (minval(pft%landunit(begp:endp)) < begl .or. maxval(pft%landunit(begp:endp)) > endl) error=.true.
    if (minval(pft%column(begp:endp)) < begc .or. maxval(pft%column(begp:endp)) > endc) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: p index ranges - ERROR'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: p index ranges - OK'

    !--- check that indices in arrays are monotonically increasing ---
    error = .false.
    do l=begl+1,endl
      if ((lun%itype(l) == lun%itype(l-1)) .and. &
           lun%gridcell(l) < lun%gridcell(l-1)) then
         ! grid cell indices should be monotonically increasing for a given landunit type
         error = .true.
      end if
      if (lun%coli(l) < lun%coli(l-1)) error = .true.
      if (lun%colf(l) < lun%colf(l-1)) error = .true.
      if (lun%pfti(l) < lun%pfti(l-1)) error = .true.
      if (lun%pftf(l) < lun%pftf(l-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: l mono increasing - ERROR'
         call endrun(decomp_index=l, clmlevel=namel, msg=errMsg(__FILE__, __LINE__))
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: l mono increasing - OK'

    error = .false.
    do c=begc+1,endc
      l = col%landunit(c)
      l_prev = col%landunit(c-1)
      if ((lun%itype(l) == lun%itype(l_prev)) .and. &
           col%gridcell(c) < col%gridcell(c-1)) then
         ! grid cell indices should be monotonically increasing for a given landunit type
         error = .true.
      end if
      if (col%landunit(c) < col%landunit(c-1)) error = .true.
      if (col%pfti(c) < col%pfti(c-1)) error = .true.
      if (col%pftf(c) < col%pftf(c-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: c mono increasing - ERROR'
         call endrun(decomp_index=c, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: c mono increasing - OK'

    error = .false.
    do p=begp+1,endp
      l = pft%landunit(p)
      l_prev = pft%landunit(p-1)
      if ((lun%itype(l) == lun%itype(l_prev)) .and. &
           pft%gridcell(p) < pft%gridcell(p-1)) then
         ! grid cell indices should be monotonically increasing for a given landunit type
         error = .true.
      end if
      if (pft%landunit(p) < pft%landunit(p-1)) error = .true.
      if (pft%column  (p) < pft%column  (p-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: p mono increasing - ERROR'
         call endrun(decomp_index=p, clmlevel=namep, msg=errMsg(__FILE__, __LINE__))
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: p mono increasing - OK'

    !--- check that the tree is internally consistent ---
    error = .false.
    do g = begg, endg
       do ltype = 1, max_lunit
          l = grc%landunit_indices(ltype, g)

          ! skip l == ispval, which implies that this landunit type doesn't exist on this grid cell
          if (l /= ispval) then
             if (lun%itype(l) /= ltype) error = .true.
             if (lun%gridcell(l) /= g) error = .true.
             if (error) then
                write(iulog,*) '   clm_ptrs_check: tree consistent - ERROR'
                call endrun(decomp_index=l, clmlevel=namel, msg=errMsg(__FILE__, __LINE__))
             endif
             do c = lun%coli(l),lun%colf(l)
                if (col%gridcell(c) /= g) error = .true.
                if (col%landunit(c) /= l) error = .true.
                if (error) then
                   write(iulog,*) '   clm_ptrs_check: tree consistent - ERROR'
                   call endrun(decomp_index=c, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
                endif
                do p = col%pfti(c),col%pftf(c)
                   if (pft%gridcell(p) /= g) error = .true.
                   if (pft%landunit(p) /= l) error = .true.
                   if (pft%column(p)   /= c) error = .true.
                   if (error) then
                      write(iulog,*) '   clm_ptrs_check: tree consistent - ERROR'
                      call endrun(decomp_index=p, clmlevel=namep, msg=errMsg(__FILE__, __LINE__))
                   endif
                enddo  ! p
             enddo  ! c
          end if  ! l /= ispval
       enddo  ! ltype
    enddo  ! g
    if (masterproc) write(iulog,*) '   clm_ptrs_check: tree consistent - OK'
    if (masterproc) write(iulog,*) ' '

    end associate
    
  end subroutine clm_ptrs_check

  !------------------------------------------------------------------------
  subroutine set_landunit_veg_compete (ltype, gi, li, ci, pi, setdata)
    !
    ! !DESCRIPTION: 
    ! Initialize vegetated landunit with competition
    !
    ! !USES
    use clmtype
    use clm_varsur, only : wt_lunit, wt_nat_pft
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varpar, only : numpft, maxpatch_pft, numcft, natpft_lb, natpft_ub
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
    !
    ! !LOCAL VARIABLES:
    integer  :: m                                ! index
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: pitype                           ! pft itype
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
          call add_patch(pi=pi, ci=ci, ptype=m, wtcol=wt_nat_pft(gi,m))
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
    use clmtype
    use clm_varsur, only : wt_lunit, wt_glc_mec
    use clm_varcon, only : istwet, istdlak, istice, istice_mec, &
         icemec_class_to_col_itype
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varpar, only : maxpatch_glcmec
    use pftvarcon , only : noveg

    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
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
    use clmtype
    use clm_varsur, only : wt_lunit, wt_cft
    use clm_varcon, only : istcrop, istsoil
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varctl, only : create_crop_landunit
    use clm_varpar, only : maxpatch_pft, numcft, crop_prog, cft_lb, cft_ub
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
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
    use clm_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, &
                              icol_road_perv, icol_road_imperv, &
                              isturb_tbd, isturb_hd, isturb_md, isturb_MIN
    use clm_varpar   , only : maxpatch_urb
    use clmtype
    use clm_varsur   , only : wt_lunit
    use subgridMod   , only : subgrid_get_gcellinfo
    use UrbanInputMod, only : urbinp
    use decompMod    , only : ldecomp
    use pftvarcon    , only : noveg
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
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

  !-----------------------------------------------------------------------
  subroutine add_landunit(li, gi, ltype, wtgcell)
    !
    ! !DESCRIPTION:
    ! Add an entry in the landunit-level arrays. li gives the index of the last landunit
    ! added; the new landunit is added at li+1, and the li argument is incremented
    ! accordingly.
    !
    ! !USES:
    use clmtype    , only : lun
    use clm_varcon , only : istsoil, istcrop, istice_mec, istdlak, isturb_MIN, isturb_MAX
    !
    ! !ARGUMENTS:
    integer  , intent(inout) :: li      ! input value is index of last landunit added; output value is index of this newly-added landunit
    integer  , intent(in)    :: gi      ! grid cell index on which this landunit should be placed
    integer  , intent(in)    :: ltype   ! landunit type
    real(r8) , intent(in)    :: wtgcell ! weight of the landunit relative to the grid cell
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'add_landunit'
    !-----------------------------------------------------------------------
    
    li = li + 1

    lun%gridcell(li) = gi
    lun%wtgcell(li) = wtgcell
    lun%itype(li) = ltype
    
    if (ltype == istsoil .or. ltype == istcrop) then
       lun%ifspecial(li) = .false.
    else
       lun%ifspecial(li) = .true.
    end if

    if (ltype == istice_mec) then
       lun%glcmecpoi(li) = .true.
    else
       lun%glcmecpoi(li) = .false.
    end if

    if (ltype == istdlak) then
       lun%lakpoi(li) = .true.
    else
       lun%lakpoi(li) = .false.
    end if

    if (ltype >= isturb_MIN .and. ltype <= isturb_MAX) then
       lun%urbpoi(li) = .true.
    else
       lun%urbpoi(li) = .false.
    end if

  end subroutine add_landunit

  !-----------------------------------------------------------------------
  subroutine add_column(ci, li, ctype, wtlunit)
    !
    ! !DESCRIPTION:
    ! Add an entry in the column-level arrays. ci gives the index of the last column
    ! added; the new column is added at ci+1, and the ci argument is incremented
    ! accordingly.
    !
    ! !USES:
    use clmtype, only : col, lun
    !
    ! !ARGUMENTS:
    integer  , intent(inout) :: ci      ! input value is index of last column added; output value is index of this newly-added column
    integer  , intent(in)    :: li      ! landunit index on which this column should be placed (assumes this landunit has already been created)
    integer  , intent(in)    :: ctype   ! column type
    real(r8) , intent(in)    :: wtlunit ! weight of the column relative to the landunit
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'add_column'
    !-----------------------------------------------------------------------

    ci = ci + 1

    col%landunit(ci) = li
    col%gridcell(ci) = lun%gridcell(li)
    col%wtlunit(ci) = wtlunit
    col%itype(ci) = ctype
    
  end subroutine add_column

  !-----------------------------------------------------------------------
  subroutine add_patch(pi, ci, ptype, wtcol)
    !
    ! !DESCRIPTION:
    ! Add an entry in the patch-level arrays. pi gives the index of the last patch added; the
    ! new patch is added at pi+1, and the pi argument is incremented accordingly.
    !
    ! !USES:
    use clmtype    , only : pft, col, lun
    use clm_varcon , only : istsoil, istcrop, ispval
    use clm_varpar , only : natpft_lb
    !
    ! !ARGUMENTS:
    integer  , intent(inout) :: pi    ! input value is index of last patch added; output value is index of this newly-added patch
    integer  , intent(in)    :: ci    ! column index on which this patch should be placed (assumes this column has already been created)
    integer  , intent(in)    :: ptype ! patch type
    real(r8) , intent(in)    :: wtcol ! weight of the patch relative to the column
    !
    ! !LOCAL VARIABLES:
    integer :: li        ! landunit index
    integer :: lb_offset ! offset between natpft_lb and 1
    
    character(len=*), parameter :: subname = 'add_patch'
    !-----------------------------------------------------------------------
    
    pi = pi + 1

    pft%column(pi) = ci
    li = col%landunit(ci)
    pft%landunit(pi) = li
    pft%gridcell(pi) = col%gridcell(ci)

    pft%wtcol(pi) = wtcol

    pft%itype(pi) = ptype

    if (lun%itype(li) == istsoil .or. lun%itype(li) == istcrop) then
       lb_offset = 1 - natpft_lb
       pft%mxy(pi) = ptype + lb_offset
    else
       pft%mxy(pi) = ispval
    end if
    

  end subroutine add_patch

  

end module initGridCellsMod
