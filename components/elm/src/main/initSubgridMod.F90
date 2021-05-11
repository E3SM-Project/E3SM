module initSubgridMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Lower-level routines for initializing the subgrid structure. This module is shared
  ! between both the production code (via initGridCellsMod) and unit testing code.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc
  use abortutils     , only : endrun
  use elm_varctl     , only : iulog
  use elm_varcon     , only : namep, namec, namel, namet
  use decompMod      , only : bounds_type
  use GridcellType   , only : grc_pp                
  Use TopounitType   , only : top_pp
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp                
  use VegetationType      , only : veg_pp                
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: elm_ptrs_compdown ! fill in data pointing down
  public :: elm_ptrs_check    ! checks and writes out a summary of subgrid data
  public :: add_topounit      ! add an entry in the topounit-level arrays
  public :: add_landunit      ! add an entry in the landunit-level arrays
  public :: add_column        ! add an entry in the column-level arrays
  public :: add_patch         ! add an entry in the patch-level arrays
  !
  !-----------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------------
  subroutine elm_ptrs_compdown(bounds)
    !
    ! !DESCRIPTION:
    ! Assumes the part of the subgrid pointing up has been set.  Fills 
    ! in the data pointing down.  Up is p_c, p_l, p_g, c_l, c_g, and l_g.
    !
    ! This algorithm assumes all indices besides grid cell are monotonically
    ! increasing.  (Note that grid cell index is NOT monotonically increasing,
    ! hence we cannot set initial & final indices at the grid cell level - 
    ! grc_pp%luni, grc_pp%lunf, etc.)
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
    use elm_varcon, only : ispval
    use topounit_varcon, only : max_topounits
    !
    ! !ARGUMENTS
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: t,l,c,p               ! loop counters
    integer :: curg,curt,curl,curc,curp ! tracks g,l,c,p indexes in arrays
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
       if (veg_pp%column(p) /= curc) then
          curc = veg_pp%column(p)
          if (curc < bounds%begc .or. curc > bounds%endc) then
             write(iulog,*) 'elm_ptrs_compdown ERROR: pcolumn ',p,curc,bounds%begc,bounds%endc
             call endrun(decomp_index=p, elmlevel=namep, msg=errMsg(__FILE__, __LINE__))
          endif
          col_pp%pfti(curc) = p
       endif
       col_pp%pftf(curc) = p
       col_pp%npfts(curc) = col_pp%pftf(curc) - col_pp%pfti(curc) + 1
       if (veg_pp%landunit(p) /= curl) then
          curl = veg_pp%landunit(p)
          if (curl < bounds%begl .or. curl > bounds%endl) then
             write(iulog,*) 'elm_ptrs_compdown ERROR: plandunit ',p,curl,bounds%begl,bounds%endl
             call endrun(decomp_index=p, elmlevel=namep, msg=errMsg(__FILE__, __LINE__))
          endif
          lun_pp%pfti(curl) = p
       endif
       lun_pp%pftf(curl) = p
       lun_pp%npfts(curl) = lun_pp%pftf(curl) - lun_pp%pfti(curl) + 1
    enddo

    curl = 0
    do c = bounds%begc,bounds%endc
       if (col_pp%landunit(c) /= curl) then
          curl = col_pp%landunit(c)
          if (curl < bounds%begl .or. curl > bounds%endl) then
             write(iulog,*) 'elm_ptrs_compdown ERROR: clandunit ',c,curl,bounds%begl,bounds%endl
             call endrun(decomp_index=c, elmlevel=namec, msg=errMsg(__FILE__, __LINE__))
          endif
          lun_pp%coli(curl) = c
       endif
       lun_pp%colf(curl) = c
       lun_pp%ncolumns(curl) = lun_pp%colf(curl) - lun_pp%coli(curl) + 1
    enddo
    
    ! Gridcell down pointers to topounits are monotonic, so those can be done like the 
    ! previous monotonic down pointers
    curg = 0
    do t = bounds%begt,bounds%endt
       if (top_pp%gridcell(t) /= curg) then
          curg = top_pp%gridcell(t)
          if (curg < bounds%begg .or. curg > bounds%endg) then
             write(iulog,*) 'elm_ptrs_compdown ERROR: tgridcell ',t,curg,bounds%begg,bounds%endg
             call endrun(decomp_index=t, elmlevel=namet, msg=errMsg(__FILE__, __LINE__))
          endif
          grc_pp%topi(curg) = t
       endif
       grc_pp%topf(curg) = t
       grc_pp%ntopounits(curg) = grc_pp%topf(curg) - grc_pp%topi(curg) + 1
    enddo

    ! Determine landunit_indices: indices into landunit-level arrays for each grid cell.
    ! Note that landunits not present in a given grid cell are set to ispval.
    ! Preliminary implementation of topounits: leave this unchanged, but will only work 
    ! for max_topounits = 1
    grc_pp%landunit_indices(:,bounds%begg:bounds%endg) = ispval
    do l = bounds%begl,bounds%endl
       ltype = lun_pp%itype(l)
       curg = lun_pp%gridcell(l)
       if (curg < bounds%begg .or. curg > bounds%endg) then
          write(iulog,*) 'elm_ptrs_compdown ERROR: gridcell landunit_indices ', l,curg,bounds%begg,bounds%endg
          call endrun(decomp_index=l, elmlevel=namel, msg=errMsg(__FILE__, __LINE__))
       end if

       if (grc_pp%landunit_indices(ltype, curg) == ispval) then
          grc_pp%landunit_indices(ltype, curg) = l
       else
          if (max_topounits == 1) then
            write(iulog,*) 'elm_ptrs_compdown ERROR: This landunit type has already been set for this gridcell'
            write(iulog,*) 'l, ltype, curg = ', l, ltype, curg
            call endrun(decomp_index=l, elmlevel=namel, msg=errMsg(__FILE__, __LINE__))
          end if
       end if
    end do

    ! Determine landunit_indices: indices into landunit-level arrays for each topounit.
    ! Note that landunits not present in a given topounit are set to ispval.
    top_pp%landunit_indices(:,bounds%begt:bounds%endt) = ispval
    do l = bounds%begl,bounds%endl
       ltype = lun_pp%itype(l)
       curt = lun_pp%topounit(l)
       if (curt < bounds%begt .or. curg > bounds%endt) then
          write(iulog,*) 'elm_ptrs_compdown ERROR: topounit landunit_indices ', l,curt,bounds%begt,bounds%endt
          call endrun(decomp_index=l, elmlevel=namel, msg=errMsg(__FILE__, __LINE__))
       end if

       if (top_pp%landunit_indices(ltype, curt) == ispval) then
          top_pp%landunit_indices(ltype, curt) = l
       else
          write(iulog,*) 'elm_ptrs_compdown ERROR: This landunit type has already been set for this topounit'
          write(iulog,*) 'l, ltype, curt = ', l, ltype, curt
          call endrun(decomp_index=l, elmlevel=namel, msg=errMsg(__FILE__, __LINE__))
       end if
    end do
    
  end subroutine elm_ptrs_compdown

  !------------------------------------------------------------------------------
  subroutine elm_ptrs_check(bounds)
    !
    ! !DESCRIPTION:
    ! Checks and writes out a summary of subgrid data
    !
    ! !USES
    use elm_varcon, only : ispval
    use landunit_varcon, only : max_lunit
    !
    ! !ARGUMENTS
    implicit none
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g,t,l,c,p     ! loop counters
    integer :: l_prev        ! l value of previous point
    integer :: ltype         ! landunit type
    logical :: error         ! error flag
    !------------------------------------------------------------------------------

    associate( &
         begg => bounds%begg, &
         endg => bounds%endg, &
         begt => bounds%begt, &
         endt => bounds%endt, &
         begl => bounds%begl, &
         endl => bounds%endl, &
         begc => bounds%begc, &
         endc => bounds%endc, &
         begp => bounds%begp, &
         endp => bounds%endp  &
         )
    
    if (masterproc) write(iulog,*) ' '
    if (masterproc) write(iulog,*) '---elm_ptrs_check:'

    !--- check index ranges ---
    error = .false.
    do g = begg, endg
       do ltype = 1, max_lunit
          l = grc_pp%landunit_indices(ltype, g)
          if (l /= ispval) then
             if (l < begl .or. l > endl) error = .true.
          end if
       end do
    end do
    if (error) then
       write(iulog,*) '   elm_ptrs_check: g index ranges - ERROR'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    if (masterproc) write(iulog,*) '   elm_ptrs_check: g index ranges - OK'

    error = .false.
    if (minval(lun_pp%gridcell(begl:endl)) < begg .or. maxval(lun_pp%gridcell(begl:endl)) > endg) error=.true.
    if (minval(lun_pp%topounit(begl:endl)) < begt .or. maxval(lun_pp%topounit(begl:endl)) > endt) error=.true.
    if (minval(lun_pp%coli(begl:endl)) < begc .or. maxval(lun_pp%coli(begl:endl)) > endc) error=.true.
    if (minval(lun_pp%colf(begl:endl)) < begc .or. maxval(lun_pp%colf(begl:endl)) > endc) error=.true.
    if (minval(lun_pp%pfti(begl:endl)) < begp .or. maxval(lun_pp%pfti(begl:endl)) > endp) error=.true.
    if (minval(lun_pp%pftf(begl:endl)) < begp .or. maxval(lun_pp%pftf(begl:endl)) > endp) error=.true.
    if (error) then
       write(iulog,*) '   elm_ptrs_check: l index ranges - ERROR'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    if (masterproc) write(iulog,*) '   elm_ptrs_check: l index ranges - OK'

    error = .false.
    if (minval(col_pp%gridcell(begc:endc)) < begg .or. maxval(col_pp%gridcell(begc:endc)) > endg) error=.true.
    if (minval(col_pp%topounit(begc:endc)) < begt .or. maxval(col_pp%topounit(begc:endc)) > endt) error=.true.
    if (minval(col_pp%landunit(begc:endc)) < begl .or. maxval(col_pp%landunit(begc:endc)) > endl) error=.true.
    if (minval(col_pp%pfti(begc:endc)) < begp .or. maxval(col_pp%pfti(begc:endc)) > endp) error=.true.
    if (minval(col_pp%pftf(begc:endc)) < begp .or. maxval(col_pp%pftf(begc:endc)) > endp) error=.true.
    if (error) then
       write(iulog,*) '   elm_ptrs_check: c index ranges - ERROR'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    if (masterproc) write(iulog,*) '   elm_ptrs_check: c index ranges - OK'

    error = .false.
    if (minval(veg_pp%gridcell(begp:endp)) < begg .or. maxval(veg_pp%gridcell(begp:endp)) > endg) error=.true.
    if (minval(veg_pp%topounit(begp:endp)) < begt .or. maxval(veg_pp%topounit(begp:endp)) > endt) error=.true.
    if (minval(veg_pp%landunit(begp:endp)) < begl .or. maxval(veg_pp%landunit(begp:endp)) > endl) error=.true.
    if (minval(veg_pp%column(begp:endp)) < begc .or. maxval(veg_pp%column(begp:endp)) > endc) error=.true.
    if (error) then
       write(iulog,*) '   elm_ptrs_check: p index ranges - ERROR'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    if (masterproc) write(iulog,*) '   elm_ptrs_check: p index ranges - OK'

    !--- check that indices in arrays are monotonically increasing ---
    error = .false.
    do l=begl+1,endl
      if ((lun_pp%itype(l) == lun_pp%itype(l-1)) .and. &
           lun_pp%gridcell(l) < lun_pp%gridcell(l-1)) then
         ! grid cell indices should be monotonically increasing for a given landunit type
         error = .true.
      end if
      if (lun_pp%coli(l) < lun_pp%coli(l-1)) error = .true.
      if (lun_pp%colf(l) < lun_pp%colf(l-1)) error = .true.
      if (lun_pp%pfti(l) < lun_pp%pfti(l-1)) error = .true.
      if (lun_pp%pftf(l) < lun_pp%pftf(l-1)) error = .true.
      if (error) then
         write(iulog,*) '   elm_ptrs_check: l mono increasing - ERROR'
         call endrun(decomp_index=l, elmlevel=namel, msg=errMsg(__FILE__, __LINE__))
      endif
    enddo
    if (masterproc) write(iulog,*) '   elm_ptrs_check: l mono increasing - OK'

    error = .false.
    do c=begc+1,endc
      l = col_pp%landunit(c)
      l_prev = col_pp%landunit(c-1)
      if ((lun_pp%itype(l) == lun_pp%itype(l_prev)) .and. &
           col_pp%gridcell(c) < col_pp%gridcell(c-1)) then
         ! grid cell indices should be monotonically increasing for a given landunit type
         error = .true.
      end if
      if (col_pp%landunit(c) < col_pp%landunit(c-1)) error = .true.
      if (col_pp%pfti(c) < col_pp%pfti(c-1)) error = .true.
      if (col_pp%pftf(c) < col_pp%pftf(c-1)) error = .true.
      if (error) then
         write(iulog,*) '   elm_ptrs_check: c mono increasing - ERROR'
         call endrun(decomp_index=c, elmlevel=namec, msg=errMsg(__FILE__, __LINE__))
      endif
    enddo
    if (masterproc) write(iulog,*) '   elm_ptrs_check: c mono increasing - OK'

    error = .false.
    do p=begp+1,endp
      l = veg_pp%landunit(p)
      l_prev = veg_pp%landunit(p-1)
      if ((lun_pp%itype(l) == lun_pp%itype(l_prev)) .and. &
           veg_pp%gridcell(p) < veg_pp%gridcell(p-1)) then
         ! grid cell indices should be monotonically increasing for a given landunit type
         error = .true.
      end if
      if (veg_pp%landunit(p) < veg_pp%landunit(p-1)) error = .true.
      if (veg_pp%column  (p) < veg_pp%column  (p-1)) error = .true.
      if (error) then
         write(iulog,*) '   elm_ptrs_check: p mono increasing - ERROR'
         call endrun(decomp_index=p, elmlevel=namep, msg=errMsg(__FILE__, __LINE__))
      endif
    enddo
    if (masterproc) write(iulog,*) '   elm_ptrs_check: p mono increasing - OK'

    !--- check that the tree is internally consistent ---
    error = .false.
    do g = begg, endg
       do ltype = 1, max_lunit
          l = grc_pp%landunit_indices(ltype, g)

          ! skip l == ispval, which implies that this landunit type doesn't exist on this grid cell
          if (l /= ispval) then
             if (lun_pp%itype(l) /= ltype) error = .true.
             if (lun_pp%gridcell(l) /= g) error = .true.
             if (error) then
                write(iulog,*) '   elm_ptrs_check: tree consistent - ERROR'
                call endrun(decomp_index=l, elmlevel=namel, msg=errMsg(__FILE__, __LINE__))
             endif
             do c = lun_pp%coli(l),lun_pp%colf(l)
                if (col_pp%gridcell(c) /= g) error = .true.
                if (col_pp%landunit(c) /= l) error = .true.
                if (error) then
                   write(iulog,*) '   elm_ptrs_check: tree consistent - ERROR'
                   call endrun(decomp_index=c, elmlevel=namec, msg=errMsg(__FILE__, __LINE__))
                endif
                do p = col_pp%pfti(c),col_pp%pftf(c)
                   if (veg_pp%gridcell(p) /= g) error = .true.
                   if (veg_pp%landunit(p) /= l) error = .true.
                   if (veg_pp%column(p)   /= c) error = .true.
                   if (error) then
                      write(iulog,*) '   elm_ptrs_check: tree consistent - ERROR'
                      call endrun(decomp_index=p, elmlevel=namep, msg=errMsg(__FILE__, __LINE__))
                   endif
                enddo  ! p
             enddo  ! c
          end if  ! l /= ispval
       enddo  ! ltype
    enddo  ! g
    if (masterproc) write(iulog,*) '   elm_ptrs_check: tree consistent - OK'
    if (masterproc) write(iulog,*) ' '

    end associate
    
  end subroutine elm_ptrs_check

  !-----------------------------------------------------------------------
  subroutine add_topounit(ti, gi, wtgcell)
    !
    ! !DESCRIPTION:
    ! Add an entry in the topounit-level arrays. ti gives the index of the last topounit
    ! added; the new topounit is added at ti+1, and the ti argument is incremented
    ! accordingly.
    !
    ! !ARGUMENTS:
    integer  , intent(inout) :: ti      ! input value is index of last topounit added; output value is index of this newly-added topounit
    integer  , intent(in)    :: gi      ! gridcell index on which this topounit should be placed 
    real(r8) , intent(in)    :: wtgcell ! weight of the topounit relative to the gridcell
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'add_topounit'
    !-----------------------------------------------------------------------

    ti = ti + 1

    top_pp%gridcell(ti) = gi
    top_pp%wtgcell(ti) = wtgcell
    
  end subroutine add_topounit

  !-----------------------------------------------------------------------
  subroutine add_landunit(li, ti, ltype, wttopounit)
    !
    ! !DESCRIPTION:
    ! Add an entry in the landunit-level arrays. li gives the index of the last landunit
    ! added; the new landunit is added at li+1, and the li argument is incremented
    ! accordingly.
    !
    ! !USES:
    use landunit_varcon , only : istsoil, istcrop, istice_mec, istdlak, isturb_MIN, isturb_MAX
    !
    ! !ARGUMENTS:
    integer  , intent(inout) :: li         ! input value is index of last landunit added; output value is index of this newly-added landunit
    integer  , intent(in)    :: ti         ! topounit index on which this landunit should be placed
    integer  , intent(in)    :: ltype      ! landunit type
    real(r8) , intent(in)    :: wttopounit ! weight of the landunit relative to the topounit
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'add_landunit'
    !-----------------------------------------------------------------------
    
    li = li + 1

    lun_pp%topounit(li) = ti
    lun_pp%gridcell(li) = top_pp%gridcell(ti)
    
    lun_pp%wttopounit(li) = wttopounit
    lun_pp%itype(li) = ltype
    
    if (ltype == istsoil .or. ltype == istcrop) then
       lun_pp%ifspecial(li) = .false.
    else
       lun_pp%ifspecial(li) = .true.
    end if

    if (ltype == istice_mec) then
       lun_pp%glcmecpoi(li) = .true.
    else
       lun_pp%glcmecpoi(li) = .false.
    end if

    if (ltype == istdlak) then
       lun_pp%lakpoi(li) = .true.
    else
       lun_pp%lakpoi(li) = .false.
    end if

    if (ltype >= isturb_MIN .and. ltype <= isturb_MAX) then
       lun_pp%urbpoi(li) = .true.
    else
       lun_pp%urbpoi(li) = .false.
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

    col_pp%landunit(ci) = li
    col_pp%topounit(ci) = lun_pp%topounit(li)
    col_pp%gridcell(ci) = lun_pp%gridcell(li)
    
    col_pp%wtlunit(ci) = wtlunit
    col_pp%itype(ci) = ctype
    
  end subroutine add_column

  !-----------------------------------------------------------------------
  subroutine add_patch(pi, ci, ptype, wtcol)
    !
    ! !DESCRIPTION:
    ! Add an entry in the patch-level arrays. pi gives the index of the last patch added; the
    ! new patch is added at pi+1, and the pi argument is incremented accordingly.
    !
    ! !USES:
    use elm_varcon      , only : ispval
    use landunit_varcon , only : istsoil, istcrop
    use elm_varpar      , only : natpft_lb
    !
    ! !ARGUMENTS:
    integer  , intent(inout) :: pi    ! input value is index of last patch added; output value is index of this newly-added patch
    integer  , intent(in)    :: ci    ! column index on which this patch should be placed (assumes this column has already been created)
    integer  , intent(in)    :: ptype ! patch type
    real(r8) , intent(in)    :: wtcol ! weight of the patch relative to the column
    !
    ! !LOCAL VARIABLES:
    integer :: li  ! landunit index, for convenience
    integer :: lb_offset ! offset between natpft_lb and 1
    
    character(len=*), parameter :: subname = 'add_patch'
    !-----------------------------------------------------------------------
    
    pi = pi + 1

    veg_pp%column(pi) = ci
    veg_pp%landunit(pi) = col_pp%landunit(ci)
    veg_pp%topounit(pi) = col_pp%topounit(ci)
    veg_pp%gridcell(pi) = col_pp%gridcell(ci)
    
    veg_pp%wtcol(pi) = wtcol
    veg_pp%itype(pi) = ptype

    li = veg_pp%landunit(pi)
    if (lun_pp%itype(li) == istsoil .or. lun_pp%itype(li) == istcrop) then
       lb_offset = 1 - natpft_lb
       veg_pp%mxy(pi) = ptype + lb_offset
    else
       veg_pp%mxy(pi) = ispval
    end if

  end subroutine add_patch


end module initSubgridMod
