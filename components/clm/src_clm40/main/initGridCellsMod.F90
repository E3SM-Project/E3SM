module initGridCellsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: initGridCellsMod
!
! !DESCRIPTION:
! Initializes sub-grid mapping for each land grid cell
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc,iam,mpicom
  use abortutils  , only : endrun
  use clm_varsur  , only : wtxy, vegxy
  use clm_varsur  , only : topoxy
  use clm_varctl  , only : iulog

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
  private clm_ptrs_compdown
  private clm_ptrs_check
  private set_landunit_veg_compete
  private set_landunit_wet_ice_lake
  private set_landunit_crop_noncompete
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !PRIVATE DATA MEMBERS: None
!EOP
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initGridcells
!
! !INTERFACE:
  subroutine initGridcells () 
!
! !DESCRIPTION: 
! Initialize sub-grid mapping and allocates space for derived type hierarchy.
! For each land gridcell determine landunit, column and pft properties.
!
! !USES
    use clmtype 
    use domainMod   , only : ldomain
    use decompMod   , only : ldecomp, get_proc_global, get_proc_bounds
    use clm_varcon  , only : istsoil, istice, istwet, istdlak, isturb, istice_mec
    use clm_varctl  , only : create_glacier_mec_landunit
    use clm_varcon  , only : istcrop
    use subgridMod  , only : subgrid_get_gcellinfo
    use shr_const_mod,only : SHR_CONST_PI
    use surfrdMod   , only : crop_prog
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: li,ci,pi,m,na,gdc,gsn,glo    ! indices
    integer :: nveg           ! number of pfts in naturally vegetated landunit
    integer :: ltype          ! landunit type
    real(r8):: wtveg          ! weight (gridcell) of naturally veg landunit
    integer :: ncrop          ! number of crop pfts in crop landunit
    real(r8):: wtcrop         ! weight (gridcell) of crop landunit
    integer :: nlake          ! number of pfts (columns) in lake landunit
    real(r8):: wtlake         ! weight (gridcell) of lake landunit
    integer :: nwetland       ! number of pfts (columns) in wetland landunit
    real(r8):: wtwetland      ! weight (gridcell) of wetland landunit
    integer :: nglacier       ! number of pfts (columns) in glacier landunit
    real(r8):: wtglacier      ! weight (gridcell) of glacier landunit
    integer :: nglacier_mec   ! number of pfts (columns) in glacier landunit
    real(r8):: wtglacier_mec  ! weight (gridcell) of glacier_mec landunit
    integer :: ier            ! error status
    integer :: numg           ! total number of gridcells across all processors
    integer :: numl           ! total number of landunits across all processors 
    integer :: numc           ! total number of columns across all processors 
    integer :: nump           ! total number of pfts across all processors 
    integer :: begg,endg      ! local beg/end gridcells gdc
    integer :: begl,endl      ! local beg/end landunits
    integer :: begc,endc      ! local beg/end columns 
    integer :: begp,endp      ! local beg/end pfts
    logical :: my_gcell       ! is gdc gridcell on my pe
    integer :: nwtxy          ! wtxy cell index

    type(gridcell_type), pointer  :: gptr ! pointer to gridcell derived subtype
    type(landunit_type), pointer  :: lptr ! pointer to landunit derived subtype
    type(column_type)  , pointer  :: cptr ! pointer to column derived subtype
    type(pft_type)     , pointer  :: pptr ! pointer to pft derived subtype
 !------------------------------------------------------------------------

    ! Set pointers into derived types for this module

    gptr => grc
    lptr => lun
    cptr => col
    pptr => pft

    ! Get total global number of grid cells, landunits, columns and pfts 
    
    call get_proc_global(numg,numl,numc,nump)
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! For each land gridcell on global grid determine landunit, column and pft properties

    li    = begl-1
    ci    = begc-1
    pi    = begp-1

    if ( crop_prog )then
       ltype = istcrop
    else
       ltype = istsoil
    end if

    !----- Set clm3 variables -----
    do gdc = begg,endg

       glo = ldecomp%gdc2glo(gdc)
       nwtxy = gdc

       my_gcell = .false.
       if (gdc >= begg .and. gdc <= endg) then
          my_gcell = .true.
       endif

       ! Determine naturally vegetated landunit

       call set_landunit_veg_compete(               &
            ltype=istsoil, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Determine crop landunit

       call set_landunit_crop_noncompete(           &
            ltype=ltype, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Determine urban landunit

       call set_landunit_urban( &
!           ltype=isturb, wtxy=wtxy, vegxy=vegxy,   &
            ltype=isturb, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       ! Determine lake, wetland and glacier landunits 

       call set_landunit_wet_ice_lake(              &
            ltype=istdlak, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       call set_landunit_wet_ice_lake(              &
            ltype=istwet, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       call set_landunit_wet_ice_lake(              &
            ltype=istice, &
            nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell)

       if (create_glacier_mec_landunit) then
          call set_landunit_wet_ice_lake(              &
               ltype=istice_mec, &
               nw=nwtxy, gi=gdc, li=li, ci=ci, pi=pi, setdata=my_gcell, &
               glcmask = ldomain%glcmask(gdc))
       endif

       ! Make ice sheet masks

       gptr%gris_mask(gdc) = 0._r8
       gptr%gris_area(gdc) = 0._r8
       gptr%aais_mask(gdc) = 0._r8
       gptr%aais_area(gdc) = 0._r8
      
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
 
            gptr%gris_mask(gdc) = 1.0_r8

      elseif (ldomain%latc(gdc) < -60.) then

            gptr%aais_mask(gdc) = 1.0_r8

       endif  ! Greenland or Antarctic grid cell

       ! Set clm3 lats/lons

       if (my_gcell) then
          gptr%gindex(gdc) = glo
          gptr%latdeg(gdc) = ldomain%latc(gdc) 
          gptr%londeg(gdc) = ldomain%lonc(gdc) 
          gptr%lat(gdc)    = gptr%latdeg(gdc) * SHR_CONST_PI/180._r8  
          gptr%lon(gdc)    = gptr%londeg(gdc) * SHR_CONST_PI/180._r8
          gptr%area(gdc)   = ldomain%area(gdc)
       endif

    enddo

    ! Fill in subgrid datatypes

    call clm_ptrs_compdown()
    call clm_ptrs_check()

  end subroutine initGridcells

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_ptrs_compdown
!
! !INTERFACE:
  subroutine clm_ptrs_compdown()
!
! !DESCRIPTION:
! Assumes the part of the subgrid pointing up has been set.  Fills 
! in the data pointing down.  Up is p_c, p_l, p_g, c_l, c_g, and l_g.
!
! This algorithm assumes all indices are monotonically increasing.
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
    use decompMod , only : get_proc_bounds

! !ARGUMENTS
    implicit none
!
! !CALLED FROM:
! subroutines initGridCellsMod
!
! !REVISION HISTORY:
! 2005.11.15  T Craig Creation
!
!
! !LOCAL VARIABLES:
    integer :: begg,endg,begl,endl,begc,endc,begp,endp ! beg/end glcp
    integer :: g,l,c,p               ! loop counters
    integer :: curg,curl,curc,curp   ! tracks g,l,c,p indexes in arrays
    type(gridcell_type), pointer  :: gptr ! pointer to gridcell derived subtype
    type(landunit_type), pointer  :: lptr ! pointer to landunit derived subtype
    type(column_type)  , pointer  :: cptr ! pointer to column derived subtype
    type(pft_type)     , pointer  :: pptr ! pointer to pft derived subtype
!EOP
!------------------------------------------------------------------------------

    gptr => grc
    lptr => lun
    cptr => col
    pptr => pft

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    !--- Set the current c,l,g (curc, curl, curg) to zero for initialization,
    !---   these indices track the current "up" index.
    !--- Take advantage of locality of g/l/c/p cells
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
    curg = 0
    do p = begp,endp
       if (pptr%column(p) /= curc) then
          curc = pptr%column(p)
          if (curc < begc .or. curc > endc) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: pcolumn ',p,curc,begc,endc
             call endrun()
          endif
          cptr%pfti(curc) = p
       endif
       cptr%pftf(curc) = p
       cptr%npfts(curc) = cptr%pftf(curc) - cptr%pfti(curc) + 1
       if (pptr%landunit(p) /= curl) then
          curl = pptr%landunit(p)
          if (curl < begl .or. curl > endl) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: plandunit ',p,curl,begl,endl
             call endrun()
          endif
          lptr%pfti(curl) = p
       endif
       lptr%pftf(curl) = p
       lptr%npfts(curl) = lptr%pftf(curl) - lptr%pfti(curl) + 1
       if (pptr%gridcell(p) /= curg) then
          curg = pptr%gridcell(p)
          if (curg < begg .or. curg > endg) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: pgridcell ',p,curg,begg,endg
             call endrun()
          endif
          gptr%pfti(curg) = p
       endif
       gptr%pftf(curg) = p
       gptr%npfts(curg) = gptr%pftf(curg) - gptr%pfti(curg) + 1
    enddo

    curg = 0
    curl = 0
    do c = begc,endc
       if (cptr%landunit(c) /= curl) then
          curl = cptr%landunit(c)
          if (curl < begl .or. curl > endl) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: clandunit ',c,curl,begl,endl
             call endrun()
          endif
          lptr%coli(curl) = c
       endif
       lptr%colf(curl) = c
       lptr%ncolumns(curl) = lptr%colf(curl) - lptr%coli(curl) + 1
       if (cptr%gridcell(c) /= curg) then
          curg = cptr%gridcell(c)
          if (curg < begg .or. curg > endg) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: cgridcell ',c,curg,begg,endg
             call endrun()
          endif
          gptr%coli(curg) = c
       endif
       gptr%colf(curg) = c
       gptr%ncolumns(curg) = gptr%colf(curg) - gptr%coli(curg) + 1
    enddo

    curg = 0
    do l = begl,endl
       if (lptr%gridcell(l) /= curg) then
          curg = lptr%gridcell(l)
          if (curg < begg .or. curg > endg) then
             write(iulog,*) 'clm_ptrs_compdown ERROR: lgridcell ',l,curg,begg,endg
             call endrun()
          endif
          gptr%luni(curg) = l
       endif
       gptr%lunf(curg) = l
       gptr%nlandunits(curg) = gptr%lunf(curg) - gptr%luni(curg) + 1
    enddo

    end subroutine clm_ptrs_compdown
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_ptrs_check
!
! !INTERFACE:
  subroutine clm_ptrs_check()
!
! !DESCRIPTION:
! Checks and writes out a summary of subgrid data
!
! !USES
    use clmtype
    use decompMod , only : get_proc_bounds

! !ARGUMENTS
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 2005.11.15  T Craig Creation
!
!
! !LOCAL VARIABLES:
    type(gridcell_type), pointer  :: gptr ! pointer to gridcell derived subtype
    type(landunit_type), pointer  :: lptr ! pointer to landunit derived subtype
    type(column_type)  , pointer  :: cptr ! pointer to column derived subtype
    type(pft_type)     , pointer  :: pptr ! pointer to pft derived subtype
    integer :: begg,endg,begl,endl,begc,endc,begp,endp   ! beg/end indices
    integer :: g,l,c,p       ! loop counters
    logical :: error         ! error flag
!EOP
!------------------------------------------------------------------------------

    gptr => grc
    lptr => lun
    cptr => col
    pptr => pft
    
    if (masterproc) write(iulog,*) ' '
    if (masterproc) write(iulog,*) '---clm_ptrs_check:'
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    !--- check index ranges ---
    error = .false.
    if (minval(gptr%luni) < begl .or. maxval(gptr%luni) > endl) error=.true.
    if (minval(gptr%lunf) < begl .or. maxval(gptr%lunf) > endl) error=.true.
    if (minval(gptr%coli) < begc .or. maxval(gptr%coli) > endc) error=.true.
    if (minval(gptr%colf) < begc .or. maxval(gptr%colf) > endc) error=.true.
    if (minval(gptr%pfti) < begp .or. maxval(gptr%pfti) > endp) error=.true.
    if (minval(gptr%pftf) < begp .or. maxval(gptr%pftf) > endp) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: g index ranges - ERROR'
       write(iulog,*)'minval,beg,maxval,end'
       write(iulog,*) minval(gptr%luni),begl,maxval(gptr%luni),endl
       write(iulog,*) minval(gptr%lunf),begl,maxval(gptr%lunf),endl
       write(iulog,*) minval(gptr%coli),begc,maxval(gptr%coli),endc
       write(iulog,*) minval(gptr%colf),begc,maxval(gptr%colf),endc
       write(iulog,*) minval(gptr%pfti),begp,maxval(gptr%pfti),endp
       write(iulog,*) minval(gptr%pftf),begp,maxval(gptr%pftf),endp
       call endrun()
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: g index ranges - OK'

    error = .false.
    if (minval(lptr%gridcell) < begg .or. maxval(lptr%gridcell) > endg) error=.true.
    if (minval(lptr%coli) < begc .or. maxval(lptr%coli) > endc) error=.true.
    if (minval(lptr%colf) < begc .or. maxval(lptr%colf) > endc) error=.true.
    if (minval(lptr%pfti) < begp .or. maxval(lptr%pfti) > endp) error=.true.
    if (minval(lptr%pftf) < begp .or. maxval(lptr%pftf) > endp) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: l index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: l index ranges - OK'

    error = .false.
    if (minval(cptr%gridcell) < begg .or. maxval(cptr%gridcell) > endg) error=.true.
    if (minval(cptr%landunit) < begl .or. maxval(cptr%landunit) > endl) error=.true.
    if (minval(cptr%pfti) < begp .or. maxval(cptr%pfti) > endp) error=.true.
    if (minval(cptr%pftf) < begp .or. maxval(cptr%pftf) > endp) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: c index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: c index ranges - OK'

    error = .false.
    if (minval(pptr%gridcell) < begg .or. maxval(pptr%gridcell) > endg) error=.true.
    if (minval(pptr%landunit) < begl .or. maxval(pptr%landunit) > endl) error=.true.
    if (minval(pptr%column) < begc .or. maxval(pptr%column) > endc) error=.true.
    if (error) then
       write(iulog,*) '   clm_ptrs_check: p index ranges - ERROR'
       call endrun()
    endif
    if (masterproc) write(iulog,*) '   clm_ptrs_check: p index ranges - OK'

    !--- check that indices in arrays are monotonically increasing ---
    error = .false.
    do g=begg+1,endg
      if (gptr%luni(g) < gptr%luni(g-1)) error = .true.
      if (gptr%lunf(g) < gptr%lunf(g-1)) error = .true.
      if (gptr%coli(g) < gptr%coli(g-1)) error = .true.
      if (gptr%colf(g) < gptr%colf(g-1)) error = .true.
      if (gptr%pfti(g) < gptr%pfti(g-1)) error = .true.
      if (gptr%pftf(g) < gptr%pftf(g-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: g mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: g mono increasing - OK'

    error = .false.
    do l=begl+1,endl
      if (lptr%gridcell(l) < lptr%gridcell(l-1)) error = .true.
      if (lptr%coli(l) < lptr%coli(l-1)) error = .true.
      if (lptr%colf(l) < lptr%colf(l-1)) error = .true.
      if (lptr%pfti(l) < lptr%pfti(l-1)) error = .true.
      if (lptr%pftf(l) < lptr%pftf(l-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: l mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: l mono increasing - OK'

    error = .false.
    do c=begc+1,endc
      if (cptr%gridcell(c) < cptr%gridcell(c-1)) error = .true.
      if (cptr%landunit(c) < cptr%landunit(c-1)) error = .true.
      if (cptr%pfti(c)     < cptr%pfti(c-1)) error = .true.
      if (cptr%pftf(c)     < cptr%pftf(c-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: c mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: c mono increasing - OK'

    error = .false.
    do p=begp+1,endp
      if (pptr%gridcell(p) < pptr%gridcell(p-1)) error = .true.
      if (pptr%landunit(p) < pptr%landunit(p-1)) error = .true.
      if (pptr%column  (p) < pptr%column  (p-1)) error = .true.
      if (error) then
         write(iulog,*) '   clm_ptrs_check: p mono increasing - ERROR'
         call endrun()
      endif
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: p mono increasing - OK'

    !--- check that the tree is internally consistent ---
    error = .false.
    do g = begg, endg
       do l = gptr%luni(g),gptr%lunf(g)
          if (lptr%gridcell(l) /= g) error = .true.
          do c = lptr%coli(l),lptr%colf(l)
             if (cptr%gridcell(c) /= g) error = .true.
             if (cptr%landunit(c) /= l) error = .true.
             do p = cptr%pfti(c),cptr%pftf(c)
                if (pptr%gridcell(p) /= g) error = .true.
                if (pptr%landunit(p) /= l) error = .true.
                if (pptr%column(p)   /= c) error = .true.
                if (error) then
                   write(iulog,*) '   clm_ptrs_check: tree consistent - ERROR'
                   call endrun()
                endif
             enddo
          enddo
       enddo
    enddo
    if (masterproc) write(iulog,*) '   clm_ptrs_check: tree consistent - OK'
    if (masterproc) write(iulog,*) ' '

end subroutine clm_ptrs_check
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_veg_compete
!
! !INTERFACE:
!  subroutine set_landunit_veg_compete (ltype, wtxy, vegxy, &
  subroutine set_landunit_veg_compete (ltype, &
                           nw, gi, li, ci, pi, setdata)
!
! !DESCRIPTION: 
! Initialize vegetated landunit with competition
!
! !USES
    use clmtype 
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varpar, only : numpft, maxpatch_pft, numcft
    use clm_varctl, only : allocate_all_vegpfts, create_crop_landunit
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
!    real(r8), intent(in)    :: wtxy(:,:)         ! subgrid patch weights
!    integer , intent(in)    :: vegxy(:,:)        ! PFT types 
    integer , intent(in)    :: nw                ! cell index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
!
! !REVISION HISTORY:
! Created by ?
! 2005.11.25 Updated by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: m                                ! m index in wtxy(nw,m)
    integer  :: n                                ! loop index
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    integer  :: pitype                           ! pft itype
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    type(landunit_type), pointer :: lptr         ! pointer to landunit
    type(column_type)  , pointer :: cptr         ! pointer to column
    type(pft_type)     , pointer :: pptr         ! pointer to pft

!------------------------------------------------------------------------

    ! Set decomposition properties

!    call subgrid_get_gcellinfo(nw, wtxy, nveg=npfts, wtveg=wtlunit2gcell)
    call subgrid_get_gcellinfo(nw, nveg=npfts, wtveg=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       lptr => lun
       cptr => col
       pptr => pft

       ncols = 1
       
       li = li + 1
       ci = ci + 1

       if (setdata) then
          ! Set landunit properties
          lptr%ifspecial(li) = .false.
          lptr%lakpoi(li)    = .false.
          lptr%urbpoi(li)    = .false.
          lptr%itype(li)     = ltype
       
          lptr%gridcell (li) = gi
          lptr%wtgcell(li) = wtlunit2gcell

          ! Set column properties for this landunit (only one column on landunit)
          cptr%itype(ci)    = 1
      
          cptr%gridcell (ci) = gi
          cptr%wtgcell(ci) = wtlunit2gcell
          cptr%landunit (ci) = li
          cptr%wtlunit(ci) = 1.0_r8
       endif ! setdata

       ! Set pft properties for this landunit

       if (create_crop_landunit) then
          do n = 1,numpft+1-numcft
             pi = pi + 1
             pitype = n-1
             if (setdata) then
                pptr%mxy(pi)      = n
                pptr%itype(pi)    = pitype
                pptr%gridcell(pi) = gi
                pptr%landunit(pi) = li
                pptr%column (pi) = ci
                pptr%wtgcell(pi) = 0.0_r8
                pptr%wtlunit(pi) = 0.0_r8
                pptr%wtcol(pi) = 0.0_r8
                do m = 1,maxpatch_pft
                   if (vegxy(nw,m) == pitype .and. wtxy(nw,m) > 0._r8) then
                      pptr%wtgcell(pi)  = pptr%wtgcell(pi) + wtxy(nw,m)
                      pptr%wtlunit(pi)  = pptr%wtlunit(pi) + wtxy(nw,m) / wtlunit2gcell
                      pptr%wtcol(pi)  = pptr%wtcol(pi) + wtxy(nw,m) / wtlunit2gcell
                   end if
                end do
             endif ! setdata
          end do
       else if (allocate_all_vegpfts) then
          do n = 1,numpft+1
             pi = pi + 1
             pitype = n-1
             if (setdata) then
                pptr%mxy(pi)      = n
                pptr%itype(pi)    = pitype
                pptr%gridcell(pi) = gi
                pptr%landunit(pi) = li
                pptr%column (pi) = ci
                pptr%wtgcell(pi) = 0.0_r8
                pptr%wtlunit(pi) = 0.0_r8
                pptr%wtcol(pi) = 0.0_r8
                do m = 1,maxpatch_pft
                   if (vegxy(nw,m) == pitype .and. wtxy(nw,m) > 0._r8) then
                      pptr%wtgcell(pi)  = pptr%wtgcell(pi) + wtxy(nw,m)
                      pptr%wtlunit(pi)  = pptr%wtlunit(pi) + wtxy(nw,m) / wtlunit2gcell
                      pptr%wtcol(pi)  = pptr%wtcol(pi) + wtxy(nw,m) / wtlunit2gcell
                   end if
                end do
             endif ! setdata
          end do
       else
          do m = 1,maxpatch_pft
             if (wtxy(nw,m) > 0._r8) then
                pi = pi + 1
                if (setdata) then
                   pptr%mxy(pi)      = m
                   pptr%itype(pi)    = vegxy(nw,m)
                   pptr%gridcell(pi) = gi
                   pptr%wtgcell(pi) = wtxy(nw,m)
                   pptr%landunit(pi) = li
                   pptr%wtlunit(pi) = wtxy(nw,m) / wtlunit2gcell
                   pptr%column (pi) = ci
                   pptr%wtcol(pi) = wtxy(nw,m) / wtlunit2gcell
                endif ! setdata
             end if
          end do
       end if

    end if

  end subroutine set_landunit_veg_compete
  
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_wet_ice_lake
!
! !INTERFACE:
!  subroutine set_landunit_wet_ice_lake (ltype, wtxy, vegxy, &
  subroutine set_landunit_wet_ice_lake (ltype, &
                           nw, gi, li, ci, pi, setdata, glcmask)
!
! !DESCRIPTION: 
! Initialize wet_ice_lake landunits that are non-urban (lake, wetland, glacier, glacier_mec)
!
! !USES
    use clmtype 
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varcon, only : istice, istwet, istdlak, istice_mec
    use clm_varpar, only : npatch_lake, npatch_glacier, npatch_wet
    use clm_varpar, only : npatch_glacier_mec

!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
!    real(r8), intent(in)    :: wtxy(:,:)         ! subgrid patch weights
!    integer , intent(in)    :: vegxy(:,:)        ! PFT types 
    integer , intent(in)    :: nw                ! cell index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
    integer , intent(in), optional :: glcmask    ! = 1 where glc requires sfc mass balance
                                                 ! = 0 otherwise
!
! !REVISION HISTORY:
! Created by Sam Levis
! 2005.11.25 Updated by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: m                                ! m index in wtxy(nw,m)
    integer  :: c                                ! column loop index
    integer  :: ctype                            ! column type
    integer  :: ier                              ! error status 
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    real(r8) :: wtcol2lunit                      ! col weight in landunit
    type(landunit_type), pointer :: lptr         ! pointer to landunit
    type(column_type)  , pointer :: cptr         ! pointer to column
    type(pft_type)     , pointer :: pptr         ! pointer to pft

!------------------------------------------------------------------------

    ! Set decomposition properties

    if (ltype == istwet) then
!       call subgrid_get_gcellinfo(nw, wtxy, nwetland=npfts, wtwetland=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nwetland=npfts, wtwetland=wtlunit2gcell)
       m = npatch_wet
    else if (ltype == istdlak) then
!       call subgrid_get_gcellinfo(nw, wtxy, nlake=npfts, wtlake=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nlake=npfts, wtlake=wtlunit2gcell)
       m = npatch_lake
    else if (ltype == istice) then 
!       call subgrid_get_gcellinfo(nw, wtxy, nglacier=npfts, wtglacier=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nglacier=npfts, wtglacier=wtlunit2gcell)
       m = npatch_glacier
    else if (ltype == istice_mec) then
!       call subgrid_get_gcellinfo(nw, wtxy, nglacier_mec=npfts, wtglacier_mec=wtlunit2gcell)
       call subgrid_get_gcellinfo(nw, nglacier_mec=npfts, wtglacier_mec=wtlunit2gcell, &
                                  glcmask = glcmask)
       ! NOTE: multiple columns per landunit, so m is not set here

    else
       write(iulog,*)' set_landunit_wet_ice_lake: ltype of ',ltype,' not valid'
       write(iulog,*)' only istwet, istdlak, istice and istice_mec ltypes are valid'
       call endrun()
    end if

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       lptr => lun
       cptr => col
       pptr => pft

       if (npfts /=1 .and. ltype /= istice_mec) then
          write(iulog,*)' set_landunit_wet_ice_lake: compete landunit must'// &
                    ' have one column and one pft '
          write(iulog,*)' current values of ncols, pfts=',ncols,npfts
          call endrun()
       end if

       if (ltype==istice_mec) then   ! multiple columns per landunit

          ! Assume that columns are of type 1 and that each column has its own pft

          ctype = 1
          li = li + 1

          if (setdata) then

             ! Determine landunit properties

             lptr%itype    (li) = ltype
             lptr%ifspecial(li) = .true.
             lptr%glcmecpoi(li) = .true.
             lptr%lakpoi   (li) = .false.
             lptr%urbpoi   (li) = .false.
             lptr%gridcell (li) = gi
             lptr%wtgcell  (li) = wtlunit2gcell

             ! Determine column and properties
             ! (Each column has its own pft)
             ! 
             ! For grid cells with glcmask = 1, make sure all the elevations classes
             !  are populated, even if some have zero fractional area.  This ensures that the 
             !  ice sheet component, glc, will receive a surface mass balance in each elevation 
             !  class wherever the SMB is needed.
             ! Columns with zero weight are referred to as "virtual" columns.
 
             do m = npatch_glacier+1, npatch_glacier_mec

                if (wtxy(nw,m) > 0._r8 .or. glcmask == 1) then

                   ci = ci + 1
                   pi = pi + 1
                   if (wtlunit2gcell > 0._r8) then
                      wtcol2lunit = wtxy(nw,m)/wtlunit2gcell
                   else   ! virtual landunit
                      wtcol2lunit = 0._r8
                   endif

                   cptr%itype    (ci) = ctype
                   cptr%gridcell (ci) = gi
                   cptr%wtgcell  (ci) = wtcol2lunit * wtlunit2gcell
                   cptr%landunit (ci) = li
                   cptr%wtlunit  (ci) = wtcol2lunit

                   ! Set sfc elevation too

                   cps%glc_topo(ci) = topoxy(nw,m)

                   ! Set pft properties

                   pptr%mxy      (pi) = m
                   pptr%itype    (pi) = vegxy(nw,m)
                   pptr%gridcell (pi) = gi
                   pptr%wtgcell  (pi) = wtcol2lunit * wtlunit2gcell
                   pptr%landunit (pi) = li
                   pptr%wtlunit  (pi) = wtcol2lunit
                   pptr%column   (pi) = ci
                   pptr%wtcol    (pi) = 1.0_r8

                endif   ! wtxy > 0 or glcmask = 1
             enddo      ! loop over columns
          endif         ! setdata

       else

          ncols = 1

          ! Currently assume that each landunit only has only one column 
          ! (of type 1) and that each column has its own pft
       
          wtcol2lunit = 1.0_r8/ncols
          ctype = 1

          li = li + 1
          ci = ci + 1
          pi = pi + 1

          if (setdata) then
       
             ! Determine landunit properties 

             lptr%itype    (li) = ltype
             lptr%ifspecial(li) = .true.
             lptr%urbpoi   (li) = .false.
             if (ltype == istdlak) then
                lptr%lakpoi(li) = .true.
             else
                lptr%lakpoi(li) = .false.
             end if
       
             lptr%gridcell (li) = gi
             lptr%wtgcell(li) = wtlunit2gcell

             ! Determine column and properties
             ! For the wet, ice or lake landunits it is assumed that each 
             ! column has its own pft
       
             cptr%itype(ci)    = ctype
       
             cptr%gridcell (ci) = gi
             cptr%wtgcell(ci) = wtcol2lunit * wtlunit2gcell
             cptr%landunit (ci) = li
             cptr%wtlunit(ci) = wtcol2lunit

             ! Set pft properties

             pptr%mxy(pi)      = m
             pptr%itype(pi)    = vegxy(nw,m)
     
             pptr%gridcell (pi) = gi
             pptr%wtgcell(pi) = wtcol2lunit * wtlunit2gcell
             pptr%landunit (pi) = li
             pptr%wtlunit(pi) = wtcol2lunit
             pptr%column (pi) = ci
             pptr%wtcol(pi) = 1.0_r8
          endif ! setdata
       end if   ! ltype = istice_mec
    endif       ! npfts > 0       

  end subroutine set_landunit_wet_ice_lake

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_crop_noncompete
!
! !INTERFACE:
!  subroutine set_landunit_crop_noncompete (ltype, wtxy, vegxy, &
  subroutine set_landunit_crop_noncompete (ltype, &
                           nw, gi, li, ci, pi, setdata)
!
! !DESCRIPTION: 
! Initialize crop landunit without competition
!
! !USES
    use clmtype 
    use subgridMod, only : subgrid_get_gcellinfo
    use clm_varctl, only : create_crop_landunit
    use clm_varpar, only : maxpatch_pft, numcft, npatch_glacier_mec
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
!    real(r8), intent(in)    :: wtxy(:,:)         ! subgrid patch weights
!    integer , intent(in)    :: vegxy(:,:)        ! PFT types 
    integer , intent(in)    :: nw                ! cell index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
!
! !REVISION HISTORY:
! Created by Sam Levis
! 2005.11.25 Updated by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: m                                ! m index in wtxy(nw,m)
    integer  :: npfts                            ! number of pfts in landunit
    integer  :: ncols                            ! number of columns in landu
    real(r8) :: wtlunit2gcell                    ! landunit weight in gridcell
    type(landunit_type), pointer :: lptr         ! pointer to landunit
    type(column_type)  , pointer :: cptr         ! pointer to column
    type(pft_type)     , pointer :: pptr         ! pointer to pft
!------------------------------------------------------------------------

    ! Set decomposition properties

!    call subgrid_get_gcellinfo(nw, wtxy, ncrop=npfts, wtcrop=wtlunit2gcell)
    call subgrid_get_gcellinfo(nw, ncrop=npfts, wtcrop=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       lptr => lun
       cptr => col
       pptr => pft
       
       ! Set landunit properties - each column has its own pft
       
       ncols = npfts
       
       li = li + 1   

       if (setdata) then
          lptr%itype(li)     = ltype
          lptr%ifspecial(li) = .false.
          lptr%lakpoi(li)    = .false.
          lptr%urbpoi(li)    = .false.
          lptr%gridcell (li) = gi
          lptr%wtgcell(li) = wtlunit2gcell
       endif ! setdata

       ! Set column and pft properties for this landunit 
       ! (each column has its own pft)

       if (create_crop_landunit) then
          do m = maxpatch_pft-numcft+1, maxpatch_pft
             ci = ci + 1
             pi = pi + 1
             
             if (setdata) then
                cptr%itype(ci)    = 1
                pptr%itype(pi)    = m - 1
                pptr%mxy(pi)      = m
          
                cptr%gridcell (ci) = gi
                cptr%wtgcell(ci) = wtxy(nw,m)
                cptr%landunit (ci) = li

                pptr%gridcell (pi) = gi
                pptr%wtgcell(pi) = wtxy(nw,m)
                pptr%landunit (pi) = li
                pptr%column (pi) = ci
                if (wtxy(nw,m) > 0._r8) then
                   cptr%wtlunit(ci) = wtxy(nw,m) / wtlunit2gcell
                   pptr%wtlunit(pi) = wtxy(nw,m) / wtlunit2gcell
                   pptr%wtcol(pi) = 1._r8
                else
                   cptr%wtlunit(ci) = 0._r8
                   pptr%wtlunit(pi) = 0._r8
                   pptr%wtcol(pi) = 0._r8
                end if
             endif ! setdata
          end do
       end if

    end if
       
  end subroutine set_landunit_crop_noncompete

!------------------------------------------------------------------------------

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_landunit_urban
!
! !INTERFACE:
!  subroutine set_landunit_urban (ltype, wtxy, vegxy, &
  subroutine set_landunit_urban (ltype, &
                                 nw, gi, li, ci, pi, setdata)
!
! !DESCRIPTION: 
! Initialize urban landunits
!
! !USES
    use clm_varcon   , only : isturb, icol_roof, icol_sunwall, icol_shadewall, &
                              icol_road_perv, icol_road_imperv
    use clm_varpar   , only : npatch_urban, maxpatch_urb
    use clmtype 
    use subgridMod   , only : subgrid_get_gcellinfo
    use UrbanInputMod, only : urbinp
    use decompMod    , only : ldecomp
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: ltype             ! landunit type
!   real(r8), intent(in)    :: wtxy(:,:)         ! subgrid patch weights
!   integer , intent(in)    :: vegxy(:,:)        ! PFT types 
    integer , intent(in)    :: nw                ! cell index
    integer , intent(in)    :: gi                ! gridcell index
    integer , intent(inout) :: li                ! landunit index
    integer , intent(inout) :: ci                ! column index
    integer , intent(inout) :: pi                ! pft index
    logical , intent(in)    :: setdata           ! set info or just compute
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: c             ! column loop index
    integer  :: m             ! m index in wtxy(nw,m)
    integer  :: ctype         ! column type
    integer  :: npfts         ! number of pfts in landunit
    integer  :: ncols         ! number of columns in landunit
    real(r8) :: wtlunit2gcell ! weight relative to gridcell of landunit
    real(r8) :: wtcol2lunit   ! weight of column with respect to landunit
    real(r8) :: wtlunit_roof  ! weight of roof with respect to landunit
    real(r8) :: wtroad_perv   ! weight of pervious road column with respect to total road
    integer  :: ier           ! error status 
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!------------------------------------------------------------------------

    ! Set decomposition properties

!   call subgrid_get_gcellinfo(nw, wtxy, nurban=npfts, wturban=wtlunit2gcell)
    call subgrid_get_gcellinfo(nw, nurban=npfts, wturban=wtlunit2gcell)

    if (npfts > 0) then

       ! Set pointers into derived types for this module

       lptr => lun
       cptr => col
       pptr => pft
       
       ! Determine landunit properties - each columns has its own pft
       
       ncols = npfts

       li = li + 1
       if (setdata) then
          lptr%itype    (li) = ltype
          lptr%ifspecial(li) = .true.
          lptr%lakpoi   (li) = .false.
          lptr%urbpoi   (li) = .true.

          lptr%gridcell (li) = gi
          lptr%wtgcell  (li) = wtlunit2gcell
       endif

       ! Loop through columns for this landunit and set the column and pft properties
       ! For the urban landunits it is assumed that each column has its own pft
       
       do m = npatch_urban, npatch_urban + maxpatch_urb - 1
          if (wtxy(nw,m) > 0._r8) then
                
             wtlunit_roof = urbinp%wtlunit_roof(nw)
             wtroad_perv  = urbinp%wtroad_perv(nw)
             
             if (m == npatch_urban  ) then
                ctype = icol_roof
                wtcol2lunit = wtlunit_roof
             else if (m == npatch_urban+1) then
                ctype = icol_sunwall
                wtcol2lunit = (1. - wtlunit_roof)/3
             else if (m == npatch_urban+2) then
                ctype = icol_shadewall
                wtcol2lunit = (1. - wtlunit_roof)/3
             else if (m == npatch_urban+3) then
                ctype = icol_road_imperv
                wtcol2lunit = ((1. - wtlunit_roof)/3) * (1.-wtroad_perv)
             else if (m == npatch_urban+4) then
                ctype = icol_road_perv
                wtcol2lunit = ((1. - wtlunit_roof)/3) * (wtroad_perv)
             end if
             
             ci = ci + 1
             pi = pi + 1 
             
             if (setdata) then
                cptr%itype(ci)     = ctype

                cptr%gridcell (ci) = gi
                cptr%wtgcell  (ci) = wtcol2lunit * wtlunit2gcell
                cptr%landunit (ci) = li
                cptr%wtlunit  (ci) = wtcol2lunit

                pptr%mxy     (pi)  = m
                pptr%itype   (pi)  = vegxy(nw,m)
                
                pptr%gridcell(pi)  = gi
                pptr%wtgcell (pi)  = wtcol2lunit * wtlunit2gcell
                pptr%landunit(pi)  = li
                pptr%wtlunit (pi)  = wtcol2lunit
                pptr%column  (pi)  = ci
                pptr%wtcol   (pi)  = 1.0_r8
             end if
             
          end if
       end do   ! end of loop through urban columns-pfts

    end if

  end subroutine set_landunit_urban

!------------------------------------------------------------------------------

end module initGridCellsMod
