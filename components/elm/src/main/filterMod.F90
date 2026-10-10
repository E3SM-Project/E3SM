module filterMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module of filters used for processing columns and pfts of particular
  ! types, including lake, non-lake, urban, soil, snow, non-snow, and
  ! naturally-vegetated patches.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use elm_varctl     , only : iulog
  use decompMod      , only : bounds_type  
  use GridcellType   , only : grc_pp
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp                
  use VegetationType , only : veg_pp  
  use TopounitType   , only : top_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type clumpfilter
     integer, pointer :: natvegp(:)      ! nat-vegetated (present) filter (pfts)
     integer :: num_natvegp              ! number of pfts in nat-vegetated filter

     integer, pointer :: pcropp(:)       ! prognostic crop filter (pfts)
     integer :: num_pcropp               ! number of pfts in prognostic crop filter
     integer, pointer :: ppercropp(:)    ! prognostic perennial crop filter (pfts)
     integer :: num_ppercropp            ! number of pfts in prognostic perennial crop filter
     integer, pointer :: soilnopcropp(:) ! soil w/o prog. crops (pfts)
     integer :: num_soilnopcropp         ! number of pfts in soil w/o prog crops

     integer, pointer :: lakep(:)        ! lake filter (pfts)
     integer :: num_lakep                ! number of pfts in lake filter
     integer, pointer :: nolakep(:)      ! non-lake filter (pfts)
     integer :: num_nolakep              ! number of pfts in non-lake filter
     integer, pointer :: lakec(:)        ! lake filter (columns)
     integer :: num_lakec                ! number of columns in lake filter
     integer, pointer :: nolakec(:)      ! non-lake filter (columns)
     integer :: num_nolakec              ! number of columns in non-lake filter

     integer, pointer :: soilc(:)        ! soil filter (columns)
     integer :: num_soilc                ! number of columns in soil filter 
     integer, pointer :: soilp(:)        ! soil filter (pfts)
     integer :: num_soilp                ! number of pfts in soil filter 

     integer, pointer :: snowc(:)        ! snow filter (columns) 
     integer :: num_snowc                ! number of columns in snow filter 
     integer, pointer :: nosnowc(:)      ! non-snow filter (columns) 
     integer :: num_nosnowc              ! number of columns in non-snow filter 

     integer, pointer :: lakesnowc(:)    ! snow filter (columns) 
     integer :: num_lakesnowc            ! number of columns in snow filter 
     integer, pointer :: lakenosnowc(:)  ! non-snow filter (columns) 
     integer :: num_lakenosnowc          ! number of columns in non-snow filter 

     integer, pointer :: hydrologyc(:)   ! hydrology filter (columns)
     integer :: num_hydrologyc           ! number of columns in hydrology filter 

     integer, pointer :: hydrononsoic(:) ! non-soil hydrology filter (columns)
     integer :: num_hydrononsoic         ! number of columns in non-soil hydrology filter

     integer, pointer :: urbanl(:)       ! urban filter (landunits)
     integer :: num_urbanl               ! number of landunits in urban filter 
     integer, pointer :: nourbanl(:)     ! non-urban filter (landunits)
     integer :: num_nourbanl             ! number of landunits in non-urban filter 

     integer, pointer :: urbanc(:)       ! urban filter (columns)
     integer :: num_urbanc               ! number of columns in urban filter
     integer, pointer :: nourbanc(:)     ! non-urban filter (columns)
     integer :: num_nourbanc             ! number of columns in non-urban filter

     integer, pointer :: urbanp(:)       ! urban filter (pfts)
     integer :: num_urbanp               ! number of pfts in urban filter
     integer, pointer :: nourbanp(:)     ! non-urban filter (pfts)
     integer :: num_nourbanp             ! number of pfts in non-urban filter

     integer, pointer :: nolakeurbanp(:) ! non-lake, non-urban filter (pfts)
     integer :: num_nolakeurbanp         ! number of pfts in non-lake, non-urban filter

     integer, pointer :: nolakeurban_barep(:) ! non-lake, non-urban, bare-ground filter (pfts)
     integer :: num_nolakeurban_barep         ! number of pfts in non-lake, non-urban, bare-ground filter
     integer, pointer :: nolakeurban_vegp(:)  ! non-lake, non-urban, vegetated filter (pfts)
     integer :: num_nolakeurban_vegp          ! number of pfts in non-lake, non-urban, vegetated filter

     integer, pointer :: icemecc(:)      ! glacier mec filter (cols)
     integer :: num_icemecc              ! number of columns in glacier mec filter
     
     integer, pointer :: do_smb_c(:)     ! glacier+bareland SMB calculations-on filter (cols)
     integer :: num_do_smb_c             ! number of columns in glacier+bareland SMB mec filter         

  end type clumpfilter
  public clumpfilter

  ! This is the standard set of filters, which should be used in most places in the code.
  ! These filters only include 'active' points.
  type(clumpfilter), allocatable, public :: filter(:)
  
  ! --- DO NOT USING THE FOLLOWING VARIABLE UNLESS YOU KNOW WHAT YOU'RE DOING! ---
  !
  ! This is a separate set of filters that contains both inactive and active points. It is
  ! rarely appropriate to use these, but they are needed in a few places, e.g., where
  ! quantities are computed before weights, active flags and filters are updated due to
  ! landuse change. Note that, for the handful of filters that are computed elsewhere
  ! (including the natvegp filter, the snow filters, and the nolakeurban_barep/
  ! nolakeurban_vegp filters), these filters are NOT
  ! included in this variable - so they can only be used from the main 'filter' variable.
  !
  ! Ideally, we would like to restructure the initialization code and driver ordering so
  ! that this version of the filters is never needed. At that point, we could remove this
  ! filter_inactive_and_active variable, and simplify filterMod to look the way it did
  ! before this variable was added (i.e., when there was only a single group of filters).
  !
  type(clumpfilter), allocatable, public :: filter_inactive_and_active(:)
  !
  public allocFilters            ! allocate memory for filters
  public setFilters              ! set filters
  public setExposedvegpFilters   ! set the dynamic, snow-dependent bare-ground/vegetated sub-filters

  private allocFiltersOneGroup  ! allocate memory for one group of filters
  private setFiltersOneGroup    ! set one group of filters
  !
  ! !REVISION HISTORY:
  ! Created by Mariana Vertenstein
  ! 11/13/03, Peter Thornton: Added soilp and num_soilp
  ! Jan/08, S. Levis: Added crop-related filters
  ! June/13, Bill Sacks: Change main filters to just work over 'active' points; 
  ! add filter_inactive_and_active
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine allocFilters()
    !
    ! !DESCRIPTION:
    ! Allocate CLM filters.
    !
    ! !REVISION HISTORY:
    ! Created by Bill Sacks
    !------------------------------------------------------------------------

    call allocFiltersOneGroup(filter)
    call allocFiltersOneGroup(filter_inactive_and_active)

  end subroutine allocFilters

  !------------------------------------------------------------------------
  subroutine allocFiltersOneGroup(this_filter)
    !
    ! !DESCRIPTION:
    ! Allocate CLM filters, for one group of filters.
    !
    ! !USES:
    use decompMod , only : get_proc_clumps, get_clump_bounds
    !
    ! !ARGUMENTS:
    type(clumpfilter), intent(inout), allocatable :: this_filter(:)  ! the filter to allocate
    !
    ! LOCAL VARAIBLES:
    integer :: nc          ! clump index
    integer :: nclumps     ! total number of clumps on this processor
    integer :: ier         ! error status
    type(bounds_type) :: bounds  
    !------------------------------------------------------------------------

    ! Determine clump variables for this processor

    nclumps = get_proc_clumps()

    ier = 0
    if( .not. allocated(this_filter)) then
       allocate(this_filter(nclumps), stat=ier)
    end if
    if (ier /= 0) then
       write(iulog,*) 'allocFiltersOneGroup(): allocation error for clumpsfilters'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Loop over clumps on this processor

!$OMP PARALLEL DO PRIVATE (nc,bounds)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds)

       allocate(this_filter(nc)%lakep(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%nolakep(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%nolakeurbanp(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%nolakeurban_barep(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%nolakeurban_vegp(bounds%endp-bounds%begp+1))

       allocate(this_filter(nc)%lakec(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%nolakec(bounds%endc-bounds%begc+1))

       allocate(this_filter(nc)%soilc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%soilp(bounds%endp-bounds%begp+1))

       allocate(this_filter(nc)%snowc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%nosnowc(bounds%endc-bounds%begc+1))

       allocate(this_filter(nc)%lakesnowc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%lakenosnowc(bounds%endc-bounds%begc+1))

       allocate(this_filter(nc)%natvegp(bounds%endp-bounds%begp+1))

       allocate(this_filter(nc)%hydrologyc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%hydrononsoic(bounds%endc-bounds%begc+1))

       allocate(this_filter(nc)%urbanp(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%nourbanp(bounds%endp-bounds%begp+1))

       allocate(this_filter(nc)%urbanc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%nourbanc(bounds%endc-bounds%begc+1))

       allocate(this_filter(nc)%urbanl(bounds%endl-bounds%begl+1))
       allocate(this_filter(nc)%nourbanl(bounds%endl-bounds%begl+1))

       allocate(this_filter(nc)%pcropp(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%ppercropp(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%soilnopcropp(bounds%endp-bounds%begp+1))

       allocate(this_filter(nc)%icemecc(bounds%endc-bounds%begc+1))      
       allocate(this_filter(nc)%do_smb_c(bounds%endc-bounds%begc+1))       
       
    end do
!$OMP END PARALLEL DO

  end subroutine allocFiltersOneGroup

  !------------------------------------------------------------------------
  subroutine setFilters(bounds, icemask_grc)
    !
    ! !DESCRIPTION:
    ! Set CLM filters.
    use decompMod , only : BOUNDS_LEVEL_CLUMP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: icemask_grc( bounds%begg: ) ! ice sheet grid coverage mask [gridcell]
    !------------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_CLUMP, errMsg(__FILE__, __LINE__))

    call setFiltersOneGroup(bounds, &
         filter, include_inactive = .false., &
         icemask_grc = icemask_grc(bounds%begg:bounds%endg))

    ! At least as of June, 2013, the 'inactive_and_active' version of the filters is
    ! static in time. Thus, we could have some logic saying whether we're in
    ! initialization, and if so, skip this call. But this is problematic for two reasons:
    ! (1) it requires that the caller of this routine (currently reweight_wrapup) know
    ! whether it is in initialization; and (2) it assumes that the filter definitions
    ! won't be changed in the future in a way that creates some variability in time. So
    ! for now, it seems cleanest and safest to just update these filters whenever the main
    ! filters are updated. But if this proves to be a performance problem, we could
    ! introduce an argument saying whether we're in initialization, and if so, skip this
    ! call.
    
    call setFiltersOneGroup(bounds, &
         filter_inactive_and_active, include_inactive = .true., &
         icemask_grc = icemask_grc(bounds%begg:bounds%endg))
    
  end subroutine setFilters


  !------------------------------------------------------------------------
  subroutine setFiltersOneGroup(bounds, this_filter, include_inactive, icemask_grc )
    !
    ! !DESCRIPTION:
    ! Set CLM filters for one group of filters.
    !
    ! "Standard" filters only include active points. However, this routine can be used to set
    ! alternative filters that also apply over inactive points, by setting include_inactive =
    ! .true.
    !
    ! !USES:
    use decompMod , only : BOUNDS_LEVEL_CLUMP
    use pftvarcon , only : iscft, crop, percrop
    use landunit_varcon, only : istsoil, istcrop, istice_mec
    use column_varcon, only : icol_road_perv
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds  
    type(clumpfilter) , intent(inout) :: this_filter(:)              ! the group of filters to set
    logical           , intent(in)    :: include_inactive            ! whether inactive points should be included in the filters
    real(r8)          , intent(in)    :: icemask_grc( bounds%begg: ) ! ice sheet grid coverage mask [gridcell]
    !
    ! LOCAL VARAIBLES:
    integer :: nc          ! clump index
    integer :: c,l,p       ! column, landunit, pft indices
    integer :: fl          ! lake filter index
    integer :: fnl,fnlu    ! non-lake filter index
    integer :: fs          ! soil filter index
    integer :: fc, fpc     ! crop and perennial crop filter index
    integer :: fnc         ! non-crop filter index
    integer :: f, fn       ! general indices
    integer :: g           !gridcell index
    integer :: t           !topounit index
    !------------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_CLUMP, errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(icemask_grc) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))

    nc = bounds%clump_index

    ! ------------------------------------------------------------------
    ! Each filter below is built in two true-parallel OpenACC passes:
    !  (1) a "reduction" pass that counts how many points satisfy each
    !      predicate (gives num_xxx, independent of execution order), and
    !  (2) a "scatter" pass that uses an atomic-capture counter to claim a
    !      unique slot in the output filter array for each point that
    !      satisfies the predicate (the array *contents* - the set of
    !      points in the filter - do not depend on scatter order, only
    !      set membership matters for how these filters are subsequently
    !      used).
    ! ------------------------------------------------------------------

    ! Create lake and non-lake filters at column-level

    fl = 0
    fnl = 0
    !$acc parallel loop independent gang vector default(present) private(t,l) reduction(+:fl,fnl)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l =col_pp%landunit(c)
             if (lun_pp%lakpoi(l)) then
                fl = fl + 1
             else
                fnl = fnl + 1
             end if
          end if
       end if
    end do
    this_filter(nc)%num_lakec = fl
    this_filter(nc)%num_nolakec = fnl

    fl = 0
    fnl = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,f) copy(fl,fnl)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l =col_pp%landunit(c)
             if (lun_pp%lakpoi(l)) then
                !$acc atomic capture
                fl = fl + 1
                f = fl
                !$acc end atomic
                this_filter(nc)%lakec(f) = c
             else
                !$acc atomic capture
                fnl = fnl + 1
                f = fnl
                !$acc end atomic
                this_filter(nc)%nolakec(f) = c
             end if
          end if
       end if
    end do

    ! Create lake and non-lake filters at pft-level

    fl = 0
    fnl = 0
    fnlu = 0
    !$acc parallel loop independent gang vector default(present) private(t,l) reduction(+:fl,fnl,fnlu)
    do p = bounds%begp,bounds%endp
       t =veg_pp%topounit(p)
       if (top_pp%active(t)) then
          if (veg_pp%active(p) .or. include_inactive) then
             l =veg_pp%landunit(p)
             if (lun_pp%lakpoi(l) ) then
                fl = fl + 1
             else
                fnl = fnl + 1
                if (.not. lun_pp%urbpoi(l)) then
                   fnlu = fnlu + 1
                end if
             end if
          end if
       end if
    end do
    this_filter(nc)%num_lakep = fl
    this_filter(nc)%num_nolakep = fnl
    this_filter(nc)%num_nolakeurbanp = fnlu

    fl = 0
    fnl = 0
    fnlu = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,f) copy(fl,fnl,fnlu)
    do p = bounds%begp,bounds%endp
       t =veg_pp%topounit(p)
       if (top_pp%active(t)) then
          if (veg_pp%active(p) .or. include_inactive) then
             l =veg_pp%landunit(p)
             if (lun_pp%lakpoi(l) ) then
                !$acc atomic capture
                fl = fl + 1
                f = fl
                !$acc end atomic
                this_filter(nc)%lakep(f) = p
             else
                !$acc atomic capture
                fnl = fnl + 1
                f = fnl
                !$acc end atomic
                this_filter(nc)%nolakep(f) = p
                if (.not. lun_pp%urbpoi(l)) then
                   !$acc atomic capture
                   fnlu = fnlu + 1
                   f = fnlu
                   !$acc end atomic
                   this_filter(nc)%nolakeurbanp(f) = p
                end if
             end if
          end if
       end if
    end do

    ! Create soil filter at column-level

    fs = 0
    !$acc parallel loop independent gang vector default(present) private(t,l) reduction(+:fs)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l =col_pp%landunit(c)
             if (col_pp%is_soil(c) .or. col_pp%is_crop(c)) then
                fs = fs + 1
             end if
          end if
       end if
    end do
    this_filter(nc)%num_soilc = fs

    fs = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,f) copy(fs)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l =col_pp%landunit(c)
             if (col_pp%is_soil(c) .or. col_pp%is_crop(c)) then
                !$acc atomic capture
                fs = fs + 1
                f = fs
                !$acc end atomic
                this_filter(nc)%soilc(f) = c
             end if
          end if
       end if
    end do

    ! Create soil filter at pft-level

    fs = 0
    !$acc parallel loop independent gang vector default(present) private(t,l) reduction(+:fs)
    do p = bounds%begp,bounds%endp
       t =veg_pp%topounit(p)
       if (top_pp%active(t)) then
          if (veg_pp%active(p) .or. include_inactive) then
             l =veg_pp%landunit(p)
             if (veg_pp%is_on_soil_col(p) .or. veg_pp%is_on_crop_col(p)) then
                fs = fs + 1
             end if
          end if
       end if
    end do
    this_filter(nc)%num_soilp = fs

    fs = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,f) copy(fs)
    do p = bounds%begp,bounds%endp
       t =veg_pp%topounit(p)
       if (top_pp%active(t)) then
          if (veg_pp%active(p) .or. include_inactive) then
             l =veg_pp%landunit(p)
             if (veg_pp%is_on_soil_col(p) .or. veg_pp%is_on_crop_col(p)) then
                !$acc atomic capture
                fs = fs + 1
                f = fs
                !$acc end atomic
                this_filter(nc)%soilp(f) = p
             end if
          end if
       end if
    end do

    ! Create column-level hydrology filter (soil and Urban pervious road cols)

    f  = 0
    fn = 0
    !$acc parallel loop independent gang vector default(present) private(t,l) reduction(+:f,fn)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l =col_pp%landunit(c)
             if (col_pp%is_soil(c) .or. col_pp%itype(c) == icol_road_perv .or. &
                  col_pp%is_crop(c)) then
                f = f + 1
                if (col_pp%itype(c) == icol_road_perv) then
                   fn = fn + 1
                end if
             end if
          end if
       end if
    end do
    this_filter(nc)%num_hydrologyc = f
    this_filter(nc)%num_hydrononsoic = fn

    f  = 0
    fn = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,fs,fc) copy(f,fn)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l =col_pp%landunit(c)
             if (col_pp%is_soil(c) .or. col_pp%itype(c) == icol_road_perv .or. &
                  col_pp%is_crop(c)) then
                !$acc atomic capture
                f = f + 1
                fs = f
                !$acc end atomic
                this_filter(nc)%hydrologyc(fs) = c

                if (col_pp%itype(c) == icol_road_perv) then
                   !$acc atomic capture
                   fn = fn + 1
                   fc = fn
                   !$acc end atomic
                   this_filter(nc)%hydrononsoic(fc) = c
                end if
             end if
          end if
       end if
    end do

    ! Create prognostic crop and soil w/o prog. crop filters at pft-level
    ! according to where the crop model should be used

    fc  = 0
    fpc = 0
    fnc = 0
    !$acc parallel loop independent gang vector default(present) private(t,l) reduction(+:fc,fpc,fnc)
    do p = bounds%begp,bounds%endp
       t =veg_pp%topounit(p)
       if (top_pp%active(t)) then
          if (veg_pp%active(p) .or. include_inactive) then
             if (.not. iscft(veg_pp%itype(p))) then
                l =veg_pp%landunit(p)
                if (veg_pp%is_on_soil_col(p) .or. veg_pp%is_on_crop_col(p)) then
                   fnc = fnc + 1
                end if
             else
                if (percrop(veg_pp%itype(p)) < 1) then
                   fc = fc + 1
                else if (percrop(veg_pp%itype(p)) >= 1) then
                   fpc = fpc + 1
                end if
             end if
          end if
       end if
    end do
    this_filter(nc)%num_pcropp   = fc
    this_filter(nc)%num_ppercropp   = fpc
    this_filter(nc)%num_soilnopcropp = fnc   ! This wasn't being set before...

    fc  = 0
    fpc = 0
    fnc = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,f) copy(fc,fpc,fnc)
    do p = bounds%begp,bounds%endp
       t =veg_pp%topounit(p)
       if (top_pp%active(t)) then
          if (veg_pp%active(p) .or. include_inactive) then
             if (.not. iscft(veg_pp%itype(p))) then
                l =veg_pp%landunit(p)
                if (veg_pp%is_on_soil_col(p) .or. veg_pp%is_on_crop_col(p)) then
                   !$acc atomic capture
                   fnc = fnc + 1
                   f = fnc
                   !$acc end atomic
                   this_filter(nc)%soilnopcropp(f) = p
                end if
             else
                if (percrop(veg_pp%itype(p)) < 1) then
                   !$acc atomic capture
                   fc = fc + 1
                   f = fc
                   !$acc end atomic
                   this_filter(nc)%pcropp(f) = p
                else if (percrop(veg_pp%itype(p)) >= 1) then
                   !$acc atomic capture
                   fpc = fpc + 1
                   f = fpc
                   !$acc end atomic
                   this_filter(nc)%ppercropp(f) = p
                end if
             end if
          end if
       end if
    end do

    ! Create landunit-level urban and non-urban filters

    f  = 0
    fn = 0
    !$acc parallel loop independent gang vector default(present) private(t) reduction(+:f,fn)
    do l = bounds%begl,bounds%endl
       t =lun_pp%topounit(l)
       if (top_pp%active(t)) then
          if (lun_pp%active(l) .or. include_inactive) then
             if (lun_pp%urbpoi(l)) then
                f = f + 1
             else
                fn = fn + 1
             end if
          end if
       end if
    end do
    this_filter(nc)%num_urbanl = f
    this_filter(nc)%num_nourbanl = fn

    f  = 0
    fn = 0
    !$acc parallel loop independent gang vector default(present) private(t,fs) copy(f,fn)
    do l = bounds%begl,bounds%endl
       t =lun_pp%topounit(l)
       if (top_pp%active(t)) then
          if (lun_pp%active(l) .or. include_inactive) then
             if (lun_pp%urbpoi(l)) then
                !$acc atomic capture
                f = f + 1
                fs = f
                !$acc end atomic
                this_filter(nc)%urbanl(fs) = l
             else
                !$acc atomic capture
                fn = fn + 1
                fs = fn
                !$acc end atomic
                this_filter(nc)%nourbanl(fs) = l
             end if
          end if
       end if
    end do

    ! Create column-level urban and non-urban filters

    f  = 0
    fn = 0
    !$acc parallel loop independent gang vector default(present) private(t,l) reduction(+:f,fn)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l = col_pp%landunit(c)
             if (lun_pp%urbpoi(l)) then
                f = f + 1
             else
                fn = fn + 1
             end if
          end if
       end if
    end do
    this_filter(nc)%num_urbanc = f
    this_filter(nc)%num_nourbanc = fn

    f  = 0
    fn = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,fs) copy(f,fn)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l = col_pp%landunit(c)
             if (lun_pp%urbpoi(l)) then
                !$acc atomic capture
                f = f + 1
                fs = f
                !$acc end atomic
                this_filter(nc)%urbanc(fs) = c
             else
                !$acc atomic capture
                fn = fn + 1
                fs = fn
                !$acc end atomic
                this_filter(nc)%nourbanc(fs) = c
             end if
          end if
       end if
    end do

    ! Create pft-level urban and non-urban filters

    f  = 0
    fn = 0
    !$acc parallel loop independent gang vector default(present) private(t,l) reduction(+:f,fn)
    do p = bounds%begp,bounds%endp
       t =veg_pp%topounit(p)
       if (top_pp%active(t)) then
          if (veg_pp%active(p) .or. include_inactive) then
             l = veg_pp%landunit(p)
             if (lun_pp%urbpoi(l)) then
                f = f + 1
             else
                fn = fn + 1
             end if
          end if
       end if
    end do
    this_filter(nc)%num_urbanp = f
    this_filter(nc)%num_nourbanp = fn

    f  = 0
    fn = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,fs) copy(f,fn)
    do p = bounds%begp,bounds%endp
       t =veg_pp%topounit(p)
       if (top_pp%active(t)) then
          if (veg_pp%active(p) .or. include_inactive) then
             l = veg_pp%landunit(p)
             if (lun_pp%urbpoi(l)) then
                !$acc atomic capture
                f = f + 1
                fs = f
                !$acc end atomic
                this_filter(nc)%urbanp(fs) = p
             else
                !$acc atomic capture
                fn = fn + 1
                fs = fn
                !$acc end atomic
                this_filter(nc)%nourbanp(fs) = p
             end if
          end if
       end if
    end do

    ! Create column-level glacier mec filter

    f = 0
    !$acc parallel loop independent gang vector default(present) private(t,l) reduction(+:f)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l = col_pp%landunit(c)
             if (lun_pp%itype(l) == istice_mec) then
                f = f + 1
             end if
          end if
       end if
    end do
    this_filter(nc)%num_icemecc = f

    f = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,fs) copy(f)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l = col_pp%landunit(c)
             if (lun_pp%itype(l) == istice_mec) then
                !$acc atomic capture
                f = f + 1
                fs = f
                !$acc end atomic
                this_filter(nc)%icemecc(fs) = c
             end if
          end if
       end if
    end do

    ! Create column-level glacier+bareland SMB filter

    f = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,g) reduction(+:f) copyin(icemask_grc)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l = col_pp%landunit(c)
             g = col_pp%gridcell(c)
             if ( lun_pp%itype(l) == istice_mec .or. &
                (col_pp%is_soil(c) .and. icemask_grc(g) > 0.)) then
                f = f + 1
             end if
          end if
       end if
    end do
    this_filter(nc)%num_do_smb_c = f

    f = 0
    !$acc parallel loop independent gang vector default(present) private(t,l,g,fs) copy(f) copyin(icemask_grc)
    do c = bounds%begc,bounds%endc
       t =col_pp%topounit(c)
       if (top_pp%active(t)) then
          if (col_pp%active(c) .or. include_inactive) then
             l = col_pp%landunit(c)
             g = col_pp%gridcell(c)
             if ( lun_pp%itype(l) == istice_mec .or. &
                (col_pp%is_soil(c) .and. icemask_grc(g) > 0.)) then
                !$acc atomic capture
                f = f + 1
                fs = f
                !$acc end atomic
                this_filter(nc)%do_smb_c(fs) = c
             end if
          end if
       end if
    end do

    ! Note: snow filters are reconstructed each time step in
    ! LakeHydrology and SnowHydrology.
    ! Note: the nolakeurban_barep/nolakeurban_vegp filters are also
    ! dynamic (they depend on the snow-dependent frac_veg_nosno flag) and
    ! are reconstructed each time step by setExposedvegpFilters, below.

  end subroutine setFiltersOneGroup

  !------------------------------------------------------------------------
  subroutine setExposedvegpFilters(bounds, this_filter, frac_veg_nosno)
    !
    ! !DESCRIPTION:
    ! Split the non-lake, non-urban pft filter (nolakeurbanp) into a
    ! bare-ground sub-filter (nolakeurban_barep) and a vegetated
    ! sub-filter (nolakeurban_vegp), based on the current, snow-dependent
    ! frac_veg_nosno flag (0 => bare ground, 1 => vegetated). Unlike the
    ! rest of the filters set up in setFiltersOneGroup, this must be
    ! recomputed every time step, since frac_veg_nosno changes with snow
    ! cover - similar to how the snow filters are reconstructed each time
    ! step in LakeHydrology and SnowHydrology.
    !
    ! !USES:
    use decompMod , only : BOUNDS_LEVEL_CLUMP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds
    type(clumpfilter)  , intent(inout) :: this_filter(:)             ! the group of filters to set
    integer            , intent(in)    :: frac_veg_nosno( bounds%begp: ) ! 0 => bare ground, 1 => vegetated [patch]
    !
    ! LOCAL VARIABLES:
    integer :: nc          ! clump index
    integer :: f, p        ! filter index, patch index
    integer :: fbare, fveg ! bare-ground / vegetated filter indices
    integer :: idx         ! captured atomic index
    !------------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_CLUMP, errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(frac_veg_nosno) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))

    nc = bounds%clump_index

    fbare = 0
    fveg  = 0
    !$acc parallel loop independent gang vector default(present) private(p) reduction(+:fbare,fveg)
    do f = 1, this_filter(nc)%num_nolakeurbanp
       p = this_filter(nc)%nolakeurbanp(f)
       if (frac_veg_nosno(p) == 0) then
          fbare = fbare + 1
       else
          fveg = fveg + 1
       end if
    end do
    this_filter(nc)%num_nolakeurban_barep = fbare
    this_filter(nc)%num_nolakeurban_vegp  = fveg

    fbare = 0
    fveg  = 0
    !$acc parallel loop independent gang vector default(present) private(p,idx) copy(fbare,fveg)
    do f = 1, this_filter(nc)%num_nolakeurbanp
       p = this_filter(nc)%nolakeurbanp(f)
       if (frac_veg_nosno(p) == 0) then
          !$acc atomic capture
          fbare = fbare + 1
          idx = fbare
          !$acc end atomic
          this_filter(nc)%nolakeurban_barep(idx) = p
       else
          !$acc atomic capture
          fveg = fveg + 1
          idx = fveg
          !$acc end atomic
          this_filter(nc)%nolakeurban_vegp(idx) = p
       end if
    end do

  end subroutine setExposedvegpFilters

end module filterMod
