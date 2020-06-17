module filterMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module of filters used for processing columns and pfts of particular
  ! types, including lake, non-lake, urban, soil, snow, non-snow, and
  ! naturally-vegetated patches.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use GridcellType   , only : grc_pp
  use LandunitType   , only : lun_pp
  use ColumnType     , only : col_pp
  use VegetationType      , only : veg_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type :: clumpfilter
     integer, pointer :: natvegp(:)   => null()   ! nat-vegetated (present) filter (pfts)
     integer, pointer :: num_natvegp  => null()   ! number of pfts in nat-vegetated filter

     integer, pointer :: pcropp(:)    => null()  ! prognostic crop filter (pfts)
     integer, pointer :: num_pcropp   => null()  ! number of pfts in prognostic crop filter
     integer, pointer :: soilnopcropp(:)  =>null()! soil w/o prog. crops (pfts)
     integer, pointer :: num_soilnopcropp =>null()        ! number of pfts in soil w/o prog crops

     integer, pointer :: lakep(:)  => null()   ! lake filter (pfts)
     integer, pointer :: num_lakep => null()   ! number of pfts in lake filter
     integer, pointer :: nolakep(:)  => null() ! non-lake filter (pfts)
     integer, pointer :: num_nolakep => null() ! number of pfts in non-lake filter
     integer, pointer :: lakec(:)  => null()   ! lake filter (columns)
     integer, pointer :: num_lakec => null()   ! number of columns in lake filter
     integer, pointer :: nolakec(:)  => null() ! non-lake filter (columns)
     integer, pointer :: num_nolakec => null() ! number of columns in non-lake filter

     integer, pointer :: soilc(:)  => null()   ! soil filter (columns)
     integer, pointer :: num_soilc => null()   ! number of columns in soil filter
     integer, pointer :: soilp(:)  => null()   ! soil filter (pfts)
     integer, pointer :: num_soilp => null()   ! number of pfts in soil filter

     integer, pointer :: snowc(:)  => null()   ! snow filter (columns)
     integer, pointer :: num_snowc => null()   ! number of columns in snow filter
     integer, pointer :: nosnowc(:)  => null() ! non-snow filter (columns)
     integer, pointer :: num_nosnowc => null() ! number of columns in non-snow filter

     integer, pointer :: lakesnowc(:)  => null()   ! snow filter (columns)
     integer, pointer :: num_lakesnowc => null()   ! number of columns in snow filter
     integer, pointer :: lakenosnowc(:)  => null() ! non-snow filter (columns)
     integer, pointer :: num_lakenosnowc => null() ! number of columns in non-snow filter

     integer, pointer :: hydrologyc(:)  =>null() ! hydrology filter (columns)
     integer, pointer :: num_hydrologyc =>null()          ! number of columns in hydrology filter

     integer, pointer :: hydrononsoic(:)   => null() ! non-soil hydrology filter (columns)
     integer, pointer :: num_hydrononsoic  => null() ! number of columns in non-soil hydrology filter

     integer, pointer :: urbanl(:)  => null()    ! urban filter (landunits)
     integer, pointer :: num_urbanl => null()             ! number of landunits in urban filter
     integer, pointer :: nourbanl(:) =>null()    ! non-urban filter (landunits)
     integer, pointer :: num_nourbanl=>null()             ! number of landunits in non-urban filter

     integer, pointer :: urbanc(:)    => null()   ! urban filter (columns)
     integer, pointer :: num_urbanc   => null()            ! number of columns in urban filter
     integer, pointer :: nourbanc(:)  => null()   ! non-urban filter (columns)
     integer, pointer :: num_nourbanc => null()            ! number of columns in non-urban filter

     integer, pointer :: urbanp(:)    => null()    ! urban filter (pfts)
     integer, pointer :: num_urbanp   => null()             ! number of pfts in urban filter
     integer, pointer :: nourbanp(:)  => null()    ! non-urban filter (pfts)
     integer, pointer :: num_nourbanp => null()             ! number of pfts in non-urban filter

     integer, pointer :: nolakeurbanp(:)  => null()! non-lake, non-urban filter (pfts)
     integer, pointer :: num_nolakeurbanp => null()         ! number of pfts in non-lake, non-urban filter

     integer, pointer :: icemecc(:)  => null()    ! glacier mec filter (cols)
     integer, pointer :: num_icemecc => null()             ! number of columns in glacier mec filter

     integer, pointer :: do_smb_c(:)  => null()  ! glacier+bareland SMB calculations-on filter (cols)
     integer, pointer :: num_do_smb_c => null()           ! number of columns in glacier+bareland SMB mec filter

  end type clumpfilter
  public clumpfilter

  ! This is the standard set of filters, which should be used in most places in the code.
  ! These filters only include 'active' points.
  type(clumpfilter), allocatable, public :: filter(:)
  !$acc declare create(filter)

  ! --- DO NOT USING THE FOLLOWING VARIABLE UNLESS YOU KNOW WHAT YOU'RE DOING! ---
  !
  ! This is a separate set of filters that contains both inactive and active points. It is
  ! rarely appropriate to use these, but they are needed in a few places, e.g., where
  ! quantities are computed before weights, active flags and filters are updated due to
  ! landuse change. Note that, for the handful of filters that are computed elsewhere
  ! (including the natvegp filter and the snow filters), these filters are NOT
  ! included in this variable - so they can only be used from the main 'filter' variable.
  !
  ! Ideally, we would like to restructure the initialization code and driver ordering so
  ! that this version of the filters is never needed. At that point, we could remove this
  ! filter_inactive_and_active variable, and simplify filterMod to look the way it did
  ! before this variable was added (i.e., when there was only a single group of filters).
  !
  type(clumpfilter), allocatable, public :: filter_inactive_and_active(:)
  !$acc declare create(filter_inactive_and_active)
  public allocFilters   ! allocate memory for filters
  public setFilters     ! set filters

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
    use decompMod , only : get_clump_bounds
    use decompMod,  only : procinfo
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

    nclumps = procinfo%nclumps

    ier = 0
    if( .not. allocated(this_filter)) then
       allocate(this_filter(nclumps), stat=ier)
    end if
    if (ier /= 0) then
       write(*,*) 'allocFiltersOneGroup(): allocation error for clumpsfilters'
       stop
    end if

    ! Loop over clumps on this processor

    !$OMP PARALLEL DO PRIVATE (nc,bounds)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds)
       allocate(this_filter(nc)%num_lakep  )
       allocate(this_filter(nc)%num_nolakep )
       allocate(this_filter(nc)%num_nolakeurbanp )
       allocate(this_filter(nc)%num_lakec  )
       allocate(this_filter(nc)%num_nolakec)
       allocate(this_filter(nc)%num_soilc)
       allocate(this_filter(nc)%num_soilp)
       allocate(this_filter(nc)%num_snowc)
       allocate(this_filter(nc)%num_nosnowc)
       allocate(this_filter(nc)%num_lakesnowc)
       allocate(this_filter(nc)%num_lakenosnowc)
       allocate(this_filter(nc)%num_natvegp)

       allocate(this_filter(nc)%lakep(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%nolakep(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%nolakeurbanp(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%lakec(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%nolakec(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%soilc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%soilp(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%snowc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%nosnowc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%lakesnowc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%lakenosnowc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%natvegp(bounds%endp-bounds%begp+1))

       allocate(this_filter(nc)%num_hydrologyc)
       allocate(this_filter(nc)%num_hydrononsoic)
       allocate(this_filter(nc)%num_urbanp)
       allocate(this_filter(nc)%num_nourbanp)
       allocate(this_filter(nc)%num_urbanc)
       allocate(this_filter(nc)%num_nourbanc)
       allocate(this_filter(nc)%num_urbanl)
       allocate(this_filter(nc)%num_nourbanl)
       allocate(this_filter(nc)%num_pcropp)
       allocate(this_filter(nc)%num_soilnopcropp)
       allocate(this_filter(nc)%num_icemecc)
       allocate(this_filter(nc)%num_do_smb_c)

       allocate(this_filter(nc)%hydrologyc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%hydrononsoic(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%urbanp(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%nourbanp(bounds%endp-bounds%begp+1))
       allocate(this_filter(nc)%urbanc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%nourbanc(bounds%endc-bounds%begc+1))
       allocate(this_filter(nc)%urbanl(bounds%endl-bounds%begl+1))
       allocate(this_filter(nc)%nourbanl(bounds%endl-bounds%begl+1))
       allocate(this_filter(nc)%pcropp(bounds%endp-bounds%begp+1))
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
      !$acc routine seq
    use decompMod , only : BOUNDS_LEVEL_CLUMP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds
    real(r8)          , intent(in) :: icemask_grc(bounds%begg: ) ! ice sheet grid coverage mask [gridcell]
    !------------------------------------------------------------------------

    call setFiltersOneGroup(bounds, &
         filter, include_inactive = .false., &
         icemask_grc = icemask_grc(bounds%begg:bounds%endg) )

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
      !$acc routine seq
    use decompMod , only : BOUNDS_LEVEL_CLUMP
    use pftvarcon , only : npcropmin
    use landunit_varcon, only : istsoil, istcrop, istice_mec
    use column_varcon, only : icol_road_perv
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds
    type(clumpfilter) , intent(inout) :: this_filter(:)              ! the group of filters to set
    logical           , intent(in)    :: include_inactive            ! whether inactive points should be included in the filters
    real(r8)          , intent(in)    :: icemask_grc(bounds%begg: ) ! ice sheet grid coverage mask [gridcell]
    !
    ! LOCAL VARAIBLES:
    integer :: nc          ! clump index
    integer :: c,l,p       ! column, landunit, pft indices
    integer :: fl          ! lake filter index
    integer :: fnl,fnlu    ! non-lake filter index
    integer :: fs          ! soil filter index
    integer :: f, fn       ! general indices
    integer :: g           ! gridcell index
    !------------------------------------------------------------------------

    nc = bounds%clump_index

    ! Create lake and non-lake filters at column-level

    fl  = 0
    fnl = 0
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .or. include_inactive) then
          l =col_pp%landunit(c)
          if (lun_pp%lakpoi(l)) then
             fl = fl + 1
             this_filter(nc)%lakec(fl) = c
          else
             fnl = fnl + 1
            ! print *, "nc, c", nc, c
             this_filter(nc)%nolakec(fnl) = c
          end if
       end if
    end do
    this_filter(nc)%num_lakec = fl
    this_filter(nc)%num_nolakec = fnl

    ! Create lake and non-lake filters at pft-level

    fl  = 0
    fnl = 0
    fnlu = 0

    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .or. include_inactive) then
          l =veg_pp%landunit(p)

          if (lun_pp%lakpoi(l) ) then
             fl = fl + 1
             this_filter(nc)%lakep(fl) = p
          else
             fnl = fnl + 1
             this_filter(nc)%nolakep(fnl) = p
             if (.not. lun_pp%urbpoi(l)) then
                fnlu = fnlu + 1
                this_filter(nc)%nolakeurbanp(fnlu) = p

             end if
          end if
       end if
    end do
    this_filter(nc)%num_lakep = fl
    this_filter(nc)%num_nolakep = fnl
    this_filter(nc)%num_nolakeurbanp = fnlu

    ! Create soil filter at column-level

    fs = 0
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .or. include_inactive) then
          l =col_pp%landunit(c)
          if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
             fs = fs + 1
             this_filter(nc)%soilc(fs) = c
          end if
       end if
    end do
    this_filter(nc)%num_soilc = fs

    ! Create soil filter at pft-level

    fs = 0
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .or. include_inactive) then
          l =veg_pp%landunit(p)
          if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
             fs = fs + 1
             this_filter(nc)%soilp(fs) = p
          end if
       end if
    end do
    this_filter(nc)%num_soilp = fs

    ! Create column-level hydrology filter (soil and Urban pervious road cols)

    f = 0
    fn= 0
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .or. include_inactive) then
          l =col_pp%landunit(c)
          if (lun_pp%itype(l) == istsoil .or. col_pp%itype(c) == icol_road_perv .or. &
               lun_pp%itype(l) == istcrop) then
             f = f + 1
             this_filter(nc)%hydrologyc(f) = c

             if (col_pp%itype(c) == icol_road_perv) then
                fn = fn + 1
                this_filter(nc)%hydrononsoic(fn) = c
             end if

          end if
       end if
    end do
    this_filter(nc)%num_hydrologyc = f
    this_filter(nc)%num_hydrononsoic = fn

    ! Create prognostic crop and soil w/o prog. crop filters at pft-level
    ! according to where the crop model should be used

    fl  = 0
    fnl = 0
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .or. include_inactive) then
          if (veg_pp%itype(p) >= npcropmin) then !skips 2 generic crop types
             fl = fl + 1
             this_filter(nc)%pcropp(fl) = p
          else
             l =veg_pp%landunit(p)
             if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
                fnl = fnl + 1
                this_filter(nc)%soilnopcropp(fnl) = p
             end if
          end if
       end if
    end do
    this_filter(nc)%num_pcropp   = fl
    this_filter(nc)%num_soilnopcropp = fnl   ! This wasn't being set before...

    ! Create landunit-level urban and non-urban filters

    f = 0
    fn = 0
    do l = bounds%begl,bounds%endl
       if (lun_pp%active(l) .or. include_inactive) then
          if (lun_pp%urbpoi(l)) then
             f = f + 1
             this_filter(nc)%urbanl(f) = l
          else
             fn = fn + 1
             this_filter(nc)%nourbanl(fn) = l
          end if
       end if
    end do
    this_filter(nc)%num_urbanl = f
    this_filter(nc)%num_nourbanl = fn

    ! Create column-level urban and non-urban filters

    f = 0
    fn = 0
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .or. include_inactive) then
          l = col_pp%landunit(c)
          if (lun_pp%urbpoi(l)) then
             f = f + 1
             this_filter(nc)%urbanc(f) = c
          else
             fn = fn + 1
             this_filter(nc)%nourbanc(fn) = c
          end if
       end if
    end do
    this_filter(nc)%num_urbanc = f
    this_filter(nc)%num_nourbanc = fn

    ! Create pft-level urban and non-urban filters

    f = 0
    fn = 0
    do p = bounds%begp,bounds%endp
       if (veg_pp%active(p) .or. include_inactive) then
          l = veg_pp%landunit(p)
          if (lun_pp%urbpoi(l)) then
             f = f + 1
             this_filter(nc)%urbanp(f) = p
          else
             fn = fn + 1
             this_filter(nc)%nourbanp(fn) = p
          end if
       end if
    end do
    this_filter(nc)%num_urbanp = f
    this_filter(nc)%num_nourbanp = fn

    f = 0
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .or. include_inactive) then
          l = col_pp%landunit(c)
          if (lun_pp%itype(l) == istice_mec) then
             f = f + 1
             this_filter(nc)%icemecc(f) = c
          end if
       end if
    end do
    this_filter(nc)%num_icemecc = f

    f = 0
    do c = bounds%begc,bounds%endc
       if (col_pp%active(c) .or. include_inactive) then
          l = col_pp%landunit(c)
          g = col_pp%gridcell(c)
          if ( lun_pp%itype(l) == istice_mec .or. &
             (lun_pp%itype(l) == istsoil .and. icemask_grc(g) > 0.)) then
             f = f + 1
             this_filter(nc)%do_smb_c(f) = c
          end if
       end if
    end do
    this_filter(nc)%num_do_smb_c = f

    ! Note: snow filters are reconstructed each time step in
    ! LakeHydrology and SnowHydrology

  end subroutine setFiltersOneGroup


end module filterMod
