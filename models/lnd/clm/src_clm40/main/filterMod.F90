module filterMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: filterMod
!
! !DESCRIPTION:
! Module of filters used for processing columns and pfts of particular
! types, including lake, non-lake, urban, soil, snow, non-snow, and
! naturally-vegetated patches.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils, only : endrun
  use clm_varctl, only : iulog, use_cndv
!
! !PUBLIC TYPES:
  implicit none
  save

  private

  type clumpfilter
     integer, pointer :: natvegp(:)      ! CNDV nat-vegetated (present) filter (pfts)
     integer :: num_natvegp              ! number of pfts in nat-vegetated filter

     integer, pointer :: pcropp(:)       ! prognostic crop filter (pfts)
     integer :: num_pcropp               ! number of pfts in prognostic crop filter
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

     integer, pointer :: hydrologyc(:)   ! hydrology filter (columns)
     integer :: num_hydrologyc           ! number of columns in hydrology filter 

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

  end type clumpfilter
  public clumpfilter

  type(clumpfilter), allocatable, public :: filter(:)
!
  public allocFilters   ! allocate memory for filters
  public setFilters     ! set filters
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 11/13/03, Peter Thornton: Added soilp and num_soilp
! Jan/08, S. Levis: Added crop-related filters
!
!EOP
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocFilters
!
! !INTERFACE:
  subroutine allocFilters()
!
! !DESCRIPTION:
! Allocate CLM filters.
!
! !USES:
    use clmtype
    use decompMod , only : get_proc_clumps, get_clump_bounds
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2004.04.27 DGVM naturally-vegetated filter added by Forrest Hoffman
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: nc          ! clump index
    integer :: nclumps     ! total number of clumps on this processor
    integer :: begp, endp  ! per-clump beginning and ending pft indices
    integer :: begc, endc  ! per-clump beginning and ending column indices
    integer :: begl, endl  ! per-clump beginning and ending landunit indices
    integer :: begg, endg  ! per-clump beginning and ending gridcell indices
    integer :: ier         ! error status
!------------------------------------------------------------------------

    ! Determine clump variables for this processor

    nclumps = get_proc_clumps()
    ier = 0
    if( .not. allocated(filter)) then
       allocate(filter(nclumps), stat=ier)
    end if
    if (ier /= 0) then
       write(iulog,*) 'allocFilters(): allocation error for clumpsfilters'
       call endrun
    end if

    ! Loop over clumps on this processor

    !$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
    do nc = 1, nclumps
       call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

       allocate(filter(nc)%lakep(endp-begp+1))
       allocate(filter(nc)%nolakep(endp-begp+1))
       allocate(filter(nc)%nolakeurbanp(endp-begp+1))

       allocate(filter(nc)%lakec(endc-begc+1))
       allocate(filter(nc)%nolakec(endc-begc+1))

       allocate(filter(nc)%soilc(endc-begc+1))
       allocate(filter(nc)%soilp(endp-begp+1))

       allocate(filter(nc)%snowc(endc-begc+1))
       allocate(filter(nc)%nosnowc(endc-begc+1))

       if (use_cndv) then
          allocate(filter(nc)%natvegp(endp-begp+1))
       end if

       allocate(filter(nc)%hydrologyc(endc-begc+1))

       allocate(filter(nc)%urbanp(endp-begp+1))
       allocate(filter(nc)%nourbanp(endp-begp+1))

       allocate(filter(nc)%urbanc(endc-begc+1))
       allocate(filter(nc)%nourbanc(endc-begc+1))

       allocate(filter(nc)%urbanl(endl-begl+1))
       allocate(filter(nc)%nourbanl(endl-begl+1))

       allocate(filter(nc)%pcropp(endp-begp+1))
       allocate(filter(nc)%soilnopcropp(endp-begp+1))
    end do
    !$OMP END PARALLEL DO

  end subroutine allocFilters

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setFilters
!
! !INTERFACE:
  subroutine setFilters( nc )
!
! !DESCRIPTION:
! Set CLM filters.
!
! !USES:
    use clmtype
    use decompMod , only : get_clump_bounds
    use pftvarcon , only : npcropmin
    use clm_varcon, only : istsoil, isturb, icol_road_perv, istice_mec
    use clm_varcon, only : istcrop
!
! !ARGUMENTS:
    implicit none
    integer, intent(IN) :: nc          ! clump index
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2004.04.27 DGVM naturally-vegetated filter added by Forrest Hoffman
! 2005.09.12 Urban related filters added by Mariana Vertenstein
!
!EOP
!
! LOCAL VARAIBLES:
    integer , pointer :: ctype(:) ! column type
    integer :: c,l,p       ! column, landunit, pft indices
    integer :: fl          ! lake filter index
    integer :: fnl,fnlu    ! non-lake filter index
    integer :: fs          ! soil filter index
    integer :: f, fn       ! general indices
    integer :: begp, endp  ! per-clump beginning and ending pft indices
    integer :: begc, endc  ! per-clump beginning and ending column indices
    integer :: begl, endl  ! per-clump beginning and ending landunit indices
    integer :: begg, endg  ! per-clump beginning and ending gridcell indices
!------------------------------------------------------------------------

    ctype => col%itype

    ! Determine clump boundaries

    call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

    ! Create lake and non-lake filters at column-level 

    fl = 0
    fnl = 0
    do c = begc,endc
       l = col%landunit(c)
       if (lun%lakpoi(l)) then
          fl = fl + 1
          filter(nc)%lakec(fl) = c
       else
          fnl = fnl + 1
          filter(nc)%nolakec(fnl) = c
       end if
    end do
    filter(nc)%num_lakec = fl
    filter(nc)%num_nolakec = fnl

    ! Create lake and non-lake filters at pft-level 
    ! Filter will only be active if weight of pft wrt gcell is nonzero

    fl = 0
    fnl = 0
    fnlu = 0
    do p = begp,endp
       l = pft%landunit(p)
       if (pft%wtgcell(p) > 0._r8    &
                   .or.                       &
           lun%itype(l)==istice_mec) then  ! some glacier_mec columns have zero weight

          l = pft%landunit(p)
          if (lun%lakpoi(l) ) then
             fl = fl + 1
             filter(nc)%lakep(fl) = p
          else
             fnl = fnl + 1
             filter(nc)%nolakep(fnl) = p
             if (lun%itype(l) /= isturb) then
                fnlu = fnlu + 1
                filter(nc)%nolakeurbanp(fnlu) = p
             end if
          end if
       end if
    end do
    filter(nc)%num_lakep = fl
    filter(nc)%num_nolakep = fnl
    filter(nc)%num_nolakeurbanp = fnlu

    ! Create soil filter at column-level

    fs = 0
    do c = begc,endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          fs = fs + 1
          filter(nc)%soilc(fs) = c
       end if
    end do
    filter(nc)%num_soilc = fs

    ! Create soil filter at pft-level
    ! Filter will only be active if weight of pft wrt gcell is nonzero

    fs = 0
    do p = begp,endp
       if (pft%wtgcell(p) > 0._r8) then
          l = pft%landunit(p)
          if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
             fs = fs + 1
             filter(nc)%soilp(fs) = p
          end if
       end if
    end do
    filter(nc)%num_soilp = fs

    ! Create column-level hydrology filter (soil and Urban pervious road cols) 

    f = 0
    do c = begc,endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. ctype(c) == icol_road_perv .or. &
           lun%itype(l) == istcrop) then
          f = f + 1
          filter(nc)%hydrologyc(f) = c
       end if
    end do
    filter(nc)%num_hydrologyc = f

    ! Create prognostic crop and soil w/o prog. crop filters at pft-level
    ! according to where the crop model should be used

    fl  = 0
    fnl = 0
    do p = begp,endp
       if (pft%wtgcell(p) > 0._r8) then
          if (pft%itype(p) >= npcropmin) then !skips 2 generic crop types
             fl = fl + 1
             filter(nc)%pcropp(fl) = p
          else
             l = pft%landunit(p)
             if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
                fnl = fnl + 1
                filter(nc)%soilnopcropp(fnl) = p
             end if
          end if
       end if
    end do
    filter(nc)%num_pcropp   = fl
    filter(nc)%soilnopcropp = fnl   ! This wasn't being set before...

    ! Create landunit-level urban and non-urban filters

    f = 0
    fn = 0
    do l = begl,endl
       if (lun%itype(l) == isturb) then
          f = f + 1
          filter(nc)%urbanl(f) = l
       else
          fn = fn + 1
          filter(nc)%nourbanl(fn) = l
       end if
    end do
    filter(nc)%num_urbanl = f
    filter(nc)%num_nourbanl = fn

    ! Create column-level urban and non-urban filters

    f = 0
    fn = 0
    do c = begc,endc
       l = col%landunit(c)
       if (lun%itype(l) == isturb) then
          f = f + 1
          filter(nc)%urbanc(f) = c
       else
          fn = fn + 1
          filter(nc)%nourbanc(fn) = c
       end if
    end do
    filter(nc)%num_urbanc = f
    filter(nc)%num_nourbanc = fn

    ! Create pft-level urban and non-urban filters

    f = 0
    fn = 0
    do p = begp,endp
       l = pft%landunit(p)
       if (lun%itype(l) == isturb .and. pft%wtgcell(p) > 0._r8) then
          f = f + 1
          filter(nc)%urbanp(f) = p
       else
          fn = fn + 1
          filter(nc)%nourbanp(fn) = p 
       end if
    end do
    filter(nc)%num_urbanp = f
    filter(nc)%num_nourbanp = fn

    ! Note: snow filters are reconstructed each time step in Hydrology2
    ! Note: CNDV "pft present" filter is reconstructed each time CNDV is run

  end subroutine setFilters

end module filterMod
