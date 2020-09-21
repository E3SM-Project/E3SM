module domainMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: domainMod
!
! !DESCRIPTION:
! Module containing 2-d global surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use spmdMod     , only : masterproc
  use elm_varctl  , only : iulog
!
! !PUBLIC TYPES:
  implicit none
  private
!
  public :: domain_type

  !--- this typically contains local domain info with arrays dim begg:endg ---
  type domain_type
     integer          :: ns         ! global size of domain
     integer          :: ni,nj      ! global axis if 2d (nj=1 if unstructured)
     logical          :: isgrid2d   ! true => global grid is lat/lon
     integer          :: nbeg,nend  ! local beg/end indices
     character(len=8) :: clmlevel   ! grid type
     integer ,pointer :: mask(:)    ! land mask: 1 = land, 0 = ocean
     real(r8),pointer :: frac(:)    ! fractional land
     real(r8),pointer :: topo(:)    ! topography
     real(r8),pointer :: latc(:)    ! latitude of grid cell (deg)
     real(r8),pointer :: lonc(:)    ! longitude of grid cell (deg)
     real(r8),pointer :: firrig(:)
     real(r8),pointer :: f_surf(:)  ! fraction of water withdraws from surfacewater
     real(r8),pointer :: f_grd(:)   ! fraction of water withdraws from groundwater
     real(r8),pointer :: xCell(:)   ! x-position of grid cell (m)
     real(r8),pointer :: yCell(:)   ! y-position of grid cell (m)
     real(r8),pointer :: area(:)    ! grid cell area (km**2)
     integer ,pointer :: pftm(:)    ! pft mask: 1=real, 0=fake, -1=notset
     integer ,pointer :: glcmask(:) ! glc mask: 1=sfc mass balance required by GLC component
                                    ! 0=SMB not required (default)
                                    ! (glcmask is just a guess at the appropriate mask, known at initialization - in contrast to icemask, which is the true mask obtained from glc)
     logical          :: set        ! flag to check if domain is set
     logical          :: decomped   ! decomposed locally or global copy

     ! pflotran:beg-----------------------------------------------------
     integer          :: nv           ! number of vertices
     real(r8),pointer :: latv(:,:)    ! latitude of grid cell's vertices (deg)
     real(r8),pointer :: lonv(:,:)    ! longitude of grid cell's vertices (deg)
     real(r8)         :: lon0         ! the origin lon/lat (Most western/southern corner, if not globally covered grids; OR -180W(360E)/-90N)
     real(r8)         :: lat0         ! the origin lon/lat (Most western/southern corner, if not globally covered grids; OR -180W(360E)/-90N)

     ! pflotran:end-----------------------------------------------------
  end type domain_type

  type(domain_type)    , public :: ldomain
  real(r8), allocatable, public :: lon1d(:), lat1d(:) ! 1d lat/lons for 2d grids
!
! !PUBLIC MEMBER FUNCTIONS:
  public domain_init          ! allocates/nans domain types
  public domain_clean         ! deallocates domain types
  public domain_check         ! write out domain info
!
! !REVISION HISTORY:
! Originally clm_varsur by Mariana Vertenstein
! Migrated from clm_varsur to domainMod by T Craig
!
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_init
!
! !INTERFACE:
  subroutine domain_init(domain,isgrid2d,ni,nj,nbeg,nend,clmlevel)
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
!
! !DESCRIPTION:
! This subroutine allocates and nans the domain type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type)   :: domain        ! domain datatype
    logical, intent(in) :: isgrid2d      ! true => global grid is lat/lon
    integer, intent(in) :: ni,nj         ! grid size, 2d
    integer         , intent(in), optional  :: nbeg,nend  ! beg/end indices
    character(len=*), intent(in), optional  :: clmlevel   ! grid type
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer ier
    integer nb,ne
!
!------------------------------------------------------------------------------

    nb = 1
    ne = ni*nj
    if (present(nbeg)) then
       if (present(nend)) then
          nb = nbeg
          ne = nend
       endif
    endif

    if (domain%set) then
       call domain_clean(domain)
    endif
    if (masterproc) then
       write(iulog,*) 'domain_init: ',ni,nj
    endif
    allocate(domain%mask(nb:ne),domain%frac(nb:ne),domain%latc(nb:ne), &
             domain%pftm(nb:ne),domain%area(nb:ne),domain%firrig(nb:ne),domain%lonc(nb:ne), &
             domain%topo(nb:ne),domain%f_surf(nb:ne),domain%f_grd(nb:ne),domain%glcmask(nb:ne), &
             domain%xCell(nb:ne),domain%yCell(nb:ne),stat=ier)
    if (ier /= 0) then
       call shr_sys_abort('domain_init ERROR: allocate mask, frac, lat, lon, area ')
    endif

    ! pflotran:beg-----------------------------------------------------
    ! 'nv' is user-defined, so it must be initialized or assigned value prior to call this subroutine
    if (domain%nv > 0 .and. domain%nv /= huge(1)) then
       if(.not.associated(domain%lonv)) then
           allocate(domain%lonv(nb:ne, 1:domain%nv), stat=ier)
           if (ier /= 0) &
           call shr_sys_abort('domain_init ERROR: allocate lonv ')
           domain%lonv     = nan
       endif
       if(.not.associated(domain%latv)) then
           allocate(domain%latv(nb:ne, 1:domain%nv))
           if (ier /= 0) &
           call shr_sys_abort('domain_init ERROR: allocate latv ')
           domain%latv     = nan
       endif
    end if
    ! pflotran:end-----------------------------------------------------

    if (present(clmlevel)) then
       domain%clmlevel = clmlevel
    endif

    domain%isgrid2d = isgrid2d
    domain%ns       = ni*nj
    domain%ni       = ni
    domain%nj       = nj
    domain%nbeg     = nb
    domain%nend     = ne
    domain%mask     = -9999
    domain%frac     = -1.0e36
    domain%topo     = 0._r8
    domain%latc     = nan
    domain%lonc     = nan
    domain%xCell    = nan
    domain%yCell    = nan
    domain%area     = nan
    domain%firrig   = 0.7_r8    
    domain%f_surf   = 1.0_r8
    domain%f_grd    = 0.0_r8

    domain%set      = .true.
    if (domain%nbeg == 1 .and. domain%nend == domain%ns) then
       domain%decomped = .false.
    else
       domain%decomped = .true.
    endif

    domain%pftm     = -9999
    domain%glcmask  = 0  

end subroutine domain_init
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_clean
!
! !INTERFACE:
  subroutine domain_clean(domain)
!
! !DESCRIPTION:
! This subroutine deallocates the domain type
!
! !ARGUMENTS:
    implicit none
    type(domain_type) :: domain        ! domain datatype
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer ier
!
!------------------------------------------------------------------------------
    if (domain%set) then
       if (masterproc) then
          write(iulog,*) 'domain_clean: cleaning ',domain%ni,domain%nj
       endif
       deallocate(domain%mask,domain%frac,domain%latc, &
                  domain%lonc,domain%area,domain%firrig,domain%pftm, &
                  domain%topo,domain%f_surf,domain%f_grd,domain%glcmask, &
                  domain%xCell,domain%yCell,stat=ier)
       if (ier /= 0) then
          call shr_sys_abort('domain_clean ERROR: deallocate mask, frac, lat, lon, area ')
       endif

       ! pflotran:beg-----------------------------------------------------
       ! 'nv' is user-defined, so it must be initialized or assigned value prior to call this subroutine
       if (domain%nv > 0 .and. domain%nv /= huge(1)) then
          if (associated(domain%lonv)) then
             deallocate(domain%lonv, stat=ier)
             if (ier /= 0) &
             call shr_sys_abort('domain_clean ERROR: deallocate lonv ')
             nullify(domain%lonv)
          endif

          if (associated(domain%latv)) then
             deallocate(domain%latv, stat=ier)
             if (ier /= 0) &
             call shr_sys_abort('domain_clean ERROR: deallocate latv ')
             nullify(domain%latv)
          endif
       endif
       ! pflotran:beg-----------------------------------------------------

    else
       if (masterproc) then
          write(iulog,*) 'domain_clean WARN: clean domain unecessary '
       endif
    endif

    domain%clmlevel   = 'NOdomain_unsetNO'
    domain%ns         = huge(1)
    domain%ni         = huge(1)
    domain%nj         = huge(1)
    domain%nbeg       = huge(1)
    domain%nend       = huge(1)
    domain%set        = .false.
    domain%decomped   = .true.

end subroutine domain_clean
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_check
!
! !INTERFACE:
  subroutine domain_check(domain)
!
! !DESCRIPTION:
! This subroutine write domain info
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(in)  :: domain        ! domain datatype
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

  if (masterproc) then
    write(iulog,*) '  domain_check set       = ',domain%set
    write(iulog,*) '  domain_check decomped  = ',domain%decomped
    write(iulog,*) '  domain_check ns        = ',domain%ns
    write(iulog,*) '  domain_check ni,nj     = ',domain%ni,domain%nj
    write(iulog,*) '  domain_check clmlevel  = ',trim(domain%clmlevel)
    write(iulog,*) '  domain_check nbeg,nend = ',domain%nbeg,domain%nend
    write(iulog,*) '  domain_check lonc      = ',minval(domain%lonc),maxval(domain%lonc)
    write(iulog,*) '  domain_check latc      = ',minval(domain%latc),maxval(domain%latc)
    write(iulog,*) '  domain_check mask      = ',minval(domain%mask),maxval(domain%mask)
    write(iulog,*) '  domain_check frac      = ',minval(domain%frac),maxval(domain%frac)
    write(iulog,*) '  domain_check topo      = ',minval(domain%topo),maxval(domain%topo)
    write(iulog,*) '  domain_check firrig    = ',minval(domain%firrig),maxval(domain%firrig)
    write(iulog,*) '  domain_check f_surf    = ',minval(domain%f_surf),maxval(domain%f_surf)
    write(iulog,*) '  domain_check f_grd     = ',minval(domain%f_grd),maxval(domain%f_grd)
    write(iulog,*) '  domain_check area      = ',minval(domain%area),maxval(domain%area)
    write(iulog,*) '  domain_check pftm      = ',minval(domain%pftm),maxval(domain%pftm)
    write(iulog,*) '  domain_check glcmask   = ',minval(domain%glcmask),maxval(domain%glcmask)
    write(iulog,*) ' '
  endif

end subroutine domain_check

!------------------------------------------------------------------------------

end module domainMod
