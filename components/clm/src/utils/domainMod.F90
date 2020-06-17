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
  use clm_varctl  , only : iulog
!
! !PUBLIC TYPES:
  implicit none
  private
!
  public :: domain_type

  !--- this typically contains local domain info with arrays dim begg:endg ---
  type domain_type
     integer         , pointer  :: ns        => null() ! global size of domain
     integer         , pointer  :: ni,nj     => null() ! global axis if 2d (nj=1 if unstructured)
     logical         , pointer  :: isgrid2d  => null() ! true => global grid is lat/lon
     integer         , pointer  :: nbeg,nend => null() ! local beg/end indices
     character(len=8), pointer  :: clmlevel  => null() ! grid type
     integer ,pointer :: mask(:)    => null() ! land mask: 1 = land, 0 = ocean
     real(r8),pointer :: frac(:)    => null() ! fractional land
     real(r8),pointer :: topo(:)    => null() ! topography
     real(r8),pointer :: latc(:)    => null() ! latitude of grid cell (deg)
     real(r8),pointer :: lonc(:)    => null() ! longitude of grid cell (deg)
     real(r8),pointer :: firrig(:)  => null()
     real(r8),pointer :: f_surf(:)  => null() ! fraction of water withdraws from surfacewater
     real(r8),pointer :: f_grd(:)   => null() ! fraction of water withdraws from groundwater
     real(r8),pointer :: xCell(:)   => null() ! x-position of grid cell (m)
     real(r8),pointer :: yCell(:)   => null() ! y-position of grid cell (m)
     real(r8),pointer :: area(:)    => null() ! grid cell area (km**2)
     integer ,pointer :: pftm(:)    => null() ! pft mask: 1=real, 0=fake, -1=notset
     integer ,pointer :: glcmask(:) => null() ! glc mask: 1=sfc mass balance required by GLC component
                                    ! 0=SMB not required (default)
                                    ! (glcmask is just a guess at the appropriate mask, known at initialization - in contrast to icemask, which is the true mask obtained from glc)
     character(len=16),pointer  :: set     => null()   ! flag to check if domain is set
     logical          ,pointer :: decomped => null()  ! decomposed locally or global copy

     ! pflotran:beg-----------------------------------------------------
     integer , pointer :: nv        => null()   ! number of vertices
     real(r8),pointer :: latv(:,:)  => null()   ! latitude of grid cell's vertices (deg)
     real(r8),pointer :: lonv(:,:)  => null()   ! longitude of grid cell's vertices (deg)
     real(r8),pointer :: lon0 => null()        ! the origin lon/lat (Most western/southern corner, if not globally covered grids; OR -180W(360E)/-90N)
     real(r8),pointer :: lat0 => null()        ! the origin lon/lat (Most western/southern corner, if not globally covered grids; OR -180W(360E)/-90N)

     ! pflotran:end-----------------------------------------------------
  end type domain_type

  type(domain_type)    , public :: ldomain
  !$acc declare create(ldomain)

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
  character*16,parameter :: set   = 'domain_set      '
  character*16,parameter :: unset = 'NOdomain_unsetNO'
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
    !#py !#py use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    !
    ! !DESCRIPTION:
    ! This subroutine allocates and nans the domain type
    !
    ! !USES:
    use clm_varcon , only : spval
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
    allocate(domain%ns       )
    allocate(domain%ni,domain%nj)
    allocate(domain%isgrid2d )
    allocate(domain%nbeg,domain%nend)
    allocate(domain%clmlevel)
    allocate(domain%set     )
    allocate(domain%decomped)
    
    if(associated(domain%nv)) then
            write(iulog,*) "already set nv"
    else
        print *, "ERROR nv should already be set"
    end if 
    
    allocate(domain%lon0, domain%lat0)
    
    print *, "done allocating pointers"
    nb = 1
    ne = ni*nj
    if (present(nbeg)) then
       if (present(nend)) then
          nb = nbeg
          ne = nend
       endif
    endif

    if (domain%set == set) then
       call domain_clean(domain)
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
           domain%lonv     = spval
       endif
       if(.not.associated(domain%latv)) then
           allocate(domain%latv(nb:ne, 1:domain%nv))
           if (ier /= 0) &
           call shr_sys_abort('domain_init ERROR: allocate latv ')
           domain%latv     = spval
       endif
    end if
    ! pflotran:end-----------------------------------------------------

    if (present(clmlevel)) then
       domain%clmlevel = clmlevel
    endif
    print *, "making assignments"
    domain%isgrid2d  = isgrid2d
    domain%ns        = ni*nj
    domain%ni        = ni
    domain%nj        = nj
    domain%nbeg      = nb
    domain%nend      = ne
    domain%mask  (:) = -9999
    domain%frac  (:) = -1.0e36
    domain%topo  (:) = 0._r8
    domain%latc  (:) = spval
    domain%lonc  (:) = spval
    domain%xCell (:) = spval
    domain%yCell (:) = spval
    domain%area  (:) = spval
    domain%firrig(:) = 0.7_r8
    domain%f_surf(:) = 1.0_r8
    domain%f_grd (:) = 0.0_r8

    domain%set      = set
    if (domain%nbeg == 1 .and. domain%nend == domain%ns) then
       domain%decomped = .false.
    else
       domain%decomped = .true.
    endif

    domain%pftm     = -9999
    domain%glcmask  = 0
    print *, "done iwth domain_init"

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
    if (domain%set == set) then
       if (masterproc) then
          write(iulog,*) 'domain_clean: cleaning ',domain%ni,domain%nj
       endif
       deallocate(domain%mask,domain%frac,domain%latc, &
                  domain%lonc,domain%area,domain%firrig,domain%pftm, &
                  domain%topo,domain%f_surf,domain%f_grd,domain%glcmask,stat=ier)
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

    domain%clmlevel   = unset
    domain%ns         = huge(1)
    domain%ni         = huge(1)
    domain%nj         = huge(1)
    domain%nbeg       = huge(1)
    domain%nend       = huge(1)
    domain%set        = unset
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
    write(*,*) '  domain_check set       = ',trim(domain%set)
    write(*,*) '  domain_check decomped  = ',domain%decomped
    write(*,*) '  domain_check ns        = ',domain%ns
    write(*,*) '  domain_check ni,nj     = ',domain%ni,domain%nj
    write(*,*) '  domain_check clmlevel  = ',trim(domain%clmlevel)
    write(*,*) '  domain_check nbeg,nend = ',domain%nbeg,domain%nend
    write(*,*) '  domain_check lonc      = ',minval(domain%lonc),maxval(domain%lonc)
    write(*,*) '  domain_check latc      = ',minval(domain%latc),maxval(domain%latc)
    write(*,*) '  domain_check mask      = ',minval(domain%mask),maxval(domain%mask)
    write(*,*) '  domain_check frac      = ',minval(domain%frac),maxval(domain%frac)
    write(*,*) '  domain_check topo      = ',minval(domain%topo),maxval(domain%topo)
    write(*,*) '  domain_check firrig    = ',minval(domain%firrig),maxval(domain%firrig)
    write(*,*) '  domain_check f_surf    = ',minval(domain%f_surf),maxval(domain%f_surf)
    write(*,*) '  domain_check f_grd     = ',minval(domain%f_grd),maxval(domain%f_grd)
    write(*,*) '  domain_check area      = ',minval(domain%area),maxval(domain%area)
    write(*,*) '  domain_check pftm      = ',minval(domain%pftm),maxval(domain%pftm)
    write(*,*) '  domain_check glcmask   = ',minval(domain%glcmask),maxval(domain%glcmask)
    write(*,*) ' '
  endif

end subroutine domain_check

!------------------------------------------------------------------------------

end module domainMod
