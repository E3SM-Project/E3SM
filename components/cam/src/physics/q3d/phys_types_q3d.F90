!-------------------------------------------------------------------------------
! Q3D physics data types module
!-------------------------------------------------------------------------------
module physics_types_q3d

#ifdef JUNG_TEST 
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use ppgrid,           only: pver
  use constituents,     only: pcnst
  use cam_logfile,      only: iulog
  use cam_abortutils,   only: endrun
  use physics_types,    only: physics_state

  implicit none
  private          ! Make default type private to the module

! Public types:
  public physics_dout

! Public interfaces
  public physics_type_sub_alloc
  public physics_dout_alloc    ! allocate individual components within dout
  public physics_dout_dealloc  ! deallocate individual components within dout

!-------------------------------------------------------------------------------
  type physics_dout
     integer   ::   psetcols=0 ! max number of columns set

     real(r8), dimension(:,:),allocatable ::  &
          dt,   & ! temperature tendency (K/s)
          du,   & ! zonal wind tendency (m/s/s)
          dv      ! meridional wind tendency (m/s/s)
          
     real(r8), dimension(:,:,:),allocatable :: &
          dqt     ! constituent mixing ratio tendency (kg/kg/s) 
     
     real(r8), dimension(:,:),allocatable ::   &
          t,    & ! temperature (K)
          u,    & ! zonal wind (m/s)
          v       ! meridional wind (m/s)
          
     real(r8), dimension(:,:,:),allocatable :: &
          qt      ! constituent mixing ratio (kg/kg moist or dry air depending on type)
     
  end type physics_dout

!===============================================================================
contains
!===============================================================================
  subroutine physics_type_sub_alloc(phys_dout, phys_state, begchunk, endchunk, psetcols)
    implicit none
    
    type(physics_dout), pointer :: phys_dout(:)
    
    type(physics_state), intent(in)  :: phys_state(:)
    integer, intent(in) :: begchunk, endchunk
    integer, intent(in) :: psetcols

    integer :: ierr=0, lchnk
    type(physics_dout), pointer :: dout

    allocate(phys_dout(begchunk:endchunk), stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'physics_dout: phys_dout allocation error = ',ierr
       call endrun('physics_types: failed to allocate physics_dout array')
    end if

    do lchnk=begchunk,endchunk
       call physics_dout_alloc(phys_dout(lchnk),phys_state(lchnk)%psetcols)
    end do

  end subroutine physics_type_sub_alloc
!===============================================================================

subroutine physics_dout_alloc(dout,psetcols)

  use infnan, only : inf, assignment(=)
! allocate the individual dout components

  type(physics_dout), intent(inout)  :: dout

  integer, intent(in) :: psetcols

  integer :: ierr = 0

  dout%psetcols = psetcols

  allocate(dout%dt(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_alloc error: allocation error for dout%dt')

  allocate(dout%du(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_alloc error: allocation error for dout%du')

  allocate(dout%dv(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_alloc error: allocation error for dout%dv')

  allocate(dout%dqt(psetcols,pver,pcnst), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%dqt')
  
  allocate(dout%t(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_alloc error: allocation error for dout%t')

  allocate(dout%u(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_alloc error: allocation error for dout%u')

  allocate(dout%v(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_alloc error: allocation error for dout%v')

  allocate(dout%qt(psetcols,pver,pcnst), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%qt')  
  
  dout%dt(:,:) = inf
  dout%du(:,:) = inf
  dout%dv(:,:) = inf
  dout%dqt(:,:,:) = inf

  dout%t(:,:) = inf
  dout%u(:,:) = inf
  dout%v(:,:) = inf
  dout%qt(:,:,:) = inf
  
end subroutine physics_dout_alloc

!===============================================================================

subroutine physics_dout_dealloc(dout)

! deallocate the individual dout components

  type(physics_dout), intent(inout)  :: dout
  integer :: psetcols
  integer :: ierr = 0

  deallocate(dout%dt, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_dealloc error: deallocation error for dout%dt')

  deallocate(dout%du, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_dealloc error: deallocation error for dout%du')

  deallocate(dout%dv, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_dealloc error: deallocation error for dout%dv')

  deallocate(dout%dqt, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_dealloc error: deallocation error for dout%dqt')

  deallocate(dout%t, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_dealloc error: deallocation error for dout%t')

  deallocate(dout%u, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_dealloc error: deallocation error for dout%u')

  deallocate(dout%v, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_dealloc error: deallocation error for dout%v')

  deallocate(dout%qt, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_dout_dealloc error: deallocation error for dout%qt')

end subroutine physics_dout_dealloc

!===============================================================================
#endif

end module physics_types_q3d
