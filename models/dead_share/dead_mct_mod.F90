module dead_mct_mod

  use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8
  use shr_sys_mod , only: shr_sys_abort
  use seq_flds_mod, only: seq_flds_dom_coord, seq_flds_dom_other
  use shr_file_mod, only: shr_file_getlogunit
  use mct_mod
  use dead_data_mod

  implicit none
  public
  save 

  public :: dead_domain_mct

  interface dead_domain_mct ; module procedure &
     dead_domain_mct_new
  end interface

!===============================================================================
contains
!===============================================================================

subroutine dead_domain_mct_new( mpicom, gbuf, gsMap, domain )

    !-------------------------------------------------------------------
    !---arguments--- 
    integer(IN)    , intent(in)   :: mpicom
    real(R8)       , intent(in)   :: gbuf(:,:)
    type(mct_gsMap), intent(in)   :: gsMap
    type(mct_ggrid), intent(out)  :: domain  

    !---local variables---
    integer(IN)  :: my_task               ! mpi task within communicator
    integer(IN)  :: lsize                 ! temporary
    integer(IN)  :: n,j,i                 ! indices	
    integer(IN)  :: ier                   ! error status
    real(R8)    , pointer  :: data(:)     ! temporary
    integer(IN) , pointer  :: idata(:)    ! temporary
    integer(IN)  :: logunit
    !-------------------------------------------------------------------

    call shr_file_getlogunit(logunit)

    !
    ! Initialize mct dead domain
    !
    call mct_gGrid_init( GGrid=domain, CoordChars=trim(seq_flds_dom_coord), &
         OtherChars=trim(seq_flds_dom_other), &
         lsize=mct_gsMap_lsize(gsMap, mpicom) )
    call mct_aVect_zero(domain%data)
    !
    ! Allocate memory
    !
    lsize = mct_gGrid_lsize(domain)
    if (size(gbuf,dim=1) /= lsize) then
       call shr_sys_abort('mct_dead_domain size error')       
    endif
    allocate(data(lsize))
    allocate(idata(lsize))
    !
    ! Initialize attribute vector with special value
    !
    call mpi_comm_rank(mpicom, my_task, ier)
    call mct_gsMap_orderedPoints(gsMap, my_task, idata)
    call mct_gGrid_importIAttr(domain,'GlobGridNum',idata,lsize)
    !
    call mct_aVect_zero(domain%data)
    data(:) = -9999.0_R8  ! generic special value 	
    call mct_gGrid_importRAttr(domain,"lat" ,data,lsize) 
    call mct_gGrid_importRAttr(domain,"lon" ,data,lsize) 
    call mct_gGrid_importRAttr(domain,"area",data,lsize) 
    call mct_gGrid_importRAttr(domain,"frac",data,lsize) 

    data(:) = 0.0_R8  ! generic special value 	
    call mct_gGrid_importRAttr(domain,"mask" ,data,lsize) 
    call mct_gGrid_importRAttr(domain,"aream",data,lsize) 
    !
    ! Fill in correct values for domain components
    !
    do n = 1,lsize
      data(n) = gbuf(n,dead_grid_lat) 
    enddo
    call mct_gGrid_importRAttr(domain,"lat",data,lsize) 

    do n = 1,lsize
      data(n) = gbuf(n,dead_grid_lon) 
    enddo
    call mct_gGrid_importRAttr(domain,"lon",data,lsize) 

    do n = 1,lsize
      data(n) = gbuf(n,dead_grid_area)
    enddo
    call mct_gGrid_importRAttr(domain,"area",data,lsize) 
    call mct_gGrid_importRAttr(domain,"aream",data,lsize) 

    do n = 1,lsize
      data(n) = gbuf(n,dead_grid_mask)
    enddo
    call mct_gGrid_importRAttr(domain,"mask"   ,data,lsize) 

    do n = 1,lsize
      data(n) = gbuf(n,dead_grid_frac)
    enddo
    call mct_gGrid_importRAttr(domain,"frac"   ,data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine dead_domain_mct_new

end module dead_mct_mod
