module spmdGathScatMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: spmdGathScatMod
!
! !DESCRIPTION:
! Perform SPMD gather and scatter operations.
!
! !USES:
  use elm_varcon, only: spval, ispval
  use decompMod, only : get_clmlevel_gsmap
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmdMod
  use mct_mod
  use abortutils, only : endrun
  use elm_varctl, only : iulog
  use perf_mod
!
! !PUBLIC TYPES:
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public  scatter_data_from_master, gather_data_to_master

  interface scatter_data_from_master
     module procedure scatter_1darray_int
     module procedure scatter_1darray_real
  end interface

  interface gather_data_to_master
     module procedure gather_1darray_int
     module procedure gather_1darray_real
  end interface
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
  integer,private,parameter :: debug = 0

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_1darray_int
!
! !INTERFACE:
  subroutine scatter_1darray_int (alocal, aglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to scatter int 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer , pointer            :: alocal(:)       ! local  data (output)
    integer , pointer            :: aglobal(:)      ! global data (input)
    character(len=*) ,intent(in) :: clmlevel    ! type of input grid
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer            :: n1,n2,lb1,ub1,lb2,ub2 ! indices
    integer            :: lsize      ! size of local array
    type(mct_aVect)    :: AVi, AVo   ! attribute vectors
    integer ,pointer   :: adata(:)   ! local data array
    character(len=256) :: rstring    ! real field list string
    character(len=256) :: istring    ! int field list string
    character(len=8)   :: fname      ! arbitrary field name
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    character(len=*),parameter :: subname = 'scatter_1darray_int'

!-----------------------------------------------------------------------

    call t_startf(trim(subname)//'_total')
    call get_clmlevel_gsmap(clmlevel,gsmap)

    lb1 = lbound(alocal,dim=1)
    ub1 = ubound(alocal,dim=1)
    lb2 = 1
    ub2 = 1

    rstring = ""
    istring = ""

    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       if (len_trim(istring) == 0) then
          istring = trim(fname)
       else
          istring = trim(istring)//":"//trim(fname)
       endif
    enddo

    if (masterproc .and. debug > 2) then
       write(iulog,*) trim(subname),' strings:',trim(rstring),' ',trim(istring)
    endif

    if (debug > 1) call t_startf(trim(subname)//'_pack')

    if (masterproc) then
       lsize = size(aglobal,dim=1)
       call mct_aVect_init(AVi,rList=trim(rstring),iList=trim(istring),lsize=lsize)
       allocate(adata(lsize))
       do n2 = lb2,ub2
          adata(1:lsize) = aglobal(1:lsize)
          write(fname,'(a1,i3.3)') 'f',n2-lb2+1
          call mct_aVect_importIattr(AVi,trim(fname),adata,lsize)
       enddo
       deallocate(adata)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_pack')
    if (debug > 1) call t_startf(trim(subname)//'_scat')

    call mct_aVect_scatter(AVi, AVo, gsmap, 0, mpicom)

    if (debug > 1) call t_stopf(trim(subname)//'_scat')
    if (debug > 1) call t_startf(trim(subname)//'_upck')

    lsize = size(alocal,dim=1)
    allocate(adata(lsize))
    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       call mct_aVect_exportIattr(AVo,trim(fname),adata,lsize)
       do n1 = lb1,ub1
          alocal(n1) = adata(n1-lb1+1)
       enddo
    enddo
    deallocate(adata)

    if (debug > 1) call t_stopf(trim(subname)//'_upck')

    if (masterproc) then
       call mct_aVect_clean(AVi)
    endif
    call mct_aVect_clean(AVo)

    call t_stopf(trim(subname)//'_total')

  end subroutine scatter_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_1darray_int
!
! !INTERFACE:
  subroutine gather_1darray_int (alocal, aglobal, clmlevel, missing)
!
! !DESCRIPTION:
! Wrapper routine to gather int 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer , pointer            :: alocal(:)       ! local  data (output)
    integer , pointer            :: aglobal(:)      ! global data (input)
    character(len=*) ,intent(in) :: clmlevel    ! type of input grid
    integer ,optional,intent(in) :: missing     ! missing value
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer            :: n1,n2,lb1,ub1,lb2,ub2 ! indices
    integer            :: lsize      ! size of local array
    type(mct_aVect)    :: AVi, AVo   ! attribute vectors
    integer ,pointer   :: adata(:)   ! temporary data array
    integer ,pointer   :: mvect(:)   ! local array for mask
    character(len=256) :: rstring    ! real field list string
    character(len=256) :: istring    ! int field list string
    character(len=8)   :: fname      ! arbitrary field name
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    character(len=*),parameter :: subname = 'gather_1darray_int'

!-----------------------------------------------------------------------

    call t_startf(trim(subname)//'_total')
    call get_clmlevel_gsmap(clmlevel,gsmap)

    lsize = size(alocal,dim=1)
    lb1 = lbound(alocal,dim=1)
    ub1 = ubound(alocal,dim=1)
    lb2 = 1
    ub2 = 1
   
    rstring = ""
    istring = ""

    if (present(missing)) then
       istring = "mask"
    endif

    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       if (len_trim(istring) == 0) then
          istring = trim(fname)
       else
          istring = trim(istring)//":"//trim(fname)
       endif
    enddo

    if (masterproc .and. debug > 2) then
       write(iulog,*) trim(subname),' strings:',trim(rstring),' ',trim(istring)
    endif

    call mct_aVect_init(AVi,rList=trim(rstring),iList=trim(istring),lsize=lsize)

    if (debug > 1) call t_startf(trim(subname)//'_pack')
    allocate(adata(lsize))
    do n2 = lb2,ub2
       do n1 = lb1,ub1
          adata(n1-lb1+1) = alocal(n1)
       enddo
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       call mct_aVect_importIattr(AVi,trim(fname),adata,lsize)
    enddo
    deallocate(adata)

    if (present(missing)) then
       allocate(mvect(lsize))
       do n1 = lb1,ub1
          mvect(n1-lb1+1) = 1
       enddo
       call mct_aVect_importIattr(AVi,"mask",mvect,lsize)
       deallocate(mvect)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_pack')
    if (debug > 1) call t_startf(trim(subname)//'_gath')

    if (present(missing)) then
! tcx wait for update in mct, then get rid of "mask"
!       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom, missing = missing)
       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
    else
       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_gath')
    if (debug > 1) call t_startf(trim(subname)//'_upck')

    if (masterproc) then
       lsize = size(aglobal,dim=1)
       allocate(adata(lsize))
       do n2 = lb2,ub2
          write(fname,'(a1,i3.3)') 'f',n2-lb2+1
          call mct_aVect_exportIattr(AVo,trim(fname),adata,lsize)
          aglobal(1:lsize) = adata(1:lsize)
       enddo
       deallocate(adata)
       if (present(missing)) then
          allocate(mvect(lsize))
          call mct_aVect_exportIattr(AVo,"mask",mvect,lsize)
          do n1 = 1,lsize
             if (mvect(n1) == 0) then
                do n2 = lb2,ub2
                   aglobal(n1) = missing
                enddo
             endif
          enddo
          deallocate(mvect)
       endif
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_upck')

    if (masterproc) then
       call mct_aVect_clean(AVo)
    endif

    call mct_aVect_clean(AVi)

    call t_stopf(trim(subname)//'_total')

  end subroutine gather_1darray_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: scatter_1darray_real
!
! !INTERFACE:
  subroutine scatter_1darray_real (alocal, aglobal, clmlevel)
!
! !DESCRIPTION:
! Wrapper routine to scatter real 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer            :: alocal(:)       ! local  data (output)
    real(r8), pointer            :: aglobal(:)      ! global data (input)
    character(len=*) ,intent(in) :: clmlevel    ! type of input grid
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer            :: n1,n2,lb1,ub1,lb2,ub2 ! indices
    integer            :: lsize      ! size of local array
    type(mct_aVect)    :: AVi, AVo   ! attribute vectors
    real(r8),pointer   :: adata(:)   ! local data array
    character(len=256) :: rstring    ! real field list string
    character(len=256) :: istring    ! int field list string
    character(len=8)   :: fname      ! arbitrary field name
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    character(len=*),parameter :: subname = 'scatter_1darray_real'

!-----------------------------------------------------------------------

    call t_startf(trim(subname)//'_total')
    call get_clmlevel_gsmap(clmlevel,gsmap)

    lb1 = lbound(alocal,dim=1)
    ub1 = ubound(alocal,dim=1)
    lb2 = 1
    ub2 = 1

    rstring = ""
    istring = ""

    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       if (len_trim(rstring) == 0) then
          rstring = trim(fname)
       else
          rstring = trim(rstring)//":"//trim(fname)
       endif
    enddo

    if (masterproc .and. debug > 2) then
       write(iulog,*) trim(subname),' strings:',trim(rstring),' ',trim(istring)
    endif

    if (debug > 1) call t_startf(trim(subname)//'_pack')

    if (masterproc) then
       lsize = size(aglobal,dim=1)
       call mct_aVect_init(AVi,rList=trim(rstring),iList=trim(istring),lsize=lsize)
       allocate(adata(lsize))
       do n2 = lb2,ub2
          adata(1:lsize) = aglobal(1:lsize)
          write(fname,'(a1,i3.3)') 'f',n2-lb2+1
          call mct_aVect_importRattr(AVi,trim(fname),adata,lsize)
       enddo
       deallocate(adata)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_pack')
    if (debug > 1) call t_startf(trim(subname)//'_scat')

    call mct_aVect_scatter(AVi, AVo, gsmap, 0, mpicom)

    if (debug > 1) call t_stopf(trim(subname)//'_scat')
    if (debug > 1) call t_startf(trim(subname)//'_upck')

    lsize = size(alocal,dim=1)
    allocate(adata(lsize))
    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       call mct_aVect_exportRattr(AVo,trim(fname),adata,lsize)
       do n1 = lb1,ub1
          alocal(n1) = adata(n1-lb1+1)
       enddo
    enddo
    deallocate(adata)

    if (debug > 1) call t_stopf(trim(subname)//'_upck')

    if (masterproc) then
       call mct_aVect_clean(AVi)
    endif
    call mct_aVect_clean(AVo)

    call t_stopf(trim(subname)//'_total')

  end subroutine scatter_1darray_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gather_1darray_real
!
! !INTERFACE:
  subroutine gather_1darray_real (alocal, aglobal, clmlevel, missing)
!
! !DESCRIPTION:
! Wrapper routine to gather real 1d array
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), pointer            :: alocal(:)       ! local  data (output)
    real(r8), pointer            :: aglobal(:)      ! global data (input)
    character(len=*) ,intent(in) :: clmlevel    ! type of input grid
    real(r8),optional,intent(in) :: missing     ! missing value
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer            :: n1,n2,lb1,ub1,lb2,ub2 ! indices
    integer            :: lsize      ! size of local array
    type(mct_aVect)    :: AVi, AVo   ! attribute vectors
    real(r8),pointer   :: adata(:)   ! temporary data array
    integer ,pointer   :: mvect(:)   ! local array for mask
    character(len=256) :: rstring    ! real field list string
    character(len=256) :: istring    ! int field list string
    character(len=8)   :: fname      ! arbitrary field name
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    character(len=*),parameter :: subname = 'gather_1darray_real'

!-----------------------------------------------------------------------

    call t_startf(trim(subname)//'_total')
    call get_clmlevel_gsmap(clmlevel,gsmap)

    lsize = size(alocal,dim=1)
    lb1 = lbound(alocal,dim=1)
    ub1 = ubound(alocal,dim=1)
    lb2 = 1
    ub2 = 1
   
    rstring = ""
    istring = ""

    if (present(missing)) then
       istring = "mask"
    endif

    do n2 = lb2,ub2
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       if (len_trim(rstring) == 0) then
          rstring = trim(fname)
       else
          rstring = trim(rstring)//":"//trim(fname)
       endif
    enddo

    if (masterproc .and. debug > 2) then
       write(iulog,*) trim(subname),' strings:',trim(rstring),' ',trim(istring)
    endif

    call mct_aVect_init(AVi,rList=trim(rstring),iList=trim(istring),lsize=lsize)

    if (debug > 1) call t_startf(trim(subname)//'_pack')
    allocate(adata(lsize))
    do n2 = lb2,ub2
       do n1 = lb1,ub1
          adata(n1-lb1+1) = alocal(n1)
       enddo
       write(fname,'(a1,i3.3)') 'f',n2-lb2+1
       call mct_aVect_importRattr(AVi,trim(fname),adata,lsize)
    enddo
    deallocate(adata)

    if (present(missing)) then
       allocate(mvect(lsize))
       do n1 = lb1,ub1
          mvect(n1-lb1+1) = 1
       enddo
       call mct_aVect_importIattr(AVi,"mask",mvect,lsize)
       deallocate(mvect)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_pack')
    if (debug > 1) call t_startf(trim(subname)//'_gath')

    if (present(missing)) then
! tcx wait for update in mct, then get rid of "mask"
!       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom, missing = missing)
       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
    else
       call mct_aVect_gather(AVi, AVo, gsmap, 0, mpicom)
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_gath')
    if (debug > 1) call t_startf(trim(subname)//'_upck')

    if (masterproc) then
       lsize = size(aglobal,dim=1)
       allocate(adata(lsize))
       do n2 = lb2,ub2
          write(fname,'(a1,i3.3)') 'f',n2-lb2+1
          call mct_aVect_exportRattr(AVo,trim(fname),adata,lsize)
          aglobal(1:lsize) = adata(1:lsize)
       enddo
       deallocate(adata)
       if (present(missing)) then
          allocate(mvect(lsize))
          call mct_aVect_exportIattr(AVo,"mask",mvect,lsize)
          do n1 = 1,lsize
             if (mvect(n1) == 0) then
                do n2 = lb2,ub2
                   aglobal(n1) = missing
                enddo
             endif
          enddo
          deallocate(mvect)
       endif
    endif

    if (debug > 1) call t_stopf(trim(subname)//'_upck')

    if (masterproc) then
       call mct_aVect_clean(AVo)
    endif

    call mct_aVect_clean(AVi)

    call t_stopf(trim(subname)//'_total')

  end subroutine gather_1darray_real

end module spmdGathScatMod
