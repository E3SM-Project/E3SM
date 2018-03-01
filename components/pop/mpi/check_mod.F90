module check_mod
  
  use kinds_mod  
#ifdef USEPIO
  use alloc_mod   ! _EXTERNAL
#endif

  implicit none
  private 

  include 'mpif.h'

  public :: checkpattern 

  interface checkpattern
      module procedure check_1D_r8, &
		       check_1D_r4, &
		       check_1D_i4
      module procedure check_3D_r8, &
		       check_3D_r4, &
		       check_3D_i4
  end interface

contains 

subroutine check_1D_r8(fname,wr_array,rd_array)
 
    character(len=*) :: fname
    real(r8) :: wr_array(:)
    real(r8) :: rd_array(:)

    real(r8),pointer :: diff(:)
    real(r8) :: lsum,gsum
    integer(i4) :: ierr,len,cbad,rank

    len = SIZE(wr_array) 

    allocate(diff(len))

    diff = wr_array - rd_array
    cbad = COUNT(diff .ne. 0.0)
    lsum = SUM(diff)
    
    call MPI_Allreduce(lsum,gsum,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(gsum .ne. 0.0) then
       call MPI_COMM_rank(MPI_COMM_WORLD,rank,ierr)
       write(*,*) 'IAM: ', rank, 'File: ',TRIM(fname),' Error dected for correctness test(R8): ',lsum,' # bad: ',cbad
    endif

end subroutine check_1D_r8

subroutine check_3D_r8(fname,wr_array,rd_array)

    character(len=*) :: fname
    real(r8) :: wr_array(:,:,:)
    real(r8) :: rd_array(:,:,:)

    real(r8), pointer :: diff(:,:,:)
    real(r8) :: lsum,gsum
    integer(i4) :: ierr,cbad,rank
    integer(i4) :: len1,len2,len3
   
    len1 = SIZE(wr_array,dim=1)
    len2 = SIZE(wr_array,dim=2)
    len3 = SIZE(wr_array,dim=3)
   
    allocate(diff(len1,len2,len3))
  
    diff = wr_array - rd_array
    cbad = COUNT(diff .ne. 0.0)
    lsum = SUM(diff)
   
    call MPI_Allreduce(lsum,gsum,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(gsum .ne. 0.0) then
       call MPI_COMM_rank(MPI_COMM_WORLD,rank,ierr)
       write(*,*) 'IAM: ', rank, 'File: ',TRIM(fname),' Error dected for correctness test(3D,R8): ',lsum,' # bad: ',cbad
    endif

end subroutine check_3D_r8

subroutine check_3D_r4(fname,wr_array,rd_array)

    character(len=*) :: fname
    real(r4) :: wr_array(:,:,:)
    real(r4) :: rd_array(:,:,:)

    real(r4), pointer :: diff(:,:,:)
    real(r4) :: lsum,gsum
    integer(i4) :: ierr,cbad,rank
    integer(i4) :: len1,len2,len3
   
    len1 = SIZE(wr_array,dim=1)
    len2 = SIZE(wr_array,dim=2)
    len3 = SIZE(wr_array,dim=3)
   
    allocate(diff(len1,len2,len3))
  
    diff = wr_array - rd_array
    cbad = COUNT(diff .ne. 0.0)
    lsum = SUM(diff)
   
    call MPI_Allreduce(lsum,gsum,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(gsum .ne. 0.0) then
       call MPI_COMM_rank(MPI_COMM_WORLD,rank,ierr)
       write(*,*) 'IAM: ', rank, 'File: ',TRIM(fname),' Error dected for correctness test(3D,R4): ',lsum,' # bad: ',cbad
    endif

end subroutine check_3D_r4

subroutine check_3D_i4(fname,wr_array,rd_array)

    character(len=*) :: fname
    integer(i4) :: wr_array(:,:,:)
    integer(i4) :: rd_array(:,:,:)

    integer(i4), pointer :: diff(:,:,:)
    integer(i4) :: lsum,gsum
    integer(i4) :: ierr,cbad,rank
    integer(i4) :: len1,len2,len3
   
    len1 = SIZE(wr_array,dim=1)
    len2 = SIZE(wr_array,dim=2)
    len3 = SIZE(wr_array,dim=3)
   
    allocate(diff(len1,len2,len3))
  
    diff = wr_array - rd_array
    cbad = COUNT(diff .ne. 0.0)
    lsum = SUM(diff)
   
    call MPI_Allreduce(lsum,gsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(gsum .ne. 0.0) then
       call MPI_COMM_rank(MPI_COMM_WORLD,rank,ierr)
       write(*,*) 'IAM: ', rank, 'File: ',TRIM(fname),' Error dected for correctness test(3D,I4): ',lsum,' # bad: ',cbad
    endif

end subroutine check_3D_i4

subroutine check_1D_r4(fname,wr_array,rd_array)
 
    character(len=*) :: fname
    real(r4) :: wr_array(:)
    real(r4) :: rd_array(:)

    real(r4),pointer :: diff(:)
    real(r4) :: lsum,gsum
    integer(i4) :: ierr,len,cbad,rank

    len = SIZE(wr_array) 

    allocate(diff(len))

    diff = wr_array - rd_array
    cbad = COUNT(diff .ne. 0.0)
    lsum = SUM(diff)
    
    call MPI_Allreduce(lsum,gsum,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(gsum .ne. 0.0) then
       call MPI_COMM_rank(MPI_COMM_WORLD,rank,ierr)
       write(*,*) 'IAM: ', rank, 'File: ',TRIM(fname),' Error dected for correctness test(R4): ',lsum,' # bad: ',cbad
    endif

end subroutine check_1D_r4

subroutine check_1D_i4(fname,wr_array,rd_array)
 
    character(len=*) :: fname
    integer(i4) :: wr_array(:)
    integer(i4) :: rd_array(:)

    integer(i4),pointer :: diff(:)
    integer(i4) :: lsum,gsum
    integer(i4) :: ierr,len,cbad,rank

    len = SIZE(wr_array) 

    allocate(diff(len))

    diff = wr_array - rd_array
    cbad = COUNT(diff .ne. 0.0)
    lsum = SUM(diff)
    
    call MPI_Allreduce(lsum,gsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(gsum .ne. 0.0) then
       call MPI_COMM_rank(MPI_COMM_WORLD,rank,ierr)
       write(*,*) 'IAM: ', rank, 'File: ',TRIM(fname),' Error dected for correctness test(I4): ',lsum,' # bad: ',cbad
    endif

end subroutine check_1D_i4

end module check_mod
