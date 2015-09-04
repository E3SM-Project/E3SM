module utils
  use filestruct, only : dim_t
contains

subroutine get_dimname_str(ndims,dimids,dims,dimname_str)
  integer, intent(in) :: ndims
  integer, intent(in) :: dimids(:)
  type(dim_t) :: dims(:)
  character(len=*),intent(out) :: dimname_str

  integer :: dlen
  integer :: j

  dimname_str = ' ' 
  
  if(ndims>0) then
     dimname_str(1:1) = '('
     dlen=2
     
     do j=1,ndims
        dimname_str(dlen:) = trim(dims(dimids(j))%name)//','
        dlen=dlen+ len_trim(dims(dimids(j))%name) + 1
     end do
     dimname_str(dlen-1:dlen-1) = ')'
  end if


end subroutine get_dimname_str

subroutine get_dim_str(ndims,loc,dim_str)
  integer, intent(in) :: ndims
  integer, intent(in) :: loc(:)
  character(len=*),intent(out) :: dim_str

  integer :: dlen
  integer :: j

  dim_str = ' ' 
  
  if(ndims>0) then
     dim_str(1:1) = '('
     dlen=2
     
     do j=1,ndims
        write(dim_str(dlen:),'(i6,a)') loc(j),','

        dlen=len_trim(dim_str)+1
     end do
     dim_str(dlen-1:dlen-1) = ')'
  end if


end subroutine get_dim_str



subroutine checknf90(ierr,returnflag)
  use netcdf, only : nf90_noerr, nf90_strerror
  integer, intent(in) :: ierr
  logical, optional, intent(in) :: returnflag

  if(ierr/=NF90_NOERR) then
     print *, trim(nf90_strerror(ierr))
     if(present(returnflag)) then
        if(returnflag) return
     end if
#ifdef AIX
     call xl__trbk()
#endif
     stop

  end if



end subroutine checknf90
  



end module utils
