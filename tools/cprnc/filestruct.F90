module filestruct
  use netcdf
  implicit none
  type dim_t
     integer :: dimsize
     integer :: start, kount  ! used for user requested dimension subsetting
     character(len=nf90_MAX_NAME) ::name = ''
  end type dim_t

  type var_t
     integer :: matchid
     integer :: ndims
     integer :: natts
     integer, pointer :: dimids(:)
     integer :: xtype
     character(len=nf90_MAX_NAME) ::name = ''
  end type var_t

  type file_t
     integer :: fh
     integer :: natts
     type(dim_t), pointer :: dim(:)
     type(var_t), pointer :: var(:)
     integer :: unlimdimid
  end type file_t

  logical :: verbose

contains
  subroutine init_file_struct( file, dimoptions )

    type(file_t) :: file
    type(dim_t), optional :: dimoptions(:)
    integer :: ndims, nvars
    integer :: dimids(NF90_MAX_DIMS)
    integer :: i, ierr, docnt, n1, n2
    integer :: j, start, kount
    character(len=NF90_MAX_NAME) :: name, dname
    ierr= nf90_inquire(file%fh, ndims, nvars, file%natts, file%unlimdimid)

    allocate(file%dim(ndims))
    allocate(file%var(nvars))


    do i=1,ndims
       ierr = nf90_inquire_dimension(file%fh, i, file%dim(i)%name, file%dim(i)%dimsize)
       file%dim(i)%start=1
       if(i==file%unlimdimid) then
          file%dim(i)%kount=1
       else
          file%dim(i)%kount=file%dim(i)%dimsize
       end if
    end do

    if(present(dimoptions)) then
       docnt = size(dimoptions)
       do j=1,docnt
          start = dimoptions(j)%start
          kount = dimoptions(j)%kount
          name = dimoptions(j)%name
          n1 = len_trim(name)
          do i=1,ndims
             dname = file%dim(i)%name
             n2 = len_trim(dname)
             if(name(1:n1).eq.dname(1:n2) ) then


                if((start > 0) .and. (start < file%dim(i)%dimsize)) then
                   file%dim(i)%start = start
                else
                   write(6,*) 'Command line start value for dim ',name(1:n1),&
                        ' out of bounds, expected 1-',file%dim(i)%dimsize,' got: ',start
                   stop
                end if
                if(kount > 0 .and. start+kount <= file%dim(i)%dimsize) then
                   file%dim(i)%kount = kount
                else if(kount == -1) then
                   file%dim(i)%kount = file%dim(i)%dimsize-file%dim(i)%start+1
                else
                   write(6,*) 'Command line count value for dim ',name(1:n1),&
                        ' out of bounds, expected 1-',file%dim(i)%dimsize-file%dim(i)%start+1,' got: ',kount
                   stop

                endif
                write(6,*) 'Setting dimension bounds for dim ',name(1:n1),file%dim(i)%start,file%dim(i)%kount

                exit
             end if
          end do
       end do
    end if

    do i=1,nvars
       file%var(i)%matchid=-1
       ierr = nf90_inquire_variable(file%fh, i, file%var(i)%name, file%var(i)%xtype, file%var(i)%ndims, dimids, &
            file%var(i)%natts)
       allocate(file%var(i)%dimids(file%var(i)%ndims))
       file%var(i)%dimids = dimids(1:file%var(i)%ndims)
    end do


  end subroutine init_file_struct


  subroutine compare_metadata(file1, file2, vid)
    type(file_t) :: file1, file2
    integer, optional, intent(in) :: vid

    integer :: id1, id2, natts1, natts2

    integer :: i, ierr
    character(len=NF90_MAX_NAME) :: attname
    integer :: atttype, attlen

    real, pointer :: attreal1(:), attreal2(:)
    double precision, pointer :: attdouble1(:),attdouble2(:)
    integer, pointer :: attint1(:),attint2(:)
    character(len=8192) :: attchar1, attchar2
    logical :: found


    if(present(vid)) then
       id1 = vid
       id2 = file1%var(id1)%matchid
       ierr = nf90_inquire_variable(file1%fh, id1, nAtts=natts1)
       ierr = nf90_inquire_variable(file2%fh, id2, nAtts=natts2)
    else
       id1 = NF90_GLOBAL
       id2 = NF90_GLOBAL
       natts1 = file1%natts
       natts2 = file2%natts
    end if

    do i=1,natts1
       found = .true.
       attname = ''
       ierr = nf90_inq_attname(file1%fh, id1, i, attname)
       ierr = nf90_inquire_attribute(file1%fh, id1, trim(attname), atttype, attlen)
       select case(atttype)
       case(nf90_char)
          attchar1=' '
          attchar2=' '

          ierr = nf90_get_att(file1%fh,id1, trim(attname), attchar1)
          ierr = nf90_get_att(file2%fh,id2, trim(attname), attchar2)
          if(ierr==NF90_NOERR) then
             if(trim(attname).ne.'case' .and. attchar1(1:attlen) .ne. attchar2(1:attlen)) then
                print *, 'Attribute ',trim(attname),' from file1: ',attchar1(1:attlen),&
                     ' does not match that found on file2: ',attchar2(1:attlen)
             end if
          else
             print *, 'Attribute ',trim(attname),' from file1: ',attchar1(1:attlen),&
                  ' not found on file2'
          end if
          if(id1==NF90_GLOBAL .and. trim(attname) .eq. 'case') then
             print *, 'CASE 1 : ',trim(attchar1)
             print *, 'CASE 2 : ',trim(attchar2)
          endif
          if(id1==NF90_GLOBAL .and. trim(attname) .eq. 'title') then
             print *, 'TITLE 1 : ',trim(attchar1)
             print *, 'TITLE 2 : ',trim(attchar2)
          end if
       case(nf90_int)
          allocate(attint1(attlen),attint2(attlen))
          ierr = nf90_get_att(file1%fh,id1, trim(attname), attint1)
          ierr = nf90_get_att(file2%fh,id2, trim(attname), attint2)

          if(ierr==NF90_NOERR) then
             if(any(attint1 /= attint2)) then
                print *, 'Attribute ',trim(attname),' from file1: ',attint1,' does not match that found on file2 ',attint2
             end if
          else
             print *, 'Attribute ',trim(attname),' from file1: ',attint1,' not found on file2'
          end if
          deallocate(attint1, attint2)


       case(nf90_float)
          allocate(attreal1(attlen),attreal2(attlen))
          ierr = nf90_get_att(file1%fh,id1, trim(attname), attreal1)
          ierr = nf90_get_att(file2%fh,id2, trim(attname), attreal2)
          if(ierr==NF90_NOERR) then
             if(any(attreal1 /= attreal2)) then
                print *, 'Attribute ',trim(attname),' from file1: ',attreal1,' does not match that found on file2 ',attreal2
             end if
          else
             print *, 'Attribute ',trim(attname),' from file1: ',attreal1,' not found on file2'
          end if
          deallocate(attreal1, attreal2)
       case(nf90_double)
          allocate(attdouble1(attlen), attdouble2(attlen))
          ierr = nf90_get_att(file1%fh,id1, trim(attname), attdouble1)
          ierr = nf90_get_att(file2%fh,id2, trim(attname), attdouble2)
          if(ierr==NF90_NOERR) then
             if(any(attdouble1 /= attdouble2)) then
                print *, 'Attribute ',trim(attname),' from file1: ',attdouble1,' does not match that found on file2 ',attdouble2
             end if
          else
             print *, 'Attribute ',trim(attname),' from file1: ',attdouble1,' not found on file2'
          end if
          deallocate(attdouble1, attdouble2)
       case default
          print *,' Did not recognize attribute type ',atttype, trim(attname), attlen
       end select
    end do

  end subroutine compare_metadata








  subroutine compare_dimensions( dimfile1, dimfile2)
    type(dim_t), intent(in) :: dimfile1(:), dimfile2(:)

    integer :: ds1, ds2
    integer :: i, j
    logical,pointer :: found(:,:)

    ds1 = size(dimfile1)
    ds2 = size(dimfile2)

    allocate(found(2,max(ds1,ds2)))

    found = .false.
    do i=1,ds1
       do j=1,ds2
          if(dimfile1(i)%name .eq. dimfile2(j)%name) then
             if(dimfile1(i)%dimsize == dimfile2(j)%dimsize) then
                print *, 'Dimension ',trim(dimfile1(i)%name), ' matches'
             else
                print *, 'Dimension ',trim(dimfile1(i)%name), ' differs ', dimfile1(i)%dimsize, ' /= ',dimfile2(j)%dimsize
             end if
             found(1,i) = .true.
             found(2,j) = .true.
          end if
       end do
    end do
    do i=1,ds1
       if(.not. found(1,i)) then
          print *, 'Could not find match for file 1 dimension ',trim(dimfile1(i)%name)
       end if
    end do
    do i=1,ds2
       if(.not. found(2,i)) then
          print *, 'Could not find match for file 2 dimension ',trim(dimfile2(i)%name)
       end if
    end do
    deallocate(found)
  end subroutine compare_dimensions


  subroutine match_vars( file1, file2 )
    type(file_t), intent(inout) :: file1, file2

    type(var_t), pointer :: varfile1(:),varfile2(:)

    integer :: vs1, vs2, i, j



    varfile1 => file1%var
    varfile2 => file2%var

    vs1 = size(varfile1)
    vs2 = size(varfile2)

    do i=1,vs1
       do j=1,vs2
          if(varfile1(i)%name .eq. varfile2(j)%name) then
             varfile1(i)%matchid=j
             varfile2(j)%matchid=i
          end if
       end do
    end do
    do i=1,vs1
       if(varfile1(i)%matchid<0) then
          print *, 'Could not find match for file1 variable ',trim(varfile1(i)%name), ' in file2'
       end if
    end do
    do i=1,vs2
       if(varfile2(i)%matchid<0) then
          print *, 'Could not find match for file2 variable ',trim(varfile2(i)%name), ' in file1'
       end if
    end do
  end subroutine match_vars


  function vdimsize(dims, dimids)
    type(dim_t), intent(in) :: dims(:)
    integer, intent(in) :: dimids(:)

    integer :: vdimsize
    integer :: i

    vdimsize=1
    do i=1,size(dimids)
       if(verbose) print *,__FILE__,__LINE__,i,dimids(i),size(dims),size(dimids)
       vdimsize = vdimsize*dims(dimids(i))%kount
    end do

  end function vdimsize





  subroutine compare_var_int(f1, f2, i1, i2, t)
    type(file_t) :: f1,f2
    integer, intent(in) :: i1, i2
    integer, optional :: t


    integer :: s1, s2, m1, m2, l1(1), l2(1), i, ierr
    integer, pointer :: v1(:), v2(:), vdiff(:)
    integer :: t1, n1
    integer :: start(NF90_MAX_DIMS), count(NF90_MAX_DIMS)

    if(present(t)) then
       t1 = t
    else
       t1 = 1
    end if

    s1 = vdimsize(f1%dim, f1%var(i1)%dimids)
    s2 = vdimsize(f2%dim, f2%var(i2)%dimids)

    if(s1 /= s2) then
       print *, 'Variable ',f1%var(i)%name,' sizes differ'
    end if

    n1 = size(f1%var(i1)%dimids)
    start = 1
    do i=1,n1
       count(i) = f1%dim(f1%var(i1)%dimids(i))%dimsize
       if(f1%var(i1)%dimids(i) == f1%unlimdimid) then
          count(i)=1
          start(i)=t1
       end if
    end do

    allocate(v1(s1), v2(s2))

    ierr = nf90_get_var(f1%fh, i1, v1, start(1:n1), count(1:n1))
    ierr = nf90_get_var(f2%fh, i2, v2, start(1:n1), count(1:n1))

    if(any(v1 /= v2)) then
       allocate(vdiff(s1))
       vdiff = abs(v1-v2)
       m1 = maxval(vdiff)
       m2 = minval(vdiff)
       l1 = maxloc(vdiff)
       l2 = minloc(vdiff)

       print *,__FILE__,__LINE__,m1,m2,l1,l2
       deallocate(vdiff)
    end if

    deallocate(v1,v2)
  end subroutine compare_var_int

  subroutine compare_var_float(f1, f2, i1, i2, t)
    type(file_t) :: f1,f2
    integer, intent(in) :: i1, i2
    integer, optional :: t


    integer :: s1, s2, m1, m2, l1(1), l2(1), i, ierr
    real, pointer :: v1(:), v2(:), vdiff(:)
    integer :: t1, n1
    integer :: start(NF90_MAX_DIMS), count(NF90_MAX_DIMS)

    if(present(t)) then
       t1 = t
    else
       t1 = 1
    end if

    s1 = vdimsize(f1%dim, f1%var(i1)%dimids)
    s2 = vdimsize(f2%dim, f2%var(i2)%dimids)

    if(s1 /= s2) then
       print *, 'Variable ',f1%var(i)%name,' sizes differ'
    end if

    n1 = size(f1%var(i1)%dimids)
    start = 1
    do i=1,n1
       count(i) =  f1%dim(f1%var(i1)%dimids(i))%dimsize
       if(f1%var(i1)%dimids(i) == f1%unlimdimid) then
          count(i)=1
          start(i)=t1
       end if
    end do

    allocate(v1(s1), v2(s2))

    ierr = nf90_get_var(f1%fh, i1, v1, start(1:n1), count(1:n1))
    ierr = nf90_get_var(f2%fh, i2, v2, start(1:n1), count(1:n1))

    if(any(v1 /= v2)) then
       allocate(vdiff(s1))
       vdiff = abs(v1-v2)
       m1 = maxval(vdiff)
       m2 = minval(vdiff)
       l1 = maxloc(vdiff)
       l2 = minloc(vdiff)

       print *,__FILE__,__LINE__,m1,m2,l1,l2
       deallocate(vdiff)
    end if

    deallocate(v1,v2)
  end subroutine compare_var_float

  subroutine compare_var_double(f1, f2, i1, i2, t)
    type(file_t) :: f1,f2
    integer, intent(in) :: i1, i2
    integer, optional :: t


    integer :: s1, s2, m1, m2, l1(1), l2(1), i, ierr
    double precision, pointer :: v1(:), v2(:), vdiff(:)
    integer :: t1, n1
    integer :: start(NF90_MAX_DIMS), count(NF90_MAX_DIMS)

    if(present(t)) then
       t1 = t
    else
       t1 = 1
    end if

    s1 = vdimsize(f1%dim, f1%var(i1)%dimids)
    s2 = vdimsize(f2%dim, f2%var(i2)%dimids)

    if(s1 /= s2) then
       print *, 'Variable ',f1%var(i)%name,' sizes differ'
    end if

    n1 = size(f1%var(i1)%dimids)
    start = 1
    do i=1,n1
       count(i) =  f1%dim(f1%var(i1)%dimids(i))%dimsize
       if(f1%var(i1)%dimids(i) == f1%unlimdimid) then
          count(i)=1
          start(i)=t1
       end if
    end do

    allocate(v1(s1), v2(s2))

    ierr = nf90_get_var(f1%fh, i1, v1, start(1:n1), count(1:n1))
    ierr = nf90_get_var(f2%fh, i2, v2, start(1:n1), count(1:n1))

    if(any(v1 /= v2)) then
       allocate(vdiff(s1))
       vdiff = abs(v1-v2)
       m1 = maxval(vdiff)
       m2 = minval(vdiff)
       l1 = maxloc(vdiff)
       l2 = minloc(vdiff)

       print *,__FILE__,__LINE__,m1,m2,l1,l2
       deallocate(vdiff)
    end if

    deallocate(v1,v2)
  end subroutine compare_var_double









end module filestruct
