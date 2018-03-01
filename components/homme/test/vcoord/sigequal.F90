      program sigeq
#ifdef NETCDF
      use netcdf_coord_file
#endif
      real*8, parameter :: p0=1000.
      real*8, dimension(:), allocatable :: A,B
      real*8, dimension(:), allocatable :: Amid,Bmid

       integer k
       integer nlev
       character(len=3)  :: charlev
       character(len=80) :: fname

       print *,"input number of vertical levels"
       read(5,*)nlev
       
       allocate(A(nlev+1))
       allocate(B(nlev+1))

       allocate(Amid(nlev))
       allocate(Bmid(nlev))

       do k=1,nlev+1
          A(k) = 0.0D0
          B(k) = REAL(k-1,kind(0.0D0)) / REAL(nlev,kind(0.0D0))
       end do

       do k=1,nlev
          Amid(k) = 0.0D0
          Bmid(k) = 0.50D0*(B(k+1) + B(k))
       end do
       write(charlev,'(i3)') nlev
#ifdef NETCDF
       fname="hyvert-"//TRIM(ADJUSTL(charlev))//".nc"
       call write_netcdf_coord_file(fname,nlev,p0,a,b,amid,bmid)
#else       

       fname="sabi-"//TRIM(ADJUSTL(charlev))//".ascii"
       open(UNIT=11,FILE=fname,form='formatted',status='unknown',&
       iostat=ios)

       fname="sabm-"//TRIM(ADJUSTL(charlev))//".ascii"
       open(UNIT=12,FILE=fname,form='formatted',status='unknown',&
       iostat=ios)

       write(11,'(i4)') nlev+1
       write(11,'(e25.16)') (A(k),k=1,nlev+1)
       write(11,'(i4)') nlev+1
       write(11,'(e25.16)') (B(k),k=1,nlev+1)

       write(12,'(i4)') nlev
       write(12,'(e25.16)') (Amid(k),k=1,nlev)
       write(12,'(i4)') nlev
       write(12,'(e25.16)') (Bmid(k),k=1,nlev)

       close(11)
       close(12)
#endif
       stop
       end
