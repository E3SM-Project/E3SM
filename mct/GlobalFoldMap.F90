!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalFoldMap - define a fold-over-the-equator
!  mapping of our grid onto the cartesian processor topology: 
!  each row of processors handles points from
!  the northern and southern hemispheres.
!
! !DESCRIPTION:
!
! !INTERFACE:

      subroutine GlobalFoldMap(GSMap, nx, ny, &
                               nxprocs, nyprocs, &
                               root, my_comm, comp_id)

!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio
      use m_GlobalSegMap,only: MCT_GSMap_init => init
      use m_GlobalSegMap,only: GlobalSegMap

      implicit none


      type(GlobalSegMap),intent(out)  :: GSMap ! Output GlobalSegMap

      integer,intent(in)              :: nx    ! number of longitudes on grid
      integer,intent(in)              :: ny    ! number of latitudes on grid
      integer,intent(in)              :: nxprocs ! # procs along lon axis   
      integer,intent(in)              :: nyprocs ! # procs along lat axis  
      integer,intent(in)              :: root    ! root on my_com
      integer,intent(in)              :: my_comm ! local communicatior
      integer,intent(in)              :: comp_id ! component model ID


! !REVISION HISTORY:
!      13Mar01 - Everest Ong <eong@mcs.anl.gov> - First working code
!EOP ___________________________________________________________________


!----------------------- MCT  model vars
      integer :: plat,plon,row,col,i,j,n,k
      integer :: ier,ierr
      integer :: myProc,mySize 

      integer,dimension(:),pointer :: starts
      integer,dimension(:),pointer :: lengths
      integer,dimension(:,:),pointer :: myglobalmap


!----------------------- MPI Checks
      call MPI_COMM_RANK (my_comm, myProc, ierr)
!     write(*,*) myProc, ' in ccm'


! check to see if there are enough processors
      call MPI_COMM_SIZE(my_comm, mySize, ierr)
      if (mySize /= nyprocs*nxprocs) then
	  write(*,*)'ERROR: found in subroutine GlobalFoldMap'
	  write(*,*)'ERROR: wrong number of processors'
	  write(*,*)'found ',mySize,' Needed',nyprocs*nxprocs
	  stop
      endif

! set local latitude and longitude size
      plat = ny / nyprocs   
      plon = nx / nxprocs

! define a Cartesian topology by assigning
! row and column indicies to each processor.
! processor with rank 0 is (0,0)
      row = myProc / nxprocs
      col = mod(myProc,nxprocs)
      

!  First, we number the grid 1 to Nax*Nay, starting
! in the South Pole and proceeding along a latitude and
! then from south to north.

      allocate(myglobalmap(nx,ny),stat=ierr)
      n=0
      do j=1,ny
       do i= 1,nx
	 n=n+1
	 myglobalmap(i,j) = n
       enddo
      enddo
!
!  For each processor, each seglength is plon
!
! the value of the global index at the start of each
! segment can be found from myglobalmap

      allocate(starts(plat),stat=ierr)
      allocate(lengths(plat),stat=ierr)

! the fist plat/2 latitudes are from the southern hemisphere
      do j=1,plat/2
       starts(j)= myglobalmap(col*plon+1, &
	(plat/2 * row) + j)
       lengths(j)=plon
      enddo

! the next plat/2 latitudes are from the northern hemisphere
      n=1
      do j=plat/2 + 1,plat
       starts(j)=myglobalmap(col*plon + &
            1,(ny - (plat/2 * (row+1))) + n)
       lengths(j)=plon
       n=n+1
      enddo

! now put all this information in a GlobalSegMap.

      call MCT_GSMap_init(GSMap,starts,lengths,root, &
			my_comm,comp_id)

      deallocate(starts,lengths,myglobalmap,stat=ierr)


      end subroutine GlobalFoldMap









