!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module io_pio

!BOP
! !MODULE: io_pio 
! !DESCRIPTION:
!  Interfaces for pio initialization
!
! !REVISION HISTORY:
! SVN:$ID: 
!

! !USES:

  use kinds_mod
  use broadcast
  use communicate
  use exit_mod
  use POP_IOUnitsMod
  use io_types
  use pio
  use shr_pio_mod,       only: shr_pio_getiosys, shr_pio_getiotype

  implicit none
  private
  save

!PUBLIC MEMBER FUNCTIONS:

  public ::  io_pio_init
  public ::  io_pio_initdecomp  

!PUBLIC DATA MEMBERS

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

  integer (i4), parameter :: nmax = 10

  type ptr_ioDesc_int_type
     type (IO_desc_t), pointer :: ioDesc(:)
  end type ptr_ioDesc_int_type
  type (ptr_ioDesc_int_type), dimension(:), allocatable :: ptr_ioDesc_i

  type ptr_ioDesc_real_type
     type (IO_desc_t), pointer :: ioDesc(:)
  end type ptr_ioDesc_real_type
  type (ptr_ioDesc_real_type), dimension(:), allocatable :: ptr_ioDesc_r

  type ptr_ioDesc_double_type
     type (IO_desc_t), pointer :: ioDesc(:)
  end type ptr_ioDesc_double_type
  type (ptr_ioDesc_double_type), dimension(:), allocatable :: ptr_ioDesc_d

  integer(i4), parameter :: iunset = -999
  integer(i4), dimension(nmax) :: nsize3d_i = iunset	
  integer(i4), dimension(nmax) :: nsize3d_r = iunset	
  integer(i4), dimension(nmax) :: nsize3d_d = iunset	

  integer(i4), dimension(nmax) :: ksize3d_i = iunset	
  integer(i4), dimension(nmax) :: ksize3d_r = iunset	
  integer(i4), dimension(nmax) :: ksize3d_d = iunset	

!EOC
!***********************************************************************

 contains

!***********************************************************************
!EOP
! !IROUTINE: io_pio_init - initialize io for input or output
! !INTERFACE: 
   subroutine io_pio_init(mode, filename, File, clobber, cdf64)
     use pio
!
! !DESCRIPTION:
!    Read the pio_inparm namelist and initialize the io subsystem
!
! !REVISION HISTORY:
!    2009-Feb-17 - J. Edwards - initial version
!
! !INPUT/OUTPUT PARAMETERS:
!
   implicit none
   character(len=*)     , intent(in)    :: mode
   character(len=*)     , intent(in)    :: filename
   type(file_desc_t)    , intent(inout) :: File
   logical,optional     , intent(in)    :: clobber
   logical,optional     , intent(in)    :: cdf64
!
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type(iosystem_desc_t), pointer :: io_pio_subsystem
   logical :: exists
   logical :: lclobber
   logical :: lcdf64
   integer :: status
   integer :: nmode
   integer :: pio_iotype
   character(*),parameter :: subName = '(io_pio_init) '


!-----------------------------------------------------------------------
!
!  read and define namelist inputs
!
!-----------------------------------------------------------------------

   io_pio_subsystem => shr_pio_getiosys(inst_name)
   pio_iotype =  shr_pio_getiotype(inst_name)

   if (trim(mode) == 'write') then
      lclobber = .false.
      if (present(clobber)) lclobber=clobber
   
      lcdf64 = .false.
      if (present(cdf64)) lcdf64=cdf64

      if (.not. pio_file_is_open(file)) then
         ! filename not open
         if (my_task == master_task) then
            inquire(file=trim(filename),exist=exists)
         end if
         call broadcast_scalar(exists, master_task)

         if (exists) then
            if (lclobber) then
               nmode = pio_clobber
               if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
               status = pio_createfile(io_pio_subsystem, File, pio_iotype, trim(filename), nmode)
               if (my_task == master_task) then
                  write(stdout,*) subname,' create file ',trim(filename)
               end if
            else
               status = pio_openfile(io_pio_subsystem, File, pio_iotype, trim(filename), pio_write)
               if (my_task == master_task) then
                  write(stdout,*) subname,' open file ',trim(filename)
               end if
            endif
         else
            nmode = pio_noclobber
            if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
            status = pio_createfile(io_pio_subsystem, File, pio_iotype, trim(filename), nmode)
            if (my_task == master_task) then
               write(stdout,*) subname,' create file ',trim(filename)
            end if
         endif
      else
         ! filename is already open, just return
      endif
   end if

   if (trim(mode) == 'read') then
      if (my_task == master_task) then
         inquire(file=trim(filename),exist=exists)
      end if
      call broadcast_scalar(exists, master_task)
      if (exists) then
         status = pio_openfile(io_pio_subsystem, File, pio_iotype, trim(filename), pio_nowrite)
      else
         if(my_task==master_task) then
            write(stdout,*) 'io_pio_ropen ERROR: file invalid ',trim(filename)
         end if
         call exit_POP(sigAbort, 'aborting in io_pio_ropen with invalid file')
      endif
   end if
      
   if (my_task==master_task) then
      call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
   end if

  end subroutine io_pio_init

!================================================================================

   subroutine io_pio_initdecomp (basetype, ndim3, kdim3, iodesc)

      use blocks, only : block, nx_block, ny_block, get_block
      use domain, only : nblocks_clinic, blocks_clinic
      use POP_DomainSizeMod, only : POP_nxGlobal, POP_nyGlobal  

      integer (i4)          , intent(in) :: basetype
      integer(kind=int_kind), intent(in) :: ndim3
      integer(kind=int_kind), intent(in) :: kdim3
      type(io_desc_t)       , pointer    :: iodesc

      integer (kind=int_kind) :: &
          iblk,ib,ie,jb,je,lon,lat,i,j,n,k,index
      
      type(iosystem_desc_t), pointer :: io_pio_subsystem
      type(block) :: this_block 

      integer(kind=int_kind), pointer :: dof3d(:)

      logical, save :: first_time = .true.

      logical :: set_iodesc

      if (first_time) then
         allocate(ptr_ioDesc_i(nmax))
         allocate(ptr_ioDesc_r(nmax))
         allocate(ptr_ioDesc_d(nmax))
         do i = 1,nmax
            allocate(ptr_ioDesc_i(i)%ioDesc(1))
            allocate(ptr_ioDesc_r(i)%ioDesc(1))
            allocate(ptr_ioDesc_d(i)%ioDesc(1))
         end do
         first_time = .false.
      end if
      
      if (basetype == PIO_INT) then
         do i = 1,nmax
            if (nsize3d_i(i) == ndim3 .and. ksize3d_i(i) == kdim3) then
               index = i
               set_ioDesc = .false.
               exit
            else if (nsize3d_i(i) == iunset .and. ksize3d_i(i) == iunset) then
               index = i
               nsize3d_i(index) = ndim3 
	       ksize3d_i(index) = kdim3
               set_ioDesc = .true.
               exit
            end if
         end do
      else if (basetype == PIO_REAL) then
         do i = 1,nmax
            if (nsize3d_r(i) == ndim3 .and. ksize3d_r(i) == kdim3) then
               index = i
               set_ioDesc = .false.
               exit
            else if (nsize3d_r(i) == iunset .and. ksize3d_r(i) == iunset) then
               index = i
               nsize3d_r(index) = ndim3 
	       ksize3d_r(index) = kdim3
               set_ioDesc = .true.
               exit
            end if
         end do
      else if (basetype == PIO_DOUBLE) then
         do i = 1,nmax  
            if (nsize3d_d(i) == ndim3 .and. ksize3d_d(i) == kdim3) then
               index = i
               set_ioDesc = .false.
               exit
            else if (nsize3d_d(i) == iunset .and. ksize3d_d(i) == iunset) then
               index = i
               nsize3d_d(index) = ndim3 
	       ksize3d_d(index) = kdim3
               set_ioDesc = .true.
               exit
            end if
         end do
      end if

      if (set_ioDesc) then

	 if ((ndim3 == 0 .and. kdim3 /= 0) .or. (ndim3 /=0 .and. kdim3 == 0)) then
            call exit_POP(sigAbort,' io_pio_initdecomp: ndim3 and kdim3 must both be zero or nonzero')
         end if

	 if (ndim3 > kdim3) then
            call exit_POP(sigAbort,' io_pio_initdecomp: ndim3 must be less than or equal to kdim3')
         end if

         if (ndim3 == 0) then
            allocate(dof3d(nx_block*ny_block*nblocks_clinic))
            n=0
            do iblk = 1, nblocks_clinic
               this_block = get_block(blocks_clinic(iblk),iblk)         
               ib = this_block%ib
               ie = this_block%ie
               jb = this_block%jb
               je = this_block%je
               
               do j=1,ny_block
               do i=1,nx_block  
                  n = n+1
                  if (j < jb .or. j>je) then
                     dof3d(n)=0
                  else if (i < ib .or. i > ie) then
                     dof3d(n) = 0
                  else
                     lon = this_block%i_glob(i)
                     lat = this_block%j_glob(j)
                     dof3d(n) = ((lat-1)*POP_nxGlobal + lon)
                  endif
               enddo !i
               enddo !j
            end do
         else
            allocate(dof3d(nx_block*ny_block*nblocks_clinic*kdim3))
            n=0
            do iblk = 1, nblocks_clinic
               this_block = get_block(blocks_clinic(iblk),iblk)         
               ib = this_block%ib
               ie = this_block%ie
               jb = this_block%jb
               je = this_block%je
               
               do k=1,kdim3
               do j=1,ny_block
               do i=1,nx_block  
                  n = n+1
                  if (j < jb .or. j>je) then
                     dof3d(n)=0
                  else if (i < ib .or. i > ie) then
                     dof3d(n) = 0
                  else
                     if (k > ndim3) then
                        dof3d(n) = 0
                     else
                        lon = this_block%i_glob(i)
                        lat = this_block%j_glob(j)
                        dof3d(n) = ((lat-1)*POP_nxGlobal + lon) + (k-1)*POP_nxGlobal*POP_nyGlobal 
                     end if
                  endif
               enddo !i
               enddo !j
               enddo !kdim3
            end do
         end if

         io_pio_subsystem => shr_pio_getiosys(inst_name)
         if (basetype == PIO_INT) then
            if (ndim3 == 0) then
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal/), &
                    dof3d, ptr_ioDesc_i(index)%ioDesc(1))
            else
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal,ndim3/), &
                    dof3d, ptr_ioDesc_i(index)%ioDesc(1))
            end if
         else if (basetype == PIO_REAL) then
            if (ndim3 == 0) then
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal/), &
                    dof3d, ptr_ioDesc_r(index)%ioDesc(1))
            else
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal,ndim3/), &
                    dof3d, ptr_ioDesc_r(index)%ioDesc(1))
            end if
         else if (basetype == PIO_DOUBLE) then
            if (ndim3 == 0) then
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal/), &
                    dof3d, ptr_ioDesc_d(index)%ioDesc(1))
            else
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal,ndim3/), &
                    dof3d, ptr_ioDesc_d(index)%ioDesc(1))
            end if
         end if

         deallocate(dof3d)
      end if

      if (basetype == PIO_INT) then
         iodesc => ptr_ioDesc_i(index)%ioDesc(1)
      elseif (basetype == PIO_REAL) then
         iodesc => ptr_ioDesc_r(index)%ioDesc(1)
      elseif (basetype == PIO_DOUBLE) then
         iodesc => ptr_ioDesc_d(index)%ioDesc(1)
      end if

    end subroutine io_pio_initdecomp

!================================================================================

end module io_pio      
