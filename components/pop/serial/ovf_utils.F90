!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ovf_utils

!SERIAL VERSION

!BOP
! !MODULE: ovf_util
! !DESCRIPTION:
!  This module contains support for summing up the contributions
!   from the following ovf regions:  
!   inflow 
!   source 
!   source adjacent 
!   entrainment 
!   entrainment adjacent 
!   product adjacent 

! !REVISION HISTORY:

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod

   use kinds_mod, only: r4, r8, i4, int_kind, char_len, log_kind, rtavg

   use constants

   use blocks, only: nx_block, ny_block, block, get_block, &
        nblocks_x, nblocks_y, nblocks_tot, get_block_ids_from_coords, &
        get_block_parameter 
   
   use global_reductions, only: global_sum

   use grid, only: DXT, DYT, dz

   use prognostic, only: RHO, TRACER

   use mpi2s_gshalo

   use communicate, only: my_task, mpi_dbl, master_task

   use pop_quicksort, only: pop_quick_sort

   use distribution

   use domain, only: nblocks_clinic, blocks_clinic, distrb_clinic, &
        ltripole_grid
   use domain_size, only: max_blocks_clinic, nt , nx_global 

   use overflow_type, only: ovf, ovf_groups, num_ovf  

   implicit none
   private



   save


! !LOCAL STUFF   
   integer (int_kind), parameter :: length_comm = nt +1

   integer(i4), dimension(:), allocatable :: ovf_print_myRequests(:)
   integer (i4), dimension(:,:), allocatable:: ovf_print_myStati(:,:)


! !PUBLIC MEMBER FUNCTIONS:
   public ::  ovf_utils_avgs   !overflows.F90
   public ::  ovf_init_groups
   public ::  ovf_print_init
   public ::  ovf_print_send
   public ::  ovf_print_get
   public ::  ovf_print_finalize
! !PUBLIC DATA MEMBERS:


!EOP


!***********************************************************************

 contains

!***********************************************************************


! !IROUTINE
! !INTERFACE:

 subroutine ovf_init_groups()

! !DESCRIPTION:
!  This routine initializes the ovf_groups structure, which lets
!  us know which processes are involved in which overflows

   ! input

   ! local
   integer (int_kind) :: n, m, region, ovf_id, i, j
   integer (int_kind) :: jmax, jmin, imax, imin
   integer (int_kind) :: count, m_end, my_pos
   integer (int_kind) :: num_ids, id, proc_count, proc_alloc_size
   integer (int_kind) :: nbors_alloc_size, min_task, cmp_task

   integer (int_kind) :: procs(1)
   integer (int_kind), pointer :: starts(:)

   logical (log_kind) :: found, comm_master_present
   real (r8), dimension(:,:,:), pointer :: g_mask !the mask

   allocate(ids(num_ovf))
   count = 0



   ! now calculate the number of groups that I am active in
   

   !serial: means active in all the groups
   count = num_ovf

   allocate(ovf_groups%groupIds(count))
   do n = 1, num_ovf
      ovf_groups%groupIds(n) = n
   end do

   ovf_groups%numTotal = num_ovf;
   ovf_groups%numMyGroups = count;
   ovf_groups%init = .true.


 
   !------------find neighbors------------------------!

   !serial = no neighbors

   !initial storage
   allocate(ovf_groups%neighborStarts(ovf_groups%numMyGroups + 1))
   ovf_groups%neighborStarts(1) = 1
   starts => ovf_groups%neighborStarts

   allocate(ovf_groups%commHandle(ovf_groups%numMyGroups))

   do n=1,ovf_groups%numMyGroups

      count = 0; !count neighbor procs for each group
      ovf_id = ovf_groups%groupIds(n) !overflow index

      !set the starts so we know where the neighbors are
      starts(n+1) = starts(n) + count

      !create the comm handle for this ovf_id
      ovf_groups%commHandle(n) = mpi2s_gshalo_init_symm_same &
           (distrb_clinic%communicator, count, &
           count , procs , .false.)

   end do ! my ovf loop


   !now for each group that I am in, determine whether I am the master task
   ! (this is used for comunicating the restart info, though it could be
   ! used for other things later.  Proc with the min id in the group is 
   ! the master of that group - EXCEPTION - the master task overall is always
   ! the mater of the groups that its in
   allocate(ovf_groups%amMaster(ovf_groups%numMyGroups))
   do n=1,ovf_groups%numMyGroups
      ovf_groups%amMaster(n) = .true.
   end do
   
 end subroutine ovf_init_groups


!***********************************************************************

! !IROUTINE
! !INTERFACE:

 subroutine ovf_utils_avgs(time_level)

! !DESCRIPTION: obtaining all the globals sums for use in ovf_reg_avgs

  include 'mpif.h'  ! MPI Fortran include file

!
! !REVISION HISTORY:
!  same as module

!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &! time indices for prognostic arrays
      time_level          ! current time level  (n)

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind)          :: &
      iblock,k,n,nn,m           ! dummy loop indices

   integer (int_kind)         :: m_end, & ! how many m loops to do	
			         region  !which region we are doing
	
   integer (int_kind)         :: g_kmin, g_kmax, ierr, ovf_id

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
         WRK                    ! temp work array

   real (r8)               :: vsum_wght  ! vertical sum (regional or adjacent) weight

   logical (log_kind)      :: use_rho, use_tracer

   real (r8)               :: local_sum_rho, global_sum_rho, &
        local_sum_wght, global_sum_wght
   real (r8), dimension(:), allocatable :: local_sum_tracer, global_sum_tracer 
   real (r8), dimension(:), allocatable :: work_array

  !! these need to be pointers so we can update original
   real (r8), dimension(:), pointer   :: g_trcr !tracer
   real (r8), pointer                 :: g_rho !density    
   real (r8), pointer                 :: g_wght
   real (r8), dimension(:,:,:), pointer :: &
      g_mask                    !the mask

   type (Schedule_t), pointer :: comm_handle

   allocate(work_array(length_comm))

   ! go thru my overflows
   do n=1, ovf_groups%numMyGroups
      ovf_id = ovf_groups%groupIds(n)
      comm_handle => ovf_groups%commHandle(n)

       !go thru each region, the 6 regions are: 1 = inflow, 2 = source, 3 = source adjacent, 
       !4 = entrainment, 5 = entrainment adjacent, 6 = product adjacent  

       ! these regions may overlap, but may or may not be on the same vertical level (k)
       do region = 1,6 

          select case (region)  
          case(1)  ! inflow
             g_wght => ovf(ovf_id)%wght_reg%inf
             g_mask => ovf(ovf_id)%mask_reg%inf
             g_trcr => ovf(ovf_id)%trcr_reg%inf
             g_rho =>  ovf(ovf_id)%rho_reg%inf
             g_kmin =  ovf(ovf_id)%reg_inf%kmin
             g_kmax =  ovf(ovf_id)%reg_inf%kmax
             use_rho = .true.	
             m_end = 1
	     use_tracer = .true.	

          case (2) !source
             g_wght => ovf(ovf_id)%wght_reg%src
             g_mask => ovf(ovf_id)%mask_reg%src
             g_trcr => ovf(ovf_id)%trcr_reg%src
             g_rho =>  ovf(ovf_id)%rho_reg%src
             g_kmin =  ovf(ovf_id)%reg_src%kmin
             g_kmax =  ovf(ovf_id)%reg_src%kmax   
             use_rho = .true. 
             m_end = 1
	     use_tracer = .true.	

          case (3) ! source adjacent
             g_wght => ovf(ovf_id)%wght_adj%src
             g_mask => ovf(ovf_id)%mask_adj%src
             g_trcr => ovf(ovf_id)%trcr_adj%src
             !g_rho = 
             g_kmin =  ovf(ovf_id)%adj_src%kmin
             g_kmax =  ovf(ovf_id)%adj_src%kmax
             use_rho = .false.
             m_end = 1
	     use_tracer = .true.	

         case (4) !entrainment region
             g_wght => ovf(ovf_id)%wght_reg%ent
             g_mask => ovf(ovf_id)%mask_reg%ent
             g_trcr => ovf(ovf_id)%trcr_reg%ent
             g_rho =>  ovf(ovf_id)%rho_reg%ent
             g_kmin =  ovf(ovf_id)%reg_ent%kmin
             g_kmax =  ovf(ovf_id)%reg_ent%kmax   
             use_rho = .true.
             m_end = 1
	     use_tracer = .true.	

         case (5) ! entrainment adjacent
             g_wght => ovf(ovf_id)%wght_adj%ent
             g_mask => ovf(ovf_id)%mask_adj%ent
             g_trcr => ovf(ovf_id)%trcr_adj%ent
             !g_rho =  
             g_kmin =  ovf(ovf_id)%adj_ent%kmin
             g_kmax =  ovf(ovf_id)%adj_ent%kmax
             use_rho = .false.
             m_end = 1
	     use_tracer = .true.	

         case (6) ! product adjacent
             !g_trcr
             use_rho = .true.
             m_end = ovf(ovf_id)%num_prd_sets
             use_tracer = .false.
             ! the following also depend on m -> set to m=1
             g_mask => ovf(ovf_id)%mask_adj%prd(:,:,:,1)
	     g_wght => ovf(ovf_id)%wght_adj%prd(1)
             g_rho =>  ovf(ovf_id)%rho_adj%prd(1)
             g_kmin = ovf(ovf_id)%adj_prd(1)%kmin
             g_kmax = ovf(ovf_id)%adj_prd(1)%kmax
            
       end select

       do m = 1, m_end

          if ( region == 6) then !prod adj. case
             g_mask => ovf(ovf_id)%mask_adj%prd(:,:,:,m)
             g_wght => ovf(ovf_id)%wght_adj%prd(m)
             g_rho =>  ovf(ovf_id)%rho_adj%prd(m)
             g_kmin =  ovf(ovf_id)%adj_prd(m)%kmin
             g_kmax =  ovf(ovf_id)%adj_prd(m)%kmax
          end if

          !calculating the area: dx*dy (only happens the first time when == c0)
           if (g_wght  .eq. c0) then
             do iblock = 1,nblocks_clinic 
               WRK(:,:,iblock) = DXT(:,:,iblock)*DYT(:,:,iblock)  &
                     *g_mask(:,:,iblock)
             end do

             local_sum_wght = ovf_local_sum(WRK, distrb_clinic, field_loc_center) 

             g_wght = local_sum_wght

          end if
       end do !end of 1st m loop

       if (use_tracer) then
          do nn = 1,nt
             g_trcr(nn) = c0
          end do
       end if

       do m=1, m_end

          if ( region == 6) then !prod adj. case
             g_mask => ovf(ovf_id)%mask_adj%prd(:,:,:,m)
             g_wght => ovf(ovf_id)%wght_adj%prd(m)
             g_rho =>  ovf(ovf_id)%rho_adj%prd(m)
             g_kmin = ovf(ovf_id)%adj_prd(m)%kmin
             g_kmax = ovf(ovf_id)%adj_prd(m)%kmax
          end if

          vsum_wght  = c0
          if (use_rho) then
             g_rho = c0
          end if
          local_sum_rho = c0

          allocate(local_sum_tracer(nt))!nt = num tracers
          local_sum_tracer = c0
          allocate(global_sum_tracer(nt))
          

!combine local sum into iblock loop for better cache use
! check perf. on higher resoltuion to see if I need to do 1 comm per
! ovf (not per region)
          do k = g_kmin, g_kmax
             vsum_wght = vsum_wght + g_wght*dz(k)
             if (use_rho) then
                do iblock = 1,nblocks_clinic
                   WRK(:,:,iblock) = RHO(:,:,k,time_level,iblock) &
                        *DXT(:,:,iblock)*DYT(:,:,iblock)*dz(k)   &
                        *g_mask(:,:,iblock)
                end do
                local_sum_rho = local_sum_rho + &
                     ovf_local_sum(WRK, distrb_clinic, field_loc_center)
             end if
             if (use_tracer) then
                do nn = 1,nt
                   do iblock = 1,nblocks_clinic
                      WRK(:,:,iblock) = TRACER(:,:,k,nn,time_level,iblock) &
                           *DXT(:,:,iblock)*DYT(:,:,iblock)*dz(k)         &
                           *g_mask(:,:,iblock)
                   end do
                   local_sum_tracer(nn) = local_sum_tracer(nn) + &
                        ovf_local_sum(WRK,distrb_clinic,field_loc_center)
                end do
             end if
          end do !end k loop

          work_array(1) = local_sum_rho
          work_array(2:nt+1) = local_sum_tracer(1:nt)

          !serial
          if (use_rho) then
             global_sum_rho = work_array(1)
             g_rho = g_rho + global_sum_rho
             g_rho = g_rho / vsum_wght
          end if
          if (use_tracer) then
             global_sum_tracer(1:nt) = work_array(2:nt+1)  
             do nn = 1, nt
                g_trcr(nn) = g_trcr(nn) + global_sum_tracer(nn)
                g_trcr(nn) = g_trcr(nn) / vsum_wght
             end do
          end if

           deallocate(local_sum_tracer)
           deallocate(global_sum_tracer)


      end do ! end of 2nd m loop
    end do !end of loop over regions
 end do !end of num_ovf loop	



deallocate(work_array)


end subroutine ovf_utils_avgs



!***********************************************************************
! !IROUTINE: ovf_local_sum

function  ovf_local_sum(X, dist, field_loc)


! !DESCRIPTION:
!  computes the local sum of the _physical domain_ of a 2-d
!  array (doubles)
!

! !INPUT PARAMS

  real (r8), dimension(:,:,:), intent(in) :: &
       X                    ! array to be summed
  
  type (distrb), intent(in) :: &
       dist                 ! block distribution for array X

  integer (int_kind), intent(in) :: &
        field_loc            ! location of field on staggered grid


! !OUTPUT PARAMETERS:

  real (r8) ::   &
       ovf_local_sum       ! resulting global sum


  
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    real (r8), dimension(:), allocatable :: &
         local_block_sum ! sum of local blocks


    real (r8) ::          &
         block_sum,         &! sum of local block
         local_sum           ! sum of all local blocks

    integer (int_kind) :: &
         i,j,n,             &! local counters
         ib,ie,jb,je,       &! beg,end of physical domain
         bid                 ! block location
    
    type (block) :: &
         this_block          ! holds local block information
    

   local_sum = c0

   if (ltripole_grid .and. (field_loc == field_loc_Nface .or. &
        field_loc == field_loc_NEcorner)) then
      !*** must remove redundant points from sum
      do n=1,nblocks_tot
         if (dist%proc(n) == my_task+1) then
            block_sum = c0
            bid = dist%local_block(n)
            this_block = get_block(n,bid)
            ib = this_block%ib
            ie = this_block%ie
            jb = this_block%jb
            je = this_block%je
            if (this_block%jblock == nblocks_y) then
               !*** for topmost row in tripole only sum
               !*** 1st half of domain - others are redundant
               do i=ib,ie
                  if (this_block%i_glob(i) <= nx_global/2) &
                       block_sum = block_sum + X(i,je,bid)
               end do
               
               je = je - 1
            endif
            do j=jb,je
               do i=ib,ie
                  block_sum = block_sum + X(i,j,bid)
               end do
            end do
            local_sum = local_sum + block_sum
         endif
      end do !block loop
   else ! regular global sum
      do n=1,nblocks_tot
         if (dist%proc(n) == my_task+1) then
            block_sum = c0
            bid = dist%local_block(n)
            call get_block_parameter(n,ib=ib,ie=ie,jb=jb,je=je)
            do j=jb,je
               do i=ib,ie
                  block_sum = block_sum + X(i,j,bid)
               end do
            end do
            local_sum = local_sum + block_sum
         endif
      end do !block loop
   endif


   !return value
   ovf_local_sum = local_sum


end function ovf_local_sum



!***********************************************************************

! !IROUTINE
! !INTERFACE:

 subroutine ovf_print_init(len, num_posts, rbuff, post_array)

! !DESCRIPTION: 

  include 'mpif.h'  ! MPI Fortran include file

!
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

  integer (int_kind), intent(in) :: len  !length of each recv

  integer (int_kind), intent(in) :: num_posts  !number of posts

  real (r8), intent(in), optional:: rbuff(:) !buffer to collect recvs

  logical (log_kind), intent(inout), &
       optional :: post_array(num_ovf)  !array of whether to post or not


!----------------------------------------------------------------------
!
!  local variables
!
!--------------------------------------------------------------------

  integer (int_kind) :: n, count, ovf_id, loc, ierr, i
  logical (log_kind) :: post

  ! no posts
  if ( my_task == master_task .and. present(post_array)) then
     do n = 1, num_ovf
        post_array(n) = .false.
     end do
  end if
  
end subroutine ovf_print_init


!***********************************************************************

! !IROUTINE
! !INTERFACE:

 subroutine ovf_print_send(len, sbuff, post_num, tag)

! !DESCRIPTION: 

  include 'mpif.h'  ! MPI Fortran include file

!
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

  integer (int_kind), intent(in) :: len  !length of send

  real (r8), intent(in), optional:: sbuff(:) !send buffer

  integer (int_kind), intent(in) :: post_num !the send number

  integer (int_kind), intent(in) :: tag ! tag for messsage (should be ovf_id)

!----------------------------------------------------------------------
!
!  local variables
!
!--------------------------------------------------------------------

  integer (int_kind) :: ierr
  logical (log_kind) :: post

 
  !nothing to do
  return

 end subroutine ovf_print_send


!***********************************************************************

! !IROUTINE
! !INTERFACE:

 subroutine ovf_print_get(post_num)

! !DESCRIPTION: 
   
  include 'mpif.h'  ! MPI Fortran include file

!
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

  integer (int_kind), intent(in) :: post_num !the recv number to get

!----------------------------------------------------------------------
!
!  local variables
!
!--------------------------------------------------------------------

  integer (int_kind) :: ierr
  
  

  !nothing to do
  return
  
end subroutine ovf_print_get


!***********************************************************************

! !IROUTINE
! !INTERFACE:

 subroutine ovf_print_finalize(num_posts)

! !DESCRIPTION: 

  include 'mpif.h'  ! MPI Fortran include file

!
!-----------------------------------------------------------------------
!  input variables
!-----------------------------------------------------------------------

  integer (int_kind), intent(in) :: num_posts !the send number


!----------------------------------------------------------------------
!
!  local variables
!
!--------------------------------------------------------------------

  integer (int_kind) :: ierr
  
  !nothing to do
  return
  
end subroutine ovf_print_finalize


!***********************************************************************

 end module ovf_utils

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
