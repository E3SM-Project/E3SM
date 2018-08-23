! Copyright (c) 2014,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_particle_list
!
!> \brief Particle framework
!> \author Phillip Wolfram
!> \date   04/10/2014
!> \details
!>  This module contains a general definition of particles which can be
!>  used in implementation of Lagrangian Particles Tracking.
!-----------------------------------------------------------------------
#define COMMA ,
#ifdef MPAS_DEBUG
#define LIGHT_DEBUG_WRITE(Z)   print *, Z
#define LIGHT_DEBUG_ALL_WRITE(Z)  print *, Z
#else
#define LIGHT_DEBUG_WRITE(Z)  ! print *, Z
#define LIGHT_DEBUG_ALL_WRITE(Z) ! print *, Z
#endif

#define LIGHT_ERROR_WRITE(Z)  print *, Z
#define LIGHT_WARNING_WRITE(Z)  print *, Z

module ocn_particle_list

  ! declare general packages used
#ifdef _MPI
  use mpi
#endif
  use mpas_kind_types
  use mpas_io_units
  use mpas_derived_types
  use mpas_dmpar
  use mpas_block_decomp
  use mpas_pool_routines
  use mpas_timer

  ! blanket statments to restrict implicit module's scope
  implicit none
  private

  ! mpi defines
#ifdef _MPI
   integer, parameter :: MPI_INTEGERKIND = MPI_INTEGER
   integer, parameter :: MPI_2INTEGERKIND = MPI_2INTEGER

#ifdef SINGLE_PRECISION
   integer, parameter :: MPI_REALKIND = MPI_REAL
   integer, parameter :: MPI_2REALKIND = MPI_2REAL
#else
   integer, parameter :: MPI_REALKIND = MPI_DOUBLE_PRECISION
   integer, parameter :: MPI_2REALKIND = MPI_2DOUBLE_PRECISION
#endif
#endif

  ! custom structures defined in mpas_grid_types

  ! define private interfaces

  ! add routines
  interface add_halo_data_to_particle_list
    !(particlelist, dataName, data)
    module procedure add_halo_data_to_particle_list_1Dreal
    module procedure add_halo_data_to_particle_list_1Dint
  end interface

  interface add_halo_data_to_particle_list_array
    module procedure add_halo_data_to_particle_list_1Dreal_array
    module procedure add_halo_data_to_particle_list_1Dint_array
  end interface

! interface add_nonhalo_data_to_particle_list
!   !(particlelist, dataName, data)
!   module procedure add_nonhalo_data_to_particle_list_1Dreal
!   module procedure add_nonhalo_data_to_particle_list_1Dint
! end interface

! interface add_nonhalo_data_to_particle_list_array
!   module procedure add_nonhalo_data_to_particle_list_1Dreal_array
!   module procedure add_nonhalo_data_to_particle_list_1Dint_array
! end interface

  ! get routines
  interface get_halo_data_from_particle_list
    !(particlelist, dataName, data)
    module procedure get_halo_data_from_particle_list_1Dreal
    module procedure get_halo_data_from_particle_list_1Dint
  end interface

  interface get_halo_data_from_particle_list_array
    module procedure get_halo_data_from_particle_list_1Dreal_array
    module procedure get_halo_data_from_particle_list_1Dint_array
  end interface

! interface get_nonhalo_data_from_particle_list
!   !(particlelist, dataName, data)
!   module procedure  get_nonhalo_data_from_particle_list_1Dreal
!   module procedure  get_nonhalo_data_from_particle_list_1Dint
! end interface

! interface get_nonhalo_data_from_particle_list_array
!   !(particlelist, dataName, data)
!   module procedure  get_nonhalo_data_from_particle_list_1Dreal_array
!   module procedure  get_nonhalo_data_from_particle_list_1Dint_array
! end interface

  !-----------------------------------------------------------------
  ! public routines and interfaces
  !-----------------------------------------------------------------
  ! define publically accessible subroutines, functions, interfaces
  public :: mpas_particle_list_build_and_assign_particle_list
  public :: mpas_particle_list_destroy_particle_list, mpas_particle_list_remove_particles_not_on_current_block
  public :: mpas_particle_list_build_computation_halos, mpas_particle_list_build_halos
  public :: mpas_particle_list_update_particle_block
  public :: mpas_particle_list_update_halos_start, mpas_particle_list_update_halos_end
  public :: mpas_particle_list_transfer_particles_from_block_to_named_block
  public :: mpas_particle_list_write_halo_data !, mpas_particle_list_write_nonhalo_data
  public :: mpas_particle_list_test_neighscalc, mpas_particle_list_test_numparticles_to_neighprocs
  public :: mpas_particle_list_test_num_current_particlelist
  public :: mpas_particle_list_self_union_halo_lists

  ! subroutine / function definition
contains

!***********************************************************************
!
!  routine mpas_particle_list_build_and_assign_particle_list
!
!> \brief   Allocates particles for initialization
!> \author  Phillip Wolfram
!> \date    07/01/2014
!> \details
!>  This routine builds and allocates particlces following initalization
!
!-----------------------------------------------------------------------
   subroutine mpas_particle_list_build_and_assign_particle_list(domain,err) !{{{
     implicit none

     !-----------------------------------------------------------------
     !
     ! input variables
     !
     !-----------------------------------------------------------------

     !-----------------------------------------------------------------
     !
     ! input/output variables
     !
     !-----------------------------------------------------------------

     type (domain_type), intent(in) :: domain

     !-----------------------------------------------------------------
     !
     ! output variables
     !
     !-----------------------------------------------------------------

     integer, intent(out) :: err !< Output: error flag

     err = 0
     ! build / allocate the listPLSend, setting ioProc (currentBlock read in as haloData)
     call build_block_particlelists(domain, err)
#ifdef MPAS_DEBUG
     LIGHT_DEBUG_ALL_WRITE('finished build_block_particlelists 1')
     call MPI_Barrier(domain % dminfo % comm, err)
#endif

     ! read in data from netCDF-injected data structures
     ! nonhalo-data is diagnostic
     call read_haloData(domain, err)
#ifdef MPAS_DEBUG
     LIGHT_DEBUG_ALL_WRITE('finished read_haloData 1')
     call MPI_Barrier(domain % dminfo % comm, err)
#endif

     ! currently there is no nonhalo data
     !! note that nonhalo data is just initialized with 0
     !! values are not imported from netCDF input file
     !call read_nonhaloData(domain, err)
#ifdef MPAS_DEBUG
     !LIGHT_DEBUG_ALL_WRITE('finished read_nonhaloData 1')
     !call MPI_Barrier(domain % dminfo % comm, err)

     !call mpas_particle_list_test_num_current_particlelist(domain)
     !LIGHT_DEBUG_ALL_WRITE('finished test in mpas_particle_list_build_and_assign_particle_list')
     !call MPI_Barrier(domain % dminfo % comm, err)
#endif
     !! test to make sure deallocation is ok before transfer...!{{{
#ifdef MPAS_DEBUG
     LIGHT_DEBUG_ALL_WRITE('    Trying to clear particlelist memory on blocks')
     call clear_block_particlelists(domain,err)
     call build_block_particlelists(domain, err)
     call read_haloData(domain, err)
     !call read_nonhaloData(domain, err)
     call mpas_particle_list_test_num_current_particlelist(domain)
     call test_currentBlock(domain)
     LIGHT_DEBUG_ALL_WRITE('    Rebuilt data structures-- ok')
     call MPI_Barrier(domain % dminfo % comm, err)
#endif
     !}}}
   end subroutine mpas_particle_list_build_and_assign_particle_list !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_particle_list_destroy_particle_list
!
!> \brief MPAS destroy particlelist
!> \author Phillip Wolfram
!> \date   06/27/2014
!> \details
!>  This routine destroys a particlelist, deallocating its memory
!>  including that of all pointers it contains
!
!-----------------------------------------------------------------------
subroutine mpas_particle_list_destroy_particle_list(particlelist) !{{{
  type (mpas_particle_list_type), pointer, intent(in) :: particlelist
  type (mpas_particle_list_type), pointer :: pLCurr, pLCurrTemp

  if(associated(particlelist)) then
    plCurr => particlelist
    do while(associated(plCurr))
      pLCurrTemp => pLCurr
      pLCurr => pLCurr % next
      ! destroy the particle too
      if(associated(pLCurrTemp % particle)) then
        call destroy_particle(pLCurrTemp % particle)
      end if
      deallocate(pLCurrTemp)
    end do
  end if

end subroutine mpas_particle_list_destroy_particle_list !}}}

!***********************************************************************
!
!  routine mpas_particle_list_remove_particles_not_on_current_block
!
!> \brief   Remove particles on block % particlelist that are not on currentBlock
!> \author  Phillip Wolfram
!> \date    07/03/2014
!> \details
!>  This routine removes particles that were transfered strictly for IO.
!>  If the particle's currentBlock is not the same as the current block,
!>  the particle is removed.
!
!-----------------------------------------------------------------------
   subroutine mpas_particle_list_remove_particles_not_on_current_block(domain, err) !{{{
     implicit none

     !-----------------------------------------------------------------
     !
     ! input variables
     !
     !-----------------------------------------------------------------

     !-----------------------------------------------------------------
     !
     ! input/output variables
     !
     !-----------------------------------------------------------------

     type (domain_type), intent(in) :: domain

     !-----------------------------------------------------------------
     !
     ! output variables
     !
     !-----------------------------------------------------------------

     integer, intent(out) :: err !< Output: error flag

     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------
     type (block_type), pointer :: block
     type (mpas_particle_list_type), pointer :: particlelist, particlelisttemp, particlelisttemp2
     integer :: thisBlock
     integer, pointer :: particleBlock
     integer :: arrayIndex
     !type (mpas_pool_type), pointer :: lagrPartTrackPool
     type (mpas_particle_type), pointer :: particle

     err = 0

     block => domain % blocklist
     do while(associated(block))
       ! for each particle on each block
       particlelist => block % particlelist
       do while(associated(particlelist))
         particle => particlelist % particle
         call mpas_pool_get_array(particle % haloDataPool, 'currentBlock', particleBlock)

         if(particleBlock /= block % blockID) then
           ! REMOVE PARTICLE FROM EXISTING PARTICLELIST
           ! remove particle and particle reference from existing particle list
           ! 3 cases: head, middle, tail

           call destroy_particle(particle)

           if(associated(particlelist % prev)) then
             particlelisttemp => particlelist % prev
             if (associated(particlelist % next)) then
               ! case of the middle
               particlelisttemp % next => particlelist % next
               particlelisttemp2 => particlelisttemp
               ! want to keep particle memory intact because their pointers were passed previously,
               ! so just empty the list, don't destroy it and its contents
               particlelisttemp => particlelist % next
               particlelisttemp % prev => particlelisttemp2
               ! just need to remove the single link, particle memory needs to stay intact
               !deallocate(particlelist % particle)
               deallocate(particlelist)
               ! get next pointer
               particlelist => particlelisttemp
             else
               ! case of tail
               nullify(particlelisttemp % next)
               !deallocate(particlelist % particle)
               deallocate(particlelist)
               ! set back to final
             end if
           else
             if(associated(particlelist % next)) then
               ! case of head, set new head (assumes more than one particle)
               particlelisttemp => particlelist % next
               nullify(particlelisttemp % prev)
               block % particlelist => particlelisttemp
               !deallocate(particlelist % particle)
               deallocate(particlelist)
               particlelist => particlelisttemp
             else
               ! case of single link / particle
               !deallocate(particlelist % particle)
               deallocate(particlelist)
               nullify(block % particlelist)
             end if
           end if
         else
           particlelist => particlelist % next
         end if
       end do

       ! this is done for each block because we want processor - processor communication
       block => block % next
     end do

   end subroutine mpas_particle_list_remove_particles_not_on_current_block !}}}

!***********************************************************************
!
!  routine  mpas_particle_list_build_computation_halos
!
!> \brief   Build up necessary info for communication of particle in
!>          halo to neighboring cell during computation step.
!> \author  Phillip Wolfram
!> \date    07/02/2014
!> \details
!>  This routine builds g_compProcNeighs which is the neighboring processor
!>  list needed to process MPI communication, assuming a list of
!>  particlelists is built up corresponding to the processors in this
!>  array.  The end result is that g_compProcNeighs is populated.
!
!-----------------------------------------------------------------------
   subroutine mpas_particle_list_build_computation_halos(domain, err, procNeighs) !{{{
     implicit none

     !-----------------------------------------------------------------
     !
     ! input variables
     !
     !-----------------------------------------------------------------

     !-----------------------------------------------------------------
     !
     ! input/output variables
     !
     !-----------------------------------------------------------------

     type (domain_type), intent(inout) :: domain

     !-----------------------------------------------------------------
     !
     ! output variables
     !
     !-----------------------------------------------------------------

     integer, intent(out) :: err !< Output: error flag
     integer, dimension(:), pointer :: procNeighs

     err = 0

     ! compute the cellOwnerBlock for cells
     ! owning processor can be obtained from
     ! mpas_get_owning_proc in mpas_block_decomp (src/framework)
     call compute_cellOwnerBlock(domain, err)

     ! in order to mediate MPI exchanges, need to determine
     ! 1. list of blockNeighs (global)
     ! stored in block % blockNeighs
     ! basically Neighs are processors who own the halo including itself
     call compute_blockNeighs(domain, err)

     ! 2. list of procNeighs (extracted from blockNeighs, also global)
     ! this information is necessary in order to know where to send data
     ! (this is like a block exchange list for particles)
     ! stored in block % procNeighs
     ! this is just the list of processors who own blockNeighs
     call compute_block_procNeighs(domain, err)

     ! need to aggregate procNeighs to be global for the processor (over each block on the
     ! processor), total number of neighboring processors to a particular processor
     call compute_procNeighs(domain, err, procNeighs)

   end subroutine mpas_particle_list_build_computation_halos !}}}

!***********************************************************************
!
!  routine mpas_particle_list_build_halos
!
!> \brief   Build the IO halo information to transmit particles from their
!>          initial host IO processor to the appropriate currentBlock
!>          processor.
!> \author  Phillip Wolfram
!> \date    07/02/2014
!> \details
!>  This routine builds the IO halo information to transmit  particles
!>  from their initial host IO processor to the appropriate currentBlock
!>  processor.  The end result is that g_ioProcNeighs is populated.
!
!-----------------------------------------------------------------------
    subroutine mpas_particle_list_build_halos(domain, err, namedBlock, ioProcNeighs) !{{{
      !{{{ initialization
      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      character(len=*), intent(in) :: namedBlock

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:), pointer, intent(out) :: ioProcNeighs
      integer, intent(out) :: err !< Output: error flag

      !}}}
      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer, dimension(:), pointer :: currentBlocks, tempInt
      integer :: i, nBlocks, nTotProcs, mpi_ierr
      logical, dimension(:), pointer :: sendNeigh, recvNeigh

      err = 0

      ! get unique list of currentBlock to determine all currentProcs
      ! that the particles must be communicated to
      ! determine complete set of ioProcs, which are assigned based on the
      ! way that PIO decomposes the nParticle dimension.
      ! at this point, processors owning the particle (at read-in) are
      ! assumed to be the ioProcessor
      call compute_all_particle_values_unique_int(domain, err, namedBlock, currentBlocks)

      nTotProcs = domain % dminfo % nprocs
      allocate(sendNeigh(nTotProcs))
      allocate(recvNeigh(nTotProcs))
      sendNeigh = .false.
      recvNeigh = .false.

      ! Get processors for currentBlocks. Note, however, this is one-sided because
      ! only the IO processors know their send location, the receivers don't know
      ! their sending processors
      if(associated(currentBlocks)) then
        nBlocks = size(currentBlocks)
        !print *,  'nBlocks = ', nBlocks, ' currentBlocks = ', currentBlocks
        allocate(ioProcNeighs(nBlocks))
        do i=1, nBlocks
          call mpas_get_owning_proc(domain % dminfo, currentBlocks(i), ioProcNeighs(i))
        end do

        ! note, each ioProc knows where it is sending data.  However, those processors
        ! do not know they should receive data (their halos are empty).  The halos
        ! are not symmetric and the communication cannot occur unless this is fixed.
        ! this is fixed via an all-to-all communication (only at initialization, otherwise
        ! parallelism can be broken)
        ! set counter for connectivity
        do i=1,size(ioProcNeighs)
          sendNeigh(ioProcNeighs(i)+1) = .true.
        end do
        deallocate(ioProcNeighs)
      end if

      LIGHT_DEBUG_WRITE('sendNeigh = ' COMMA sendNeigh)
      LIGHT_DEBUG_WRITE('recvNeigh= ' COMMA recvNeigh)
#ifdef _MPI
      !call MPI_Barrier(domain % dminfo % comm, mpi_ierr)
#endif
      ! send with MPI all to update in recvNeigh (should only have to be done once)
#ifdef _MPI
      call MPI_Alltoall(sendNeigh, 1, MPI_LOGICAL, recvNeigh, 1, MPI_LOGICAL, domain % dminfo % comm, mpi_ierr)
#endif
#ifdef _MPI
      !call MPI_Barrier(domain % dminfo % comm, err)
#endif
      LIGHT_DEBUG_WRITE('Finished all to all with mpi_ierr = ' COMMA mpi_ierr)
      LIGHT_DEBUG_WRITE('sendNeigh = ' COMMA sendNeigh)
      LIGHT_DEBUG_WRITE('recvNeigh= ' COMMA recvNeigh)
      LIGHT_DEBUG_WRITE('communicate = ' COMMA sendNeigh .or. recvNeigh)

      ! could possibly optimize here by keeping track of send / recv lists separately
      ! however, if there is nothing to be sent the only message that is sent
      ! is the number of particles to be transferred...
      ! update ioProcNeighs from logical lists
      ! "add" the lists
      recvNeigh = recvNeigh .or. sendNeigh
      allocate(tempInt(nTotProcs))
      ! this can artificially create a problem if there isn't a single block to a processor
      tempInt = domain % dminfo % my_proc_id
      do i=1,nTotProcs
        if (recvNeigh(i))   tempInt(i) = i-1
      end do
      deallocate(sendNeigh)
      deallocate(recvNeigh)

      ! get a complete list of the processors (including itself)
      call uniqueIntegerList(tempInt,ioProcNeighs)
      call removeValueFromIntList(ioProcNeighs, domain % dminfo % my_proc_id)
      LIGHT_DEBUG_WRITE('ioProcNeighs = ' COMMA ioProcNeighs)
      deallocate(tempInt)

    end subroutine mpas_particle_list_build_halos !}}}

!***********************************************************************
!
!  routine  mpas_particle_list_update_halos_end
!
!> \brief   Updates halo processors for communication, noting that
!>          the receiving processors must be informed of changes
!> \author  Phillip Wolfram
!> \date    07/08/2014
!> \details
!>  This routine transmits a logical list of processors that will
!>  transmit data for each compProc or ioProc.  On the compProcs or
!>  ioProcs, these lists must be aggregated to build out the full list
!>  of processors from which data will be received.
!
!-----------------------------------------------------------------------
   subroutine mpas_particle_list_update_halos_end(domain, err, destinationName, sendProcNeighs, sendProcSendList, & !{{{
                                                  sendProcRecvList)
      implicit none
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      logical, dimension(:,:), intent(in) :: sendProcRecvList !< x: sendProcNeighs for send. y: each receiving
                                                              !<    processors on x denoted by true
      ! 'ioBlock' and 'currentBlock' are options for destinationName
      character(len=*), intent(in) :: destinationName
      logical, dimension(:), pointer, intent(inout) :: sendProcSendList   !< location of true indicates processors to send data to
      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:), pointer, intent(inout)   :: sendProcNeighs !< list of io halo processors

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i, nsendProcNeighs, nProcs
      logical, dimension(:), pointer :: completeList
      logical, dimension(:,:), pointer :: recvList
      integer, dimension(:), pointer :: intArray
      integer, dimension(:), pointer :: sendRequestID, recvRequestID
      integer :: mpi_ierr
      logical :: firstTime

      err = 0

      ! compute send list now that all particles reside on correct block (processor) 'currentBlock'
      ! didn't show up with serial IO because all computational processors sent data to proc 0
      ! note that this could be missing communication where before transfer particle on proc A
      ! has sendProc of A and is sent to B (must have previously kept a record that B is in A's halo list
      ! note, may be slightly redundant because we could update sendProcSendList once particles are transfered
      call compute_additional_particle_send_list(domain, destinationName, sendProcSendList)

      ! need to update IO processors as to the change also so that they know where to get data from!!!
      LIGHT_DEBUG_WRITE('destinationName= ' COMMA destinationName)
      LIGHT_DEBUG_WRITE('sendProcSendList = ' COMMA sendProcSendList)
      LIGHT_DEBUG_WRITE('sendProcRecvList = ' COMMA sendProcRecvList)
      LIGHT_DEBUG_WRITE('sendProcNeighs before = ' COMMA sendProcNeighs)

      ! proceed to update the halo
      nProcs = domain % dminfo % nprocs
      nsendProcNeighs = size(sendProcNeighs)
      allocate(completeList(nProcs), recvList(nProcs,nsendProcNeighs))
      allocate(sendRequestID(nsendProcNeighs), recvRequestID(nsendProcNeighs))

      completeList = .False.
      recvList = .False.
      ! for each sendProc, send logical array information
      do i = 1, nsendProcNeighs
#ifdef _MPI
        call MPI_ISend(sendProcRecvList(:,i), nProcs, MPI_LOGICAL, sendProcNeighs(i), domain % dminfo % my_proc_id, &
          domain % dminfo % comm, sendRequestID(i), mpi_ierr)
#endif
        LIGHT_DEBUG_WRITE('sendProcNeigh= ' COMMA sendProcNeighs(i))
        LIGHT_DEBUG_WRITE('send data = ' COMMA sendProcRecvList(i,:))
#ifdef _MPI
        call MPI_IRecv(recvList(:,i), nProcs, MPI_LOGICAL, sendProcNeighs(i), sendProcNeighs(i), &
          domain % dminfo % comm, recvRequestID(i), mpi_ierr)
#endif
      end do

#ifdef _MPI
      ! wait until the data is in the buffer
      call MPI_WaitAll(nsendProcNeighs,recvRequestID, MPI_STATUSES_IGNORE, mpi_ierr)
#endif

      do i = 1, nsendProcNeighs
        ! aggregate results after wait, making sure that we have the most
        ! comprehensive list of sendProcs for receiving
        LIGHT_DEBUG_WRITE('sendProcNeigh= ' COMMA sendProcNeighs(i))
        LIGHT_DEBUG_WRITE('recvList before = ' COMMA recvList(:,i))
        LIGHT_DEBUG_WRITE('completeList before = ' COMMA completeList)
        completeList = completeList .or. recvList(:,i)
        LIGHT_DEBUG_WRITE('recvList after = ' COMMA recvList(:,i))
        LIGHT_DEBUG_WRITE('completeList after = ' COMMA completeList)
      end do

      ! wait to make sure (just in case) that all sends have completed
#ifdef _MPI
      call MPI_WaitAll(nsendProcNeighs, sendRequestID, MPI_STATUSES_IGNORE, mpi_ierr)
#endif

      ! "add" receiving and sending lists into a comple list
      completeList = completeList .or. sendProcSendList
      !print *,  'completeList = ', completeList

      ! convert complete list into a unique list of processor numbers
      allocate(intArray(nProcs))
      firstTime = .True.
      do i = 1, nProcs
        if (completeList(i)) then
          if (firstTime) then
            intArray = i - 1
            firstTime = .False.
          else
            intArray(i) = i - 1
          end if
        end if
      end do

      ! now get the desired integer halo list
      deallocate(sendProcNeighs)
      call uniqueIntegerList(intArray, sendProcNeighs)
      call removeValueFromIntList(sendProcNeighs, domain % dminfo % my_proc_id)

      deallocate(intArray, completeList, sendRequestID, recvRequestID)
      ! need to update IO processors as to the change also so that they know where to get data from!!!
      ! should uncomment for testing when multiple sendProcs are utilized (parallel IO)
      LIGHT_DEBUG_WRITE('sendProcSendList after = ' COMMA sendProcSendList)
      LIGHT_DEBUG_WRITE('sendProcRecvList after = ' COMMA sendProcRecvList)
      LIGHT_DEBUG_WRITE('sendProcNeighs after = ' COMMA sendProcNeighs)

   end subroutine mpas_particle_list_update_halos_end !}}}

!***********************************************************************
!
!  routine  mpas_particle_list_transfer_particles_from_block_to_named_block
!
!> \brief   Move particles to the appropriate currentBlock
!> \author  Phillip Wolfram
!> \date    07/01/2014
!> \details
!>  This routine uses MPI communication to ensure particles end up
!>  on their appropriate currentBlock.
!
!-----------------------------------------------------------------------
   subroutine mpas_particle_list_transfer_particles_from_block_to_named_block(domain, err, & !{{{
       haloOnly, copyOnly, namedBlock, procNeighs)
     implicit none
     !-----------------------------------------------------------------
     !
     ! input variables
     !
     !-----------------------------------------------------------------
     integer, dimension(:), pointer, intent(in) :: procNeighs
     character(len=*), intent(in) :: namedBlock
     logical, intent(in) :: copyOnly, haloOnly

     !-----------------------------------------------------------------
     !
     ! input/output variables
     !
     !-----------------------------------------------------------------

     type (domain_type), intent(in) :: domain

     !-----------------------------------------------------------------
     !
     ! output variables
     !
     !-----------------------------------------------------------------

     integer, intent(out) :: err !< Output: error flag
     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------
     integer, dimension(:), pointer :: nPartSend => NULL(), nPartRecv => NULL()
     type (block_type), pointer :: block
     type (mpas_pool_type), pointer :: lagrPartTrackPool
     type (mpas_list_of_particle_list_type), dimension(:), pointer :: listPLSend => NULL(), listPLRecv => NULL()

     err = 0

     ! allocate particle list in terms of its communication dimesion to other processors
     allocate(listPLSend(size(procNeighs)),listPLRecv(size(procNeighs)))
     allocate(nPartSend(size(procNeighs)),nPartRecv(size(procNeighs)))

     ! if statement inside each of the group's subroutines needs removed for parallel IO which
     ! assumes a high-level nParticle decomposition which will allocate IO blocks via the
     ! decomposition


     ! this is the communication list and it needs to be a bi-directional graph so that
     ! sending processors have their receiving processor listening, even if there is no
     ! data transfer required

     ! presently don't distribute particles from each block to correct, owning blocks
     !1. communicating from block to block (basically like make_proc_to_proc_particlelists but for blocks)
     ! this is all done locally and all it will do is affect the particlelists on each block, forming temporary lists,
     ! and then appending particles on these temporary lists back to block % particlelist
     !this doesn't have to be tested unless there is more than one block per processor,
     !which presently is not how things are done
     ! this takes the block % particlelist on the first block and makes sure
     ! it is appropriately distributed on other blocks on the same processor
     ! call make_block_to_block_particlelists(domain)


     ! 1. Make temporary particle lists for transfers based on currentBock of cells in halo.
     !    These lists live on each block and correspond to other
     !    blocks that the list must be moved to.  The particle must end up on its currentBlock specified.
     ! convention on particles lists is that g_compProcNeighs specifies the processor neighbor numbers
     ! corresponding to each index in the particlelists pointer array.  Should use linked-list
     ! because processor neighbors are typically going to be somewhere less than 10
     ! (perfect partitioning of plane gives hexagons with 6 cell neighbors, for instance).
     ! strategy is to make linked list corresponding to each gProcNeigh and then append
     ! particles beloning to the list on the list
     LIGHT_DEBUG_WRITE('before make_proc_to_proc_particlelists')
#ifdef MPAS_DEBUG
     call MPI_Barrier(domain % dminfo % comm, err)
     call mpas_particle_list_test_numparticles_to_neighprocs(domain % dminfo % my_proc_id, procNeighs, procNeighs)
     call MPI_Barrier(domain % dminfo % comm, err)
#endif
     call make_proc_to_proc_particlelists(domain, copyOnly, namedBlock, listPLSend, procNeighs, err)
     LIGHT_DEBUG_WRITE('after make_proc_to_proc_particlelists')
#ifdef MPAS_DEBUG
     call MPI_Barrier(domain % dminfo % comm, err)
     call mpas_particle_list_test_numparticles_to_neighprocs(domain % dminfo % my_proc_id, procNeighs, procNeighs)
     call MPI_Barrier(domain % dminfo % comm, err)
#endif
     LIGHT_DEBUG_WRITE('finished make_proc_to_proc_particlelists')
#ifdef MPAS_DEBUG
     call MPI_Barrier(domain % dminfo % comm, err)
     call mpas_particle_list_test_num_current_particlelist(domain)
     call MPI_Barrier(domain % dminfo % comm, err)
     call test_num_particles_on_particlelist(listPLSend, size(procNeighs))
     call MPI_Barrier(domain % dminfo % comm, err)
#endif

     ! now move data from one proc to another, assuming that data is move to 1st block
     ! on foreign processor
     ! tell other processors how many particles are going to be communicated to them
     call get_num_particlelists(listPLSend, size(procNeighs), nPartSend)
#ifdef MPAS_DEBUG
     LIGHT_DEBUG_WRITE('finished get_num_particlelists')
     LIGHT_DEBUG_WRITE('proc id=' COMMA domain % dminfo % my_proc_id COMMA ' comm=' COMMA domain % dminfo % comm)
     LIGHT_DEBUG_WRITE('before 1st barrier communicate_num_particles_send_recv')
     call MPI_Barrier(domain % dminfo % comm, err)
     LIGHT_DEBUG_WRITE('after 1st barrier communicate_num_particles_send_recv')
#endif
     call communicate_num_particles_send_recv(domain, procNeighs, nPartSend, nPartRecv)
#ifdef MPAS_DEBUG
     ! Error is above of here
     LIGHT_DEBUG_WRITE('before 2nd barrier communicate_num_particles_send_recv')
     call MPI_Barrier(domain % dminfo % comm, err)
     LIGHT_DEBUG_WRITE('after 2nd barrier communicate_num_particles_send_recv')
     LIGHT_DEBUG_WRITE('finished communicate_num_particles_send_recv')
#endif

     !2. communicating from processor to processor
     ! make appropriate list
     call allocate_list_particlelists(nPartRecv, listPLRecv)

     ! communicate data (assumes that each processor has knowledge about construction of the pool lagrPartTrackPoolHalo,
     ! from the registry.  If this changes, this will break this member...  It also assumes the code is deterministic
     ! and that pools are built and computed the exact same way on each processor.
     call mpas_timer_start("communicate_data_halo_LPT")
     !call MPI_Barrier(domain % dminfo % comm, err)
     call communicate_particle_halo_data(domain, procNeighs, nPartSend, nPartRecv, listPLSend, listPLRecv)
     !call MPI_Barrier(domain % dminfo % comm, err)
     call mpas_timer_stop("communicate_data_halo_LPT")
     LIGHT_DEBUG_WRITE('finished communicate_particle_halo_data')

     ! note that at present there is no halo data
     LIGHT_DEBUG_WRITE('there is no halo data to transfer currently')
     !call mpas_timer_start("communicate_data_nonhalo_LPT")
     !if(haloOnly) then
     !  ! need to also allocate the nonHalo data somehow, it just needs initialized so that it can be "filled in" when
     !  ! necessary for output
     !  ! can utilize empty 0'd fields for the nonhalo data portion to initialize the field for all the processors
     !  call allocate_list_nonHalo_data(domain, listPLRecv)
     !else
     !  !call communicate_particle_nonhalo_data(domain, procNeighs, nPartSend, nPartRecv, listPLSend, listPLRecv)
     !end if
     !call mpas_timer_stop("communicate_data_nonhalo_LPT")

     ! now there should be the complete particles on listPLRecv.  These, however, need moved to the blocks of the processor
     ! for use in calculations
     !3. communicating from block to block
     call distribute_particlelist_to_blocks(domain, namedBlock, listPLRecv)

     LIGHT_DEBUG_WRITE('finished distributing particlelist')
     LIGHT_DEBUG_WRITE('halo procs = ' COMMA procNeighs)
     LIGHT_DEBUG_WRITE('nSend = ' COMMA nPartSend COMMA ' nRecv = ' COMMA nPartRecv)

     ! deallocate listPLSend and listPLRecv
     if (copyOnly) then
       ! just empty the list without destroying the particle data
       !call mpas_log_write( 'just emptying the particlelist')
       call empty_list_particlelists(listPLSend)
     else
       ! remove list of particles as well as particle data because it was just sent to the other processors
       call destory_list_particlelists(listPLSend)
     endif
     ! just empty the list without destroying the particle data (because we need particles that were just transferred!)
     call empty_list_particlelists(listPLRecv)

     deallocate(nPartSend, nPartRecv)

     LIGHT_DEBUG_WRITE('finished mpas_particle_list_test_num_current_particlelist')

   end subroutine mpas_particle_list_transfer_particles_from_block_to_named_block !}}}

!***********************************************************************
!
!  routine mpas_particle_list_write_halo_data
!
!> \brief   Writes haloData to struct arrays for output
!> \author  Phillip Wolfram
!> \date    06/03/2014
!> \details
!>  This routine writes haloData output for this MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------
   subroutine mpas_particle_list_write_halo_data(domain, err)!{{{

      implicit none

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_particle_list_type), pointer :: particlelist
      type (mpas_pool_type), pointer :: lagrPartTrackPool
      type (mpas_pool_iterator_type) :: dimItr
      type (field1DReal), pointer :: field1DRealPointer
      type (field1DInteger), pointer :: field1DIntPointer
      real (kind=RKIND), dimension(:), pointer :: Array1DRealPointer => NULL()
      integer, dimension(:), pointer :: Array1DIntPointer => NULL()
      integer, dimension(:), pointer :: indexToParticleIDOriginal => NULL(), &
                                        indexToParticleIDNew => NULL(), &
                                        orderingVector => NULL()
      character (len=StrKIND) :: message

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! particle related pointers
        particlelist => block % particlelist
        ! iterate over each member of the pool and make the relevant assignment
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackHalo', lagrPartTrackPool)
        ! need to compute the ordering matrices
        call mpas_pool_get_array(lagrPartTrackPool, 'indexToParticleID', indexToParticleIDOriginal)
        call get_halo_data_from_particle_list_array(particlelist, 'indexToParticleID', indexToParticleIDNew)
        ! note: orderingVector can be a subset of indexToParticleIDNew because this index can include compute as well
        !       as IO particles however, it must be of the same size as indexToParticleIDOriginal
        call compute_ordering_vector(indexToParticleIDOriginal, indexToParticleIDNew, orderingVector)
              LIGHT_DEBUG_WRITE('write halo data')
              LIGHT_DEBUG_WRITE('indexToParticleIDOriginal =' COMMA indexToParticleIDOriginal)
              LIGHT_DEBUG_WRITE('indexToParticleIDNew =' COMMA indexToParticleIDNew)
              LIGHT_DEBUG_WRITE('ordering vector=' COMMA orderingVector)

        allocate(Array1DRealPointer(count_particlelist(particlelist)))
        allocate(Array1DIntPointer(count_particlelist(particlelist)))

        ! iterate over contents of pool and transfer
        call mpas_pool_begin_iteration(lagrPartTrackPool)
        do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
          ! determine the type of data
          if (dimItr % memberType == MPAS_POOL_FIELD) then
            if (dimItr % dataType == MPAS_POOL_REAL) then
              ! get data and place it in appropriate array
              call mpas_pool_get_field(lagrPartTrackPool, dimItr % memberName, field1DRealPointer)
                !{{{
              LIGHT_DEBUG_WRITE('write halo data')
              LIGHT_DEBUG_WRITE('member name =' COMMA trim(dimItr % memberName))
              LIGHT_DEBUG_WRITE('particlelistSize= ' COMMA  count_particlelist(particlelist))
              LIGHT_DEBUG_WRITE('memory arraysize= ' COMMA  size(field1DRealPointer % array))
              LIGHT_DEBUG_WRITE(trim(dimItr % memberName) COMMA ' = ' COMMA field1DRealPointer % array)
              !}}}
              call get_halo_data_from_particle_list_array(particlelist, dimItr % memberName, Array1DRealPointer)
              ! reorder
              field1DRealPointer % array = Array1DRealPointer(orderingVector)
            elseif (dimItr % dataType == MPAS_POOL_INTEGER) then
              ! get data and place it in appropriate array
              call mpas_pool_get_field(lagrPartTrackPool, dimItr % memberName, field1DIntPointer)
              call get_halo_data_from_particle_list_array(particlelist, dimItr % memberName, Array1DIntPointer)
              ! reorder
              field1DIntPointer % array = Array1DIntPointer(orderingVector)
            else
              LIGHT_DEBUG_ALL_WRITE("Different field type than implemented in Halo write!")
            end if
          elseif (dimItr % memberType == MPAS_POOL_DIMENSION) then
            ! ignore dimensions for now and have this code so they aren't printed as an error message
          else
            write(message, *) "Different type expected in registry for key ", trim(dimItr % memberName), &
                              " in Halo data for write-- don't know what to do!"
            LIGHT_DEBUG_ALL_WRITE(message)
          end if
        end do

        ! free memory for the next loop
        deallocate(indexToParticleIDNew)
        deallocate(orderingVector)
        deallocate(Array1DRealPointer)
        deallocate(Array1DIntPointer)

        block => block % next
      end do

   end subroutine mpas_particle_list_write_halo_data!}}}

!***********************************************************************
!
!  routine mpas_particle_list_write_nonhalo_data
!
!> \brief   Writes nonhaloData to struct arrays for output
!> \author  Phillip Wolfram
!> \date    06/03/2014
!> \details
!>  This routine writes nonhaloData output for this MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------
!  subroutine mpas_particle_list_write_nonhalo_data(domain, err)!{{{

!     implicit none

!     !-----------------------------------------------------------------
!     ! input variables
!     !-----------------------------------------------------------------
!     type (domain_type), intent(in) :: domain

!     !-----------------------------------------------------------------
!     ! output variables
!     !-----------------------------------------------------------------
!     integer, intent(out) :: err !< Output: error flag

!     !-----------------------------------------------------------------
!     ! local variables
!     !-----------------------------------------------------------------
!     type (block_type), pointer :: block
!     type (mpas_particle_list_type), pointer :: particlelist
!     type (mpas_pool_type), pointer :: lagrPartTrackPool
!     type (mpas_pool_iterator_type) :: dimItr
!     type (field1DReal), pointer :: field1DRealPointer
!     type (field1DInteger), pointer :: field1DIntPointer
!     real (kind=RKIND), dimension(:), pointer :: Array1DRealPointer => NULL()
!     integer, dimension(:), pointer :: Array1DIntPointer => NULL()
!     integer, dimension(:), pointer :: indexToParticleIDOriginal => NULL(), &
!                                       indexToParticleIDNew => NULL(), &
!                                       orderingVector =>NULL()
!     character (len=StrKIND) :: message

!     err = 0

!     block => domain % blocklist
!     do while (associated(block))
!       ! particle related pointers
!       particlelist => block % particlelist
!       ! iterate over each member of the pool and make the relevant assignment
!       call mpas_pool_get_subpool(block % structs, 'lagrPartTrackHalo', lagrPartTrackPool)
!       ! need to compute the ordering matrices
!       call mpas_pool_get_array(lagrPartTrackPool, 'indexToParticleID', indexToParticleIDOriginal)
!       call get_halo_data_from_particle_list_array(particlelist, 'indexToParticleID', indexToParticleIDNew)
!       ! note: orderingVector can be a subset of indexToParticleIDNew because this index can include compute as well as
!       !       IO particles however, it must be of the same size as indexToParticleIDOriginal
!       call compute_ordering_vector(indexToParticleIDOriginal, indexToParticleIDNew, orderingVector)

!       ! iterate over each member of the pool and make the relevant assignment
!       call mpas_pool_get_subpool(block % structs, 'lagrPartTrackNonHalo', lagrPartTrackPool)
!       call mpas_pool_begin_iteration(lagrPartTrackPool)
!       do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
!         ! determine the type of data
!         if (dimItr % memberType == MPAS_POOL_FIELD) then
!           if (dimItr % dataType == MPAS_POOL_REAL) then
!             ! get data and place it in appropriate array
!             call mpas_pool_get_field(lagrPartTrackPool, dimItr % memberName, field1DRealPointer)
!               !{{{
!             LIGHT_DEBUG_WRITE('write nonhalo data')
!             LIGHT_DEBUG_WRITE('member name =' COMMA  dimItr % memberName)
!             LIGHT_DEBUG_WRITE('particlelistSize= ' COMMA  count_particlelist(particlelist))
!             LIGHT_DEBUG_WRITE('memory arraysize= ' COMMA  size(field1DRealPointer % array))
!             !}}}
!             allocate(Array1DRealPointer(count_particlelist(particlelist)))
!             call get_nonhalo_data_from_particle_list_array(particlelist, dimItr % memberName, Array1DRealPointer)
!             ! reorder
!             field1DRealPointer % array = Array1DRealPointer(orderingVector)
!             deallocate(Array1DRealPointer)
!           elseif (dimItr % dataType == MPAS_POOL_INTEGER) then
!             write(message, *) "Integer type in registry for key ", dimItr % memberName, &
!                               " in nonHalo data for write, not yet tested!"
!             LIGHT_DEBUG_ALL_WRITE(message)
!             ! get data and place it in appropriate array
!             call mpas_pool_get_field(lagrPartTrackPool, dimItr % memberName, field1DIntPointer)
!             allocate(Array1DIntPointer(count_particlelist(particlelist)))
!             call get_nonhalo_data_from_particle_list_array(particlelist, dimItr % memberName, Array1DIntPointer)
!             ! reorder
!             field1DIntPointer % array = Array1DIntPointer(orderingVector)
!             deallocate(Array1DIntPointer)
!           else
!             LIGHT_DEBUG_ALL_WRITE("Different field type than implemented in nonHalo write!")
!           end if
!         elseif (dimItr % memberType == MPAS_POOL_DIMENSION) then
!           ! ignore dimensions for now and have this code so they aren't printed as an error message
!         else
!           write(message, *) "Different type expected in registry for key ", trim(dimItr % memberName), &
!                             " in nonHalo data for write-- don't know what to do!"
!           LIGHT_DEBUG_ALL_WRITE(message)
!         end if
!       end do

!       ! free memory for the next loop
!       deallocate(indexToParticleIDNew)
!       deallocate(orderingVector)

!       block => block % next
!       end do

!  end subroutine mpas_particle_list_write_nonhalo_data!}}}

!-----------------------------------------------------------------------
!
! PRIVATE SUBROUTINES
!
!-----------------------------------------------------------------------
!{{{

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine append_particle_to_particlelist
!
!> \brief MPAS add particle to particlelist, starting new list if
!>        not allocated
!> \author Phillip Wolfram
!> \date   06/26/2014
!> \details
!>  This routine takes a particle and places it on an existing particlelist,
!>  or in the event of no list, initializes the list with the particle
!
!-----------------------------------------------------------------------
subroutine append_particle_to_particlelist(particle, particlelist)!{{{
  implicit none
  type (mpas_particle_type), pointer, intent(in) :: particle
  type (mpas_particle_list_type), pointer, intent(inout) :: particlelist
  type (mpas_particle_list_type), pointer :: headPL=>NULL(), tempPL=>NULL()

  ! case where particeList needs to be created and populated by particle
  if(.not.associated(particlelist)) then
    allocate(particlelist)
    nullify(particlelist % next)
    nullify(particlelist % prev)
    particlelist % particle => particle
  else
    ! case where list exists get the head
    headPL => particlelist
    if (.not.associated(headPL % particle)) then
      ! populate empty link
      headPL % particle => particle
    else
      ! build new link
      allocate(tempPL)
      tempPL % particle => particle
      nullify(tempPL % prev)
      tempPL % next => headPL
      ! connect to the list
      headPL % prev => tempPL
      particlelist => tempPL
    end if
  end if

end subroutine append_particle_to_particlelist !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_get_halodata_to_particlelist_1Dreal_array
!
!> \brief MPAS get halodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   07/01/2014
!> \details
!>  This routine reads 0D real arrays in each particle of the particlelist
!>  and places the data into a 1D real field output.
!
!-----------------------------------------------------------------------
subroutine get_halo_data_from_particle_list_1Dreal_array &  !{{{
    (particlelist, dataName, array1DRealPointer)
    ! input data
    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
    character(len=*), intent(in) :: dataName
    real (kind=RKIND), dimension(:), pointer, intent(out) :: array1DRealPointer

    ! subroutine data
    integer :: dataNumber
    type (mpas_particle_list_type), pointer :: particlelistCurr
    type (field0DReal), pointer :: field0DRealPointer

    ! allocate the array if it isn't allocated
    if(.not.associated(array1DRealPointer)) then
      allocate(array1DRealPointer(count_particlelist(particlelist)))
    end if

    ! loop over all elements of the list and insert the data
    dataNumber = 1
    particlelistCurr => particlelist
    ! while we have a real link
    do while(associated(particlelistCurr))

      call mpas_pool_get_field(particlelistCurr % particle % haloDataPool, dataName, field0DRealPointer)
      array1DRealPointer(dataNumber) = field0DRealPointer % scalar

      ! increment for new dataNumber
      dataNumber = dataNumber + 1
      ! get next link
      particlelistCurr => particlelistCurr % next
    end do

end subroutine get_halo_data_from_particle_list_1Dreal_array !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_get_halodata_to_particlelist_1Dreal
!
!> \brief MPAS get halodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   06/03/2014
!> \details
!>  This routine reads 0D real arrays in each particle of the particlelist
!>  and places the data into a 1D real field output.
!
!-----------------------------------------------------------------------
subroutine get_halo_data_from_particle_list_1Dreal &  !{{{
    (particlelist, dataName, field1DRealPointer)
    ! input data
    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
    character(len=*), intent(in) :: dataName
    type (field1DReal), pointer, intent(out) :: field1DRealPointer

    ! subroutine data
    integer :: dataNumber
    type (mpas_particle_list_type), pointer :: particlelistCurr
    type (field0DReal), pointer :: field0DRealPointer

    ! loop over all elements of the list and insert the data
    dataNumber = 1
    particlelistCurr => particlelist
    ! while we have a real link
    do while(associated(particlelistCurr))

      call mpas_pool_get_field(particlelistCurr % particle % haloDataPool, dataName, field0DRealPointer)
      field1DRealPointer % array(dataNumber) = field0DRealPointer % scalar

      ! increment for new dataNumber
      dataNumber = dataNumber + 1
      ! get next link
      particlelistCurr => particlelistCurr % next
    end do

end subroutine get_halo_data_from_particle_list_1Dreal !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_get_nonhalodata_to_particlelist_1Dint
!
!> \brief MPAS get nonhalodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   07/07/2014
!> \details
!
!-----------------------------------------------------------------------
!subroutine get_nonhalo_data_from_particle_list_1Dint&  !{{{
!    (particlelist, dataName, field1DIntPointer)
!    ! input data
!    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
!    character(len=*), intent(in) :: dataName
!    type (field1DInteger), pointer, intent(out) :: field1DIntPointer
!
!    ! subroutine data
!    integer :: dataNumber
!    type (mpas_particle_list_type), pointer :: particlelistCurr
!    type (field0DInteger), pointer :: field0DIntPointer
!
!    ! loop over all elements of the list and insert the data
!    dataNumber = 1
!    particlelistCurr => particlelist
!    ! while we have a real link
!    do while(associated(particlelistCurr))
!      call mpas_pool_get_field(particlelistCurr % particle % nonhaloDataPool, dataName, field0DIntPointer)
!      field1DIntPointer % array(dataNumber) = field0DIntPointer % scalar
!
!      ! increment for new dataNumber
!      dataNumber = dataNumber + 1
!      ! get next link
!      particlelistCurr => particlelistCurr % next
!    end do
!
!end subroutine get_nonhalo_data_from_particle_list_1Dint !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_get_nonhalodata_to_particlelist_1Dint_array
!
!> \brief MPAS get nonhalodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   07/07/2014
!> \details
!
!-----------------------------------------------------------------------
!subroutine get_nonhalo_data_from_particle_list_1Dint_array &  !{{{
!    (particlelist, dataName, array1DIntPointer)
!    ! input data
!    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
!    character(len=*), intent(in) :: dataName
!    integer, dimension(:), pointer, intent(out) :: array1DIntPointer
!
!    ! subroutine data
!    integer :: dataNumber
!    type (mpas_particle_list_type), pointer :: particlelistCurr
!    type (field0DInteger), pointer :: field0DIntPointer
!
!    ! allocate the array if it isn't allocated
!    if(.not.associated(array1DIntPointer)) then
!      allocate(array1DIntPointer(count_particlelist(particlelist)))
!    end if
!
!    ! loop over all elements of the list and insert the data
!    dataNumber = 1
!    particlelistCurr => particlelist
!    ! while we have a real link
!    do while(associated(particlelistCurr))
!      call mpas_pool_get_field(particlelistCurr % particle % nonhaloDataPool, dataName, field0DIntPointer)
!      array1DIntPointer(dataNumber) = field0DIntPointer % scalar
!
!      ! increment for new dataNumber
!      dataNumber = dataNumber + 1
!      ! get next link
!      particlelistCurr => particlelistCurr % next
!    end do
!
!end subroutine get_nonhalo_data_from_particle_list_1Dint_array !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_get_nonhalodata_to_particlelist_1Dreal
!
!> \brief MPAS get nonhalodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   06/04/2014
!> \details
!
!-----------------------------------------------------------------------
!subroutine get_nonhalo_data_from_particle_list_1Dreal &  !{{{
!    (particlelist, dataName, field1DRealPointer)
!    ! input data
!    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
!    character(len=*), intent(in) :: dataName
!    type (field1DReal), pointer, intent(out) :: field1DRealPointer
!
!    ! subroutine data
!    integer :: dataNumber
!    type (mpas_particle_list_type), pointer :: particlelistCurr
!    type (field0DReal), pointer :: field0DRealPointer
!
!    ! loop over all elements of the list and insert the data
!    dataNumber = 1
!    particlelistCurr => particlelist
!    ! while we have a real link
!    do while(associated(particlelistCurr))
!      call mpas_pool_get_field(particlelistCurr % particle % nonhaloDataPool, dataName, field0DRealPointer)
!      field1DRealPointer % array(dataNumber) = field0DRealPointer % scalar
!
!      ! increment for new dataNumber
!      dataNumber = dataNumber + 1
!      ! get next link
!      particlelistCurr => particlelistCurr % next
!    end do
!
!end subroutine get_nonhalo_data_from_particle_list_1Dreal !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_get_nonhalodata_to_particlelist_1Dreal_array
!
!> \brief MPAS get nonhalodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   07/07/2014
!> \details
!
!-----------------------------------------------------------------------
!subroutine get_nonhalo_data_from_particle_list_1Dreal_array &  !{{{
!    (particlelist, dataName, array1DRealPointer)
!    ! input data
!    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
!    character(len=*), intent(in) :: dataName
!    real (kind=RKIND), dimension(:), pointer, intent(out) :: array1DRealPointer
!
!    ! subroutine data
!    integer :: dataNumber
!    type (mpas_particle_list_type), pointer :: particlelistCurr
!    type (field0DReal), pointer :: field0DRealPointer
!
!    ! allocate the array if it isn't allocated
!    if(.not.associated(array1DRealPointer)) then
!      allocate(array1DRealPointer(count_particlelist(particlelist)))
!    end if
!
!    ! loop over all elements of the list and insert the data
!    dataNumber = 1
!    particlelistCurr => particlelist
!    ! while we have a real link
!    do while(associated(particlelistCurr))
!      call mpas_pool_get_field(particlelistCurr % particle % nonhaloDataPool, dataName, field0DRealPointer)
!      array1DRealPointer(dataNumber) = field0DRealPointer % scalar
!
!      ! increment for new dataNumber
!      dataNumber = dataNumber + 1
!      ! get next link
!      particlelistCurr => particlelistCurr % next
!    end do
!
!end subroutine get_nonhalo_data_from_particle_list_1Dreal_array !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_get_nonhalodata_to_particlelist_2Dreal
!
!> \brief MPAS get nonhalodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   04/10/2014
!> \details
!>  This routine takes a an array of 1D real arrays and places the data into
!>  the nonhaloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
!subroutine get_nonhalo_data_from_particle_list_2Dreal &  !{{{
!    (particlelist, dataName, field2DRealPointer)
!    ! input data
!    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
!    character(len=*), intent(in) :: dataName
!    type (field2DReal), pointer, intent(out) :: field2DRealPointer
!
!    ! subroutine data
!    integer :: dataNumber
!    type (mpas_particle_list_type), pointer :: particlelistCurr
!    type (field1DReal), pointer :: field1DRealPointer
!
!    ! loop over all elements of the list and insert the data
!    dataNumber = 1
!    particlelistCurr => particlelist
!    ! while we have a real link
!    do while(associated(particlelistCurr))
!      call mpas_pool_get_field(particlelistCurr % particle % nonhaloDataPool, dataName, field1DRealPointer)
!      field2DRealPointer % array(dataNumber, :) = field1DRealPointer % array(:)
!
!      ! increment for new dataNumber
!      dataNumber = dataNumber + 1
!      ! get next link
!      particlelistCurr => particlelistCurr % next
!    end do
!
!end subroutine get_nonhalo_data_from_particle_list_2Dreal !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_get_halodata_to_particlelist_1Dint_array
!
!> \brief MPAS get halodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   04/10/2014
!> \details
!>  This routine takes a an array of 1D int arrays and places the data into
!>  the haloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
subroutine get_halo_data_from_particle_list_1Dint_array &  !{{{
    (particlelist, dataName, array1DIntPointer)
  ! input data
  type (mpas_particle_list_type), pointer, intent(in) :: particlelist
  character(len=*), intent(in) :: dataName
  integer, dimension(:), pointer, intent(out) :: array1DIntPointer

  ! subroutine data
  integer :: dataNumber
  type (mpas_particle_list_type), pointer :: particlelistCurr
  type (field0DInteger), pointer :: field0DIntPointer

  ! allocate the array if it isn't allocated
  if(.not.associated(array1DIntPointer)) then
    allocate(array1DIntPointer(count_particlelist(particlelist)))
  end if

  ! loop over all elements of the list and insert the data
  dataNumber = 1
  particlelistCurr => particlelist
  ! while we have a int link
  do while(associated(particlelistCurr))
    call mpas_pool_get_field(particlelistCurr % particle % haloDataPool, dataName, field0DIntPointer)
    array1DIntPointer(dataNumber) = field0DIntPointer % scalar

    ! increment for new dataNumber
    dataNumber = dataNumber + 1
    ! get next link
    particlelistCurr => particlelistCurr % next
  end do

end subroutine get_halo_data_from_particle_list_1Dint_array !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_get_halodata_to_particlelist_1Dint
!
!> \brief MPAS get halodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   04/10/2014
!> \details
!>  This routine takes a an array of 1D int arrays and places the data into
!>  the haloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
subroutine get_halo_data_from_particle_list_1Dint &  !{{{
    (particlelist, dataName, field1DIntPointer)
  ! input data
  type (mpas_particle_list_type), pointer, intent(in) :: particlelist
  character(len=*), intent(in) :: dataName
  type (field1DInteger), pointer, intent(out) :: field1DIntPointer

  ! subroutine data
  integer :: dataNumber
  type (mpas_particle_list_type), pointer :: particlelistCurr
  type (field0DInteger), pointer :: field0DIntPointer

  ! loop over all elements of the list and insert the data
  dataNumber = 1
  particlelistCurr => particlelist
  ! while we have a int link
  do while(associated(particlelistCurr))
    call mpas_pool_get_field(particlelistCurr % particle % haloDataPool, dataName, field0DIntPointer)
    field1DIntPointer % array(dataNumber) = field0DIntPointer % scalar

    ! increment for new dataNumber
    dataNumber = dataNumber + 1
    ! get next link
    particlelistCurr => particlelistCurr % next
  end do

end subroutine get_halo_data_from_particle_list_1Dint !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_add_halodata_to_particlelist_1Dreal_array
!
!> \brief MPAS add halodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   07/01/2014
!> \details
!>  This routine takes a an array of real scalars and places the data into
!>  the haloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
subroutine add_halo_data_to_particle_list_1Dreal_array & !{{{
    (particlelist, dataName, array1DRealPointer)
    ! input data
    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
    character(len=*), intent(in) :: dataName
    real (kind=RKIND), dimension(:), intent(in) :: array1DRealPointer

    ! subroutine data
    integer :: dataNumber
    type (mpas_particle_list_type), pointer :: particlelistCurr
    type (field0DReal), pointer :: field0DRealPointer

    ! loop over all elements of the list and insert the data
    dataNumber = 1
    particlelistCurr => particlelist
    ! while we have a real link
    do while(associated(particlelistCurr))

      allocate(field0DRealPointer)
      field0DRealPointer % scalar = array1DRealPointer(dataNumber)
      call mpas_pool_add_field(particlelistCurr % particle % haloDataPool, dataName, field0DRealPointer)

      ! increment for new dataNumber
      dataNumber = dataNumber + 1
      ! get next link
      particlelistCurr => particlelistCurr % next
    end do

end subroutine add_halo_data_to_particle_list_1Dreal_array !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_add_halodata_to_particlelist_1Dreal
!
!> \brief MPAS add halodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   06/03/2014
!> \details
!>  This routine takes a an array of real scalars and places the data into
!>  the haloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
subroutine add_halo_data_to_particle_list_1Dreal & !{{{
    (particlelist, dataName, field1DRealPointer)
    ! input data
    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
    character(len=*), intent(in) :: dataName
    type (field1DReal), pointer, intent(in) :: field1DRealPointer

    ! subroutine data
    integer :: dataNumber
    type (mpas_particle_list_type), pointer :: particlelistCurr
    type (field0DReal), pointer :: field0DRealPointer

    ! loop over all elements of the list and insert the data
    dataNumber = 1
    particlelistCurr => particlelist
    ! while we have a real link
    do while(associated(particlelistCurr))

      allocate(field0DRealPointer)
      field0DRealPointer % scalar = field1DRealPointer % array(dataNumber)
      call mpas_pool_add_field(particlelistCurr % particle % haloDataPool, dataName, field0DRealPointer)

      ! increment for new dataNumber
      dataNumber = dataNumber + 1
      ! get next link
      particlelistCurr => particlelistCurr % next
    end do

end subroutine add_halo_data_to_particle_list_1Dreal !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_add_nonhalodata_to_particlelist_1Dreal_array
!
!> \brief MPAS add nonhalodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   07/01/2014
!> \details
!>  This routine takes a an array of real scalars and places the data into
!>  the nonhaloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
!subroutine add_nonhalo_data_to_particle_list_1Dreal_array & !{{{
!    (particlelist, dataName, array1DRealPointer)
!    ! input data
!    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
!    character(len=*), intent(in) :: dataName
!    real (kind=RKIND), dimension(:), pointer :: array1DRealPointer
!
!    ! subroutine data
!    integer :: dataNumber
!    type (mpas_particle_list_type), pointer :: particlelistCurr
!    type (field0DReal), pointer :: field0DRealPointer
!
!    ! loop over all elements of the list and insert the data
!    dataNumber = 1
!    particlelistCurr => particlelist
!    ! while we have a real link
!    do while(associated(particlelistCurr))
!
!      allocate(field0DRealPointer)
!      field0DRealPointer % scalar = array1DRealPointer(dataNumber)
!      call mpas_pool_add_field(particlelistCurr % particle % nonhaloDataPool, dataName, field0DRealPointer)
!
!      ! increment for new dataNumber
!      dataNumber = dataNumber + 1
!      ! get next link
!      particlelistCurr => particlelistCurr % next
!    end do
!
!end subroutine add_nonhalo_data_to_particle_list_1Dreal_array !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_add_nonhalodata_to_particlelist_1Dreal
!
!> \brief MPAS add nonhalodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   06/03/2014
!> \details
!>  This routine takes a an array of real scalars and places the data into
!>  the nonhaloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
!subroutine add_nonhalo_data_to_particle_list_1Dreal & !{{{
!    (particlelist, dataName, field1DRealPointer)
!    ! input data
!    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
!    character(len=*), intent(in) :: dataName
!    type (field1DReal), pointer :: field1DRealPointer
!
!    ! subroutine data
!    integer :: dataNumber
!    type (mpas_particle_list_type), pointer :: particlelistCurr
!    type (field0DReal), pointer :: field0DRealPointer
!
!    ! loop over all elements of the list and insert the data
!    dataNumber = 1
!    particlelistCurr => particlelist
!    ! while we have a real link
!    do while(associated(particlelistCurr))
!
!      allocate(field0DRealPointer)
!      field0DRealPointer % scalar = field1DRealPointer % array(dataNumber)
!      call mpas_pool_add_field(particlelistCurr % particle % nonhaloDataPool, dataName, field0DRealPointer)
!
!      ! increment for new dataNumber
!      dataNumber = dataNumber + 1
!      ! get next link
!      particlelistCurr => particlelistCurr % next
!    end do
!
!end subroutine add_nonhalo_data_to_particle_list_1Dreal !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_add_halodata_to_particlelist_1Dint_array
!
!> \brief MPAS add halodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   07/01/2014
!> \details
!>  This routine takes a an array of int scalars and places the data into
!>  the haloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
subroutine add_halo_data_to_particle_list_1Dint_array & !{{{
    (particlelist, dataName, array1DIntPointer)
    ! input data
    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
    character(len=*), intent(in) :: dataName
    integer, dimension(:), intent(in) :: array1DIntPointer

    ! subroutine data
    integer :: dataNumber
    type (mpas_particle_list_type), pointer :: particlelistCurr
    type (field0DInteger), pointer :: field0DIntPointer

    ! loop over all elements of the list and insert the data
    dataNumber = 1
    particlelistCurr => particlelist
    ! while we have a int link
    do while(associated(particlelistCurr))

      allocate(field0DIntPointer)
      field0DIntPointer % scalar = array1DIntPointer(dataNumber)
      call mpas_pool_add_field(particlelistCurr % particle % haloDataPool, dataName, field0DIntPointer)

      ! increment for new dataNumber
      dataNumber = dataNumber + 1
      ! get next link
      particlelistCurr => particlelistCurr % next
    end do

end subroutine add_halo_data_to_particle_list_1Dint_array !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_add_halodata_to_particlelist_1Dint
!
!> \brief MPAS add halodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   06/03/2014
!> \details
!>  This routine takes a an array of int scalars and places the data into
!>  the haloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
subroutine add_halo_data_to_particle_list_1Dint & !{{{
    (particlelist, dataName, field1DIntPointer)
    ! input data
    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
    character(len=*), intent(in) :: dataName
    type (field1DInteger), pointer, intent(in) :: field1DIntPointer

    ! subroutine data
    integer :: dataNumber
    type (mpas_particle_list_type), pointer :: particlelistCurr
    type (field0DInteger), pointer :: field0DIntPointer

    ! loop over all elements of the list and insert the data
    dataNumber = 1
    particlelistCurr => particlelist
    ! while we have a int link
    do while(associated(particlelistCurr))

      allocate(field0DIntPointer)
      field0DIntPointer % scalar = field1DIntPointer % array(dataNumber)
      call mpas_pool_add_field(particlelistCurr % particle % haloDataPool, dataName, field0DIntPointer)

      ! increment for new dataNumber
      dataNumber = dataNumber + 1
      ! get next link
      particlelistCurr => particlelistCurr % next
    end do

end subroutine add_halo_data_to_particle_list_1Dint !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_add_nonhalodata_to_particlelist_1Dint_array
!
!> \brief MPAS add nonhalodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   07/01/2014
!> \details
!>  This routine takes a an array of int scalars and places the data into
!>  the nonhaloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
!subroutine add_nonhalo_data_to_particle_list_1Dint_array & !{{{
!    (particlelist, dataName, array1DIntPointer)
!    ! input data
!    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
!    character(len=*), intent(in) :: dataName
!    integer, dimension(:), pointer, intent(in) :: array1DIntPointer
!
!    ! subroutine data
!    integer :: dataNumber
!    type (mpas_particle_list_type), pointer :: particlelistCurr
!    type (field0DInteger), pointer :: field0DIntPointer
!
!    ! loop over all elements of the list and insert the data
!    dataNumber = 1
!    particlelistCurr => particlelist
!    ! while we have a int link
!    do while(associated(particlelistCurr))
!
!      allocate(field0DIntPointer)
!      field0DIntPointer % scalar = array1DIntPointer(dataNumber)
!      call mpas_pool_add_field(particlelistCurr % particle % nonhaloDataPool, dataName, field0DIntPointer)
!
!      ! increment for new dataNumber
!      dataNumber = dataNumber + 1
!      ! get next link
!      particlelistCurr => particlelistCurr % next
!    end do
!
!end subroutine add_nonhalo_data_to_particle_list_1Dint_array !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_add_nonhalodata_to_particlelist_1Dint
!
!> \brief MPAS add nonhalodata to particle contained in a list
!> \author Phillip Wolfram
!> \date   06/03/2014
!> \details
!>  This routine takes a an array of int scalars and places the data into
!>  the nonhaloData pool of the particles in the list.  It is assumed
!>  that the dimensionality of the field is reduced in order by one
!>  on each particle.
!
!-----------------------------------------------------------------------
!subroutine add_nonhalo_data_to_particle_list_1Dint & !{{{
!    (particlelist, dataName, field1DIntPointer)
!    ! input data
!    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
!    character(len=*), intent(in) :: dataName
!    type (field1DInteger), pointer, intent(in) :: field1DIntPointer
!
!    ! subroutine data
!    integer :: dataNumber
!    type (mpas_particle_list_type), pointer :: particlelistCurr
!    type (field0DInteger), pointer :: field0DIntPointer
!
!    ! loop over all elements of the list and insert the data
!    dataNumber = 1
!    particlelistCurr => particlelist
!    ! while we have a int link
!    do while(associated(particlelistCurr))
!
!      allocate(field0DIntPointer)
!      field0DIntPointer % scalar = field1DIntPointer % array(dataNumber)
!      call mpas_pool_add_field(particlelistCurr % particle % nonhaloDataPool, dataName, field0DIntPointer)
!
!      ! increment for new dataNumber
!      dataNumber = dataNumber + 1
!      ! get next link
!      particlelistCurr => particlelistCurr % next
!    end do
!
!end subroutine add_nonhalo_data_to_particle_list_1Dint !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine count_particlelist_particles
!
!> \brief MPAS count particles in particlelist (must be allocated)
!> \author Phillip Wolfram
!> \date   04/14/2014
!> \details
!>  This routine counts number of particles in particlelist.
!
!-----------------------------------------------------------------------
integer function count_particlelist_particles(particlelist) !{{{
  implicit none
    ! input/output data
    type (mpas_particle_list_type), pointer, intent(in) :: particlelist

    ! subroutine data
    type (mpas_particle_list_type), pointer :: ParticleLinkCurr

    count_particlelist_particles = 0
    !call mpas_log_write( 'count_particlelist')
    !if(.not.associated(particlelist)) then
    !  call mpas_log_write( 'particleLinkCurr not associated')
    !  return
    !end if

    ! current link
    particleLinkCurr => particlelist

    do while (associated(particleLinkCurr))
      ! increment the list
      if (associated(particleLinkCurr % particle)) then
        count_particlelist_particles = count_particlelist_particles + 1
      end if
      ! get next item on the list
      particleLinkCurr => particleLinkCurr % next
    end do

    return

end function count_particlelist_particles !}}}


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine build_new_particlelist
!
!> \brief MPAS build list of particles
!> \author Phillip Wolfram
!> \date   06/03/2014
!> \details
!>  This routine builds up a list of empty particles ready to be populated.
!>  This is just a raw constructor for a particle list.  This assumes
!>  that a list is more than 1 link, such that nParticles >= 2.
!
!-----------------------------------------------------------------------
subroutine build_new_particlelist(nParticles, particlelist, ioBlock) !{{{
    ! input/output data
    ! number of particles
    integer, intent(in) :: nParticles
    integer, intent(in), optional :: ioBlock
    type (field0dInteger), pointer :: ioBlockfield
    type (mpas_particle_list_type), pointer, intent(inout) :: particlelist
    ! subroutine data
    integer aParticle
    type (mpas_particle_type), pointer :: particle
    type (mpas_particle_list_type), pointer :: newParticleLink
    type (mpas_particle_list_type), pointer :: ParticleLinkCurr

    integer :: counter

    if(nParticles == 0) then
      return
    end if

    ! instantiate a list of empty particles of dimension nParticles
    if(.not.associated(particlelist)) then
      allocate(particlelist)
    end if

    ! allocate memory for the new particle
    allocate(particle)
    call mpas_pool_create_pool(particle % haloDataPool)
    !call mpas_pool_create_pool(particle % nonhaloDataPool)
    if (present(ioBlock)) then
      allocate(ioBlockfield)
      ioBlockfield % scalar = ioBlock
      call mpas_pool_add_field(particle % haloDataPool, 'ioBlock', ioBlockfield)
    end if

    !! allocate start of list link (this must have already been done!
    ! assign allocated particle memory to link
    particlelist % particle => particle

    ! current link
    particleLinkCurr => particlelist

    ! add more links
    do aParticle = 2, nParticles
      ! allocate memory for the  new particle
      allocate(particle)
      call mpas_pool_create_pool(particle % haloDataPool)
      !call mpas_pool_create_pool(particle % nonhaloDataPool)
      if(present(ioBlock)) then
        allocate(ioBlockfield)
        ioBlockfield % scalar = ioBlock
        call mpas_pool_add_field(particle % haloDataPool, 'ioBlock', ioBlockfield)
      end if
      ! we already have one link so make a new one
      allocate(newParticleLink)
      ! place the particle in the list link
      newParticleLink % particle => particle
      nullify(newParticleLink % next)
      ! connect new link to current link
      newParticleLink % prev => particleLinkCurr
      ! next  link is the new link
      particleLinkCurr % next => newParticleLink
      ! reset link to last link
      particleLinkCurr => newParticleLink
    end do

end subroutine build_new_particlelist !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine count_particlelist
!
!> \brief MPAS count particles in particlelist
!> \author Phillip Wolfram
!> \date   04/14/2014
!> \details
!>  This routine counts number of particles in particlelist.
!
!-----------------------------------------------------------------------
integer function count_particlelist(particlelist) !{{{
  implicit none
    ! input/output data
    type (mpas_particle_list_type), pointer, intent(in) :: particlelist

    ! subroutine data
    type (mpas_particle_list_type), pointer :: ParticleLinkCurr

    count_particlelist = 0
    !call mpas_log_write( 'count_particlelist')
    !if(.not.associated(particlelist)) then
    !  call mpas_log_write( 'particleLinkCurr not associated')
    !  return
    !end if

    ! current link
    particleLinkCurr => particlelist

    do while (associated(particleLinkCurr))
      ! increment the list
      count_particlelist = count_particlelist + 1
      !print *,  'count_particlelist = ', count_particlelist
      ! get next item on the list
      particleLinkCurr => particleLinkCurr % next
    end do

    return

end function count_particlelist !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine destroy_particle
!
!> \brief MPAS destroy particle
!> \author Phillip Wolfram
!> \date   06/03/2014
!> \details
!>  This routine destroys a particle, deallocating its memory
!
!-----------------------------------------------------------------------
subroutine destroy_particle(particle) !{{{
  implicit none

  type (mpas_particle_type), pointer, intent(inout) :: particle

  if(associated(particle)) then
    if(associated(particle % haloDataPool)) then
      call mpas_pool_destroy_pool(particle % haloDataPool)
    end if
    !if(associated(particle % nonhaloDataPool)) then
    !  call mpas_pool_destroy_pool(particle % nonhaloDataPool)
    !end if
    deallocate(particle)
  end if

end subroutine destroy_particle !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine empty_particlelist
!
!> \brief MPAS empty particlelist
!> \author Phillip Wolfram
!> \date   06/27/2014
!> \details
!>  This routine emptys a particlelist, deallocating its memory
!>  but keeping memory of particles it contains intact
!
!-----------------------------------------------------------------------
subroutine empty_particlelist(particlelist) !{{{
  type (mpas_particle_list_type), pointer, intent(in) :: particlelist
  type (mpas_particle_list_type), pointer :: pLCurr, pLCurrTemp

  plCurr => particlelist
  do while(associated(plCurr))
    pLCurrTemp => pLCurr
    pLCurr => pLCurr % next
    deallocate(pLCurrTemp)
    ! N.B., particle contained in pLCurrTemp is not deallocated!!!
  end do

end subroutine empty_particlelist !}}}

subroutine destory_list_particlelists(listParticleLists) !{{{
  type (mpas_list_of_particle_list_type), dimension(:), pointer, intent(inout) :: listParticleLists
  integer :: i, sizeList

  if(associated(listParticleLists)) then
    sizeList = size(listParticleLists)
    do i=1,sizeList
      call mpas_particle_list_destroy_particle_list(listParticleLists(i) % list)
    end do
    deallocate(listParticleLists)
  end if

end subroutine destory_list_particlelists !}}}

subroutine empty_list_particlelists(listParticleLists) !{{{
  type (mpas_list_of_particle_list_type), dimension(:), pointer, intent(inout) :: listParticleLists
  integer :: i, sizeList

  if(associated(listParticleLists)) then
    sizeList = size(listParticleLists)
    do i=1,sizeList
      call empty_particlelist(listParticleLists(i) % list)
    end do
    deallocate(listParticleLists)
  end if

end subroutine empty_list_particlelists!}}}

!***********************************************************************
!
!  routine compute_cellOwnerBlock(domain, err)
!
!> \brief   Compute owner block arrays to specify halo ownership
!> \author  Phillip Wolfram
!> \date    06/25/2014
!> \details
!>  This routine computes the cellOwnerBlock for all cells on a block,
!>  diagnosing the block which owns each cell.  Could potentially be
!>  generalized and put in framework (probably need package variables
!>  if particles are used).
!
!-----------------------------------------------------------------------
   subroutine compute_cellOwnerBlock(domain, err) !{{{
      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: lagrPartTrackPool
      type (field1DInteger), pointer :: cellOwnerBlock

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! get pool to access cellOwnerBlock
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackCells', lagrPartTrackPool)
        ! prepare cellOwnerBlock for use in determining blockID to
        ! be passed from one processor to another
        call mpas_pool_get_field(lagrPartTrackPool, 'cellOwnerBlock', cellOwnerBlock)
        ! set cellOwnerBlock to be current block
        cellOwnerBlock % array(:) = block % blockID

        block => block % next
      end do
      ! exchange halos
      call mpas_dmpar_exch_halo_field(cellOwnerBlock)

   end subroutine compute_cellOwnerBlock !}}}

!***********************************************************************
!
!  routine compute_blockNeighs(domain, err)
!
!> \brief   Compute neighboring blocks from owner block,
!>          parsing the halo to get unique values
!> \author  Phillip Wolfram
!> \date    06/26/2014
!> \details
!>  This routine computes the neighboring blocks block % blockNeighs
!>  for each block.  Note that the size of this is not predetermined
!>  and depends on the partitioning (typically in graph.info.part.#)
!
!-----------------------------------------------------------------------
   subroutine compute_blockNeighs(domain, err) !{{{
      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: lagrPartTrackPool
      type (field1DInteger), pointer :: cellOwnerBlock

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! get pool to access cellOwnerBlock
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackCells', lagrPartTrackPool)
        ! get cellOwnerBlock
        call mpas_pool_get_field(lagrPartTrackPool, 'cellOwnerBlock', cellOwnerBlock)

        ! compute unique list obtained from cellOwnerBlock
        call uniqueIntegerList(cellOwnerBlock % array, block % blockNeighs)

        block => block % next
      end do

   end subroutine compute_blockNeighs !}}}

!***********************************************************************
!
!  routine compute_block_procNeighs(domain, err)
!
!> \brief   Compute neighboring processors from blockNeighs,
!>          getting unique values on a block
!> \author  Phillip Wolfram
!> \date    06/26/2014
!> \details
!>  This routine computes the neighboring processors corresponding to
!>  block % blockNeighs for each block.
!
!-----------------------------------------------------------------------
   subroutine compute_block_procNeighs(domain, err) !{{{
      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      !type (mpas_pool_type), pointer :: lagrPartTrackPool
      integer, dimension(:), pointer :: array
      integer :: numBlockNeighs, i

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! convert blockNeighs to procNeighs array
        numBlockNeighs = size(block % blockNeighs)
        allocate(array(numBlockNeighs))
        do i=1, numBlockNeighs
          call mpas_get_owning_proc(domain % dminfo, block % blockNeighs(i), array(i))
        end do
        ! compute unique list for processor neighbors
        call uniqueIntegerList(array, block % procNeighs)

        ! free up temporary memory
        deallocate(array)

        ! get the next block
        block => block % next
      end do

   end subroutine compute_block_procNeighs !}}}

!***********************************************************************
!
!  routine compute_procNeighs(domain, err)
!
!> \brief   Compute neighboring processors from blockNeighs,
!>          getting unique values across all blocks
!> \author  Phillip Wolfram
!> \date    06/26/2014
!> \details
!>  This routine computes the neighboring processors corresponding to
!>  block % procNeighs for each block.
!
!-----------------------------------------------------------------------
   subroutine compute_procNeighs(domain, err, procNeighs) !{{{
      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag
      integer, dimension(:), pointer, intent(out) :: procNeighs

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      integer, dimension(:), pointer :: tempIntegerArray
      integer :: localSize, totalSize, iStart

      err = 0

      totalSize = 0

      block => domain % blocklist
      do while(associated(block))
        totalSize = totalSize + size(block % procNeighs)
        block => block % next
      end do

      allocate(tempIntegerArray(totalSize))

      block => domain % blocklist
      iStart = 1
      do while(associated(block))
        localSize = size(block % procNeighs)
        tempIntegerArray(iStart:iStart+localSize-1) = block % procNeighs
        iStart = iStart + localSize
        block => block % next
      end do

      ! get unique array for complete list
      call uniqueIntegerList(tempIntegerArray, procNeighs)
      call removeValueFromIntList(procNeighs, domain % dminfo % my_proc_id)

      deallocate(tempIntegerArray)

   end subroutine compute_procNeighs !}}}


!***********************************************************************
!
!  routine removeValueFromIntList(list, removeval)
!
!> \brief   Remove value removeval from list vector
!> \author  Phillip Wolfram
!> \date    10/22/2015
!> \details
!>  This routine remotes the removeval value from the list in-place.
!
!-----------------------------------------------------------------------
   subroutine removeValueFromIntList(list, removeval) !{{{
     implicit none
     integer, dimension(:), pointer, intent(inout) :: list
     integer, intent(in) :: removeval

     integer, dimension(:), pointer :: tmparray
     integer :: nsize, i

     ! get size of new array
     nsize = 0
     do i = 1, size(list)
       if (list(i) /= removeval) then
         nsize = nsize + 1
       end if
     end do

     allocate(tmparray(nsize))
     nsize = 0
     do i = 1, size(list)
       if (list(i) /= removeval) then
         nsize = nsize + 1
         tmparray(nsize) = list(i)
       end if
     end do

     ! resize the list
     deallocate(list)
     allocate(list(nsize))

     ! transfer contents back to list
     do i = 1, nsize
       list(i) = tmparray(i)
     end do

     deallocate(tmparray)

   end subroutine removeValueFromIntList

!***********************************************************************
!
!  routine mpas_particle_list_self_union_halo_lists(self, list)
!
!> \brief   Taken union of self and list and store in self
!> \author  Phillip J. Wolfram
!> \date    10/29/2015
!> \details
!>  This routine computes the unique entries in an array using a
!>  linked list for dynamic memory storage
!
!-----------------------------------------------------------------------
   subroutine mpas_particle_list_self_union_halo_lists(self, list, nprocs, procid) !{{{
     implicit none
     integer, dimension(:), pointer, intent(inout) :: self
     integer, dimension(:), pointer, intent(in) :: list
     integer, intent(in) :: nprocs
     integer, intent(in) :: procid

     logical, dimension(:), pointer :: theunion
     integer :: i, thesum

     allocate(theunion(nprocs))
     theunion = .False.

     ! use Fortran integer indexing
     theunion(self+1) = .True.
     theunion(list+1) = .True.

     ! compute number of processors in halo
     thesum = 0
     do i=1,nprocs
       if (theunion(i)) then
         thesum = thesum + 1
       end if
     end do

     ! rebuild up self
     deallocate(self)
     allocate(self(thesum))

     ! store processor number in new list
     thesum = 0
     do i=1,nprocs
       if (theunion(i)) then
         thesum = thesum + 1
         self(thesum) = i-1
       end if
     end do

     deallocate(theunion)
     call removeValueFromIntList(self, procid)

   end subroutine mpas_particle_list_self_union_halo_lists !}}}

!***********************************************************************
!
!  routine uniqueIntegerList(array, uniqueList)
!
!> \brief   Compute unique entries in array with output in uniqueList
!> \author  Phillip Wolfram
!> \date    06/26/2014
!> \details
!>  This routine computes the unique entries in an array using a
!>  linked list for dynamic memory storage
!
!-----------------------------------------------------------------------
   subroutine uniqueIntegerList(array, uniqueList) !{{{
     implicit none
     integer, dimension(:), pointer, intent(out) :: uniqueList
     integer, dimension(:), pointer, intent(in) :: array

      type simplelist
        type (simplelist), pointer :: next => null()
        integer :: num
      end type simplelist

      type (simplelist), pointer :: listhead, templist, currlist

      integer :: nlist, i

      ! parse the simple list, looking for unique values
      ! algorithm will scale like Nblocks*NCells

      if(.not.associated(array)) then
        uniqueList => null()
        return
      end if

      ! first entry
      allocate(listhead)
      listhead % num = array(1)
      nlist = 1

      do i = 2, size(array)
        currlist => listhead
        ! check to see if value is on list
        do while(associated(currlist))
          if(currlist % num == array(i)) then
            exit
          else
            if(.not.associated(currlist % next)) then
              ! we are on the end of the list, so we should add the entry
              allocate(templist)
              templist % num = array(i)
              nlist = nlist + 1
              currlist % next => templist
            end if
            currlist => currlist % next
          end if
        end do
      end do

      ! now we have a complete, unique list so store it
      if(associated(uniqueList)) LIGHT_ERROR_WRITE('Trying to allocate uniqueList, which is already allocated!')
      allocate(uniqueList(nlist))

      currlist => listhead
      i = 1
      do while(associated(currlist))
        uniqueList(i) = currlist % num
        currlist => currlist % next
        i = i + 1
      end do

      ! deallocate linked list
      currlist => listhead
      do while(associated(currlist))
        templist => currlist % next
        deallocate(currlist)
        currlist => templist
      end do

   end subroutine uniqueIntegerList !}}}

!***********************************************************************
!
!  routine  compute_all_particle_values_unique_int
!
!> \brief   Get a unique list of values across particlelists on all blocks
!> \author  Phillip Wolfram
!> \date    07/02/2014
!> \details
!>  This routine computes a unique list of values for all the particles
!>  on the processor.
!
!-----------------------------------------------------------------------
     subroutine compute_all_particle_values_unique_int(domain, err, attrName, attrData) !{{{
      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      character(len=*) :: attrName

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag
      integer, dimension(:), pointer, intent(out) :: attrData

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      integer, dimension(:), pointer :: tempIntegerArray
      integer :: localSize, totalSize, iStart, nPart

      err = 0

      ! need to get full set of attrData from each processor and form unique list

      totalSize = 0

      block => domain % blocklist
      do while(associated(block))
        nPart = count_particlelist(block % particlelist)
        if (nPart > 0) then
          allocate(attrData(nPart))
          call get_halo_data_from_particle_list_array(block % particlelist, trim(attrName), attrData)
          totalSize = totalSize + size(attrData)
          deallocate(attrData)
        end if
        block => block % next
      end do

      if (totalSize > 0) allocate(tempIntegerArray(totalSize))

      block => domain % blocklist
      iStart = 1
      do while(associated(block))
        nPart = count_particlelist(block % particlelist)
        if (nPart > 0) then
          allocate(attrData(nPart))
          call get_halo_data_from_particle_list_array(block % particlelist, trim(attrName), attrData)
          localSize = size(attrData)
          tempIntegerArray(iStart:iStart+localSize-1) = attrData
          iStart = iStart + localSize
          deallocate(attrData)
        end if
        block => block % next
      end do

      ! get unique array for complete list
      !if(associated(tempIntegerArray)) print *,  'tempIntegerArray = ', tempIntegerArray
      call uniqueIntegerList(tempIntegerArray, attrData)

      if (associated(tempIntegerArray)) deallocate(tempIntegerArray)

     end subroutine compute_all_particle_values_unique_int !}}}

!***********************************************************************
!
!  routine make_proc_to_proc_particlelist(domain, err)
!
!> \brief   Compute neighboring processors from blockNeighs,
!>          getting unique values across all blocks and forming the lists
!> \author  Phillip Wolfram
!> \date    06/26/2014
!> \details
!>  This routine computes the lists of neighboring processors
!
!-----------------------------------------------------------------------
   subroutine make_proc_to_proc_particlelists(domain, copyOnly, blockSendToName, interProcPLArray, procNeighs, err) !{{{
      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      character(len=*), intent(in) :: blockSendToName
      logical, intent(in) :: copyOnly

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain
      integer, dimension(:), pointer, intent(in) :: procNeighs
      type (mpas_list_of_particle_list_type), dimension(:), &
                            pointer, intent(inout) :: interProcPLArray

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_particle_list_type), pointer :: particlelist, particlelisttemp, particlelisttemp2
      integer :: thisBlock
      integer, pointer :: particleBlock
      integer :: particleProc, arrayIndex
      !type (mpas_pool_type), pointer :: lagrPartTrackPool
      type (mpas_particle_type), pointer :: particle
      character (len=StrKIND) :: message

      err = 0

      LIGHT_DEBUG_WRITE('make_proc_to_proc_particlelists')

      block => domain % blocklist
      do while(associated(block))
        ! for each particle on each block
        particlelist => block % particlelist
        do while(associated(particlelist))
          LIGHT_DEBUG_ALL_WRITE('link on particlelist ' COMMA  loc(particlelist))
#ifdef MPAS_DEBUG
          if(.not.associated(particlelist % particle)) then
            LIGHT_ERROR_WRITE('particle not associated!')
          end if
#endif
          particle => particlelist % particle
#ifdef MPAS_DEBUG
          if(.not.associated(particle % haloDataPool)) then
            LIGHT_ERROR_WRITE('particle haloDataPool not associated!')
          end if
#endif
          call mpas_pool_get_array(particle % haloDataPool, blockSendToName, particleBlock)
          ! For example:
          !call mpas_pool_get_array(particle % haloDataPool, 'currentBlock', particleBlock)

          LIGHT_DEBUG_ALL_WRITE('particleBlock ' COMMA  particleBlock COMMA  ' for ' COMMA  trim(blockSendToName))
          LIGHT_DEBUG_ALL_WRITE('before mpas call')
          call mpas_get_owning_proc(domain % dminfo, particleBlock, particleProc)
#ifdef MPAS_DEBUG
          if(particleBlock /= particleProc) then
            write(message, *) 'Error if one block per proc: currentBlock = ', particleBlock, ' currentProc = ', particleProc
            LIGHT_ERROR_WRITE(message)
          end if
#endif
          LIGHT_DEBUG_WRITE('after mpas call')
          ! determine whether it belongs on thisBlock
          LIGHT_DEBUG_WRITE('myproc= ' COMMA  domain % dminfo % my_proc_id COMMA  ' particleProc= ' COMMA  particleProc)
          ! eventual support for multiple blocks
          !if(particleBlock /= block % blockID) then
          if(particleProc /= domain % dminfo % my_proc_id) then
            ! we need to move the particle to a list for export to particleProc
            ! get index for array and place particle in list at that index location
            !call mpas_log_write( 'get array index')
            arrayIndex = find_index(procNeighs, particleProc)
            if(arrayIndex == -1) LIGHT_ERROR_WRITE('Found processor is not on list of "halo" processors!')
            LIGHT_DEBUG_WRITE('ammending particle to processor index  ' COMMA  arrayIndex)
            call append_particle_to_particlelist(particle, interProcPLArray(arrayIndex)%list)
            write(message, *) 'added particle ', loc(particle), ' on link ', loc(particlelist), ' to list'
            LIGHT_DEBUG_ALL_WRITE(message)
            if (.not.copyOnly) then
              ! REMOVE PARTICLE FROM EXISTING PARTICLELIST
              ! remove particle reference from existing particle list to prevent bug / hitting null
              ! particle if particle is deallocated when interProcPLArray particlelists are destroyed
              ! next two lines mean (particlelist % prev) % next => particlelist % next, which
              ! fortran doesn't allow even though syntactically this makes perfect sense.
              ! 3 cases: head, middle, tail
              if(associated(particlelist % prev)) then
                LIGHT_DEBUG_ALL_WRITE('have previous particlelist link')
                particlelisttemp => particlelist % prev
                if (associated(particlelist % next)) then
                  LIGHT_DEBUG_ALL_WRITE('in middle of list')
                  ! case of the middle
                  particlelisttemp % next => particlelist % next
                  particlelisttemp2 => particlelisttemp
                  ! want to keep particle memory intact because their pointers were passed previously,
                  ! so just empty the list, don't destroy it and its contents
                  particlelisttemp => particlelist % next
                  particlelisttemp % prev => particlelisttemp2

                  particlelisttemp => particlelisttemp % prev
                  ! just need to remove the single link, particle memory needs to stay intact
                  !print *,  'deallocating link ', loc(particlelist)
                  ! deallocate link and get next link
                  deallocate(particlelist)
                  particlelist => particlelisttemp
                else
                  LIGHT_DEBUG_ALL_WRITE('in tail of list')
                  ! case of tail
                  nullify(particlelisttemp % next)
                  !print *,  'deallocating link ', loc(particlelist)
                  ! deallocate link and get next link
                  deallocate(particlelist)
                  ! no other links to process
                end if
              else
                if(associated(particlelist % next)) then
                  LIGHT_DEBUG_ALL_WRITE('in head of list')
                  ! case of head, set new head (assumes more than one particle)
                  particlelisttemp => particlelist % next
                  nullify(particlelisttemp % prev)
                  deallocate(particlelist)
                  block % particlelist => particlelisttemp
                  !print *,  'deallocating link ', loc(particlelist)
                  ! deallocate link and get next link
                  particlelist => particlelisttemp
                else
                  LIGHT_DEBUG_ALL_WRITE('single link ' COMMA  loc(particlelist))
                  LIGHT_DEBUG_ALL_WRITE('block % particlelist' COMMA  loc(block % particlelist))
                  ! case of single link / particle
                  LIGHT_DEBUG_ALL_WRITE('deallocating link ' COMMA  loc(particlelist))
                  deallocate(particlelist)
                  nullify(block % particlelist) !deallocate(block % particlelist)
                  ! there is no other particle in the list (we now have 0!)
                end if
              end if
            else
              particlelist => particlelist % next
            end if
          else
            LIGHT_DEBUG_ALL_WRITE('just go to the next particle')
            particlelist => particlelist % next
          end if
        end do

        ! this is done for each block because we want processor - processor communication
        block => block % next
      end do

   end subroutine make_proc_to_proc_particlelists!}}}

   subroutine get_num_particlelists(particlelists, numLists, nPartList) !{{{
      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_list_of_particle_list_type), dimension(:), &
                            pointer, intent(in) :: particlelists
      integer, intent(in) :: numLists

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------
      integer, dimension(:), pointer, intent(inout) :: nPartList

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer :: i, allerr=0

      ! now get number of particles in each list
      ! for each entry in list of particlelist
      nPartList = -1
      do i=1, numLists
        nPartList(i) = count_particlelist(particlelists(i)%list)
        LIGHT_DEBUG_WRITE('nPartList(i)=' COMMA  nPartList(i))
      end do

   end subroutine get_num_particlelists !}}}

!***********************************************************************
!
!  routine communicate_num_particles_send_recv
!
!> \brief   Communicate the number of particles to be sent/recv'd
!>          from adjacent processors
!> \author  Phillip Wolfram
!> \date    06/27/2014
!> \details
!>  This routine transmitts nPartSend to be stored in nPartRecv of
!>  adjacent processors in procNeighs.  Not that matrix with
!>  procNeighs in columns stored in rows for each processor
!>  must be symmetric for this to work.
!
!-----------------------------------------------------------------------
   subroutine communicate_num_particles_send_recv(domain, procNeighs, nPartSend, nPartRecv) !{{{
     implicit none

     type (domain_type), intent(in) :: domain
     integer, dimension(:), pointer, intent(in) :: procNeighs
     integer, dimension(:), pointer, intent(in) :: nPartSend
     integer, dimension(:), pointer, intent(inout) :: nPartRecv

     integer :: i, j, numProcs
     integer, dimension(:), pointer :: requestID
     integer :: mpi_ierr
     character (len=StrKIND) :: message

     numProcs = size(procNeighs)

     ! set to be -1 for error catching
     nPartRecv = -1

     ! want to send individual values in nPartSend to each entry in procNeighs (paired data)
     ! values obtained after communication are to be stored in nPartRecv

#ifdef MPAS_DEBUG
     !i = 990 + domain % dminfo % my_proc_id
     !write(message,*) 'proc_debug_message_', domain % dminfo % my_proc_id
     !open(unit=i, file=message)
     !write(i,*) 'procNeighs=', procNeighs
     !write(i,*) '\n'
     !write(i,*) 'nPartSend=', nPartSend
     !close(unit=i)

#ifdef _MPI
     call MPI_Barrier(domain % dminfo % comm, mpi_ierr)
#endif
#endif

     allocate(requestID(2*numProcs))
#ifdef MPAS_DEBUG
     LIGHT_DEBUG_WRITE('before 1st barrier')
#ifdef _MPI
     call MPI_Barrier(domain % dminfo % comm, mpi_ierr)
#endif
     LIGHT_DEBUG_WRITE('after 1st barrier')
     LIGHT_DEBUG_WRITE('receiving data')
     LIGHT_DEBUG_WRITE('numProcs =' COMMA  numProcs)
#endif
     do i=1,numProcs
       LIGHT_DEBUG_WRITE('receiving data i=' COMMA i COMMA  ' procNeighs(i)=' COMMA  procNeighs(i))
#ifdef _MPI
         call MPI_IRecv(nPartRecv(i), 1, MPI_INTEGERKIND, procNeighs(i), procNeighs(i), &
           domain % dminfo % comm, requestID(i), mpi_ierr)
#endif
     end do
     do i=1,numProcs
       LIGHT_DEBUG_WRITE('procNeighs=' COMMA  procNeighs)
       write(message, *) 'sending data i=', i,  ' procNeighs(i)=',  procNeighs(i),  ' nPartSend(i)=',  nPartSend(i)
       LIGHT_DEBUG_WRITE(message)
#ifdef _MPI
         call MPI_ISend(nPartSend(i), 1, MPI_INTEGERKIND, procNeighs(i), domain % dminfo % my_proc_id, &
           domain % dminfo % comm, requestID(numProcs + i), mpi_ierr)
#endif
     end do
     LIGHT_DEBUG_WRITE('done with send receive calls')
#ifdef _MPI
     call MPI_WaitAll(2*numProcs, requestID, MPI_STATUSES_IGNORE, mpi_ierr)
#endif
#ifdef MPAS_DEBUG
     LIGHT_DEBUG_WRITE('before 2nd barrier')
#ifdef _MPI
     call MPI_Barrier(domain % dminfo % comm, mpi_ierr)
#endif
     LIGHT_DEBUG_WRITE('after 2nd barrier')
     LIGHT_DEBUG_WRITE('done with wait')
#endif
     deallocate(requestID)

   end subroutine communicate_num_particles_send_recv !}}}

!***********************************************************************
!
!  routine compute_additional_particle_send_list
!
!> \brief   Update send list
!> \author  Phillip Wolfram
!> \date    06/24/2015
!> \details
!>  Compute send list for all particles residing on correct block 'currentBlock'
!-----------------------------------------------------------------------
   subroutine compute_additional_particle_send_list(domain, destinationName, ioProcSendList) !{{{
     implicit none
     type (domain_type), intent(in) :: domain
     character(len=*), intent(in) :: destinationName
     logical, dimension(:), pointer, intent(inout) :: ioProcSendList

     ! local variables
     type (block_type), pointer :: block
     type (mpas_particle_list_type), pointer :: particlelist
     type (mpas_particle_type), pointer :: particle
     integer, pointer :: sendBlock
     integer :: ioProc

     block => domain % blocklist
     do while (associated(block)) !{{{
       particlelist => block % particlelist
       do while(associated(particlelist)) !{{{
         particle => particlelist % particle
         call mpas_pool_get_array(particle % haloDataPool, destinationName, sendBlock)
         call mpas_get_owning_proc(domain % dminfo, sendBlock, ioProc)
         ioProcSendList(ioProc+1) = .True.
         ! get next particle to process on the list
         particlelist => particlelist % next
       end do !}}}
       ! get next block
       block => block % next
     end do !}}}

   end subroutine compute_additional_particle_send_list !}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine allocate_list_particlelists
!
!> \brief MPAS allocate list of particlelists
!> \author Phillip Wolfram
!> \date   07/01/2014
!> \details
!>  This routine allocates a list of particlelists.  Useful in preparation
!>  to receive data from MPI communication
!
!-----------------------------------------------------------------------
subroutine allocate_list_particlelists(nPart, listPL) !{{{
  implicit none

  type (mpas_list_of_particle_list_type), dimension(:), intent(inout), pointer :: listPL
  integer, dimension(:), pointer, intent(in) :: nPart

  integer :: nProcs, i

  nProcs = size(nPart)

  ! allocate size of particeLists
  do i=1,nProcs
    call build_new_particlelist(nPart(i), listPL(i)%list)
  end do

end subroutine allocate_list_particlelists !}}}

!***********************************************************************
!
!  routine distribute_particlelist_to_blocks
!
!> \brief   Take a list of particlelists and move them onto the appropriate
!>          block on the processor
!> \author  Phillip Wolfram
!> \date    07/01/2014
!> \details
!>  This routine places particles in a list of particlelist on the appropriate
!>  block (currentBlock) on the processor
!
!-----------------------------------------------------------------------
   subroutine distribute_particlelist_to_blocks(domain, blockSendToName, listPL) !{{{
     implicit none

     type (domain_type), intent(in) :: domain
     type (mpas_list_of_particle_list_type), dimension(:), pointer, intent(in) :: listPL
     character(len=*), intent(in) :: blockSendToName

     integer :: i, numLists
     type (block_type), pointer :: block
     type (mpas_particle_list_type), pointer :: tempPL
     type (mpas_particle_type), pointer :: particle
     type (mpas_pool_type), pointer :: lagrPartTrackPool
     integer, pointer :: blockNum

     ! need to copy pointers to particles from each entry of the listPL and move the
     ! particle to the appropriate block

     ! loop over each part of the particle list
     numLists = size(listPL)
     do i=1,numLists
       ! get a single list of particles
       tempPL => listPL(i) % list
       ! iterate through the list and assign particles
       do while(associated(tempPL))
         ! get the particle blockNum
         particle => tempPL % particle
         call mpas_pool_get_array(particle % haloDataPool, blockSendToName, blockNum)

         ! determine the block the particle should be moved to
         block => domain % blocklist
         LIGHT_DEBUG_ALL_WRITE('blockNum = '  COMMA  blockNum)
         blocksearch: do while(associated(block))
           !print *,  'blockID = ', block % blockID, 'blockNum = ', blockNum
           if (block % blockID == blockNum) then
             ! we have the correct block
             exit blocksearch
           end if
           block => block % next
         end do blocksearch

         ! check to make sure block is correct
         LIGHT_DEBUG_ALL_WRITE('blockNum = ' COMMA  blockNum)
         if (blockNum /= block % blockID) then
           LIGHT_DEBUG_WRITE('block is not correct! for blockNum =' COMMA blockNum COMMA ' and blockID = ' COMMA block % blockID)
         end if

         ! move the particle to the block
         call append_particle_to_particlelist(particle, block % particlelist)

         ! get the next particle in the list
         tempPL => tempPL % next
       end do

     end do

   end subroutine distribute_particlelist_to_blocks !}}}

!***********************************************************************
!
!  routine compute_ordering_vector
!
!> \brief   Compute the orderingVector for restructuring of mixed up particles
!>          into original order
!> \author  Phillip Wolfram
!> \date    07/03/2014
!> \details
!>  This routine ensures that the particles are ordered consistently with
!>  the ordering of the original data.
!
!-----------------------------------------------------------------------
   subroutine compute_ordering_vector(arrayOrig, arrayNew, orderingVector) !{{{
     implicit none

     integer, dimension(:), pointer, intent(in) :: arrayOrig, arrayNew
     integer, dimension(:), pointer, intent(out) :: orderingVector

     integer :: i, idx

     allocate(orderingVector(size(arrayOrig)))
     ! allocate to -1 to that if a value isn't found it isn't random data, assuming arrays
     ! are all indexed starting from 1
     orderingVector = -1
     ! loop over each element and figure out index
     do i=1,size(arrayOrig)
       !N.B. should in principle check to make sure index is found
       idx = find_index(arrayNew, arrayOrig(i))
#ifdef MPAS_DEBUG
       if ( idx  > 0 ) then
#endif
         orderingVector(i) = idx
#ifdef MPAS_DEBUG
       else
         LIGHT_ERROR_WRITE("ERROR! Didn't find correct ordering index")
       end if
#endif
     end do

   end subroutine compute_ordering_vector !}}}

!***********************************************************************
!
!  routine  find_index
!
!> \brief   find the index corresponding to num in array
!> \author  Phillip Wolfram
!> \date    07/03/2014
!> \details
!>  This routine returns the index such that array(find_index) = num.
!>  If the index doesn't exist it returns -1.
!
!-----------------------------------------------------------------------
   integer function find_index(array, num) !{{{
     implicit none

     integer, dimension(:), pointer, intent(in) :: array
     integer, intent(in) :: num

     integer :: i

     ! allocate to negative number to make sure it breaks if
     ! an index is not found

     find_index = -1
     ! an error here could be caused by running with different
     ! processors (blocks) specified in the input file
     ! than run with MPI
     do i=1,size(array)
       if(num == array(i)) then
         find_index = i
         return
       end if
     end do

#ifdef MPAS_DEBUG
     if(find_index == -1) LIGHT_WARNING_WRITE('Error: Index number ' COMMA num COMMA ' not found in array!')
#endif

   end function find_index !}}}

!***********************************************************************
!
!  routine communicate_particle_nonhalo_data
!
!> \brief   Communicate particle nonhalo data to relevant processors
!>          based on particlelists
!> \author  Phillip Wolfram
!> \date    07/07/2014
!> \details
!>  This routine transmitts nonHaloData from particlelists
!
!-----------------------------------------------------------------------
!  subroutine communicate_particle_nonhalo_data(domain, procNeighs, nPartSend, nPartRecv, listSend, listRecv) !{{{
!    implicit none

!    type (domain_type), intent(in) :: domain
!    integer, dimension(:), pointer, intent(in) :: procNeighs
!    integer, dimension(:), pointer, intent(in) :: nPartSend
!    integer, dimension(:), pointer, intent(in) :: nPartRecv
!    type (mpas_list_of_particle_list_type), dimension(:), pointer :: listSend, listRecv

!    integer :: i, j, numProcs, numFields, numRecv, numSends
!    integer, dimension(:), pointer :: recvRequestID, sendRequestID
!    integer :: mpi_ierr
!    type (mpas_pool_type), pointer :: lagrPartTrackPool
!    type (mpas_pool_iterator_type) :: dimItr
!    type array1DReal_list
!      real (kind=RKIND), dimension(:), pointer :: val
!    end type
!    type array1DInt_list
!      integer, dimension(:), pointer :: val
!    end type
!    type (array1DInt_list), dimension(:), pointer :: array1DIntSend, array1DIntRecv
!    type (array1DReal_list), dimension(:), pointer :: array1DRealSend, array1DRealRecv
!    character (len=StrKIND) :: message

!    ! for each entry in the halo pool, want to send and recv the data

!ifdef _MPI
!    !call MPI_Barrier(domain % dminfo % comm)
!endif

!    numProcs = size(procNeighs)
!    allocate(array1DRealSend(numProcs), array1DRealRecv(numProcs))
!    allocate(array1DIntSend(numProcs),  array1DIntRecv(numProcs))

!    numSends = 0
!    do i = 1, numProcs
!      if (nPartSend(i) > 0) numSends = numSends + 1
!    end do
!    allocate(sendRequestID(numSends))

!    numRecv = 0
!    do i = 1, numProcs
!      if (nPartRecv(i) > 0) numRecv = numRecv + 1
!    end do
!    allocate(recvRequestID(numRecv))

!    !Notes !{{{
!    !! get number of items that need transfered from halo pool, numFields which is a constant
!    !call mpas_pool_get_subpool(domain % blocklist % structs, 'lagrPartTrackHalo', lagrPartTrackPool)
!    !call mpas_pool_begin_iteration(lagrPartTrackPool)
!    !numFields = 0
!    !do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
!    !  ! only need to transfer pool
!    !  if (dimItr % memberType == MPAS_POOL_FIELD) then
!    !    numFields = numFields + 1
!    !  end if
!    !end do
!    !! assume, for now, that this will be constant accross processors.  If not, it would need to be sent to other
!    !! processors too.  This also presumes that properties will be fixed accross the processesors.
!    !}}}

!    ! on each list, transmit relevant fields to associated processors (note using the var struct since it has the names
!    ! required and this information is on each processor, even if the pool's fields are empty their names and types
!    ! are there from the registry).
!    call mpas_pool_get_subpool(domain % blocklist % structs, 'lagrPartTrackNonHalo', lagrPartTrackPool)
!    call mpas_pool_begin_iteration(lagrPartTrackPool)
!    do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
!      if (dimItr % memberType == MPAS_POOL_FIELD) then
!        !print *,  'transfering ', trim(dimItr % memberName)
!        if (dimItr % dataType == MPAS_POOL_REAL) then
!          ! recv
!          j = 1
!          do i=1,numProcs
!            if(nPartRecv(i) > 0) then
!              allocate(array1DRealRecv(i)%val(nPartRecv(i)))
!             !print *,  'receiving real ', trim(dimItr % memberName), ' from ', procNeighs(i)
!              ! receive communicated data
!ifdef _MPI
!              call MPI_IRecv(array1DRealRecv(i)%val, nPartRecv(i), MPI_REALKIND, procNeighs(i), procNeighs(i), &
!                domain % dminfo % comm, recvRequestID(j), mpi_ierr)
!endif
!              j = j + 1
!            end if
!          end do
!          ! send
!          j = 1
!          do i=1,numProcs
!            if (nPartSend(i) > 0) then
!              allocate(array1DRealSend(i)%val(nPartSend(i)))
!              !print *,  'sending real ', trim(dimItr % memberName), ' to ', procNeighs(i)
!              call get_nonhalo_data_from_particle_list_array(listSend(i)%list, dimItr % memberName, &
!                array1DRealSend(i)%val)
!ifdef _MPI
!              call MPI_ISend(array1DRealSend(i)%val, nPartSend(i), MPI_REALKIND, procNeighs(i), &
!                domain % dminfo % my_proc_id, domain % dminfo % comm, sendRequestID(j), mpi_ierr)
!endif
!              j = j + 1
!            end if
!          end do

!ifdef _MPI
!         call MPI_WaitAll(numRecv, recvRequestID, MPI_STATUSES_IGNORE, mpi_ierr)
!endif
!         if (mpi_ierr /= 0) call mpas_log_write('recv: mpi_ierr = ' COMMA mpi_ierr)
!ifdef _MPI
!         call MPI_WaitAll(numSends, sendRequestID, MPI_STATUSES_IGNORE, mpi_ierr)
!endif
!         if (mpi_ierr /= 0) call mpas_log_write('send: mpi_ierr = ' COMMA mpi_ierr)

!         ! store values
!         j = 1
!         do i=1,numProcs
!           if(nPartRecv(i) > 0) then
!             ! place it in particle list
!             call add_nonhalo_data_to_particle_list_array(listRecv(i)%list, dimItr % memberName, &
!               array1DRealRecv(i)%val)
!             j = j + 1
!           end if
!         end do

!          do i=1,numProcs
!            if(nPartSend(i) > 0) deallocate(array1DRealSend(i)%val)
!          end do
!          do i=1,numProcs
!            if(nPartRecv(i) > 0) deallocate(array1DRealRecv(i)%val)
!          end do
!         !call mpas_log_write( 'finished'

!       elseif (dimItr % dataType == MPAS_POOL_INTEGER) then
!         ! recv
!         j = 1
!         do i=1,numProcs
!           if(nPartRecv(i) > 0) then
!             allocate(array1DIntRecv(i)%val(nPartRecv(i)))
!             ! receive communicated data
!             !print *,  'receiving int ', trim(dimItr % memberName), ' from ', procNeighs(i)
!ifdef _MPI
!             call MPI_IRecv(array1DIntRecv(i)%val, nPartRecv(i), MPI_INTEGERKIND, procNeighs(i), procNeighs(i), &
!               domain % dminfo % comm, recvRequestID(j), mpi_ierr)
!endif
!             !if( trim(dimItr % memberName) == 'currentBlock') print *,  'currentBlock received ', &
!             !                                                 nPartRecv(i), ' from', procNeighs(i), ' = ', array1DIntRecv(i)%val
!             j = j + 1
!           end if
!         end do
!         ! send
!         j = 1
!         do i=1,numProcs
!           if(nPartSend(i) > 0) then
!             allocate(array1DIntSend(i)%val(nPartSend(i)))
!             !print *,  'sending int ', trim(dimItr % memberName), ' to ', procNeighs(i)
!             call get_nonhalo_data_from_particle_list_array(listSend(i)%list, dimItr % memberName, array1DIntSend(i)%val)
!             !if( trim(dimItr % memberName) == 'currentBlock') print *,  'currentBlock sent ',nPartSend(i), &
!             !                                                 ' to ', procNeighs(i), ' = ', array1DIntSend(i)%val
!ifdef _MPI
!             call MPI_ISend(array1DIntSend(i)%val, nPartSend(i), MPI_INTEGERKIND, procNeighs(i), &
!               domain % dminfo % my_proc_id, domain % dminfo % comm, sendRequestID(j), mpi_ierr)
!endif
!             j = j + 1
!           end if
!         end do

!ifdef _MPI
!         call MPI_WaitAll(numRecv, recvRequestID, MPI_STATUSES_IGNORE, mpi_ierr)
!endif
!         if (mpi_ierr /= 0) call mpas_log_write('mpi_ierr = ' COMMA mpi_ierr)
!ifdef _MPI
!         call MPI_WaitAll(numSends, sendRequestID, MPI_STATUSES_IGNORE, mpi_ierr)
!endif
!         if (mpi_ierr /= 0) print *,  'mpi_ierr = ', mpi_ierr

!         !do i=1,numProcs
!         !  if(nPartRecv(i) > 0) print *,  'Received ', trim(dimItr % memberName), ' = ', array1DIntRecv(i) % val
!         !end do

!         ! store values
!         do i=1,numProcs
!           if(nPartRecv(i) > 0) then
!             ! place it in particle list
!             call add_nonhalo_data_to_particle_list_array(listRecv(i)%list, dimItr % memberName, array1DIntRecv(i)%val)
!             j = j + 1
!           end if
!         end do

!         do i=1,numProcs
!           if(nPartSend(i) > 0) deallocate(array1DIntSend(i)%val)
!         end do
!         do i=1,numProcs
!           if(nPartRecv(i) > 0) deallocate(array1DIntRecv(i)%val)
!         end do
!         !call mpas_log_write( 'finished')
!       else
!          !call mpas_log_write( "Different field type than implemented during nonHalo communication!")
!        end if
!      elseif (dimItr % memberType == MPAS_POOL_DIMENSION) then
!        ! ignore dimensions for now and have this code so they aren't printed as an error message
!      else
!        write(message, *) "Different type expected in registry for key ", trim(dimItr % memberName), &
!                          " in nonHalo data for communication-- don't know what to do!"
!        LIGHT_DEBUG_ALL_WRITE(message)
!      end if
!    end do

!    deallocate(array1DIntSend, array1DIntRecv, array1DRealSend, array1DRealRecv, recvRequestID, sendRequestID)
!    LIGHT_DEBUG_WRITE('Finished primary MPI communication for nonhalo')

!  end subroutine communicate_particle_nonhalo_data!}}}

!***********************************************************************
!
!  routine communicate_particle_halo_data
!
!> \brief   Communicate particle halo data to relevant processors
!>          based on particlelists
!> \author  Phillip Wolfram
!> \date    07/01/2014
!> \details
!>  This routine transmitts haloData from particlelists
!
!-----------------------------------------------------------------------
   subroutine communicate_particle_halo_data(domain, procNeighs, nPartSend, nPartRecv, listSend, listRecv) !{{{
     implicit none

     type (domain_type), intent(in) :: domain
     integer, dimension(:), pointer, intent(in) :: procNeighs
     integer, dimension(:), pointer, intent(in) :: nPartSend
     integer, dimension(:), pointer, intent(in) :: nPartRecv
     integer :: MPI_TYPE
     type (mpas_list_of_particle_list_type), dimension(:), pointer :: listSend, listRecv

     integer :: i, j, ii, numProcs, numFields, numRecv, numSends, numReals, numInts
     integer, dimension(:), pointer :: recvRequestID, sendRequestID
     integer :: mpi_ierr
     type (mpas_pool_type), pointer :: lagrPartTrackPool
     type (mpas_pool_iterator_type) :: dimItr
     type array1DReal_list
       real (kind=RKIND), dimension(:), pointer :: val
     end type
     real (kind=RKIND), dimension(:), pointer :: realVal
     type array1DInt_list
       integer, dimension(:), pointer :: val
     end type
     integer, dimension(:), pointer :: intVal
     type (array1DInt_list), dimension(:), pointer :: array1DIntSend, array1DIntRecv
     type (array1DReal_list), dimension(:), pointer :: array1DRealSend, array1DRealRecv
     character (len=StrKIND) :: message

     ! for each entry in the halo pool, want to send and recv the data

#ifdef _MPI
     !call MPI_Barrier(domain % dminfo % comm)
#endif
     
     numReals = 0
     numInts = 0
     call mpas_pool_get_subpool(domain % blocklist % structs, 'lagrPartTrackHalo', lagrPartTrackPool)
     call mpas_pool_begin_iteration(lagrPartTrackPool)
     do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
       if (dimItr % memberType == MPAS_POOL_FIELD) then
         !print *,  'transfering ', trim(dimItr % memberName)
         if (dimItr % dataType == MPAS_POOL_REAL) then 
           numReals = numReals + 1
         elseif (dimItr % dataType == MPAS_POOL_INTEGER) then 
           numInts = numInts + 1
         end if
       else
         ! only need to print this warning once in this function to simplify if-else statements
         write(message, *) "Different type expected in registry for key ", trim(dimItr % memberName), &
                           " of type", dimItr % memberType, "in halo data for communication-- don't know what to do!"
         LIGHT_DEBUG_ALL_WRITE(message)
       end if 
     end do

! NOTE, code assumes that memberName is consistent between processors in terms of pool interation
!#ifdef mpas_debug
!     ! verify that membername order is the same between processors
!     call mpas_log_write('printing particle pool contents (verify that membername consistent between processors):') 
!     call mpas_pool_print_summary(lagrPartTrackPool, mpas_pool_field, .true.)
!#endif

     numProcs = size(procNeighs)

     numSends = 0
     do i = 1, numProcs
       if (nPartSend(i) > 0) numSends = numSends + 1
     end do
     allocate(sendRequestID(numSends))

     numRecv = 0
     do i = 1, numProcs
       if (nPartRecv(i) > 0) numRecv = numRecv + 1
     end do
     allocate(recvRequestID(numRecv))

     !!!!! begin REALS !!!! ! {{{
     allocate(array1DRealSend(numProcs), array1DRealRecv(numProcs))
     ! recv
     j = 1
     do i=1,numProcs ! {{{ 
       if(nPartRecv(i) > 0) then
         allocate(array1DRealRecv(i)%val(nPartRecv(i)*numReals))
#ifdef _MPI
         call MPI_IRecv(array1DRealRecv(i)%val, nPartRecv(i)*numReals, MPI_REALKIND, &
           procNeighs(i), procNeighs(i), domain % dminfo % comm, recvRequestID(j), mpi_ierr)
#endif
         j = j + 1
       end if
     end do ! }}}

     ! aggregate data to send
     do i=1,numProcs ! {{{
       if (nPartSend(i) > 0) then
         allocate(array1DRealSend(i)%val(nPartSend(i)*numReals))
         allocate(realVal(nPartSend(i)))
         j = 1
         ! perform the parallel communication
         call mpas_pool_begin_iteration(lagrPartTrackPool)
         do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
           if (dimItr % memberType == MPAS_POOL_FIELD) then
             if (dimItr % dataType == MPAS_POOL_REAL) then ! {{{
               call get_halo_data_from_particle_list_array(listSend(i)%list, dimItr % memberName, &
                 realVal)
               do ii=1,nPartSend(i)
                 array1DRealSend(i)%val((ii+(j-1)*nPartSend(i))) = realval(ii)
               end do
               j = j + 1
             end if ! }}}
           end if
         end do
         deallocate(realVal)
       end if
     end do !}}}

     ! send
     j = 1
     do i=1,numProcs ! {{{
       if (nPartSend(i) > 0) then
#ifdef _MPI
         call MPI_ISend(array1DRealSend(i)%val, nPartSend(i)*numReals, MPI_REALKIND, procNeighs(i), &
           domain % dminfo % my_proc_id, domain % dminfo % comm, sendRequestID(j), mpi_ierr)
#endif
         j = j + 1
       end if
     end do !}}}

#ifdef _MPI
     call MPI_WaitAll(numRecv, recvRequestID, MPI_STATUSES_IGNORE, mpi_ierr)
#endif
     if (mpi_ierr /= 0) call mpas_log_write('recv: mpi_ierr = %i', intArgs=(/ mpi_ierr /) )
#ifdef _MPI
     call MPI_WaitAll(numSends, sendRequestID, MPI_STATUSES_IGNORE, mpi_ierr)
#endif
     if (mpi_ierr /= 0) call mpas_log_write('send: mpi_ierr = %i', intArgs=(/ mpi_ierr /) )

     ! store values
     do i=1,numProcs ! {{{{
       if(nPartRecv(i) > 0) then
         j = 1
         call mpas_pool_begin_iteration(lagrPartTrackPool)
         do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
           if (dimItr % memberType == MPAS_POOL_FIELD) then
             if (dimItr % dataType == MPAS_POOL_REAL) then ! {{{
               ! place it in particle list
               call add_halo_data_to_particle_list_array(listRecv(i)%list, dimItr % memberName, &
                       array1DRealRecv(i)%val((1+(j-1)*nPartRecv(i)) : (j*nPartRecv(i))))
               j = j + 1
             end if ! }}}
           end if
         end do
       end if
     end do ! }}}

     do i=1,numProcs
       if(nPartSend(i) > 0) deallocate(array1DRealSend(i)%val)
     end do
     do i=1,numProcs
       if(nPartRecv(i) > 0) deallocate(array1DRealRecv(i)%val)
     end do

     deallocate(array1DRealSend, array1DRealRecv)
     !!!! end REALS !!!! ! }}}
     
     !!!!! begin INTS !!!! ! {{{
     allocate(array1DIntSend(numProcs), array1DIntRecv(numProcs))
     ! recv
     j = 1
     do i=1,numProcs ! {{{ 
       if(nPartRecv(i) > 0) then
         allocate(array1DIntRecv(i)%val(nPartRecv(i)*numInts))
#ifdef _MPI
         call MPI_IRecv(array1DIntRecv(i)%val, nPartRecv(i)*numInts, MPI_INTEGERKIND, &
           procNeighs(i), procNeighs(i), domain % dminfo % comm, recvRequestID(j), mpi_ierr)
#endif
         j = j + 1
       end if
     end do ! }}}

     ! aggregate data to send
     do i=1,numProcs ! {{{
       if (nPartSend(i) > 0) then
         allocate(array1DIntSend(i)%val(nPartSend(i)*numInts))
         allocate(intVal(nPartSend(i)))
         j = 1
         ! perform the parallel communication
         call mpas_pool_begin_iteration(lagrPartTrackPool)
         do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
           if (dimItr % memberType == MPAS_POOL_FIELD) then
             if (dimItr % dataType == MPAS_POOL_INTEGER) then ! {{{
               call get_halo_data_from_particle_list_array(listSend(i)%list, dimItr % memberName, &
                 intVal)
               do ii=1,nPartSend(i)
                 array1DIntSend(i)%val((ii+(j-1)*nPartSend(i))) = intval(ii)
               end do
               j = j + 1
             end if ! }}}
           end if
         end do
         deallocate(intVal)
       end if
     end do !}}}

     ! send
     j = 1
     do i=1,numProcs ! {{{
       if (nPartSend(i) > 0) then
#ifdef _MPI
         call MPI_ISend(array1DIntSend(i)%val, nPartSend(i)*numInts, MPI_INTEGERKIND, procNeighs(i), &
           domain % dminfo % my_proc_id, domain % dminfo % comm, sendRequestID(j), mpi_ierr)
#endif
         j = j + 1
       end if
     end do !}}}

#ifdef _MPI
     call MPI_WaitAll(numRecv, recvRequestID, MPI_STATUSES_IGNORE, mpi_ierr)
#endif
     if (mpi_ierr /= 0) call mpas_log_write('recv: mpi_ierr = %i', intArgs=(/ mpi_ierr /) )
#ifdef _MPI
     call MPI_WaitAll(numSends, sendRequestID, MPI_STATUSES_IGNORE, mpi_ierr)
#endif
     if (mpi_ierr /= 0) call mpas_log_write('send: mpi_ierr = %i', intArgs=(/ mpi_ierr /) )

     ! store values
     do i=1,numProcs ! {{{{
       if(nPartRecv(i) > 0) then
         j = 1
         call mpas_pool_begin_iteration(lagrPartTrackPool)
         do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
           if (dimItr % memberType == MPAS_POOL_FIELD) then
             if (dimItr % dataType == MPAS_POOL_INTEGER) then ! {{{
               ! place it in particle list
               call add_halo_data_to_particle_list_array(listRecv(i)%list, dimItr % memberName, &
                       array1DIntRecv(i)%val((1+(j-1)*nPartRecv(i)) : (j*nPartRecv(i))))
               j = j + 1
             end if ! }}}
           end if
         end do
       end if
     end do ! }}}

     do i=1,numProcs
       if(nPartSend(i) > 0) deallocate(array1DIntSend(i)%val)
     end do
     do i=1,numProcs
       if(nPartRecv(i) > 0) deallocate(array1DIntRecv(i)%val)
     end do

     deallocate(array1DIntSend, array1DIntRecv)
     !!!! end INTS !!!! ! }}}
     
     deallocate(recvRequestID, sendRequestID)
     LIGHT_DEBUG_WRITE('Finished primary MPI communication for halo')

   end subroutine communicate_particle_halo_data!}}}

!***********************************************************************
!
!  routine allocate_nonHalo_data
!
!> \brief   Allocate space for nonHalo data for diagnostic output
!>          on particlelists
!> \author  Phillip Wolfram
!> \date    07/01/2014
!> \details
!>  This routine allocates space for nonHaloData on the particlelist
!
!-----------------------------------------------------------------------
!  subroutine allocate_nonHalo_data(domain, particlelist) !{{{
!    implicit none

!    type (domain_type), intent(in) :: domain
!    type (mpas_particle_list_type), pointer, intent(in) :: particlelist
!    type (mpas_pool_type), pointer :: lagrPartTrackPool
!    type (mpas_pool_iterator_type) :: dimItr

!    integer :: i, nPart
!    real (kind=RKIND), dimension(:), pointer :: array1DRealPointer
!    integer, dimension(:), pointer :: array1DIntPointer
!    character (len=StrKIND) :: message


!    ! get number of particle on list
!    nPart = count_particlelist(particlelist)

!    ! allocate zero arrays
!    allocate(array1DRealPointer(nPart))
!    allocate(array1DIntPointer(nPart))
!    array1DRealPointer = 0.0_RKIND
!    array1DIntPointer = 0

!    ! on each list, transmit relevant fields to associated processors
!    call mpas_pool_get_subpool(domain % blocklist % structs, 'lagrPartTrackNonHalo', lagrPartTrackPool)
!    call mpas_pool_begin_iteration(lagrPartTrackPool)
!    do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
!      if (dimItr % memberType == MPAS_POOL_FIELD) then
!        if (dimItr % dataType == MPAS_POOL_REAL) then
!            call add_nonhalo_data_to_particle_list_array(particlelist, dimItr % memberName, array1DRealPointer)
!        elseif (dimItr % dataType == MPAS_POOL_INTEGER) then
!            call add_nonhalo_data_to_particle_list_array(particlelist, dimItr % memberName, array1DIntPointer)
!        else
!          LIGHT_DEBUG_ALL_WRITE("Different field type than implemented during halo communication!")
!        end if
!      elseif (dimItr % memberType == MPAS_POOL_DIMENSION) then
!        ! ignore dimensions for now and have this code so they aren't printed as an error message
!      else
!        write(message, *) "Different type expected in registry for key ", trim(dimItr % memberName), &
!                          " in halo data for communication-- don't know what to do!"
!        LIGHT_DEBUG_ALL_WRITE(message)
!      end if
!    end do

!    ! deallocate arrays
!    deallocate(array1DRealPointer)
!    deallocate(array1DIntPointer)

!  end subroutine allocate_nonHalo_data !}}}

!  subroutine allocate_list_nonHalo_data(domain, listPL) !{{{
!    implicit none

!    type (domain_type), intent(in) :: domain
!    type (mpas_list_of_particle_list_type), dimension(:), pointer, intent(inout) :: listPL

!    integer :: i, numList

!    numList = size(listPL)
!    do i=1, numList
!      call allocate_nonHalo_data(domain, listPL(i)%list)
!    end do

!  end subroutine allocate_list_nonHalo_data !}}}

!***********************************************************************
!
!  routine build_block_particlelists
!
!> \brief   Allocates empty particlelist
!> \author  Phillip Wolfram
!> \date    04/15/2014
!> \details
!>  This routine allocates empty particlelist data structures
!
!-----------------------------------------------------------------------
   subroutine build_block_particlelists(domain, err)!{{{

      implicit none

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_particle_list_type), pointer :: particlelist
      integer, pointer :: nParticles

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! allocate pointers
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_dimension(meshPool, 'nParticles', nParticles)

        ! allocate the memory in its location for use
        LIGHT_DEBUG_WRITE('in build_block_particlelists: nParticles = ' COMMA  nParticles)

        if (nParticles > 0) then
          allocate(block % particlelist)
          particlelist => block % particlelist

          !-----------------------------------------------------------------
          ! populate list of particles from input data structures
          !-----------------------------------------------------------------

          ! initialize the particlelist for population
          call build_new_particlelist(nParticles, particlelist, block % blockID)
        end if

        block => block % next
      end do

   end subroutine build_block_particlelists!}}}

   subroutine clear_block_particlelists(domain, err) !{{{
     implicit none

     !-----------------------------------------------------------------
     ! input/output variables
     !-----------------------------------------------------------------
     type (domain_type), intent(in) :: domain

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------
     integer, intent(out) :: err !< Output: error flag

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     type (block_type), pointer :: block

     err = 0

     block => domain % blocklist
     do while (associated(block))

       call mpas_particle_list_destroy_particle_list(block % particlelist)

       block => block % next
     end do

   end subroutine clear_block_particlelists !}}}

!***********************************************************************
!
!  routine read_haloData
!
!> \brief   Reads haloData from netCDF-injected struct arrays
!> \author  Phillip Wolfram
!> \date    04/15/2014
!> \details
!>  This routine reads haloData input for this MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------
   subroutine read_haloData(domain, err)!{{{

      implicit none

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: lagrPartTrackPool
      type (mpas_pool_iterator_type) :: dimItr
      type (mpas_particle_list_type), pointer :: particlelist
      type (field1DReal), pointer :: field1DRealPointer
      !type (field2DReal), pointer :: field2DRealPointer
      type (field1DInteger), pointer :: field1DIntPointer
      integer, dimension(:), pointer :: array1DInt
      character (len=StrKIND) :: message

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! allocate pointers
        particlelist => block % particlelist
        call mpas_pool_get_subpool(block % structs, 'lagrPartTrackHalo', lagrPartTrackPool)

        ! iterate over each member of the pool and make the relevant assignment
        call mpas_pool_begin_iteration(lagrPartTrackPool)
        do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
          ! determine the type of data
          if (dimItr % memberType == MPAS_POOL_FIELD) then
            if (dimItr % dataType == MPAS_POOL_REAL) then
              call mpas_pool_get_field(lagrPartTrackPool, dimItr % memberName, field1DRealPointer)
              LIGHT_DEBUG_WRITE(dimItr % memberName COMMA  ' = ')
              LIGHT_DEBUG_WRITE(field1DRealPointer % array)
              call add_halo_data_to_particle_list(particlelist, dimItr % memberName, field1DRealPointer)
            elseif (dimItr % dataType == MPAS_POOL_INTEGER) then
              call mpas_pool_get_field(lagrPartTrackPool, dimItr % memberName, field1DIntPointer)
              ! assign ioBlock explicitly during initial read
              if (dimItr % memberName == 'ioBlock') then
                field1DIntPointer % array = domain % dminfo % my_proc_id
                LIGHT_DEBUG_WRITE('ioBlock = ' COMMA  field1DIntPointer % array)
              end if
              LIGHT_DEBUG_WRITE(dimItr % memberName COMMA ' = ')
              LIGHT_DEBUG_WRITE(field1DIntPointer % array)
              call add_halo_data_to_particle_list(particlelist, dimItr % memberName, field1DIntPointer)
            else
              LIGHT_DEBUG_WRITE("Different field type than implemented during halo read!")
            end if
          elseif (dimItr % memberType == MPAS_POOL_DIMENSION) then
            ! ignore dimensions for now and have this code so they aren't printed as an error message
          else
            write(message, *) "Different type expected in registry for key ",  trim(dimItr % memberName), &
                              " in halo data for read-- don't know what to do!"
            LIGHT_DEBUG_WRITE(message)
            ! false warning for
            !Different type expected in registry for key on_a_sphere in Halo data for read, don't know what to do!
            !Different type expected in registry for key sphere_radius in Halo data for read, don't know what to do!
            !Different type expected in registry for key is_periodic in Halo data for read, don't know what to do!
            !Different type expected in registry for key x_period in Halo data for read, don't know what to do!
            !Different type expected in registry for key y_period in Halo data for read, don't know what to do!
          end if
        end do

         block => block % next
      end do
      LIGHT_DEBUG_WRITE('Finished reading halo data')

   end subroutine read_haloData!}}}

!***********************************************************************
!
!  routine read_nonhaloData
!
!> \brief   Reads nonhaloData from netCDF-injected struct arrays
!> \author  Phillip Wolfram
!> \date    04/15/2014
!> \details
!>  This routine reads nonhaloData input for this MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------
!  subroutine read_nonhaloData(domain, err)!{{{

!     implicit none

!     !-----------------------------------------------------------------
!     ! input/output variables
!     !-----------------------------------------------------------------
!     type (domain_type), intent(in) :: domain

!     !-----------------------------------------------------------------
!     ! output variables
!     !-----------------------------------------------------------------
!     integer, intent(out) :: err !< Output: error flag

!     !-----------------------------------------------------------------
!     ! local variables
!     !-----------------------------------------------------------------
!     type (block_type), pointer :: block
!     type (mpas_pool_type), pointer :: lagrPartTrackPool
!     type (mpas_pool_iterator_type) :: dimItr
!     type (mpas_particle_list_type), pointer :: particlelist
!     type (field1DReal), pointer :: field1DRealPointer
!     type (field1DInteger), pointer :: field1DIntPointer
!     character (len=StrKIND) :: message

!     err = 0

!     block => domain % blocklist
!     do while (associated(block))
!       ! allocate pointers
!       particlelist => block % particlelist
!       call mpas_pool_get_subpool(block % structs, 'lagrPartTrackNonHalo', lagrPartTrackPool)

!       ! iterate over each member of the pool and make the relevant assignment
!       call mpas_pool_begin_iteration(lagrPartTrackPool)
!       do while(mpas_pool_get_next_member(lagrPartTrackPool, dimItr))
!         ! determine the type of data
!         if (dimItr % memberType == MPAS_POOL_FIELD) then
!           if (dimItr % dataType == MPAS_POOL_REAL) then
!             call mpas_pool_get_field(lagrPartTrackPool, dimItr % memberName, field1DRealPointer)
!             call add_nonhalo_data_to_particle_list(particlelist, dimItr % memberName, field1DRealPointer)
!           elseif (dimItr % dataType == MPAS_POOL_INTEGER) then
!             call mpas_pool_get_field(lagrPartTrackPool, dimItr % memberName, field1DIntPointer)
!             call add_nonhalo_data_to_particle_list(particlelist, dimItr % memberName, field1DIntPointer)
!           else
!             LIGHT_DEBUG_WRITE("Different field type than implemented in nonHalo read!")
!           end if
!         elseif (dimItr % memberType == MPAS_POOL_DIMENSION) then
!           ! ignore dimensions for now and have this code so they aren't printed as an error message
!         else
!           write(message, *) "Different type expected in registry for key ", trim(dimItr % memberName),  &
!                             " in nonHalo data for read-- don't know what to do!"
!           LIGHT_DEBUG_WRITE(message)
!         end if
!       end do

!       ! alternatively, could initialize these fields or just make sure that they exist!

!        block => block % next
!     end do
!     LIGHT_DEBUG_WRITE('Finished reading non-halo data')

!  end subroutine read_nonhaloData!}}}

!***********************************************************************
!
!  routine mpas_particle_list_update_particle_block
!
!> \brief   Update particle block
!> \author  Phillip Wolfram
!> \date    10/28/2015
!> \details
!>  This routine updates the currentBlock for particles within the
!>  particlelist loop.
!
!-----------------------------------------------------------------------
   subroutine mpas_particle_list_update_particle_block(domain, block, particle, poolname, iCell) !{{{
     implicit none
     type (domain_type), intent(inout) :: domain
     type (block_type), intent(inout), pointer :: block
     type (mpas_particle_type), intent(inout), pointer :: particle
     character(len=*), intent(in) :: poolname
     integer, intent(inout) :: iCell

     ! local variables
     type (mpas_pool_type), pointer :: lagrPartTrackCellsPool
     integer, dimension(:), pointer :: cellOwnerBlock
     integer, pointer :: currentBlock, transfered

     call mpas_pool_get_subpool(block % structs, trim(poolname), lagrPartTrackCellsPool)
     call mpas_pool_get_array(lagrPartTrackCellsPool, 'cellOwnerBlock', cellOwnerBlock)
     call mpas_pool_get_array(particle % haloDataPool, 'currentBlock', currentBlock)
     if(cellOwnerBlock(iCell) /= currentBlock) then
       ! increment transfer counter
       call mpas_pool_get_array(particle % haloDataPool, 'transfered', transfered)
       transfered = transfered + 1
       ! set new current block
       currentBlock = cellOwnerBlock(iCell)
       ! reset cell_id to be brute force computed on new block after trasnfer
       iCell = -1
     end if

   end subroutine mpas_particle_list_update_particle_block !}}}

!***********************************************************************
!
!  routine mpas_particle_list_update_halos_start
!
!> \brief   Update halo
!> \author  Phillip Wolfram
!> \date    10/28/2015
!> \details
!>  This routine updates the halos for particles within the particlelist loop.
!>  This facilitates a computational transfer of a particle from one domain
!>  to another.
!
!-----------------------------------------------------------------------
   subroutine mpas_particle_list_update_halos_start(domain, block, particle, poolname, iCell, &
       arrayIndex, destinationName, sendProcRecvList, sendProcSendList, gsendProcNeighs ) !{{{
     implicit none
     type (domain_type), intent(inout) :: domain
     type (block_type), intent(inout), pointer :: block
     type (mpas_particle_type), intent(inout), pointer :: particle
     character(len=*), intent(in) :: poolname
     ! 'ioBlock' and 'currentBlock' are options for destinationName
     character(len=*), intent(in) :: destinationName
     integer, intent(inout) :: iCell, arrayIndex
     integer, dimension(:), pointer, intent(inout) :: gsendProcNeighs
     logical, dimension(:,:), pointer, intent(inout) :: sendProcRecvList
     logical, dimension(:), pointer, intent(inout) :: sendProcSendList

     ! local variables
     integer :: currentProc, sendProc
     type (mpas_pool_type), pointer :: lagrPartTrackCellsPool
     integer, pointer :: currentBlock, sendBlock

     call mpas_pool_get_subpool(block % structs, trim(poolname), lagrPartTrackCellsPool)
     call mpas_pool_get_array(particle % haloDataPool, 'currentBlock', currentBlock)
     call mpas_pool_get_array(particle % haloDataPool, destinationName, sendBlock)
     call mpas_get_owning_proc(domain % dminfo, currentBlock, currentProc)
     call mpas_get_owning_proc(domain % dminfo, sendBlock, sendProc)

     ! increment data for receiving processors (sum should be total number of particles on processor)
     LIGHT_DEBUG_WRITE('destinationName=' COMMA destinationName)
     LIGHT_DEBUG_WRITE('g_sendProcNeighs=' COMMA gsendProcNeighs)
     LIGHT_DEBUG_WRITE('sendProc=' COMMA sendProc)
     ! do not need to transfer information for particles on-processor, this is computed from the send list
     ! which is dependent upon current, on-processor particles
     if (sendProc /= domain % dminfo % my_proc_id) then
       arrayIndex = find_index(gsendProcNeighs, sendProc)
       sendProcRecvList(currentProc+1, arrayIndex) = .True.
     else
       ! consider the case where a particle on A has an sendProc of A and is sent to B (need to have B in A's halo).
       sendProcSendList(currentProc+1) = .True.
     end if
     ! must be computed after computational particles are transferred (this was a bug left-over from serial IO)
   end subroutine mpas_particle_list_update_halos_start!}}}

!}}}

!-----------------------------------------------------------------------
!
! TESTING SUBROUTINES
!
!-----------------------------------------------------------------------
!{{{
   subroutine mpas_particle_list_test_neighscalc(domain, err) !{{{
      implicit none

      type (domain_type), intent(in) :: domain
      integer, intent(out) :: err !< Output: error flag
      type (block_type), pointer :: block

      err = 0
      block => domain % blocklist
      do while (associated(block))

        ! write out all blockNeighs and procNeighs
        call mpas_log_write( 'blockID = $i', intArgs=(/  block % blockID /) )
        call mpas_log_write( 'blockNeighs = $i', intArgs=(/ block % blockNeighs /) )
        call mpas_log_write( 'procNeighs = $i', intArgs=(/ block % procNeighs /) )

        block => block % next
      end do

   end subroutine mpas_particle_list_test_neighscalc !}}}

   subroutine mpas_particle_list_test_numparticles_to_neighprocs(myproc, procNeighs, ioProcNeighs) !{{{
     implicit none
     integer, intent(in) :: myproc
     integer, dimension(:), pointer, intent(in) :: procNeighs, ioProcNeighs

     integer :: i, numNeighs

     numNeighs = size(procNeighs)

     LIGHT_DEBUG_WRITE( 'myproc, procNeighs')
     do i=1,numNeighs
       LIGHT_DEBUG_WRITE(myproc COMMA procNeighs(i))
     end do
     if(associated(ioProcNeighs)) then
       LIGHT_DEBUG_WRITE( 'myproc, ioProcNeighs')
       numNeighs = size(ioProcNeighs)
       do i=1,numNeighs
         LIGHT_DEBUG_WRITE(myproc COMMA ioProcNeighs(i))
       end do
     end if

   end subroutine mpas_particle_list_test_numparticles_to_neighprocs !}}}

   subroutine mpas_particle_list_test_num_current_particlelist(domain) !{{{
     implicit none
     type (domain_type), intent(in) :: domain
     type (block_type), pointer :: block
     integer :: nPartList, nPartList_particles

     block => domain % blocklist
     do while(associated(block))
       ! THIS LINE CAUSED A VERY, VERY, VERY NASTY BUG-- BE CAREFUL ABOUT GETTING THE LOCATION OF POTENTIALLY NULLS!
       nPartList = count_particlelist(block % particlelist)
       nPartList_particles = count_particlelist_particles(block % particlelist)
       LIGHT_DEBUG_WRITE('block = ' COMMA block % blockID COMMA ' nPartList= ' COMMA  nPartList COMMA ' nparticles = ' COMMA nPartList_particles)
       if (nPartList > nPartList_particles) then
         LIGHT_DEBUG_WRITE('Possible error! ' COMMA nPartList - nPartList_particles COMMA ' particles on particle list is not allocated!')
       end if

       block => block % next
     end do

   end subroutine mpas_particle_list_test_num_current_particlelist!}}}

   subroutine test_currentBlock(domain) !{{{
     implicit none
     type (domain_type), intent(in) :: domain
     type (block_type), pointer :: block
     type (mpas_pool_type), pointer :: lagrPartTrackPool
     integer, dimension(:), pointer :: field1DInt
     integer :: currentBlock, countOnProc, i
     countOnProc = 0
     block => domain % blocklist
     do while(associated(block))
       call mpas_pool_get_subpool(block % structs, 'lagrPartTrackHalo', lagrPartTrackPool)
       call mpas_pool_get_array(lagrPartTrackPool, 'currentBlock', field1DInt)
       do i=1,size(field1DInt(:))
         if (field1DInt(i) == block % blockID) countOnProc = countOnProc + 1
       end do

       block => block % next
     end do
     LIGHT_DEBUG_WRITE('number of particles on processor = ' COMMA countOnProc)

   end subroutine test_currentBlock !}}}

   subroutine test_num_particles_on_particlelist(listPL, nlist) !{{{
     implicit none
     integer, intent(in) :: nlist
     type (mpas_list_of_particle_list_type), dimension(:), pointer, intent(in) :: listPL

     integer :: i, sumtot, listparticles

     sumtot = 0
     do i = 1, nlist
       listparticles = count_particlelist(listPL(i)%list)
       LIGHT_DEBUG_WRITE('list i=' COMMA i COMMA ' nparticles=' COMMA listparticles)
       sumtot = sumtot + listparticles
     end do
     LIGHT_DEBUG_WRITE('total_on_list=' COMMA sumtot)
   end subroutine test_num_particles_on_particlelist !}}}
!}}}


end module ocn_particle_list
! vim: foldmethod=marker
