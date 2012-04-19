#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!=======================================================================================================!
program dg3d_main
!=======================================================================================================!
        use init_mod, only : init
! -----------------------------------------------
        use dg3d_primeq_mod, only : primeq_dg
! -----------------------------------------------
!      use kinds
! -----------------------------------------------
      use parallel_mod, only : parallel_t, initmp, syncmp, haltmp
! -----------------------------------------------
      use thread_mod, only : nthreads, omp_get_thread_num, omp_set_num_threads
! -----------------------------------------------
!      use time_mod
! -----------------------------------------------
      use dimensions_mod, only : nelemd, np, ne, nlev
! -----------------------------------------------
! -----------------------------------------------
      use domain_mod, only : domain1d_t, decompose
! -----------------------------------------------
      use element_mod, only : element_t
! -----------------------------------------------
!      use state_mod
! -----------------------------------------------
      use edge_mod, only : EdgeBuffer_t
! -----------------------------------------------
      use reduction_mod, only : ReductionBuffer_ordered_1d_t  
! -----------------------------------------------
      use hybvcoord_mod, only : hvcoord_t, hvcoord_init
! -----------------------------------------------
      use control_mod, only : vfile_mid, vfile_int
! -----------------------------------------------
      use hybrid_mod, only : hybrid_t
      
!=======================================================================================================!
      implicit none
      type (Element_t), pointer :: elem(:)
      type (EdgeBuffer_t)  :: edge1            ! 1 component edge buffer (1, 3d scalar field)
      type (EdgeBuffer_t)  :: edge2            ! 2 component edge buffer (1, 3d vector field)
      type (EdgeBuffer_t)  :: edge3            ! 3 component edge buffer (1, 3d vector + 1 3d scalar field)
      type (ReductionBuffer_ordered_1d_t)  :: red    ! reduction buffer for cg
      type (parallel_t)    :: par              ! parallel structure for distributed memory programming
      type (domain1d_t), allocatable :: dom_mt(:)
      type (hybrid_t)       :: hybrid ! parallel structure for shared memory/distributed memory
      type (hvcoord_t)     :: hvcoord         ! hybrid vertical coordinate struct

      integer nets,nete
      integer ithr
      integer ierr
!=======================================================================================================!

!=======================================================================================================!
!	Begin executable code set distributed memory world...						!
!=======================================================================================================!      
      par=initmp()
      call init(elem,edge1,edge2,edge3,red,par)
  
!=======================================================================================================!
!	Allocate state variables									!
!=======================================================================================================!
      if(par%masterproc) print *,'allocating state variables...'
 
!=======================================================================================================!
!	Set number of threads...									!
!=======================================================================================================!
      if(par%masterproc) print *,'Main:NThreads=',NThreads

      call omp_set_num_threads(NThreads)
     
      allocate(dom_mt(0:NThreads-1))
      do ithr=0,NThreads-1
         dom_mt(ithr)=decompose(1,nelemd,NThreads,ithr)
      end do
!=======================================================================================================!
!	Sync-up to make sure timing is clean								!
!=======================================================================================================!
      call syncmp(par)

!=======================================================================================================!
!	Begin threaded region...									!
!=======================================================================================================!
      ithr=omp_get_thread_num()
      nets=dom_mt(ithr)%start
      nete=dom_mt(ithr)%end
!=======================================================================================================!
!	Initialize thread decomposition									!
!=======================================================================================================!
      write(6,9) par%rank,ithr,nets,nete 
  9   format('process: ',i5,1x,'thread: ',i2,1x,'element limits: ',i4," - ",i4)

! ==================================
! Initialize the vertical coordinate  (cam initializes hvcoord externally)
! ==================================

      hvcoord = hvcoord_init(vfile_mid, vfile_int, .true., hybrid%par%masterproc, ierr)
      if (ierr /= 0) then
         call haltmp("error in hvcoord_init")
      end if
      
      call primeq_dg(elem,edge1,edge2,edge3,red,par,ithr,nets,nete, hvcoord)

      if(par%masterproc) then 
         print *,'Finishing DG 3D Computation...'
      endif
!=======================================================================================================!
!	Sync-up to make sure timing is clean								!
!=======================================================================================================!
      call syncmp(par)

!=======================================================================================================!

10    format(1x,'#process:',i5,', nlev:',i3,', np:',i2,', nv:',i2,', ne:',i2,', nelem:',i5,&
           ', model total time (sec)=',f10.2)

!=======================================================================================================!
!	End distributed memory region									!
!=======================================================================================================!
      call haltmp('exiting program...')
      
!=======================================================================================================!
end program dg3d_main
