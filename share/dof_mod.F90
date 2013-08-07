#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dof_mod
  use kinds, only : real_kind,int_kind,long_kind
  use dimensions_mod, only : np, npsq, nelem, nelemd
  use quadrature_mod, only : quadrature_t
  use element_mod, only : element_t,index_t
  use parallel_mod, only : parallel_t, mpiinteger_t
  use edge_mod, only : longedgebuffer_t,initlongedgebuffer,freelongedgebuffer, &
		       longedgevpack, longedgevunpackmin
  use bndry_mod, only : bndry_exchangev
implicit none
private
  ! public data
  ! public subroutines
  public :: global_dof
  public :: genLocalDof
  public :: PrintDofP
  public :: UniquePoints
  public :: PutUniquePoints
  public :: UniqueNcolsP
  public :: UniqueCoords
  public :: CreateUniqueIndex
  public :: SetElemOffset
  public :: CreateMetaData

  interface UniquePoints
     module procedure UniquePoints2D
     module procedure UniquePoints3D
     module procedure UniquePoints4D
  end interface
  interface PutUniquePoints
     module procedure PutUniquePoints2D
     module procedure PutUniquePoints3D
     module procedure PutUniquePoints4D
  end interface


contains

  subroutine genLocalDof(ig,npts,ldof)

    integer(kind=int_kind), intent(in) :: ig
    integer(kind=int_kind), intent(in) :: npts
    integer(kind=int_kind), intent(inout) :: ldof(:,:)

    integer(kind=int_kind) :: i,j,npts2
  
   
     npts2=npts*npts
     do j=1,npts
       do i=1,npts
           ldof(i,j) = (ig-1)*npts2 + (j-1)*npts + i
       enddo
     enddo

  end subroutine genLocalDOF

! ===========================================
! global_dof
!
! Compute the global degree of freedom for each element...
! ===========================================

  subroutine global_dof(par,elem)

    type (parallel_t),intent(in) :: par
    type (element_t)             :: elem(:)

    type (LongEdgeBuffer_t)    :: edge

    real(kind=real_kind)  da                     ! area element

    type (quadrature_t) :: gp

    integer (kind=int_kind) :: ldofP(np,np,nelemd)

    integer ii
    integer i,j,ig,ie
    integer kptr
    integer iptr

    ! ===================
    ! begin code
    ! ===================
    call initLongEdgeBuffer(edge,1)

    ! =================================================
    ! mass matrix on the velocity grid
    ! =================================================    

 
    do ie=1,nelemd
       ig = elem(ie)%vertex%number
       call genLocalDOF(ig,np,ldofP(:,:,ie))
	 
       kptr=0
       call LongEdgeVpack(edge,ldofP(:,:,ie),1,kptr,elem(ie)%desc)
    end do

    ! ==============================
    ! Insert boundary exchange here
    ! ==============================

    call bndry_exchangeV(par,edge)

    do ie=1,nelemd
       ! we should unpack directly into elem(ie)%gdofV, but we dont have
       ! a VunpackMIN that takes integer*8.  gdofV integer*8 means  
       ! more than 2G grid points.
       kptr=0
       call LongEdgeVunpackMIN(edge,ldofP(:,:,ie),1,kptr,elem(ie)%desc)
       elem(ie)%gdofP(:,:)=ldofP(:,:,ie)
    end do
#if (! defined VERT_OPENMP)
!$OMP BARRIER
#endif
    call FreeLongEdgeBuffer(edge)
       
  end subroutine global_dof


  subroutine UniquePoints2D(idxUnique,src,dest)
    type (index_t) :: idxUnique
    real (kind=real_kind) :: src(:,:)
    real (kind=real_kind) :: dest(:)

    integer(kind=int_kind) :: i,j,ii
    

    do ii=1,idxUnique%NumUniquePts
       i=idxUnique%ia(ii)
       j=idxUnique%ja(ii)
       dest(ii)=src(i,j)
    enddo

  end subroutine UniquePoints2D

! putUniquePoints first zeros out the destination array, then fills the unique points of the 
! array with values from src.  A boundary communication should then be called to fill in the 
! redundent points of the array

  subroutine putUniquePoints2D(idxUnique,src,dest)
    type (index_t) :: idxUnique
    real (kind=real_kind),intent(in) :: src(:)
    real (kind=real_kind),intent(out) :: dest(:,:)

    integer(kind=int_kind) :: i,j,ii
    
    dest=0.0D0
    do ii=1,idxUnique%NumUniquePts
       i=idxUnique%ia(ii)
       j=idxUnique%ja(ii)
       dest(i,j)=src(ii)
    enddo

  end subroutine putUniquePoints2D

  subroutine UniqueNcolsP(elem,idxUnique,cid)    
    use element_mod, only : GetColumnIdP, element_t
    type (element_t), intent(in) :: elem
    type (index_t), intent(in) :: idxUnique
    integer,intent(out) :: cid(:)
    integer(kind=int_kind) :: i,j,ii


    do ii=1,idxUnique%NumUniquePts
       i=idxUnique%ia(ii)
       j=idxUnique%ja(ii)
       cid(ii)=GetColumnIdP(elem,i,j)
    enddo
    
  end subroutine UniqueNcolsP


  subroutine UniqueCoords(idxUnique,src,lat,lon)

    use coordinate_systems_mod, only  : spherical_polar_t
    type (index_t), intent(in) :: idxUnique

    type (spherical_polar_t) :: src(:,:)
    real (kind=real_kind), intent(out) :: lat(:)
    real (kind=real_kind), intent(out) :: lon(:)

    integer(kind=int_kind) :: i,j,ii

    do ii=1,idxUnique%NumUniquePts
       i=idxUnique%ia(ii)
       j=idxUnique%ja(ii)
       lat(ii)=src(i,j)%lat
       lon(ii)=src(i,j)%lon
    enddo

  end subroutine UniqueCoords

  subroutine UniquePoints3D(idxUnique,nlyr,src,dest)
    type (index_t) :: idxUnique
    integer(kind=int_kind) :: nlyr
    real (kind=real_kind) :: src(:,:,:)
    real (kind=real_kind) :: dest(:,:)
    
    integer(kind=int_kind) :: i,j,k,ii

    do ii=1,idxUnique%NumUniquePts
       i=idxUnique%ia(ii)
       j=idxUnique%ja(ii)
       do k=1,nlyr
          dest(ii,k)=src(i,j,k)
       enddo
    enddo

  end subroutine UniquePoints3D
  subroutine UniquePoints4D(idxUnique,d3,d4,src,dest)
    type (index_t) :: idxUnique
    integer(kind=int_kind) :: d3,d4
    real (kind=real_kind) :: src(:,:,:,:)
    real (kind=real_kind) :: dest(:,:,:)
    
    integer(kind=int_kind) :: i,j,k,n,ii

    do n=1,d4
       do k=1,d3
          do ii=1,idxUnique%NumUniquePts
             i=idxUnique%ia(ii)
             j=idxUnique%ja(ii)
             dest(ii,k,n)=src(i,j,k,n)
          enddo
       end do
    enddo

  end subroutine UniquePoints4D

! putUniquePoints first zeros out the destination array, then fills the unique points of the 
! array with values from src.  A boundary communication should then be called to fill in the 
! redundent points of the array

  subroutine putUniquePoints3D(idxUnique,nlyr,src,dest)
    type (index_t) :: idxUnique
    integer(kind=int_kind) :: nlyr
    real (kind=real_kind),intent(in) :: src(:,:)
    real (kind=real_kind),intent(out) :: dest(:,:,:)
    
    integer(kind=int_kind) :: i,j,k,ii

    dest=0.0D0
    do k=1,nlyr
       do ii=1,idxUnique%NumUniquePts
          i=idxUnique%ia(ii)
          j=idxUnique%ja(ii)
          dest(i,j,k)=src(ii,k)
       enddo
    enddo

  end subroutine putUniquePoints3D

  subroutine putUniquePoints4D(idxUnique,d3,d4,src,dest)
    type (index_t) :: idxUnique
    integer(kind=int_kind) :: d3,d4
    real (kind=real_kind),intent(in) :: src(:,:,:)
    real (kind=real_kind),intent(out) :: dest(:,:,:,:)
    
    integer(kind=int_kind) :: i,j,k,n,ii

    dest=0.0D0
    do n=1,d4
       do k=1,d3
          do ii=1,idxunique%NumUniquePts
             i=idxUnique%ia(ii)
             j=idxUnique%ja(ii)
             dest(i,j,k,n)=src(ii,k,n)
          enddo
       enddo
    end do
  end subroutine putUniquePoints4D

  subroutine SetElemOffset(par,elem,GlobalUniqueColsP)
#ifdef _MPI
     use parallel_mod, only : mpi_sum
#endif
     type (parallel_t) :: par
     type (element_t) :: elem(:)
     integer, intent(out) :: GlobalUniqueColsP

     integer(kind=int_kind), allocatable :: numElemP(:),numElem2P(:)
     integer(kind=int_kind), allocatable :: numElemV(:),numElem2V(:)
     integer(kind=int_kind), allocatable :: gOffset(:)
    
     integer(kind=int_kind) :: ie,ig,nprocs,ierr

     logical,parameter :: Debug = .FALSE.

     nprocs = par%nprocs
     allocate(numElemP(nelem))
     allocate(numElem2P(nelem))
     allocate(gOffset(nelem))
     numElemP=0;numElem2P=0;gOffset=0

     do ie=1,nelemd
	ig = elem(ie)%GlobalId
	numElemP(ig) = elem(ie)%idxP%NumUniquePts
     enddo
#ifdef _MPI
     call MPI_Allreduce(numElemP,numElem2P,nelem,MPIinteger_t,MPI_SUM,par%comm,ierr) 
#else
     numElem2P=numElemP
#endif

     gOffset(1)=1
     do ig=2,nelem
	gOffset(ig) = gOffset(ig-1)+numElem2P(ig-1)
     enddo
     do ie=1,nelemd
        ig = elem(ie)%GlobalId
        elem(ie)%idxP%UniquePtOffset=gOffset(ig)
     enddo
     GlobalUniqueColsP = gOffset(nelem)+numElem2P(nelem)-1

     deallocate(numElemP)
     deallocate(numElem2P)
     deallocate(gOffset)
  end subroutine SetElemOffset

  subroutine CreateUniqueIndex(ig,gdof,idx)

    integer(kind=int_kind) :: ig
    type (index_t) :: idx 
    integer(kind=long_kind) :: gdof(:,:)
    
    integer, allocatable :: ldof(:,:)
    integer :: i,j,ii,npts


    npts = size(gdof,dim=1)
    allocate(ldof(npts,npts))
    ! ====================
    ! Form the local DOF
    ! ====================
    call genLocalDOF(ig,npts,ldof)
    
    ii=1
    
    do j=1,npts
       do i=1,npts
          ! ==========================
          ! check for point ownership
          ! ==========================
          if(gdof(i,j) .eq. ldof(i,j)) then
             idx%ia(ii) = i
             idx%ja(ii) = j
             ii=ii+1
          endif
       enddo
    enddo
    
    idx%NumUniquePts=ii-1
    deallocate(ldof)

  end subroutine CreateUniqueIndex


  subroutine CreateMetaData(par,elem,subelement_corners, fdofp)
    type (parallel_t),intent(in) :: par
    type (element_t), target    :: elem(:)

    integer, intent(out),optional         :: subelement_corners((np-1)*(np-1)*nelemd,4)
    integer(kind=int_kind), optional :: fdofp(np,np,nelemd)

    type (index_t), pointer  :: idx 
    type (LongEdgeBuffer_t)    :: edge
    integer :: i, j, ii, ie, base
    integer(kind=long_kind), pointer :: gdof(:,:)
    integer :: fdofp_local(np,np,nelemd)

    call initLongEdgeBuffer(edge,1)
    fdofp_local=0
    
    do ie=1,nelemd
       idx => elem(ie)%idxP
       do ii=1,idx%NumUniquePts
          i=idx%ia(ii)
          j=idx%ja(ii)
          
          fdofp_local(i,j,ie) = -(idx%UniquePtoffset+ii-1)
       end do
       call LongEdgeVpack(edge,fdofp_local(:,:,ie),1,0,elem(ie)%desc)
    end do
    call bndry_exchangeV(par,edge)
    do ie=1,nelemd
       base = (ie-1)*(np-1)*(np-1)
       call LongEdgeVunpackMIN(edge,fdofp_local(:,:,ie),1,0,elem(ie)%desc)
       if(present(subelement_corners)) then
          ii=0       
          do j=1,np-1
             do i=1,np-1
                ii=ii+1
                subelement_corners(base+ii,1) = -fdofp_local(i,j,ie)
                subelement_corners(base+ii,2) = -fdofp_local(i,j+1,ie)
                subelement_corners(base+ii,3) = -fdofp_local(i+1,j+1,ie)
                subelement_corners(base+ii,4) = -fdofp_local(i+1,j,ie)
             end do
          end do
       end if
    end do
    if(present(fdofp)) then
       fdofp=-fdofp_local
    end if
    


  end subroutine CreateMetaData


! ==========================================
!  PrintDofP
!
!   Prints the degree of freedom 
! ==========================================
  subroutine PrintDofP(elem)

   implicit none
   type (element_t), intent(in) :: elem(:)

   integer :: ie,nse,i,j
   

   nse = SIZE(elem)
 
   do ie=1,nse
      print *,'Element # ',elem(ie)%vertex%number
      do j=np,1,-1
         write(6,*) (elem(ie)%gdofP(i,j), i=1,np)
      enddo
   enddo
 10 format('I5')

 end subroutine PrintDofP

end module dof_mod

