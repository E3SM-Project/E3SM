!------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS
!------------------------------------------------------------------------
      MODULE unstructured
!BOP
!
! !MODULE: unstructured
!
! !USES:
#include "pilgrim.h"
      IMPLICIT NONE

!
! !DESCRIPTION:



! !REVISION HISTORY:
!   01.10.30   Sawyer     Creation
!   02.08.14   Sawyer     Added explicit precisions from pilgrim.h
!
! !PUBLIC TYPES:
      PUBLIC Pump, GetPoint, AddTriangle

!EOP
      CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Pump --- Increase numbers of points and triangles 
!
! !INTERFACE:
      subroutine Pump( X, Y, Z, Tri )
! !USES:
      implicit none
 
! !INPUT/OUTPUT PARAMETERS:
      real(CPP_REAL8), pointer     :: X(:)          ! X
      real(CPP_REAL8), pointer     :: Y(:)          ! Y
      real(CPP_REAL8), pointer     :: Z(:)          ! Z
      integer, pointer      :: Tri(:,:)      ! Corner vertices

! !DESCRIPTION:
!
!     This routine puts a 2-D ghost region at the end of a buffer, first in
!     X then in Y.
!
! !LOCAL VARIABLES:
      integer  :: NVert, NTri, NVertNew, NTriNew
      real(CPP_REAL8), pointer :: Xnew(:), Ynew(:), Znew(:)
      integer, pointer  :: TriNew(:,:)
      integer, allocatable :: DB(:,:)

!
! !REVISION HISTORY:
!   01.10.30   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
       integer avert, bvert, cvert, abvert, bcvert, acvert
       integer i, j, k, count


       NVert    = SIZE(X)
       NTri     = SIZE(Tri, 2)

       NVertNew = 3*NTri/2 + NVert
       NTriNew  = 4*NTri

       ALLOCATE( DB(NVert,12) )   ! The vertex database 
       DB(:,:) = 0

       ALLOCATE( Xnew(NVertNew) )
       ALLOCATE( Ynew(NVertNew) )
       ALLOCATE( Znew(NVertNew) )
       ALLOCATE( TriNew(3,NTriNew) )

!
! Copy over old coordinates (adopt the vertex numbering)
!
       Xnew(1:NVert) = X(1:NVert)
       Ynew(1:NVert) = Y(1:NVert)
       Znew(1:NVert) = Z(1:NVert)

       Count = 0
!
! For each triangle
       DO i=1,NTri
         avert=min(min(Tri(1,i),Tri(2,i)),Tri(3,i))
         cvert=max(max(Tri(1,i),Tri(2,i)),Tri(3,i))
         DO k=1,3
           IF (Tri(k,i) /= avert .AND. Tri(k,i) /= cvert) bvert=Tri(k,i)
         ENDDO
         CALL GetPoint(avert,bvert,X,Y,Z,DB,NVert,Xnew,Ynew,Znew,abvert)
         CALL GetPoint(bvert,cvert,X,Y,Z,DB,NVert,Xnew,Ynew,Znew,bcvert)
         CALL GetPoint(avert,cvert,X,Y,Z,DB,NVert,Xnew,Ynew,Znew,acvert)


!  Add new triangles
         
         CALL AddTriangle(avert,abvert,acvert,Count,TriNew)
         CALL AddTriangle(bvert,abvert,bcvert,Count,TriNew)
         CALL AddTriangle(cvert,bcvert,acvert,Count,TriNew)
         CALL AddTriangle(abvert,bcvert,acvert,Count,TriNew)
       ENDDO

       IF (NVert /= NVertNew) print *, "Error in Pump", NVert,"!=",NVertNew
       IF (Count /= NTriNew) print *, "Error in Pump", Count,"!=",NTriNew

!
! Swap the old and new data sets
!
       DEALLOCATE(X)
       DEALLOCATE(Y)
       DEALLOCATE(Z)
       DEALLOCATE(Tri)
       X => Xnew
       Y => Ynew
       Z => Znew
       Tri => TriNew

       DEALLOCATE(DB)
!EOC
      END SUBROUTINE Pump
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: GetPoint --- Get the index of the new point
!
! !INTERFACE:
      SUBROUTINE GetPoint(avert,bvert,X,Y,Z,DB,NVert,Xnew,Ynew,Znew,abvert)

! !USES:
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in)  :: avert
      integer, intent(in)  :: bvert
      real(CPP_REAL8), pointer    :: X(:)
      real(CPP_REAL8), pointer    :: Y(:)
      real(CPP_REAL8), pointer    :: Z(:)

! !INPUT/OUTPUT PARAMETERS:
      integer              :: DB(:,:)
      real(CPP_REAL8), pointer    :: Xnew(:)
      real(CPP_REAL8), pointer    :: Ynew(:)
      real(CPP_REAL8), pointer    :: Znew(:)
      integer, intent(inout)             :: NVert

! !OUTPUT PARAMETERS:
      integer, intent(out)               :: abvert

! !DESCRIPTION:
!
!     Determine the midpoint if it has not already been done
!
! !REVISION HISTORY:
!   01.10.30   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC

      INTEGER k
      LOGICAL Found
      REAL (CPP_REAL8) :: xab, yab, zab, norm

         k = 0
         Found = .FALSE.
         DO WHILE ( .NOT. Found  .AND. k < 6 )
           k = k+1
           IF ( DB(avert,k) == bvert ) THEN
             abvert = DB(avert,k+6)
             Found = .TRUE.
           ELSEIF ( DB(avert,k) == 0 ) THEN
             NVert = NVert+1
             abvert = NVert
             DB(avert,k) = bvert
             DB(avert,k+6) = NVert
             xab = x(avert)+x(bvert)
             yab = y(avert)+y(bvert)
             zab = z(avert)+z(bvert)
             norm = sqrt(xab*xab + yab*yab + zab*zab)
             xnew(abvert) = xab / norm
             ynew(abvert) = yab / norm
             Znew(abvert) = zab / norm
             Found = .TRUE.
           ENDIF
         ENDDO

!EOC
      END SUBROUTINE GetPoint
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: AddTriangle --- Add a triangle to the list
      SUBROUTINE AddTriangle(vert1,vert2,vert3,NTri,TriNew)

! !USES:
      implicit none
        integer, intent(in)    :: vert1, vert2, vert3
        integer, intent(inout) :: NTri
        integer, pointer       :: TriNew(:,:)

        NTri = NTri+1
        TriNew(1,NTri) = vert1
        TriNew(2,NTri) = vert2
        TriNew(3,NTri) = vert3

!EOC
      END SUBROUTINE AddTriangle
!-------------------------------------------------------------------------

      END MODULE unstructured

