MODULE z_coord
! Prepares the (stretched) vertical coordinates and map factors

      USE shr_kind_mod, only: r8 => shr_kind_r8
      
      USE parmsld, only: nk3
      USE constld, only: zz,zt,fnz,fnt

IMPLICIT NONE
PRIVATE

PUBLIC :: coords_2d,coords_2d_mixed

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE coords_2d
! SUBROUTINE coords_2d_mixed
!-----------------------------------------------------------------------

!=======================================================================
   SUBROUTINE COORDS_2D ( CZ1, CZ2, DZ, ZB )
!======================================================================= 
!  OUTPUT: ZZ,ZT,FNZ,FNT
  
      REAL (KIND=r8), INTENT(IN) :: CZ1,CZ2,DZ,ZB
        
      INTEGER :: K

      ZZ(1) = ZB

      DO K = 2, NK3
       ZZ(K) = ZZ(K-1) + DZ
      ENDDO

      ZT(1) = ZZ(1)
      ZT(2) = ZZ(1) + DZ / 2.0_r8

      DO K = 3, NK3
       ZT(K) = ZT(K-1) + DZ
      ENDDO
      
!     DEFINE TRANSFORMATION FUNCTIONS ( KLEMP & CHEN ,1982 )
!     AND PHYSICAL COORDINATES

      DO K = 1, NK3
       FNZ(K) = 1.0_r8 / ( CZ1 + 2.0_r8 * CZ2 * ZZ(K) )
       ZZ(K) = ZZ(K) * ( CZ1 + CZ2 * ZZ(K) )
      ENDDO

      DO K = 1, NK3
       FNT(K) = 1.0_r8 / ( CZ1 + 2.0_r8 * CZ2 * ZT(K) )
       ZT(K) = ZT(K) * ( CZ1 + CZ2 * ZT(K) )
      ENDDO

   END SUBROUTINE COORDS_2D

!=======================================================================
   SUBROUTINE COORDS_2D_MIXED ( CZ1, CZ2, DZ, ZB, DZCONST )
!=======================================================================   
!  mixed vertical grids: constant below 4 km and stretched above
!  OUTPUT: ZZ,ZT,FNZ,FNT

      REAL (KIND=r8), INTENT(IN) :: CZ1,CZ2,DZ,ZB,DZCONST
         
      REAL (KIND=r8), PARAMETER :: ZFIX_TOP=4000.0_r8
      INTEGER :: LFIX_TOP,K  
      
      LFIX_TOP = INT(ZFIX_TOP/DZCONST)
   
!-------------------------------------      
!     Grids with a stretching depth
!------------------------------------- 
      ZZ(LFIX_TOP+1) = ZB
      
      DO K = LFIX_TOP+2,NK3
       ZZ(K) = ZZ(K-1) + DZ
      ENDDO

      ZT(LFIX_TOP+2) = ZZ(LFIX_TOP+1) + DZ/2.0_r8
      DO K = LFIX_TOP+3,NK3
       ZT(K) = ZT(K-1) + DZ
      ENDDO

!     DEFINE TRANSFORMATION FUNCTIONS ( KLEMP & CHEN ,1982 )
!     AND PHYSICAL COORDINATES

      DO K = LFIX_TOP+1, NK3
       FNZ(K) = 1.0_r8 / ( CZ1 + 2.0_r8 * CZ2 * ZZ(K) )
       ZZ(K) = ZZ(K) * ( CZ1 + CZ2 * ZZ(K) ) + ZFIX_TOP
      ENDDO

      DO K = LFIX_TOP+2, NK3
       FNT(K) = 1.0_r8 / ( CZ1 + 2.0_r8 * CZ2 * ZT(K) )
       ZT(K) = ZT(K) * ( CZ1 + CZ2 * ZT(K) ) + ZFIX_TOP
      ENDDO
      
!-------------------------------------      
!     Grids with a constant depth
!------------------------------------- 
      
      ZZ(1) = 0.0_r8
      DO K=2,LFIX_TOP
       ZZ(K) = ZZ(K-1) + DZCONST
      ENDDO
      DO K=1,LFIX_TOP
       FNZ(K) = DZ/DZCONST
      ENDDO      
      
      ZT(1) = ZZ(1)
      ZT(2) = ZZ(1) + DZCONST / 2.0_r8
      DO K = 3, LFIX_TOP+1
       ZT(K) = ZT(K-1) + DZCONST
      ENDDO
      DO K = 1, LFIX_TOP+1
       FNT(K) = DZ/DZCONST
      ENDDO

!-------------------------------------      
!     Correction of m factor 
!------------------------------------- 
      FNZ(LFIX_TOP+1)=DZ/(ZT(LFIX_TOP+2)-ZT(LFIX_TOP+1))
      
   END SUBROUTINE COORDS_2D_MIXED

END MODULE z_coord
