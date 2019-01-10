MODULE bound_channel_module

! Fill the halo points (along the channel) at the cube surface edges.
! Since the whole channel is in the same node, there is no mpi-communication.
! Mapping of halo points is done by separate procedures (halo_q, halo_vort, and halo_z).

! All variables are handled as if they were at q-point, except for z3dz & psi.
! For z3dz & psi, the possibility of overlap at the Edge is considered (data shift).
! (In the stand-alone version, bound_sync & data shift are done in halo_z.)

USE shr_kind_mod,   only: r8 => shr_kind_r8
USE vvm_data_types, only: channel_t

USE parmsld,        only: ntracer,nk1,nk2

IMPLICIT NONE
PRIVATE

! public member functions
PUBLIC ::  bound_channel

CONTAINS

! Local Subroutines:
!-----------------------------------------------------------------------
! SUBROUTINE bound_channel
! SUBROUTINE halo_group1 (& halo_group1_2)
! SUBROUTINE halo_group2 (& halo_group2_2)
! SUBROUTINE halo_group3 (& halo_group3_2)
!-----------------------------------------------------------------------

!=========================================================================================
      SUBROUTINE BOUND_CHANNEL (num_halo,channel,TH3D,QV3D,QT3D,  &
                                QC3D,QI3D,QR3D,QS3D,QG3D,         &
                                W3D,Z3DX,Z3DY,Z3DZ,               &
                                PSI,CHI,TOPOZ,ZROUGH,GWET,TG)
!=========================================================================================
!     Number of halo points specified can be different (NUM_HALO),
!     but, the horizontal array size of the variables is fixed, i.e., (mim:mip,mjm:mjp)
!-----------------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: num_halo
      type(channel_t), INTENT(INOUT) :: channel   ! channel data

      LOGICAL, INTENT(IN), OPTIONAL :: TH3D,QV3D,QT3D,QC3D,QI3D,QR3D,QS3D,QG3D, &
                                       W3D,Z3DX,Z3DY,Z3DZ, &
                                       PSI,CHI,TOPOZ,ZROUGH,GWET,TG

      ! Local
      INTEGER :: nt,mi1(4),mim(4),mip(4),mj1(4),mjm(4),mjp(4)

      ! Size of segment 1
      mi1(1) = channel%seg(1)%mi1
      mim(1) = channel%seg(1)%mim
      mip(1) = channel%seg(1)%mip
      mj1(1) = channel%seg(1)%mj1
      mjm(1) = channel%seg(1)%mjm
      mjp(1) = channel%seg(1)%mjp

!---------------------------------
! TH3D
!---------------------------------
      IF (PRESENT(TH3D)) THEN
        if (TH3D) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%TH3D(:,:,2:NK2),            &
                            channel%seg(2)%TH3D(:,:,2:NK2),            &
                            channel%seg(3)%TH3D(:,:,2:NK2),            &
                            channel%seg(4)%TH3D(:,:,2:NK2))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%TH3D(:,:,2:NK2),            &
                            channel%seg(2)%TH3D(:,:,2:NK2),            &
                            channel%seg(3)%TH3D(:,:,2:NK2),            &
                            channel%seg(4)%TH3D(:,:,2:NK2))

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%TH3D(:,:,2:NK2),            &
                            channel%seg(2)%TH3D(:,:,2:NK2),            &
                            channel%seg(3)%TH3D(:,:,2:NK2),            &
                            channel%seg(4)%TH3D(:,:,2:NK2))
        END SELECT
        endif
      ENDIF
!---------------------------------
! QV3D
!---------------------------------
      IF (PRESENT(QV3D)) THEN
        if (QV3D) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QV3D(:,:,2:NK2),            &
                            channel%seg(2)%QV3D(:,:,2:NK2),            &
                            channel%seg(3)%QV3D(:,:,2:NK2),            &
                            channel%seg(4)%QV3D(:,:,2:NK2))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QV3D(:,:,2:NK2),            &
                            channel%seg(2)%QV3D(:,:,2:NK2),            &
                            channel%seg(3)%QV3D(:,:,2:NK2),            &
                            channel%seg(4)%QV3D(:,:,2:NK2))

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QV3D(:,:,2:NK2),            &
                            channel%seg(2)%QV3D(:,:,2:NK2),            &
                            channel%seg(3)%QV3D(:,:,2:NK2),            &
                            channel%seg(4)%QV3D(:,:,2:NK2))
        END SELECT
        endif
      ENDIF
!---------------------------------
! QC3D
!---------------------------------
      IF (PRESENT(QC3D)) THEN
        if (QC3D) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QC3D(:,:,2:NK2),            &
                            channel%seg(2)%QC3D(:,:,2:NK2),            &
                            channel%seg(3)%QC3D(:,:,2:NK2),            &
                            channel%seg(4)%QC3D(:,:,2:NK2))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QC3D(:,:,2:NK2),            &
                            channel%seg(2)%QC3D(:,:,2:NK2),            &
                            channel%seg(3)%QC3D(:,:,2:NK2),            &
                            channel%seg(4)%QC3D(:,:,2:NK2))

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QC3D(:,:,2:NK2),            &
                            channel%seg(2)%QC3D(:,:,2:NK2),            &
                            channel%seg(3)%QC3D(:,:,2:NK2),            &
                            channel%seg(4)%QC3D(:,:,2:NK2))
        END SELECT
        endif
      ENDIF
!---------------------------------
! QI3D
!---------------------------------
      IF (PRESENT(QI3D)) THEN
        if (QI3D) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QI3D(:,:,2:NK2),            &
                            channel%seg(2)%QI3D(:,:,2:NK2),            &
                            channel%seg(3)%QI3D(:,:,2:NK2),            &
                            channel%seg(4)%QI3D(:,:,2:NK2))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QI3D(:,:,2:NK2),            &
                            channel%seg(2)%QI3D(:,:,2:NK2),            &
                            channel%seg(3)%QI3D(:,:,2:NK2),            &
                            channel%seg(4)%QI3D(:,:,2:NK2))

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QI3D(:,:,2:NK2),            &
                            channel%seg(2)%QI3D(:,:,2:NK2),            &
                            channel%seg(3)%QI3D(:,:,2:NK2),            &
                            channel%seg(4)%QI3D(:,:,2:NK2))
        END SELECT
        endif
      ENDIF
!---------------------------------
! QR3D
!---------------------------------
      IF (PRESENT(QR3D)) THEN
        if (QR3D) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QR3D(:,:,2:NK2),            &
                            channel%seg(2)%QR3D(:,:,2:NK2),            &
                            channel%seg(3)%QR3D(:,:,2:NK2),            &
                            channel%seg(4)%QR3D(:,:,2:NK2))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QR3D(:,:,2:NK2),            &
                            channel%seg(2)%QR3D(:,:,2:NK2),            &
                            channel%seg(3)%QR3D(:,:,2:NK2),            &
                            channel%seg(4)%QR3D(:,:,2:NK2))

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QR3D(:,:,2:NK2),            &
                            channel%seg(2)%QR3D(:,:,2:NK2),            &
                            channel%seg(3)%QR3D(:,:,2:NK2),            &
                            channel%seg(4)%QR3D(:,:,2:NK2))
        END SELECT
        endif
      ENDIF
!---------------------------------
! QS3D
!---------------------------------
      IF (PRESENT(QS3D)) THEN
        if (QS3D) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QS3D(:,:,2:NK2),            &
                            channel%seg(2)%QS3D(:,:,2:NK2),            &
                            channel%seg(3)%QS3D(:,:,2:NK2),            &
                            channel%seg(4)%QS3D(:,:,2:NK2))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QS3D(:,:,2:NK2),            &
                            channel%seg(2)%QS3D(:,:,2:NK2),            &
                            channel%seg(3)%QS3D(:,:,2:NK2),            &
                            channel%seg(4)%QS3D(:,:,2:NK2))

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QS3D(:,:,2:NK2),            &
                            channel%seg(2)%QS3D(:,:,2:NK2),            &
                            channel%seg(3)%QS3D(:,:,2:NK2),            &
                            channel%seg(4)%QS3D(:,:,2:NK2))
        END SELECT
        endif
      ENDIF
!---------------------------------
! QG3D
!---------------------------------
      IF (PRESENT(QG3D)) THEN
        if (QG3D) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QG3D(:,:,2:NK2),            &
                            channel%seg(2)%QG3D(:,:,2:NK2),            &
                            channel%seg(3)%QG3D(:,:,2:NK2),            &
                            channel%seg(4)%QG3D(:,:,2:NK2))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QG3D(:,:,2:NK2),            &
                            channel%seg(2)%QG3D(:,:,2:NK2),            &
                            channel%seg(3)%QG3D(:,:,2:NK2),            &
                            channel%seg(4)%QG3D(:,:,2:NK2))

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QG3D(:,:,2:NK2),            &
                            channel%seg(2)%QG3D(:,:,2:NK2),            &
                            channel%seg(3)%QG3D(:,:,2:NK2),            &
                            channel%seg(4)%QG3D(:,:,2:NK2))
        END SELECT
        endif
      ENDIF
!---------------------------------
! QT3D
!---------------------------------
      IF (PRESENT(QT3D)) THEN
        if (QT3D) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          DO nt = 1,ntracer
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QT3D(:,:,2:NK2,nt),         &
                            channel%seg(2)%QT3D(:,:,2:NK2,nt),         &
                            channel%seg(3)%QT3D(:,:,2:NK2,nt),         &
                            channel%seg(4)%QT3D(:,:,2:NK2,nt))
          ENDDO
        CASE(2)
          DO nt = 1,ntracer
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QT3D(:,:,2:NK2,nt),         &
                            channel%seg(2)%QT3D(:,:,2:NK2,nt),         &
                            channel%seg(3)%QT3D(:,:,2:NK2,nt),         &
                            channel%seg(4)%QT3D(:,:,2:NK2,nt))
          ENDDO
        CASE(3)
          DO nt = 1,ntracer
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1,                              &
                            channel%seg(1)%QT3D(:,:,2:NK2,nt),         &
                            channel%seg(2)%QT3D(:,:,2:NK2,nt),         &
                            channel%seg(3)%QT3D(:,:,2:NK2,nt),         &
                            channel%seg(4)%QT3D(:,:,2:NK2,nt))
          ENDDO
        END SELECT
        endif
      ENDIF
!---------------------------------
! Z3DX
!---------------------------------
      IF (PRESENT(Z3DX)) THEN
        if (Z3DX) then

        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1-1,                            &
                            channel%seg(1)%Z3DX(:,:,2:NK1),            &
                            channel%seg(2)%Z3DX(:,:,2:NK1),            &
                            channel%seg(3)%Z3DX(:,:,2:NK1),            &
                            channel%seg(4)%Z3DX(:,:,2:NK1))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1-1,                            &
                            channel%seg(1)%Z3DX(:,:,2:NK1),            &
                            channel%seg(2)%Z3DX(:,:,2:NK1),            &
                            channel%seg(3)%Z3DX(:,:,2:NK1),            &
                            channel%seg(4)%Z3DX(:,:,2:NK1))

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1-1,                            &
                            channel%seg(1)%Z3DX(:,:,2:NK1),            &
                            channel%seg(2)%Z3DX(:,:,2:NK1),            &
                            channel%seg(3)%Z3DX(:,:,2:NK1),            &
                            channel%seg(4)%Z3DX(:,:,2:NK1))
        END SELECT

        endif
      ENDIF
!---------------------------------
! Z3DY
!---------------------------------
      IF (PRESENT(Z3DY)) THEN
        if (Z3DY) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1-1,                            &
                            channel%seg(1)%Z3DY(:,:,2:NK1),            &
                            channel%seg(2)%Z3DY(:,:,2:NK1),            &
                            channel%seg(3)%Z3DY(:,:,2:NK1),            &
                            channel%seg(4)%Z3DY(:,:,2:NK1))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1-1,                            &
                            channel%seg(1)%Z3DY(:,:,2:NK1),            &
                            channel%seg(2)%Z3DY(:,:,2:NK1),            &
                            channel%seg(3)%Z3DY(:,:,2:NK1),            &
                            channel%seg(4)%Z3DY(:,:,2:NK1))

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1-1,                            &
                            channel%seg(1)%Z3DY(:,:,2:NK1),            &
                            channel%seg(2)%Z3DY(:,:,2:NK1),            &
                            channel%seg(3)%Z3DY(:,:,2:NK1),            &
                            channel%seg(4)%Z3DY(:,:,2:NK1))
        END SELECT
        endif
      ENDIF
!---------------------------------
! Z3DZ
!---------------------------------
      IF (PRESENT(Z3DZ)) THEN
        if (Z3DZ) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,1,                                &
                            channel%seg(1)%Z3DZ(:,:,NK2),              &
                            channel%seg(2)%Z3DZ(:,:,NK2),              &
                            channel%seg(3)%Z3DZ(:,:,NK2),              &
                            channel%seg(4)%Z3DZ(:,:,NK2))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,1,                                &
                            channel%seg(1)%Z3DZ(:,:,NK2),              &
                            channel%seg(2)%Z3DZ(:,:,NK2),              &
                            channel%seg(3)%Z3DZ(:,:,NK2),              &
                            channel%seg(4)%Z3DZ(:,:,NK2),OVERLAP=.TRUE.)

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,1,                                &
                            channel%seg(1)%Z3DZ(:,:,NK2),              &
                            channel%seg(2)%Z3DZ(:,:,NK2),              &
                            channel%seg(3)%Z3DZ(:,:,NK2),              &
                            channel%seg(4)%Z3DZ(:,:,NK2),OVERLAP=.TRUE.)
        END SELECT
        endif
      ENDIF

!---------------------------------
! W3D
!---------------------------------
      IF (PRESENT(W3D)) THEN
        if (W3D) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1-1,                            &
                            channel%seg(1)%W3D(:,:,2:NK1),             &
                            channel%seg(2)%W3D(:,:,2:NK1),             &
                            channel%seg(3)%W3D(:,:,2:NK1),             &
                            channel%seg(4)%W3D(:,:,2:NK1))
        CASE(2)
          CALL HALO_GROUP2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1-1,                            &
                            channel%seg(1)%W3D(:,:,2:NK1),             &
                            channel%seg(2)%W3D(:,:,2:NK1),             &
                            channel%seg(3)%W3D(:,:,2:NK1),             &
                            channel%seg(4)%W3D(:,:,2:NK1))

        CASE(3)
          CALL HALO_GROUP3 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1), &
                            num_halo,NK1-1,                            &
                            channel%seg(1)%W3D(:,:,2:NK1),             &
                            channel%seg(2)%W3D(:,:,2:NK1),             &
                            channel%seg(3)%W3D(:,:,2:NK1),             &
                            channel%seg(4)%W3D(:,:,2:NK1))
        END SELECT
        endif
      ENDIF

!---------------------------------
! PSI
!---------------------------------
      IF (PRESENT(PSI)) THEN
        if (PSI) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%PSI,channel%seg(2)%PSI,              &
                              channel%seg(3)%PSI,channel%seg(4)%PSI)
        CASE(2)
          CALL HALO_GROUP2_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%PSI,channel%seg(2)%PSI,              &
                              channel%seg(3)%PSI,channel%seg(4)%PSI,OVERLAP=.TRUE.)

        CASE(3)
          CALL HALO_GROUP3_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%PSI,channel%seg(2)%PSI,              &
                              channel%seg(3)%PSI,channel%seg(4)%PSI,OVERLAP=.TRUE.)
        END SELECT
        endif
      ENDIF
!---------------------------------
! CHI
!---------------------------------
      IF (PRESENT(CHI)) THEN
        if (CHI) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%CHI,channel%seg(2)%CHI,              &
                              channel%seg(3)%CHI,channel%seg(4)%CHI)
        CASE(2)
          CALL HALO_GROUP2_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%CHI,channel%seg(2)%CHI,              &
                              channel%seg(3)%CHI,channel%seg(4)%CHI)

        CASE(3)
          CALL HALO_GROUP3_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%CHI,channel%seg(2)%CHI,              &
                              channel%seg(3)%CHI,channel%seg(4)%CHI)
        END SELECT
        endif
      ENDIF
!---------------------------------
! TOPOZ
!---------------------------------
      IF (PRESENT(TOPOZ)) THEN
        if (TOPOZ) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%TOPOZ,channel%seg(2)%TOPOZ,          &
                              channel%seg(3)%TOPOZ,channel%seg(4)%TOPOZ)
        CASE(2)
          CALL HALO_GROUP2_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%TOPOZ,channel%seg(2)%TOPOZ,          &
                              channel%seg(3)%TOPOZ,channel%seg(4)%TOPOZ)

        CASE(3)
          CALL HALO_GROUP3_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%TOPOZ,channel%seg(2)%TOPOZ,          &
                              channel%seg(3)%TOPOZ,channel%seg(4)%TOPOZ)
        END SELECT
        endif
      ENDIF
!---------------------------------
! ZROUGH
!---------------------------------
      IF (PRESENT(ZROUGH)) THEN
        if (ZROUGH) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%ZROUGH,channel%seg(2)%ZROUGH,        &
                              channel%seg(3)%ZROUGH,channel%seg(4)%ZROUGH)
        CASE(2)
          CALL HALO_GROUP2_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%ZROUGH,channel%seg(2)%ZROUGH,        &
                              channel%seg(3)%ZROUGH,channel%seg(4)%ZROUGH)

        CASE(3)
          CALL HALO_GROUP3_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%ZROUGH,channel%seg(2)%ZROUGH,        &
                              channel%seg(3)%ZROUGH,channel%seg(4)%ZROUGH)
        END SELECT
        endif
      ENDIF
!---------------------------------
! GWET
!---------------------------------
      IF (PRESENT(GWET)) THEN
        if (GWET) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%GWET,channel%seg(2)%GWET,            &
                              channel%seg(3)%GWET,channel%seg(4)%GWET)
        CASE(2)
          CALL HALO_GROUP2_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%GWET,channel%seg(2)%GWET,            &
                              channel%seg(3)%GWET,channel%seg(4)%GWET)

        CASE(3)
          CALL HALO_GROUP3_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%GWET,channel%seg(2)%GWET,            &
                              channel%seg(3)%GWET,channel%seg(4)%GWET)
        END SELECT
        endif
      ENDIF
!---------------------------------
! TG
!---------------------------------
      IF (PRESENT(TG)) THEN
        if (TG) then
        SELECT CASE (channel%num_chg)
        CASE(1)
          CALL HALO_GROUP1_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%TG,channel%seg(2)%TG,                &
                              channel%seg(3)%TG,channel%seg(4)%TG)
        CASE(2)
          CALL HALO_GROUP2_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%TG,channel%seg(2)%TG,                &
                              channel%seg(3)%TG,channel%seg(4)%TG)

        CASE(3)
          CALL HALO_GROUP3_2 (mi1(1),mim(1),mip(1),mj1(1),mjm(1),mjp(1),num_halo, &
                              channel%seg(1)%TG,channel%seg(2)%TG,                &
                              channel%seg(3)%TG,channel%seg(4)%TG)
        END SELECT
        endif
      ENDIF

      END SUBROUTINE bound_channel

!=========================================================================================
      SUBROUTINE HALO_GROUP1 (mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1, &
                              num_halo,ksize,a1,a2,a3,a4)
!=========================================================================================
!     specify the halo values of GROUP 1 channel
!     (SEG1: Face 4, SEG2: Face 1, SEG3: Face 2, SEG4: Face 3)

      INTEGER, INTENT(IN) :: mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1,num_halo,ksize
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1,ksize), INTENT(INOUT) :: A1
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1,ksize), INTENT(INOUT) :: A2
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1,ksize), INTENT(INOUT) :: A3
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1,ksize), INTENT(INOUT) :: A4

      INTEGER :: chl,chl_p,chl_m,chl_rev,chw,chwm,chwp
      INTEGER :: nn,chn,k

!--------------------------------------------------------
! Seg 1, 2, 3, and 4 (same domain shape and size, x-array)
!--------------------------------------------------------
      chl  = mi1_1
      chw  = mj1_1
      chwm = mjm_1
      chwp = mjp_1

      DO nn = 1,num_halo
        chl_p   = chl + nn
        chl_m   = 1 - nn
        chl_rev = chl - nn + 1
        DO K = 1, ksize
        DO chn = chwm, chwp
          A1(chl_p,chn,K) = A2(nn,chn,K)
          A1(chl_m,chn,K) = A4(chl_rev,chn,K)

          A2(chl_p,chn,K) = A3(nn,chn,K)
          A2(chl_m,chn,K) = A1(chl_rev,chn,K)

          A3(chl_p,chn,K) = A4(nn,chn,K)
          A3(chl_m,chn,K) = A2(chl_rev,chn,K)

          A4(chl_p,chn,K) = A1(nn,chn,K)
          A4(chl_m,chn,K) = A3(chl_rev,chn,K)
        ENDDO
        ENDDO
      ENDDO

      end subroutine halo_group1

!=========================================================================================
      SUBROUTINE HALO_GROUP1_2 (mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1, &
                                num_halo,a1,a2,a3,a4)
!=========================================================================================
!     Same as HALO_GROUP1, but used for a 2D variable.

      INTEGER, INTENT(IN) :: mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1,num_halo
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1), INTENT(INOUT) :: A1
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1), INTENT(INOUT) :: A2
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1), INTENT(INOUT) :: A3
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1), INTENT(INOUT) :: A4

      INTEGER :: chl,chl_p,chl_m,chl_rev,chw,chwm,chwp
      INTEGER :: nn,chn

!--------------------------------------------------------
! Seg 1, 2, 3, and 4 (same domain shape and size, x-array)
!--------------------------------------------------------
      chl  = mi1_1
      chw  = mj1_1
      chwm = mjm_1
      chwp = mjp_1

      DO nn = 1,num_halo
        chl_p   = chl + nn
        chl_m   = 1 - nn
        chl_rev = chl - nn + 1
        DO chn = chwm, chwp
          A1(chl_p,chn) = A2(nn,chn)
          A1(chl_m,chn) = A4(chl_rev,chn)

          A2(chl_p,chn) = A3(nn,chn)
          A2(chl_m,chn) = A1(chl_rev,chn)

          A3(chl_p,chn) = A4(nn,chn)
          A3(chl_m,chn) = A2(chl_rev,chn)

          A4(chl_p,chn) = A1(nn,chn)
          A4(chl_m,chn) = A3(chl_rev,chn)
        ENDDO
      ENDDO

      end subroutine halo_group1_2

!=========================================================================================
      SUBROUTINE HALO_GROUP2 (mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1, &
                              num_halo,ksize,a1,a2,a3,a4,OVERLAP)
!=========================================================================================
!     specify the halo values of GROUP 2 channel
!     (SEG1: Face 6, SEG2: Face 1, SEG3: Face 5, SEG4: Face 3)

      INTEGER, INTENT(IN) :: mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1,num_halo,ksize
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1,ksize), INTENT(INOUT) :: A1
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1,ksize), INTENT(INOUT) :: A2
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1,ksize), INTENT(INOUT) :: A3
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1,ksize), INTENT(INOUT) :: A4

      LOGICAL, INTENT(IN), OPTIONAL :: OVERLAP

      ! Local
      INTEGER :: chl,chl_p,chl_m,chl_rev,chw,chwm,chwp
      INTEGER :: nn,chn,chn_rev,k

!--------------------------------------------------------
! Seg 1, 2, 3, and 4 (same domain shape and size, y-array)
!--------------------------------------------------------
      chl  = mj1_1
      chw  = mi1_1
      chwm = mim_1
      chwp = mip_1

      DO nn = 1,num_halo
        chl_p   = chl + nn
        chl_m   = 1 - nn
        chl_rev = chl - nn + 1

        DO K = 1, ksize
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A1(chn,chl_p,K) = A2(chn,nn,K)
          A1(chn,chl_m,K) = A4(chn_rev,nn,K)

          A2(chn,chl_p,K) = A3(chn,nn,K)
          A2(chn,chl_m,K) = A1(chn,chl_rev,K)

          A3(chn,chl_m,K) = A2(chn,chl_rev,K)
          A4(chn,chl_m,K) = A1(chn_rev,nn,K)
        ENDDO
        ENDDO

      IF (PRESENT(OVERLAP)) THEN
        if (OVERLAP) then
        DO K = 1, ksize
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn
          A3(chn,chl_p,K) = A4(chn_rev,chl_rev-1,K)   ! For zeta & psi
          A4(chn,chl_p,K) = A3(chn_rev,chl_rev-1,K)   ! For zeta & psi
        ENDDO
        ENDDO
        endif
      ELSE
        DO K = 1, ksize
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn
          A3(chn,chl_p,K) = A4(chn_rev,chl_rev,K)
          A4(chn,chl_p,K) = A3(chn_rev,chl_rev,K)
        ENDDO
        ENDDO
      ENDIF

      ENDDO   ! nn-loop

      end subroutine halo_group2

!=========================================================================================
      SUBROUTINE HALO_GROUP2_2 (mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1, &
                                num_halo,a1,a2,a3,a4,OVERLAP)
!=========================================================================================
!     Same as HALO_GROUP2, but used for a 2D variable.

      INTEGER, INTENT(IN) :: mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1,num_halo
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1), INTENT(INOUT) :: A1
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1), INTENT(INOUT) :: A2
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1), INTENT(INOUT) :: A3
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1), INTENT(INOUT) :: A4

      LOGICAL, INTENT(IN), OPTIONAL :: OVERLAP

      ! Local
      INTEGER :: chl,chl_p,chl_m,chl_rev,chw,chwm,chwp
      INTEGER :: nn,chn,chn_rev

!--------------------------------------------------------
! Seg 1, 2, 3, and 4 (same domain shape and size, y-array)
!--------------------------------------------------------
      chl  = mj1_1
      chw  = mi1_1
      chwm = mim_1
      chwp = mip_1

      DO nn = 1,num_halo
        chl_p   = chl + nn
        chl_m   = 1 - nn
        chl_rev = chl - nn + 1

        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A1(chn,chl_p) = A2(chn,nn)
          A1(chn,chl_m) = A4(chn_rev,nn)

          A2(chn,chl_p) = A3(chn,nn)
          A2(chn,chl_m) = A1(chn,chl_rev)

          A3(chn,chl_m) = A2(chn,chl_rev)
          A4(chn,chl_m) = A1(chn_rev,nn)
        ENDDO

      IF (PRESENT(OVERLAP)) THEN
        if (OVERLAP) then
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn
          A3(chn,chl_p) = A4(chn_rev,chl_rev-1)   ! For zeta & psi
          A4(chn,chl_p) = A3(chn_rev,chl_rev-1)   ! For zeta & psi
        ENDDO
        endif
      ELSE
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn
          A3(chn,chl_p) = A4(chn_rev,chl_rev)
          A4(chn,chl_p) = A3(chn_rev,chl_rev)
        ENDDO
      ENDIF

      ENDDO   ! nn-loop

      end subroutine halo_group2_2

!=========================================================================================
      SUBROUTINE HALO_GROUP3 (mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1, &
                              num_halo,ksize,a1,a2,a3,a4,OVERLAP)
!=========================================================================================
!     specify the halo values of GROUP 3 channel
!     (SEG1: Face 4, SEG2: Face 5, SEG3: Face 2, SEG4: Face 6)

      INTEGER, INTENT(IN) :: mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1,num_halo,ksize
      
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1,ksize), INTENT(INOUT) :: A1
      REAL (kind=r8), DIMENSION(mjm_1:mjp_1,mim_1:mip_1,ksize), INTENT(INOUT) :: A2
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1,ksize), INTENT(INOUT) :: A3
      REAL (kind=r8), DIMENSION(mjm_1:mjp_1,mim_1:mip_1,ksize), INTENT(INOUT) :: A4

      LOGICAL, INTENT(IN), OPTIONAL :: OVERLAP

      ! Local
      INTEGER :: chl,chl_p,chl_m,chl_rev,chw,chwm,chwp
      INTEGER :: nn,chn,chn_rev,k

!----------------------------
      ! Seg 1 (y-array)
!----------------------------
      chl = mj1_1

      chw  = mi1_1
      chwm = mim_1
      chwp = mip_1
      DO nn = 1,num_halo
        chl_p = chl + nn
        chl_m = 1 - nn

        DO K = 1, ksize
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A1(chn,chl_p,K) = A2(nn,chn_rev,K)
          A1(chn,chl_m,K) = A4(nn,chn,K)

        ENDDO
        ENDDO
      ENDDO
!----------------------------
      ! Seg 3 (y-array)
!----------------------------
      chl = mj1_1

      chw  = mi1_1
      chwm = mim_1
      chwp = mip_1
      DO nn = 1,num_halo
        chl_p   = chl + nn
        chl_m   = 1 - nn
        
        chl_rev = chl - nn + 1

      IF (PRESENT(OVERLAP)) THEN
        if (OVERLAP) then
        DO K = 1, ksize
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A3(chn,chl_p,K) = A2(chl_rev-1,chn,K)      ! for zeta & psi
          A3(chn,chl_m,K) = A4(chl_rev,chn_rev,K)

        ENDDO
        ENDDO
        endif
      ELSE
        DO K = 1, ksize
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A3(chn,chl_p,K) = A2(chl_rev,chn,K)
          A3(chn,chl_m,K) = A4(chl_rev,chn_rev,K)

        ENDDO
        ENDDO
      ENDIF

      ENDDO  ! nn-loop
!----------------------------
      ! Seg 2 (x-array)
!----------------------------
      chl = mj1_1

      chw  = mi1_1
      chwm = mim_1
      chwp = mip_1
      DO nn = 1,num_halo
        chl_p   = chl + nn
        chl_m   = 1 - nn
        chl_rev = chl - nn + 1

      IF (PRESENT(OVERLAP)) THEN
        if (OVERLAP) then
        DO K = 1, ksize
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A2(chl_p,chn,K) = A3(chn,chl_rev-1,K)     ! For zeta & psi
          A2(chl_m,chn,K) = A1(chn_rev,chl_rev,K)

        ENDDO
        ENDDO
        endif
      ELSE
        DO K = 1, ksize
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A2(chl_p,chn,K) = A3(chn,chl_rev,K)
          A2(chl_m,chn,K) = A1(chn_rev,chl_rev,K)

        ENDDO
        ENDDO
      ENDIF

      ENDDO  ! nn-loop
!----------------------------
      ! Seg 4 (x-array)
!----------------------------
      chl = mj1_1

      chw  = mi1_1
      chwm = mim_1
      chwp = mip_1
      DO nn = 1,num_halo
        chl_p   = chl + nn
        chl_m   = 1 - nn

        DO K = 1, ksize
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A4(chl_p,chn,K) = A3(chn_rev,nn,K)
          A4(chl_m,chn,K) = A1(chn,nn,K)

        ENDDO
        ENDDO
      ENDDO  ! nn-loop

      end subroutine halo_group3

!=========================================================================================
      SUBROUTINE HALO_GROUP3_2 (mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1, &
                                num_halo,a1,a2,a3,a4,OVERLAP)
!=========================================================================================
!     Same as HALO_GROUP3, but used for a 2D variable.

!     specify the halo values of GROUP 3 channel
!     (SEG1: Face 4, SEG2: Face 5, SEG3: Face 2, SEG4: Face 6)

      INTEGER, INTENT(IN) :: mi1_1,mim_1,mip_1,mj1_1,mjm_1,mjp_1,num_halo
      
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1), INTENT(INOUT) :: A1
      REAL (kind=r8), DIMENSION(mjm_1:mjp_1,mim_1:mip_1), INTENT(INOUT) :: A2
      REAL (kind=r8), DIMENSION(mim_1:mip_1,mjm_1:mjp_1), INTENT(INOUT) :: A3
      REAL (kind=r8), DIMENSION(mjm_1:mjp_1,mim_1:mip_1), INTENT(INOUT) :: A4

      LOGICAL, INTENT(IN), OPTIONAL :: OVERLAP

      ! Local
      INTEGER :: chl,chl_p,chl_m,chl_rev,chw,chwm,chwp
      INTEGER :: nn,chn,chn_rev

!----------------------------
      ! Seg 1 (y-array)
!----------------------------
      chl = mj1_1

      chw  = mi1_1
      chwm = mim_1
      chwp = mip_1
      DO nn = 1,num_halo
        chl_p = chl + nn
        chl_m = 1 - nn

        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A1(chn,chl_p) = A2(nn,chn_rev)
          A1(chn,chl_m) = A4(nn,chn)

        ENDDO

      ENDDO
!----------------------------
      ! Seg 3 (y-array)
!----------------------------
      chl = mj1_1

      chw  = mi1_1
      chwm = mim_1
      chwp = mip_1
      DO nn = 1,num_halo
        chl_p   = chl + nn
        chl_m   = 1 - nn
        
        chl_rev = chl - nn + 1

      IF (PRESENT(OVERLAP)) THEN
        if (OVERLAP) then

        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A3(chn,chl_p) = A2(chl_rev-1,chn)      ! for zeta & psi
          A3(chn,chl_m) = A4(chl_rev,chn_rev)

        ENDDO

        endif
      ELSE
        
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A3(chn,chl_p) = A2(chl_rev,chn)
          A3(chn,chl_m) = A4(chl_rev,chn_rev)

        ENDDO
        
      ENDIF

      ENDDO  ! nn-loop
!----------------------------
      ! Seg 2 (x-array)
!----------------------------
      chl = mj1_1

      chw  = mi1_1
      chwm = mim_1
      chwp = mip_1
      DO nn = 1,num_halo
        chl_p   = chl + nn
        chl_m   = 1 - nn
        chl_rev = chl - nn + 1

      IF (PRESENT(OVERLAP)) THEN
        if (OVERLAP) then
        
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A2(chl_p,chn) = A3(chn,chl_rev-1)     ! For zeta & psi
          A2(chl_m,chn) = A1(chn_rev,chl_rev)

        ENDDO
        
        endif
      ELSE
        
        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A2(chl_p,chn) = A3(chn,chl_rev)
          A2(chl_m,chn) = A1(chn_rev,chl_rev)

        ENDDO
        
      ENDIF

      ENDDO  ! nn-loop
!----------------------------
      ! Seg 4 (x-array)
!----------------------------
      chl = mj1_1

      chw  = mi1_1
      chwm = mim_1
      chwp = mip_1
      DO nn = 1,num_halo
        chl_p   = chl + nn
        chl_m   = 1 - nn

        DO chn = chwm, chwp
          chn_rev = chw + 1 - chn

          A4(chl_p,chn) = A3(chn_rev,nn)
          A4(chl_m,chn) = A1(chn,nn)

        ENDDO
        
      ENDDO  ! nn-loop

      end subroutine halo_group3_2

END MODULE bound_channel_module
