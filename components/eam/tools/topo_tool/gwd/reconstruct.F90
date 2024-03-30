!-----------------------------------------------------------------------------
! MODULE GECORE_Interp_ABP
!
! Purpose:
!   Provides functions for performing conservative interpolation between
!   cubed sphere and lat lon grids.
!
!      Date       Programmer       Affiliation          Description of change
!      ====       ==========       ===========          =====================
!    06/25/08    P.A.Ullrich       CGD,NCAR             Original Code
!    07/11/08    P.A.Ullrich       CGD,NCAR             First order interp.
!    07/18/08    P.A.Ullrich       CGD,NCAR             Renamed to GECORE
!
!-----------------------------------------------------------------------------

MODULE reconstruct
  USE remap
!  INTEGER, PARAMETER ::                           &
!       int_kind  = KIND(1),                       &
!       real_kind = SELECTED_REAL_KIND(p=14,r=100),&
!       dbl_kind  = selected_real_kind(13)        

!  LOGICAL, PARAMETER:: ldbg_r = .FALSE.


  INTEGER, PRIVATE :: ncube_reconstruct
  REAL(kind=real_kind), PARAMETER, PRIVATE ::          &
!       one = 1.0                             ,&
!       aa  = 1.0                             ,&
!       tiny= 1.0E-10                         ,&
       pi  = 3.14159265358979323846264338327 ,&
       piq = 0.25*pi                         ,&
       pih = 0.50*pi                         ,&
       pi2 = 2.0*pi 

          
           
  REAL    (KIND=dbl_kind), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: abp_centroid
  REAL    (KIND=dbl_kind), DIMENSION(:), ALLOCATABLE, PRIVATE :: gp
  REAL    (KIND=dbl_kind), PARAMETER, PRIVATE :: lammax = 60_dbl_kind !for selective limiting
!  REAL    (KIND=dbl_kind),PARAMETER,PRIVATE :: bignum = 1.e+20_dbl_kind



CONTAINS
  !
  ! wrapper subroutine to use the CSLAM code
  !
  ! note that ncube_reconstruct is ncube_reconstruct+1 in topo software
  !
  SUBROUTINE get_reconstruction(fcube_in,order, kmono, recons, kpd,DAcube,ncube,nreconstruction,centroids)
    IMPLICIT NONE

    REAL (KIND=dbl_kind), &
            DIMENSION(6*(ncube-1)*(ncube-1)), INTENT(IN) :: fcube_in

    INTEGER (KIND=int_kind), INTENT(IN) :: order, kmono, kpd,ncube,nreconstruction

    REAL (KIND=dbl_kind), DIMENSION(nreconstruction,6*(ncube-1)*(ncube-1)), INTENT(OUT) :: recons
    REAL (KIND=dbl_kind), DIMENSION(ncube-1,ncube-1), INTENT(IN) :: dAcube
    REAL (KIND=dbl_kind), DIMENSION(nreconstruction,6*(ncube-1)*(ncube-1)), INTENT(OUT) :: centroids
    ! Local variables
    REAL (KIND=dbl_kind), &
            DIMENSION(1:ncube-1, 1:ncube-1, 6) :: fcube
    INTEGER (KIND=int_kind) :: status,k,ip,ii,ix,iy

    REAL (KIND=dbl_kind), DIMENSION(-1:ncube+1, -1:ncube+1, 6) :: fcubehalo
    REAL (KIND=dbl_kind), DIMENSION(nreconstruction,ncube-1,ncube-1,6) :: recons_4D

    IF (order<2) THEN
      WRITE(*,*) "order out of range",order
      STOP
    END IF

    ncube_reconstruct = ncube

    DO ip=1,6
      DO iy=1,ncube-1
        DO ix=1,ncube-1
          ii = (ip-1)*(ncube-1)*(ncube-1)+(iy-1)*(ncube-1)+ix
          fcube(ix,iy,ip) = fcube_in(ii)
        END DO
      END DO
    END DO


    ALLOCATE(gp (ncube_reconstruct), STAT=status)
    ! Equi-angular grid in [-piq, piq]
    DO k = 1, ncube_reconstruct
      gp(k) = -piq + (pi / DBLE(2 * (ncube_reconstruct-1))) * DBLE(k-1)
    ENDDO
    
    WRITE(*,*) "Compute centroids"
    CALL ComputeABPElementCentroids(order,DAcube,ncube_reconstruct)
!    WRITE(*,*) "abp_centroid(1,1,1)",abp_centroid(1,1,1)
    WRITE(*,*) "Compute reconstruction coefficients"
    CALL ReconstructABPGradient(fcube, 3, 2, order, kmono, recons_4D, kpd, 2)
    DEALLOCATE(gp)
    WRITE(*,*) "map to 1D arrays"
    recons = 0.0
    DO ip=1,6
      DO iy=1,ncube-1
        DO ix=1,ncube-1
          ii = (ip-1)*(ncube-1)*(ncube-1)+(iy-1)*(ncube-1)+ix
          recons(:,ii) = recons_4D(:,ix,iy,ip)     
          centroids(:,ii) = abp_centroid(:,ix,iy)
        END DO
      END DO
    END DO
    DEALLOCATE(abp_centroid)
  END SUBROUTINE get_reconstruction

!------------------------------------------------------------------------------
! SUBROUTINE ComputeABPElementCentroids
!
! Description:
!   Compute the centroid coordinates of ABP elements.
!
! Note:
!   ComputeABPW8.0s_1 must be called prior to calling this function.
!
! Parameters:
!   order - Order of the method being applied
!------------------------------------------------------------------------------
  SUBROUTINE ComputeABPElementCentroids(order,DAcube,ncube_reconstruct)
    IMPLICIT NONE

    INTEGER (KIND=int_kind)                         , INTENT(IN) :: order,ncube_reconstruct
    REAL (KIND=dbl_kind), DIMENSION(ncube_reconstruct-1,ncube_reconstruct-1), INTENT(IN) :: dAcube

    INTEGER (KIND=int_kind) :: i, j, k
    REAL    (KIND=dbl_kind) :: alpha, beta, alpha_next, beta_next, area, lint1, lint2

    ! Equi-angular grid in [-piq, piq]
    DO k = 1, ncube_reconstruct
      gp(k) = -piq + (pi / DBLE(2 * (ncube_reconstruct-1))) * DBLE(k-1)
    ENDDO

    ! Allocate space for centroid calculations
    IF (order == 1) THEN
      RETURN
    ELSEIF (order == 2) THEN
      WRITE(*,*) "allocate abp_centroid"
      ALLOCATE(abp_centroid(2, -1:ncube_reconstruct+1, -1:ncube_reconstruct+1))
    ELSEIF (order == 3) THEN
      WRITE(*,*) "allocate abp_centroid"
      ALLOCATE(abp_centroid(5, -1:ncube_reconstruct+1, -1:ncube_reconstruct+1))
    ELSE
      WRITE(*,*) 'Fatal Error: In ComputeABPElementCentroids'
      WRITE(*,*) 'order out of range [1-3], given: ', order
      STOP
    ENDIF

    WRITE(*,*) 'Begin computing ABP element centroids'

    ! Compute centroids via line integrals
    abp_centroid = 0.0

    DO j = -1, ncube_reconstruct+1
      IF ((j > 0) .AND. (j < ncube_reconstruct)) THEN
        beta = gp(j)
        beta_next = gp(j+1)
      ELSEIF (j == -1) THEN
        beta = -piq - (gp(3) + piq)
        beta_next = -piq - (gp(2) + piq)
      ELSEIF (j == 0) THEN
        beta = -piq - (gp(2) + piq)
        beta_next = -piq
      ELSEIF (j == ncube_reconstruct) THEN
        beta = piq
        beta_next = piq + (piq - gp(ncube_reconstruct-1))
      ELSEIF (j == ncube_reconstruct+1) THEN
        beta = piq + (piq - gp(ncube_reconstruct-1))
        beta_next = piq + (piq - gp(ncube_reconstruct-2))
      ENDIF

      DO i = -1, ncube_reconstruct+1
        IF ((i > 0) .AND. (i < ncube_reconstruct)) THEN
          alpha = gp(i)
          alpha_next = gp(i+1)
        ELSEIF (i == -1) THEN
          alpha = -piq - (gp(3) + piq)
          alpha_next = -piq - (gp(2) + piq)
        ELSEIF (i == 0) THEN
          alpha = -piq - (gp(2) + piq)
          alpha_next = -piq
        ELSEIF (i == ncube_reconstruct) THEN
          alpha = piq
          alpha_next = piq + (piq - gp(ncube_reconstruct-1))
        ELSEIF (i == ncube_reconstruct+1) THEN
          alpha = piq + (piq - gp(ncube_reconstruct-1))
          alpha_next = piq + (piq - gp(ncube_reconstruct-2))
        ENDIF
        abp_centroid(1,i,j) =                             &
               I_10_ab(alpha_next,beta_next)-I_10_ab(alpha     ,beta_next)+&
               I_10_ab(alpha     ,beta     )-I_10_ab(alpha_next,beta     )
!          - ASINH(COS(alpha_next) * TAN(beta_next)) &
!          + ASINH(COS(alpha_next) * TAN(beta))         &
!          + ASINH(COS(alpha) * TAN(beta_next))         &
!          - ASINH(COS(alpha) * TAN(beta))

        abp_centroid(2,i,j) =                             &
               I_01_ab(alpha_next,beta_next)-I_01_ab(alpha     ,beta_next)+&
               I_01_ab(alpha     ,beta     )-I_01_ab(alpha_next,beta     )
!          - ASINH(TAN(alpha_next) * COS(beta_next)) &
!          + ASINH(TAN(alpha_next) * COS(beta))         &
!          + ASINH(TAN(alpha) * COS(beta_next))         &
!          - ASINH(TAN(alpha) * COS(beta))

        !ADD PHL START
        IF (order>2) THEN
           ! TAN(alpha)^2 component
           abp_centroid(3,i,j)  =&
                I_20_ab(alpha_next,beta_next)-I_20_ab(alpha     ,beta_next)+&
                I_20_ab(alpha     ,beta     )-I_20_ab(alpha_next,beta     )

           ! TAN(beta)^2 component
           abp_centroid(4,i,j) = &
                I_02_ab(alpha_next,beta_next)-I_02_ab(alpha     ,beta_next)+&
                I_02_ab(alpha     ,beta     )-I_02_ab(alpha_next,beta     )

           ! TAN(alpha) TAN(beta) component
           abp_centroid(5,i,j) = &
                I_11_ab(alpha_next,beta_next)-I_11_ab(alpha     ,beta_next)+&
                I_11_ab(alpha     ,beta     )-I_11_ab(alpha_next,beta     )
        ENDIF
        !ADD PHL END
     ENDDO
  ENDDO

!
! PHL outcommented below
!
    ! High order calculations
!    IF (order > 2) THEN
!      DO k = 1, nlon
!        DO i = 1, int_nx(nlat,k)-1
!          IF ((int_itype(i,k) > 4) .AND. (int_np(1,i,k) == 1)) THEN
!            abp_centroid(3, int_a(i,k), int_b(i,k)) =                  &
!              abp_centroid(3, int_a(i,k), int_b(i,k)) + int_wt_2a(i,k)
!            abp_centroid(4, int_a(i,k), int_b(i,k)) =                  &
!              abp_centroid(4, int_a(i,k), int_b(i,k)) + int_wt_2b(i,k)
!            abp_centroid(5, int_a(i,k), int_b(i,k)) =                  &
!              abp_centroid(5, int_a(i,k), int_b(i,k)) + int_wt_2c(i,k)
!          ENDIF
!        ENDDO
!      ENDDO
!    ENDIF

    ! Normalize with element areas
    DO j = -1, ncube_reconstruct+1
      IF ((j > 0) .AND. (j < ncube_reconstruct)) THEN
        beta = gp(j)
        beta_next = gp(j+1)
      ELSEIF (j == -1) THEN
        beta = -piq - (gp(3) + piq)
        beta_next = -piq - (gp(2) + piq)
      ELSEIF (j == 0) THEN
        beta = -piq - (gp(2) + piq)
        beta_next = -piq
      ELSEIF (j == ncube_reconstruct) THEN
        beta = piq
        beta_next = piq + (piq - gp(ncube_reconstruct-1))
      ELSEIF (j == ncube_reconstruct+1) THEN
        beta = piq + (piq - gp(ncube_reconstruct-1))
        beta_next = piq + (piq - gp(ncube_reconstruct-2))
      ENDIF
      DO i = -1, ncube_reconstruct+1
        IF ((i > 0) .AND. (i < ncube_reconstruct)) THEN
          alpha = gp(i)
          alpha_next = gp(i+1)
        ELSEIF (i == -1) THEN
          alpha = -piq - (gp(3) + piq)
          alpha_next = -piq - (gp(2) + piq)
        ELSEIF (i == 0) THEN
          alpha = -piq - (gp(2) + piq)
          alpha_next = -piq
        ELSEIF (i == ncube_reconstruct) THEN
          alpha = piq
          alpha_next = piq + (piq - gp(ncube_reconstruct-1))
        ELSEIF (i == ncube_reconstruct+1) THEN
          alpha = piq + (piq - gp(ncube_reconstruct-1))
          alpha_next = piq + (piq - gp(ncube_reconstruct-2))
        ENDIF

        IF ((i > 0) .AND. (i < ncube_reconstruct) .AND. (j > 0) .AND. (j < ncube_reconstruct)) THEN
          area = DAcube(i,j)
        ELSE
          area = EquiangularElementArea(alpha, alpha_next - alpha,  &
                                        beta,  beta_next - beta)
        ENDIF

        abp_centroid(1,i,j) = abp_centroid(1,i,j) / area
        abp_centroid(2,i,j) = abp_centroid(2,i,j) / area

        IF (order > 2) THEN
          IF ((i > 0) .AND. (i < ncube_reconstruct) .AND. (j > 0) .AND. (j < ncube_reconstruct)) THEN
            abp_centroid(3,i,j) = abp_centroid(3,i,j) / area
            abp_centroid(4,i,j) = abp_centroid(4,i,j) / area
            abp_centroid(5,i,j) = abp_centroid(5,i,j) / area
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    WRITE(*,*) '...Done computing ABP element centroids'

  END SUBROUTINE ComputeABPElementCentroids

!------------------------------------------------------------------------------
! FUNCTION EvaluateABPReconstruction
!
! Description:
!   Evaluate the sub-grid scale reconstruction at the given point.
!
! Parameters:
!   fcubehalo - Array of element values
!   recons - Array of reconstruction coefficients
!   a - Index of element in alpha direction (1 <= a <= ncube_reconstruct-1)
!   b - Index of element in beta direction (1 <= b <= ncube_reconstruct-1)
!   p - Panel index of element
!   alpha - Alpha coordinate of evaluation point
!   beta - Beta coordinate of evaluation point
!   order - Order of the reconstruction
!   value (OUT) - Result of function evaluation at given point
!------------------------------------------------------------------------------
  SUBROUTINE EvaluateABPReconstruction( &
               fcubehalo, recons, a, b, p, alpha, beta, order, value)
    IMPLICIT NONE

    ! Dummy variables
    REAL (KIND=dbl_kind), DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6), &
                          INTENT(IN) :: fcubehalo

    REAL    (KIND=dbl_kind), DIMENSION(:,:,:,:), INTENT(IN) :: recons
    INTEGER (KIND=int_kind), INTENT(IN) :: a, b, p
    REAL    (KIND=dbl_kind), INTENT(IN) :: alpha, beta
    INTEGER (KIND=int_kind), INTENT(IN) :: order

    REAL    (KIND=dbl_kind), INTENT(OUT) :: value

    ! Evaluate constant order terms
    value = fcubehalo(a,b,p)

    ! Evaluate linear order terms
    IF (order > 1) THEN
      value = value + recons(1,a,b,p) * (TAN(alpha) - abp_centroid(1,a,b))
      value = value + recons(2,a,b,p) * (TAN(beta) - abp_centroid(2,a,b))
    ENDIF

    ! Evaluate second order terms
    IF (order > 2) THEN
      value = value + recons(3,a,b,p) *                              &
                      (abp_centroid(1,a,b)**2 - abp_centroid(3,a,b))
      value = value + recons(4,a,b,p) *                              &
                      (abp_centroid(2,a,b)**2 - abp_centroid(4,a,b))
      value = value + recons(5,a,b,p) *                              &
                      (abp_centroid(1,a,b) * abp_centroid(2,a,b) -   &
                       abp_centroid(5,a,b))

      value = value + recons(3,a,b,p) * (TAN(alpha) - abp_centroid(1,a,b))**2
      value = value + recons(4,a,b,p) * (TAN(beta) - abp_centroid(2,a,b))**2
      value = value + recons(5,a,b,p) * (TAN(alpha) - abp_centroid(1,a,b)) &
                                      * (TAN(beta) - abp_centroid(2,a,b))
    ENDIF

  END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE ABPHaloMinMax
!
! Description:
!   Calculate the minimum and maximum values of the cell-averaged function
!   around the given element.
!
! Parameters:
!   fcubehalo - Cell-averages for the cubed sphere
!   a - Local element alpha index
!   b - Local element beta index
!   p - Local element panel index
!   min_val (OUT) - Minimum value in the halo
!   max_val (OUT) - Maximum value in the halo
!   nomiddle - whether to not include the middle cell (index a,b) in the search.
!
! NOTE: Since this routine is not vectorized, it will likely be called MANY times.
! To speed things up, make sure to pass the first argument as the ENTIRE original 
!   array, not as a subset of it, since repeatedly cutting up that array and creating
!   an array temporary (on some compilers) is VERY slow.
! ex: 
! CALL APBHaloMinMax(zarg, a, ...) !YES
! CALL ABPHaloMinMax(zarg(-1:ncube_reconstruct+1,-1:ncube_reconstruct+1,:)) !NO -- slow
!------------------------------------------------------------------------------
  SUBROUTINE ABPHaloMinMax(fcubehalo, a, b, p, min_val, max_val, nomiddle)    
    IMPLICIT NONE

    REAL (KIND=dbl_kind), DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6), &
                          INTENT(IN) :: fcubehalo

    INTEGER (KIND=int_kind), INTENT(IN)  :: a, b, p
    REAL    (KIND=dbl_kind), INTENT(OUT) :: min_val, max_val
    LOGICAL, INTENT(IN) :: nomiddle

    ! Local variables
    INTEGER (KIND=int_kind) :: i, j, il, jl, inew, jnew
    REAL    (KIND=dbl_kind) :: value

    min_val = fcubehalo(a,b,p)
    max_val = fcubehalo(a,b,p)
    value   = fcubehalo(a,b,p)

    DO il = a-1,a+1
      DO jl = b-1,b+1

         i = il
         j = jl

         inew = i
         jnew = j

         IF (nomiddle .AND. i==a .AND. j==b) CYCLE

         !Interior
        IF ((i > 0) .AND. (i < ncube_reconstruct) .AND. (j > 0) .AND. (j < ncube_reconstruct)) THEN
           value = fcubehalo(i,j,p)

        ELSE


          !The next 4.0 regions are cases in which a,b themselves lie in the panel's halo, and the cell's "halo" (in this usage the 8.0 cells surrounding it) might wrap around into another part of the halo. This happens for (a,b) = {(1,:0),(ncube_reconstruct-1,:0),(1,ncube_reconstruct:),(ncube_reconstruct-1,ncube_reconstruct:)} and for the transposes thereof ({(:0,1), etc.}). In these cases (i,j) could lie in the "Corners" where nothing should lie. We correct this by moving i,j to its appropriate position on the "facing" halo, and then the remainder of the routine then moves it onto the correct face.
         
101      FORMAT("ERROR cannot find (i,j) = (", I4, ", ", I4, ") for (a,b,p) = ", I4, ",", I4, ",", I4, ")")
102      FORMAT("i,j,p = ", 3I4, " moved to " 2I4, " (CASE ", I1, ")")
         !NOTE: we need the general case to be able to properly handle (0,0), (ncube_reconstruct,0), etc. Note that we don't need to bother with (0,0), etc. when a, b lie in the interior, since both sides of the (0,0) cell are already accounted for by this routine.
         !LOWER LEFT
         IF (i < 1 .AND. j < 1) THEN
            IF (a < 1) THEN !(a,b) centered on left halo, cross to lower halo
               inew = 1-j
               jnew = i
            ELSE IF (b < 1) THEN !(a,b) centered on lower halo, cross to left halo
               jnew = 1-i
               inew = j
            END IF
!            WRITE(*,102) i, j, p, inew, jnew, 1
         !LOWER RIGHT
         ELSE IF (i > ncube_reconstruct-1 .AND. j < 1) THEN
            IF (a > ncube_reconstruct-1) THEN !(a,b) centered on right halo, cross to lower halo
               inew = ncube_reconstruct-1+j
               jnew = ncube_reconstruct-i
            ELSE IF (b < 1) THEN !(a,b) centered on lower halo, cross to right halo
               jnew = 1+(i-ncube_reconstruct)
               inew = ncube_reconstruct-j
            END IF
!            WRITE(*,102) i, j, p, inew, jnew, 2
         !UPPER LEFT
         ELSE IF (i < 1 .AND. j > ncube_reconstruct-1) THEN
            IF (a < 1) THEN! (a,b) centered on left halo, cross to upper halo
               inew = 1-(j-ncube_reconstruct)
               jnew = ncube_reconstruct-i
            ELSE IF (b > ncube_reconstruct-1) THEN !(a,b) centered on upper halo, cross to left halo
               inew = ncube_reconstruct-j
               jnew = ncube_reconstruct-1-i
            END IF
!            WRITE(*,102) i, j, p, inew, jnew, 3
         !UPPER RIGHT
         ELSE IF (i > ncube_reconstruct-1 .AND. j > ncube_reconstruct-1) THEN
            IF (a > ncube_reconstruct-1) THEN !(a,b) centered on right halo, cross to upper halo
               inew = ncube_reconstruct-1-(ncube_reconstruct-j)
               jnew = i
            ELSE IF (b > ncube_reconstruct-1) THEN !(a,b) centered on upper halo, cross to right halo
               inew = j
               jnew = ncube_reconstruct-1-(ncube_reconstruct-i)
            END IF  
!            WRITE(*,102) i, j, p, inew, jnew, 4
         END IF

         i = inew
         j = jnew


          !Lower halo ("halo" meaning the panel's halo, not the nine-cell halo
        IF ((i < 1) .AND. (j > 0) .AND. (j < ncube_reconstruct)) THEN
          IF (p == 1) THEN
            value = fcubehalo(ncube_reconstruct-1+i,j,4)
          ELSEIF (p == 2) THEN
            value = fcubehalo(ncube_reconstruct-1+i,j,1)
          ELSEIF (p == 3) THEN
            value = fcubehalo(ncube_reconstruct-1+i,j,2)
          ELSEIF (p == 4) THEN
            value = fcubehalo(ncube_reconstruct-1+i,j,3)
          ELSEIF (p == 5) THEN
            value = fcubehalo(j,1-i,4)
          ELSEIF (p == 6) THEN
            value = fcubehalo(ncube_reconstruct-j,ncube_reconstruct-1+i,4)
          ENDIF

          !Upper halo
        ELSEIF ((i > ncube_reconstruct-1) .AND. (j > 0) .AND. (j < ncube_reconstruct)) THEN
          IF (p == 1) THEN
            value = fcubehalo(i-ncube_reconstruct+1,j,2)
          ELSEIF (p == 2) THEN
            value = fcubehalo(i-ncube_reconstruct+1,j,3)
          ELSEIF (p == 3) THEN
            value = fcubehalo(i-ncube_reconstruct+1,j,4)
          ELSEIF (p == 4) THEN
            value = fcubehalo(i-ncube_reconstruct+1,j,1)
          ELSEIF (p == 5) THEN
            value = fcubehalo(ncube_reconstruct-j,i-ncube_reconstruct+1,2)
          ELSEIF (p == 6) THEN
            value = fcubehalo(j,2*ncube_reconstruct-i-1,2)
          ENDIF

          !Left halo
        ELSEIF ((j < 1) .AND. (i > 0) .AND. (i < ncube_reconstruct)) THEN
          IF (p == 1) THEN
            value = fcubehalo(i,ncube_reconstruct-1+j,5)
          ELSEIF (p == 2) THEN
            value = fcubehalo(ncube_reconstruct-1+j,ncube_reconstruct-i,5)
          ELSEIF (p == 3) THEN
            value = fcubehalo(ncube_reconstruct-i,1-j,5)
          ELSEIF (p == 4) THEN
            value = fcubehalo(1-j,i,5)
          ELSEIF (p == 5) THEN
            value = fcubehalo(ncube_reconstruct-i,1-j,3)
          ELSEIF (p == 6) THEN
            value = fcubehalo(i,ncube_reconstruct-1+j,1)
          ENDIF

          !Right halo
        ELSEIF ((j > ncube_reconstruct-1) .AND. (i > 0) .AND. (i < ncube_reconstruct)) THEN
          IF (p == 1) THEN
            value = fcubehalo(i,j-ncube_reconstruct+1,6)
          ELSEIF (p == 2) THEN
            value = fcubehalo(2*ncube_reconstruct-j-1,i,6)
          ELSEIF (p == 3) THEN
            value = fcubehalo(ncube_reconstruct-i, 2*ncube_reconstruct-j-1,6)
          ELSEIF (p == 4) THEN
            value = fcubehalo(j-ncube_reconstruct+1,ncube_reconstruct-i,6)
          ELSEIF (p == 5) THEN
            value = fcubehalo(i,j-ncube_reconstruct+1,1)
          ELSEIF (p == 6) THEN
            value = fcubehalo(ncube_reconstruct-i, 2*ncube_reconstruct-j-1,3)
          ENDIF

        ENDIF

        END IF
        min_val = MIN(min_val, value)
        max_val = MAX(max_val, value)
      ENDDO
    ENDDO
  END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE MonotonizeABPGradient
!
! Description:
!   Apply a monotonic filter to the calculated ABP gradient.
!
! Parameters:
!   fcubehalo - Scalar field on the cubed sphere to use in reconstruction
!   order - Order of the reconstruction
!   recons (INOUT) - Array of reconstructed coefficients
!   selective - whether to apply a simple form of selective limiting, 
  !which assumes that if a point is larger/smaller than ALL of its 
  !surrounding points, that the extremum is physical, and that 
  !filtering should not be applied to it.
!
! Remarks:
!   This monotonizing scheme is based on the monotone scheme for unstructured
!   grids of Barth and Jespersen (1989).
!------------------------------------------------------------------------------
  SUBROUTINE MonotonizeABPGradient(fcubehalo, order, recons, selective)

!    USE selective_limiting

    IMPLICIT NONE

    REAL (KIND=dbl_kind), DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6), &
                          INTENT(IN) :: fcubehalo

    INTEGER (KIND=int_kind), INTENT(IN) :: order

    LOGICAL, INTENT(IN) :: selective

    REAL (KIND=dbl_kind), DIMENSION(:,:,:,:), INTENT(INOUT) :: recons

    ! Local variables
    INTEGER (KIND=int_kind) :: i, j, k, m, n, skip

    REAL (KIND=dbl_kind) :: local_min, local_max, value, phi, min_phi
    REAL (KIND=dbl_kind) :: disc, mx, my, lam, gamma_min, gamma_max
    REAL (KIND=dbl_kind), DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6) :: &
         gamma

    ! The first-order piecewise constant scheme is monotone by construction
    IF (order == 1) THEN
      RETURN
    ENDIF

!
! xxxxx
!
!    IF (selective) THEN
!       CALL smoothness2D(fcubehalo, gamma, 2)
!       WRITE(*,*) 'gamma range: max ', MAXVAL(gamma), " min ", MINVAL(gamma)
!       DO i=1,ncube_reconstruct-1
!          WRITE(*,*) gamma(i, i, 3)
!       ENDDO
!       skip = 0
!    END IF
       

    ! Apply monotone limiting
    DO k = 1, 6
      DO j = 1, ncube_reconstruct-1
        DO i = 1, ncube_reconstruct-1


           IF (selective) THEN

              CALL ABPHaloMinMax(gamma, i, j, k, gamma_min, gamma_max, .FALSE.)

              IF (gamma_max/(gamma_min + tiny) < lammax) THEN
                 skip = skip + 1
                 CYCLE
              END IF

           END IF
              
           CALL ABPHaloMinMax(fcubehalo, i, j, k, local_min, local_max,.FALSE.)


          ! Initialize the limiter
          min_phi = one

          ! For the second-order calculation, the minima and maxima will occur
          ! at the corner points of the element
          DO m = i, i+1
            DO n = j, j+1

              ! Evaluate the function at each corner point
              CALL EvaluateABPReconstruction(                                &
                     fcubehalo, recons, i, j, k, gp(m), gp(n), order, value)

              CALL AdjustLimiter(                                            &
                     value, fcubehalo(i,j,k), local_min, local_max, min_phi)
            ENDDO
          ENDDO

          ! For the third order method, the minima and maxima may occur along
          ! the line segments given by du/dx = 0 and du/dy = 0.  Also check
          ! for the presence of a maxima / minima of the quadratic within
          ! the domain.
          IF (order == 3) THEN
            disc = recons(5,i,j,k)**2 - 4.0 * recons(4,i,j,k) * recons(3,i,j,k)

            ! Check if the quadratic is minimized within the element
            IF (ABS(disc) > tiny) THEN
              mx = - recons(5,i,j,k) * recons(2,i,j,k)        &
                   + 2.0 * recons(4,i,j,k) * recons(1,i,j,k)
              my = - recons(5,i,j,k) * recons(1,i,j,k)        &
                   + 2.0 * recons(3,i,j,k) * recons(2,i,j,k)

              mx = mx / disc + abp_centroid(1,i,j)
              my = my / disc + abp_centroid(2,i,j)

              IF ((mx - TAN(gp(i))   > -tiny) .AND. &
                  (mx - TAN(gp(i+1)) <  tiny) .AND. &
                  (my - TAN(gp(j))   > -tiny) .AND. &
                  (my - TAN(gp(j+1)) <  tiny)       &
              ) THEN
                CALL EvaluateABPReconstruction(                         &
                       fcubehalo, recons, i, j, k, ATAN(mx), ATAN(my),  &
                       order, value)

                CALL AdjustLimiter(                         &
                       value, fcubehalo(i,j,k),             &
                       local_min, local_max, min_phi)
              ENDIF
            ENDIF

            ! Check all potential minimizer points along element boundaries
            IF (ABS(recons(5,i,j,k)) > tiny) THEN

              ! Left/right edge, intercept with du/dx = 0
              DO m = i, i+1
                my = - recons(1,i,j,k) - 2.0 * recons(3,i,j,k) * &
                     (TAN(gp(m)) - abp_centroid(1,i,j))

                my = my / recons(5,i,j,k) + abp_centroid(2,i,j)

                IF ((my < TAN(gp(j))) .OR. (my > TAN(gp(j+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, gp(m), ATAN(my),  &
                      order, value)

                CALL AdjustLimiter(                   &
                      value, fcubehalo(i,j,k),        &
                      local_min, local_max, min_phi)
              ENDDO

              ! Top/bottom edge, intercept with du/dy = 0
              DO n = j, j+1
                mx = - recons(2,i,j,k) - 2.0 * recons(4,i,j,k) * &
                     (TAN(gp(n)) - abp_centroid(2,i,j))

                mx = mx / recons(5,i,j,k) + abp_centroid(1,i,j)

                IF ((mx < TAN(gp(i))) .OR. (mx > TAN(gp(i+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, ATAN(mx), gp(n),  &
                      order, value)

                CALL AdjustLimiter(                   &
                      value, fcubehalo(i,j,k),        &
                      local_min, local_max, min_phi)
              ENDDO
            ENDIF

            ! Top/bottom edge, intercept with du/dx = 0
            IF (ABS(recons(3,i,j,k)) > tiny) THEN
              DO n = j, j+1
                mx = - recons(1,i,j,k) - recons(5,i,j,k) * &
                     (TAN(gp(n)) - abp_centroid(2,i,j))

                mx = mx / (2.0 * recons(3,i,j,k)) + abp_centroid(1,i,j)

                IF ((mx < TAN(gp(i))) .OR. (mx > TAN(gp(i+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, ATAN(mx), gp(n),  &
                      order, value)

                CALL AdjustLimiter(                   &
                      value, fcubehalo(i,j,k),        &
                      local_min, local_max, min_phi)
              ENDDO
            ENDIF

            ! Left/right edge, intercept with du/dy = 0
            IF (ABS(recons(4,i,j,k)) > tiny) THEN
              DO m = i, i+1
                my = - recons(2,i,j,k) - recons(5,i,j,k) * &
                     (TAN(gp(m)) - abp_centroid(1,i,j))

                my = my / (2.0 * recons(4,i,j,k)) + abp_centroid(2,i,j)

                IF ((my < TAN(gp(j))) .OR. (my > TAN(gp(j+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, gp(m), ATAN(my),  &
                      order, value)

                CALL AdjustLimiter(                   &
                      value, fcubehalo(i,j,k),        &
                      local_min, local_max, min_phi)
              ENDDO
            ENDIF
          ENDIF

          IF ((min_phi < -tiny) .OR. (min_phi > one + tiny)) THEN
            WRITE (*,*) 'Fatal Error: In MonotonizeABPGradient'
            WRITE (*,*) 'Slope limiter out of range: ', min_phi
            STOP
          ENDIF

          ! Apply monotone limiter to all reconstruction coefficients
          recons(1,i,j,k) = min_phi * recons(1,i,j,k)
          recons(2,i,j,k) = min_phi * recons(2,i,j,k)

          IF (order > 2) THEN
            recons(3,i,j,k) = min_phi * recons(3,i,j,k)
            recons(4,i,j,k) = min_phi * recons(4,i,j,k)
            recons(5,i,j,k) = min_phi * recons(5,i,j,k)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    IF (selective) WRITE(*,*) 'skipped ', skip, ' points out of ', 6*(ncube_reconstruct-1)**2

  END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE PosDefABPGradient
!
! Description:
!   Scale the reconstructions so they are positive definite
!
! Parameters:
!   fcubehalo - Scalar field on the cubed sphere to use in reconstruction
!   order - Order of the reconstruction
!   recons (INOUT) - Array of reconstructed coefficients
!
! Remarks:
!   This monotonizing scheme is based on the monotone scheme for unstructured
!   grids of Barth and Jespersen (1989), but simpler. This simply finds the 
!   minimum and then scales the reconstruction so that it is 0. 
!------------------------------------------------------------------------------
  SUBROUTINE PosDefABPGradient(fcubehalo, order, recons)

    IMPLICIT NONE

    REAL (KIND=dbl_kind), DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6), &
                          INTENT(IN) :: fcubehalo

    INTEGER (KIND=int_kind), INTENT(IN) :: order

    REAL (KIND=dbl_kind), DIMENSION(:,:,:,:), INTENT(INOUT) :: recons

    ! Local variables
    INTEGER (KIND=int_kind) :: i, j, k, m, n

    REAL (KIND=dbl_kind) :: local_min, local_max, value, phi, min_phi
    REAL (KIND=dbl_kind) :: disc, mx, my, lam, gamma_min, gamma_max
    REAL (KIND=dbl_kind), DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6) :: &
         gamma

    ! The first-order piecewise constant scheme is monotone by construction
    IF (order == 1) THEN
      RETURN
    ENDIF


    ! Apply monotone limiting
    DO k = 1, 6
      DO j = 1, ncube_reconstruct-1
        DO i = 1, ncube_reconstruct-1

           !If the average value in the cell is 0.0, then we should skip 
           !all of the scaling and just set the reconstruction to 0.0
!           IF (ABS(fcubehalo(i,j,k)) < tiny) THEN
!              recons(:,i,j,k) = 0.0
!              CYCLE
!           END IF
           CALL ABPHaloMinMax(fcubehalo, i, j, k, local_min, local_max,.FALSE.)

           
           !This allowance for miniscule negative values appearing around the cell being 
           !filtered/limited. Before this, negative values would be caught in adjust_limiter
           !and would stop the model. Doing this only causes minor negative values; no blowing
           !up is observed. The rationale is the same as for the monotone filter, which does
           !allow miniscule negative values due to roundoff error --- of the order E-10 ---
           !in flux-form methods (and E-17 in the s-L method, indicating that roundoff error
           !is more severe in the flux-form method, as we expect since we are often subtracting
           !2.0 values which are very close together. 
           local_min = MIN(0.0,local_min) 
           local_max = bignum !prevents scaling upward; for positive 
                              !definite limiting we don't care about the upper bound

          ! Initialize the limiter
          min_phi = one

          ! For the second-order calculation, the minima and maxima will occur
          ! at the corner points of the element
          DO m = i, i+1
            DO n = j, j+1

              ! Evaluate the function at each corner point
              CALL EvaluateABPReconstruction(                                &
                     fcubehalo, recons, i, j, k, gp(m), gp(n), order, value)

              CALL AdjustLimiter(                                            &
                     value, fcubehalo(i,j,k), local_min, local_max, min_phi)
            ENDDO
          ENDDO

          ! For the third order method, the minima and maxima may occur along
          ! the line segments given by du/dx = 0 and du/dy = 0.  Also check
          ! for the presence of a maxima / minima of the quadratic within
          ! the domain.
          IF (order == 3) THEN
            disc = recons(5,i,j,k)**2 - 4.0 * recons(4,i,j,k) * recons(3,i,j,k)

            ! Check if the quadratic is minimized within the element
            IF (ABS(disc) > tiny) THEN
              mx = - recons(5,i,j,k) * recons(2,i,j,k)        &
                   + 2.0 * recons(4,i,j,k) * recons(1,i,j,k)
              my = - recons(5,i,j,k) * recons(1,i,j,k)        &
                   + 2.0 * recons(3,i,j,k) * recons(2,i,j,k)

              mx = mx / disc + abp_centroid(1,i,j)
              my = my / disc + abp_centroid(2,i,j)

              IF ((mx - TAN(gp(i))   > -tiny) .AND. &
                  (mx - TAN(gp(i+1)) <  tiny) .AND. &
                  (my - TAN(gp(j))   > -tiny) .AND. &
                  (my - TAN(gp(j+1)) <  tiny)       &
              ) THEN
                CALL EvaluateABPReconstruction(                         &
                       fcubehalo, recons, i, j, k, ATAN(mx), ATAN(my),  &
                       order, value)

                CALL AdjustLimiter(                         &
                       value, fcubehalo(i,j,k),             &
                       local_min, local_max, min_phi)
              ENDIF
            ENDIF

            ! Check all potential minimizer points along element boundaries
            IF (ABS(recons(5,i,j,k)) > tiny) THEN

              ! Left/right edge, intercept with du/dx = 0
              DO m = i, i+1
                my = - recons(1,i,j,k) - 2.0 * recons(3,i,j,k) * &
                     (TAN(gp(m)) - abp_centroid(1,i,j))

                my = my / recons(5,i,j,k) + abp_centroid(2,i,j)

                IF ((my < TAN(gp(j))) .OR. (my > TAN(gp(j+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, gp(m), ATAN(my),  &
                      order, value)

                CALL AdjustLimiter(                   &
                      value, fcubehalo(i,j,k),        &
                      local_min, local_max, min_phi)
              ENDDO

              ! Top/bottom edge, intercept with du/dy = 0
              DO n = j, j+1
                mx = - recons(2,i,j,k) - 2.0 * recons(4,i,j,k) * &
                     (TAN(gp(n)) - abp_centroid(2,i,j))

                mx = mx / recons(5,i,j,k) + abp_centroid(1,i,j)

                IF ((mx < TAN(gp(i))) .OR. (mx > TAN(gp(i+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, ATAN(mx), gp(n),  &
                      order, value)

                CALL AdjustLimiter(                   &
                      value, fcubehalo(i,j,k),        &
                      local_min, local_max, min_phi)
              ENDDO
            ENDIF

            ! Top/bottom edge, intercept with du/dx = 0
            IF (ABS(recons(3,i,j,k)) > tiny) THEN
              DO n = j, j+1
                mx = - recons(1,i,j,k) - recons(5,i,j,k) * &
                     (TAN(gp(n)) - abp_centroid(2,i,j))

                mx = mx / (2.0 * recons(3,i,j,k)) + abp_centroid(1,i,j)

                IF ((mx < TAN(gp(i))) .OR. (mx > TAN(gp(i+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, ATAN(mx), gp(n),  &
                      order, value)

                CALL AdjustLimiter(                   &
                      value, fcubehalo(i,j,k),        &
                      local_min, local_max, min_phi)
              ENDDO
            ENDIF

            ! Left/right edge, intercept with du/dy = 0
            IF (ABS(recons(4,i,j,k)) > tiny) THEN
              DO m = i, i+1
                my = - recons(2,i,j,k) - recons(5,i,j,k) * &
                     (TAN(gp(m)) - abp_centroid(1,i,j))

                my = my / (2.0 * recons(4,i,j,k)) + abp_centroid(2,i,j)

                IF ((my < TAN(gp(j))) .OR. (my > TAN(gp(j+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, gp(m), ATAN(my),  &
                      order, value)

                CALL AdjustLimiter(                   &
                      value, fcubehalo(i,j,k),        &
                      local_min, local_max, min_phi)
              ENDDO
            ENDIF
          ENDIF

          IF ((min_phi < -tiny) .OR. (min_phi > one + tiny)) THEN
            WRITE (*,*) 'Fatal Error: In MonotonizeABPGradient'
            WRITE (*,*) 'Slope limiter out of range: ', min_phi
            STOP
          ENDIF

          ! Apply monotone limiter to all reconstruction coefficients
          recons(1,i,j,k) = min_phi * recons(1,i,j,k)
          recons(2,i,j,k) = min_phi * recons(2,i,j,k)

          IF (order > 2) THEN
            recons(3,i,j,k) = min_phi * recons(3,i,j,k)
            recons(4,i,j,k) = min_phi * recons(4,i,j,k)
            recons(5,i,j,k) = min_phi * recons(5,i,j,k)
          ENDIF

        ENDDO
      ENDDO
    ENDDO


  END SUBROUTINE PosDefABPGradient

!------------------------------------------------------------------------------
! SUBROUTINE MonotonizeABPGradient_New
!
! Description:
!   Apply a monotonic filter to the calculated ABP gradient.
!
! Parameters:
!   fcubehalo - Scalar field on the cubed sphere to use in reconstruction
!   order - Order of the reconstruction
!   recons (INOUT) - Array of reconstructed coefficients
!
! Remarks:
!   This monotonizing scheme is similar to the one in MonotonizeABPGradient,
!   except the second order derivatives are limited after the first order
!   derivatives.
!------------------------------------------------------------------------------
  SUBROUTINE MonotonizeABPGradient_New(fcubehalo, order, recons)

    IMPLICIT NONE

    REAL (KIND=dbl_kind), DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6), &
                          INTENT(IN) :: fcubehalo

    INTEGER (KIND=int_kind), INTENT(IN) :: order

    REAL (KIND=dbl_kind), DIMENSION(:,:,:,:), INTENT(INOUT) :: recons

    ! Local variables
    INTEGER (KIND=int_kind) :: i, j, k, m, n

    REAL (KIND=dbl_kind) :: local_min, local_max, value, phi, min_phi, linval
    REAL (KIND=dbl_kind) :: disc, mx, my

    ! The first-order piecewise constant scheme is monotone by construction
    IF (order == 1) THEN
      RETURN
    ENDIF

    ! Apply monotone limiting
    DO k = 1, 6
      DO j = 1, ncube_reconstruct-1
        DO i = 1, ncube_reconstruct-1
          CALL ABPHaloMinMax(fcubehalo, i, j, k, local_min, local_max, .FALSE.)

          ! Initialize the limiter
          min_phi = one

          ! For the second-order calculation, the minima and maxima will occur
          ! at the corner points of the element
          DO m = i, i+1
            DO n = j, j+1

              ! Evaluate the function at each corner point, only taking into
              ! account the linear component of the reconstruction.
              value =                                                        &
                fcubehalo(i,j,k)                                             &
                + recons(1,i,j,k) * (TAN(gp(m)) - abp_centroid(1,i,j))       &
                + recons(2,i,j,k) * (TAN(gp(n)) - abp_centroid(2,i,j))

              CALL AdjustLimiter(                                            &
                     value, fcubehalo(i,j,k), local_min, local_max, min_phi)
            ENDDO
          ENDDO

          ! Apply monotone limiter to all reconstruction coefficients
          IF ((min_phi < -tiny) .OR. (min_phi > one + tiny)) THEN
            WRITE (*,*) 'Fatal Error: In MonotonizeABPGradient'
            WRITE (*,*) 'Slope limiter out of range: ', min_phi
            STOP
          ENDIF

          recons(1,i,j,k) = min_phi * recons(1,i,j,k)
          recons(2,i,j,k) = min_phi * recons(2,i,j,k)

          ! For the third order method, the minima and maxima may occur along
          ! the line segments given by du/dx = 0 and du/dy = 0.  Also check
          ! for the presence of a maxima / minima of the quadratic within
          ! the domain.
          IF (order == 3) THEN
            ! Reset the limiter
            min_phi = one

            ! Calculate discriminant, which we use to determine the absolute
            ! minima/maxima of the paraboloid
            disc = recons(5,i,j,k)**2 - 4.0 * recons(4,i,j,k) * recons(3,i,j,k)

            ! Check if the quadratic is minimized within the element
            IF (ABS(disc) > tiny) THEN
              mx = - recons(5,i,j,k) * recons(2,i,j,k)        &
                   + 2.0 * recons(4,i,j,k) * recons(1,i,j,k)
              my = - recons(5,i,j,k) * recons(1,i,j,k)        &
                   + 2.0 * recons(3,i,j,k) * recons(2,i,j,k)

              mx = mx / disc + abp_centroid(1,i,j)
              my = my / disc + abp_centroid(2,i,j)

              IF ((mx - TAN(gp(i))   > -tiny) .AND. &
                  (mx - TAN(gp(i+1)) <  tiny) .AND. &
                  (my - TAN(gp(j))   > -tiny) .AND. &
                  (my - TAN(gp(j+1)) <  tiny)       &
              ) THEN
                CALL EvaluateABPReconstruction(                         &
                       fcubehalo, recons, i, j, k, ATAN(mx), ATAN(my),  &
                       order, value)

                linval =                                             &
                  fcubehalo(i,j,k)                                   &
                  + recons(1,i,j,k) * (mx - abp_centroid(1,i,j))     &
                  + recons(2,i,j,k) * (my - abp_centroid(2,i,j))

                IF (linval < local_min) THEN
                  linval = local_min
                ENDIF
                IF (linval > local_max) THEN
                  linval = local_max
                ENDIF

                CALL AdjustLimiter(                                  &
                       value, linval, local_min, local_max, min_phi)
              ENDIF
            ENDIF

            ! Check all potential minimizer points along element boundaries
            IF (ABS(recons(5,i,j,k)) > tiny) THEN

              ! Left/right edge, intercept with du/dx = 0
              DO m = i, i+1
                my = - recons(1,i,j,k) - 2.0 * recons(3,i,j,k) * &
                     (TAN(gp(m)) - abp_centroid(1,i,j))

                my = my / recons(5,i,j,k) + abp_centroid(2,i,j)

                IF ((my < TAN(gp(j))) .OR. (my > TAN(gp(j+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, gp(m), ATAN(my),  &
                      order, value)

                linval =                                                     &
                  fcubehalo(i,j,k)                                           &
                  + recons(1,i,j,k) * (TAN(gp(m)) - abp_centroid(1,i,j))     &
                  + recons(2,i,j,k) * (my - abp_centroid(2,i,j))

                IF (linval < local_min) THEN
                  linval = local_min
                ENDIF
                IF (linval > local_max) THEN
                  linval = local_max
                ENDIF

                CALL AdjustLimiter(                                 &
                      value, linval, local_min, local_max, min_phi)
              ENDDO

              ! Top/bottom edge, intercept with du/dy = 0
              DO n = j, j+1
                mx = - recons(2,i,j,k) - 2.0 * recons(4,i,j,k) * &
                     (TAN(gp(n)) - abp_centroid(2,i,j))

                mx = mx / recons(5,i,j,k) + abp_centroid(1,i,j)

                IF ((mx < TAN(gp(i))) .OR. (mx > TAN(gp(i+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, ATAN(mx), gp(n),  &
                      order, value)

                linval =                                                     &
                  fcubehalo(i,j,k)                                           &
                  + recons(1,i,j,k) * (mx - abp_centroid(1,i,j))             &
                  + recons(2,i,j,k) * (TAN(gp(n)) - abp_centroid(2,i,j))

                IF (linval < local_min) THEN
                  linval = local_min
                ENDIF
                IF (linval > local_max) THEN
                  linval = local_max
                ENDIF

                CALL AdjustLimiter(                   &
                      value, linval, local_min, local_max, min_phi)
              ENDDO
            ENDIF

            ! Top/bottom edge, intercept with du/dx = 0
            IF (ABS(recons(3,i,j,k)) > tiny) THEN
              DO n = j, j+1
                mx = - recons(1,i,j,k) - recons(5,i,j,k) * &
                     (TAN(gp(n)) - abp_centroid(2,i,j))

                mx = mx / (2.0 * recons(3,i,j,k)) + abp_centroid(1,i,j)

                IF ((mx < TAN(gp(i))) .OR. (mx > TAN(gp(i+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, ATAN(mx), gp(n),  &
                      order, value)

                linval =                                                     &
                  fcubehalo(i,j,k)                                           &
                  + recons(1,i,j,k) * (mx - abp_centroid(1,i,j))             &
                  + recons(2,i,j,k) * (TAN(gp(n)) - abp_centroid(2,i,j))

                IF (linval < local_min) THEN
                  linval = local_min
                ENDIF
                IF (linval > local_max) THEN
                  linval = local_max
                ENDIF

                CALL AdjustLimiter(                   &
                      value, linval, local_min, local_max, min_phi)
              ENDDO
            ENDIF

            ! Left/right edge, intercept with du/dy = 0
            IF (ABS(recons(4,i,j,k)) > tiny) THEN
              DO m = i, i+1
                my = - recons(2,i,j,k) - recons(5,i,j,k) * &
                     (TAN(gp(m)) - abp_centroid(1,i,j))

                my = my / (2.0 * recons(4,i,j,k)) + abp_centroid(2,i,j)

                IF ((my < TAN(gp(j))) .OR. (my > TAN(gp(j+1)))) THEN
                  CYCLE
                ENDIF

                CALL EvaluateABPReconstruction(                     &
                      fcubehalo, recons, i, j, k, gp(m), ATAN(my),  &
                      order, value)

                linval =                                                     &
                  fcubehalo(i,j,k)                                           &
                  + recons(1,i,j,k) * (TAN(gp(m)) - abp_centroid(1,i,j))     &
                  + recons(2,i,j,k) * (my - abp_centroid(2,i,j))

                IF (linval < local_min) THEN
                  linval = local_min
                ENDIF
                IF (linval > local_max) THEN
                  linval = local_max
                ENDIF

                CALL AdjustLimiter(                   &
                      value, linval, local_min, local_max, min_phi)
              ENDDO
            ENDIF

            ! For the second-order calculation, the minima and maxima will occur
            ! at the corner points of the element
            DO m = i, i+1
              DO n = j, j+1

                ! Evaluate the function at each corner point
                CALL EvaluateABPReconstruction(                              &
                       fcubehalo, recons, i, j, k, gp(m), gp(n),             &
                       order, value)

                linval =                                                   &
                  fcubehalo(i,j,k)                                         &
                  + recons(1,i,j,k) * (TAN(gp(m)) - abp_centroid(1,i,j))   &
                  + recons(2,i,j,k) * (TAN(gp(n)) - abp_centroid(2,i,j))

                IF (linval < local_min) THEN
                  linval = local_min
                ENDIF
                IF (linval > local_max) THEN
                  linval = local_max
                ENDIF

                CALL AdjustLimiter(                                          &
                       value, linval, local_min, local_max, min_phi)
              ENDDO
            ENDDO

            IF ((min_phi < -tiny) .OR. (min_phi > one + tiny)) THEN
              WRITE (*,*) 'Fatal Error: In MonotonizeABPGradient'
              WRITE (*,*) 'Slope limiter out of range: ', min_phi
              STOP
            ENDIF

            WRITE (*,*) '2: ', min_phi

            recons(1,i,j,k) = min_phi * recons(1,i,j,k)
            recons(2,i,j,k) = min_phi * recons(2,i,j,k)
            recons(3,i,j,k) = min_phi * recons(3,i,j,k)
            recons(4,i,j,k) = min_phi * recons(4,i,j,k)
            recons(5,i,j,k) = min_phi * recons(5,i,j,k)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE ReconstructABPGradient_NEL
!
! Description:
!   Construct a non-equidistant linear reconstruction of the gradient
!   within each element on an ABP grid.
!
! Parameters:
!   fcubehalo - Scalar field on the ABP grid to use in reconstruction
!   recons (OUT) - Array of reconstructed coefficients for total elements
!   order - Order of the scheme (2 or 3)
!------------------------------------------------------------------------------
  SUBROUTINE ReconstructABPGradient_NEL(fcubehalo, recons, order)

!    USE CubedSphereTrans
!    USE InterpolateCSLL_Utils

    IMPLICIT NONE

    REAL (KIND=dbl_kind),   &
          DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6), INTENT(IN) :: fcubehalo

    REAL (KIND=dbl_kind), DIMENSION(:,:,:,:), INTENT(OUT) :: recons

    INTEGER (KIND=int_kind), INTENT(IN) :: order

    ! Local variables
    INTEGER (KIND=int_kind) :: i, j, p

    REAL (KIND=dbl_kind) :: alpha1, alpha2, beta1, beta2
    REAL (KIND=dbl_kind) :: dx_left, dx_right, top_value, bot_value

    DO p = 1, 6
      DO j = 1, ncube_reconstruct-1
        DO i = 1, ncube_reconstruct-1
          dx_left = abp_centroid(1,i-1,j) - abp_centroid(1,i,j)
          dx_right = abp_centroid(1,i+1,j) - abp_centroid(1,i,j)

          recons(1,i,j,p) =                                     &
            (+ fcubehalo(i-1,j,p) * dx_right**2                 &
             - fcubehalo(i+1,j,p) * dx_left**2                  &
             - fcubehalo(i,j,p) * (dx_right**2 - dx_left**2)) / &
            (dx_right * dx_left * (dx_right - dx_left))

          dx_left = abp_centroid(2,i,j-1) - abp_centroid(2,i,j)
          dx_right = abp_centroid(2,i,j+1) - abp_centroid(2,i,j)

          recons(2,i,j,p) =                                     &
            (+ fcubehalo(i,j-1,p) * dx_right**2                 &
             - fcubehalo(i,j+1,p) * dx_left**2                  &
             - fcubehalo(i,j,p) * (dx_right**2 - dx_left**2)) / &
            (dx_right * dx_left * (dx_right - dx_left))

          IF (order > 2) THEN
            dx_left = abp_centroid(1,i-1,j) - abp_centroid(1,i,j)
            dx_right = abp_centroid(1,i+1,j) - abp_centroid(1,i,j)

            recons(3,i,j,p) =                               &
              (+ fcubehalo(i-1,j,p) * dx_right              &
               - fcubehalo(i+1,j,p) * dx_left               &
               - fcubehalo(i,j,p) * (dx_right - dx_left)) / &
              (dx_right * dx_left * (dx_left - dx_right))

            dx_left = abp_centroid(2,i,j-1) - abp_centroid(2,i,j)
            dx_right = abp_centroid(2,i,j+1) - abp_centroid(2,i,j)

            recons(4,i,j,p) =                               &
              (+ fcubehalo(i,j-1,p) * dx_right              &
               - fcubehalo(i,j+1,p) * dx_left               &
               - fcubehalo(i,j,p) * (dx_right - dx_left)) / &
              (dx_right * dx_left * (dx_left - dx_right))
          ENDIF
        ENDDO
      ENDDO

      IF (order > 2) THEN
        DO j = 1, ncube_reconstruct-1
          DO i = 1, ncube_reconstruct-1
            dx_left = abp_centroid(1,i-1,j+1) - abp_centroid(1,i,j+1)
            dx_right = abp_centroid(1,i+1,j+1) - abp_centroid(1,i,j+1)
  
            top_value =                                             &
              (+ fcubehalo(i-1,j+1,p) * dx_right**2                 &
               - fcubehalo(i+1,j+1,p) * dx_left**2                  &
               - fcubehalo(i,j+1,p) * (dx_right**2 - dx_left**2)) / &
              (dx_right * dx_left * (dx_right - dx_left))
  
            dx_left = abp_centroid(1,i-1,j-1) - abp_centroid(1,i,j-1)
            dx_right = abp_centroid(1,i+1,j-1) - abp_centroid(1,i,j-1)
  
            bot_value =                                             &
              (+ fcubehalo(i-1,j-1,p) * dx_right**2                 &
               - fcubehalo(i+1,j-1,p) * dx_left**2                  &
               - fcubehalo(i,j-1,p) * (dx_right**2 - dx_left**2)) / &
              (dx_right * dx_left * (dx_right - dx_left))
  
            dx_left = abp_centroid(2,i,j-1) - abp_centroid(2,i,j)
            dx_right = abp_centroid(2,i,j+1) - abp_centroid(2,i,j)
  
            recons(5,i,j,p) =                                    &
              (+ bot_value * dx_right**2                         &
               - top_value * dx_left**2                          &
               - recons(1,i,j,p) * (dx_right**2 - dx_left**2)) / &
              (dx_right * dx_left * (dx_right - dx_left))

          ENDDO
        ENDDO
      ENDIF
    ENDDO

  END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE ReconstructABPGradient_NEP
!
! Description:
!   Construct a non-equidistant parabolic reconstruction of the gradient
!   within each element on an ABP grid.
!
! Parameters:
!   fcubehalo - Scalar field on the ABP grid to use in reconstruction
!   recons (OUT) - Array of reconstructed coefficients for total elements
!   order - Order of the scheme (2 or 3)
!------------------------------------------------------------------------------
  SUBROUTINE ReconstructABPGradient_NEP(fcubehalo, recons, order)


!    USE CubedSphereTrans
!    USE InterpolateCSLL_Utils

    IMPLICIT NONE

    REAL (KIND=dbl_kind),   &
          DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6), INTENT(IN) :: fcubehalo

    REAL (KIND=dbl_kind), DIMENSION(:,:,:,:), INTENT(OUT) :: recons

    INTEGER (KIND=int_kind), INTENT(IN) :: order

    ! Local variables
    INTEGER (KIND=int_kind) :: i, j, p

    REAL (KIND=dbl_kind) :: x1, x2, x4, x5, y1, y2, y3, y4, y5

    REAL (KIND=dbl_kind), DIMENSION(5) :: t, pa, denom

    DO p = 1, 6
      DO j = 1, ncube_reconstruct-1
        DO i = 1, ncube_reconstruct-1
          ! X-direction reconstruction
          x1 = abp_centroid(1,i-2,j) - abp_centroid(1,i,j)
          x2 = abp_centroid(1,i-1,j) - abp_centroid(1,i,j)
          x4 = abp_centroid(1,i+1,j) - abp_centroid(1,i,j)
          x5 = abp_centroid(1,i+2,j) - abp_centroid(1,i,j)

          !IF (i == 1) THEN
          !  x1 = piq
          !ELSEIF (i == ncube_reconstruct-1) THEN
          !  x5 = -piq
          !ENDIF

          y1 = fcubehalo(i-2,j,p)
          y2 = fcubehalo(i-1,j,p)
          y3 = fcubehalo(i,j,p)
          y4 = fcubehalo(i+1,j,p)
          y5 = fcubehalo(i+2,j,p)

          denom(1) = (x2 - x1) * (x4 - x1) * (x5 - x1) * x1
          denom(2) = (x1 - x2) * (x4 - x2) * (x5 - x2) * x2
          denom(4) = (x1 - x4) * (x2 - x4) * (x5 - x4) * x4
          denom(5) = (x1 - x5) * (x2 - x5) * (x4 - x5) * x5

          t(1) = x5 * x4 * x2
          t(2) = x5 * x4 * x1
          t(4) = x5 * x2 * x1
          t(5) = x4 * x2 * x1
          t(3) = (t(1) + t(2) + t(4) + t(5)) / (x1 * x2 * x4 * x5) 

          pa(1) = x2 * x4 + x2 * x5 + x4 * x5
          pa(2) = x1 * x4 + x1 * x5 + x4 * x5
          pa(4) = x1 * x2 + x1 * x5 + x2 * x5
          pa(5) = x1 * x2 + x1 * x4 + x2 * x4
          pa(3) = (pa(1) + pa(2) + pa(4) + pa(5)) / (2.0 * x1 * x2 * x4 * x5)

          recons(1,i,j,p) =                                           &
            + y1 * t(1) / denom(1)                                    &
            + y2 * t(2) / denom(2)                                    &
            - y3 * t(3)                                               &
            + y4 * t(4) / denom(4)                                    &
            + y5 * t(5) / denom(5)

          IF (order > 2) THEN
            recons(3,i,j,p) =                                         &
              - y1 * pa(1) / denom(1)                                 &
              - y2 * pa(2) / denom(2)                                 &
              + y3 * pa(3)                                            &
              - y4 * pa(4) / denom(4)                                 &
              - y5 * pa(5) / denom(5)
          ENDIF

          ! Y-direction reconstruction
          x1 = abp_centroid(2,i,j-2) - abp_centroid(2,i,j)
          x2 = abp_centroid(2,i,j-1) - abp_centroid(2,i,j)
          x4 = abp_centroid(2,i,j+1) - abp_centroid(2,i,j)
          x5 = abp_centroid(2,i,j+2) - abp_centroid(2,i,j)

          !IF (j == 1) THEN
          !  x1 = piq
          !ELSEIF (j == ncube_reconstruct-1) THEN
          !  x5 = -piq
          !ENDIF

          y1 = fcubehalo(i,j-2,p)
          y2 = fcubehalo(i,j-1,p)
          y3 = fcubehalo(i,j,p)
          y4 = fcubehalo(i,j+1,p)
          y5 = fcubehalo(i,j+2,p)

          denom(1) = (x2 - x1) * (x4 - x1) * (x5 - x1) * x1
          denom(2) = (x1 - x2) * (x4 - x2) * (x5 - x2) * x2
          denom(4) = (x1 - x4) * (x2 - x4) * (x5 - x4) * x4
          denom(5) = (x1 - x5) * (x2 - x5) * (x4 - x5) * x5

          t(1) = x5 * x4 * x2
          t(2) = x5 * x4 * x1
          t(4) = x5 * x2 * x1
          t(5) = x4 * x2 * x1
          t(3) = (t(1) + t(2) + t(4) + t(5)) / (x1 * x2 * x4 * x5) 

          pa(1) = x2 * x4 + x2 * x5 + x4 * x5
          pa(2) = x1 * x4 + x1 * x5 + x4 * x5
          pa(4) = x1 * x2 + x1 * x5 + x2 * x5
          pa(5) = x1 * x2 + x1 * x4 + x2 * x4
          pa(3) = (pa(1) + pa(2) + pa(4) + pa(5)) / (2.0 * x1 * x2 * x4 * x5)

          recons(2,i,j,p) =                                           &
            + y1 * t(1) / denom(1)                                    &
            + y2 * t(2) / denom(2)                                    &
            - y3 * t(3)                                               &
            + y4 * t(4) / denom(4)                                    &
            + y5 * t(5) / denom(5)

          IF (order > 2) THEN
            recons(4,i,j,p) =                                         &
              - y1 * pa(1) / denom(1)                                 &
              - y2 * pa(2) / denom(2)                                 &
              + y3 * pa(3)                                            &
              - y4 * pa(4) / denom(4)                                 &
              - y5 * pa(5) / denom(5)
            recons(5,i,j,p) = 0.0
          ENDIF

        ENDDO
      ENDDO
      IF (order > 2) THEN
        DO j = 1, ncube_reconstruct-1
          DO i = 1, ncube_reconstruct-1
            x1 = abp_centroid(1,i-1,j+1) - abp_centroid(1,i,j+1)
            x2 = abp_centroid(1,i+1,j+1) - abp_centroid(1,i,j+1)
  
            y2 = (+ fcubehalo(i-1,j+1,p) * x2**2            &
                  - fcubehalo(i+1,j+1,p) * x1**2            &
                  - fcubehalo(i,j+1,p) * (x2**2 - x1**2)) / &
                 (x2 * x1 * (x2 - x1))

            x1 = abp_centroid(1,i-1,j-1) - abp_centroid(1,i,j-1)
            x2 = abp_centroid(1,i+1,j-1) - abp_centroid(1,i,j-1)

            y1 = (+ fcubehalo(i-1,j-1,p) * x2**2            &
                  - fcubehalo(i+1,j-1,p) * x1**2            &
                  - fcubehalo(i,j-1,p) * (x2**2 - x1**2)) / &
                 (x2 * x1 * (x2 - x1))
  
            x1 = abp_centroid(2,i,j-1) - abp_centroid(2,i,j)
            x2 = abp_centroid(2,i,j+1) - abp_centroid(2,i,j)
  
            recons(5,i,j,p) =                         &
              (+ y1 * x2**2                           &
               - y2 * x1**2                           &
               - recons(1,i,j,p) * (x2**2 - x1**2)) / &
              (x2 * x1 * (x2 - x1))

          ENDDO
        ENDDO
      ENDIF
    ENDDO

  END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE ReconstructABPGradient_PLM
!
! Description:
!   Construct a piecewise linear reconstruction of the gradient within
!   each element on an ABP grid.
!
! Parameters:
!   fcubehalo - Scalar field on the ABP grid to use in reconstruction
!   recons (OUT) - Array of reconstructed coefficients for total elements
!   order - Order of the scheme (2 or 3)
!------------------------------------------------------------------------------
  SUBROUTINE ReconstructABPGradient_PLM(fcubehalo, recons, order)

!    USE CubedSphereTrans
!    USE InterpolateCSLL_Utils

    IMPLICIT NONE

    REAL (KIND=dbl_kind),   &
          DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6), INTENT(IN) :: fcubehalo

    REAL (KIND=dbl_kind), DIMENSION(:,:,:,:), INTENT(OUT) :: recons

    INTEGER (KIND=int_kind), INTENT(IN) :: order

    ! Local variables
    INTEGER (KIND=int_kind) :: i, j, p

    REAL (KIND=dbl_kind) :: width

    ! ABP width between elements
    width = pih / DBLE(ncube_reconstruct-1)

    DO p = 1, 6
      DO j = 1, ncube_reconstruct-1
        DO i = 1, ncube_reconstruct-1
          ! df/dx
          recons(1,i,j,p) = (fcubehalo(i+1,j,p) - fcubehalo(i-1,j,p)) / &
                            (2.0 * width)

          ! df/dy
          recons(2,i,j,p) = (fcubehalo(i,j+1,p) - fcubehalo(i,j-1,p)) / &
                            (2.0 * width)

          ! Stretching
          recons(1,i,j,p) = recons(1,i,j,p) / (one + abp_centroid(1,i,j)**2)
          recons(2,i,j,p) = recons(2,i,j,p) / (one + abp_centroid(2,i,j)**2)

          ! Third order scheme
          IF (order > 2) THEN
            ! d^2f/dx^2
            recons(3,i,j,p) =                              &
              (fcubehalo(i+1,j,p) - 2.0 * fcubehalo(i,j,p) &
               + fcubehalo(i-1,j,p)) / (width * width)

            ! d^2f/dy^2
            recons(4,i,j,p) =                              &
              (fcubehalo(i,j+1,p) - 2.0 * fcubehalo(i,j,p) &
               + fcubehalo(i,j-1,p)) / (width * width)

            ! d^2f/dxdy
            recons(5,i,j,p) =                                 &
              (+ fcubehalo(i+1,j+1,p) - fcubehalo(i-1,j+1,p)  &
               - fcubehalo(i+1,j-1,p) + fcubehalo(i-1,j-1,p)  &
              ) / (4.0 * width * width)

            ! Stretching
            recons(3,i,j,p) =                                                 &
              (- 2.0 * abp_centroid(1,i,j) * (one + abp_centroid(1,i,j)**2) * recons(1,i,j,p)                  &
               + recons(3,i,j,p)) / (one + abp_centroid(1,i,j)**2)**2

            recons(4,i,j,p) =                                                 &
              (- 2.0 * abp_centroid(2,i,j) * (one + abp_centroid(2,i,j)**2) * recons(2,i,j,p)                  &
               + recons(4,i,j,p)) / (one + abp_centroid(2,i,j)**2)**2

            recons(5,i,j,p) = recons(5,i,j,p) /                               &
              ((one + abp_centroid(1,i,j)**2) * (one + abp_centroid(2,i,j)**2))

            ! Scaling
            recons(3,i,j,p) = 0.5 * recons(3,i,j,p)
            recons(4,i,j,p) = 0.5 * recons(4,i,j,p)

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE ReconstructABPGradient_PPM
!
! Description:
!   Construct a piecewise parabolic reconstruction of the gradient within
!   each element on an ABP grid.
!
! Parameters:
!   fcubehalo - Scalar field on the ABP grid to use in reconstruction
!   recons (OUT) - Array of reconstructed coefficients for total elements
!   order - Order of the scheme (2 or 3)
!------------------------------------------------------------------------------
  SUBROUTINE ReconstructABPGradient_PPM(fcubehalo, recons, order)


!    USE CubedSphereTrans
!    USE InterpolateCSLL_Utils

    IMPLICIT NONE

    REAL (KIND=dbl_kind),   &
          DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6), INTENT(IN) :: fcubehalo

    REAL (KIND=dbl_kind), DIMENSION(:,:,:,:), INTENT(OUT) :: recons

    INTEGER (KIND=int_kind), INTENT(IN) :: order

    ! Local variables
    INTEGER (KIND=int_kind) :: i, j, p

    REAL (KIND=dbl_kind) :: width

    ! ABP width between elements
    width = pih / DBLE(ncube_reconstruct-1)

    DO p = 1, 6
      DO j = 1, ncube_reconstruct-1
        DO i = 1, ncube_reconstruct-1
          ! df/dalfa
          recons(1,i,j,p) =                                       &
            (+ fcubehalo(i+2,j,p) - 8.0 * fcubehalo(i+1,j,p)    &
             + 8.0 * fcubehalo(i-1,j,p) - fcubehalo(i-2,j,p)) / &
            (- 12.0 * width)

          ! df/dbeta
          recons(2,i,j,p) =                                       &
            (+ fcubehalo(i,j+2,p) - 8.0 * fcubehalo(i,j+1,p)    &
             + 8.0 * fcubehalo(i,j-1,p) - fcubehalo(i,j-2,p)) / &
            (- 12.0 * width)

          ! Stretching
          recons(1,i,j,p) = recons(1,i,j,p) / (one + abp_centroid(1,i,j)**2)
          recons(2,i,j,p) = recons(2,i,j,p) / (one + abp_centroid(2,i,j)**2)

          ! Third order scheme
          IF (order > 2) THEN
            ! d^2f/dx^2
            recons(3,i,j,p) = (- fcubehalo(i+2,j,p)                &
                               + 16_dbl_kind * fcubehalo(i+1,j,p)  &
                               - 30_dbl_kind * fcubehalo(i,j,p)    &
                               + 16_dbl_kind * fcubehalo(i-1,j,p)  &
                               - fcubehalo(i-2,j,p)                &
                              ) / (12_dbl_kind * width**2)

            ! d^2f/dy^2
            recons(4,i,j,p) = (- fcubehalo(i,j+2,p)                &
                               + 16_dbl_kind * fcubehalo(i,j+1,p)  &
                               - 30_dbl_kind * fcubehalo(i,j,p)    &
                               + 16_dbl_kind * fcubehalo(i,j-1,p)  &
                               - fcubehalo(i,j-2,p)                &
                              ) / (12_dbl_kind * width**2)

            ! d^2f/dxdy
            recons(5,i,j,p) =                                 &
              (+ fcubehalo(i+1,j+1,p) - fcubehalo(i-1,j+1,p)  &
               - fcubehalo(i+1,j-1,p) + fcubehalo(i-1,j-1,p)  &
              ) / (4.0 * width * width)

            ! Stretching
            recons(3,i,j,p) =                                                 &
              (- 2.0 * abp_centroid(1,i,j) * (one + abp_centroid(1,i,j)**2) * recons(1,i,j,p)                  &
               + recons(3,i,j,p)) / (one + abp_centroid(1,i,j)**2)**2

            recons(4,i,j,p) =                                                 &
              (- 2.0 * abp_centroid(2,i,j) * (one + abp_centroid(2,i,j)**2) * recons(2,i,j,p)                  &
               + recons(4,i,j,p)) / (one + abp_centroid(2,i,j)**2)**2

            recons(5,i,j,p) = recons(5,i,j,p) /                               &
              ((one + abp_centroid(1,i,j)**2) * (one + abp_centroid(2,i,j)**2))

            ! Scaling
            recons(3,i,j,p) = 0.5 * recons(3,i,j,p)
            recons(4,i,j,p) = 0.5 * recons(4,i,j,p)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE ReconstructABPGradient
!
! Description:
!   Compute the reconstructed gradient in gnomonic coordinates for each
!   ABP element.
!
! Parameters:
!   fcube - Scalar field on the cubed sphere to use in reconstruction
!   halomethod - Method for computing halo elements
!                (0) Piecewise constant
!                (1) Piecewise linear
!                (3) Piecewise cubic
!   recons_method - Method for computing the sub-grid scale gradient
!                   (0) Non-equidistant linear reconstruction
!                   (1) Non-equidistant parabolic reconstruction
!                   (2) Piecewise linear reconstruction with stretching
!                   (3) Piecewise parabolic reconstruction with stretching
!   order - Order of the method being applied
!   kmono - Apply monotone limiting (1) or not (0)
!   recons (INOUT) - Array of reconstructed coefficients
!------------------------------------------------------------------------------
  SUBROUTINE ReconstructABPGradient(                                   &
               fcube, halomethod, recons_method, order, kmono, recons, kpd, kscheme)

!    USE InterpolateCSLL_Utils

    IMPLICIT NONE

    REAL (KIND=dbl_kind), &
            DIMENSION(1:ncube_reconstruct-1, 1:ncube_reconstruct-1, 6), INTENT(IN) :: fcube

    INTEGER (KIND=int_kind), INTENT(IN) :: halomethod, recons_method
    INTEGER (KIND=int_kind), INTENT(IN) :: order, kmono, kpd, kscheme

    REAL (KIND=dbl_kind), DIMENSION(:,:,:,:), INTENT(INOUT) :: recons

    ! Local variables
    INTEGER (KIND=int_kind) :: i, j, p

    REAL (KIND=dbl_kind), DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6) :: fcubehalo

    ! Report status
    WRITE (*,*) '...Performing sub-grid scale reconstruction on ABP grid'

    ! Compute element haloes
    WRITE(*,*) "fill cubed-sphere halo for reconstruction"
    DO p = 1, 6
    IF (halomethod == 0) THEN
       CALL CubedSphereFillHalo(fcube, fcubehalo, p, ncube_reconstruct, 2)

      ELSEIF (halomethod == 1) THEN
        CALL CubedSphereFillHalo_Linear(fcube, fcubehalo, p, ncube_reconstruct)

      ELSEIF (halomethod == 3) THEN
        !halomethod is always 3 in the standard CSLAM setup
        CALL CubedSphereFillHalo_Cubic(fcube, fcubehalo, p, ncube_reconstruct)
      ELSE
        WRITE (*,*) 'Fatal Error: In ReconstructABPGradient'
        WRITE (*,*) 'Invalid halo method: ', halomethod
        WRITE (*,*) 'Halo method must be 0, 1 or 3.'
        STOP
      ENDIF
    ENDDO

    ! Nonequidistant linear reconstruction
    IF (recons_method == 1) THEN
      CALL ReconstructABPGradient_NEL(fcubehalo, recons, order)

    ! Nonequidistant parabolic reconstruction (JCP paper)
    ELSEIF (recons_method == 2) THEN
      WRITE(*,*) "Nonequidistant parabolic reconstruction"
      CALL ReconstructABPGradient_NEP(fcubehalo, recons, order)

    ! Piecewise linear reconstruction with rotation
    ELSEIF (recons_method == 3) THEN
      CALL ReconstructABPGradient_PLM(fcubehalo, recons, order)

    ! Piecewise parabolic reconstruction with rotation
    ELSEIF (recons_method == 4) THEN
      CALL ReconstructABPGradient_PPM(fcubehalo, recons, order)

    ELSE
      WRITE(*,*) 'Fatal Error: In ReconstructABPGradient'
      WRITE(*,*) 'Specified recons_method out of range. Given: ', recons_method
      WRITE(*,*) 'Valid values: 1, 2, 3, 4'
      STOP
    ENDIF

     ! Apply monotone filtering
    SELECT CASE (kmono)
    CASE (0) !Do nothing
      WRITE(*,*) "no filter applied to the reconstruction"
    CASE (1)

       !Simplest filter: just scales the recon so it's extreme value 
       !is no bigger than the original values of this point and its neighbors
      CALL MonotonizeABPGradient(fcubehalo, order, recons, .FALSE.)

    CASE (2)
       
       !Applies a more sophisticated Van Leer limiter (or, to be consistent, a filter)
       CALL VanLeerLimit(fcubehalo, order, recons)

    CASE (3)

       !Applies a selective filter
       CALL MonotonizeABPGradient(fcubehalo, order, recons, .TRUE.)

    CASE (4)

       !A filter that filters the linear part first
       CALL MonotonizeABPGradient_New(fcubehalo, order, recons)

    CASE DEFAULT
       WRITE(*,*) "Limiter kmono = ", kmono, " does not exist."
       STOP 1201

    END SELECT

    !Apply positive-definite filtering, if desired. This should 
    !ONLY be applied to the S-L method, since the flux-form 
    !method needs something different done. (In particular, using 
    !positive-definite reconstructions does not ensure that a flux-
    !form scheme is positive definite, since we could get negatives 
    !when subtracting the resulting fluxes.)
    !HOWEVER...we will allow this to be enabled, for testing purposes
    IF ( (kpd > 0 .AND. kscheme == 2) .OR. (kpd == 2 .AND. kscheme == 4) ) THEN
      WRITE(*,*) "applying positive deifnite constraint"
       CALL PosDefABPGradient(fcubehalo, order, recons)
    END IF


  END SUBROUTINE



!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! SUBROUTINE AdjustLimiter
!
! Description:
!   Adjust the slope limiter based on new point values.
!
! Parameters:
!   value - Point value
!   element_value - Value at the center of the element
!   local_max - Local maximum value of the function (from neighbours)
!   local_min - Local minimum value of the function (to neighbours)
!   min_phi (INOUT) - Slope limiter
!------------------------------------------------------------------------------
  SUBROUTINE AdjustLimiter(value, element_value, local_min, local_max, min_phi)

    IMPLICIT NONE

    REAL (KIND=dbl_kind), INTENT(IN)    :: value, element_value
    REAL (KIND=dbl_kind), INTENT(IN)    :: local_min, local_max
    REAL (KIND=dbl_kind), INTENT(INOUT) :: min_phi

    ! Local variables
    REAL (KIND=dbl_kind) :: phi = 0.0

    IF ((local_min > element_value ) .OR. (local_max < element_value )) THEN
      WRITE (*,*) 'Fatal Error: In AdjustLimiter'
      WRITE (*,*) 'Local min: ', local_min, ' max: ', local_max
      WRITE (*,*) 'Elemn: ', element_value
      STOP
    ENDIF

    ! Check against the minimum bound on the reconstruction
    IF (value - element_value > tiny * value) THEN
      phi = (local_max - element_value) / &
            (value - element_value)

      min_phi = MIN(min_phi, phi)

    ! Check against the maximum bound on the reconstruction
    ELSEIF (value - element_value < -tiny * value) THEN
      phi = (local_min - element_value) / &
            (value - element_value)

      min_phi = MIN(min_phi, phi)

    ENDIF

    IF (min_phi < 0.0) THEN
      WRITE (*,*) 'Fatal Error: In AdjustLimiter'
      WRITE (*,*) 'Min_Phi: ', min_phi
      WRITE (*,*) 'Phi: ', phi
      WRITE (*,*) 'Value: ', value
      WRITE (*,*) 'Elemn: ', element_value
      WRITE (*,*) 'Val-E: ', value - element_value
      STOP
    ENDIF

  END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE VanLeerLimit
!
! Description:
!   Apply a 2D Van Leer-type limiter to a reconstruction. This acts ONLY 
!   on the linear part of the reconstruction , if any. If passed a PCoM 
!   reconstruction, this just returns without altering the recon.
!
! Parameters:
!   fcubehalo - Scalar field on the cubed sphere to use in reconstruction
!   order - Order of the reconstruction
!   recons (INOUT) - Array of reconstructed coefficients
!
! Remarks:
!   The Van Leer Limiter described here is given on pages 328--329 
!   of Dukowicz and Baumgardner (2000). There are no guarantees 
!   on what it will do to PPM.
!------------------------------------------------------------------------------
  SUBROUTINE VanLeerLimit(fcubehalo, order, recons)


    IMPLICIT NONE

    REAL (KIND=dbl_kind), DIMENSION(-1:ncube_reconstruct+1, -1:ncube_reconstruct+1, 6), &
                          INTENT(IN) :: fcubehalo

    INTEGER (KIND=int_kind), INTENT(IN) :: order

    REAL (KIND=dbl_kind), DIMENSION(:,:,:,:), INTENT(INOUT) :: recons

    ! Local variables
    INTEGER (KIND=int_kind) :: i, j, k, m, n

    REAL (KIND=dbl_kind) :: local_min, local_max, value, phi, min_phi, &
         recon_min, recon_max

    ! The first-order piecewise constant scheme is monotone by construction
    IF (order == 1) THEN
      RETURN
    ENDIF

    ! Apply monotone limiting
    DO k = 1, 6
    DO j = 1, ncube_reconstruct-1
    DO i = 1, ncube_reconstruct-1
       CALL ABPHaloMinMax(fcubehalo, i, j, k, local_min, local_max,.FALSE.)

       ! Initialize the limiter
       min_phi = one

       ! For the second-order calculation, the minima and maxima will occur
       ! at the corner points of the element. For the Van Leer limiter, we 
       !wish to find BOTH of the reconstruction extrema.
       recon_min = bignum
       recon_max = -bignum

       DO m = i, i+1
       DO n = j, j+1

          ! Evaluate the function at each corner point
          CALL EvaluateABPReconstruction(                                &
               fcubehalo, recons, i, j, k, gp(m), gp(n), order, value)
          recon_min = MIN(recon_min, value)
          recon_max = MAX(recon_max, value)
          
       ENDDO
       ENDDO

       !This is equation 27 in Dukowicz and Baumgardner 2000
       min_phi = MIN(one, MAX(0.0, (local_min - fcubehalo(i,j,k))/(recon_min - fcubehalo(i,j,k))), &
            MAX(0.0, (local_max - fcubehalo(i,j,k))/(recon_max - fcubehalo(i,j,k))) )

       IF ((min_phi < -tiny) .OR. (min_phi > one + tiny)) THEN
          WRITE (*,*) 'Fatal Error: In MonotonizeABPGradient'
          WRITE (*,*) 'Slope limiter out of range: ', min_phi
          STOP
       ENDIF

       ! Apply monotone limiter to all reconstruction coefficients
       recons(1,i,j,k) = min_phi * recons(1,i,j,k)
       recons(2,i,j,k) = min_phi * recons(2,i,j,k)
          
    END DO
    END DO
    END DO




  END SUBROUTINE VanLeerLimit

  !------------------------------------------------------------------------------
  ! SUBROUTINE EquiangularElementArea
  !
  ! Description:
  !   Compute the area of a single equiangular cubed sphere grid cell.
  !
  ! Parameters: 
  !   alpha - Alpha coordinate of lower-left corner of grid cell
  !   da - Delta alpha
  !   beta - Beta coordinate of lower-left corner of grid cell
  !   db - Delta beta
  !------------------------------------------------------------------------------
  REAL(KIND=dbl_kind) FUNCTION EquiangularElementArea(alpha, da, beta, db)

    IMPLICIT NONE

!    REAL (kind=dbl_kind) :: EquiangularElementArea
    REAL (kind=dbl_kind) :: alpha, da, beta, db
    REAL (kind=dbl_kind) :: a1, a2, a3, a4

    ! Calculate interior grid angles
    a1 =      EquiangularGridAngle(alpha   , beta   )
    a2 = pi - EquiangularGridAngle(alpha+da, beta   )
    a3 = pi - EquiangularGridAngle(alpha   , beta+db)
    a4 =      EquiangularGridAngle(alpha+da, beta+db)

    ! Area = r*r*(-2*pi+sum(interior angles))
    EquiangularElementArea = -pi2 + a1 + a2 + a3 + a4

  END FUNCTION EquiangularElementArea

  !------------------------------------------------------------------------------
  ! FUNCTION EquiangularGridAngle
  !
  ! Description:
  !   Compute the angle between equiangular cubed sphere projection grid lines.
  !
  ! Parameters: 
  !   alpha - Alpha coordinate of evaluation point
  !   beta - Beta coordinate of evaluation point
  !------------------------------------------------------------------------------
  REAL(KIND=dbl_kind) FUNCTION EquiangularGridAngle(alpha, beta)
    IMPLICIT NONE
    REAL (kind=dbl_kind) :: alpha, beta
    EquiangularGridAngle = ACOS(-SIN(alpha) * SIN(beta))
  END FUNCTION EquiangularGridAngle

!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereFillHalo
!
! Description:
!   Recompute the cubed sphere data storage array, with the addition of a
!   halo region around the specified panel.
!
! Parameters:
!   parg - Current panel values
!   zarg (OUT) - Calculated panel values with halo/ghost region
!   np - Panel number
!   ncube - Dimension of the cubed sphere (# of grid lines)
!   nhalo - Number of halo/ghost elements around each panel
!------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo(parg, zarg, np, ncube, nhalo)

    IMPLICIT NONE

    REAL (KIND=dbl_kind), DIMENSION(ncube-1, ncube-1, 6), INTENT(IN) :: parg

    REAL (KIND=dbl_kind),                                            &
         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6), &
         INTENT(OUT) :: zarg

    INTEGER (KIND=int_kind), INTENT(IN) :: np, ncube,nhalo

    ! Local variables
    INTEGER (KIND=int_kind)               :: jh,jhy

    !zarg = 0.0 !DBG
    zarg(1:ncube-1,1:ncube-1,np) = parg(1:ncube-1,1:ncube-1,np)

    zarg(1-nhalo:0,1-nhalo:0,np) = 0.0
    zarg(1-nhalo:0,ncube:ncube+nhalo-1,np) = 0.0
    zarg(ncube:ncube+nhalo-1,1-nhalo:0,np) = 0.0
    zarg(ncube:ncube+nhalo-1,ncube:ncube+nhalo-1,np) = 0.0

    ! Equatorial panels
    IF (np==1) THEN
       DO jh=1,nhalo
          zarg(ncube+jh-1,1:ncube-1 ,1) = parg(jh        ,1:ncube-1 ,2)  !exchange right
          zarg(1-jh      ,1:ncube-1 ,1) = parg(ncube-jh  ,1:ncube-1 ,4)  !exchange left
          zarg(1:ncube-1 ,1-jh      ,1) = parg(1:ncube-1 ,ncube-jh  ,5)  !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,1) = parg(1:ncube-1 ,jh        ,6)  !exchange over
       ENDDO

    ELSE IF (np==2) THEN
       DO jh=1,nhalo
          zarg(1-jh      ,1:ncube-1 ,2) = parg(ncube-jh,1:ncube-1   ,1)  !exchange left
          zarg(ncube+jh-1,1:ncube-1 ,2) = parg(jh      ,1:ncube-1   ,3)  !exchange right
          zarg(1:ncube-1 ,1-jh      ,2) = parg(ncube-jh,ncube-1:1:-1,5)  !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,2) = parg(ncube-jh,1:ncube-1   ,6)  !exchange over
       ENDDO

    ELSE IF (np==3) THEN
       DO jh=1,nhalo
          zarg(ncube+jh-1,1:ncube-1 ,3) = parg(jh          ,1:ncube-1,4)  !exchange right
          zarg(1-jh      ,1:ncube-1 ,3) = parg(ncube-jh    ,1:ncube-1,2)  !exchange left
          zarg(1:ncube-1 ,1-jh      ,3) = parg(ncube-1:1:-1,jh       ,5)  !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,3) = parg(ncube-1:1:-1,ncube-jh ,6)  !exchange over
       ENDDO

    ELSE IF (np==4) THEN
       DO jh=1,nhalo
          zarg(1-jh      ,1:ncube-1 ,4) = parg(ncube-jh,1:ncube-1   ,3) !exchange left
          zarg(ncube+jh-1,1:ncube-1 ,4) = parg(jh      ,1:ncube-1   ,1) !exchange right
          zarg(1:ncube-1 ,1-jh      ,4) = parg(jh      ,1:ncube-1   ,5) !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,4) = parg(jh      ,ncube-1:1:-1,6) !exchange over
       ENDDO

    ! Bottom panel
    ELSE IF (np==5) THEN
       DO jh=1,nhalo
          zarg(1-jh      ,1:ncube-1 ,5) = parg(1:ncube-1   ,jh      ,4) !exchange left
          zarg(ncube+jh-1,1:ncube-1 ,5) = parg(ncube-1:1:-1,jh      ,2) !exchange right
          zarg(1:ncube-1 ,1-jh      ,5) = parg(ncube-1:1:-1,jh      ,3) !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,5) = parg(1:ncube-1   ,jh      ,1) !exchange over
       ENDDO

    ! Top panel
    ELSE IF (np==6) THEN
       DO jh=1,nhalo
          zarg(1-jh      ,1:ncube-1 ,6) = parg(ncube-1:1:-1,ncube-jh,4) !exchange left
          zarg(ncube+jh-1,1:ncube-1 ,6) = parg(1:ncube-1   ,ncube-jh,2) !exchange right
          zarg(1:ncube-1 ,1-jh      ,6) = parg(1:ncube-1   ,ncube-jh,1) !exchange below
          zarg(1:ncube-1 ,ncube+jh-1,6) = parg(ncube-1:1:-1,ncube-jh,3) !exchange over
       ENDDO

    ELSE
       WRITE (*,*) 'Fatal error: In CubedSphereFillHalo'
       WRITE (*,*) 'Invalid panel id ', np
       STOP
    ENDIF

  END SUBROUTINE CubedSphereFillHalo

!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereFillHalo_Linear
!
! Description:
!   Recompute the cubed sphere data storage array, with the addition of a
!   2-element halo region around the specified panel.  Use linear order
!   interpolation to translate between panels.
!
! Parameters:
!   parg - Current panel values
!   zarg (OUT) - Calculated panel values with halo/ghost region
!   np - Panel number
!   ncube - Dimension of the cubed sphere (# of grid lines)
!------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo_Linear(parg, zarg, np, ncube)

!    USE CubedSphereTrans  ! Cubed sphere transforms

    IMPLICIT NONE

    INTEGER (KIND=int_kind), PARAMETER :: nhalo = 2

    REAL (KIND=dbl_kind), DIMENSION(ncube-1, ncube-1, 6), INTENT(IN) :: parg

    REAL (KIND=dbl_kind),                                            &
         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6), &
         INTENT(OUT) :: zarg

    INTEGER (KIND=int_kind), INTENT(IN) :: np, ncube

    ! Local variables
    INTEGER (KIND=int_kind) :: ii, iref, jj, ipanel, imin, imax
    REAL    (KIND=dbl_kind) :: width, lon, lat, beta, a, newbeta

    REAL    (KIND=dbl_kind), DIMENSION(0:ncube, nhalo) :: prealpha
    REAL    (KIND=dbl_kind), DIMENSION(0:ncube, nhalo) :: newalpha

    REAL (KIND=dbl_kind), &
         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6) :: yarg

    ! Use 0.0 order interpolation to begin
    CALL CubedSphereFillHalo(parg, yarg, np, ncube, nhalo)

    zarg(:,:,np) = yarg(:,:,np)

    ! Calculate the overlapping alpha coordinates
    width = pih / DBLE(ncube-1)

    DO jj = 1, nhalo
      DO ii = 0, ncube
        prealpha(ii, jj) = width * (DBLE(ii-1) + 0.5) - piq
        beta = - width * (DBLE(jj-1) + 0.5) - piq

        CALL CubedSphereABPFromABP(prealpha(ii,jj), beta, 1, 5, &
                                   newalpha(ii,jj), newbeta)
      ENDDO
    ENDDO

    ! Now apply linear interpolation to obtain edge components
    DO jj = 1, nhalo
      ! Reset the reference index
      iref = 2 

      ! Interpolation can be applied to more elements after first band
      IF (jj == 1) THEN
        imin = 1
        imax = ncube-1
      ELSE
        imin = 0
        imax = ncube
      ENDIF

      ! Apply linear interpolation
      DO ii = imin, imax
        DO WHILE ((iref .NE. ncube-1) .AND. &
                  (newalpha(ii,jj) > prealpha(iref,jj)))
          iref = iref + 1
        ENDDO

        IF ((newalpha(ii,jj) > prealpha(iref-1,jj)) .AND.    &
            (newalpha(ii,jj) .LE. prealpha(iref  ,jj)))      &
        THEN
          a = (newalpha(ii,jj)   - prealpha(iref-1,jj)) / &
              (prealpha(iref,jj) - prealpha(iref-1,jj))

          IF ((a < 0.0) .OR. (a > one)) THEN
            WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
            WRITE (*,*) 'a out of bounds'
            STOP
          ENDIF

          ! Bottom edge of panel
          zarg(ii, 1-jj, np) =                   &
            (one - a) * yarg(iref-1, 1-jj, np) + &
                   a  * yarg(iref, 1-jj, np)

          ! Left edge of panel
          zarg(1-jj, ii, np) =                   &
            (one - a) * yarg(1-jj, iref-1, np) + &
                   a  * yarg(1-jj, iref, np)

          ! Top edge of panel
          zarg(ii, ncube+jj-1, np) =                   &
            (one - a) * yarg(iref-1, ncube+jj-1, np) + &
                   a  * yarg(iref, ncube+jj-1, np)

          ! Right edge of panel
          zarg(ncube+jj-1, ii, np) =                   &
            (one - a) * yarg(ncube+jj-1, iref-1, np) + &
                   a  * yarg(ncube+jj-1, iref, np)

        ELSE
          WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
          WRITE (*,*) 'ii: ', ii, ' jj: ', jj
          WRITE (*,*) 'newalpha: ', newalpha(ii,jj)
          WRITE (*,*) 'prealpha: ', prealpha(iref-1,jj), '-', prealpha(iref,jj)
          STOP
        ENDIF
      ENDDO
    ENDDO

    ! Fill in corner bits
    zarg(0, 0, np) =                         &
      0.25 * (zarg(1,0,np) + zarg(0,1,np) + &
               zarg(-1,0,np) + zarg(0,-1,np))
    zarg(0, ncube, np) =                                 &
      0.25 * (zarg(0,ncube-1,np) + zarg(0,ncube+1,np) + &
               zarg(-1,ncube,np)  + zarg(1,ncube,np))
    zarg(ncube, 0, np) =                                 &
      0.25 * (zarg(ncube-1,0,np) + zarg(ncube+1,0,np) + &
               zarg(ncube,-1,np)  + zarg(ncube,1,np))
    zarg(ncube, ncube, np) =                                     &
      0.25 * (zarg(ncube-1,ncube,np) + zarg(ncube+1,ncube,np) + &
               zarg(ncube,ncube-1,np) + zarg(ncube,ncube+1,np))

  END SUBROUTINE CubedSphereFillHalo_Linear

!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereFillHalo_Cubic
!
! Description:
!   Recompute the cubed sphere data storage array, with the addition of a
!   2-element halo region around the specified panel.  Use higher order 
!   interpolation to translate between panels.
!
! Parameters:
!   parg - Current panel values
!   zarg (OUT) - Calculated panel values with halo/ghost region
!   np - Panel number
!   ncube - Dimension of the cubed sphere (# of grid lines)
!------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo_Cubic(parg, zarg, np, ncube)

!    USE CubedSphereTrans  ! Cubed sphere transforms
!    USE MathUtils         ! Has function for 1D cubic interpolation

    IMPLICIT NONE

    INTEGER (KIND=int_kind), PARAMETER :: nhalo = 2

    REAL (KIND=dbl_kind), DIMENSION(ncube-1, ncube-1, 6), INTENT(IN) :: parg

    REAL (KIND=dbl_kind),                                            &
         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6), &
         INTENT(OUT) :: zarg

    INTEGER (KIND=int_kind), INTENT(IN) :: np, ncube

    ! Local variables
    INTEGER (KIND=int_kind) :: ii, iref, ibaseref, jj, ipanel, imin, imax
    REAL    (KIND=dbl_kind) :: width, lon, lat, beta, a, newbeta

    REAL    (KIND=dbl_kind), DIMENSION(0:ncube, nhalo) :: prealpha
    REAL    (KIND=dbl_kind), DIMENSION(0:ncube, nhalo) :: newalpha
    REAL    (KIND=dbl_kind), DIMENSION(1:4) :: C, D, X

    REAL (KIND=dbl_kind), &
         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6) :: yarg

    ! Use 0.0 order interpolation to begin
    CALL CubedSphereFillHalo(parg, yarg, np, ncube, nhalo)

    zarg(:,:,np) = yarg(:,:,np)

    ! Calculate the overlapping alpha coordinates
    width = pih / DBLE(ncube-1)

    DO jj = 1, nhalo
      DO ii = 0, ncube
        !
        ! alpha,beta for the cell center (extending the panel)
        !
        prealpha(ii, jj) = width * (DBLE(ii-1) + 0.5) - piq
        beta = - width * (DBLE(jj-1) + 0.5) - piq

        CALL CubedSphereABPFromABP(prealpha(ii,jj), beta, 1, 5, &
                                   newalpha(ii,jj), newbeta)
      ENDDO
    ENDDO

    ! Now apply cubic interpolation to obtain edge components
    DO jj = 1, nhalo
      ! Reset the reference index, which gives the element in newalpha that
      ! is closest to ii, looking towards larger values of alpha.
      iref = 2 

      ! Interpolation can be applied to more elements after first band
!      IF (jj == 1) THEN
!        imin = 1
!        imax = ncube-1
!      ELSE
        imin = 0
        imax = ncube
!      ENDIF

      ! Apply cubic interpolation
      DO ii = imin, imax
        DO WHILE ((iref .NE. ncube-1) .AND. &
                  (newalpha(ii,jj) > prealpha(iref,jj)))
          iref = iref + 1
        ENDDO

        ! Smallest index for cubic interpolation - apply special consideration
        IF (iref == 2) THEN
          ibaseref = iref-1

        ! Largest index for cubic interpolation - apply special consideration
        ELSEIF (iref == ncube-1) THEN
          ibaseref = iref-3

        ! Normal range
        ELSE
          ibaseref = iref-2
        ENDIF

        ! Bottom edge of panel
        zarg(ii, 1-jj, np) =                                &
          CUBIC_EQUISPACE_INTERP(                           &
            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
            yarg(ibaseref:ibaseref+3, 1-jj, np))

        ! Left edge of panel
        zarg(1-jj, ii, np) =                                &
          CUBIC_EQUISPACE_INTERP(                           &
            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
            yarg(1-jj, ibaseref:ibaseref+3, np))

        ! Top edge of panel
        zarg(ii, ncube+jj-1, np) =                          &
          CUBIC_EQUISPACE_INTERP(                           &
            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
            yarg(ibaseref:ibaseref+3, ncube+jj-1, np))

        ! Right edge of panel
        zarg(ncube+jj-1, ii, np) =                          &
          CUBIC_EQUISPACE_INTERP(                           &
            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
            yarg(ncube+jj-1, ibaseref:ibaseref+3, np))

      ENDDO
    ENDDO

    ! Fill in corner bits
    zarg(0, 0, np) =                         &
      0.25 * (zarg(1,0,np) + zarg(0,1,np) + &
               zarg(-1,0,np) + zarg(0,-1,np))
    zarg(0, ncube, np) =                                 &
      0.25 * (zarg(0,ncube-1,np) + zarg(0,ncube+1,np) + &
               zarg(-1,ncube,np)  + zarg(1,ncube,np))
    zarg(ncube, 0, np) =                                 &
      0.25 * (zarg(ncube-1,0,np) + zarg(ncube+1,0,np) + &
               zarg(ncube,-1,np)  + zarg(ncube,1,np))
    zarg(ncube, ncube, np) =                                     &
      0.25 * (zarg(ncube-1,ncube,np) + zarg(ncube+1,ncube,np) + &
               zarg(ncube,ncube-1,np) + zarg(ncube,ncube+1,np))

  END SUBROUTINE CubedSphereFillHalo_Cubic

!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereABPFromABP
!
! Description:
!   Determine the (alpha,beta,idest) coordinate of a source point on
!   panel isource.
!
! Parameters:
!   alpha_in - Alpha coordinate in
!   beta_in - Beta coordinate in
!   isource - Source panel
!   idest - Destination panel
!   alpha_out (OUT) - Alpha coordinate out
!   beta_out (OUT) - Beta coordiante out
!------------------------------------------------------------------------------
  SUBROUTINE CubedSphereABPFromABP(alpha_in,  beta_in, isource, idest, &
                                   alpha_out, beta_out)

    IMPLICIT NONE

    REAL    (KIND=dbl_kind), INTENT(IN)  :: alpha_in, beta_in
    INTEGER (KIND=int_kind), INTENT(IN)  :: isource, idest
    REAL    (KIND=dbl_kind), INTENT(OUT) :: alpha_out, beta_out

    ! Local variables
    REAL    (KIND=dbl_kind) :: a1, b1
    REAL    (KIND=dbl_kind) :: xx, yy, zz
    REAL    (KIND=dbl_kind) :: sx, sy, sz

    ! Convert to relative Cartesian coordinates
    a1 = TAN(alpha_in)
    b1 = TAN(beta_in)

    sz = (one + a1 * a1 + b1 * b1)**(-0.5)
    sx = sz * a1
    sy = sz * b1

    ! Convert to full Cartesian coordinates
    IF (isource == 6) THEN
      yy = sx; xx = -sy; zz = sz

    ELSEIF (isource == 5) THEN
      yy = sx; xx = sy; zz = -sz

    ELSEIF (isource == 1) THEN
      yy = sx; zz = sy; xx = sz

    ELSEIF (isource == 3) THEN
      yy = -sx; zz = sy; xx = -sz

    ELSEIF (isource == 2) THEN
      xx = -sx; zz = sy; yy = sz

    ELSEIF (isource == 4) THEN
      xx = sx; zz = sy; yy = -sz

    ELSE
      WRITE(*,*) 'Fatal Error: Source panel invalid in CubedSphereABPFromABP'
      WRITE(*,*) 'panel = ', isource
      STOP
    ENDIF

    ! Convert to relative Cartesian coordinates on destination panel
    IF (idest == 6) THEN
      sx = yy; sy = -xx; sz = zz

    ELSEIF (idest == 5) THEN
      sx = yy; sy = xx; sz = -zz

    ELSEIF (idest == 1) THEN
      sx = yy; sy = zz; sz = xx

    ELSEIF (idest == 3) THEN
      sx = -yy; sy = zz; sz = -xx

    ELSEIF (idest == 2) THEN
      sx = -xx; sy = zz; sz = yy

    ELSEIF (idest == 4) THEN
      sx = xx; sy = zz; sz = -yy

    ELSE
      WRITE(*,*) 'Fatal Error: Dest panel invalid in CubedSphereABPFromABP'
      WRITE(*,*) 'panel = ', idest
      STOP
    ENDIF
    IF (sz < 0) THEN
      WRITE(*,*) 'Fatal Error: In CubedSphereABPFromABP'
      WRITE(*,*) 'Invalid relative Z coordinate'
      STOP
    ENDIF

    ! Use panel information to calculate (alpha, beta) coords
    alpha_out = ATAN(sx / sz)
    beta_out = ATAN(sy / sz)

  END SUBROUTINE


!------------------------------------------------------------------------------
! FUNCTION CUBIC_EQUISPACE_INTERP
!
! Description:
!   Apply cubic interpolation on the specified array of values, where all
!   points are equally spaced.
!
! Parameters:
!   dx - Spacing of points
!   x - X coordinate where interpolation is to be applied
!   y - Array of 4 values = f(x + k * dx) where k = 0,1,2,3
!------------------------------------------------------------------------------
  FUNCTION CUBIC_EQUISPACE_INTERP(dx, x, y)
    
    IMPLICIT NONE
    
    REAL (KIND=dbl_kind) :: CUBIC_EQUISPACE_INTERP
    REAL (KIND=dbl_kind) :: dx, x
    REAL (KIND=dbl_kind), DIMENSION(1:4) :: y
    
    CUBIC_EQUISPACE_INTERP =                                                   &
         (-y(1) / (6.0 * dx**3)) * (x - dx) * (x - 2.0 * dx) * (x - 3.0 * dx) + &
         ( y(2) / (2.0 * dx**3)) * (x) * (x - 2.0 * dx) * (x - 3.0 * dx) +      &
         (-y(3) / (2.0 * dx**3)) * (x) * (x - dx) * (x - 3.0 * dx) +            &
         ( y(4) / (6.0 * dx**3)) * (x) * (x - dx) * (x - 2.0 * dx)
    
  END FUNCTION CUBIC_EQUISPACE_INTERP
  
!  FUNCTION I_10_ab(alpha,beta)
!    IMPLICIT NONE
!    REAL (KIND=dbl_kind) :: I_10_AB
!    REAL (KIND=dbl_kind), INTENT(IN) :: alpha,beta
!    I_10_ab = -ASINH(COS(alpha) * TAN(beta))
!  END FUNCTION I_10_AB
!!
!
!  REAL (KIND=dbl_kind) FUNCTION I_01_ab(alpha,beta)
!    IMPLICIT NONE
!    REAL (KIND=dbl_kind), INTENT(IN) :: alpha,beta
!    I_01_ab = -ASINH(COS(beta) * TAN(alpha))
!  END FUNCTION I_01_AB
!
!  REAL (KIND=dbl_kind) FUNCTION I_20_ab(alpha,beta)
!    IMPLICIT NONE
!    REAL (KIND=dbl_kind), INTENT(IN) :: alpha,beta
!
!    I_20_ab = TAN(beta)*ASINH(COS(beta)*TAN(alpha))+ACOS(SIN(alpha)*SIN(beta))
!  END FUNCTION I_20_AB
!
!  REAL (KIND=dbl_kind) FUNCTION I_02_ab(alpha,beta)
!    IMPLICIT NONE
!    REAL (KIND=dbl_kind), INTENT(IN) :: alpha,beta
!    
!    I_02_ab = TAN(alpha)*ASINH(TAN(beta)*COS(alpha))+ACOS(SIN(alpha)*SIN(beta))
!  END FUNCTION I_02_AB
! 
!  REAL (KIND=dbl_kind) FUNCTION I_11_ab(alpha,beta)
!    IMPLICIT NONE
!    REAL (KIND=dbl_kind), INTENT(IN) :: alpha,beta
!    
!    I_11_ab = -SQRT(1.0+TAN(alpha)**2+TAN(beta)**2)
!  END FUNCTION I_11_AB
!


END MODULE reconstruct

