module KineticsMod
  ! !DESCRIPTION:
  ! Subroutines to do substrate kinetics
  ! Created by Jinyun Tang, Apr 11, 2013
  ! !USES:

  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun
  use clm_varctl,   only: iulog
  implicit none
  real(r8),public, parameter :: kd_infty = 1.e40_r8      !internal parameter

  interface mmcomplex   !the m-m kinetics
     module procedure mmcomplex_v1s,mmcomplex_v1e, mmcomplex_m
  end interface mmcomplex

  interface ecacomplex  !the eca kinetics
     module procedure ecacomplex_v1s,ecacomplex_v1e, ecacomplex_m
  end interface ecacomplex

  interface ecacomplex_cell_norm  !the eca kinetics
     module procedure ecacomplex_cell_norm_v1s,ecacomplex_cell_norm_v1e, ecacomplex_cell_norm_m
  end interface ecacomplex_cell_norm

contains
  !-------------------------------------------------------------------------------
  subroutine mmcomplex_v1s(kd,ee,ss,siej)

    ! !DESCRIPTION:
    ! Compute concentrations of the enzyme substrate complexes
    ! many microbes vs single substrate
    ! using the traditional M-M kinetics

    ! !USES:

    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in)  :: kd
    real(r8), dimension(:), intent(in)  :: ee
    real(r8), intent(in)                :: ss
    real(r8), dimension(:), intent(out) :: siej

    ! !LOCAL VARIABLES:
    integer :: jj, j
    real(r8) :: dS
    jj = size(ee)
    siej = 0._r8

    do j = 1, jj
       if(kd(j)>0._r8 .and. (kd(j)<.9*kd_infty))then
          siej(j) = ss * ee(j) / (kd(j) + ss)
       endif
    enddo
    ds = sum(siej)
    if(ds>ss)then
       do j = 1, jj
          siej(j) = siej(j) * ss / ds
       enddo
    endif

  end subroutine mmcomplex_v1s
  !-------------------------------------------------------------------------------
  subroutine mmcomplex_v1e(kd,ee,ss,siej)
    ! !DESCRIPTION:
    !compute concentrations of the enzyme substrate complexes
    !using the traditional M-M kinetics
    !many substrates vs single microbe

    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in)  :: kd
    real(r8), dimension(:), intent(in)  :: ss
    real(r8), intent(in)                :: ee
    real(r8), dimension(:), intent(out) :: siej

    ! !LOCAL VARIABLES
    integer :: ii
    integer :: i
    real(r8) :: dE
    ii = size(ss)
    siej = 0._r8

    do i = 1, ii
       if(kd(i)>0._r8 .and. (kd(i)<.9*kd_infty))then
          siej(i) = ss(i) * ee / (kd(i) + ss(i))
       endif
    enddo
    dE = sum(siej)
    if(dE>ee)then
       do i = 1, ii
          siej(i) = siej(i) * ee / dE
       enddo
    endif

  end subroutine mmcomplex_v1e
  !-------------------------------------------------------------------------------
  subroutine mmcomplex_m(kd,ee,ss,siej)
    ! !DESCRIPTION:
    !compute concentrations of the enzyme substrate complexes
    !using the traditional M-M kinetics
    !many substrates vs many microbes

    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:,:), intent(in)  :: kd
    real(r8), dimension(:), intent(in)    :: ee, ss
    real(r8), dimension(:,:), intent(out) :: siej

    ! !LOCAL VARIABLES:
    integer :: ii,jj
    integer :: i, j
    real(r8) :: dS, dE

    ii = size(ss)
    jj = size(ee)
    siej = 0._r8
    do i = 1, ii
       do j = 1, jj
          if(kd(i,j)>0._r8 .and. (kd(i,j)<.9*kd_infty))then
             siej(i,j) = ss(i) * ee(j) / (kd(i,j) + ss(i))
          endif
       enddo
       ds = sum(siej(i,:))
       if(ds>ss(i))then
          do j = 1, jj
             siej(i,j) = siej(i,j) * ss(i) / ds
          enddo
       endif
    enddo

    do j = 1, jj
       dE = sum(siej(:,j))
       if(dE>ee(j))then
          do i = 1, ii
             siej(i,j) = siej(i,j) * ee(j) / dE
          enddo
       endif
    enddo
  end subroutine mmcomplex_m
  !-------------------------------------------------------------------------------
  subroutine ecacomplex_v1s(kd,ss,ee,siej)
    ! !DESCRIPTION:
    !compute concentrations of the enzyme substrate complexes
    !using the first order accurate ECA kinetics
    !many microbes vs one substrate
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in)  :: kd
    real(r8), dimension(:), intent(in)  :: ee
    real(r8), intent(in)                :: ss
    real(r8), dimension(:), intent(out) :: siej

    ! !LOCAL VARIABLES:
    integer :: jj
    integer :: j
    real(r8) :: dnm2

    jj = size(ee)
    siej = 0._r8

    dnm2=1._r8
    do j = 1, jj
       if(kd(j)>0._r8 .and. (kd(j)<.9*kd_infty))then
          dnm2=dnm2 + ee(j)/kd(j)
       endif
    enddo
    do j = 1, jj
       if(kd(j)>0._r8 .and. (kd(j)<.9*kd_infty))then
          siej(j) = ss*ee(j)/(kd(j)*(dnm2+ss/kd(j)))
       endif
    enddo
  end subroutine ecacomplex_v1s
!-------------------------------------------------------------------------------
   subroutine ecacomplex_v1e(kd,ss,ee,siej)
   ! !DESCRIPTION:
   !compute concentrations of the enzyme substrate complexes
   !using the first order accurate ECA kinetics
   !many substrate vs single microbe
   implicit none
   ! !ARGUMENTS:
   real(r8), dimension(:), intent(in)  :: kd
   real(r8), dimension(:), intent(in)  :: ss
   real(r8), intent(in)                :: ee
   real(r8), dimension(:), intent(out) :: siej

   ! !LOCAL VARIABLES:
   integer :: ii
   integer :: i
   real(r8) :: dnm1

   ii = size(ss)
   siej = 0._r8
   dnm1=1._r8
   do i = 1, ii
      if(kd(i)>0._r8 .and. (kd(i)<.9*kd_infty))then
         dnm1 = dnm1 + ss(i)/kd(i)
      endif
   enddo
   do i = 1, ii
      if(kd(i)>0._r8 .and. (kd(i)<.9*kd_infty))then
         siej(i) = ss(i)*ee/(kd(i)*(dnm1+ee/kd(i)))
      endif
   enddo
   end subroutine ecacomplex_v1e
   !-------------------------------------------------------------------------------
   subroutine ecacomplex_m(kd,ss,ee,siej)
     ! !DESCRIPTION:
     !compute concentrations of the enzyme substrate complexes
     !using the first order accurate ECA kinetics
     ! many substrate vs many enzymes
     implicit none
     ! !ARGUMENTS:
     real(r8), dimension(:,:), intent(in)  :: kd
     real(r8), dimension(:), intent(in)    :: ee, ss
     real(r8), dimension(:,:), intent(out) :: siej

     ! !LOCAL VARIABLES:
     integer :: ii,jj
     integer :: i, j, k
     real(r8) :: dnm1, dnm2

     ii = size(ss)       !number of substrates, dim 1
     jj = size(ee)       !number of enzymes, dim2
     if(ii/=size(siej,1) .or. jj/=size(siej,2))then
        write(iulog,*)'wrong matrix shape in ecacomplex_m'
        write(iulog,*)'clm model is stopping'
        call endrun()
     endif
     siej = 0._r8
     do i = 1, ii
        dnm1 = 0._r8
        do k = 1, jj
           if(kd(i,k)>0._r8 .and. (kd(i,k)<.9*kd_infty))then
              dnm1 = dnm1 + ee(k)/kd(i,k)
           endif
        enddo
        do j = 1, jj
           dnm2 = 0._r8
           if(kd(i,j)>0._r8 .and. (kd(i,j)<.9*kd_infty) )then
              do k = 1, ii
                 if(kd(k,j)>0._r8 .and. (kd(k,j)<.9*kd_infty))then
                    dnm2=dnm2 + ss(k)/kd(k,j)
                 endif
              enddo
              siej(i,j) = ss(i)*ee(j)/(kd(i,j)*(1._r8+dnm1+dnm2))
           endif
        enddo
     enddo
   end subroutine ecacomplex_m
   !-------------------------------------------------------------------------------
   subroutine ecacomplex_cell_norm_m(kd,ss,ee,siej)
     ! !DESCRIPTION:
     ! compute concentrations of the enzyme substrate complexes
     ! using the first order accurate ECA kinetics
     ! and noramlize the return value with cell abundance
     ! many substrates vs many enzymes
     implicit none
     ! !ARGUMENTS:
     real(r8), dimension(:,:), intent(in)  :: kd
     real(r8), dimension(:), intent(in)    :: ee, ss
     real(r8), dimension(:,:), intent(out) :: siej
     ! !LOCAL VARIABLES:
     integer :: ii,jj
     integer :: i, j, k
     real(r8) :: dnm1, dnm2

     ii = size(ss)       !number of substrates, dim 1
     jj = size(ee)       !number of enzymes, dim2
     if(ii/=size(siej,1) .or. jj/=size(siej,2))then
        write(iulog,*)'wrong matrix shape in ecacomplex_m'
        write(iulog,*)'clm model is stopping'
        call endrun()
     endif
     siej = 0._r8
     do i = 1, ii
        dnm1 = 0._r8
        do k = 1, jj
           if(kd(i,k)>0._r8 .and. (kd(i,k)<.9*kd_infty))then
              dnm1 = dnm1 + ee(k)/kd(i,k)
           endif
        enddo
        do j = 1, jj
           dnm2 = 0._r8
           if(kd(i,j)>0._r8 .and. (kd(i,j)<.9*kd_infty))then
              do k = 1, ii
                 if(kd(k,j)>0._r8)then
                    dnm2=dnm2 + ss(k)/kd(k,j)
                 endif
              enddo
              siej(i,j) = ss(i)/(kd(i,j)*(1._r8+dnm1+dnm2))
           endif
        enddo
     enddo
   end subroutine ecacomplex_cell_norm_m

   !-------------------------------------------------------------------------------
   subroutine ecacomplex_cell_norm_v1s(kd,ss,ee,siej)
     ! !DESCRIPTION:
     !compute concentrations of the enzyme substrate complexes
     !using the first order accurate ECA kinetics
     !many microbes vs one substrate
     !and normalize return value with cell abundance

     implicit none
     ! !ARGUMENTS:
     real(r8), dimension(:), intent(in) :: kd
     real(r8), dimension(:), intent(in) :: ee
     real(r8), intent(in) :: ss
     real(r8), dimension(:), intent(out) :: siej
     ! !LOCAL VARIABLES:
     integer :: jj
     integer :: j
     real(r8) :: dnm2

     jj = size(ee)
     siej = 0._r8

     dnm2=1._r8
     do j = 1, jj
        if(kd(j)>0._r8 .and. (kd(j)<.9*kd_infty))then
           dnm2=dnm2 + ee(j)/kd(j)
        endif
     enddo
     do j = 1, jj
        if(kd(j)>0._r8 .and. (kd(j)<.9*kd_infty))then
           siej(j) = ss/(kd(j)*(dnm2+ss/kd(j)))
        endif
     enddo
   end subroutine ecacomplex_cell_norm_v1s
   !-------------------------------------------------------------------------------
   subroutine ecacomplex_cell_norm_v1e(kd,ss,ee,siej)
     ! !DESCRIPTION:
     ! compute concentrations of the enzyme substrate complexes
     ! using the first order accurate ECA kinetics
     ! many substrate vs single microbe
     ! and normalize the return value with cell abundance
     implicit none
     ! !ARGUMENTS:
     real(r8), dimension(:), intent(in)  :: kd
     real(r8), dimension(:), intent(in)  :: ss
     real(r8), intent(in)                :: ee
     real(r8), dimension(:), intent(out) :: siej
     ! !LOCAL VARIABLES:
     integer :: ii
     integer :: i
     real(r8) :: dnm1

     ii = size(ss)
     siej = 0._r8
     dnm1=1._r8
     do i = 1, ii
        if(kd(i)>0._r8  .and. (kd(i)<.9*kd_infty))then
           dnm1 = dnm1 + ss(i)/kd(i)
        endif
     enddo
     do i = 1, ii
        if(kd(i)>0._r8 .and. (kd(i)<.9*kd_infty))then
           siej(i) = ss(i)/(kd(i)*(dnm1+ee/kd(i)))
        endif
     enddo
   end subroutine ecacomplex_cell_norm_v1e
end module KineticsMod
