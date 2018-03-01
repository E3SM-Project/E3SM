module KineticsMod
  ! !DESCRIPTION:
  ! Subroutines to do substrate kinetics
  ! Created by Jinyun Tang, Apr 11, 2013
  ! !USES:
#include "bshr_assert.h"
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
  implicit none
  private
  character(len=*), parameter :: mod_filename = &
       __FILE__

  real(r8),public, parameter :: kd_infty = 1.e40_r8      !internal parameter
  public :: mmcomplex, ecacomplex, ecacomplex_cell_norm

  interface mmcomplex   !the m-m kinetics
     module procedure mmcomplex_v1s,mmcomplex_v1e, mmcomplex_m
  end interface mmcomplex

  interface ecacomplex  !the eca kinetics
     module procedure ecacomplex_v1s,ecacomplex_v1e, ecacomplex_m
  end interface ecacomplex

  interface ecacomplex_cell_norm  !the eca kinetics
     module procedure ecacomplex_cell_norm_v1s,ecacomplex_cell_norm_v1e, ecacomplex_cell_norm_m
  end interface ecacomplex_cell_norm


  public :: ECA_normflux
  public :: supeca_normflux
  public :: supecaijk
  public :: ecaij
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
   subroutine ecacomplex_m(kd,ss,ee,siej, bstatus)
     ! !DESCRIPTION:
     !compute concentrations of the enzyme substrate complexes
     !using the first order accurate ECA kinetics
     ! many substrate vs many enzymes
     use BetrStatusType  , only : betr_status_type
     implicit none
     ! !ARGUMENTS:
     real(r8), dimension(:,:), intent(in)  :: kd
     real(r8), dimension(:), intent(in)    :: ee, ss
     real(r8), dimension(:,:), intent(out) :: siej
     type(betr_status_type), intent(out)   :: bstatus

     ! !LOCAL VARIABLES:
     integer :: ii,jj
     integer :: i, j, k
     real(r8) :: dnm1, dnm2
     real(r8), dimension(:), allocatable :: Fr, Fc

     call bstatus%reset()
     ii = size(ss)       !number of substrates, dim 1
     jj = size(ee)       !number of enzymes, dim2
     if(ii/=size(siej,1) .or. jj/=size(siej,2))then
        call bstatus%set_msg(msg='wrong matrix shape in ecacomplex_m ' &
           //errMsg(mod_filename, __LINE__), err=-1)
        if(bstatus%check_status())return
     endif
     siej = 0._r8
     allocate(Fr(1:ii))
     allocate(Fc(1:jj))
     call ECA_normflux(ii, jj, kd, ss,ee, Fc,Fr)
     do j = 1, jj
       do i = 1, ii
         siej(i,j) = ecaij(i, j, ii, jj, kd, ss, ee, Fc, Fr)
       enddo
     enddo
     if(allocated(Fr))deallocate(Fr)
     if(allocated(Fc))deallocate(Fc)
   end subroutine ecacomplex_m

   !-------------------------------------------------------------------------------
   function ecaij(i, j, ii, jj, kd, ss, ee, Fc, Fr) result(cij)

   implicit none
   integer, intent(in)  :: i, j
   integer, intent(in)  :: ii, jj
   real(r8), intent(in) :: kd(1:ii,1:jj)
   real(r8), intent(in) :: ss(1:ii)
   real(r8), intent(in) :: ee(1:jj)
   real(r8), intent(in) :: Fc(1:jj)
   real(r8), intent(in) :: Fr(1:ii)

   real(r8) :: cij

   if(is_active_kd(kd(i,j)))then
     cij = ee(j) * ss(i)/kd(i,j)/(1._r8 + Fc(j) + Fr(i))
   else
     cij = 0._r8
   endif
   end function ecaij
   !-------------------------------------------------------------------------------
   subroutine ECA_normflux(ii, jj, kd, ss,ee, Fc,Fr)
   !
   !DESCRIPTION
   !compute the normalized ECA row and column fluxes
   implicit none
   !ARGUMENTS
   integer , intent(in)  :: ii, jj
   real(r8), intent(in)  :: kd(1:ii,1:jj)
   real(r8), intent(in)  :: ss(1:ii)
   real(r8), intent(in)  :: ee(1:jj)
   real(r8), intent(out) :: Fc(1:jj)
   real(r8), intent(out) :: Fr(1:ii)

   integer :: i, j
   Fc = 0._r8; Fr=0._r8
   do j = 1, jj
     do i = 1, ii
       if(is_active_kd(kd(i,j)))then
         Fc(j) = Fc(j) + ss(i)/kd(i,j)
         Fr(i) = Fr(i) + ee(j)/kd(i,j)
       endif
     enddo
   enddo

   end subroutine ECA_normflux
   !-------------------------------------------------------------------------------
   function is_active_kd(kd)result(ans)
   implicit none
   real(r8), intent(in) :: kd

   logical :: ans

   ans = kd >0._r8 .and. (kd<.9*kd_infty)
   end function is_active_kd
   !-------------------------------------------------------------------------------

   function supecaijk(i,j,k, ii,jj,kk, kdae, kdbe, aa, bb, ee, Fcak, Fcbk, &
      Frai, Frbj, Gaik, Gbjk) result(cijk)

   !compute the supeca complex
   implicit none
   integer, intent(in) :: i, j, k
   integer, intent(in) :: ii, jj, kk
   real(r8), intent(in) :: kdae(1:ii,1:kk)
   real(r8), intent(in) :: kdbe(1:jj,1:kk)
   real(r8), intent(in) :: aa(1:ii)
   real(r8), intent(in) :: bb(1:jj)
   real(r8), intent(in) :: ee(1:kk)
   real(r8), intent(in) :: Fcak(1:kk)
   real(r8), intent(in) :: Fcbk(1:kk)
   real(r8), intent(in) :: Frai(1:ii)
   real(r8), intent(in) :: Frbj(1:jj)
   real(r8), intent(in) :: Gaik(1:ii,1:kk)
   real(r8), intent(in) :: Gbjk(1:jj,1:kk)

   real(r8) :: cijk
   real(r8) :: Gabijk, Fcabk, denorm

   cijk = 0._r8
   if(is_active_kd(kdae(i,k)))then
      if(is_active_kd(kdbe(j,k)))then
        !double substrate
        !compute the denorminator
        Fcabk = Fcak(k) + Fcbk(k)
        Gabijk = Gaik(i,k) + Gbjk(j,k)
        denorm = Gaik(i,k) * Gbjk(j,k) * Fcabk / Gabijk + Fcabk - &
           (Fcak(k) * Gbjk(j,k) + Gaik(i,k)*Fcbk(k) - Gaik(i,k) * Gbjk(j,k))/Gabijk
        cijk =ee(k) * aa(i)/Kdae(i,k) * bb(j) / Kdbe(j,k) /denorm
      endif
   endif
   end function supecaijk

   !-------------------------------------------------------------------------------
   subroutine supeca_normflux(ii,jj, kk, kdae,kdbe, aa, bb, ee, Fcak, Fcbk, &
      Frai, Frbj, Gaik, Gbjk)
   !compute normalized supeca flux
   implicit none
   integer, intent(in)  :: ii
   integer, intent(in)  :: jj
   integer, intent(in)  :: kk
   real(r8), intent(in) :: kdae(1:ii,1:kk)
   real(r8), intent(in) :: kdbe(1:jj,1:kk)
   real(r8), intent(in) :: aa(1:ii)
   real(r8), intent(in) :: bb(1:jj)
   real(r8), intent(in) :: ee(1:kk)
   real(r8), intent(out):: Fcak(1:kk)
   real(r8), intent(out):: Fcbk(1:kk)
   real(r8), intent(out):: Frai(1:ii)
   real(r8), intent(out):: Frbj(1:jj)
   real(r8), intent(out):: Gaik(1:ii,1:kk)
   real(r8), intent(out):: Gbjk(1:jj,1:kk)

   integer :: i, j, k

   call  ECA_normflux(ii, kk, kdae, aa, ee, Fcak, Frai)

   call  ECA_normflux(jj, kk, kdbe, bb, ee, Fcbk, Frbj)

   do k = 1, kk
     do i = 1, ii
       Gaik(i,k) = Fcak(k) + Frai(i)
     enddo
     do j = 1, jj
       Gbjk(j,k) = Fcbk(k) + Frbj(j)
     enddo
   enddo

   end subroutine supeca_normflux
   !-------------------------------------------------------------------------------
   subroutine ecacomplex_cell_norm_m(kd,ss,ee,siej, bstatus)
     ! !DESCRIPTION:
     ! compute concentrations of the enzyme substrate complexes
     ! using the first order accurate ECA kinetics
     ! and noramlize the return value with cell abundance
     ! many substrates vs many enzymes
     use BetrStatusType  , only : betr_status_type
     implicit none
     ! !ARGUMENTS:
     real(r8), dimension(:,:), intent(in)  :: kd
     real(r8), dimension(:), intent(in)    :: ee, ss
     real(r8), dimension(:,:), intent(out) :: siej
     type(betr_status_type), intent(out)   :: bstatus
     ! !LOCAL VARIABLES:
     integer :: ii,jj
     integer :: i, j, k
     real(r8) :: dnm1, dnm2

     call bstatus%reset()
     ii = size(ss)       !number of substrates, dim 1
     jj = size(ee)       !number of enzymes, dim2
     if(ii/=size(siej,1) .or. jj/=size(siej,2))then
        call bstatus%set_msg(msg='wrong matrix shape in ecacomplex_m ' &
          //errMsg(mod_filename, __LINE__), err=-1)
        if(bstatus%check_status())return
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
     real(r8),               intent(in)  :: ee
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
