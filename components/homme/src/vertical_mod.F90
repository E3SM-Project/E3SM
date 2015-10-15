#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module vertical_mod
  ! ------------------------------
  use kinds, only : real_kind
  ! ------------------------------
  use dimensions_mod, only : nlev
  ! ------------------------------
implicit none
private

  type, public :: vcoord_t
     real (kind=real_kind) :: A(nlev+1)         ! hybrid pressure coordinate coefficients on interfaces
     real (kind=real_kind) :: B(nlev+1)         ! hybrid sigma coordinate (dp/deta) coefficient
     real (kind=real_kind) :: Amid(nlev)        ! hybrid pressure coordinate coefficients on mid levels
     real (kind=real_kind) :: Bmid(nlev)        ! hybrid sigma coordinate (dp/deta) coefficient on mid levels
     real (kind=real_kind) :: dB(nlev)          ! delta B on mid levels
     real (kind=real_kind) :: Hmat(nlev,nlev)   ! hydrostatic integral matrix        (sigma coordinates)
     real (kind=real_kind) :: Cmat(nlev,nlev)   ! energy convergence integral matrix (sigma coordinates)
  end type 

  public :: hmat_init
  public :: cmat_init

  public :: vcoord_init
  public :: vcoord_print

  public :: hmat_print

  private :: ecmwf_hmat_init
  private :: ecmwf_cmat_init

  private :: ccm_hmat_init
  private :: ccm_cmat_init


contains
  
  ! ====================================================
  ! subroutine ecmwf_Bmid encapsulates the
  ! ECMWF method for computing sigma mid levels.
  !
  ! reference: Ritchie,et.al. and Burridge and Simmons.
  ! ====================================================

  subroutine ecmwf_Bmid(sigma,dsigma,sigmam)

    real (kind=real_kind) :: sigma(nlev+1)  ! interface sigma level
    real (kind=real_kind) :: dsigma(nlev)   ! delta sigma level
    real (kind=real_kind) :: sigmam(nlev)   ! mid point sigma level

    ! Local variables

    integer  k
    real (kind=real_kind) :: C

    C = 1.0D0
    sigmam(1) = 0.5D0*dsigma(1)
    do k=2,nlev
       sigmam(k) = EXP((1.0D0/dsigma(k))*(sigma(k+1)*LOG(sigma(k+1))   &
                                        - sigma(k  )*LOG(sigma(k  )) ) - C)
    end do

  end subroutine ecmwf_Bmid

  ! ====================================================
  ! subroutine ccm1_sigma_mid encapsulates the
  ! CCM1 method for computing sigma mid levels.
  !
  ! reference: ???
  ! ====================================================

  subroutine ccm1_Bmid(sigma,dsigma,sigmam)

    real (kind=real_kind) :: sigma(nlev+1)  ! interface sigma level
    real (kind=real_kind) :: dsigma(nlev)   ! delta sigma level
    real (kind=real_kind) :: sigmam(nlev)   ! mid point sigma level
    integer  k

    do k=1,nlev
       sigmam(k) = 0.5D0*(sigma(k+1) + sigma(k))
    end do

  end subroutine ccm1_Bmid

!  ================================================
!  cmat_init:
!
!  ================================================

  function cmat_init(p,k,l,npts) result(Cmat)

    integer, intent(in) :: npts
    real(kind=real_kind), intent(in) :: p(npts,npts,nlev)
    integer, intent(in) :: k,l

    real(kind=real_kind) :: Cmat(npts,npts)
    integer              :: i,j

    if (l<k) then
       Cmat(:,:)=1.0D0/p(:,:,k)
    else if (l==k) then
       Cmat(:,:)=0.5D0/p(:,:,k)
    else
       Cmat(:,:)=0.0D0
    end if

  end function cmat_init



!  ================================================
!  hmat_init:
!
!  Compute the Hydrostatic Matrix Element H(k,l)
!  using the relationships in eq. (3.a.109) of the
!  CCM-2 Description (June 1993), NCAR/TN-382+STR
!
!  ================================================

  function hmat_init(logps,logp) result(Hmat)
    use dimensions_mod, only : np
    real(kind=real_kind), intent(in) :: logps(np,np)
    real(kind=real_kind), intent(in) :: logp(np,np,nlev)
    integer :: k,m

    real(kind=real_kind) :: Hmat(np,np,nlev,nlev)
    integer              :: i,j


    do k=2,nlev
!      if(m<k)
       do m=1,k-1
          hmat(:,:,m,k)=0.0_real_kind
       end do
    end do
    do k=1,nlev-1
!      if(m==k)
       do j=1,np
          do i=1,np
             Hmat(i,j,k,k)  = 0.50D0*(logp(i,j,k+1)-logp(i,j,k))
          end do
       end do
!      m>k<nlev       
       do m=k+1,nlev-1
          do j=1,np
             do i=1,np
                Hmat(i,j,m,k)  = 0.50D0*(logp(i,j,m+1)-logp(i,j,m-1))
             end do
          end do

       end do
!      m==nlev!=k
       do j=1,np
          do i=1,np
             Hmat(i,j,nlev,k) = logps(i,j) - 0.50D0*(logp(i,j,nlev) + logp(i,j,nlev-1))   
          end do
       end do
    end do
!   m==k==nlev    
    do j=1,np
       do i=1,np
          Hmat(i,j,nlev,nlev)  = logps(i,j) - logp(i,j,nlev)
       end do
    end do
 
   end function hmat_init

  ! ======================================================
  ! vcoord_print:
  ! 
  ! print hybrid coordinate values (A and B coefficients)
  ! ======================================================

  subroutine vcoord_print(vert,formulation)

    type (vcoord_t)  :: vert
    character(len=*) formulation   ! formulation choice (ecmwf,ccm1)

    integer k

    write(6,10)formulation
 10 format(a6," hybrid vertical coordinate formulation")

    write(6,*)"---------------------------------------"
    write(6,*)"interfaces"
    do k=1,nlev+1
       write(6,11)k-1,vert%A(k),k-1,vert%B(k)
 11    format("A(",i4,"+1/2)=",e22.15," B(",i4,"+1/2)=",e22.15)
    end do

    write(6,*)""
    write(6,*)"---------------------------------------"
    write(6,*)"mid levels:"
    do k=1,nlev
       write(6,12)k,vert%Amid(k),k-1,vert%Bmid(k)
 12    format("A(",i4,")=",e22.15," B(",i4,")=",e22.15)
    end do
    write(6,*)""

  end subroutine vcoord_print

  ! =============================================================
  ! vcoord_init: Initialize vertical coordinate system...
  ! returns: 0 for success
  ! =============================================================

  function vcoord_init(vert, formulation, vfile_mid, vfile_int) result(ierr)

    integer :: ierr

    type (vcoord_t)  :: vert

    character(len=*) formulation   ! formulation choice (ecmwf,ccm)
    character(len=*) vfile_mid     ! mid levels vertical coordinate system 
    character(len=*) vfile_int     ! interface levels vertical coordinate system 
    
    ! Local variables

    real (kind=real_kind) :: s(nlev+1)
    integer ierr11,ierr12

    integer k

    ierr=0

!$OMP CRITICAL
    open(UNIT=11,FILE=vfile_int, form='unformatted',status='old',iostat=ierr11)

    if (ierr11 /= 0) then
       ierr=ierr11
    end if

    open(UNIT=12,FILE=vfile_mid, form='unformatted',status='old',iostat=ierr12)

    if (ierr12 /= 0) then
       ierr=ierr12
    end if

    do k=1,nlev+1
       read(11)vert%A(k)
       read(11)vert%B(k)
    end do

    do k=1,nlev
       read(12)vert%Amid(k)
       read(12)vert%Bmid(k)
    end do

    close(11)
    close(12)
!$OMP END CRITICAL

    print *,"successfully read vertical coordinates"

    do k=1,nlev
       vert%dB(k)   = vert%B(k+1) - vert%B(k)
    end do

  end function vcoord_init

  ! ===============================================
  ! hmat_print prints the entries of the 
  ! CCM1 Hydrostatic matrix (Hmat) in the case of 
  ! sigma coordinates.
  ! ===============================================

  subroutine hmat_print(vert)

    type (vcoord_t) :: vert

    ! Local variables
  
    integer k,l

    print *
    print *,"CCM Hydrostatic Matrix"
    do k=1,nlev
       do l=1,nlev
          write(6,10)k,l,vert%Hmat(k,l)
 10       format("B(",i4,",",i4,")=",e22.15)
       end do
       print *
    end do

  end subroutine hmat_print

  ! ===============================================
  ! ccm_hmat_init fills in entries of CCM1 H matrix.
  ! This applies only to the sigma coordinate
  ! implementation of CCM-1.
  !
  ! see Description of CCM1 (NCAR/TN-285+STR), 
  ! page 40-41 for details.
  ! ===============================================

  subroutine ccm_hmat_init(vert)

    type (vcoord_t) :: vert

    ! Local variables
  
    integer k,l

    ! =================================== 
    ! eq (3.b.14) of Description of CCM1 
    !                (NCAR/TN-285+STR)
    ! ===================================

    vert%Hmat(nlev,nlev)= -LOG(vert%Bmid(nlev))    

    ! ===================================
    ! eq (3.b.15), ibid.
    ! ===================================

    do k=1,nlev-1             
       vert%Hmat(nlev,k)=0.0D0
    end do
    ! ===================================
    ! eq (3.b.16), ibid.
    ! ===================================

    do k=1,nlev-1
       vert%Hmat(k,k)=0.50D0*LOG(vert%Bmid(k+1)/vert%Bmid(k))
    end do
   
    ! ===================================
    ! eq (3.b.17), ibid.
    ! ===================================

    do k=1,nlev-1
       vert%Hmat(k,k+1)= vert%Hmat(k+1,k+1) + 0.50D0*LOG(vert%Bmid(k+1)/vert%Bmid(k))
    end do

    ! ===================================
    ! eq (3.b.18), ibid.
    ! ===================================

    do l=3,nlev
       do k=l-2,1,-1
          vert%Hmat(k,l)=vert%Hmat(k+1,l)
       end do
    end do

    ! ===================================
    ! eq (3.b.19), ibid.
    ! ===================================

    do k=1,nlev
       do l=1,k-1
          vert%Hmat(k,l)=0.0D0
       end do
    end do

  end subroutine ccm_hmat_init

  ! ===============================================
  ! ccm1_cmat_init fills in entries of CCM1 C matrix.
  ! see Description of CCM1 (NCAR/TN-285+STR), 
  ! page 41, eq. (3.b.21) for details.
  !
  ! ===============================================

  subroutine ccm_cmat_init(vert)

    type (vcoord_t) :: vert

    ! Local variables

    integer j,k
     
    do k=1,nlev
       do j=1,nlev
          vert%Cmat(k,j) = vert%Hmat(j,k)*(vert%dB(j)/vert%dB(k))
       end do
    end do

  end subroutine ccm_cmat_init

  ! ===============================================
  ! ecmwf_hmat_init fills in entries of ECMWF hydrostatic
  ! integral or B matrix.
  !  (see RDL nodes deriving this from
  !   eqs. (2.21) and (2.22) of Ritchie, et. al.)
  ! ===============================================

  subroutine ecmwf_hmat_init(vert)

    type (vcoord_t) :: vert

    ! Local variables
  
    integer j,k

    do k=1,nlev
       do j=1,k-1
          vert%Hmat(k,j)=0.0D0
       end do

       vert%Hmat(k,k)=1.0D0 - (vert%B(k)/vert%dB(k))*LOG(vert%B(k+1)/vert%B(k))

       do j=k+1,nlev
          vert%Hmat(k,j)=LOG(vert%B(j+1)/vert%B(j))
       end do
    end do        
     
    ! ======================================================
    ! Special case (see Ritchie discussion of alpha(1))
    ! ======================================================

    vert%Hmat(1,1) = LOG(2.0D0)

  end subroutine ecmwf_hmat_init

  ! ===============================================
  ! ecmwf_cmat_init fills in entries of ECMWF energy
  ! conversion integral, or C matrix.
  !  (see RDL nodes deriving this from
  !   eqs. (2.25) of Ritchie, et. al.)
  ! ===============================================

  subroutine ecmwf_cmat_init(vert)

    type (vcoord_t) :: vert

    ! Local variables
  
    integer j,k

    do k=1,nlev

       do j=1,k-1
          vert%Cmat(k,j)= (vert%dB(j)/vert%dB(k))*LOG(vert%B(k+1)/vert%B(k))
       end do

       vert%Cmat(k,k)=1.0D0 - (vert%B(k)/vert%dB(k))*LOG(vert%B(k+1)/vert%B(k))

       do j=k+1,nlev
          vert%Cmat(k,j)=0.0D0
       end do
    end do        
     
    ! ======================================================
    ! Special case (see Ritchie discussion of alpha(1))
    ! ======================================================

    vert%Cmat(1,1) = LOG(2.0D0)

  end subroutine ecmwf_cmat_init

end module vertical_mod
