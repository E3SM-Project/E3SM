#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!MODULE FVM_RECONSTRUCTION_MOD--------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! This module contains everything  to do (ONLY) a CUBIC (3rd order) reconstruction  ! 
!                                                                                   !
! IMPORTANT: the implementation is done for a ncfl > 1, which is not working        !
!            but it works for ncfl=1                                                !
!
! This module has been recoded for multi-tracer efficiency (May, 2014)
!
!-----------------------------------------------------------------------------------!
module fvm_reconstruction_mod

  use kinds, only                  : int_kind, real_kind
  use dimensions_mod, only         : nc,nhe,nhr,nht,ns,nhc
  use coordinate_systems_mod, only : cartesian2D_t,cartesian3D_t
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  use fvm_control_volume_mod, only: fvm_struct
  use parallel_mod, only: abortmp
  implicit none
  private
  integer, parameter, private:: nh = nhr+(nhe-1) ! = 2 (nhr=2; nhe=1)
                                                 ! = 3 (nhr=2; nhe=2)
  logical, private, parameter :: lplot = .false.
  public :: reconstruction
contains
! ----------------------------------------------------------------------------------!
!SUBROUTINE RECONSTRUCTION------------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: controls the cubic (3rd order) reconstructions:                      !
!                                                                                   !
! CALLS: fillhalo_cubic, reconstruction_cubic                                       !
! INPUT: fcube    ...  tracer values incl. the halo zone                            !
!        fvm    ...  structure incl. tracer values aso                            !                                   ! 
! OUTPUT:recons   ...  has the reconstruction coefficients (5) for the 3rd order    !
!                      reconstruction: dx, dy, dx^2, dy^2, dxdy                     !
!-----------------------------------------------------------------------------------!
subroutine reconstruction(fcube,fvm,recons)
  use fvm_control_volume_mod, only: fvm_struct
  implicit none
  type (fvm_struct), intent(in)                                  :: fvm
  !
  ! dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc)
  !
  real (kind=real_kind), dimension(1-nhc:, 1-nhc:), intent(in)   :: fcube
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe), intent(out)  :: recons 
                                                    

  real (kind=real_kind), dimension(1-nht:nc+nht, 1-nht:nc+nht)    :: fpanel
  real   (kind=real_kind), dimension(1-nht:nc+nht,1-nht:nc+nht,2) :: fotherpanel

  call fill_halo(fcube,fvm,fpanel,fotherpanel)
  call get_reconstruction(fvm, fpanel, fotherpanel,recons)
  
!  call debug_halo(fvm,fcubenew,fpanel)
!  call debug_halo_recons(fvm,recons,recons_trunk)
!  call print_which_case(fvm)
!  
!  call debug_halo_neighbor       (fvm,fotherface,fotherpanel)
!  call debug_halo_neighbor_recons(fvm,recons,recons_trunk)
end subroutine reconstruction
!END SUBROUTINE RECONSTRUCTION--------------------------------------------CE-for FVM!

function matmul_w(w,f)
  implicit none
  real (kind=real_kind)                          :: matmul_w
  real (kind=real_kind),dimension(:), intent(in) :: w,f      !dimension(ns)
  real (kind=real_kind)                          :: tmp
  integer                                        :: k
  matmul_w = 0.0D0
  do k=1,ns
     matmul_w = matmul_w+w(k)*f(k)
  end do
end function matmul_w

subroutine fill_halo(fcube,fvm,fpanel,fotherpanel)!dbg
  implicit none
  real (kind=real_kind),   &
       dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc), intent(in)          :: fcube
  type (fvm_struct), intent(in)                                   :: fvm
        

!  real (kind=real_kind)  , dimension(1-nht:nc+nht, 1-nht:nc+nht ), intent(out) :: fpanel
!  real   (kind=real_kind), dimension(1-nht:nc+nht,1-nht:nc+nht,2), intent(out) :: fotherpanel
  real (kind=real_kind)  , dimension(1-nht:, 1-nht:),  intent(out) :: fpanel
  real   (kind=real_kind), dimension(1-nht:,1-nht:,:), intent(out) :: fotherpanel


  integer (kind=int_kind)                                         :: i, j, halo,ibaseref
  real (kind=real_kind), dimension(1-nh:nc+nh,1:nhr,1:ns) :: w
  !
  fpanel = 1.0E19 !dbg
  !
  ! 
  ! Stencil for reconstruction is:
  !
  !     ---------------------
  !     |   |   | i |   |   |
  !     ---------------------
  !     |   | i | i | i |   |
  !     ---------------------
  !     | i | i | R | i | i |
  !     ---------------------
  !     |   | i | i | i |   |
  !     ---------------------
  !     |   |   | i |   |   |
  !     ---------------------
  !
  ! where
  !
  !   "R" is cell for which we whish to do the reconstruction
  !   "i" is the stencil 
  !
  !
  ! If one or more point in the stencil is on another panel(s) then we need to interpolate
  ! to a stencil that is an extension of the panel on which R is located
  ! (this is done using one dimensional cubic Lagrange interpolation along panel side)
  !
  ! Example: say that southern most "s" on Figure above is on another panels projection then the stencil becomes
  !
  !
  !   ---------------------------------
  !   |   |   |   |   |   | i |   |   |
  !   ----------------|----------------
  !   |   |   |   |   | i | i | i |   |
  !   ----------------|----------------
  !   |   |   |   | i | i | R | i | i |
  !   ----------------|----------------
  !   |   |   |   |   | i | i | i |   |
  !   ---------------------------------
  !   /   /   /   /   / S /S&i/ S / S /
  !  /---/---/---/---/---/---/---/---/
  ! /   /   /   /   /   /   /   /   /
  !/---/---/---/---/---/---/---/---/
  !
  !
  ! where "S" are the cell average values used for the cubic interpolation (located on the South panel)
  !
  !
  if (fvm%cubeboundary==0) then 
     fpanel(1-nht:nc+nht,1-nht:nc+nht)=fcube(1-nht:nc+nht,1-nht:nc+nht)
  else if (fvm%cubeboundary==west) then
     !                                                       !
     !                                                       ! Case shown below: nhr=2, nhe=1, nht=nhr+nhe
     !                                                       ! (nhr = reconstruction width along x and y)
     !                                                       ! (nhe = max. Courant number)
     !                                                       !  
     !                                                       !
     !    Figure below shows the element in question         ! In terms of data structure:
     !    (center element) and the surrounding elements      !     
     !    on the element in question's projection            !     * "H" is on same panel average value
     !                                                       !     * "w" is west panel average values that need 
     !    Notation: "0" marks the element boundaries         !       to be interpolated to main element 
     !                                                       !       projection
     !    Elements to the west are on a different projection !     * "i" is extra halo required by the cubic
     !                                                       !       interpolation
     !     0                                                 !
     !     |0000                                             !
     !     |   |00000                                        !
     !     |\--|   |000000000000000000000000000000000000     !    -x---x---x---x---x---x---x---x---x---x---x---x
     !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   |   | i |   |   |   |   |   |   |   |   |
     !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x 
     !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | i | H | H | H | H | H |   |   |   |
     !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x 
     !     0   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | w | H | H | H | H | H | H |   |   |
     !     |0000   |\--0---------------0---------------0     !    -------------x---------------x---------------x
     !     |   |0000   0   |   |   |   0   |   |   |   0     !     |   | w | w | r | r | r | r | r | H | H |   |
     !     |\--|   |000000000000000000000000000000000000     !    -x---x---x---00000000000000000---x---x---x---x
     !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
     !     |\--|   \---0---------------0---------------0     !    -------------0---------------0---------------x
     !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
     !     |\--|   \---0---------------0---------------0     !    -------------0---------------0---------------x
     !     0   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
     !     |0000   |\--0---------------0---------------0     !    -------------0---------------0---------------x
     !     |   |0000   0   |   |   |   0   |   |   |   0     !     |   | w | w 0 r | r | r | r 0 r | H | H |   |
     !     |\--|   |000000000000000000000000000000000000     !    -x---x---x---00000000000000000---x---x---x---x
     !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | w | w | r | r | r | r | r | H | H |   |
     !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x
     !     |   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | w | H | H | H | H | H | H |   |   |
     !     |\--|   |\--0---------------0---------------0     !    -------------x---------------x---------------x
     !     0   |\--|   0   |   |   |   0   |   |   |   0     !     |   | i | i | H | H | H | H | H |   |   |   |
     !      0000   |\--0---------------0---------------0     !    -------------x---------------x---------------x
     !          0000   0   |   |   |   0   |   |   |   0     !     |   |   | i |   |   |   |   |   |   |   |   |
     !              000000000000000000000000000000000000     !    -x---x---x---x---x---x---x---x---x---x---x---x
     !  
     !      
     !      -2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8       
     !
     !
     !
     ! fill in values (incl. halo) that are on the "main" panels projection
     ! 
     fpanel(1:nc+nht,1-nht:nc+nht)=fcube(1:nc+nht,1-nht:nc+nht)
     !
     ! fill in values that are on the west panels projection
     ! 
     fotherpanel (1-nht:0,1-nht:nc+nht,1)=fcube(1-nht:0,1-nht:nc+nht)
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
           ibaseref=fvm%ibase(i,halo,1)                                              
!           ibaseref = ibase(i,halo,1)
           fpanel(1-halo ,i) = matmul_w(w(i,halo,:),fcube(1-halo ,ibaseref:ibaseref+ns-1))
           !
           ! Exploit symmetry in interpolation weights
           !
           fotherpanel(halo,i,1)     = matmul_w(w(i,halo,:),fcube(halo   ,ibaseref:ibaseref+ns-1))
        end do
     end do
  else if (fvm%cubeboundary==east) then
     !
     ! north part is on different panel
     !
     ! stencil
     !
     ! CN<1 case                                             !
     !                                                       !
     !                                                  
     !                                                  
     !                                                 0     !                                      
     !                                             0000|     !                                      
     !                                         0000|   |     !                                      
     !     000000000000000000000000000000000000|   |--/|     !     x---x---x---x---x---x---x---x---x---x---x---x-
     !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   |   |   |   |   |   |   |   | i |   |   |
     !     0---------------0---------------0--/    |--/|     !     x---------------x---------------x---x---x---x-
     !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   |   |   | H | H | H | H | H | i | i |   |
     !     0---------------0---------------0--/|   |--/|     !     x---------------x---------------x---x---x---x-
     !     0   |   |   |   0   |   |   |   0   |--/|   0     !     |   |   | H | H | H | H | H | H | e | i |   |
     !     0---------------0---------------0--/|   0000|     !     x---------------x---------------x---x---x---x-
     !     0   |   |   |   0   |   |   |   0   0000|   |     !     |   | H | H | r | r | r | r | r | e | e |   |
     !     000000000000000000000000000000000000|   |--/|     !     x---x---x---x---00000000000000000---x---x---x-
     !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
     !     0---------------0---------------0--/|   |--/|     !     x---------------0---------------0---x---x---x-
     !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
     !     0---------------0---------------0--/|   |--/|     !     x---------------0---------------0---x---x---x-
     !     0   |   |   |   0   |   |   |   0   |--/|   0     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
     !     0---------------0---------------0--/|   0000|     !     x---------------0---------------0---x---x---x-
     !     0   |   |   |   0   |   |   |   0   0000|   |     !     |   | H | H | r 0 r | r | r | r 0 e | e |   |
     !     000000000000000000000000000000000000|   |--/|     !     x---x---x---x---00000000000000000---x---x---x-
     !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   | H | H | r | r | r | r | r | e | e |   |
     !     0---------------0---------------0--/|   |--/|     !     ----------------x---------------x---x---x---x-
     !     0   |   |   |   0   |   |   |   0   |--/|   |     !     |   |   | H | H | H | H | H | H | e | i |   |
     !     0---------------0---------------0--/|   |--/|     !     ----------------x---------------x---x---x---x-
     !     0   |   |   |   0   |   |   |   0   |--/|   0     !     |   |   |   | H | H | H | H | H | i | i |   |
     !     0---------------0---------------0--/|   0000      !     ----------------x---------------x---x---x---x-
     !     0   |   |   |   0   |   |   |   0   0000          !     |   |   |   |   |   |   |   |   | i |   |   |
     !     000000000000000000000000000000000000              !     x---x---x---x---x---x---x---x---x---x---x---x-
     !  
     !      
     !      -3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 
     !
     fpanel      (1-nht:nc     ,1-nht:nc+nht  )=fcube(1-nht:nc     ,1-nht:nc+nht)
     fotherpanel (nc+1 :nc+nht ,1-nht:nc+nht,1)=fcube(nc+1 :nc+nht ,1-nht:nc+nht) !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
!           ibaseref=fvm%ibase(i,halo,1 )                                              
           ibaseref = fvm%ibase(i,halo,1)
           fpanel      (nc+halo   ,i  ) = matmul_w(w(i,halo,:),fcube(nc  +halo,ibaseref:ibaseref+ns-1))
           fotherpanel (nc+1-halo ,i,1) = matmul_w(w(i,halo,:),fcube(nc+1-halo,ibaseref:ibaseref+ns-1))
        end do
     end do
  else if (fvm%cubeboundary==north) then
     !
     ! north part is on different panel
     !
     ! stencil
     !
     ! CN<1 case                                            
     !                                                      !   x---------------x---------------x---------------x
     !                                                      !   |   |   |   |   |   |   |   |   |   |   |   |   |
     !0---\---\---\---0---\---\---\---0---\---\---\---0     !   x---------------x---------------x---------------x
     ! 0   \   \   \   0   \   \   \   0   \   \   \   0    !   |   | i | i | n | n | n | n | n | n | i | i |   |
     !  0---\---\---\---0---\---\---\---0---\---\---\---0   !   x---------------x---------------x---------------x
     !   0   \   \   \   0   \   \   \   0   \   \   \   0  !   | i | i | n | n | n | n | n | n | n | n | i | i |
     !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
     !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
     !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
     !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
     !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r | r | r | r | r | r | H | H |   |
     !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   | H | H | H | H | H | H | H | H |   |   |
     !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   | H | H | H | H | H | H |   |   |   |
     !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   |   |   |   |   |   |   |   |   |   |
     !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---x---x---x---x---x---x---x---x---x
     !        
     !    -3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8    !    -3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 
     !
     ! fill in values that are on the same projection as "main" element
     fpanel      (1-nht:nc+nht ,1-nht:nc)=fcube(1-nht:nc+nht ,1-nht:nc)
     ! fill in halo for north element
     fotherpanel (1-nht:nc+nht ,nc+1:nc+nht,1)=fcube(1-nht:nc+nht ,nc+1:nc+nht)
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
           ibaseref = fvm%ibase(i,halo,1)
           fpanel      (i,nc+halo    ) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,nc+halo  )) !north
           fotherpanel (i,nc+1-halo,1) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,nc+1-halo))
        end do
     end do

  else if (fvm%cubeboundary==south) then
     !
     ! south part is on different panel
     !
     ! stencil
     !
     !                                                      !  
     !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---x---x---x---x---x---x---x---x---x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   |   |   |   |   |   |   |   |   |   |
     !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   |   | H | H | H | H | H | H |   |   |   |
     !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   |   | H | H | H | H | H | H | H | H |   |   |
     !   0---------------0---------------0---------------0  !   x---------------x---------------x---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r | r | r | r | r | r | H | H |   |
     !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
     !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
     !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
     !   0---------------0---------------0---------------0  !   x---------------0---------------0---------------x
     !   0   |   |   |   0   |   |   |   0   |   |   |   0  !   |   | H | H | r 0 r | r | r | r 0 r | H | H |   |
     !   0000000000000000000000000000000000000000000000000  !   x---x---x---x---00000000000000000---x---x---x---x
     !   0   /   /   /   0   /   /   /   0   /   /   /   0  !   | i | i | s | s | s | s | s | s | s | s | i | i |
     !  0---/---/---/---0---/---/---/---0---/---/---/---0   !   x---------------x---------------x---------------x
     ! 0   /   /   /   0   /   /   /   0   /   /   /   0    !   |   | i | i | s | s | s | s | s | s | i | i |   |
     !0---/---/---/---0---/---/---/---0---/---/---/---0     !   x---------------x---------------x---------------x
     !                                                      !   |   |   |   |   |   |   |   |   |   |   |   |   |
     !
     !     0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9           !     0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
     !
     ! fill in values that are on the same projection as "main" element (marked with "i" in Figure above)
     !
     fpanel      (1-nht:nc+nht,1:nc+nht  )=fcube(1-nht:nc+nht,1:nc+nht)
     fotherpanel (1-nht:nc+nht,1-nht:0 ,1)=fcube(1-nht:nc+nht,1-nht:0 )
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
           ibaseref=fvm%ibase(i,halo,1)!fvm%ibase(i,halo,2) 
           fpanel      (i,1-halo ) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,1-halo))  !south
           fotherpanel (i,  halo,1) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,  halo))
        end do
     end do
  else if (fvm%cubeboundary==swest) then
     !
     ! south and west neighboring cells are on different panel
     !
     ! stencil
     !
     !
     ! CN<1 case                                                
     !
     !
     !
     !     |000000000000000000000000000000000000   !   x---x---x---x---x---x---x---x---x---x---x---x---x
     !  0000   0   |   |   |   0   |   |   |   0   !   |   |   |   |   |   |   |   |   |   |   |   |   |
     ! 0   |/--0---------------0---------------0   !   x---------------x---------------x---------------x
     ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   |   | w | H | H | H | H | H |   |   |   |
     ! |   |/--0---------------0---------------0   !   x---------------x---------------x---------------x
     ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w | H | H | H | H | H | H |   |   |
     ! |   |/--0---------------0---------------0   !   x---------------x---------------x---------------x
     ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w | r | r | r | r | r | H | H |   |
     ! |   |000000000000000000000000000000000000   !   x---x---x---x---00000000000000000---x---x---x---x
     ! |0000   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
     ! 0   |/--0---------------0---------------0   !   x---------------0---------------0---------------x
     ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
     ! |   |/--0---------------0---------------0   !   x---------------0---------------0---------------x
     ! |/--|   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
     ! |   |/--0---------------0---------------0   !   x---------------0---------------0---------------x
     ! |   |   0   |   |   |   0   |   |   |   0   !   |   |   | w | w 0 r | r | r | r 0 r | H | H |   |
     ! | -/|   000000000000000000000000000000000   !   x---x---x---x---00000000000000000---x---x---x---x
     ! |/  | 0    /   /   /   0   /   /   /   0    !   |   |   | w |   | s | s | s | s | s | s |   |   |
     ! |   0-----/---/---/---0---/---/---/---0     !   x---------------x---------------x---------------x
     ! | 0      /   /   /   0   /   /   /   0      !   |   |   |   | s | s | s | s | s | s |   |   |   |
     ! 0-------/---/---/---0---/---/---/---0       !   x---------------x---------------x---------------x
     !
     !
     !  -1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |       !   |-3 |-2 |-1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
     !
     ! fill in values that are on the same projection as "main" element (marked with "i" in Figure above)
     !
     fpanel(1:nc+nht,1:nc+nht)=fcube(1:nc+nht,1:nc+nht)
     !
     ! fill in west part (marked with "w" on Figure above) and south part (marked with "s")
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=max(halo-nh,0),nc+nh-(halo-1)
           ibaseref=fvm%ibase(i,halo,1)!fvm%ibase(i,halo,1)      
           fpanel(1-halo ,i) = matmul_w(w(i,halo,:),fcube(1-halo ,ibaseref:ibaseref+ns-1)) !west
           fpanel(i,1-halo ) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,1-halo))  !south
        end do
     end do
     !
     ! corner value
     !
     fpanel(0,0)=0.25D0*(fpanel(0,1)+fpanel(1,0)+fpanel(-1,0)+fpanel(0,-1))
     !
     ! ****************************************************************
     !
     ! fill halo for reconstruction on south neighbor panel projection
     !
     ! ****************************************************************
     !
     ! On the south panel projection the neighbors are arragened as follow (nwest case):
     !
     !
     ! \
     !  \    p
     !   \  
     !    \-----
     !    |
     !  w |  s
     !    |  
     !
     !
     ! x---x---x---x---00000000000000000---x---x---x---x
     ! |   |   |   |   0   |   |   |   0   |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   |   |   0   |   |   |   0   |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   |   | p 0 p | p | p | p 0 p |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   | w | wp0 p | p | p | p 0 p | p |   |   |
     ! x---x---x---x---00000000000000000---x---x---x---x
     ! |   |   | w | w | r | r | r | r | r | i | i |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   | w | i | i | i | i | i | i |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   | i | i | i | i | i |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---  
     !
     !
     ! fill values on same panel projection ("r" and "i" on Figure above)
     !
     fotherpanel(1:nc+nht,1-nht:0,1)  = fcube(1:nc+nht,1-nht:0)
     !
     ! compute interpolated cell average values in "p" cells on Figure on above
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=max(halo-nh,0),nc+nh-(halo-1)
           ibaseref=fvm%ibase(i,halo,1)
           !
           ! use same weights as interpolation south from main panel (symmetric)
           !
           fotherpanel(i,halo,1)  = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,halo)) 
        end do
     end do
     !
     ! compute interpolated cell average values in "w" cells on Figure on above
     !
     w = fvm%halo_interp_weight(:,:,:,2)
     do halo=1,nhr
        do i=nc+halo-nhr,nc+1
           ibaseref=fvm%ibase(i,halo,2)-nc
           !
           ! fotherpanel indexing follows main panel indexing
           ! fcube indexing most be "rotated":
           !
           ! ===============================
           ! |              |              |
           ! |  W      ^    |   S          |    
           ! |         |    |              |
           ! |       x |    |              |
           ! |         |    |              |
           ! !              |              |
           ! !   <-----     |              |
           ! !      y       |              |
           ! !              |              |
           ! ===============================
           !
           fotherpanel(1-halo,i-nc,1)  = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,halo)) 
        end do
     end do
     fotherpanel(0,1,1) = 0.25D0*(fotherpanel(-1,1,1)+fotherpanel(1,1,1)+fotherpanel(0,2,1)+fotherpanel(0,0,1))
     !
     ! ****************************************************************
     !
     ! fill halo for reconstruction on west neighbor panel projection
     !
     ! ****************************************************************
     !
     ! On the west panel projection the neighbors are arragened as follow (seast case):
     !
     !   --------
     !   |      |
     !   |  w   |    p
     !   |      |
     !   -------\
     !           \
     !       s    \
     !       
     !
     !
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   | i |   |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   | i | i | e |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   | i | i | r | e | e |   |   |   |   |   |   |
     ! x---x---x---x---00000000000000000---x---x---x---x
     ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   | i | i | r 0 e | e |   |   0   |   |   |   |
     ! x---x---x---x---00000000000000000---x---x---x---x
     ! |   |   | s | s | se| e |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   | s | s |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---   
     !
     !
     ! fill values on same panel projection ("r" and "i" on Figure above)
     !
     fotherpanel(1-nht:nc,1:nc+nht,2)  = fcube(1-nht:nc,1:nc+nht)
     !
     ! compute interpolated cell average values in "p" cells on Figure on above
     !
     w = fvm%halo_interp_weight(:,:,:,1) ! symmetry
     do halo=1,nhr
        do i=max(halo-nh,0),nc+nh-(halo-1)
           ibaseref=fvm%ibase(i,halo,1)
           !
           ! use same weights as interpolation south from main panel (symmetric)
           !
           fotherpanel(halo,i,2)  = matmul_w(w(i,halo,:),fcube(halo,ibaseref:ibaseref+ns-1)) 
        end do
     end do
     !
     ! compute interpolated cell average values in "s" cells on Figure on above
     !
     w = fvm%halo_interp_weight(:,:,:,2)
     do halo=1,nhr
        do i=nc+halo-nhr,nc+1
           ibaseref=fvm%ibase(i,halo,2)-nc
           !
           ! fotherpanel indexing follows main panel indexing
           ! fcube indexing most be "rotated":
           !
           ! ===============================
           ! |              |              |
           ! |  W      ^    |   S          |    
           ! |         |    |              |
           ! |       x |    |              |
           ! |         |    |              |
           ! !              |              |
           ! !   <-----     |              |
           ! !      y       |              |
           ! !              |              |
           ! ===============================
           !
           fotherpanel(i-nc,1-halo,2)  = matmul_w(w(i,halo,:),fcube(halo,ibaseref:ibaseref+ns-1)) 
        end do
     end do
     fotherpanel(1,0,2) = 0.25D0*(fotherpanel(0,0,2)+fotherpanel(2,0,2)+fotherpanel(1,-1,2)+fotherpanel(1,1,2))
  else if (fvm%cubeboundary==seast) then
     !
     ! south and east neighboring cells are on different panel
     !
     !
     !
     ! 000000000000000000000000000000000000|       
     ! 0   |   |   |   0   |   |   |   0   0000    !   |   |   |   |   |   |   |   |   |   |   |   |   |
     ! 0---------------0---------------0--\|   0   !   x---------------x---------------x---------------x   
     ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   |   |   |   | H | H | H | H |   |   |   |   |
     ! 0---------------0---------------0--\|   |   !   x---------------x---------------x---------------x   
     ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   |   |   | H | H | H | H | H | e |   |   |   |
     ! 0---------------0---------------0--\|   |   !   x---------------x---------------x---------------x   
     ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   | H | H | r | r | r | r | r | e | e |   |   |
     ! 000000000000000000000000000000000000|   |   !   x---x---x---x---00000000000000000---x---x---x---x
     ! 0   |   |   |   0   |   |   |   0   0000|   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
     ! 0---------------0---------------0--\|   0   !   x---------------0---------------0---------------x   
     ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
     ! 0---------------0---------------0--\|   |   !   x---------------0---------------0---------------x   
     ! 0   |   |   |   0   |   |   |   0   |--\|   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
     ! 0---------------0---------------0--\|   |   !   x---------------0---------------0---------------x   
     ! 0   |   |   |   0   |   |   |   0   |   |   !   |   | H | H | r 0 r | r | r | r 0 e | e |   |   |
     ! 000000000000000000000000000000000   |\- |   !   x---x---x---x---00000000000000000---x---x---x---x
     !  0   \   \   \   0   \   \   \    0 |  \|   !   |   |   | s | s | s | s | s | s |s/e| e |   |   |
     !   0---\---\---\---0---\---\---\-----0   |   !   x---------------x---------------x---------------x   
     !    0   \   \   \   0   \   \   \      0 |   !   |   |   |   | s | s | s | s | s | s |   |   |   |
     !     0---\---\---\---0---\---\---\-------0   !   x---------------x---------------x---------------x   
     !
     !
     fpanel       (1-nht:nc,1:nc+nht)=fcube(1-nht:nc,1:nc+nht)
     !
     ! east
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=max(halo-nh,0),nc+nh-(halo-1)
           ibaseref = fvm%ibase(i,halo,1)
           fpanel(nc+halo,i) = matmul_w(w(i,halo,:),fcube(nc  +halo,ibaseref:ibaseref+ns-1))
        end do
     end do
     !
     ! south
     !
     w = fvm%halo_interp_weight(:,:,:,2)
     do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref = fvm%ibase(i,halo,2)
           fpanel(i,1-halo ) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,1-halo))  !south
        end do
     end do
     fpanel(nc+1,0   )=0.25D0*(&
          fpanel(nc+1,1)+fpanel(nc,0)+fpanel(nc+2,0)+fpanel(nc+1,-1)) 
     !
     ! ****************************************************************
     !
     ! fill halo for reconstruction on south neighbor panel projection
     !
     ! ****************************************************************
     !
     ! On the south panel projection the neighbors are arragened as follow (neast case):
     !
     !
     !             /
     !       P    /
     !           /
     !    ------/
     !    |     |  E
     !    |  S  | 
     !    |     |
     ! 
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   |   |   |   |   |   |
     ! x---x---x---x---00000000000000000---x---x---x---x
     ! |   |   |   |   0   |   |   |   0   |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   |   |   0   |   |   |   0   |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   |   | n 0 n | n | n | n 0 n |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   | n | n 0 n | n | n | n 0 ne| e |   |   |
     ! x---x---x---x---00000000000000000---x---x---x---x
     ! |   | i | i | r | r | r | r | r | e | e |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   | i | i | i | i | i | i | e |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   | i | i | i | i | i |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x   
     !     
     !
     !
     ! fill values on same panel projection ("r" and "i" on Figure above)
     !
     fotherpanel(1-nht:nc,1-nht:0,1)  = fcube(1-nht:nc,1-nht:0)
     !
     w = fvm%halo_interp_weight(:,:,:,2)
     !
     ! fill in "n" on Figure above
     !
     do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref = fvm%ibase(i,halo,2)
           fotherpanel (i,halo,1) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,  halo))
        end do
     end do
     !
     ! fill in "e" on Figure above
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=0,nht-halo!nc+nh-(halo-1)
           ibaseref = fvm%ibase(i,halo,1)
           !
           ! fother panel follows indexing on main panel
           !
           ! use symmetry for weights (same weights as East from main panel but for south panel
           ! projection the indecies are rotated)
           !
           fotherpanel (nc+halo ,1-i,1) = matmul_w(w(i,halo,:),fcube(nc+ibaseref:nc+ibaseref+ns-1,halo))
        end do
     end do
     fotherpanel(nc+1,1,1) = 0.25D0*(fotherpanel(nc+2,1,1)+fotherpanel(nc,1,1)&
          +fotherpanel(nc+1,2,1)+fotherpanel(nc+1,0,1))

     !
     ! ****************************************************************
     !
     ! fill halo for reconstruction on east neighbor panel projection
     !
     ! ****************************************************************
     !
     ! On the south panel projection the neighbors are arragened as follow (neast case):
     !
     !
     !             |     |  
     !         P   |  E  | 
     !             |-----|
     !            /
     !           /    S
     !          /
     ! 
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   |   | i |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   | w | i | i |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   | w | w | r | i | i |   |
     ! x---x---x---x---00000000000000000---x---x---x---x
     ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   |   |   0   |   | w | w 0 r | i | i |   |
     ! x---x---x---x---00000000000000000---x---x---x---x
     ! |   |   |   |   |   |   | w | ws| s | s |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   | s | s |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     !
     !
     !
     ! fill values on same panel projection ("r" and "i" on Figure above)
     !
     fotherpanel(nc+1:nc+nht,1:nc+nht,2)  = fcube(nc+1:nc+nht,1:nc+nht)
     !
     !
     ! fill in "w" on Figure above
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=0,nc+nh-(halo-1)
           ibaseref = fvm%ibase(i,halo,1)
           fotherpanel(nc+1-halo,i,2) = matmul_w(w(i,halo,:),fcube(nc+1-halo,ibaseref:ibaseref+ns-1))
        end do
     end do
     !
     ! fill in "s" on Figure above
     !
     w = fvm%halo_interp_weight(:,:,:,2)
     do halo=1,nhr
        do i=nc+1-nht+halo,nc+1
           !
           !
           ! !  P           |  E
           ! !              |
           ! !              |
           ! ================
           ! |              |
           ! |  S      |    |  <----- y 
           ! |         |    |           ^
           ! |       x |    |           |
           ! |         v    |           |
           ! !              |           |
           ! !    ----->    |           x
           ! !      y       |
           ! !              |
           ! ================
           !
           !
           ! shift (since we are using south weights from main panel interpolation
           !
           ibaseref = fvm%ibase(i,halo,2)-nc 
           !
           ! fotherpanel index: reverse
           !
           ! fcube index: due to rotation (see Figure above)
           !
           fotherpanel(nc+(nc+1-i),1-halo,2) = matmul_w(w(i,halo,:),fcube(nc+1-halo,ibaseref:ibaseref+ns-1))
        end do
     end do
     fotherpanel(nc,0,2) = 0.25D0*(fotherpanel(nc+1,0,2)+fotherpanel(nc-1,0,2)&
          +fotherpanel(nc,1,2)+fotherpanel(nc,-1,2))
  else if (fvm%cubeboundary==nwest) then
     !
     !
     ! 0-------\---\---\---0---\---\---\---0       !   --------x---------------x---------------x   
     ! | 0      \   \   \   0   \   \   \   0      !   |   | n | n | n | n | n | n |   |   |   |
     ! |   0-----\---\---\---0---\---\---\---0     !   --------x---------------x---------------x   
     ! |   | 0    \   \   \   0   \   \   \   0    !   | w | a | n | n | n | n | n | n |   |   |
     ! |\  |   000000000000000000000000000000000   !   --------00000000000000000---------------x   
     ! | -\|   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
     ! |   |\--0---------------0---------------0   !   --------0---------------0---------------x   
     ! |\--|   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
     ! |   |\--0---------------0---------------0   !   --------0---------------0---------------x   
     ! |\--|   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
     ! 0   |\--0---------------0---------------0   !   --------0---------------0---------------x   
     ! |0000   0   |   |   |   0   |   |   |   0   !   | w | w 0 r | r | r | r 0 r | H | H |   |
     ! |   |000000000000000000000000000000000000   !   --------00000000000000000---------------x   
     ! |\--|   0   |   |   |   0   |   |   |   0   !   | w | w | r | r | r | r | r | H | H |   |
     ! |   |\--0---------------0---------------0   !   --------x---------------x---------------x   
     ! |\--|   0   |   |   |   0   |   |   |   0   !   |   | w | H | H | H | H | H | H |   |   |
     ! |   |\--0---------------0---------------0   !   --------x---------------x---------------x   
     ! |\--|   0   |   |   |   0   |   |   |   0   !   |   |   | H | H | H | H | H |   |   |   |
     ! 0   |\--0---------------0---------------0   !   --------x---------------x---------------x   
     !  0000   0   |   |   |   0   |   |   |   0   !   |   |   |   |   |   |   |   |   |   |   |
     !      000000000000000000000000000000000000   !   --------x---------------x---------------x   
     !
     !
     !
     fpanel(1:nc+nht,1-nht:nc)=fcube(1:nc+nht,1-nht:nc)
     !
     ! west
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref=fvm%ibase(i,halo,1)                                              
           fpanel(1-halo ,i) = matmul_w(w(i,halo,:),fcube(1-halo ,ibaseref:ibaseref+ns-1))
        end do
     end do
     ! 
     ! north
     !
     w = fvm%halo_interp_weight(:,:,:,2)
     do halo=1,nhr
        do i=max(halo-nh,0),nc+nh-(halo-1)
           ibaseref = fvm%ibase(i,halo,2)
           fpanel(i,nc+halo) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,nc+halo  )) !north
        end do
     end do
     fpanel(0   ,nc+1)=0.25D0*(&
          fpanel(0,nc)+fpanel(1,nc+1)+fpanel(-1,nc+1)+fpanel(0,nc+2))

     !
     ! ****************************************************************
     !
     ! fill halo for reconstruction on north neighbor panel projection
     !
     ! ****************************************************************
     !
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   | i | i | i | i | i |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   | w | i | i | i | i | i | i |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   | w | w | r | r | r | r | r | i | i |   |
     !x---x---x---x---00000000000000000---x---x---x---x
     !|   |   | w | ws0 s | s | s | s 0 s | s |   |   |
     !x---x---x---x---0---x---x---x---0---x---x---x---x
     !|   |   |   | s 0 s | s | s | s 0 s |   |   |   |
     !x---x---x---x---0---x---x---x---0---x---x---x---x
     !|   |   |   |   0   |   |   |   0   |   |   |   |
     !x---x---x---x---0---x---x---x---0---x---x---x---x
     !|   |   |   |   0   |   |   |   0   |   |   |   |
     !x---x---x---x---00000000000000000---x---x---x---x
     !
     !
     ! fill values on same panel projection ("r" and "i" on Figure above)
     !
     fotherpanel(1:nc+nht,nc+1:nc+nht,1)  = fcube(1:nc+nht,nc+1:nc+nht)
     !
     !
     ! fill in "s" on Figure above
     !
     ! (use code from north above)
     !
     w = fvm%halo_interp_weight(:,:,:,2)
     do halo=1,nhr
        do i=max(halo-nh,0),nc+nh-(halo-1)
           ibaseref = fvm%ibase(i,halo,2)
           fotherpanel(i,nc+1-halo,1) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,nc+1-halo  ))
        end do
     end do
     !
     ! fill in "w" on Figure above
     !
     ! (use code from west above)
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=nc+1-nht+halo,nc+1
           ibaseref=fvm%ibase(i,halo,1)-nc
           fotherpanel(1-halo,nc-(i-(nc+1)),1) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,nc+1-halo))
        end do
     end do
     fotherpanel(0,nc,1)=0.25D0*(&
          fotherpanel(1,nc,1)+fotherpanel(-1,nc,1)+fotherpanel(0,nc+1,1)+fotherpanel(0,nc-1,1))

     !
     ! ****************************************************************
     !
     ! fill halo for reconstruction on west neighbor panel projection
     !
     ! ****************************************************************
     !
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   | n | n |   |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   | n | n | ne| e |   |   |   |   |   |   |
     !x---x---x---x---00000000000000000---x---x---x---x
     !|   | i | i | r 0 e | e |   |   0   |   |   |   |
     !x---x---x---x---0---x---x---x---0---x---x---x---x
     !|   | i | i | r 0 e | e |   |   0   |   |   |   |
     !x---x---x---x---0---x---x---x---0---x---x---x---x
     !|   | i | i | r 0 e | e |   |   0   |   |   |   |
     !x---x---x---x---0---x---x---x---0---x---x---x---x
     !|   | i | i | r 0 e | e |   |   0   |   |   |   |
     !x---x---x---x---00000000000000000---x---x---x---x
     !|   | i | i | r | e | e |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   | i | i | e |   |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   | i |   |   |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---
     !
     !
     ! fill values on same panel projection ("r" and "i" on Figure above)
     !
     fotherpanel(1-nht:nc,1-nht:nc,2)  = fcube(1-nht:nc,1-nht:nc)
     !
     !
     ! fill in "e" on Figure above
     !
     ! (use code from west above)
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref=fvm%ibase(i,halo,1)                                              
           fotherpanel(halo ,i,2) = matmul_w(w(i,halo,:),fcube(halo ,ibaseref:ibaseref+ns-1))
        end do
     end do
     !
     !
     ! fill in "n" on Figure above
     !
     ! (use code from north above)
     !
     w = fvm%halo_interp_weight(:,:,:,2)
     do halo=1,nhr
        do i=0,nht-halo
           ibaseref = fvm%ibase(i,halo,2)+nc
           fotherpanel(1-i,nc+halo,2) = matmul_w(w(i,halo,:),fcube(halo,ibaseref:ibaseref+ns-1)) !north
        end do
     end do
     fotherpanel(1,nc+1,2)=0.25D0*(&
          fotherpanel(2,nc+1,2)+fotherpanel(0,nc+1,2)+fotherpanel(1,nc+2,2)+fotherpanel(1,nc,2))

  else if (fvm%cubeboundary==neast) then
     !
     !
     !     0---/---/---/---0---/---/---/-------0     !   x---------------x---------------x--------
     !    0   /   /   /   0   /   /   /      0 |     !   |   |   |   |   | n | n | n | n | n |   |
     !   0---/---/---/---0---/---/---/-----0   |     !   x---------------x---------------x--------
     !  0   /   /   /   0   /   /   /    0 |   |     !   |   |   |   | n | n | n | n | n | a | e |
     ! 000000000000000000000000000000000   |   |     !   x---------------00000000000000000--------
     ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
     ! 0---------------0---------------0--/|   |     !   x---------------0---------------0--------
     ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
     ! 0---------------0---------------0--/|   |     !   x---------------0---------------0--------
     ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
     ! 0---------------0---------------0--/|   0     !   x---------------0---------------0--------
     ! 0   |   |   |   0   |   |   |   0   0000|     !   |   |   | H | H 0 r | r | r | r 0 e | e |
     ! 000000000000000000000000000000000000|   |     !   x---------------00000000000000000--------
     ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   | H | H | r | r | r | r | e | e |
     ! 0---------------0---------------0--/|   |     !   x---------------x---------------x--------
     ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   |   | H | H | H | H | H | e |   |
     ! 0---------------0---------------0--/|   |     !   x---------------x---------------x--------
     ! 0   |   |   |   0   |   |   |   0   |--/|     !   |   |   |   |   | H | H | H | H |   |   |
     ! 0---------------0---------------0--/|   0     !   x---------------x---------------x--------
     ! 0   |   |   |   0   |   |   |   0   0000      !   |   |   |   |   |   |   |   |   |   |   |
     ! 000000000000000000000000000000000000          !   x---------------x---------------x--------
     !
     !
     !
     fpanel(1-nht:nc,1-nht:nc)=fcube(1-nht:nc,1-nht:nc)
!     fotherpanel (nc+1 :nc+nht ,1-nht:nc+nht)=fcube(nc+1 :nc+nht ,1-nht:nc+nht)
     !
     ! east
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref=fvm%ibase(i,halo,1 )                                              
           fpanel(nc+halo,i) = matmul_w(w(i,halo,:),fcube(nc  +halo,ibaseref:ibaseref+ns-1))
        end do
     end do
     !
     ! north
     !
!     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref=fvm%ibase(i,halo,1) 
           fpanel(i,nc+halo) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,nc+halo  )) !north
        end do
     end do
     fpanel(nc+1,nc+1)=0.25D0*(&
          fpanel(nc,nc+1)+fpanel(nc+1,nc)+fpanel(nc+1,nc+2)+fpanel(nc+2,nc+1))
     !
     ! ****************************************************************
     !
     ! fill halo for reconstruction on north neighbor panel projection
     !
     ! ****************************************************************
     !
     ! On the north panel projection the neighbors are arragened as follow (seast case):
     !
     !
     !             |     |  
     !             |  N  | E
     !             |-----|
     !                   \
     !                S   \
     !                     \
     ! 
     ! 
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   |   |   |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   |   | i | i | i | i | i |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   |   | i | i | i | i | i | i | e |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     ! |   | i | i | r | r | r | r | r | e | e |   |   |
     ! x---x---x---x---00000000000000000---x---x---x---x
     ! |   |   | s | s 0 s | s | s | s 0 se| e |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   |   | s 0 s | s | s | s 0 s |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   |   |   0   |   |   |   0   |   |   |   |
     ! x---x---x---x---0---x---x---x---0---x---x---x---x
     ! |   |   |   |   0   |   |   |   0   |   |   |   |
     ! x---x---x---x---00000000000000000---x---x---x---x
     ! |   |   |   |   |   |   |   |   |   |   |   |   |
     ! x---x---x---x---x---x---x---x---x---x---x---x---x
     !
     !
     ! fill values on same panel projection ("r" and "i" on Figure above)
     !
     fotherpanel(1-nht:nc,nc+1:nc+nht,1)  = fcube(1-nht:nc,nc+1:nc+nht)
     !
     ! fill in "s" on Figure above
     !
     ! (use north case from above and shift/reverse j-index
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref=fvm%ibase(i,halo,1) 
           fotherpanel (i,nc+1-halo,1) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,nc+1-halo))
        end do
     end do
     !
     ! fill in "e" on Figure above
     !
     w = fvm%halo_interp_weight(:,:,:,2)
     do halo=1,nhr
        do i=max(halo-nh,0),nht-halo
           ibaseref=fvm%ibase(i,halo,2) +nc
           !
           ! fotherpanel uses indexing of main panel's projection
           ! fcube: rotated indexing
           !
           fotherpanel (nc+halo,nc+i,1) = matmul_w(w(i,halo,:),fcube(ibaseref:ibaseref+ns-1,nc+1-halo))
        end do
     end do
     fotherpanel(nc+1,nc,1)=0.25D0*(&
          fotherpanel(nc+2,nc,1)+fotherpanel(nc,nc,1)+fotherpanel(nc+1,nc+1,1)+fotherpanel(nc+1,nc-1,1))
     !
     ! ****************************************************************
     !
     ! fill halo for reconstruction on east neighbor panel projection
     !
     ! ****************************************************************
     !
     ! On the north panel projection the neighbors are arragened as follow (seast case):
     !
     !
     !           \    N
     !            \  
     !             \------
     !             |     |
     !         P   |  E  |
     !             |     |
     !             -------
     !
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   |   | n | n |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   | w | wn| n | n |   |   |
     !x---x---x---x---00000000000000000---x---x---x---x
     !|   |   |   |   0   |   | w | w 0 r | i | i |   |
     !x---x---x---x---0---x---x---x---0---x---x---x---x
     !|   |   |   |   0   |   | w | w 0 r | i | i |   |
     !x---x---x---x---0---x---x---x---0---x---x---x---x
     !|   |   |   |   0   |   | w | w 0 r | i | i |   |
     !x---x---x---x---0---x---x---x---0---x---x---x---x
     !|   |   |   |   0   |   | w | w 0 r | i | i |   |
     !x---x---x---x---00000000000000000---x---x---x---x
     !|   |   |   |   |   |   | w | w | r | i | i |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   |   | w | i | i |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   |   |   | i |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---x
     !|   |   |   |   |   |   |   |   |   |   |   |   |
     !x---x---x---x---x---x---x---x---x---x---x---x---
     !
     !
     !
     ! fill values on same panel projection ("r" and "i" on Figure above)
     !
     fotherpanel(nc+1:nc+nht,1-nht:nc,2)  = fcube(nc+1:nc+nht,1-nht:nc)
     !
     ! fill in "w" on Figure above
     !
     ! (use east case from above and shift/reverse j-index
     !
     w = fvm%halo_interp_weight(:,:,:,1)
     do halo=1,nhr
        do i=halo-nh,min(nc+nh-(halo-1),nc+1)
           ibaseref=fvm%ibase(i,halo,1 )                                              
           fotherpanel(nc+1-halo,i,2) = matmul_w(w(i,halo,:),fcube(nc+1-halo,ibaseref:ibaseref+ns-1))
        end do
     end do
     !
     ! fill in "n" on Figure above
     !
     w = fvm%halo_interp_weight(:,:,:,2)
     do halo=1,nhr
        do i=max(halo-nh,0),nht-halo
           ibaseref=fvm%ibase(i,halo,2) +nc
           !
           ! fotherpanel uses indexing of main panel's projection
           ! fcube: rotated indexing
           !
           fotherpanel (nc+i,nc+halo,2) = matmul_w(w(i,halo,:),fcube(nc+1-halo,ibaseref:ibaseref+ns-1))
        end do
     end do
     fotherpanel(nc,nc+1,2)=0.25D0*(&
          fotherpanel(nc+1,nc+1,2)+fotherpanel(nc-1,nc+1,2)+fotherpanel(nc,nc+2,2)+fotherpanel(nc,nc,2))

  end if
end subroutine fill_halo


subroutine get_reconstruction(fvm, fpanel, fotherpanel, recons)
  implicit none

  type (fvm_struct), intent(in)                                     :: fvm
!  real (kind=real_kind), dimension(    1-nht:nc+nht,1-nht:nc+nht  ), intent(in)   :: fpanel
!  real (kind=real_kind), dimension(    1-nht:nc+nht,1-nht:nc+nht,2), intent(in)   :: fotherpanel
!  real (kind=real_kind), dimension(1:5,1-nhe:nc+nhe,1-nhe:nc+nhe  ), intent(inout):: recons
!  real (kind=real_kind), dimension(1-nht:nc+nht,1-nht:nc+nht,1:5,1:5)   :: fcube_matrix

  real (kind=real_kind), dimension(  1-nht:,1-nht:  ), intent(in)   :: fpanel
  real (kind=real_kind), dimension(  1-nht:,1-nht:,:), intent(in)   :: fotherpanel
  real (kind=real_kind), dimension(:,1-nhe:,1-nhe:  ), intent(inout):: recons

  real (kind=real_kind), dimension(1-nht:nc+nht,1-nht:nc+nht,1:5,1:5)   :: fcube_matrix
  integer  :: i, j, k, h
  !
  ! on panel reconstruction
  !
  recons       = 0.0D0
  fcube_matrix = 0.0D0 
  !
  call fill_fcube_matrix(fvm,fpanel,fotherpanel,fcube_matrix) 
  do j = 1-nhe,nc+nhe
     do i = 1-nhe,nc+nhe
        !
        ! loop over reconstruction coefficients
        !
        do k=1,5
           do h=1,5
              !
              ! loop over finite-difference (and metric terms) matrix
              !
              recons(k,i,j) = recons(k,i,j) + fvm%recons_matrix(h,k,i,j)*fcube_matrix(i,j,k,h)
           end do
        end do
        !
        ! stretching of d^2f/dx^2 and d^2f/dy^2 (can not be cassed in matrix form)
        !
        recons(3,i,j) = (recons(3,i,j) + fvm%recons_matrix(1,6,i,j)*recons(1,i,j))*fvm%recons_matrix(2,6,i,j)
        recons(4,i,j) = (recons(4,i,j) + fvm%recons_matrix(3,6,i,j)*recons(2,i,j))*fvm%recons_matrix(4,6,i,j)
     end do
  end do
end subroutine get_reconstruction


subroutine fill_fcube_matrix(fvm,fpanel,fotherpanel,fcube_matrix)
  implicit none
  type (fvm_struct), intent(in) :: fvm
!  real (kind=real_kind), dimension(1-nht:nc+nht,1-nht:nc+nht        ), intent(in)  :: fpanel
!  real (kind=real_kind), dimension(1-nht:nc+nht,1-nht:nc+nht,2      ), intent(in)  :: fotherpanel
!  real (kind=real_kind), dimension(1-nht:nc+nht,1-nht:nc+nht,1:5,1:5), intent(out) :: fcube_matrix

  real (kind=real_kind), dimension(1-nht:,1-nht:    ), intent(in)  :: fpanel
  real (kind=real_kind), dimension(1-nht:,1-nht:,:  ), intent(in)  :: fotherpanel
  real (kind=real_kind), dimension(1-nht:,1-nht:,:,:), intent(out) :: fcube_matrix


  integer  :: i,j, k, kmax

  integer, dimension(2)  :: imin,imax,jmin,jmax,invx,invy
  logical, dimension(2)  :: swap

  do j = fvm%jy_min, fvm%jy_max-1
     do i = fvm%jx_min, fvm%jx_max-1
        call get_fcube_matrix(fcube_matrix(i,j,:,:),fpanel(i-2:i+2,j-2:j+2))
     end do
  end do

  if (fvm%cubeboundary>0) then
     imin(1) = fvm%jx_min1; imax(1) = fvm%jx_max1-1; jmin(1) = fvm%jy_min1; jmax(1) = fvm%jy_max1-1
     swap(1) = fvm%swap1  ; invx(1) = fvm%invx1    ; invy(1) = fvm%invy1
     if (fvm%cubeboundary>4) then
        kmax = 2
        imin(2) = fvm%jx_min2; imax(2) = fvm%jx_max2-1
        jmin(2) = fvm%jy_min2; jmax(2) = fvm%jy_max2-1
        swap(2) = fvm%swap2
        invx(2) = fvm%invx2  ; invy(2) = fvm%invy2
     else
        kmax = 1
     end if

     do k=1,kmax
        if (swap(k)) then
           do j=jmin(k),jmax(k)
              do i=imin(k),imax(k)
                 call get_fcube_matrix_otherpanel_swap(&
                      fcube_matrix(i,j,:,:),fotherpanel(i-2:i+2,j-2:j+2,k),invx(k),invy(k))
              end do
           end do
        else
           do j=jmin(k),jmax(k)
              do i=imin(k),imax(k)
                 call get_fcube_matrix_otherpanel(&
                      fcube_matrix(i,j,:,:),fotherpanel(i-2:i+2,j-2:j+2,k),invx(k),invy(k))
              end do
           end do
        end if
     end do
  end if
end subroutine fill_fcube_matrix

subroutine get_fcube_matrix(fcube_matrix,f)
  implicit none
!  real (kind=real_kind), dimension(1:5,1:5),intent(out) :: fcube_matrix
  real (kind=real_kind), dimension(:,:),intent(out) :: fcube_matrix
  real (kind=real_kind), dimension(-2:,-2:), intent(in)   :: f
  ! for x-derivative
  fcube_matrix(1,:) = f(:,0); fcube_matrix(3,:) = fcube_matrix(1,:)
  ! for y-derivative
  fcube_matrix(2,:) = f(0,:); fcube_matrix(4,:) = fcube_matrix(2,:)
  ! cross terms
  fcube_matrix(5,1) = f(-1,-1); fcube_matrix(5,2) = f(-1, 1)
  fcube_matrix(5,3) = f( 0, 0); fcube_matrix(5,4) = f( 1,-1)
  fcube_matrix(5,5) = f( 1, 1)
end subroutine get_fcube_matrix

subroutine get_fcube_matrix_otherpanel(fcube_matrix,f,invx,invy)
  implicit none
!  real (kind=real_kind), dimension(1:5,1:5)  , intent(out) :: fcube_matrix
  real (kind=real_kind), dimension(:,:)  , intent(out) :: fcube_matrix
  real (kind=real_kind), dimension(-2:,-2:), intent(in)  :: f
  integer                                    , intent(in)  :: invx, invy

  integer                                                  :: inv

  inv = invx*invy

  ! for x-derivative
  fcube_matrix(1,:) = invx*f(:,0); fcube_matrix(3,:) = f(:,0)
  ! for y-derivative
  fcube_matrix(2,:) = invy*f(0,:); fcube_matrix(4,:) = f(0,:)
  ! cross terms
  fcube_matrix(5,1) = inv*f(-1,-1); fcube_matrix(5,2) = inv*f(-1, 1)
  fcube_matrix(5,3) = inv*f( 0, 0); fcube_matrix(5,4) = inv*f( 1,-1)
  fcube_matrix(5,5) = inv*f( 1, 1)
end subroutine get_fcube_matrix_otherpanel

subroutine get_fcube_matrix_otherpanel_swap(fcube_matrix,f,invx,invy)
  implicit none
!  real (kind=real_kind), dimension(1:5,1:5)  , intent(out) :: fcube_matrix
  real (kind=real_kind), dimension(:,:)  , intent(out) :: fcube_matrix
  real (kind=real_kind), dimension(-2:,-2:), intent(in)  :: f
  integer                                    , intent(in)  :: invx, invy

  integer                                                  :: inv

  inv = invx*invy

  ! for x-derivative
  fcube_matrix(2,:) = invx*f(:,0); fcube_matrix(4,:) = f(:,0)
  ! for y-derivative
  fcube_matrix(1,:) = invy*f(0,:); fcube_matrix(3,:) = f(0,:)
  ! cross terms
  fcube_matrix(5,1) = inv*f(-1,-1); fcube_matrix(5,2) = inv*f(-1, 1)
  fcube_matrix(5,3) = inv*f( 0, 0); fcube_matrix(5,4) = inv*f( 1,-1)
  fcube_matrix(5,5) = inv*f( 1, 1)
end subroutine get_fcube_matrix_otherpanel_swap
 


!
! the subroutines below are just for debugging
!
!
subroutine debug_halo(fvm,fcubenew,fpanel)!
  implicit none
  type (fvm_struct), intent(in)                                 :: fvm
        
  real (kind=real_kind),   &
        dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc), intent(in)            :: fcubenew !dbg
  real (kind=real_kind), dimension(1-nht:nc+nht, 1-nht:nc+nht) :: fpanel

  integer (kind=int_kind)                                         :: i, j

  integer (kind=int_kind)                                         :: i1,i2, j1,j2,istart,iend, count

  logical, dimension(1-nhc:nc+nhc,1-nhc:nc+nhc) :: lhalo

  call indicator_fct_recons(fvm,lhalo,.false.)
  
  do j=1-nhc,nc+nhc
     do i=1-nhc,nc+nhc
        if (lhalo(i,j)) then
           if (ABS(fcubenew(i,j)-fpanel(i,j))>1.0E-12 &
                   .or.fpanel(i,j).ne.fpanel(i,j).or.&
                   fcubenew(i,j).ne.fcubenew(i,j)) then
              write(*,*) "difference",i,j,fcubenew(i,j),fpanel(i,j),fcubenew(i,j)-fpanel(i,j)
              stop
           else
!              write(*,*) "pass test",i,j,fcubenew(i,j),fpanel(i,j),fcubenew(i,j)-fpanel(i,j)
           end if
        end if
     end do
  end do
end subroutine debug_halo


subroutine debug_halo_recons(fvm,recons,recons_trunk)!
  implicit none
  type (fvm_struct), intent(in)                                 :: fvm

  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)   :: recons, recons_trunk

  integer (kind=int_kind)                                         :: i, j, k

  integer (kind=int_kind)                                         :: i1,i2, j1,j2,istart,iend, count

  logical, dimension(1-nhc:nc+nhc,1-nhc:nc+nhc) :: lhalo

  call indicator_fct_recons(fvm,lhalo,.true.)
  
  call print_which_case(fvm)

  do j=1-nhe,nc+nhe
     do i=1-nhe,nc+nhe
        if (lhalo(i,j)) then
           do k=1,5
              if (ABS(recons(k,i,j)-recons_trunk(k,i,j))>1.0E-12&
                   .or.recons(k,i,j).ne.recons(k,i,j).or.&
                       recons(k,i,j).ne.recons(k,i,j)) then
              write(*,*) "recons difference",k,i,j,recons(k,i,j)-recons_trunk(k,i,j),recons(k,i,j),recons_trunk(k,i,j)
              stop
           else
!              write(*,*) "recons pass test",k,i,j,recons(k,i,j)-recons_trunk(k,i,j),recons(k,i,j),recons_trunk(k,i,j)
!              write(*,*) "recons pass test",i,j,k,fcubenew(i,j),fpanel(i,j),fcubenew(i,j)-fpanel(i,j)
           end if
        end do
        end if
     end do
  end do
write(*,*) "recons pass test for "
call print_which_case(fvm)
end subroutine debug_halo_recons


subroutine indicator_fct_recons(fvm,lhalo,lrecons)!
  implicit none
  type (fvm_struct), intent(in)                                 :: fvm
        
  integer (kind=int_kind)                                         :: i, j

  integer (kind=int_kind)                                         :: i1,i2, j1,j2,istart,iend, count

  integer (kind=int_kind)                                         :: imin,imax,jmin,jmax

  logical, dimension(1-nhc:nc+nhc,1-nhc:nc+nhc), intent(out) :: lhalo
  CHARACTER(len=2), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)     :: chalo
  logical, intent(in) :: lrecons
  !
  lhalo = .false.
  chalo = "  "

  imin = 1-nhe; imax = nc+nhe; jmin = 1-nhe; jmax = nc+nhe

  select case (fvm%cubeboundary)
  case(west) 
     imin = 1
  case(east) 
     imax = nc
  case(south) 
     jmin = 1
  case(north) 
     jmax = nc
  case(nwest) 
     jmax = nc; imin = 1
  case(neast) 
     jmax = nc; imax = nc
  case(seast) 
     jmin = 1; imax = nc
  case(swest)
     jmin = 1; imin = 1
  case default
     imin = 1-nhe; imax = nc+nhe; jmin = 1-nhe; jmax = nc+nhe
  end select

  do j=max(jmin,1-nhe),min(jmax,nc+nhe)
     do i=max(imin,1-nhe),min(imax,nc+nhe)
        call get_stencil (lhalo,chalo,i,j,fvm%cubeboundary,lrecons)!        
     end do
  end do
  
  if (lplot) call plot_stencil(lhalo,chalo)
end subroutine indicator_fct_recons


subroutine debug_halo_neighbor(fvm,fcubenew,fotherpanel)!
  implicit none
  type (fvm_struct), intent(in)                                 :: fvm
        
!  real (kind=real_kind),   &
!phl        dimension(1-nhc:nc+nhc, 1-nhc:nc+nhc)                     :: fcubenew
 
  real (kind=real_kind), dimension(1-nht:nc+nht, 1-nht:nc+nht,2) :: fotherpanel,fcubenew

  integer (kind=int_kind)                                         :: i, j

  integer (kind=int_kind)                                         :: i1,i2, j1,j2,istart,iend, count,k

  logical, dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,2) :: lhalo
  CHARACTER(len=2), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,2)   :: chalo

  call indicator_fct_recons_neighbor(fvm,lhalo,chalo,.false.)
  
  !
  ! plot halo to screen
  !
  call print_which_case(fvm)
  if (fvm%cubeboundary>0) then
     do k=1,2
        write(*,*) "k=",k
        do j=1-nhc,nc+nhc
           do i=1-nhc,nc+nhc
              if (lhalo(i,j,k)) then
                 if (ABS(fcubenew(i,j,k)-fotherpanel(i,j,k))>1.0E-12&
                      .or.fotherpanel(i,j,k).ne.fotherpanel(i,j,k).or.&
                      fcubenew(i,j,k).ne.fcubenew(i,j,k)) then
                    write(*,*) "difference",i,j,fcubenew(i,j,k),fotherpanel(i,j,k),fcubenew(i,j,k)-fotherpanel(i,j,k)
                    stop
                 else
!                    write(*,*) "pass test",i,j,fcubenew(i,j,k),fotherpanel(i,j,k),fcubenew(i,j,k)-fotherpanel(i,j,k)
                 end if
              end if
           end do
        end do
     end do
  end if
end subroutine debug_halo_neighbor


subroutine debug_halo_neighbor_recons(fvm,recons,recons_trunk)
  implicit none
  type (fvm_struct), intent(in)                                 :: fvm
        
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)   :: recons, recons_trunk
  integer (kind=int_kind)                                         :: i, j, h, k


  logical, dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,2) :: lhalo
  CHARACTER(len=2), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,2)   :: chalo

  call indicator_fct_recons_neighbor(fvm,lhalo,chalo,.true.)
  
  !
  ! plot halo to screen
  !

  do j=1-nhe,nc+nhe
     do i=1-nhe,nc+nhe
        do h=1,2
           if (lhalo(i,j,h)) then
              do k=1,5
                 if (ABS(recons(k,i,j)-recons_trunk(k,i,j))>1.0E-12&
                      .or.recons(k,i,j).ne.recons(k,i,j).or.&
                      recons(k,i,j).ne.recons(k,i,j)) then
                    write(*,*) "recons difference neighbor h,k,i,j",h,k,i,j,recons(k,i,j)-recons_trunk(k,i,j),recons(k,i,j),recons_trunk(k,i,j)
                    stop
                 else
!                    write(*,*) "recons pass test neighbor h,k,i,j",h,k,i,j,recons(k,i,j)-recons_trunk(k,i,j),recons(k,i,j),recons_trunk(k,i,j)
                    !              write(*,*) "recons pass test",i,j,k,fcubenew(i,j),fpanel(i,j),fcubenew(i,j)-fpanel(i,j)
                 end if
              end do
           end if
        end do
     end do
  end do
write(*,*) "recons neighbor pass test for "
call print_which_case(fvm)


end subroutine debug_halo_neighbor_recons

subroutine indicator_fct_recons_neighbor(fvm,lhalo,chalo,lrecons)!
  implicit none
  type (fvm_struct), intent(in)                                 :: fvm
        
  integer (kind=int_kind)                                         :: i, j

  integer (kind=int_kind)                                         :: i1,i2, j1,j2,istart,iend, count, halo

  logical, dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,2), intent(out) :: lhalo
  CHARACTER(len=2), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc,2)  , intent(out) :: chalo

  logical, dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)          :: lh
  CHARACTER(len=2), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc) :: ch

  logical, intent(in) :: lrecons

  !
  lhalo = .false.
  chalo = "  "

  if (fvm%cubeboundary==west) then
     lh    = .false.
     ch    = "  "
     !
     ! index here follows "as if we were on the south panel
     !
     !
     ! should be recoded with get_stencil
     !
     do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
           lhalo(halo,i,1) = .true.
        end do
     end do
     do halo=1,nht
        do i=halo-nht,nc+nht-(halo-1)
           lhalo(1-halo,i,1) = .true.
        end do
     end do
     
  else  if (fvm%cubeboundary==south) then
     do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
           lhalo(i,halo,1) = .true.
        end do
     end do
     do halo=1,nht
        do i=halo-nht,nc+nht-(halo-1)
           lhalo(i,1-halo,1) = .true.  
        end do
     end do
  else if (fvm%cubeboundary==east) then
     do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
           lhalo (nc+1-halo ,i,1) = .true.
        end do
     end do
     do halo=1,nht
        do i=halo-nht,nc+nht-(halo-1)
           lhalo (nc+halo ,i,1) = .true.
        end do
     end do
  else if (fvm%cubeboundary==north) then
     do halo=1,nhr
        do i=halo-nh,nc+nh-(halo-1)
           lhalo(i,nc+1-halo,1) = .true.
        end do
     end do
     do halo=1,nht
        do i=halo-nht,nc+nht-(halo-1)
           lhalo(i,nc+halo,1) = .true.
        end do
     end do
  else if (fvm%cubeboundary==swest) then
     !
     !*******************************************************************
     !
     !  panel reconstruction (which is nwest for south panel projection)
     !
     !*******************************************************************
     !
     write(*,*) "reconstruction is nwest for south panel projection"
     lh    = .false.
     ch    = "  "
     !
     ! index here follows "as if we were on the south panel
     !
     do j=nc+1-nhe,nc
        do i=1,nc+nhe
           call get_stencil (lh,ch,i,j,nwest,lrecons)!
        end do
     end do
     !
     ! shift back to main panel indexing
     !
     lhalo(:,(nc-nh)-nc:(nc+nhr)-nc,1) = lh(:,(nc-nh):(nc+nhr))
     chalo(:,(nc-nh)-nc:(nc+nhr)-nc,1) = ch(:,(nc-nh):(nc+nhr))
     !
     if (lplot) call plot_stencil(lhalo(:,:,1),chalo(:,:,1))
     !
     !*********************************************************************
     !
     ! West panel reconstruction (which is seast for west panel projection)
     !
     !*********************************************************************
     !
     lh    = .false.
     ch    = "  "
     write(*,*) "reconstruction is seast for west panel projection"
     do j=1,nc+nhe
        do i=nc+1-nhe,nc
           call get_stencil (lh(:,:),ch(:,:),i,j,seast,lrecons)!
        end do
     end do
     !
     ! shift back to main panel indexing
     !
     lhalo(nc+1-nht-nc:nc+nhr-nc,:,2) = lh(nc+1-nht:nc+nhr,:)
     chalo(nc+1-nht-nc:nc+nhr-nc,:,2) = ch(nc+1-nht:nc+nhr,:)
     if (lplot) call plot_stencil(lhalo(:,:,2),chalo(:,:,2))
  else if (fvm%cubeboundary==seast) then
     !
     !*******************************************************************
     !
     !  panel reconstruction (which is neast for south panel projection)
     !
     !*******************************************************************
     !
     write(*,*) "case seast"
     write(*,*) "reconstruction is neast for south panel projection"
     lh    = .false.
     ch    = "  "
     !
     ! index here follows "as if we were on the south panel"
     !
     do j=nc+1-nhe,nc
        do i=1-nhe,nc
           call get_stencil (lh,ch,i,j,neast,lrecons)!
        end do
     end do
     !
     ! shift back to main panel indexing
     !
     lhalo(:,(nc+1-nht)-nc:(nc+nhr)-nc,1) = lh(:,(nc+1-nht):(nc+nhr))
     chalo(:,(nc+1-nht)-nc:(nc+nhr)-nc,1) = ch(:,(nc+1-nht):(nc+nhr))
     if (lplot) call plot_stencil(lhalo(:,:,1),chalo(:,:,1))
     !
     !*******************************************************************
     !
     !  panel reconstruction (which is swest for east panel projection)
     !
     !*******************************************************************
     !
     write(*,*) "case seast"
     write(*,*) "reconstruction is seast for east panel projection"
     lh    = .false.
     ch    = "  "
     !
     ! index here follows "as if we were on the east panel"
     !
     do j=1,nc+nhe
        do i=1,nhe
           call get_stencil (lh,ch,i,j,swest,lrecons)!
        end do
     end do
     !
     ! shift back to main panel indexing
     !
     lhalo(nc+1-nhr:nc+nht,:,2) = lh(1-nhr:nht,:)
     chalo(nc+1-nhr:nc+nht,:,2) = ch(1-nhr:nht,:)

     !
     if (lplot) call plot_stencil(lhalo(:,:,2),chalo(:,:,2))
  else if (fvm%cubeboundary==neast) then
     !
     !*******************************************************************
     !
     !  panel reconstruction (which is seast for north panel projection)
     !
     !*******************************************************************
     !
     write(*,*) "case neast"
     write(*,*) "reconstruction is seast for north panel projection"
     lh    = .false.
     ch    = "  "
     !
     ! index here follows "as if we were on the north panel"
     !
     do j=1,nhe
        do i=1-nhe,nc
           call get_stencil (lh,ch,i,j,seast,lrecons)!
        end do
     end do
     !
     ! shift back to main panel indexing
     !
     lhalo(:,1-nhr+nc:nht+nc,1) = lh(:,1-nhr:nht)
     chalo(:,1-nhr+nc:nht+nc,1) = ch(:,1-nhr:nht)
     if (lplot) call plot_stencil(lhalo(:,:,1),chalo(:,:,1))
     !
     !*******************************************************************
     !
     !  panel reconstruction (which is seast for north panel projection)
     !
     !*******************************************************************
     !
     write(*,*) "case neast"
     write(*,*) "reconstruction is nwest for east panel projection"
     lh    = .false.
     ch    = "  "
     !
     ! index here follows "as if we were on the east panel"
     !
     do j=1-nhe,nc
        do i=1,nhe
           call get_stencil (lh,ch,i,j,nwest,lrecons)!
        end do
     end do
     !
     ! shift back to main panel indexing
     !
     lhalo(1-nhr+nc:nht+nc,:,2) = lh(1-nhr:nht,:)
     chalo(1-nhr+nc:nht+nc,:,2) = ch(1-nhr:nht,:)
     if (lplot) call plot_stencil(lhalo(:,:,2),chalo(:,:,2))
  else if (fvm%cubeboundary==nwest) then
     !
     !*******************************************************************
     !
     !  panel reconstruction (which is swest for north panel projection)
     !
     !*******************************************************************
     !
     write(*,*) "case nwest"
     write(*,*) "reconstruction is swest for north panel projection"
     lh    = .false.
     ch    = "  "
     !
     ! index here follows "as if we were on the north panel"
     !
     do j=1,nhe
        do i=1,nc+nhe
           call get_stencil (lh,ch,i,j,swest,lrecons)!
        end do
     end do
     !
     ! shift back to main panel indexing
     !
     lhalo(:,1-nhr+nc:nht+nc,1) = lh(:,1-nhr:nht)
     chalo(:,1-nhr+nc:nht+nc,1) = ch(:,1-nhr:nht)
     if (lplot) call plot_stencil(lhalo(:,:,1),chalo(:,:,1))
     !
     !*******************************************************************
     !
     !  panel reconstruction (which is swest for north panel projection)
     !
     !*******************************************************************
     !
     write(*,*) "case nwest"
     write(*,*) "reconstruction is neast for west panel projection"
     lh    = .false.
     ch    = "  "
     !
     ! index here follows "as if we were on the west panel"
     !
     do j=1-nhe,nc
        do i=nc+1-nhe,nc
           call get_stencil (lh,ch,i,j,neast,lrecons)!
        end do
     end do
     !
     ! shift back to main panel indexing
     !
     lhalo(nc+1-nht-nc:nc+nhr-nc,:,2) = lh(nc+1-nht:nc+nhr,:)
     chalo(nc+1-nht-nc:nc+nhr-nc,:,2) = ch(nc+1-nht:nc+nhr,:)
     if (lplot) call plot_stencil(lhalo(:,:,2),chalo(:,:,2))
  else if (fvm%cubeboundary==0) then
     lhalo = .false.
  else
     lhalo=.true.
  end if
end subroutine indicator_fct_recons_neighbor


subroutine print_which_case(fvm)
  implicit none
  type (fvm_struct), intent(in) :: fvm

  select case (fvm%cubeboundary)
  case(west) 
     write(*,*) "case west"
  case(east) 
     write(*,*) "case east"
  case(north) 
     write(*,*) "case north"
  case(south) 
     write(*,*) "case south"
  case(nwest) 
     write(*,*) "case nwest"
  case(neast) 
     write(*,*) "case neast"
  case(seast) 
     write(*,*) "case seast"
  case(swest)
     write(*,*) "case swest" 
  case (0)
     write(*,*) "enterior element"
  case default
     write(*,*) "cubeboundary ill-defined in print_which_case in fvm_reconstruction_mod.F90"
  end select

end subroutine print_which_case


subroutine plot_stencil(lhalo,chalo)
  implicit none
  integer (kind=int_kind)                                         :: i, j

  logical, dimension(1-nhc:nc+nhc,1-nhc:nc+nhc), intent(in) :: lhalo
  CHARACTER(len=2), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc), intent(in)   :: chalo
  !
  ! plot stencil to screen
  !
  write(*,*) "   "
  write(*,*) "   "
  do j=nc+nhc,1-nhc,-1
     do i=1-nhc,nc+nhc
        if ((j==nc.or.j==0).and.i>0.and.i<nc+2) then
           if (i==nc+1) then
              write(*,'(A4)',advance='no') "0---"
           else
              write(*,'(A4)',advance='no') "0000"
           end if
        else                
           if (j>0.and.j<nc.and.(i==1.or.i==nc+1)) then
              write(*,'(A4)',advance='no') "0---"
           else
              write(*,'(A4)',advance='no') "x---"
           end if
        end if
     end do
     write(*,'(A1)') "x"
     do i=1-nhc,nc+nhc
        if ((i==1.or.i==nc+1).and.(j>0.and.j<nc+1)) then
           write(*,'(A2)',advance='no') "0 "
        else
           write(*,'(A2)',advance='no') "| "
        end if
        write(*,'(A2)',advance='no') chalo(i,j)
        !
     end do
     write(*,'(A1)') "|"
  end do
  do i=1-nhc,nc+nhc
     write(*,'(A4)',advance='no') "x---"
  end do
  write(*,*) "   "
  write(*,*) "   "
  write(*,*) "   "
end subroutine plot_stencil


subroutine get_stencil (lhalo,chalo,i,j,cubeboundary,lrecons)!
  implicit none
  integer (kind=int_kind)                               , intent(in)    :: i, j
  integer (kind=int_kind)                               , intent(in)    :: cubeboundary
  logical         , dimension(1-nhc:nc+nhc,1-nhc:nc+nhc), intent(inout) :: lhalo
  CHARACTER(len=2), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc), intent(inout) :: chalo
  logical, intent(in) :: lrecons
  !
  integer (kind=int_kind)   :: ii, jj
  !
  if (lrecons) then
     lhalo(i  ,j  ) = .true.;                             chalo(i  ,j  ) = "r "
  else
     lhalo(i  ,j  ) = .true.;                             chalo(i  ,j  ) = "r "
     lhalo(i+1,j  ) = .true.; if (chalo(i+1,j  ).ne."r ") chalo(i+1,j  ) = "i "
     lhalo(i-1,j  ) = .true.; if (chalo(i-1,j  ).ne."r ") chalo(i-1,j  ) = "i "
     lhalo(i  ,j+1) = .true.; if (chalo(i  ,j+1).ne."r ") chalo(i  ,j+1) = "i "
     lhalo(i  ,j-1) = .true.; if (chalo(i  ,j-1).ne."r ") chalo(i  ,j-1) = "i "
     lhalo(i+1,j+1) = .true.; if (chalo(i+1,j+1).ne."r ") chalo(i+1,j+1) = "i "
     lhalo(i+1,j-1) = .true.; if (chalo(i+1,j-1).ne."r ") chalo(i+1,j-1) = "i "
     lhalo(i-1,j+1) = .true.; if (chalo(i-1,j+1).ne."r ") chalo(i-1,j+1) = "i " 
     lhalo(i-1,j-1) = .true.; if (chalo(i-1,j-1).ne."r ") chalo(i-1,j-1) = "i "
     
     lhalo(i+2,j  ) = .true.; if (chalo(i+2,j  ).ne."r ") chalo(i+2,j  ) = "i "
     lhalo(i-2,j  ) = .true.; if (chalo(i-2,j  ).ne."r ") chalo(i-2,j  ) = "i "
     lhalo(i  ,j+2) = .true.; if (chalo(i  ,j+2).ne."r ") chalo(i  ,j+2) = "i "
     lhalo(i  ,j-2) = .true.; if (chalo(i  ,j-2).ne."r ") chalo(i  ,j-2) = "i "
     
     do jj=j-2,j+2
        do ii=i-2,i+2
           !
           ! these points are not part of stencil
           !
           if ((ii==i+2.and.jj.ne.j).or.(ii==i+1.and.jj==j+2).or.(ii==i+1.and.jj==j-2).or.&
                (ii==i-2.and.jj.ne.j).or.(ii==i-1.and.jj==j+2).or.(ii==i-1.and.jj==j-2)) cycle
           !
           !
           !
           select case (cubeboundary)
           case(west) 
              if (ii<1) chalo(ii,jj) = "w "
           case(east) 
              if (ii>nc) chalo(ii,jj) = "e "
           case(north) 
              if (jj>nc) chalo(ii,jj) = "n "
           case(south) 
              if (jj<1) chalo(ii,jj) = "s "
           case(nwest) 
              if (jj>nc.and.ii>0  ) chalo(ii,jj) = "n "
              if (ii<1.and.jj<nc+1) chalo(ii,jj) = "w "
           case(neast) 
              if (jj>nc.and.ii<nc+1) chalo(ii,jj) = "n "
              if (ii>nc.and.jj<nc+1) chalo(ii,jj) = "e "
           case(seast) 
              if (jj<1 .and.jj<nc+1) chalo(ii,jj) = "s "
              if (ii>nc.and.jj>0   ) chalo(ii,jj) = "e "
           case(swest)
              if (jj<1.and.ii>0) chalo(ii,jj) = "s "
              if (ii<1.and.jj>0) chalo(ii,jj) = "w "
           end select
        end do
     end do
     !
     ! corner points have special stencil
     ! 
     select case (cubeboundary)     
     case(nwest) 
        if (i==1.and.j==nc) then
           lhalo( 0,nc+1) = .true.; chalo( 0,nc+1) = "wn";
           lhalo(-1,nc+1) = .true.; chalo(-1,nc+1) = "w ";
           lhalo( 0,nc+2) = .true.; chalo( 0,nc+2) = "n ";
        end if
     case(neast) 
        if (i==nc.and.j==nc) then
           lhalo(nc+1,nc+1) = .true.; chalo(nc+1,nc+1) = "ne"
           lhalo(nc+2,nc+1) = .true.; chalo(nc+2,nc+1) = "e "
           lhalo(nc+1,nc+2) = .true.; chalo(nc+1,nc+2) = "n "
        end if
     case(seast) 
        if (i==nc.and.j==1) then
           lhalo(nc+1, 0) = .true.; chalo(nc+1, 0) = "se"
           lhalo(nc+2, 0) = .true.; chalo(nc+2, 0) = "e "
           lhalo(nc+1,-1) = .true.; chalo(nc+1,-1) = "s"
        end if
     case(swest)
        if (i==1.and.j==1) then
           lhalo(0 , 0) = .true.; chalo(0 , 0) = "ws"
           lhalo(-1, 0) = .true.; chalo(-1, 0) = "w "
           lhalo(0 ,-1) = .true.; chalo(0 ,-1) = "s "
        end if
     end select
  end if
end subroutine get_stencil


end module fvm_reconstruction_mod
