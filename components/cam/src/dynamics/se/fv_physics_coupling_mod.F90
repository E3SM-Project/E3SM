!---------------------------------------------------------------------------------------------------
! This module contains routines for mapping between dynamics and physics columns
! when dynamics is on the spectral element grid (i.e. GLL) and the physics is on
! a quasi-equal area finite volume grid that evenly divides the dyn element
! 
! The current mapping method utilizes a simple piece-wise linear approach where
! the physics tendencies are uniformly copied to the underlying GLL nodes. This
! comes with the restriction of only using 1 or 2x2 phsyics cells per element, 
! but it also carries the advantage of not needing any information outside of 
! the element. Higher order mapping method can be implemented, but they will 
! require special consideration of the model's vertical coordinate, which is 
! currently based on moist pressure. In either case, dyn_to_fv_phys() would stay
! the same, with simple weighted averaging over element subcells.
! 
! Author: Walter Hannah (LLNL)
!---------------------------------------------------------------------------------------------------
module fv_physics_coupling_mod
  use element_mod,    only: element_t
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use kinds,          only: real_kind, int_kind
  use constituents,   only: pcnst, cnst_name
  use dimensions_mod, only: np, npsq, nelemd, nlev, fv_nphys
  use ppgrid,         only: pcols, pver, pverp
    
  private
  ! These method encapsulate the coupling method for the fv phys grid,
  ! they are used by d_p_coupling() and p_d_coupling() in dp_coupling.F90 
  public :: fv_phys_to_dyn
  public :: dyn_to_fv_phys

contains
  !=================================================================================================
  !=================================================================================================
  subroutine fv_phys_to_dyn(elem,T_tmp,uv_tmp,Q_tmp)
    implicit none
    !---------------------------------------------------------------------------
    ! interface arguments
    type(element_t),      intent(inout) :: elem(:)          ! dynamics element structure
    real(kind=real_kind), intent(inout) :: T_tmp (:,:,:)    ! temp array to hold T
    real(kind=real_kind), intent(inout) :: uv_tmp(:,:,:,:)  ! temp array to hold u and v
    real(kind=real_kind), intent(inout) :: q_tmp (:,:,:,:)  ! temp to hold advected constituents
    ! local variables
    integer(kind=int_kind)   :: ie, m, i, j, icol, ilyr    ! loop iterators
    integer(kind=int_kind)   :: gi(2), gj(2)               ! index list used to simplify pg2 case
    integer(kind=int_kind)   :: di, dj
    !---------------------------------------------------------------------------
    ! Copy tendencies on the physics grid over to the dynamics grid (GLL)
    !---------------------------------------------------------------------------
    do ie = 1,nelemd
      icol = 0
      do j = 1,fv_nphys
        do i = 1,fv_nphys 
          icol = icol + 1
          do ilyr = 1,pver
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
            ! the pg1 case is simple, not sure how to generalize with pg2
            if (fv_nphys == 1) then
              elem(ie)%derived%FT(:,:,  ilyr)   = T_tmp (icol,  ilyr,ie)
              elem(ie)%derived%FM(:,:,1,ilyr)   = uv_tmp(icol,1,ilyr,ie)
              elem(ie)%derived%FM(:,:,2,ilyr)   = uv_tmp(icol,2,ilyr,ie)
              do m = 1,pcnst
                elem(ie)%derived%FQ(:,:,ilyr,m) = q_tmp (icol,  ilyr,m,ie)
              end do
            end if ! fv_nphys == 1
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
            ! for pg2 we need to copy the FV state to quadrants of GLL grid
            if (fv_nphys == 2) then
              ! define indices of GLL points in quadrant
              di = (i-1)+(i-2)  ! either i=1 & di=-1 or i=2 & di=+1
              dj = (j-1)+(j-2)  ! either j=1 & dj=-1 or j=2 & dj=+1
              ! gi & gj are indices of the 4 GLL nodes in each element quadrant
              gi(1:2) = (/i+1,i+1+di/)
              gj(1:2) = (/j+1,j+1+dj/)
              ! copy physics column values to GLL nodes
              elem(ie)%derived%FT(gi,gj,ilyr)   = T_tmp(icol,   ilyr,ie)
              elem(ie)%derived%FM(gi,gj,1,ilyr) = uv_tmp(icol,1,ilyr,ie)
              elem(ie)%derived%FM(gi,gj,2,ilyr) = uv_tmp(icol,2,ilyr,ie)
              do m = 1,pcnst
                elem(ie)%derived%FQ(gi,gj,ilyr,m) = q_tmp(icol,ilyr,m,ie)
              end do
            end if ! fv_nphys == 2
            !-------------------------------------------------------------------
            !-------------------------------------------------------------------
          end do ! ilyr
        end do ! i
      end do ! j
    end do ! ie

  end subroutine fv_phys_to_dyn
  !================================================================================================
  !================================================================================================
  subroutine dyn_to_fv_phys(elem,ps_tmp,zs_tmp,T_tmp,uv_tmp,w_tmp,Q_tmp)
    use derivative_mod,     only: subcell_integration
    use dyn_comp,           only: TimeLevel
    implicit none
    !---------------------------------------------------------------------------
    ! interface arguments
    type(element_t),      intent(inout) :: elem(:)          ! dynamics element structure
    real(kind=real_kind), intent(inout) :: ps_tmp(:,:)      ! temp array to hold ps
    real(kind=real_kind), intent(inout) :: zs_tmp(:,:)      ! temp array to hold phis  
    real(kind=real_kind), intent(inout) :: T_tmp (:,:,:)    ! temp array to hold T
    real(kind=real_kind), intent(inout) :: uv_tmp(:,:,:,:)  ! temp array to hold u and v
    real(kind=real_kind), intent(inout) :: w_tmp (:,:,:)    ! temp array to hold omega
    real(kind=real_kind), intent(inout) :: q_tmp (:,:,:,:)  ! temp to hold advected constituents
    ! local variables
    integer(kind=int_kind) :: ie, m, icol, ilyr             ! loop iterators
    integer                :: tl_f                          ! time level
    integer                :: nphys_sq
    real(r8), dimension(np,np)             :: tmp_area
    real(r8), dimension(fv_nphys,fv_nphys) :: inv_area
    real(r8), dimension(fv_nphys*fv_nphys) :: inv_area_reshape
    real(r8), dimension(np,np)             :: dp_gll
    real(r8), dimension(fv_nphys,fv_nphys) :: dp_fvm
    real(r8), dimension(fv_nphys*fv_nphys) :: inv_dp_fvm_reshape
    !---------------------------------------------------------------------------
    ! Integrate dynamics field with appropriate weighting 
    ! to get average state in each physics cell
    !---------------------------------------------------------------------------
    tl_f = TimeLevel%n0
    nphys_sq = fv_nphys*fv_nphys
    tmp_area(:,:) = 1.0_r8

    do ie = 1,nelemd

      inv_area(:,:) = 1.0_r8/subcell_integration(tmp_area,np,fv_nphys,elem(ie)%metdet(:,:))
      inv_area_reshape(:) = RESHAPE( inv_area(:,:), (/nphys_sq/) )

      ps_tmp(:,ie) = RESHAPE( subcell_integration(                  &
                     elem(ie)%state%ps_v(:,:,tl_f),                 &
                     np, fv_nphys, elem(ie)%metdet(:,:) ) ,         &
                     (/nphys_sq/)  ) * inv_area_reshape
      zs_tmp(:,ie) = RESHAPE( subcell_integration(                  &
                     elem(ie)%state%phis(:,:),                      &
                     np, fv_nphys, elem(ie)%metdet(:,:) ) ,         &
                     (/nphys_sq/) ) * inv_area_reshape

      do ilyr = 1,pver

        dp_gll = elem(ie)%state%dp3d(:,:,ilyr,tl_f)
        dp_fvm = subcell_integration(dp_gll,np,fv_nphys,elem(ie)%metdet(:,:))
        inv_dp_fvm_reshape = 1.0 / RESHAPE( dp_fvm , (/nphys_sq/) )

        T_tmp(:,ilyr,ie)      = RESHAPE( subcell_integration(             &
                                elem(ie)%state%T(:,:,ilyr,tl_f)*dp_gll,   &
                                np, fv_nphys, elem(ie)%metdet(:,:) ) ,    &
                                (/nphys_sq/) ) *inv_dp_fvm_reshape

        w_tmp(:,ilyr,ie)      = RESHAPE( subcell_integration(             &
                                elem(ie)%derived%omega_p(:,:,ilyr),       &
                                np, fv_nphys, elem(ie)%metdet(:,:) ) ,    &
                                (/nphys_sq/) ) *inv_area_reshape
        do m = 1,2
          uv_tmp(:,m,ilyr,ie) = RESHAPE( subcell_integration(             &
                                elem(ie)%state%V(:,:,m,ilyr,tl_f),        &
                                np, fv_nphys, elem(ie)%metdet(:,:) ) ,    &
                                (/nphys_sq/) ) *inv_area_reshape
        end do
        do m = 1,pcnst
          Q_tmp(:,ilyr,m,ie)  = RESHAPE( subcell_integration(             &
                                elem(ie)%state%Q(:,:,ilyr,m)*dp_gll,      &
                                np, fv_nphys, elem(ie)%metdet(:,:) ) ,    &
                                (/nphys_sq/) ) *inv_dp_fvm_reshape
        end do

      end do ! ilyr
    end do ! ie
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
  end subroutine dyn_to_fv_phys
  !=================================================================================================
  !=================================================================================================
end module fv_physics_coupling_mod