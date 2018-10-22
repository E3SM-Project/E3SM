! Module: mo_cloud_sampling

! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015-2016,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Cloud state sampling for Monte Carlo Independent Column Approximation. 

module mo_cloud_sampling 
  use mo_rte_kind,      only: wp
  use mo_rng,              only: ty_rng
  implicit none 
  
  interface sample_clouds
    module procedure sample_one_cloudfrac
  end interface sample_clouds 
  
  private 
  public :: set_overlap, sample_clouds, RAN_OVERLAP, MAX_OVERLAP, MAX_RAN_OVERLAP
  
  ! Overlap kinds -- should this be an ENUM? 
  integer, parameter :: RAN_OVERLAP = 3, MAX_OVERLAP = 2, MAX_RAN_OVERLAP = 1
  integer            :: overlap = MAX_RAN_OVERLAP
contains 
  !--------------------------------------------------------------------------------------------------------------------
  function set_overlap(overlap_in) result(error_msg) 
    integer, intent(in) :: overlap_in
    character(len=128)  :: error_msg 
  
    if(.not. any(overlap_in == [RAN_OVERLAP,MAX_OVERLAP,MAX_RAN_OVERLAP])) then 
      error_msg = "cloud sampling: invalid overlap"
    else
      overlap = overlap_in
      error_msg = ""
    end if
  end function set_overlap
  !--------------------------------------------------------------------------------------------------------------------
  function sample_one_cloudfrac(ncol, nlay, ngpt, rngs, cld_frac, top_at_1, cld_mask) & 
    result(error_msg) 
    ! Input is a single cloud fraction 
    ! Output is a cloud mask
    ! Inputs and outputs have different dimensions because cloud optics is more efficient 
    !   with spectral index first
    ! See, for example, Raisanen et al., 2004, https://dx.doi.org/10.1256/qj.03.99
    
    integer,                             intent(in   ) :: ncol, nlay, ngpt
    class(ty_rng), &
              dimension(:),              intent(inout) :: rngs      !< random number generator states
                                                                    ! SHould be length ncol
    real(wp), dimension(ncol,nlay),      intent(in  )  :: cld_frac
    logical,                             intent(in   ) :: top_at_1
    ! Consider making rank a variable without column dimension 
    !   Save space but wouldn't let operations be atomic 
    logical,  dimension(ngpt,nlay,ncol), intent(  out) :: cld_mask
    character(len=128)                                 :: error_msg 
    !---------------------------------
    ! Local variables
    integer :: icol, ilay 
    real(wp), dimension(ngpt,nlay,ncol) :: rank 
    
    !---------------------------------
    error_msg = ""
    select case(overlap) 
      case default 
        error_msg = 'sample_clouds: unknown overlap specified'
      case(MAX_OVERLAP) 
        !
        ! For maximum overlap the same random number is used through the whole atmosphere 
        ! Raisanen et al., 2004, https://dx.doi.org/10.1256/qj.03.99, Eq 4. 
        !
        do icol = 1, ncol
          rank(1:ngpt,1,icol) = rngs(icol)%get_random(ngpt) 
          do ilay = 2, nlay
            rank(1:ngpt,ilay,icol) = rank(1:ngpt,1,icol)
          end do 
        end do 
      case(RAN_OVERLAP, MAX_RAN_OVERLAP) 
        ! 
        ! For random overlap each layer gets a different random number, but there's no 
        !   reason to generate the random numbers if there are no clouds 
        !   (efficiency may overrule this) 
        ! Raisanen et al., 2004, https://dx.doi.org/10.1256/qj.03.99, Eq 4. 
        ! Max-random is a perturbation on this -- see below. 
        !
        do icol = 1, ncol
          do ilay = 1, nlay
            if(cld_frac(icol,ilay) > 0._wp) then 
              rank(1:ngpt,ilay,icol) = rngs(icol)%get_random(ngpt) 
            else
              rank(1:ngpt,ilay,icol) = 0._wp
            end if 
          end do 
        end do 
    end select 
    
    if(overlap == MAX_RAN_OVERLAP) then 
      ! Walk from top down and enforce max/random
        ! Raisanen et al., 2004, https://dx.doi.org/10.1256/qj.03.99, Eq 14. 
      if(top_at_1) then
        do icol = 1, ncol
          do ilay = 2, nlay
            if(cld_frac(icol,ilay-1) > 0._wp) & 
              rank(1:ngpt,ilay,icol) = MERGE(rank(1:ngpt,ilay-1,icol), & 
                                             rank(1:ngpt,ilay  ,icol) * (1._wp - cld_frac(icol,ilay-1)), & 
                                             rank(1:ngpt,ilay-1,icol) > (1._wp - cld_frac(icol,ilay-1)))
          end do  
        end do
      else
        do icol = 1, ncol
          do ilay = nlay-1, 1, -1
            if(cld_frac(icol,ilay+1) > 0._wp) & 
              rank(1:ngpt,ilay,icol) = MERGE(rank(1:ngpt,ilay+1,icol), & 
                                             rank(1:ngpt,ilay  ,icol) * (1._wp - cld_frac(icol,ilay+1)), & 
                                             rank(1:ngpt,ilay+1,icol) > (1._wp - cld_frac(icol,ilay+1)))
          end do  
        end do
      end if 
    end if 
    
    !
    ! Convert from rank/random deviate to cloud mask
    !
    do icol = 1, ncol
      do ilay = 1, nlay
        ! Could test for cloud fraction > 1; this would save lots of comparisons. 
        cld_mask(1:ngpt,ilay,icol) = (rank(1:ngpt,ilay,icol) >= 1._wp - cld_frac(icol,ilay)) 
      end do 
    end do 
  end function sample_one_cloudfrac
  !--------------------------------------------------------------------------------------------------------------------
end module mo_cloud_sampling 

