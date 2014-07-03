module CNVerticalProfileMod
#ifdef CN
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNVerticalProfileMod
!
! !DESCRIPTION:
! Module holding routines for vertical discretization of C and N inputs into deocmposing pools
!
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_TKFRZ
  use clm_varctl  , only: iulog
  use clm_varcon, only: dzsoi_decomp
  
  implicit none
  save
  private
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: decomp_vertprofiles
  
#ifdef VERTSOILC
  logical, public :: exponential_rooting_profile = .true.
  logical, public :: pftspecific_rootingprofile = .true.
  real(r8), public :: rootprof_exp  = 3.       ! parameter for how steep the profile is for root C inputs (1/ e-folding depth) (1/m)
  real(r8), public :: surfprof_exp  = 10.      ! parameter for how steep the profile is for surface components (1/ e_folding depth) (1/m)
#endif
  
!
! !REVISION HISTORY:
! 6/27/2011 Created by C. Koven
!
!EOP
!-----------------------------------------------------------------------
contains
  
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: decomp_vertprofiles
!
! !INTERFACE:
  subroutine decomp_vertprofiles(lbp, ubp, lbc,ubc,num_soilc,filter_soilc,num_soilp,filter_soilp)
!
! !DESCRIPTION:
!
! !USES:
    use clmtype
    use clm_time_manager, only: get_step_size
    use pft2colMod, only: p2c
    use clm_varcon, only: zsoi, dzsoi, zisoi
    use clm_varpar   , only: nlevdecomp, nlevgrnd, nlevdecomp_full, maxpatch_pft
    use pftvarcon, only : rootprof_beta, noveg
    use abortutils   , only : endrun

    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp        ! pft-index bounds
    integer, intent(in) :: lbc, ubc        ! column-index bounds
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:) ! filter for soil columns
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
    !
    ! !CALLED FROM:
    ! subroutine CNDecompAlloc in module CNDecompMod.F90
    !
    ! !REVISION HISTORY:
    ! 10/5/2010: created by C. Koven to calculate vertical profiles for distributing soil and litter C and N
    !
    ! !LOCAL VARIABLES:
    ! local pointers to implicit in scalars
    !
    ! column level
    real(r8), pointer :: t_soisno(:,:)           ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: nfixation_prof(:,:)     ! (1/m) profile for N fixation additions
    real(r8), pointer :: ndep_prof(:,:)          ! (1/m) profile for N fixation additions
    integer, pointer :: altmax_lastyear_indx(:)  ! frost table depth (m)
    integer , pointer :: npfts(:)                ! number of pfts for each column
    integer , pointer :: pfti(:)                 ! beginning pft index for each column
    
    ! pft level
    integer , pointer :: ivt(:)                  ! pft vegetation type
    real(r8), pointer :: rootfr(:,:)             ! fraction of roots in each soil layer  (nlevgrnd)
    integer , pointer :: pcolumn(:)              ! pft's column index
    real(r8), pointer :: leaf_prof(:,:)          ! (1/m) profile of leaves
    real(r8), pointer :: froot_prof(:,:)         ! (1/m) profile of fine roots
    real(r8), pointer :: croot_prof(:,:)         ! (1/m) profile of coarse roots
    real(r8), pointer :: stem_prof(:,:)          ! (1/m) profile of stems
    real(r8), pointer :: wtcol(:)                ! pft weight relative to column (0-1)
    logical , pointer :: pactive(:)              ! true=>do computations on this pft (see reweightMod for details)

    ! local variables
    real(r8) :: surface_prof(1:nlevdecomp)
    real(r8) :: surface_prof_tot
    real(r8) :: rootfr_tot
    real(r8) :: cinput_rootfr(lbp:ubp, 1:nlevdecomp_full)      ! pft-native root fraction used for calculating inputs
    real(r8) :: col_cinput_rootfr(lbc:ubc, 1:nlevdecomp_full)  ! col-native root fraction used for calculating inputs
    integer  :: c, j, fc, p, fp, pi
    integer  :: alt_ind
    
    ! debugging temp variables
    real(r8) :: froot_prof_sum
    real(r8) :: croot_prof_sum
    real(r8) :: leaf_prof_sum
    real(r8) :: stem_prof_sum
    real(r8) :: ndep_prof_sum
    real(r8) :: nfixation_prof_sum
    real(r8) :: delta = 1.e-10
    character(len=32) :: subname = 'decomp_vertprofiles'

    ! assign local pointers at the column level    
    nfixation_prof                    => clm3%g%l%c%cps%nfixation_prof
    ndep_prof                         => clm3%g%l%c%cps%ndep_prof
    altmax_lastyear_indx              => clm3%g%l%c%cps%altmax_lastyear_indx
    npfts                             => clm3%g%l%c%npfts
    pfti                              => clm3%g%l%c%pfti
    
    ! assign local pointers at the pft level
    ivt                               => clm3%g%l%c%p%itype
    leaf_prof                         => clm3%g%l%c%p%pps%leaf_prof
    froot_prof                        => clm3%g%l%c%p%pps%froot_prof
    croot_prof                        => clm3%g%l%c%p%pps%croot_prof
    stem_prof                         => clm3%g%l%c%p%pps%stem_prof
    pcolumn                           => clm3%g%l%c%p%column
    rootfr                            => clm3%g%l%c%p%pps%rootfr
    wtcol                             => clm3%g%l%c%p%wtcol
    pactive                           => clm3%g%l%c%p%active


#ifdef VERTSOILC
    ! define a single shallow surface profile for surface additions (leaves, stems, and N deposition)
    surface_prof(:) = 0._r8
    do j = 1, nlevdecomp
       surface_prof(j) = exp(-surfprof_exp * zsoi(j)) / dzsoi_decomp(j)
    end do
    
    ! initialize profiles to zero
    leaf_prof(:,:) = 0._r8
    froot_prof(:,:) = 0._r8
    croot_prof(:,:) = 0._r8
    stem_prof(:,:) = 0._r8
    nfixation_prof(:,:) = 0._r8
    ndep_prof(:,:) = 0._r8

    cinput_rootfr(:,:) = 0._r8
    col_cinput_rootfr(:,:) = 0._r8
    
    if ( exponential_rooting_profile ) then
       if ( .not. pftspecific_rootingprofile ) then
          ! define rooting profile from exponential parameters
          do j = 1, nlevdecomp
             do fp = 1,num_soilp
                p = filter_soilp(fp)
                cinput_rootfr(p,j) = exp(-rootprof_exp * zsoi(j)) / dzsoi_decomp(j)
             end do
          end do
       else
          ! use beta distribution parameter from Jackson et al., 1996
          do p = lbp, ubp
             c = pcolumn(p)
             if (ivt(p) /= noveg) then
                do j = 1, nlevdecomp
                   cinput_rootfr(p,j) = ( rootprof_beta(ivt(p)) ** (zisoi(j-1)*100._r8) - rootprof_beta(ivt(p)) ** (zisoi(j)*100._r8) ) &
                        / dzsoi_decomp(j)
                end do
             else
                cinput_rootfr(p,1) = 1._r8 / dzsoi_decomp(1)
             endif
          end do
       endif
    else
       do j = 1, nlevdecomp
          ! use standard CLM root fraction profiles
          do fp = 1,num_soilp
             p = filter_soilp(fp)
             cinput_rootfr(p,j) = rootfr(p,j) / dzsoi_decomp(j)
          end do
       end do
    endif
    
    do fp = 1,num_soilp
       p = filter_soilp(fp)
       c = pcolumn(p)
       ! integrate rootfr over active layer of soil column
       rootfr_tot = 0._r8
       surface_prof_tot = 0._r8
       do j = 1, min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
          rootfr_tot = rootfr_tot + cinput_rootfr(p,j) * dzsoi_decomp(j)
          surface_prof_tot = surface_prof_tot + surface_prof(j)  * dzsoi_decomp(j)
       end do
       if ( (altmax_lastyear_indx(c) .gt. 0) .and. (rootfr_tot .gt. 0._r8) .and. (surface_prof_tot .gt. 0._r8) ) then
          ! where there is not permafrost extending to the surface, integrate the profiles over the active layer
          ! this is equivalnet to integrating over all soil layers outside of permafrost regions
          do j = 1, min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
             froot_prof(p,j) = cinput_rootfr(p,j) / rootfr_tot
             croot_prof(p,j) = cinput_rootfr(p,j) / rootfr_tot
             ! set all surface processes to shallower profile
             leaf_prof(p,j) = surface_prof(j)/ surface_prof_tot
             stem_prof(p,j) = surface_prof(j)/ surface_prof_tot
          end do
       else
          ! if fully frozen, or no roots, put everything in the top layer
          froot_prof(p,1) = 1./dzsoi_decomp(1)
          croot_prof(p,1) = 1./dzsoi_decomp(1)
          leaf_prof(p,1) = 1./dzsoi_decomp(1)
          stem_prof(p,1) = 1./dzsoi_decomp(1)
       endif
       
    end do

    !! aggregate root profile to column
    ! call p2c (lbp, ubp, lbc, ubc, nlevdecomp_full, cinput_rootfr, col_cinput_rootfr, 'unity')
    do pi = 1,maxpatch_pft
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          if (pi <=  npfts(c)) then
             p = pfti(c) + pi - 1
             if (pactive(p)) then
                do j = 1,nlevdecomp
                   col_cinput_rootfr(c,j) = col_cinput_rootfr(c,j) + cinput_rootfr(p,j) * wtcol(p)
                end do
             end if
          end if
       end do
    end do

    ! repeat for column-native profiles: Ndep and Nfix
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       rootfr_tot = 0._r8
       surface_prof_tot = 0._r8
       ! redo column ntegration over active layer for column-native profiles
       do j = 1, min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
          rootfr_tot = rootfr_tot + col_cinput_rootfr(c,j) * dzsoi_decomp(j)
          surface_prof_tot = surface_prof_tot + surface_prof(j) * dzsoi_decomp(j)
       end do
       if ( (altmax_lastyear_indx(c) .gt. 0) .and. (rootfr_tot .gt. 0._r8) .and. (surface_prof_tot .gt. 0._r8) ) then
          do j = 1,  min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
             nfixation_prof(c,j) = col_cinput_rootfr(c,j) / rootfr_tot
             ndep_prof(c,j) = surface_prof(j)/ surface_prof_tot
          end do
       else
          nfixation_prof(c,1) = 1./dzsoi_decomp(1)
          ndep_prof(c,1) = 1./dzsoi_decomp(1)
       endif
    end do
    
#else
    
    ! for one layer decomposition model, set profiles to unity
    leaf_prof(:,:) = 1._r8
    froot_prof(:,:) = 1._r8
    croot_prof(:,:) = 1._r8
    stem_prof(:,:) = 1._r8
    nfixation_prof(:,:) = 1._r8
    ndep_prof(:,:) = 1._r8
    
#endif
    

    ! check to make sure integral of all profiles = 1.
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       ndep_prof_sum = 0.
       nfixation_prof_sum = 0.
       do j = 1, nlevdecomp
          ndep_prof_sum = ndep_prof_sum + ndep_prof(c,j) *  dzsoi_decomp(j)
          nfixation_prof_sum = nfixation_prof_sum + nfixation_prof(c,j) *  dzsoi_decomp(j)
       end do
       if ( ( abs(ndep_prof_sum - 1._r8) .gt. delta ) .or.  ( abs(nfixation_prof_sum - 1._r8) .gt. delta ) ) then
          write(iulog, *) 'profile sums: ', ndep_prof_sum, nfixation_prof_sum
          write(iulog, *) 'c: ', c
          write(iulog, *) 'altmax_lastyear_indx: ', altmax_lastyear_indx(c)
          write(iulog, *) 'nfixation_prof: ', nfixation_prof(c,:)
          write(iulog, *) 'ndep_prof: ', ndep_prof(c,:)
          write(iulog, *) 'cinput_rootfr: ', cinput_rootfr(c,:)
          write(iulog, *) 'dzsoi_decomp: ', dzsoi_decomp(:)
          write(iulog, *) 'surface_prof: ', surface_prof(:)
          write(iulog, *) 'npfts(c): ', npfts(c)
          do p = pfti(c), pfti(c) + npfts(c) -1
             write(iulog, *) 'p, ivt(p), wtcol(p): ', p, ivt(p), wtcol(p)
             write(iulog, *) 'cinput_rootfr(p,:): ', cinput_rootfr(p,:)
          end do
          call endrun( trim(subname)//" ERROR: _prof_sum-1>delta" )
       endif
    end do

    do fp = 1,num_soilp
       p = filter_soilp(fp)
       froot_prof_sum = 0.
       croot_prof_sum = 0.
       leaf_prof_sum = 0.
       stem_prof_sum = 0.
       do j = 1, nlevdecomp
          froot_prof_sum = froot_prof_sum + froot_prof(p,j) *  dzsoi_decomp(j)
          croot_prof_sum = croot_prof_sum + croot_prof(p,j) *  dzsoi_decomp(j)
          leaf_prof_sum = leaf_prof_sum + leaf_prof(p,j) *  dzsoi_decomp(j)
          stem_prof_sum = stem_prof_sum + stem_prof(p,j) *  dzsoi_decomp(j)
       end do
       if ( ( abs(froot_prof_sum - 1._r8) .gt. delta ) .or.  ( abs(croot_prof_sum - 1._r8) .gt. delta ) .or. &
            ( abs(stem_prof_sum - 1._r8) .gt. delta ) .or.  ( abs(leaf_prof_sum - 1._r8) .gt. delta ) ) then
          write(iulog, *) 'profile sums: ', froot_prof_sum, croot_prof_sum, leaf_prof_sum, stem_prof_sum
          call endrun( trim(subname)//' ERROR: sum-1 > delta' )
       endif
    end do


  end subroutine decomp_vertprofiles
  
    
#endif
  
end module CNVerticalProfileMod
