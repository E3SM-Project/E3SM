module VerticalProfileMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines for vertical discretization of C and N inputs into deocmposing pools
  !
  ! !USES:
  use shr_kind_mod    , only: r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use subgridAveMod   , only : p2c
  use SoilStateType   , only : soilstate_type
  use CanopyStateType , only : canopystate_type
  use CNStateType     , only : cnstate_type
  use ColumnType      , only : col_pp
  use VegetationType  , only : veg_pp
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: decomp_vertprofiles
  public :: decomp_vertprofiles_gpu
  !
  logical , public :: exponential_rooting_profile = .true.
  logical , public :: pftspecific_rootingprofile  = .true.
  ! how steep profile is for root C inputs (1/ e-folding depth) (1/m)
  real(r8), public :: rootprof_exp  = 3.
  ! how steep profile is for surface components (1/ e_folding depth) (1/m)
  real(r8), public :: surfprof_exp  = 10.
  !$acc declare copyin(rootprof_exp,surfprof_exp)
  !$acc declare copyin(exponential_rooting_profile, pftspecific_rootingprofile)
  !-----------------------------------------------------------------------

contains

  subroutine decomp_vertprofiles(bounds, &
       num_soilc,filter_soilc,num_soilp,filter_soilp, &
       soilstate_vars, canopystate_vars, cnstate_vars)
    !
    ! !DESCRIPTION:
    !  calculate vertical profiles for distributing soil and litter C and N
    !
    !  Note (WJS, 6-12-13): Because of this routine's placement in the driver sequence (it is
    !  called very early in each timestep, before weights are adjusted and filters are
    !  updated), it may be necessary for this routine to compute values over inactive as well
    !  as active points (since some inactive points may soon become active) - so that's what
    !  is done now. Currently, it seems to be okay to do this, because the variables computed
    !  here seem to only depend on quantities that are valid over inactive as well as active
    !  points. However, note that this routine is (mistakenly) called from two places
    !  currently - the above note applies to its call from the driver, but its call from
    !  SoilLittDecompMod uses the standard filters that just apply over active points
    !
    ! !USES:
    use elm_varcon  , only : zsoi, dzsoi, zisoi, dzsoi_decomp
    use elm_varpar  , only : nlevdecomp, nlevgrnd, nlevdecomp_full, maxpatch_pft
    use elm_varctl  , only : use_vertsoilc, iulog, use_dynroot
    use pftvarcon   , only : rootprof_beta, noveg
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilstate_type)   , intent(in)    :: soilstate_vars
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(cnstate_type)     , intent(inout) :: cnstate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8) :: surface_prof(1:nlevdecomp)!, filter_wtpatch(1:num_soilp)
    real(r8) :: surface_prof_tot
    real(r8) :: rootfr_tot
    real(r8) :: cinput_rootfr(bounds%begp:bounds%endp, 1:nlevdecomp_full)      ! pft-native root fraction used for calculating inputs
    real(r8) :: col_cinput_rootfr(1:num_soilc, 1:nlevdecomp_full)  ! col-native root fraction used for calculating inputs
    integer  :: c, j, fc, p, fp, pi
    integer  :: alt_ind
    ! debugging temp variables
    real(r8) :: froot_prof_sum
    real(r8) :: croot_prof_sum
    real(r8) :: leaf_prof_sum
    real(r8) :: stem_prof_sum
    real(r8) :: ndep_prof_sum
    real(r8) :: nfixation_prof_sum
    real(r8) :: pdep_prof_sum
    real(r8) :: temp_sum
    real(r8) :: delta = 1.e-10
    logical, parameter :: debug=.false.
    real :: startt, stopt
    character(len=32) :: subname = 'decomp_vertprofiles_gpu'
    !-----------------------------------------------------------------------
    associate(                                                               &
         rootfr               => soilstate_vars%rootfr_patch               , & ! Input:  [real(r8)  (:,:) ]  fraction of roots in each soil layer  (nlevgrnd)

         altmax_lastyear_indx => canopystate_vars%altmax_lastyear_indx_col , & ! Input:  [integer   (:)   ]  frost table depth (m)

         nfixation_prof       => cnstate_vars%nfixation_prof_col           , & ! Input:  [real(r8)  (:,:) ]  (1/m) profile for N fixation additions
         ndep_prof            => cnstate_vars%ndep_prof_col                , & ! Input:  [real(r8)  (:,:) ]  (1/m) profile for N fixation additions
         pdep_prof            => cnstate_vars%pdep_prof_col                , & ! Input:  [real(r8)  (:,:) ]  (1/m) profile for P depostition additions

         leaf_prof            => cnstate_vars%leaf_prof_patch              , & ! Output:  [real(r8) (:,:) ]  (1/m) profile of leaves
         froot_prof           => cnstate_vars%froot_prof_patch             , & ! Output:  [real(r8) (:,:) ]  (1/m) profile of fine roots
         croot_prof           => cnstate_vars%croot_prof_patch             , & ! Output:  [real(r8) (:,:) ]  (1/m) profile of coarse roots
         stem_prof            => cnstate_vars%stem_prof_patch              , & ! Output:  [real(r8) (:,:) ]  (1/m) profile of stems

         begp                 => bounds%begp                               , &
         endp                 => bounds%endp                               , &
         begc                 => bounds%begc                               , &
         endc                 => bounds%endc                                 &
         )

      if (use_vertsoilc) then

         call cpu_time(startt)
         do j = 1, nlevdecomp
            surface_prof(j) = exp(-surfprof_exp * zsoi(j))/dzsoi_decomp(j)
         end do


         call cpu_time(stopt)
         print *, "VerticalProfile::init",(stopt-startt)*1.e+3,"ms"
         call cpu_time(startt)
         if ( exponential_rooting_profile ) then
            if ( .not. pftspecific_rootingprofile ) then
               ! define rooting profile from exponential parameters
               do j = 1, nlevdecomp
                  do fp = 1,num_soilp
                     cinput_rootfr(p,j) = exp(-rootprof_exp * zsoi(j)) / dzsoi_decomp(j)
                  end do
               end do
            else
               ! use beta distribution parameter from Jackson et al., 1996
               do j = 1, nlevdecomp
                  do fp = 1,num_soilp
                     p = filter_soilp(fp)
                     if (veg_pp%itype(p) /= noveg) then
                        cinput_rootfr(p,j) = ( rootprof_beta(veg_pp%itype(p)) ** (zisoi(j-1)*100._r8) - &
                             rootprof_beta(veg_pp%itype(p)) ** (zisoi(j)*100._r8) ) &
                             / dzsoi_decomp(j)
                     else
                        cinput_rootfr(p,1) = 1._r8 / dzsoi_decomp(1)
                     endif
                  end do
               end do
            endif
            call cpu_time(stopt)
            print *, "VerticleProfile::cinput_rootfr",(stopt-startt)*1.e+3,"ms"
         else
            call cpu_time(startt)
            do j = 1, nlevdecomp
               ! use standard CLM root fraction profiles
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  cinput_rootfr(p,j) = rootfr(p,j) / dzsoi_decomp(j)
               end do
            end do
            call cpu_time(stopt)
            print *, "VerticleProfile::cinput_rootfr",(stopt-startt)*1.e+3,"ms"
         endif
         call cpu_time(startt)
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = veg_pp%column(p)
            ! integrate rootfr over active layer of soil column
            rootfr_tot = 0._r8
            surface_prof_tot = 0._r8
            do j = 1, min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
               rootfr_tot = rootfr_tot + cinput_rootfr(p,j) * dzsoi_decomp(j)
               surface_prof_tot = surface_prof_tot + surface_prof(j)  * dzsoi_decomp(j)
            end do

            if ( (altmax_lastyear_indx(c) > 0) .and. (rootfr_tot > 0._r8) .and. (surface_prof_tot > 0._r8) ) then
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
               leaf_prof (p,1) = 1./dzsoi_decomp(1)
               stem_prof (p,1) = 1./dzsoi_decomp(1)
            endif

         end do
         call cpu_time(stopt)
         print *, "VerticleProfile::Surface reduction",(stopt-startt)*1.e+3,"ms"

         !! aggregate root profile to column
         ! call p2c (decomp, nlevdecomp_full, &
         !      cinput_rootfr(bounds%begp:bounds%endp, :), &
         !      col_cinput_rootfr(bounds%begc:bounds%endc, :), &
         !      'unity')

         call cpu_time(startt)
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               temp_sum = 0._r8
               do p = col_pp%pfti(c),col_pp%pftf(c)
                  temp_sum = temp_sum + cinput_rootfr(p,j) * veg_pp%wtcol(p)
               end do
               col_cinput_rootfr(c,j) = temp_sum
            end do
         end do
         call cpu_time(stopt)
         print *, "VerticleProfile::col_rootfr",(stopt-startt)*1.e+3,"ms"

         call cpu_time(startt)
         ! repeat for column-native profiles: Ndep and Nfix
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            rootfr_tot = 0._r8
            surface_prof_tot = 0._r8
            ! redo column ntegration over active layer for column-native profiles
            do j = 1, min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
               rootfr_tot = rootfr_tot + col_cinput_rootfr(fc,j) * dzsoi_decomp(j)
               surface_prof_tot = surface_prof_tot + surface_prof(j) * dzsoi_decomp(j)
            end do
            if ( (altmax_lastyear_indx(c) > 0) .and. (rootfr_tot > 0._r8) .and. (surface_prof_tot > 0._r8) ) then
               do j = 1,  min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
                  nfixation_prof(c,j) = col_cinput_rootfr(fc,j) / rootfr_tot
                  ndep_prof(c,j) = surface_prof(j)/ surface_prof_tot
                  pdep_prof(c,j) = surface_prof(j)/ surface_prof_tot
               end do
            else
               nfixation_prof(c,1) = 1./dzsoi_decomp(1)
               ndep_prof(c,1) = 1./dzsoi_decomp(1)
               pdep_prof(c,1) = 1./dzsoi_decomp(1)
            endif
         end do
         call cpu_time(stopt)
         print *, "VerticleProfile::NdepNfix",(stopt-startt)*1.E+3,"ms"

      else

         ! for one layer decomposition model, set profiles to unity
         do j = 1, nlevdecomp
            do fp = 1, num_soilp
               p = filter_soilp(fp)
               leaf_prof(p, j) = 1._r8
               froot_prof(p,j) = 1._r8
               croot_prof(p,j) = 1._r8
               stem_prof(p,j) = 1._r8
            end do
         end do

         do j = 1, nlevdecomp
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               nfixation_prof(c,j) = 1._r8
               ndep_prof(c,j) = 1._r8
               pdep_prof(c,j) = 1._r8
            end do
         end do

      end if

      call cpu_time(startt)
      ! check to make sure integral of all profiles = 1.
      if(debug) then
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         ndep_prof_sum = 0._r8
         nfixation_prof_sum = 0._r8
         pdep_prof_sum = 0._r8
         do j = 1, nlevdecomp
            ndep_prof_sum = ndep_prof_sum + ndep_prof(c,j) *  dzsoi_decomp(j)
            nfixation_prof_sum = nfixation_prof_sum + nfixation_prof(c,j) *  dzsoi_decomp(j)
            pdep_prof_sum = pdep_prof_sum + pdep_prof(c,j) *  dzsoi_decomp(j)
         end do
         print *, "ndep_prof_sum",ndep_prof_sum
         print *, "nfixation_prof_sum",nfixation_prof_sum
         print *, "pdep_prof_sum", pdep_prof_sum

      end do

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         froot_prof_sum = 0._r8
         croot_prof_sum = 0._r8
         leaf_prof_sum = 0._r8
         stem_prof_sum = 0._r8
         do j = 1, nlevdecomp
            froot_prof_sum = froot_prof_sum + froot_prof(p,j) *  dzsoi_decomp(j)
            croot_prof_sum = croot_prof_sum + croot_prof(p,j) *  dzsoi_decomp(j)
            leaf_prof_sum  = leaf_prof_sum + leaf_prof(p,j) *  dzsoi_decomp(j)
            stem_prof_sum  = stem_prof_sum + stem_prof(p,j) *  dzsoi_decomp(j)
         end do
         print *, "froot_prof_sum",froot_prof_sum
         print *, "croot_prof_sum",croot_prof_sum
         print *, "leaf_prof_sum ",leaf_prof_sum
         print *, "stem_prof_sum ",stem_prof_sum
      end do

      call cpu_time(stopt)
      print *, "VerticleProfile::Debug",(stopt-startt)*1.e+3,"ms"
      end if
    end associate

  end subroutine decomp_vertprofiles


  !-----------------------------------------------------------------------
  subroutine decomp_vertprofiles_gpu(bounds, &
       num_soilc,filter_soilc,num_soilp,filter_soilp, &
       soilstate_vars, canopystate_vars, cnstate_vars)
    !
    ! !DESCRIPTION:
    !  calculate vertical profiles for distributing soil and litter C and N
    !
    !  Note (WJS, 6-12-13): Because of this routine's placement in the driver sequence (it is
    !  called very early in each timestep, before weights are adjusted and filters are
    !  updated), it may be necessary for this routine to compute values over inactive as well
    !  as active points (since some inactive points may soon become active) - so that's what
    !  is done now. Currently, it seems to be okay to do this, because the variables computed
    !  here seem to only depend on quantities that are valid over inactive as well as active
    !  points. However, note that this routine is (mistakenly) called from two places
    !  currently - the above note applies to its call from the driver, but its call from
    !  SoilLittDecompMod uses the standard filters that just apply over active points
    !
    ! !USES:
    use elm_varcon  , only : zsoi, dzsoi, zisoi, dzsoi_decomp
    use elm_varpar  , only : nlevdecomp, nlevgrnd, nlevdecomp_full, maxpatch_pft
    use elm_varctl  , only : use_vertsoilc, iulog, use_dynroot
    use pftvarcon   , only : rootprof_beta, noveg
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilstate_type)   , intent(in)    :: soilstate_vars
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(cnstate_type)     , intent(inout) :: cnstate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8) :: surface_prof(1:nlevdecomp)!, filter_wtpatch(1:num_soilp)
    real(r8) :: surface_prof_tot
    real(r8) :: rootfr_tot
    real(r8) :: cinput_rootfr(bounds%begp:bounds%endp, 1:nlevdecomp_full)      ! pft-native root fraction used for calculating inputs
    real(r8) :: col_cinput_rootfr(1:num_soilc, 1:nlevdecomp_full)  ! col-native root fraction used for calculating inputs
    integer  :: c, j, fc, p, fp, pi
    integer  :: alt_ind
    ! debugging temp variables
    real(r8) :: froot_prof_sum
    real(r8) :: croot_prof_sum
    real(r8) :: leaf_prof_sum
    real(r8) :: stem_prof_sum
    real(r8) :: ndep_prof_sum
    real(r8) :: nfixation_prof_sum
    real(r8) :: pdep_prof_sum
    real(r8) :: temp_sum
    real(r8) :: delta = 1.e-10
    logical, parameter :: debug=.false.
    real :: startt, stopt
    character(len=32) :: subname = 'decomp_vertprofiles_gpu'
    !-----------------------------------------------------------------------
    associate(                                                               &
         rootfr               => soilstate_vars%rootfr_patch               , & ! Input:  [real(r8)  (:,:) ]  fraction of roots in each soil layer  (nlevgrnd)

         altmax_lastyear_indx => canopystate_vars%altmax_lastyear_indx_col , & ! Input:  [integer   (:)   ]  frost table depth (m)

         nfixation_prof       => cnstate_vars%nfixation_prof_col           , & ! Input:  [real(r8)  (:,:) ]  (1/m) profile for N fixation additions
         ndep_prof            => cnstate_vars%ndep_prof_col                , & ! Input:  [real(r8)  (:,:) ]  (1/m) profile for N fixation additions
         pdep_prof            => cnstate_vars%pdep_prof_col                , & ! Input:  [real(r8)  (:,:) ]  (1/m) profile for P depostition additions

         leaf_prof            => cnstate_vars%leaf_prof_patch              , & ! Output:  [real(r8) (:,:) ]  (1/m) profile of leaves
         froot_prof           => cnstate_vars%froot_prof_patch             , & ! Output:  [real(r8) (:,:) ]  (1/m) profile of fine roots
         croot_prof           => cnstate_vars%croot_prof_patch             , & ! Output:  [real(r8) (:,:) ]  (1/m) profile of coarse roots
         stem_prof            => cnstate_vars%stem_prof_patch              , & ! Output:  [real(r8) (:,:) ]  (1/m) profile of stems

         begp                 => bounds%begp                               , &
         endp                 => bounds%endp                               , &
         begc                 => bounds%begc                               , &
         endc                 => bounds%endc                                 &
         )

      if (use_vertsoilc) then

         call cpu_time(startt)
         do j = 1, nlevdecomp
            surface_prof(j) = exp(-surfprof_exp * zsoi(j))/dzsoi_decomp(j)
         end do
         !$acc enter data copyin(surface_prof(:)) create(cinput_rootfr(begp:endp,1:nlevdecomp_full),&
         !$acc                  col_cinput_rootfr(1:num_soilc,1:nlevdecomp_full))


         call cpu_time(stopt)
         print *, "VerticalProfile::init",(stopt-startt)*1.e+3,"ms"
         ! ! initialize profiles to zero
         ! leaf_prof(begp:endp, :)      = 0._r8
         ! froot_prof(begp:endp, :)     = 0._r8
         ! croot_prof(begp:endp, :)     = 0._r8
         ! stem_prof(begp:endp, :)      = 0._r8
         ! nfixation_prof(begc:endc, :) = 0._r8
         ! ndep_prof(begc:endc, :)      = 0._r8
         ! pdep_prof(begc:endc, :)      = 0._r8
         !
         ! cinput_rootfr(begp:endp, :)     = 0._r8
         ! col_cinput_rootfr(begc:endc, :) = 0._r8
         call cpu_time(startt)
         if ( exponential_rooting_profile ) then
            if ( .not. pftspecific_rootingprofile ) then
               ! define rooting profile from exponential parameters
               !$acc parallel loop independent gang default(present)
               do j = 1, nlevdecomp
                  !$acc loop vector independent private(p)
                  do fp = 1,num_soilp
                     cinput_rootfr(p,j) = exp(-rootprof_exp * zsoi(j)) / dzsoi_decomp(j)
                  end do
               end do
            else
               ! use beta distribution parameter from Jackson et al., 1996
               !$acc parallel loop independent gang default(present)
               do j = 1, nlevdecomp
                  !$acc loop vector independent private(p)
                  do fp = 1,num_soilp
                     p = filter_soilp(fp)
                     if (veg_pp%itype(p) /= noveg) then
                        cinput_rootfr(p,j) = ( rootprof_beta(veg_pp%itype(p)) ** (zisoi(j-1)*100._r8) - &
                             rootprof_beta(veg_pp%itype(p)) ** (zisoi(j)*100._r8) ) &
                             / dzsoi_decomp(j)
                     else
                        cinput_rootfr(p,1) = 1._r8 / dzsoi_decomp(1)
                     endif
                  end do
               end do
            endif
            call cpu_time(stopt)
            print *, "VerticleProfile::cinput_rootfr",(stopt-startt)*1.e+3,"ms"
         else
            call cpu_time(startt)
            !$acc parallel loop independent gang default(present)
            do j = 1, nlevdecomp
               ! use standard CLM root fraction profiles
               !$acc loop vector independent private(p)
               do fp = 1,num_soilp
                  p = filter_soilp(fp)
                  cinput_rootfr(p,j) = rootfr(p,j) / dzsoi_decomp(j)
               end do
            end do
            call cpu_time(stopt)
            print *, "VerticleProfile::cinput_rootfr",(stopt-startt)*1.e+3,"ms"
         endif
         call cpu_time(startt)
         !$acc parallel loop independent gang worker default(present) private(p,c,rootfr_tot,surface_prof_tot)
         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = veg_pp%column(p)
            ! integrate rootfr over active layer of soil column
            rootfr_tot = 0._r8
            surface_prof_tot = 0._r8
            !$acc loop vector reduction(+:rootfr_tot,surface_prof_tot)
            do j = 1, min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
               rootfr_tot = rootfr_tot + cinput_rootfr(p,j) * dzsoi_decomp(j)
               surface_prof_tot = surface_prof_tot + surface_prof(j)  * dzsoi_decomp(j)
            end do

            if ( (altmax_lastyear_indx(c) > 0) .and. (rootfr_tot > 0._r8) .and. (surface_prof_tot > 0._r8) ) then
               ! where there is not permafrost extending to the surface, integrate the profiles over the active layer
               ! this is equivalnet to integrating over all soil layers outside of permafrost regions
               !$acc loop vector independent
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
               leaf_prof (p,1) = 1./dzsoi_decomp(1)
               stem_prof (p,1) = 1./dzsoi_decomp(1)
            endif

         end do
         call cpu_time(stopt)
         print *, "VerticleProfile::Surface reduction",(stopt-startt)*1.e+3,"ms"

         !! aggregate root profile to column
         ! call p2c (decomp, nlevdecomp_full, &
         !      cinput_rootfr(bounds%begp:bounds%endp, :), &
         !      col_cinput_rootfr(bounds%begc:bounds%endc, :), &
         !      'unity')

         call cpu_time(startt)
         !$acc parallel loop collapse(2) independent default(present) private(temp_sum,c)
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               temp_sum = 0._r8
               !$acc loop vector reduction(+:temp_sum)
               do p = col_pp%pfti(c),col_pp%pftf(c)
                  temp_sum = temp_sum + cinput_rootfr(p,j) * veg_pp%wtcol(p)
               end do
               col_cinput_rootfr(c,j) = temp_sum
            end do
         end do
         call cpu_time(stopt)
         print *, "VerticleProfile::col_rootfr",(stopt-startt)*1.e+3,"ms"

         call cpu_time(startt)
         ! repeat for column-native profiles: Ndep and Nfix
         !$acc parallel loop independent gang worker default(present) private(c,rootfr_tot,surface_prof_tot)
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            rootfr_tot = 0._r8
            surface_prof_tot = 0._r8
            ! redo column ntegration over active layer for column-native profiles
            !$acc loop vector reduction(+:rootfr_tot,surface_prof_tot)
            do j = 1, min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
               rootfr_tot = rootfr_tot + col_cinput_rootfr(fc,j) * dzsoi_decomp(j)
               surface_prof_tot = surface_prof_tot + surface_prof(j) * dzsoi_decomp(j)
            end do
            if ( (altmax_lastyear_indx(c) > 0) .and. (rootfr_tot > 0._r8) .and. (surface_prof_tot > 0._r8) ) then
               !$acc loop vector independent
               do j = 1,  min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
                  nfixation_prof(c,j) = col_cinput_rootfr(fc,j) / rootfr_tot
                  ndep_prof(c,j) = surface_prof(j)/ surface_prof_tot
                  pdep_prof(c,j) = surface_prof(j)/ surface_prof_tot
               end do
            else
               nfixation_prof(c,1) = 1./dzsoi_decomp(1)
               ndep_prof(c,1) = 1./dzsoi_decomp(1)
               pdep_prof(c,1) = 1./dzsoi_decomp(1)
            endif
         end do
         call cpu_time(stopt)
         print *, "VerticleProfile::NdepNfix",(stopt-startt)*1.E+3,"ms"

      else

         ! for one layer decomposition model, set profiles to unity
         !$acc parallel loop independent gang worker default(present)
         do j = 1, nlevdecomp
            !$acc loop vector independent private(p)
            do fp = 1, num_soilp
               p = filter_soilp(fp)
               leaf_prof(p, j) = 1._r8
               froot_prof(p,j) = 1._r8
               croot_prof(p,j) = 1._r8
               stem_prof(p,j) = 1._r8
            end do
         end do

         !$acc parallel loop independent gang worker default(present)
         do j = 1, nlevdecomp
            !$acc loop vector independent private(c)
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               nfixation_prof(c,j) = 1._r8
               ndep_prof(c,j) = 1._r8
               pdep_prof(c,j) = 1._r8
            end do
         end do

      end if

      call cpu_time(startt)
      ! check to make sure integral of all profiles = 1.
      if(debug) then
      !$acc parallel loop independent gang worker default(present) private(c,ndep_prof_sum,nfixation_prof_sum,pdep_prof_sum)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         ndep_prof_sum = 0._r8
         nfixation_prof_sum = 0._r8
         pdep_prof_sum = 0._r8
         !$acc loop vector reduction(+:ndep_prof_sum,nfixation_prof_sum,pdep_prof_sum)
         do j = 1, nlevdecomp
            ndep_prof_sum = ndep_prof_sum + ndep_prof(c,j) *  dzsoi_decomp(j)
            nfixation_prof_sum = nfixation_prof_sum + nfixation_prof(c,j) *  dzsoi_decomp(j)
            pdep_prof_sum = pdep_prof_sum + pdep_prof(c,j) *  dzsoi_decomp(j)
         end do
         print *, "ndep_prof_sum",ndep_prof_sum
         print *, "nfixation_prof_sum",nfixation_prof_sum
         print *, "pdep_prof_sum", pdep_prof_sum

      end do

      !$acc parallel loop gang worker independent default(present) &
      !$acc private(p, froot_prof_sum,croot_prof_sum,leaf_prof_sum,stem_prof_sum)
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         froot_prof_sum = 0._r8
         croot_prof_sum = 0._r8
         leaf_prof_sum = 0._r8
         stem_prof_sum = 0._r8
         !$acc loop vector reduction(+:froot_prof_sum,croot_prof_sum,leaf_prof_sum,stem_prof_sum)
         do j = 1, nlevdecomp
            froot_prof_sum = froot_prof_sum + froot_prof(p,j) *  dzsoi_decomp(j)
            croot_prof_sum = croot_prof_sum + croot_prof(p,j) *  dzsoi_decomp(j)
            leaf_prof_sum  = leaf_prof_sum + leaf_prof(p,j) *  dzsoi_decomp(j)
            stem_prof_sum  = stem_prof_sum + stem_prof(p,j) *  dzsoi_decomp(j)
         end do
         print *, "froot_prof_sum",froot_prof_sum
         print *, "croot_prof_sum",croot_prof_sum
         print *, "leaf_prof_sum ",leaf_prof_sum
         print *, "stem_prof_sum ",stem_prof_sum
      end do

      call cpu_time(stopt)
      print *, "VerticleProfile::Debug",(stopt-startt)*1.e+3,"ms"
      end if
         !$acc exit data delete(surface_prof(:),cinput_rootfr(begp:endp,1:nlevdecomp_full),&
         !$acc                  col_cinput_rootfr(1:num_soilc,1:nlevdecomp_full))
    end associate

  end subroutine decomp_vertprofiles_gpu

end module VerticalProfileMod
