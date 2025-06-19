module SoilLittVertTranspMod

  !-----------------------------------------------------------------------
  ! calculate vertical mixing of all decomposing C and N pools
  !
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use elm_varctl             , only : iulog, use_c13, use_c14, spinup_state, use_vertsoilc
  use elm_varcon             , only : secspday
  use decompMod              , only : bounds_type
  use abortutils             , only : endrun
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CanopyStateType        , only : canopystate_type
  use CNStateType            , only : cnstate_type
  use elm_varctl             , only : nu_com
  use ColumnDataType         , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType         , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType         , only : col_ns, col_nf, col_ps, col_pf
  use timeinfoMod
  !
  implicit none
  save
  !
  public :: SoilLittVertTransp
  public :: createLitterTransportList
  public :: readSoilLittVertTranspParams
  private :: calc_diffus_advflux

  type, public :: SoilLittVertTranspParamsType
     real(r8)  :: som_diffus                  ! Soil organic matter diffusion
     real(r8)  :: cryoturb_diffusion_k        ! The cryoturbation diffusive constant cryoturbation to the active layer thickness
     real(r8)  :: max_altdepth_cryoturbation  ! (m) maximum active layer thickness for cryoturbation to occur
  end type SoilLittVertTranspParamsType

  type(SoilLittVertTranspParamsType), public ::  SoilLittVertTranspParamsInst
  !$acc declare create(SoilLittVertTranspParamsInst)

  type, public :: ConcTransportType
     !! Type that points to decomposition pools for vertical transport calculations

     real(r8), pointer :: conc_ptr(:,:,:) => null()
     real(r8), pointer :: src_ptr(:,:,:)  => null()
     real(r8), pointer :: trcr_tend_ptr(:,:,:) => null()
  end type ConcTransportType
  type(ConcTransportType), public, allocatable :: transport_ptr_list(:)
  !$acc declare create(transport_ptr_list(:))

  real(r8), public :: som_adv_flux =  0._r8
  !$acc declare create(som_adv_flux)
  real(r8), public :: max_depth_cryoturb = 3._r8   ! (m) this is the maximum depth of cryoturbation
  !$acc declare create(max_depth_cryoturb)
  !-----------------------------------------------------------------------

contains

   subroutine createLitterTransportList()
      ! This subroutine creates a list that will point to the
      ! litter/som fields needed for the vertical transport
      ! calculations.

      implicit none

      integer :: ntype

      ntype = 3
      if ( use_c13 ) then
         ntype = ntype+1
      endif
      if ( use_c14 ) then
         ntype = ntype+1
      endif

      allocate(transport_ptr_list(ntype))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! C
      transport_ptr_list(1)%conc_ptr        => col_cs%decomp_cpools_vr
      transport_ptr_list(1)%src_ptr         => col_cf%decomp_cpools_sourcesink
      transport_ptr_list(1)%trcr_tend_ptr   => col_cf%decomp_cpools_transport_tendency
      ! N
      transport_ptr_list(2)%conc_ptr        => col_ns%decomp_npools_vr
      transport_ptr_list(2)%src_ptr         => col_nf%decomp_npools_sourcesink
      transport_ptr_list(2)%trcr_tend_ptr   => col_nf%decomp_npools_transport_tendency
      ! P
      transport_ptr_list(3)%conc_ptr        => col_ps%decomp_ppools_vr
      transport_ptr_list(3)%src_ptr         => col_pf%decomp_ppools_sourcesink
      transport_ptr_list(3)%trcr_tend_ptr   => col_pf%decomp_ppools_transport_tendency
      ! c13 and c14 if there
      if(use_c14 .and. use_c13) then
         !
         transport_ptr_list(4)%conc_ptr       => c13_col_cs%decomp_cpools_vr
         transport_ptr_list(4)%src_ptr        => c13_col_cf%decomp_cpools_sourcesink
         transport_ptr_list(4)%trcr_tend_ptr  => c13_col_cf%decomp_cpools_transport_tendency
         !
         transport_ptr_list(5)%conc_ptr       => c14_col_cs%decomp_cpools_vr
         transport_ptr_list(5)%src_ptr        => c14_col_cf%decomp_cpools_sourcesink
         transport_ptr_list(5)%trcr_tend_ptr  => c14_col_cf%decomp_cpools_transport_tendency
      else
         if(use_c13) then
            transport_ptr_list(4)%conc_ptr       => c13_col_cs%decomp_cpools_vr
            transport_ptr_list(4)%src_ptr        => c13_col_cf%decomp_cpools_sourcesink
            transport_ptr_list(4)%trcr_tend_ptr  => c13_col_cf%decomp_cpools_transport_tendency
         end if
         if (use_c14) then
            transport_ptr_list(4)%conc_ptr       => c14_col_cs%decomp_cpools_vr
            transport_ptr_list(4)%src_ptr        => c14_col_cf%decomp_cpools_sourcesink
            transport_ptr_list(4)%trcr_tend_ptr  => c14_col_cf%decomp_cpools_transport_tendency
         end if
      end if

   end subroutine createLitterTransportList

   subroutine cleanupLitterTransportList()
      !! Nullifies pointers and frees allocated memory
      integer :: i

      if (allocated(transport_ptr_list)) then
         ! Loop through each element and nullify pointers
         do i = 1, size(transport_ptr_list)
            if (associated(transport_ptr_list(i)%conc_ptr)) then
               nullify(transport_ptr_list(i)%conc_ptr)
            endif
            if (associated(transport_ptr_list(i)%src_ptr)) then
               nullify(transport_ptr_list(i)%src_ptr)
            endif
            if (associated(transport_ptr_list(i)%trcr_tend_ptr)) then
               nullify(transport_ptr_list(i)%trcr_tend_ptr)
            endif
         end do

         ! Deallocate the transport_ptr_list array itself
         deallocate(transport_ptr_list)
      endif

   end subroutine cleanupLitterTransportList


  !-----------------------------------------------------------------------
  subroutine readSoilLittVertTranspParams ( ncid )
    !
    use ncdio_pio   , only : file_desc_t,ncd_io
    !
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    character(len=32)  :: subname = 'SoilLittVertTranspType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    !
    ! read in parameters
    !
    ! REMOVE THESE?
    tString='som_diffus'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    !SoilLittVertTranspParamsInst%som_diffus=tempr
    ! FIX(SPM,032414) - can't be pulled out since division makes things not bfb
    SoilLittVertTranspParamsInst%som_diffus = 1e-4_r8 / (secspday * 365._r8)

    tString='cryoturb_diffusion_k'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    !SoilLittVertTranspParamsInst%cryoturb_diffusion_k=tempr
    !FIX(SPM,032414) Todo.  This constant cannot be on file since the divide makes things
    !SPM Todo.  This constant cannot be on file since the divide makes things
    !not bfb

    SoilLittVertTranspParamsInst%cryoturb_diffusion_k = 5e-4_r8 / (secspday * 365._r8)  ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr

    tString='max_altdepth_cryoturbation'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    SoilLittVertTranspParamsInst%max_altdepth_cryoturbation=tempr

    !$acc enter data copyin(SoilLittVertTranspParamsInst)
  end subroutine readSoilLittVertTranspParams

  function aaa(pe) result(res)
     !$acc routine seq
     implicit none
     real(r8) :: res
     real(r8) :: pe
     res =  max (0._r8, (1._r8 - 0.1_r8 * abs(pe))**5)
  end function

  !-----------------------------------------------------------------------
  subroutine SoilLittVertTransp(num_soilc, filter_soilc,   &
       canopystate_vars, cnstate_vars )
    !
    ! !DESCRIPTION:
    ! Calculate vertical mixing of soil and litter pools.  Also reconcile sources and sinks of these pools
    ! calculated in the CarbonStateUpdate1 and NStateUpdate1 subroutines.
    ! Advection-diffusion code based on algorithm in Patankar (1980)
    ! Initial code by C. Koven and W. Riley
    !
    ! !USES:
    use elm_varpar       , only : nlevdecomp, ndecomp_pools, nlevdecomp_full
    use elm_varcon       , only : zsoi, dzsoi_decomp, zisoi
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8) :: pe          ! Pe for "A" function in Patankar
    real(r8) :: w_m1, w_p1  ! Weights for calculating harmonic mean of diffusivity
    real(r8) :: d_m1, d_p1  ! Harmonic mean of diffusivity
    real(r8) :: d_p1_zp1    ! diffusivity/delta_z for next j  (set to zero for no diffusion)
    real(r8) :: d_m1_zm1    ! diffusivity/delta_z for previous j (set to zero for no diffusion)
    real(r8) :: pe_p1       ! Peclet # for next j
    real(r8) :: pe_m1       ! Peclet # for previous j
    real(r8) :: dz_node,dz_nodep1          ! difference between nodes
    real(r8) :: a_p_0
    integer  :: ntype
    integer  :: i_type,s,fc,c,j,l  ! indices
    integer  :: jtop(num_soilc)    ! top level at each column
    real(r8) :: spinup_term                  ! spinup accelerated decomposition factor, used to accelerate transport as well
    real(r8), parameter :: epsilon=1.e-30     ! small number
    !!added to remove arrays:
    real(r8) :: diffus_j, diffus_jm1, diffus_jp1 ! diffusivity (m2/s)  (includes spinup correction, if any)
    real(r8) :: adv_flux_j,adv_flux_jm1, adv_flux_jp1        ! advective flux (m/s)  (includes spinup correction, if any)
    real(r8) :: a_tri(num_soilc,0:nlevdecomp+1,ndecomp_pools)      ! "a" vector for tridiagonal matrix
    real(r8) :: b_tri(num_soilc,0:nlevdecomp+1,ndecomp_pools)      ! "b" vector for tridiagonal matrix
    real(r8) :: c_tri(num_soilc,0:nlevdecomp+1,ndecomp_pools)      ! "c" vector for tridiagonal matrix
    real(r8) :: r_tri(num_soilc,0:nlevdecomp+1,ndecomp_pools)      ! "r" vector for tridiagonal solution
    real(r8) :: conc_trcr(num_soilc,0:nlevdecomp+1,ndecomp_pools)                  !
    real(r8) :: bet
    real(r8) :: gam(0:nlevdecomp+1)
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------

    ! Set statement functions
    associate(                                                      &
         is_cwd           => decomp_cascade_con%is_cwd            , & ! Input:  [logical (:)    ]  TRUE => pool is a cwd pool
         spinup_factor    => decomp_cascade_con%spinup_factor     , & ! Input:  [real(r8) (:)   ]  spinup accelerated decomposition factor, used to accelerate transport as well

         altmax           => canopystate_vars%altmax_col          , & ! Input:  [real(r8) (:)   ]  maximum annual depth of thaw
         altmax_lastyear  => canopystate_vars%altmax_lastyear_col , & ! Input:  [real(r8) (:)   ]  prior year maximum annual depth of thaw

         som_adv_coef     => cnstate_vars%som_adv_coef_col       , & ! Output: [real(r8) (:,:) ]  SOM advective flux (m/s)
         som_diffus_coef  => cnstate_vars%som_diffus_coef_col    ,  & ! Output: [real(r8) (:,:) ]  SOM diffusivity due to bio/cryo-turbation (m2/s)
         ! !Set parameters of vertical mixing of SOM
          som_diffus                 => SoilLittVertTranspParamsInst%som_diffus   , &
          cryoturb_diffusion_k       => SoilLittVertTranspParamsInst%cryoturb_diffusion_k  , &
          max_altdepth_cryoturbation => SoilLittVertTranspParamsInst%max_altdepth_cryoturbation &
         )
      
      !$acc enter data create(a_tri(:,:,:),b_tri(:,:,:),&
      !$acc     c_tri(:,:,:),r_tri(:,:,:), &
      !$acc     conc_trcr(:,:,:), gam(:) )
      ntype = 3
      if ( use_c13 ) then
         ntype = ntype+1
      endif
      if ( use_c14 ) then
         ntype = ntype+1
      endif
   
      !$acc enter data create(spinup_term, i_type) 
      spinup_term = 1._r8
      !$acc update device(spinup_term)

      if (use_vertsoilc) then
         !------ first get diffusivity / advection terms -------!
         ! use different mixing rates for bioturbation and cryoturbation, with fixed bioturbation and cryoturbation set to a maximum depth
         !$acc parallel loop independent gang default(present)
         do j = 1,nlevdecomp+1
            !$acc loop vector independent private(c)
            do fc = 1, num_soilc
               c = filter_soilc (fc)
               if  ( ( max(altmax(c), altmax_lastyear(c)) <= max_altdepth_cryoturbation ) .and. &
                  ( max(altmax(c), altmax_lastyear(c)) > 0._r8) ) then
                  ! use mixing profile modified slightly from Koven et al. (2009): constant through active layer, linear decrease from base of active layer to zero at a fixed depth
                  if ( zisoi(j) < max(altmax(c), altmax_lastyear(c)) ) then
                     som_diffus_coef(c,j) = cryoturb_diffusion_k
                     som_adv_coef(c,j) = 0._r8
                  else
                     som_diffus_coef(c,j) = max(cryoturb_diffusion_k * &
                          ( 1._r8 - ( zisoi(j) - max(altmax(c), altmax_lastyear(c)) ) / &
                          ( max_depth_cryoturb - max(altmax(c), altmax_lastyear(c)) ) ), 0._r8)  ! go linearly to zero between ALT and max_depth_cryoturb
                     som_adv_coef(c,j) = 0._r8
                  endif
               elseif (  max(altmax(c), altmax_lastyear(c)) > 0._r8 ) then
                  ! constant advection, constant diffusion
                  som_adv_coef(c,j) = som_adv_flux
                  som_diffus_coef(c,j) = som_diffus
               else
                  ! completely frozen soils--no mixing
                  som_adv_coef(c,j) = 0._r8
                  som_diffus_coef(c,j) = 0._r8
               endif
            end do
         end do
      endif
   
      !------ loop over litter/som types
      do i_type = 1, ntype

         !$acc update device(i_type)

         if (use_vertsoilc) then
            ! Set Pe (Peclet #) and D/dz throughout column

            !$acc parallel loop independent gang default(present)
            do s = 1, ndecomp_pools
               if ( .not. is_cwd(s) ) then
                  !$acc loop independent worker vector private(c)
                  do fc = 1, num_soilc ! dummy terms here
                     c = filter_soilc (fc)
                     conc_trcr(fc,0,s) = 0._r8
                     conc_trcr(fc,nlevdecomp+1,s) = 0._r8

                     a_tri(fc,0,s) = 0._r8
                     b_tri(fc,0,s) = 1._r8
                     c_tri(fc,0,s) = -1._r8
                     r_tri(fc,0,s) = 0._r8

                     conc_trcr(fc,nlevdecomp+1,s) = transport_ptr_list(i_type)%conc_ptr(c,nlevdecomp+1,s)
                     a_tri(fc,nlevdecomp+1,s) = -1._r8
                     b_tri(fc,nlevdecomp+1,s) = 1._r8
                     c_tri(fc,nlevdecomp+1,s) = 0._r8
                     r_tri(fc,nlevdecomp+1,s) = 0._r8
                  end do
               end if
            end do

            !$acc parallel loop independent gang worker vector collapse(3) default(present) 
            do s = 1, ndecomp_pools
               do j = 1,nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     if(.not. is_cwd(s)) then

                        if ( spinup_state .eq. 1 ) then
                           ! increase transport (both advection and diffusion) by the same factor as accelerated decomposition for a given pool
                           spinup_term = spinup_factor(s)
                        else
                           spinup_term = 1.
                        endif
                        conc_trcr(fc,j,s) = transport_ptr_list(i_type)%conc_ptr(c,j,s)
                        ! dz_tracer below is the difference between gridcell edges  (dzsoi_decomp)
                        ! dz_node_tracer is difference between cell centers
                        call calc_diffus_advflux(spinup_term,year_curr, som_diffus_coef(c,j), som_adv_coef(c,j), &
                                                 cnstate_vars%scalaravg_col(c,j),adv_flux_j, diffus_j)

                        ! Calculate the D and F terms in the Patankar algorithm
                        if (j == 1) then
                          call calc_diffus_advflux(spinup_term,year_curr, som_diffus_coef(c,j+1), som_adv_coef(c,j+1), &
                                                   cnstate_vars%scalaravg_col(c,j+1),adv_flux_jp1, diffus_jp1)
                           dz_nodep1 =  zsoi(j+1) - zsoi(j)
                           d_m1_zm1 = 0._r8
                           w_p1 = (zsoi(j+1) - zisoi(j)) / dz_nodep1
                           if (diffus_jp1 > 0._r8 .and. diffus_j > 0._r8) then
                             d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus_j + w_p1 / diffus_jp1) ! Harmonic mean of diffus
                           else
                              d_p1 = 0._r8
                           endif

                           d_p1_zp1 = d_p1 / dz_nodep1
                           pe_m1 = 0._r8
                           pe_p1 = adv_flux_jp1 / d_p1_zp1 ! Peclet #

                           a_p_0 =  dzsoi_decomp(j) / dtime_mod
                           a_tri(fc,j,s) = -(d_m1_zm1 * aaa(pe_m1) + max( adv_flux_j, 0._r8)) ! Eqn 5.47 Patankar
                           c_tri(fc,j,s) = -(d_p1_zp1 * aaa(pe_p1) + max(-adv_flux_jp1, 0._r8))
                           b_tri(fc,j,s) = -a_tri(fc,j,s) - c_tri(fc,j,s) + a_p_0
                           r_tri(fc,j,s) = transport_ptr_list(i_type)%src_ptr(c,j,s) * dzsoi_decomp(j) /dtime_mod + (a_p_0 - adv_flux_j) * conc_trcr(fc,j,s)
                        else
                          ! Use distance from j-1 node to interface with j divided by distance between nodes
                          call calc_diffus_advflux(spinup_term,year_curr, som_diffus_coef(c,j-1), som_adv_coef(c,j-1), &
                                                   cnstate_vars%scalaravg_col(c,j-1),adv_flux_jm1, diffus_jm1)

                          call calc_diffus_advflux(spinup_term,year_curr, som_diffus_coef(c,j+1), som_adv_coef(c,j+1), &
                                                   cnstate_vars%scalaravg_col(c,j+1),adv_flux_jp1, diffus_jp1)
                           ! Use distance from j-1 node to interface with j divided by distance between nodes
                           dz_node = zsoi(j) - zsoi(j-1)
                           w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node

                           if ( diffus_jm1 > 0._r8 .and. diffus_j > 0._r8) then
                              d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus_j + w_m1 / diffus_jm1) ! Harmonic mean of diffus
                           else
                              d_m1 = 0._r8
                           endif

                           dz_nodep1 = zsoi(j+1) - zsoi(j)
                           w_p1 = (zsoi(j+1) - zisoi(j)) / dz_nodep1

                           if ( diffus_jp1 > 0._r8 .and. diffus_j > 0._r8) then
                              d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus_j + w_p1 / diffus_jp1) ! Harmonic mean of diffus
                           else
                              d_p1 = (1._r8 - w_m1) * diffus_j + w_p1 * diffus_jp1  ! Arithmetic mean of diffus
                           endif
                           d_m1_zm1 = d_m1 / dz_node
                           d_p1_zp1 = d_p1 / dz_nodep1
                           pe_m1 = adv_flux_j / d_m1_zm1 ! Peclet #
                           pe_p1 = adv_flux_jp1 / d_p1_zp1 ! Peclet #

                           a_p_0 =  dzsoi_decomp(j) / dtime_mod
                           a_tri(fc,j,s) = -(d_m1_zm1 * aaa(pe_m1) + max( adv_flux_j, 0._r8)) ! Eqn 5.47 Patankar
                           c_tri(fc,j,s) = -(d_p1_zp1 * aaa(pe_p1) + max(-adv_flux_jp1, 0._r8))
                           b_tri(fc,j,s) = -a_tri(fc,j,s) - c_tri(fc,j,s) + a_p_0
                           r_tri(fc,j,s) = transport_ptr_list(i_type)%src_ptr(c,j,s) * dzsoi_decomp(j) /dtime_mod + a_p_0 * conc_trcr(fc,j,s)
                        end if
                     end if
                  enddo ! fc
               enddo ! j; nlevdecomp
            end do ! s: ndecomp_pools

            ! subtract initial concentration and source terms for tendency calculation
            !$acc parallel loop independent collapse(3) gang vector default(present)
            do s = 1, ndecomp_pools
               do j = 1, nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     if(.not. is_cwd(s)) then
                        transport_ptr_list(i_type)%trcr_tend_ptr(c,j,s) = 0._r8 - (conc_trcr(fc,j,s) + transport_ptr_list(i_type)%src_ptr(c,j,s))
                     end if
                  end do
               end do
            end do

            ! Solve for the concentration profile for this time step

            !$acc parallel loop independent gang worker vector collapse(2) default(present) private(bet, gam(0:nlevdecomp+1))
            do s = 1, ndecomp_pools
               do fc = 1,num_soilc
                  if(.not. is_cwd(s)) then
                     bet = b_tri(fc,0,s)

                     !$acc loop seq
                     do j = 0, nlevdecomp+1
                        if (j == 0) then
                           conc_trcr(fc,j,s) = r_tri(fc,j,s) / bet
                        else
                           gam(j) = c_tri(fc,j-1,s) / bet
                           bet = b_tri(fc,j,s) - a_tri(fc,j,s) * gam(j)
                           conc_trcr(fc,j,s) = (r_tri(fc,j,s) - a_tri(fc,j,s)*conc_trcr(fc,j-1,s)) / bet
                        end if
                     end do

                     !$acc loop seq
                     do j = nlevdecomp,1,-1
                       conc_trcr(fc,j,s) = conc_trcr(fc,j,s) - gam(j+1) * conc_trcr(fc,j+1,s)
                     end do
                  end if
               end do
            end do

            ! add post-transport concentration to calculate tendency term
            !$acc parallel loop independent gang collapse(2) default(present)
            do s = 1, ndecomp_pools
               do j = 1, nlevdecomp
                  if(.not. is_cwd(s)) then
                     !$acc loop vector independent private(c)
                     do fc = 1, num_soilc
                        c = filter_soilc (fc)
                        transport_ptr_list(i_type)%trcr_tend_ptr(c,j,s) = (transport_ptr_list(i_type)%trcr_tend_ptr(c,j,s) + conc_trcr(fc,j,s))/dtime_mod
                     end do
                  end if
                  !
               end do
            end do

            ! for CWD pools, just add
            !$acc parallel loop independent gang default(present)
            do s = 1, ndecomp_pools
               if(is_cwd(s)) then
                  !$acc loop worker vector collapse(2) independent private(c)
                  do j = 1,nlevdecomp
                     do fc = 1, num_soilc
                        c = filter_soilc (fc)
                        conc_trcr(fc,j,s) = transport_ptr_list(i_type)%conc_ptr(c,j,s) + transport_ptr_list(i_type)%src_ptr(c,j,s)
                     end do
                  end do
               end if
            end do


            !$acc parallel loop independent gang collapse(2) default(present)
            do s = 1, ndecomp_pools
               do j = 1,nlevdecomp
                  !$acc loop vector independent private(c)
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     transport_ptr_list(i_type)%conc_ptr(c,j,s) = conc_trcr(fc,j,s)
                  end do
               end do

            end do ! s (pool loop)

         else !use_vertsoilc?

            !! for single level case, no transport; just update the fluxes calculated in the StateUpdate1 subroutines
            !$acc parallel loop independent collapse(2) default(present)
            do l = 1, ndecomp_pools
               do j = 1,nlevdecomp
                  !$acc loop vector independent private(c)
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     transport_ptr_list(i_type)%conc_ptr(c,j,l) = transport_ptr_list(i_type)%conc_ptr(c,j,l) &
                                                                  +transport_ptr_list(i_type)%src_ptr(c,j,l)
                     transport_ptr_list(i_type)%trcr_tend_ptr(c,j,l) = 0._r8
                  end do
               end do
            end do

         endif

      end do  ! i_type
   
      !$acc exit data delete(a_tri(:,:,:),b_tri(:,:,:),&
      !$acc     c_tri(:,:,:),r_tri(:,:,:), gam(:), &
      !$acc     conc_trcr(:,:,:), spinup_term, i_type)
    end associate

  end subroutine SoilLittVertTransp

  subroutine calc_diffus_advflux(spinup_term, year, som_diffus_coef, som_adv_coef,&
                                   cnscalaravg_col, adv_fluxj, diffusj)
    !$acc routine seq
    real(r8), intent(in)  ::  spinup_term
    integer ,  intent(in) :: year
    real(r8), intent(in)  ::  som_diffus_coef, som_adv_coef, cnscalaravg_col
    real(r8), intent(out) ::  adv_fluxj, diffusj

    real(r8), parameter :: eps = 1d-30
    logical :: use_scaling

    ! Determine if we should scale by cnscalaravg_col
    use_scaling = (spinup_term > 1 .and. year >= 40 .and. spinup_state .eq. 1)

    ! Compute adv_fluxj
    if (abs(som_adv_coef) * spinup_term < eps) then
        adv_fluxj = eps
    else
        adv_fluxj = som_adv_coef * spinup_term
        if (use_scaling) adv_fluxj = adv_fluxj / cnscalaravg_col
    endif

    ! Compute diffusj
    if (abs(som_diffus_coef) * spinup_term < eps) then
        diffusj = eps
    else
        diffusj = som_diffus_coef * spinup_term
        if (use_scaling) diffusj = diffusj / cnscalaravg_col
    endif

 end subroutine calc_diffus_advflux


end module SoilLittVertTranspMod
