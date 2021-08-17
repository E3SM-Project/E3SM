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
  public :: readSoilLittVertTranspParams

  type, private :: SoilLittVertTranspParamsType
     real(r8),pointer  :: som_diffus                => null() ! Soil organic matter diffusion
     real(r8),pointer  :: cryoturb_diffusion_k      => null() ! The cryoturbation diffusive constant cryoturbation to the active layer thickness
     real(r8),pointer  :: max_altdepth_cryoturbation => null() ! (m) maximum active layer thickness for cryoturbation to occur
  end type SoilLittVertTranspParamsType

  type(SoilLittVertTranspParamsType), public ::  SoilLittVertTranspParamsInst
  !$acc declare create(SoilLittVertTranspParamsInst)

  !
  real(r8), public :: som_adv_flux =  0._r8
  !$acc declare copyin(som_adv_flux)
  real(r8), public :: max_depth_cryoturb = 3._r8   ! (m) this is the maximum depth of cryoturbation
  !$acc declare copyin(max_depth_cryoturb)
  real(r8) :: som_diffus                   ! [m^2/sec] = 1 cm^2 / yr
  real(r8) :: cryoturb_diffusion_k         ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr
  real(r8) :: max_altdepth_cryoturbation   ! (m) maximum active layer thickness for cryoturbation to occur
  !$acc declare create(som_diffus                )
  !$acc declare create(cryoturb_diffusion_k      )
  !$acc declare create(max_altdepth_cryoturbation)
  !-----------------------------------------------------------------------

contains

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
    allocate(SoilLittVertTranspParamsInst%som_diffus                )
    allocate(SoilLittVertTranspParamsInst%cryoturb_diffusion_k      )
    allocate(SoilLittVertTranspParamsInst%max_altdepth_cryoturbation)
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

  !-----------------------------------------------------------------------
  subroutine SoilLittVertTransp( num_soilc, filter_soilc,   &
       canopystate_vars, cnstate_vars )
    !
    ! !DESCRIPTION:
    ! Calculate vertical mixing of soil and litter pools.  Also reconcile sources and sinks of these pools
    ! calculated in the CarbonStateUpdate1 and NStateUpdate1 subroutines.
    ! Advection-diffusion code based on algorithm in Patankar (1980)
    ! Initial code by C. Koven and W. Riley
    !
    ! !USES:
      !$acc routine seq
    use elm_varpar       , only : nlevdecomp, ndecomp_pools, nlevdecomp_full
    use elm_varcon       , only : zsoi, dzsoi_decomp, zisoi
    use TridiagonalMod   , only : Tridiagonal_SoilLittVertTransp
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    real(r8) :: dtime  ! land model time step (sec)
    integer  :: year, mon, day, secs

    !
    ! !LOCAL VARIABLES:
    real(r8) :: aaa                                                ! "A" function in Patankar
    real(r8) :: pe                                                 ! Pe for "A" function in Patankar
    real(r8) :: w_m1, w_p1                                         ! Weights for calculating harmonic mean of diffusivity
    real(r8) :: d_m1, d_p1                                         ! Harmonic mean of diffusivity
    real(r8) :: a_tri(num_soilc,0:nlevdecomp+1)      ! "a" vector for tridiagonal matrix
    real(r8) :: b_tri(num_soilc,0:nlevdecomp+1)      ! "b" vector for tridiagonal matrix
    real(r8) :: c_tri(num_soilc,0:nlevdecomp+1)      ! "c" vector for tridiagonal matrix
    real(r8) :: r_tri(num_soilc,0:nlevdecomp+1)      ! "r" vector for tridiagonal solution
    real(r8) :: d_p1_zp1   ! diffusivity/delta_z for next j  (set to zero for no diffusion)
    real(r8) :: d_m1_zm1   ! diffusivity/delta_z for previous j (set to zero for no diffusion)
    real(r8) :: pe_p1       ! Peclet # for next j
    real(r8) :: pe_m1       ! Peclet # for previous j
    real(r8) :: dz_node,dz_nodep1          ! difference between nodes
    real(r8) :: conc_trcr(num_soilc,0:nlevdecomp+1)                  !
    real(r8), pointer :: conc_ptr(:,:,:)                           ! pointer, concentration state variable being transported
    real(r8), pointer :: source(:,:,:)                             ! pointer, source term
    real(r8), pointer :: trcr_tendency_ptr(:,:,:)                  ! poiner, store the vertical tendency (gain/loss due to vertical transport)
    real(r8) :: a_p_0
    real(r8) :: deficit
    integer  :: ntype
    integer  :: i_type,s,fc,c,j,l             ! indices
    integer  :: jtop(num_soilc) ! top level at each column
    integer  :: zerolev_diffus
    real(r8) :: spinup_term                   ! spinup accelerated decomposition factor, used to accelerate transport as well
    real(r8) :: epsilon                       ! small number
    !!added to remove arrays:
    real(r8) :: diffus_j, diffus_jm1, diffus_jp1 ! diffusivity (m2/s)  (includes spinup correction, if any)
    real(r8) :: adv_flux_j,adv_flux_jm1, adv_flux_jp1        ! advective flux (m/s)  (includes spinup correction, if any)
    !-----------------------------------------------------------------------


    ! Set statement functions
    aaa (pe) = max (0._r8, (1._r8 - 0.1_r8 * abs(pe))**5)  ! A function from Patankar, Table 5.2, pg 95

    associate(                                                      &
         is_cwd           => decomp_cascade_con%is_cwd            , & ! Input:  [logical (:)    ]  TRUE => pool is a cwd pool
         spinup_factor    => decomp_cascade_con%spinup_factor     , & ! Input:  [real(r8) (:)   ]  spinup accelerated decomposition factor, used to accelerate transport as well

         altmax           => canopystate_vars%altmax_col          , & ! Input:  [real(r8) (:)   ]  maximum annual depth of thaw
         altmax_lastyear  => canopystate_vars%altmax_lastyear_col , & ! Input:  [real(r8) (:)   ]  prior year maximum annual depth of thaw

         som_adv_coef     => cnstate_vars%som_adv_coef_col        , & ! Output: [real(r8) (:,:) ]  SOM advective flux (m/s)
         som_diffus_coef  => cnstate_vars%som_diffus_coef_col       & ! Output: [real(r8) (:,:) ]  SOM diffusivity due to bio/cryo-turbation (m2/s)
         )

      !Set parameters of vertical mixing of SOM
      som_diffus                 = SoilLittVertTranspParamsInst%som_diffus
      cryoturb_diffusion_k       = SoilLittVertTranspParamsInst%cryoturb_diffusion_k
      max_altdepth_cryoturbation = SoilLittVertTranspParamsInst%max_altdepth_cryoturbation

      dtime = dtime_mod
      year = year_curr
      mon = mon_curr
      day = day_curr


      ntype = 3
      if ( use_c13 ) then
         ntype = ntype+1
      endif
      if ( use_c14 ) then
         ntype = ntype+1
      endif

      spinup_term = 1._r8
      epsilon = 1.e-30

      if (use_vertsoilc) then
         !------ first get diffusivity / advection terms -------!
         ! use different mixing rates for bioturbation and cryoturbation, with fixed bioturbation and cryoturbation set to a maximum depth
         do fc = 1, num_soilc
            c = filter_soilc (fc)
            if  ( ( max(altmax(c), altmax_lastyear(c)) <= max_altdepth_cryoturbation ) .and. &
                 ( max(altmax(c), altmax_lastyear(c)) > 0._r8) ) then
               ! use mixing profile modified slightly from Koven et al. (2009): constant through active layer, linear decrease from base of active layer to zero at a fixed depth
               do j = 1,nlevdecomp+1
                  if ( zisoi(j) < max(altmax(c), altmax_lastyear(c)) ) then
                     som_diffus_coef(c,j) = cryoturb_diffusion_k
                     som_adv_coef(c,j) = 0._r8
                  else
                     som_diffus_coef(c,j) = max(cryoturb_diffusion_k * &
                          ( 1._r8 - ( zisoi(j) - max(altmax(c), altmax_lastyear(c)) ) / &
                          ( max_depth_cryoturb - max(altmax(c), altmax_lastyear(c)) ) ), 0._r8)  ! go linearly to zero between ALT and max_depth_cryoturb
                     som_adv_coef(c,j) = 0._r8
                  endif
               end do
            elseif (  max(altmax(c), altmax_lastyear(c)) > 0._r8 ) then
               ! constant advection, constant diffusion
               do j = 1,nlevdecomp+1
                  som_adv_coef(c,j) = som_adv_flux
                  som_diffus_coef(c,j) = som_diffus
               end do
            else
               ! completely frozen soils--no mixing
               do j = 1,nlevdecomp+1
                  som_adv_coef(c,j) = 0._r8
                  som_diffus_coef(c,j) = 0._r8
               end do
            endif
         end do

         ! ! Set the distance between the node and the one ABOVE it
         ! dz_node(1) = zsoi(1)
         ! do j = 2,nlevdecomp+1
         !    dz_node(j)= zsoi(j) - zsoi(j-1)
         ! enddo

      endif

      !------ loop over litter/som types
      do i_type = 1, ntype

         select case (i_type)
         case (1)  ! C
            conc_ptr          => col_cs%decomp_cpools_vr
            source            => col_cf%decomp_cpools_sourcesink
            trcr_tendency_ptr => col_cf%decomp_cpools_transport_tendency
         case (2)  ! N
            conc_ptr          => col_ns%decomp_npools_vr
            source            => col_nf%decomp_npools_sourcesink
            trcr_tendency_ptr => col_nf%decomp_npools_transport_tendency
         case (3)  ! P
            conc_ptr          => col_ps%decomp_ppools_vr
            source            => col_pf%decomp_ppools_sourcesink
            trcr_tendency_ptr => col_pf%decomp_ppools_transport_tendency
         case (4)
            if ( use_c13 ) then
               ! C13
               conc_ptr          => c13_col_cs%decomp_cpools_vr
               source            => c13_col_cf%decomp_cpools_sourcesink
               trcr_tendency_ptr => c13_col_cf%decomp_cpools_transport_tendency
            else
               ! C14
               conc_ptr          => c14_col_cs%decomp_cpools_vr
               source            => c14_col_cf%decomp_cpools_sourcesink
               trcr_tendency_ptr => c14_col_cf%decomp_cpools_transport_tendency
            endif
         case (5)
            if ( use_c14 .and. use_c13 ) then
               ! C14
               conc_ptr          => c14_col_cs%decomp_cpools_vr
               source            => c14_col_cf%decomp_cpools_sourcesink
               trcr_tendency_ptr => c14_col_cf%decomp_cpools_transport_tendency
            else
#ifndef _OPENACC
               write(iulog,*) 'error.  ncase = 5, but c13 and c14 not both enabled.'
               call endrun(msg=errMsg(__FILE__, __LINE__))
#endif
            endif
         end select

         if (use_vertsoilc) then

            do s = 1, ndecomp_pools

               if ( spinup_state .eq. 1 ) then
                  ! increase transport (both advection and diffusion) by the same factor as accelerated decomposition for a given pool
                  spinup_term = spinup_factor(s)
               else
                  spinup_term = 1.
               endif

               if ( .not. is_cwd(s) ) then

                  ! Set Pe (Peclet #) and D/dz throughout column
                  !
                  do fc = 1, num_soilc ! dummy terms here
                     c = filter_soilc (fc)
                     conc_trcr(fc,0) = 0._r8
                     conc_trcr(fc,nlevdecomp+1) = 0._r8

                     a_tri(fc,0) = 0._r8
                     b_tri(fc,0) = 1._r8
                     c_tri(fc,0) = -1._r8
                     r_tri(fc,0) = 0._r8

                     conc_trcr(fc,nlevdecomp+1) = conc_ptr(c,nlevdecomp+1,s)
                     a_tri(fc,nlevdecomp+1) = -1._r8
                     b_tri(fc,nlevdecomp+1) = 1._r8
                     c_tri(fc,nlevdecomp+1) = 0._r8
                     r_tri(fc,nlevdecomp+1) = 0._r8

                  end do


                  do j = 1,nlevdecomp
                     do fc = 1, num_soilc
                        c = filter_soilc (fc)

                        conc_trcr(fc,j) = conc_ptr(c,j,s)
                        ! dz_tracer below is the difference between gridcell edges  (dzsoi_decomp)
                        ! dz_node_tracer is difference between cell centers
                        call calc_diffus_advflux(spinup_term,year, som_diffus_coef(c,j), som_adv_coef(c,j), &
                                                 cnstate_vars%scalaravg_col(c,j),adv_flux_j, diffus_j)

                        ! Use distance from j-1 node to interface with j divided by distance between nodes
                        call calc_diffus_advflux(spinup_term,year, som_diffus_coef(c,j-1), som_adv_coef(c,j-1), &
                                                 cnstate_vars%scalaravg_col(c,j-1),adv_flux_jm1, diffus_jm1)

                        call calc_diffus_advflux(spinup_term,year, som_diffus_coef(c,j+1), som_adv_coef(c,j+1), &
                                                 cnstate_vars%scalaravg_col(c,j+1),adv_flux_jp1, diffus_jp1)

                        ! Calculate the D and F terms in the Patankar algorithm
                        if (j == 1) then
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

                           a_p_0 =  dzsoi_decomp(j) / dtime
                           a_tri(fc,j) = -(d_m1_zm1 * aaa(pe_m1) + max( adv_flux_j, 0._r8)) ! Eqn 5.47 Patankar
                           c_tri(fc,j) = -(d_p1_zp1 * aaa(pe_p1) + max(-adv_flux_jp1, 0._r8))
                           b_tri(fc,j) = -a_tri(fc,j) - c_tri(fc,j) + a_p_0
                           r_tri(fc,j) = source(c,j,s) * dzsoi_decomp(j) /dtime + (a_p_0 - adv_flux_j) * conc_trcr(fc,j)
                        else
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
                             d_p1 = (1._r8 - w_m1) * diffus_j + w_p1 * diffus_jp1 ! Arithmetic mean of diffus
                           endif
                           d_m1_zm1 = d_m1 / dz_node
                           d_p1_zp1 = d_p1 / dz_nodep1
                           pe_m1 = adv_flux_j / d_m1_zm1 ! Peclet #
                           pe_p1 = adv_flux_jp1 / d_p1_zp1 ! Peclet #

                           a_p_0 =  dzsoi_decomp(j) / dtime
                           a_tri(fc,j) = -(d_m1_zm1 * aaa(pe_m1) + max( adv_flux_j, 0._r8)) ! Eqn 5.47 Patankar
                           c_tri(fc,j) = -(d_p1_zp1 * aaa(pe_p1) + max(-adv_flux_jp1, 0._r8))
                           b_tri(fc,j) = -a_tri(fc,j) - c_tri(fc,j) + a_p_0
                           r_tri(fc,j) = source(c,j,s) * dzsoi_decomp(j) /dtime + a_p_0 * conc_trcr(fc,j)
                        end if
                     enddo ! fc
                  enddo ! j; nlevdecomp


                  do fc = 1, num_soilc
                     jtop(fc) = 0
                  enddo

                  ! subtract initial concentration and source terms for tendency calculation
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     do j = 1, nlevdecomp
                        trcr_tendency_ptr(c,j,s) = 0.-(conc_trcr(fc,j) + source(c,j,s))
                     end do
                  end do

                  ! Solve for the concentration profile for this time step
                  call Tridiagonal_SoilLittVertTransp(0, nlevdecomp+1, &
                       jtop(1:num_soilc), &
                       num_soilc, filter_soilc, &
                       a_tri(1:num_soilc, :), &
                       b_tri(1:num_soilc, :), &
                       c_tri(1:num_soilc, :), &
                       r_tri(1:num_soilc, :), &
                       conc_trcr(1:num_soilc,0:nlevdecomp+1))

                  ! add post-transport concentration to calculate tendency term
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     do j = 1, nlevdecomp
                        trcr_tendency_ptr(c,j,s) = trcr_tendency_ptr(c,j,s) + conc_trcr(fc,j)
                        trcr_tendency_ptr(c,j,s) = trcr_tendency_ptr(c,j,s) / dtime
                     end do
                  end do

               else
                  ! for CWD pools, just add
                  do j = 1,nlevdecomp
                     do fc = 1, num_soilc
                        c = filter_soilc (fc)
                        conc_trcr(fc,j) = conc_ptr(c,j,s) + source(c,j,s)
                     end do
                  end do

               end if ! not CWD

               do j = 1,nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     conc_ptr(c,j,s) = conc_trcr(fc,j)
                  end do
               end do

            end do ! s (pool loop)

         else

            !! for single level case, no transport; just update the fluxes calculated in the StateUpdate1 subroutines
            do l = 1, ndecomp_pools
               do j = 1,nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)

                     conc_ptr(c,j,l) = conc_ptr(c,j,l) + source(c,j,l)

                     trcr_tendency_ptr(c,j,l) = 0._r8

                  end do
               end do
            end do

         endif

      end do  ! i_type

    end associate

  end subroutine SoilLittVertTransp

  subroutine calc_diffus_advflux(spinup_term,year, som_diffus_coef, som_adv_coef, cnscalaravg_col, &
                          adv_fluxj, diffusj)
    !$acc routine seq
    real(r8), intent(in)  ::  spinup_term
    integer , value, intent(in) :: year
    real(r8), intent(in)  ::  som_diffus_coef, som_adv_coef, cnscalaravg_col
    real(r8), intent(out) ::  adv_fluxj, diffusj

    real(r8) , parameter :: eps = 1d-30

    if ( abs(som_adv_coef) * spinup_term < eps ) then
       adv_fluxj = eps
    else
         if (spinup_term > 1 .and. year >= 40 .and. spinup_state .eq. 1) then
           adv_fluxj = som_adv_coef * spinup_term / cnscalaravg_col
         else
           adv_fluxj = som_adv_coef * spinup_term
         end if
    endif
    !
    if ( abs(som_diffus_coef) * spinup_term < eps ) then
       diffusj = eps
    else
         if (spinup_term > 1 .and. year >= 40 .and. spinup_state .eq. 1) then
           diffusj = som_diffus_coef * spinup_term / cnscalaravg_col
         else
           diffusj = som_diffus_coef * spinup_term
         end if
    endif
  end subroutine calc_diffus_advflux

end module SoilLittVertTranspMod
