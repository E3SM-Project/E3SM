module CNSoilLittVertTranspMod

  !-----------------------------------------------------------------------
  ! calculate vertical mixing of all decomposing C and N pools
  !
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varctl             , only : iulog, use_c13, use_c14, spinup_state, use_vertsoilc
  use clm_varcon             , only : secspday
  use decompMod              , only : bounds_type
  use abortutils             , only : endrun
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CanopyStateType        , only : canopystate_type
  use CNCarbonFluxType       , only : carbonflux_type
  use CNCarbonStateType      , only : carbonstate_type
  use CNStateType            , only : cnstate_type
  use CNNitrogenFluxType     , only : nitrogenflux_type
  use CNNitrogenStateType    , only : nitrogenstate_type
  use PhosphorusFluxType     , only : phosphorusflux_type
  use PhosphorusStateType    , only : phosphorusstate_type
  !
  implicit none
  save
  !
  public :: CNSoilLittVertTransp
  public :: readCNSoilLittVertTranspParams

  type, private :: CNSoilLittVertTranspParamsType
     real(r8) :: som_diffus                 ! Soil organic matter diffusion
     real(r8) :: cryoturb_diffusion_k       ! The cryoturbation diffusive constant
                                            ! cryoturbation to the active layer thickness
     real(r8) :: max_altdepth_cryoturbation ! (m) maximum active layer thickness for cryoturbation to occur
  end type CNSoilLittVertTranspParamsType

  type(CNSoilLittVertTranspParamsType),     private ::  CNSoilLittVertTranspParamsInst

  !
  real(r8), public :: som_adv_flux =  0._r8
  real(r8), public :: max_depth_cryoturb = 3._r8   ! (m) this is the maximum depth of cryoturbation
  real(r8) :: som_diffus                   ! [m^2/sec] = 1 cm^2 / yr
  real(r8) :: cryoturb_diffusion_k         ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr
  real(r8) :: max_altdepth_cryoturbation   ! (m) maximum active layer thickness for cryoturbation to occur
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------  
  subroutine readCNSoilLittVertTranspParams ( ncid )
    !
    use ncdio_pio   , only : file_desc_t,ncd_io
    !
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    character(len=32)  :: subname = 'CNSoilLittVertTranspType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    !
    ! read in parameters
    !
     tString='som_diffus'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     !CNSoilLittVertTranspParamsInst%som_diffus=tempr
     ! FIX(SPM,032414) - can't be pulled out since division makes things not bfb
     CNSoilLittVertTranspParamsInst%som_diffus = 1e-4_r8 / (secspday * 365._r8)  

     tString='cryoturb_diffusion_k'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     !CNSoilLittVertTranspParamsInst%cryoturb_diffusion_k=tempr
     !FIX(SPM,032414) Todo.  This constant cannot be on file since the divide makes things
     !SPM Todo.  This constant cannot be on file since the divide makes things
     !not bfb
     CNSoilLittVertTranspParamsInst%cryoturb_diffusion_k = 5e-4_r8 / (secspday * 365._r8)  ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr

     tString='max_altdepth_cryoturbation'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     CNSoilLittVertTranspParamsInst%max_altdepth_cryoturbation=tempr
    
  end subroutine readCNSoilLittVertTranspParams

  !-----------------------------------------------------------------------
  subroutine CNSoilLittVertTransp(bounds, num_soilc, filter_soilc,   &
       canopystate_vars, cnstate_vars,                               &
       carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, &
       carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,    &
       nitrogenstate_vars,nitrogenflux_vars,&
       phosphorusstate_vars,phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! Calculate vertical mixing of soil and litter pools.  Also reconcile sources and sinks of these pools 
    ! calculated in the CStateUpdate1 and NStateUpdate1 subroutines.
    ! Advection-diffusion code based on algorithm in Patankar (1980)
    ! Initial code by C. Koven and W. Riley
    !
    ! !USES:
    use clm_time_manager , only : get_step_size, get_curr_date
    use clm_varpar       , only : nlevdecomp, ndecomp_pools, nlevdecomp_full
    use clm_varcon       , only : zsoi, dzsoi_decomp, zisoi
    use TridiagonalMod   , only : Tridiagonal
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds 
    integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars

    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    real(r8) :: diffus (bounds%begc:bounds%endc,1:nlevdecomp+1)    ! diffusivity (m2/s)  (includes spinup correction, if any)
    real(r8) :: adv_flux(bounds%begc:bounds%endc,1:nlevdecomp+1)   ! advective flux (m/s)  (includes spinup correction, if any)
    real(r8) :: aaa                                                ! "A" function in Patankar
    real(r8) :: pe                                                 ! Pe for "A" function in Patankar
    real(r8) :: w_m1, w_p1                                         ! Weights for calculating harmonic mean of diffusivity
    real(r8) :: d_m1, d_p1                                         ! Harmonic mean of diffusivity
    real(r8) :: a_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)      ! "a" vector for tridiagonal matrix
    real(r8) :: b_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)      ! "b" vector for tridiagonal matrix
    real(r8) :: c_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)      ! "c" vector for tridiagonal matrix
    real(r8) :: r_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)      ! "r" vector for tridiagonal solution
    real(r8) :: d_p1_zp1(bounds%begc:bounds%endc,1:nlevdecomp+1)   ! diffusivity/delta_z for next j  (set to zero for no diffusion)
    real(r8) :: d_m1_zm1(bounds%begc:bounds%endc,1:nlevdecomp+1)   ! diffusivity/delta_z for previous j (set to zero for no diffusion)
    real(r8) :: f_p1(bounds%begc:bounds%endc,1:nlevdecomp+1)       ! water flux for next j
    real(r8) :: f_m1(bounds%begc:bounds%endc,1:nlevdecomp+1)       ! water flux for previous j
    real(r8) :: pe_p1(bounds%begc:bounds%endc,1:nlevdecomp+1)      ! Peclet # for next j
    real(r8) :: pe_m1(bounds%begc:bounds%endc,1:nlevdecomp+1)      ! Peclet # for previous j
    real(r8) :: dz_node(1:nlevdecomp+1)                            ! difference between nodes
    real(r8) :: epsilon_t (bounds%begc:bounds%endc,1:nlevdecomp+1,1:ndecomp_pools) !
    real(r8) :: conc_trcr(bounds%begc:bounds%endc,0:nlevdecomp+1)                  !
    real(r8), pointer :: conc_ptr(:,:,:)                           ! pointer, concentration state variable being transported
    real(r8), pointer :: source(:,:,:)                             ! pointer, source term
    real(r8), pointer :: trcr_tendency_ptr(:,:,:)                  ! poiner, store the vertical tendency (gain/loss due to vertical transport)
    real(r8) :: a_p_0
    real(r8) :: deficit
    integer  :: ntype
    integer  :: i_type,s,fc,c,j,l             ! indices
    integer  :: jtop(bounds%begc:bounds%endc) ! top level at each column
    real(r8) :: dtime                         ! land model time step (sec)
    integer  :: zerolev_diffus
    real(r8) :: spinup_term                   ! spinup accelerated decomposition factor, used to accelerate transport as well
    real(r8) :: epsilon                       ! small number
    integer  :: year, mon, day, sec

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
      som_diffus                 = CNSoilLittVertTranspParamsInst%som_diffus 
      cryoturb_diffusion_k       = CNSoilLittVertTranspParamsInst%cryoturb_diffusion_k 
      max_altdepth_cryoturbation = CNSoilLittVertTranspParamsInst%max_altdepth_cryoturbation 

      dtime = get_step_size()

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

         ! Set the distance between the node and the one ABOVE it   
         dz_node(1) = zsoi(1)
         do j = 2,nlevdecomp+1
            dz_node(j)= zsoi(j) - zsoi(j-1)
         enddo

      endif

      !------ loop over litter/som types
      do i_type = 1, ntype

         select case (i_type)
         case (1)  ! C
            conc_ptr          => carbonstate_vars%decomp_cpools_vr_col
            source            => carbonflux_vars%decomp_cpools_sourcesink_col
            trcr_tendency_ptr => carbonflux_vars%decomp_cpools_transport_tendency_col
         case (2)  ! N
            conc_ptr          => nitrogenstate_vars%decomp_npools_vr_col
            source            => nitrogenflux_vars%decomp_npools_sourcesink_col
            trcr_tendency_ptr => nitrogenflux_vars%decomp_npools_transport_tendency_col
         case (3)  ! P
            conc_ptr          => phosphorusstate_vars%decomp_ppools_vr_col
            source            => phosphorusflux_vars%decomp_ppools_sourcesink_col
            trcr_tendency_ptr => phosphorusflux_vars%decomp_ppools_transport_tendency_col
         case (4)
            if ( use_c13 ) then
               ! C13
               conc_ptr          => c13_carbonstate_vars%decomp_cpools_vr_col
               source            => c13_carbonflux_vars%decomp_cpools_sourcesink_col
               trcr_tendency_ptr => c13_carbonflux_vars%decomp_cpools_transport_tendency_col
            else
               ! C14
               conc_ptr          => c14_carbonstate_vars%decomp_cpools_vr_col
               source            => c14_carbonflux_vars%decomp_cpools_sourcesink_col
               trcr_tendency_ptr => c14_carbonflux_vars%decomp_cpools_transport_tendency_col
            endif
         case (5)
            if ( use_c14 .and. use_c13 ) then
               ! C14
               conc_ptr          => c14_carbonstate_vars%decomp_cpools_vr_col
               source            => c14_carbonflux_vars%decomp_cpools_sourcesink_col
               trcr_tendency_ptr => c14_carbonflux_vars%decomp_cpools_transport_tendency_col
            else
               write(iulog,*) 'error.  ncase = 5, but c13 and c14 not both enabled.'
               call endrun(msg=errMsg(__FILE__, __LINE__))
            endif
         end select

	 call get_curr_date(year, mon, day, sec)

         if (use_vertsoilc) then

            do s = 1, ndecomp_pools

               if ( spinup_state .eq. 1 ) then
                  ! increase transport (both advection and diffusion) by the same factor as accelerated decomposition for a given pool
                  spinup_term = spinup_factor(s)
               else
                  spinup_term = 1.
               endif

               if ( .not. is_cwd(s) ) then

                  do j = 1,nlevdecomp+1
                     do fc = 1, num_soilc
                        c = filter_soilc (fc)
                        !
                        if ( abs(som_adv_coef(c,j)) * spinup_term < epsilon ) then
                           adv_flux(c,j) = epsilon
                        else
			   if (spinup_term > 1 .and. year >= 40 .and. spinup_state .eq. 1) then 
                             adv_flux(c,j) = som_adv_coef(c,j) * spinup_term / cnstate_vars%scalaravg_col(c)
 			   else
                             adv_flux(c,j) = som_adv_coef(c,j) * spinup_term
			   end if			     
                        endif
                        !
                        if ( abs(som_diffus_coef(c,j)) * spinup_term < epsilon ) then
                           diffus(c,j) = epsilon
                        else
			   if (spinup_term > 1 .and. year >= 40 .and. spinup_state .eq. 1) then 
                             diffus(c,j) = som_diffus_coef(c,j) * spinup_term / cnstate_vars%scalaravg_col(c)
                           else
                             diffus(c,j) = som_diffus_coef(c,j) * spinup_term
			   end if
                        endif
                        !
                     end do
                  end do

                  ! Set Pe (Peclet #) and D/dz throughout column

                  do fc = 1, num_soilc ! dummy terms here
                     c = filter_soilc (fc)
                     conc_trcr(c,0) = 0._r8
                     conc_trcr(c,nlevdecomp+1) = 0._r8
                  end do


                  do j = 1,nlevdecomp+1
                     do fc = 1, num_soilc
                        c = filter_soilc (fc)

                        conc_trcr(c,j) = conc_ptr(c,j,s)
                        ! dz_tracer below is the difference between gridcell edges  (dzsoi_decomp)
                        ! dz_node_tracer is difference between cell centers 

                        ! Calculate the D and F terms in the Patankar algorithm
                        if (j == 1) then
                           d_m1_zm1(c,j) = 0._r8
                           w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
                           if ( diffus(c,j+1) > 0._r8 .and. diffus(c,j) > 0._r8) then
                              d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus(c,j) + w_p1 / diffus(c,j+1)) ! Harmonic mean of diffus
                           else
                              d_p1 = 0._r8
                           endif
                           d_p1_zp1(c,j) = d_p1 / dz_node(j+1)
                           f_m1(c,j) = adv_flux(c,j)  ! Include infiltration here
                           f_p1(c,j) = adv_flux(c,j+1)
                           pe_m1(c,j) = 0._r8
                           pe_p1(c,j) = f_p1(c,j) / d_p1_zp1(c,j) ! Peclet #
                        elseif (j == nlevdecomp+1) then
                           ! At the bottom, assume no gradient in d_z (i.e., they're the same)
                           w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
                           if ( diffus(c,j) > 0._r8 .and. diffus(c,j-1) > 0._r8) then
                              d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(c,j) + w_m1 / diffus(c,j-1)) ! Harmonic mean of diffus
                           else
                              d_m1 = 0._r8
                           endif
                           d_m1_zm1(c,j) = d_m1 / dz_node(j)
                           d_p1_zp1(c,j) = d_m1_zm1(c,j) ! Set to be the same
                           f_m1(c,j) = adv_flux(c,j)
                           !f_p1(c,j) = adv_flux(c,j+1)
                           f_p1(c,j) = 0._r8
                           pe_m1(c,j) = f_m1(c,j) / d_m1_zm1(c,j) ! Peclet #
                           pe_p1(c,j) = f_p1(c,j) / d_p1_zp1(c,j) ! Peclet #
                        else
                           ! Use distance from j-1 node to interface with j divided by distance between nodes
                           w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
                           if ( diffus(c,j-1) > 0._r8 .and. diffus(c,j) > 0._r8) then
                              d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(c,j) + w_m1 / diffus(c,j-1)) ! Harmonic mean of diffus
                           else
                              d_m1 = 0._r8
                           endif
                           w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
                           if ( diffus(c,j+1) > 0._r8 .and. diffus(c,j) > 0._r8) then
                              d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus(c,j) + w_p1 / diffus(c,j+1)) ! Harmonic mean of diffus
                           else
                              d_p1 = (1._r8 - w_m1) * diffus(c,j) + w_p1 * diffus(c,j+1) ! Arithmetic mean of diffus
                           endif
                           d_m1_zm1(c,j) = d_m1 / dz_node(j)
                           d_p1_zp1(c,j) = d_p1 / dz_node(j+1)
                           f_m1(c,j) = adv_flux(c,j)
                           f_p1(c,j) = adv_flux(c,j+1)
                           pe_m1(c,j) = f_m1(c,j) / d_m1_zm1(c,j) ! Peclet #
                           pe_p1(c,j) = f_p1(c,j) / d_p1_zp1(c,j) ! Peclet #
                        end if
                     enddo ! fc
                  enddo ! j; nlevdecomp


                  ! Calculate the tridiagonal coefficients
                  do j = 0,nlevdecomp +1
                     do fc = 1, num_soilc
                        c = filter_soilc (fc)
                        ! g = cgridcell(c)

                        if (j > 0 .and. j < nlevdecomp+1) then
                           a_p_0 =  dzsoi_decomp(j) / dtime
                        endif

                        if (j == 0) then ! top layer (atmosphere)
                           a_tri(c,j) = 0._r8
                           b_tri(c,j) = 1._r8
                           c_tri(c,j) = -1._r8
                           r_tri(c,j) = 0._r8
                        elseif (j == 1) then
                           a_tri(c,j) = -(d_m1_zm1(c,j) * aaa(pe_m1(c,j)) + max( f_m1(c,j), 0._r8)) ! Eqn 5.47 Patankar
                           c_tri(c,j) = -(d_p1_zp1(c,j) * aaa(pe_p1(c,j)) + max(-f_p1(c,j), 0._r8))
                           b_tri(c,j) = -a_tri(c,j) - c_tri(c,j) + a_p_0
                           r_tri(c,j) = source(c,j,s) * dzsoi_decomp(j) /dtime + (a_p_0 - adv_flux(c,j)) * conc_trcr(c,j) 
                        elseif (j < nlevdecomp+1) then
                           a_tri(c,j) = -(d_m1_zm1(c,j) * aaa(pe_m1(c,j)) + max( f_m1(c,j), 0._r8)) ! Eqn 5.47 Patankar
                           c_tri(c,j) = -(d_p1_zp1(c,j) * aaa(pe_p1(c,j)) + max(-f_p1(c,j), 0._r8))
                           b_tri(c,j) = -a_tri(c,j) - c_tri(c,j) + a_p_0
                           r_tri(c,j) = source(c,j,s) * dzsoi_decomp(j) /dtime + a_p_0 * conc_trcr(c,j)
                        else ! j==nlevdecomp+1; 0 concentration gradient at bottom
                           a_tri(c,j) = -1._r8
                           b_tri(c,j) = 1._r8
                           c_tri(c,j) = 0._r8 
                           r_tri(c,j) = 0._r8
                        endif
                     enddo ! fc; column
                  enddo ! j; nlevdecomp

                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     jtop(c) = 0
                  enddo

                  ! subtract initial concentration and source terms for tendency calculation
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     do j = 1, nlevdecomp
                        trcr_tendency_ptr(c,j,s) = 0.-(conc_trcr(c,j) + source(c,j,s))
                     end do
                  end do

                  ! Solve for the concentration profile for this time step
                  call Tridiagonal(bounds, 0, nlevdecomp+1, &
                       jtop(bounds%begc:bounds%endc), &
                       num_soilc, filter_soilc, &
                       a_tri(bounds%begc:bounds%endc, :), &
                       b_tri(bounds%begc:bounds%endc, :), &
                       c_tri(bounds%begc:bounds%endc, :), &
                       r_tri(bounds%begc:bounds%endc, :), &
                       conc_trcr(bounds%begc:bounds%endc,0:nlevdecomp+1))


                  ! add post-transport concentration to calculate tendency term
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     do j = 1, nlevdecomp
                        trcr_tendency_ptr(c,j,s) = trcr_tendency_ptr(c,j,s) + conc_trcr(c,j)
                        trcr_tendency_ptr(c,j,s) = trcr_tendency_ptr(c,j,s) / dtime
                     end do
                  end do

               else
                  ! for CWD pools, just add
                  do j = 1,nlevdecomp
                     do fc = 1, num_soilc
                        c = filter_soilc (fc)
                        conc_trcr(c,j) = conc_ptr(c,j,s) + source(c,j,s)
                     end do
                  end do

               end if ! not CWD

               do j = 1,nlevdecomp
                  do fc = 1, num_soilc
                     c = filter_soilc (fc)
                     conc_ptr(c,j,s) = conc_trcr(c,j) 
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

  end subroutine CNSoilLittVertTransp
 
end module CNSoilLittVertTranspMod
