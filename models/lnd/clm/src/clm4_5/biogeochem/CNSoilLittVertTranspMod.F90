
module CNSoilLittVertTranspMod
#ifdef CN
  
  
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNSoilLittVertTranspMod
!
! !DESCRIPTION:
! calculate vertical mixing of all decomposing C and N pools
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl  , only : iulog, use_c13, use_c14, spinup_state
  use clm_varcon  , only : secspday
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNSoilLittVertTransp   
!
! !PUBLIC DATA MEMBERS:
   real(r8), public :: som_diffus = 1e-4_r8 / (secspday * 365._r8)  ! [m^2/sec] = 1 cm^2 / yr
   real(r8), public :: som_adv_flux =  0._r8
   real(r8), public :: cryoturb_diffusion_k =  5e-4_r8 / (secspday * 365._r8)  ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr
   real(r8), public :: max_depth_cryoturb = 3._r8          ! (m) this is the maximum depth of cryoturbation
   real(r8), public :: max_altdepth_cryoturbation = 2._r8  ! (m) maximum active layer thickness for cryoturbation to occur
! !REVISION HISTORY:
!
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !PRIVATE TYPES:
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!

!
! !IROUTINE: CNSoilLittVertTransp
!
! !INTERFACE:
subroutine CNSoilLittVertTransp(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
!
! Calculate vertical mixing of soil and litter pools.  Also reconcile sources and sinks of these pools 
! calculated in the CStateUpdate1 and NStateUpdate1 subroutines.
! Advection-diffusion code based on algorithm in Patankar (1980)
! Initial code by C. Koven and W. Riley
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varpar   , only: nlevdecomp, ndecomp_pools, nlevdecomp_full
   use clm_varcon   , only: zsoi, dzsoi_decomp, zisoi
   use TridiagonalMod  , only : Tridiagonal
   use abortutils   , only : endrun

! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: lbc, ubc        ! column-index bounds

! !LOCAL VARIABLES:
!EOP
   ! real(r8) :: som_diffus_coef (lbc:ubc,1:nlevdecomp+1)     ! diffusivity (m2/s)  
   ! real(r8) :: som_adv_coef(lbc:ubc,1:nlevdecomp+1)         ! advective flux (m/s)
   real(r8), pointer :: som_adv_coef(:,:)                     ! SOM advective flux (m/s)
   real(r8), pointer :: som_diffus_coef(:,:)                  ! SOM diffusivity due to bio/cryo-turbation (m2/s)
   real(r8) :: diffus (lbc:ubc,1:nlevdecomp+1)                ! diffusivity (m2/s)  (includes spinup correction, if any)
   real(r8) :: adv_flux(lbc:ubc,1:nlevdecomp+1)               ! advective flux (m/s)  (includes spinup correction, if any)
   real(r8) :: aaa                                            ! "A" function in Patankar
   real(r8) :: pe                                             ! Pe for "A" function in Patankar
   real(r8) :: w_m1, w_p1                                     ! Weights for calculating harmonic mean of diffusivity
   real(r8) :: d_m1, d_p1                                     ! Harmonic mean of diffusivity
   real(r8) :: a_tri (lbc:ubc,0:nlevdecomp+1)                 ! "a" vector for tridiagonal matrix
   real(r8) :: b_tri (lbc:ubc,0:nlevdecomp+1)                 ! "b" vector for tridiagonal matrix
   real(r8) :: c_tri (lbc:ubc,0:nlevdecomp+1)                 ! "c" vector for tridiagonal matrix
   real(r8) :: r_tri (lbc:ubc,0:nlevdecomp+1)                 ! "r" vector for tridiagonal solution
   real(r8) :: d_p1_zp1 (lbc:ubc,1:nlevdecomp+1)              ! diffusivity/delta_z for next j  
                                                              ! (set to zero for no diffusion)
   real(r8) :: d_m1_zm1 (lbc:ubc,1:nlevdecomp+1)              ! diffusivity/delta_z for previous j 
                                                              ! (set to zero for no diffusion)
   real(r8) :: f_p1     (lbc:ubc,1:nlevdecomp+1)              ! water flux for next j
   real(r8) :: f_m1     (lbc:ubc,1:nlevdecomp+1)              ! water flux for previous j
   real(r8) :: pe_p1     (lbc:ubc,1:nlevdecomp+1)             ! Peclet # for next j
   real(r8) :: pe_m1     (lbc:ubc,1:nlevdecomp+1)             ! Peclet # for previous j
   real(r8) :: dz_node(1:nlevdecomp+1)                        ! difference between nodes
   real(r8) :: epsilon_t (lbc:ubc,1:nlevdecomp+1,1:ndecomp_pools) !
   real(r8) :: conc_trcr(lbc:ubc,0:nlevdecomp+1)  
   real(r8), pointer :: conc_ptr(:,:,:)                       ! pointer for concentration state variable being transported
   real(r8), pointer :: source(:,:,:)                         ! pointer for source term
   logical, pointer  :: is_cwd(:)                             ! TRUE => pool is a cwd pool
   real(r8) :: a_p_0
   real(r8) :: deficit
   real(r8), pointer :: trcr_tendency_ptr(:,:,:)              ! pointer to store the vertical tendency 
                                                              ! (gain/loss due to vertical transport)
   real(r8), pointer :: altmax(:)                             ! maximum annual depth of thaw
   real(r8), pointer :: altmax_lastyear(:)                    ! prior year maximum annual depth of thaw


   integer :: ntype
   integer :: i_type,s,fc,c,j,l          ! indices
   integer :: jtop(lbc:ubc)              ! top level at each column
   real(r8) :: dtime                     ! land model time step (sec)
   integer :: zerolev_diffus
   real(r8), pointer :: spinup_factor(:) ! spinup accelerated decomposition factor, used to accelerate transport as well
   real(r8) :: spinup_term               ! spinup accelerated decomposition factor, used to accelerate transport as well
   real(r8) :: epsilon                   ! small number

!EOP
!-----------------------------------------------------------------------
   



   ! Set statement functions
   aaa (pe) = max (0._r8, (1._r8 - 0.1_r8 * abs(pe))**5)  ! A function from Patankar, Table 5.2, pg 95
   
   is_cwd                                  => decomp_cascade_con%is_cwd
   spinup_factor                           => decomp_cascade_con%spinup_factor
   altmax                          => clm3%g%l%c%cps%altmax
   altmax_lastyear                 => clm3%g%l%c%cps%altmax_lastyear
   som_adv_coef                    => clm3%g%l%c%cps%som_adv_coef
   som_diffus_coef                 => clm3%g%l%c%cps%som_diffus_coef

   
   dtime = get_step_size()
   
   ntype = 2
   if ( use_c13 ) then
      ntype = ntype+1
   endif
   if ( use_c14 ) then
      ntype = ntype+1
   endif
   spinup_term = 1._r8
   epsilon = 1.e-30
   
#ifdef VERTSOILC
   !------ first get diffusivity / advection terms -------!
   ! use different mixing rates for bioturbation and cryoturbation, with fixed bioturbation and cryoturbation set to a maximum depth
   do fc = 1, num_soilc
      c = filter_soilc (fc)
      if  ( ( max(altmax(c), altmax_lastyear(c)) .le. max_altdepth_cryoturbation ) .and. ( max(altmax(c), altmax_lastyear(c)) .gt. 0._r8) ) then
         ! use mixing profile modified slightly from Koven et al. (2009): constant through active layer, linear decrease from base of active layer to zero at a fixed depth
         do j = 1,nlevdecomp+1
            if ( zisoi(j) .lt. max(altmax(c), altmax_lastyear(c)) ) then
               som_diffus_coef(c,j) = cryoturb_diffusion_k 
               som_adv_coef(c,j) = 0._r8
            else
               som_diffus_coef(c,j) = max(cryoturb_diffusion_k * ( 1._r8 - ( zisoi(j) - max(altmax(c), altmax_lastyear(c)) ) / &
                    ( max_depth_cryoturb - max(altmax(c), altmax_lastyear(c)) ) ), 0._r8)  ! go linearly to zero between ALT and max_depth_cryoturb
               som_adv_coef(c,j) = 0._r8
            endif
         end do
      elseif (  max(altmax(c), altmax_lastyear(c)) .gt. 0._r8 ) then
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

#endif
   
   !------ loop over litter/som types
   do i_type = 1, ntype
            
      select case (i_type)
      case (1)  ! C
         conc_ptr => clm3%g%l%c%ccs%decomp_cpools_vr
         source    => clm3%g%l%c%ccf%decomp_cpools_sourcesink
         trcr_tendency_ptr => clm3%g%l%c%ccf%decomp_cpools_transport_tendency
      case (2)  ! N
         conc_ptr => clm3%g%l%c%cns%decomp_npools_vr
         source    => clm3%g%l%c%cnf%decomp_npools_sourcesink
         trcr_tendency_ptr => clm3%g%l%c%cnf%decomp_npools_transport_tendency
      case (3)
         if ( use_c13 ) then
            ! C13
            conc_ptr => clm3%g%l%c%cc13s%decomp_cpools_vr
            source    => clm3%g%l%c%cc13f%decomp_cpools_sourcesink
            trcr_tendency_ptr => clm3%g%l%c%cc13f%decomp_cpools_transport_tendency
         else
            ! C14
            conc_ptr => clm3%g%l%c%cc14s%decomp_cpools_vr
            source    => clm3%g%l%c%cc14f%decomp_cpools_sourcesink
            trcr_tendency_ptr => clm3%g%l%c%cc14f%decomp_cpools_transport_tendency
         endif
      case (4)
         if ( use_c14 .and. use_c13 ) then
            ! C14
            conc_ptr => clm3%g%l%c%cc14s%decomp_cpools_vr
            source    => clm3%g%l%c%cc14f%decomp_cpools_sourcesink
            trcr_tendency_ptr => clm3%g%l%c%cc14f%decomp_cpools_transport_tendency
         else
            write(iulog,*) 'error.  ncase = 4, but c13 and c14 not both enabled.'
            call endrun()
         endif
      end select

#ifdef VERTSOILC
      
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
                  if ( abs(som_adv_coef(c,j)) * spinup_term .lt. epsilon ) then
                     adv_flux(c,j) = epsilon
                  else
                     adv_flux(c,j) = som_adv_coef(c,j) * spinup_term
                  endif
                  !
                  if ( abs(som_diffus_coef(c,j)) * spinup_term .lt. epsilon ) then
                     diffus(c,j) = epsilon
                  else
                     diffus(c,j) = som_diffus_coef(c,j) * spinup_term
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
                     if ( diffus(c,j+1) .gt. 0._r8 .and. diffus(c,j) .gt. 0._r8) then
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
                     if ( diffus(c,j) .gt. 0._r8 .and. diffus(c,j-1) .gt. 0._r8) then
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
                     if ( diffus(c,j-1) .gt. 0._r8 .and. diffus(c,j) .gt. 0._r8) then
                        d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(c,j) + w_m1 / diffus(c,j-1)) ! Harmonic mean of diffus
                     else
                        d_m1 = 0._r8
                     endif
                     w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
                     if ( diffus(c,j+1) .gt. 0._r8 .and. diffus(c,j) .gt. 0._r8) then
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
                  elseif (j .lt. nlevdecomp+1) then
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
            call Tridiagonal(lbc, ubc, 0, nlevdecomp+1, jtop, num_soilc, filter_soilc, &
                 a_tri, b_tri, c_tri, r_tri, conc_trcr(lbc:ubc,0:nlevdecomp+1))
            
            
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
      
#else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      
#endif
      
   end do  ! i_type
   
end subroutine CNSoilLittVertTransp
 
#endif
end module CNSoilLittVertTranspMod

