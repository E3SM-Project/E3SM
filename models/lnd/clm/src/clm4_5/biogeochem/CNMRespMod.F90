module CNMRespMod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNMRespMod
!
! !DESCRIPTION:
! Module holding maintenance respiration routines for coupled carbon
! nitrogen code.
!
! !USES:
   use shr_kind_mod , only: r8 => shr_kind_r8
   use clm_varpar   , only: nlevgrnd
   use shr_const_mod, only: SHR_CONST_TKFRZ
   implicit none
   save
   private
! !PUBLIC MEMBER FUNCTIONS:
   public :: CNMResp
!
! !REVISION HISTORY:
! 8/14/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: Migrated all subroutines to vector data structures.
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNMResp
!
! !INTERFACE:
subroutine CNMResp(lbc, ubc, num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
!
! !USES:
   use clmtype
   use pftvarcon    , only : npcropmin
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc                    ! column-index bounds
   integer, intent(in) :: num_soilc                 ! number of soil points in column filter
   integer, intent(in) :: filter_soilc(:)   ! column filter for soil points
   integer, intent(in) :: num_soilp                 ! number of soil points in pft filter
   integer, intent(in) :: filter_soilp(:)   ! pft filter for soil points
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 8/14/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   ! column level
   real(r8), pointer :: t_soisno(:,:) ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   ! pft level
   real(r8), pointer :: t_ref2m(:)    ! 2 m height surface air temperature (Kelvin)
   real(r8), pointer :: leafn(:)      ! (gN/m2) leaf N
   real(r8), pointer :: frootn(:)     ! (gN/m2) fine root N
   real(r8), pointer :: livestemn(:)  ! (gN/m2) live stem N
   real(r8), pointer :: livecrootn(:) ! (gN/m2) live coarse root N
   real(r8), pointer :: grainn(:)     ! (kgN/m2) grain N
   real(r8), pointer :: rootfr(:,:)   ! fraction of roots in each soil layer  (nlevgrnd)
   integer , pointer :: ivt(:)        ! pft vegetation type
   integer , pointer :: pcolumn(:)    ! index into column level quantities
   integer , pointer :: plandunit(:)  ! index into landunit level quantities
   integer , pointer :: clandunit(:)  ! index into landunit level quantities
   integer , pointer :: itypelun(:)   ! landunit type
   ! ecophysiological constants
   real(r8), pointer :: woody(:)      ! binary flag for woody lifeform (1=woody, 0=not woody)
   logical , pointer :: croplive(:)   ! Flag, true if planted, not harvested
!
! local pointers to implicit in/out arrays
!
   ! pft level
   real(r8), pointer :: leaf_mr(:)
   real(r8), pointer :: froot_mr(:)
   real(r8), pointer :: livestem_mr(:)
   real(r8), pointer :: livecroot_mr(:)
   real(r8), pointer :: grain_mr(:)
   real(r8), pointer :: lmrsun(:)         ! sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
   real(r8), pointer :: lmrsha(:)         ! shaded leaf maintenance respiration rate (umol CO2/m**2/s)
   real(r8), pointer :: laisun(:)         ! sunlit projected leaf area index
   real(r8), pointer :: laisha(:)         ! shaded projected leaf area index
   integer , pointer :: frac_veg_nosno(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p,j          ! indices
   integer :: fp             ! soil filter pft index
   integer :: fc             ! soil filter column index
   real(r8):: mr             ! maintenance respiration (gC/m2/s)
   real(r8):: br             ! base rate (gC/gN/s)
   real(r8):: q10            ! temperature dependence
   real(r8):: tc             ! temperature correction, 2m air temp (unitless)
   real(r8):: tcsoi(lbc:ubc,nlevgrnd) ! temperature correction by soil layer (unitless)
!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays
   t_soisno       => clm3%g%l%c%ces%t_soisno
   t_ref2m        => clm3%g%l%c%p%pes%t_ref2m
   leafn          => clm3%g%l%c%p%pns%leafn
   frootn         => clm3%g%l%c%p%pns%frootn
   livestemn      => clm3%g%l%c%p%pns%livestemn
   livecrootn     => clm3%g%l%c%p%pns%livecrootn
   grainn         => clm3%g%l%c%p%pns%grainn
   rootfr         => clm3%g%l%c%p%pps%rootfr
   leaf_mr        => clm3%g%l%c%p%pcf%leaf_mr
   froot_mr       => clm3%g%l%c%p%pcf%froot_mr
   livestem_mr    => clm3%g%l%c%p%pcf%livestem_mr
   livecroot_mr   => clm3%g%l%c%p%pcf%livecroot_mr
   grain_mr       => clm3%g%l%c%p%pcf%grain_mr
   lmrsun          => clm3%g%l%c%p%pcf%lmrsun
   lmrsha          => clm3%g%l%c%p%pcf%lmrsha
   laisun         => clm3%g%l%c%p%pps%laisun
   laisha         => clm3%g%l%c%p%pps%laisha
   frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
   ivt            => clm3%g%l%c%p%itype
   pcolumn        => clm3%g%l%c%p%column
   plandunit      => clm3%g%l%c%p%landunit
   clandunit      => clm3%g%l%c%landunit
   itypelun       => clm3%g%l%itype
   woody          => pftcon%woody
   croplive       => clm3%g%l%c%p%pps%croplive

   ! base rate for maintenance respiration is from:
   ! M. Ryan, 1991. Effects of climate change on plant respiration.
   ! Ecological Applications, 1(2), 157-167.
   ! Original expression is br = 0.0106 molC/(molN h)
   ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
   br = 2.525e-6_r8
   ! Peter Thornton: 3/13/09 
   ! Q10 was originally set to 2.0, an arbitrary choice, but reduced to 1.5 as part of the tuning
   ! to improve seasonal cycle of atmospheric CO2 concentration in global
   ! simulatoins
   q10 = 1.5_r8

   ! column loop to calculate temperature factors in each soil layer
   do j=1,nlevgrnd
      do fc = 1, num_soilc
         c = filter_soilc(fc)

         ! calculate temperature corrections for each soil layer, for use in
         ! estimating fine root maintenance respiration with depth

         tcsoi(c,j) = q10**((t_soisno(c,j)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
      end do
   end do

   ! pft loop for leaves and live wood
   do fp = 1, num_soilp
      p = filter_soilp(fp)

      ! calculate maintenance respiration fluxes in
      ! gC/m2/s for each of the live plant tissues.
      ! Leaf and live wood MR

      tc = q10**((t_ref2m(p)-SHR_CONST_TKFRZ - 20.0_r8)/10.0_r8)
      if (frac_veg_nosno(p) == 1) then
         leaf_mr(p) = lmrsun(p) * laisun(p) * 12.011e-6_r8 + &
                      lmrsha(p) * laisha(p) * 12.011e-6_r8
      else
         leaf_mr(p) = 0._r8
      end if

      if (woody(ivt(p)) == 1) then
         livestem_mr(p) = livestemn(p)*br*tc
         livecroot_mr(p) = livecrootn(p)*br*tc
      else if (ivt(p) >= npcropmin) then
         livestem_mr(p) = livestemn(p)*br*tc
         grain_mr(p) = grainn(p)*br*tc
      end if
   end do

   ! soil and pft loop for fine root
   do j = 1,nlevgrnd
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = pcolumn(p)

         ! Fine root MR
         ! rootfr(j) sums to 1.0 over all soil layers, and
         ! describes the fraction of root mass that is in each
         ! layer.  This is used with the layer temperature correction
         ! to estimate the total fine root maintenance respiration as a
         ! function of temperature and N content.

         froot_mr(p) = froot_mr(p) + frootn(p)*br*tcsoi(c,j)*rootfr(p,j)
      end do
   end do

end subroutine CNMResp

#endif

end module CNMRespMod
