module initSoilParVICMod

!-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs mapping between VIC and CLM layers 
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initSoilParVIC      ! map clm soil parameters to vic parameters
  !
  ! !REVISION HISTORY:
  ! Created by Maoyi Huang
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
 subroutine initSoilParVIC(c, claycol, sandcol, om_fraccol)
   !
   ! !DESCRIPTION:
   ! This subroutine converts default CLM soil properties to VIC parameters to be used for runoff simulations
   ! added by M. Huang
   !
   ! !USES:
   use clmtype     , only : cps 
   use clm_varcon  , only : denh2o, denice, pondmx
   use clm_varpar  , only : nlevsoi, nlayer, nlayert, nlevgrnd 

   ! !ARGUMENTS:
   implicit none
   integer , intent(in) :: c               ! column index
   real(r8),pointer     :: sandcol(:,:)    ! read in - soil texture: percent sand
   real(r8),pointer     :: claycol(:,:)    ! read in - soil texture: percent clay
   real(r8),pointer     :: om_fraccol(:,:) ! read in - organic matter: kg/m3

   ! !LOCAL VARIABLES:
   real(r8) :: om_watsat    = 0.9_r8             ! porosity of organic soil
   real(r8) :: om_hksat     = 0.1_r8             ! saturated hydraulic conductivity of organic soil [mm/s]
   real(r8) :: om_tkm       = 0.25_r8            ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
   real(r8) :: om_sucsat    = 10.3_r8            ! saturated suction for organic matter (Letts, 2000)
   real(r8) :: om_csol      = 2.5_r8             ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
   real(r8) :: om_tkd       = 0.05_r8            ! thermal conductivity of dry organic soil (Farouki, 1981)
   real(r8) :: om_b         = 2.7_r8             ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
   real(r8) :: om_expt      = 3._r8+2._r8*2.7_r8 ! soil expt for VIC        
   real(r8) :: csol_bedrock = 2.0e6_r8           ! vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
   real(r8) :: pc           = 0.5_r8             ! percolation threshold
   real(r8) :: pcbeta       = 0.139_r8           ! percolation exponent
   real(r8) :: xksat                             ! maximum hydraulic conductivity of soil [mm/s]
   real(r8) :: perc_frac                         ! "percolating" fraction of organic soil
   real(r8) :: perc_norm                         ! normalize to 1 when 100% organic soil
   real(r8) :: uncon_hksat                       ! series conductivity of mineral/organic soil
   real(r8) :: uncon_frac                        ! fraction of "unconnected" soil
   real(r8) :: temp_sum_frac                     ! sum of node fractions in each VIC layer
   real(r8) :: sandvic(1:nlayert)                ! temporary, weighted averaged sand% for VIC layers
   real(r8) :: clayvic(1:nlayert)                ! temporary, weighted averaged clay% for VIC layers
   real(r8) :: om_fracvic(1:nlayert)             ! temporary, weighted averaged organic matter fract for VIC layers
   integer  :: i, j                              ! indices
   !-------------------------------------------------------------------------------------------

   ! cps%depth(:,:)           Output: layer depth of upper layer(m) 
   ! cps%vic_clm_fract(:,:,:) Output: fraction of VIC layers in CLM layers
   ! cps%c_param(:)           Output: baseflow exponent (Qb)
   ! cps%expt(:,:)            Output: pore-size distribution related paramter(Q12)
   ! cps%ksat(:,:)            Output: Saturated hydrologic conductivity (mm/s)
   ! cps%phi_s(:,:)           Output: soil moisture dissusion parameter
   ! cps%porosity(:,:)        Output: soil porosity
   ! cps%max_moist(:,:)       Output: maximum soil moisture (ice + liq)

   !  map parameters between VIC layers and CLM layers
   cps%c_param(c) = 2._r8

   ! map the CLM layers to VIC layers 
   ! There might have better way to do this process

   do i = 1, nlayer      

      sandvic(i)    = 0._r8
      clayvic(i)    = 0._r8   
      om_fracvic(i) = 0._r8  
      temp_sum_frac = 0._r8     
      do j = 1, nlevsoi
         sandvic(i)    = sandvic(i) + sandcol(c,j) * cps%vic_clm_fract(c,i,j)
         clayvic(i)    = clayvic(i) + claycol(c,j) * cps%vic_clm_fract(c,i,j)
         om_fracvic(i) = om_fracvic(i) + om_fraccol(c,j) * cps%vic_clm_fract(c,i,j) 
         temp_sum_frac = temp_sum_frac + cps%vic_clm_fract(c,i,j)
      end do

      !average soil properties, M.Huang, 08/11/2010
      sandvic(i) = sandvic(i)/temp_sum_frac
      clayvic(i) = clayvic(i)/temp_sum_frac
      om_fracvic(i) = om_fracvic(i)/temp_sum_frac

      !make sure sand, clay and om fractions are between 0 and 100% 
      sandvic(i)    = min(100._r8 , sandvic(i))
      clayvic(i)    = min(100._r8 , clayvic(i))
      om_fracvic(i) = min(100._r8 , om_fracvic(i))
      sandvic(i)    = max(0._r8   , sandvic(i))
      clayvic(i)    = max(0._r8   , clayvic(i))
      om_fracvic(i) = max(0._r8   , om_fracvic(i))

      !calculate other parameters based on teh percentages
      cps%porosity(c, i) = 0.489_r8 - 0.00126_r8*sandvic(i)
      cps%expt(c, i) = 3._r8+ 2._r8*(2.91_r8 + 0.159_r8*clayvic(i))
      xksat = 0.0070556 *( 10.**(-0.884+0.0153*sandvic(i)) )

      !consider organic matter, M.Huang 
      cps%expt(c, i)    = (1._r8 - om_fracvic(i))*cps%expt(c, i)    + om_fracvic(i)*om_expt 
      cps%porosity(c,i) = (1._r8 - om_fracvic(i))*cps%porosity(c,i) + om_watsat*om_fracvic(i) 

      ! perc_frac is zero unless perf_frac greater than percolation threshold
      if (om_fracvic(i) > pc) then
          perc_norm=(1._r8 - pc)**(-pcbeta)
          perc_frac=perc_norm*(om_fracvic(i) - pc)**pcbeta
      else
          perc_frac=0._r8
      endif
      ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
      uncon_frac=(1._r8-om_fracvic(i))+(1._r8-perc_frac)*om_fracvic(i)

      ! uncon_hksat is series addition of mineral/organic conductivites
      if (om_fracvic(i) .lt. 1._r8) then
          uncon_hksat=uncon_frac/((1._r8-om_fracvic(i))/xksat &
                    +((1._r8-perc_frac)*om_fracvic(i))/om_hksat)
      else
          uncon_hksat = 0._r8
      end if

      cps%ksat(c,i)  = uncon_frac*uncon_hksat + (perc_frac*om_fracvic(i))*om_hksat
      cps%max_moist(c,i) = cps%porosity(c,i) * cps%depth(c,i) * 1000._r8 !in mm!

      cps%phi_s(c,i)=-(exp((1.54_r8 - 0.0095_r8*sandvic(i) + &
                 0.0063_r8*(100.0_r8-sandvic(i)-clayvic(i)))*log(10.0_r8))*9.8e-5_r8)
     end do

  end subroutine initSoilParVIC

end module initSoilParVICMod
