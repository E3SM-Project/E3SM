module CNC14DecayMod

  !-----------------------------------------------------------------------
  ! Module for 14-carbon flux variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar   , only: ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: C14Decay
  public:: C14BombSpike
  public:: C14_init_BombSpike

  ! !PUBLIC TYPES:
  logical, public :: use_c14_bombspike = .false.         ! do we use time-varying atmospheric C14?
  character(len=256), public :: atm_c14_filename = ' '   ! file name of C14 input data

  ! !PRIVATE TYPES:
  real(r8), allocatable, private :: atm_c14file_time(:)
  real(r8), allocatable, private :: atm_delta_c14(:)
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the radioactive decay of C14
    !
    ! !USES:
    use clmtype
    use clm_time_manager, only: get_step_size, get_days_per_year
    use clm_varcon, only: secspday
    use clm_varctl, only : spinup_state
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_soilc       ! number of soil columns filter
    integer, intent(in) :: filter_soilc(:) ! filter for soil columns
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
    !
    ! !LOCAL VARIABLES:
    integer :: fp,j,l,p,fc,c,i
    real(r8) :: dt                             ! radiation time step (seconds)
    real(r8) :: half_life
    real(r8) :: decay_const
    real(r8) :: days_per_year                  ! days per year
    real(r8) :: spinup_term               ! spinup accelerated decomposition factor, used to accelerate transport as well
    !-----------------------------------------------------------------------

   associate(& 
   decomp_cpools_vr               =>    cc14s%decomp_cpools_vr       , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   seedc                          =>    cc14s%seedc                  , & ! InOut:  [real(r8) (:)]                                          
   cpool                          =>    pc14s%cpool                  , & ! InOut:  [real(r8) (:)]  (gC/m2) temporary photosynthate C pool  
   xsmrpool                       =>    pc14s%xsmrpool               , & ! InOut:  [real(r8) (:)]  (gC/m2) execss maint resp C pool        
   deadcrootc                     =>    pc14s%deadcrootc             , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C              
   deadcrootc_storage             =>    pc14s%deadcrootc_storage     , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C storage      
   deadcrootc_xfer                =>    pc14s%deadcrootc_xfer        , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer     
   deadstemc                      =>    pc14s%deadstemc              , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C                     
   deadstemc_storage              =>    pc14s%deadstemc_storage      , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C storage             
   deadstemc_xfer                 =>    pc14s%deadstemc_xfer         , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C transfer            
   frootc                         =>    pc14s%frootc                 , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C                     
   frootc_storage                 =>    pc14s%frootc_storage         , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C storage             
   frootc_xfer                    =>    pc14s%frootc_xfer            , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C transfer            
   gresp_storage                  =>    pc14s%gresp_storage          , & ! InOut:  [real(r8) (:)]  (gC/m2) growth respiration storage      
   gresp_xfer                     =>    pc14s%gresp_xfer             , & ! InOut:  [real(r8) (:)]  (gC/m2) growth respiration transfer     
   leafc                          =>    pc14s%leafc                  , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C                          
   leafc_storage                  =>    pc14s%leafc_storage          , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C storage                  
   leafc_xfer                     =>    pc14s%leafc_xfer             , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C transfer                 
   livecrootc                     =>    pc14s%livecrootc             , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C              
   livecrootc_storage             =>    pc14s%livecrootc_storage     , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C storage      
   livecrootc_xfer                =>    pc14s%livecrootc_xfer        , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C transfer     
   livestemc                      =>    pc14s%livestemc              , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C                     
   livestemc_storage              =>    pc14s%livestemc_storage      , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C storage             
   livestemc_xfer                 =>    pc14s%livestemc_xfer         , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C transfer            
   pft_ctrunc                     =>    pc14s%pft_ctrunc             , & ! InOut:  [real(r8) (:)]  (gC/m2) pft-level sink for C truncation 
   spinup_factor                  =>    decomp_cascade_con%spinup_factor    & ! InOut:  [real(r8) (:)]  factor for AD spinup associated with each pool
   )

    ! set time steps
    dt = real( get_step_size(), r8 )
    days_per_year = get_days_per_year()

    half_life = 5568._r8 * secspday * days_per_year  !! libby half-life value, for comparison against ages calculated with this value
    ! half_life = 5730._r8 * secspday * days_per_year  !! recent half-life value
    decay_const = - log(0.5_r8) / half_life


    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       seedc(c) = seedc(c) *  (1._r8 - decay_const * dt)
    end do ! end of columns loop

    do l = 1, ndecomp_pools
       if ( spinup_state .eq. 1) then
          ! speed up radioactive decay by the same factor as decomposition so tat SOM ages prematurely in all respects
          spinup_term = spinup_factor(l) 
       else
          spinup_term = 1.
       endif
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             decomp_cpools_vr(c,j,l) = decomp_cpools_vr(c,j,l) * (1._r8 - decay_const * spinup_term * dt)
          end do
       end do
    end do ! end of columns loop

 
    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       cpool(p)              = cpool(p)               * (1._r8 - decay_const * dt)
       xsmrpool(p)           = xsmrpool(p)            * (1._r8 - decay_const * dt)
       deadcrootc(p)         = deadcrootc(p)          * (1._r8 - decay_const * dt)
       deadcrootc_storage(p) = deadcrootc_storage(p)  * (1._r8 - decay_const * dt)
       deadcrootc_xfer(p)    = deadcrootc_xfer(p)     * (1._r8 - decay_const * dt)
       deadstemc(p)          = deadstemc(p)           * (1._r8 - decay_const * dt)
       deadstemc_storage(p)  = deadstemc_storage(p)   * (1._r8 - decay_const * dt)
       deadstemc_xfer(p)     = deadstemc_xfer(p)      * (1._r8 - decay_const * dt)
       frootc(p)             = frootc(p)              * (1._r8 - decay_const * dt)
       frootc_storage(p)     = frootc_storage(p)      * (1._r8 - decay_const * dt)
       frootc_xfer(p)        = frootc_xfer(p)         * (1._r8 - decay_const * dt)
       gresp_storage(p)      = gresp_storage(p)       * (1._r8 - decay_const * dt)
       gresp_xfer(p)         = gresp_xfer(p)          * (1._r8 - decay_const * dt)
       leafc(p)              = leafc(p)               * (1._r8 - decay_const * dt)
       leafc_storage(p)      = leafc_storage(p)       * (1._r8 - decay_const * dt)
       leafc_xfer(p)         = leafc_xfer(p)          * (1._r8 - decay_const * dt)
       livecrootc(p)         = livecrootc(p)          * (1._r8 - decay_const * dt)
       livecrootc_storage(p) = livecrootc_storage(p)  * (1._r8 - decay_const * dt)
       livecrootc_xfer(p)    = livecrootc_xfer(p)     * (1._r8 - decay_const * dt)
       livestemc(p)          = livestemc(p)           * (1._r8 - decay_const * dt)
       livestemc_storage(p)  = livestemc_storage(p)   * (1._r8 - decay_const * dt)
       livestemc_xfer(p)     = livestemc_xfer(p)      * (1._r8 - decay_const * dt)
       pft_ctrunc(p)         = pft_ctrunc(p)          * (1._r8 - decay_const * dt)
    end do

    end associate 
  end subroutine C14Decay

  !-----------------------------------------------------------------------
  subroutine C14BombSpike(num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! for transient pulse simulation, impose a simplified bomb spike
    !
    ! !USES:
    use clmtype
    use clm_time_manager, only: get_curr_date,get_days_per_year
    use clm_varcon  , only : c14ratio,secspday
    use ncdio_pio
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
    !
    ! !LOCAL VARIABLES:
    integer :: yr, mon, day, tod, offset
    real(r8) :: dateyear
    real(r8) :: delc14o2_atm 
    real(r8) :: days_per_year                    ! days per year
    integer :: fp, p, nt
    integer :: ind_below
    integer :: ntim_atm_ts
    real(r8) :: twt_1, twt_2                     ! weighting fractions for interpolating
    real(r8) :: min, max
    !-----------------------------------------------------------------------

   associate(& 
   rc14_atm  =>  pepv%rc14_atm                  & ! InOut:  [real(r8) (:)] C14O2/C12O2 in atmosphere                
   )

   if ( use_c14_bombspike ) then

      ! get current date
      call get_curr_date(yr, mon, day, tod, offset)
      days_per_year = get_days_per_year()
      dateyear = real(yr) + real(mon)/12._r8 + real(day)/days_per_year + real(tod)/(secspday*days_per_year)

      ! find points in atm timeseries to interpolate between
      ntim_atm_ts = size(atm_c14file_time)
      ind_below = 0
      do nt = 1, ntim_atm_ts
         if (dateyear .ge. atm_c14file_time(nt) ) then
            ind_below = ind_below+1
         endif
      end do

      ! interpolate between nearest two points in atm c14 timeseries
      if (ind_below .eq. 0 ) then 
         delc14o2_atm = atm_delta_c14(1)
      elseif (ind_below .eq. ntim_atm_ts ) then
         delc14o2_atm = atm_delta_c14(ntim_atm_ts)
      else
         twt_2 = min(1._r8, max(0._r8,(dateyear-atm_c14file_time(ind_below)) &
             / (atm_c14file_time(ind_below+1)-atm_c14file_time(ind_below))))
         twt_1 = 1._r8 - twt_2
         delc14o2_atm = atm_delta_c14(ind_below) * twt_1 +  atm_delta_c14(ind_below+1) * twt_2
      endif

      ! change delta units to ratio, put on pft loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)   
         rc14_atm(p) = (delc14o2_atm * 1.e-3_r8 + 1._r8) * c14ratio
      end do
      
   else
      ! for constant 14c concentration
      ! pft loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)   
         rc14_atm(p) = c14ratio
      end do
   endif

    end associate 
 end subroutine C14BombSpike

 !-----------------------------------------------------------------------
 subroutine C14_init_BombSpike()
   !
   ! !DESCRIPTION:
   ! read netcdf file containing a timeseries of atmospheric delta C14 values; save in module-level array 
   !
   ! !USES:
   use ncdio_pio
   use fileutils   , only : getfil
   use abortutils  , only : endrun
   use clm_varctl  , only : iulog
   use spmdMod     , only : masterproc
   use shr_log_mod , only : errMsg => shr_log_errMsg
   !
   ! !ARGUMENTS:
   implicit none
   !
   ! !LOCAL VARIABLES:
   character(len=256) :: locfn           ! local file name
   type(file_desc_t)  :: ncid            ! netcdf id
   integer :: dimid,varid                ! input netCDF id's
   integer :: ntim                       ! number of input data time samples
   integer :: t
   !-----------------------------------------------------------------------
   
   if ( use_c14_bombspike ) then
      
      if ( masterproc ) then
         write(iulog, *) 'C14_init_BombSpike: preparing to open file:'
         write(iulog, *) trim(locfn)
      endif

      call getfil(atm_c14_filename, locfn, 0)
      
      call ncd_pio_openfile (ncid, trim(locfn), 0)
      
      call ncd_inqdlen(ncid,dimid,ntim,'time')

      !! allocate arrays based on size of netcdf timeseries
      allocate(atm_c14file_time(ntim))
      allocate(atm_delta_c14(ntim))

      call ncd_io(ncid=ncid, varname='time', flag='read', data=atm_c14file_time)
      
      call ncd_io(ncid=ncid, varname='atm_delta_c14', flag='read', data=atm_delta_c14)
      
      call ncd_pio_closefile(ncid)

      ! check to make sure that time dimension is well behaved
      do t = 2, ntim
         if ( atm_c14file_time(t) - atm_c14file_time(t-1) .le. 0._r8 ) then
            write(iulog, *) 'C14_init_BombSpike: error.  time axis must be monotonically increasing'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         endif
      end do
   endif
   
end subroutine C14_init_BombSpike

end module CNC14DecayMod
