module accFldsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains subroutines that initialize, update and extract
  ! the user-specified fields over user-defined intervals. Each interval
  ! and accumulation type is unique to each field processed.
  ! Subroutine [initAccumFlds] defines the fields to be processed
  ! and the type of accumulation. Subroutine [updateAccumFlds] does
  ! the actual accumulation for a given field. Fields are accumulated
  ! by calls to subroutine [update_accum_field]. To accumulate a field,
  ! it must first be defined in subroutine [initAccumFlds] and then
  ! accumulated by calls to [updateAccumFlds].
  ! Four types of accumulations are possible:
  !   o average over time interval
  !   o running mean over time interval
  !   o running accumulation over time interval
  ! Time average fields are only valid at the end of the averaging interval.
  ! Running means are valid once the length of the simulation exceeds the
  ! averaging interval. Accumulated fields are continuously accumulated.
  ! The trigger value "-99999." resets the accumulation to zero.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_log_mod , only: errMsg => shr_log_errMsg
  use abortutils,   only: endrun
  use clm_varctl,   only: iulog, use_cndv
  use clm_varpar,   only: crop_prog
  use decompMod ,   only: bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initAccFlds     ! Initialization accumulator fields
  public :: initAccClmtype  ! Initialize clmtype variables obtained from accum fields
  public :: updateAccFlds   ! Update accumulator fields
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine initAccFlds(bounds)
    !
    ! !DESCRIPTION:
    ! Initializes accumulator and sets up array of accumulated fields
    !
    ! !USES:
    use accumulMod       , only : init_accum_field, print_accum_fields 
    use clm_time_manager , only : get_step_size
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: dtime                     !time step size
    integer, parameter :: not_used = huge(1)
    !------------------------------------------------------------------------

    ! Hourly average of 2m temperature.
    dtime = get_step_size()
    call init_accum_field(name='TREFAV', units='K', &
         desc='average over an hour of 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! Hourly average of Urban 2m temperature.
    call init_accum_field(name='TREFAV_U', units='K', &
         desc='average over an hour of urban 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! Hourly average of Rural 2m temperature.
    call init_accum_field(name='TREFAV_R', units='K', &
         desc='average over an hour of rural 2-m temperature', &
         accum_type='timeavg', accum_period=nint(3600._r8/dtime), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 24hr average of vegetation temperature 
    call init_accum_field (name='T_VEG24', units='K', &
         desc='24hr average of vegetation temperature', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 240hr average of vegetation temperature
    call init_accum_field (name='T_VEG240', units='K', &
         desc='240hr average of vegetation temperature', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 24hr average of direct solar radiation 
    call init_accum_field (name='FSD24', units='W/m2', &
         desc='24hr average of direct solar radiation', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 240hr average of direct solar radiation
    call init_accum_field (name='FSD240', units='W/m2', &
         desc='240hr average of direct solar radiation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 24hr average of diffuse solar radiation
    call init_accum_field (name='FSI24', units='W/m2', &
         desc='24hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 240hr average of diffuse solar radiation
    call init_accum_field (name='FSI240', units='W/m2', &
         desc='240hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 24hr average of fraction of canopy that is sunlit 
    call init_accum_field (name='FSUN24', units='fraction', &
         desc='24hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! 240hr average of fraction of canopy that is sunlit 
    call init_accum_field (name='FSUN240', units='fraction', &
         desc='240hr average of diffuse solar radiation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! Average of LAI from previous and current timestep      
    call init_accum_field (name='LAIP', units='m2/m2', &
         desc='leaf area index average over timestep', &
         accum_type='runmean', accum_period=1, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    ! The following is a running mean.
    ! The accumulation period is set to -10 for a 10-day running mean.

    call init_accum_field (name='T10', units='K', &
         desc='10-day running mean of 2-m temperature', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1,init_value=SHR_CONST_TKFRZ+20._r8)

    if (use_cndv) then
       ! 30-day average of 2m temperature.
       
       call init_accum_field (name='TDA', units='K', &
            desc='30-day average of 2-m temperature', &
            accum_type='timeavg', accum_period=-30, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       ! The following are running means.
       ! The accumulation period is set to -365 for a 365-day running mean.

       call init_accum_field (name='PREC365', units='MM H2O/S', &
            desc='365-day running mean of total precipitation', &
            accum_type='runmean', accum_period=-365, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       ! The following are accumulated fields.
       ! These types of fields are accumulated until a trigger value resets
       ! the accumulation to zero (see subroutine update_accum_field).
       ! Hence, [accper] is not valid.

       call init_accum_field (name='AGDDTW', units='K', &
            desc='growing degree-days base twmax', &
            accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='AGDD', units='K', &
            desc='growing degree-days base 5C', &
            accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0._r8)
    end if

    call init_accum_field (name='PREC60', units='MM H2O/S', &
         desc='60-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-60, &
         subgrid_type='pft', numlev=1, init_value=0._r8)
    
    call init_accum_field (name='PREC10', units='MM H2O/S', &
         desc='10-day running mean of total precipitation', &
         accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    if ( crop_prog )then
       ! 10-day average of min 2m temperature.

       call init_accum_field (name='TDM10', units='K', &
            desc='10-day running mean of min 2-m temperature', &
            accum_type='runmean', accum_period=-10, &
            subgrid_type='pft', numlev=1, init_value=SHR_CONST_TKFRZ)

       ! 5-day average of min 2m temperature.

       call init_accum_field (name='TDM5', units='K', &
            desc='5-day running mean of min 2-m temperature', &
            accum_type='runmean', accum_period=-5, &
            subgrid_type='pft', numlev=1, init_value=SHR_CONST_TKFRZ)

       ! All GDD summations are relative to the planting date
       ! (Kucharik & Brye 2003)

       call init_accum_field (name='GDD0', units='K', &
            desc='growing degree-days base 0C from planting', &
            accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDD8', units='K', &
            desc='growing degree-days base 8C from planting', &
            accum_type='runaccum', accum_period=not_used, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDD10', units='K', &
            desc='growing degree-days base 10C from planting', &
            accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDDPLANT', units='K', &
            desc='growing degree-days from planting', &
            accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='GDDTSOI', units='K', &
            desc='growing degree-days from planting (top two soil layers)', &
            accum_type='runaccum', accum_period=not_used,  &
            subgrid_type='pft', numlev=1, init_value=0._r8)
    end if

    ! Print output of accumulated fields

    call print_accum_fields()

  end subroutine initAccFlds

  !------------------------------------------------------------------------
  subroutine updateAccFlds(bounds)
    !
    ! !DESCRIPTION:
    ! Update and/or extract accumulated fields
    !
    ! !USES:
    use clmtype          , only : dgv_pftcon, pft, col, lun, grc
    use clmtype          , only : cps, ces, cws
    use clmtype          , only : pps, pes, pvs, pdgvs 
    use clm_atmlnd       , only : clm_a2l, a2l_downscaled_col
    use clm_varcon       , only : spval
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep, is_end_curr_day, get_curr_date
    use accumulMod       , only : update_accum_field, extract_accum_field, accumResetVal
    use pftvarcon        , only : nwcereal, nwcerealirrig, mxtmp, baset
    use pftvarcon        , only : ndllf_dcd_brl_tree
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: g,l,c,p                   ! indices
    integer :: itypveg                   ! vegetation type
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: year                      ! year (0, ...) for nstep
    integer :: month                     ! month (1, ..., 12) for nstep
    integer :: day                       ! day of month (1, ..., 31) for nstep
    integer :: secs                      ! seconds into current date for nstep
    logical :: end_cd                    ! temporary for is_end_curr_day() value
    integer :: ier                       ! error status
    real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level
    
!------------------------------------------------------------------------

    !    a2l_downscaled_col%forc_t     Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)                  
    !    a2l_downscaled_col%forc_rain  Input:  [real(r8) (:)]  rain rate [mm/s]                                  
    !    a2l_downscaled_col%forc_snow  Input:  [real(r8) (:)]  snow rate [mm/s]                                  
    !    clm_a2l%forc_solad 	       Input:  [real(r8) (:,:)]  direct beam radiation (visible only)            
    !    clm_a2l%forc_solai 	       Input:  [real(r8) (:,:)]  diffuse radiation     (visible only)            
    !    pps%croplive                  Input:  [logical (:)]  Flag, true if planted, not harvested               
    !    pps%vf                        Input:  [real(r8) (:)]  vernalization factor                              
    !    ces%t_soisno                  Input:  [real(r8) (:,:)]  soil temperature (K)                            
    !    cws%h2osoi_liq                Input:  [real(r8) (:,:)]  liquid water (kg/m2)                            
    !    cps%watsat                    Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
    !    cps%dz                        Input:  [real(r8) (:,:)]  layer thickness depth (m)                       
    !    pes%t_ref2m                   Input:  [real(r8) (:)]  2 m height surface air temperature (Kelvin)       
    !    pes%t_ref2m_u                 Input:  [real(r8) (:)]  Urban 2 m height surface air temperature (Kelvin) 
    !    pes%t_ref2m_r                 Input:  [real(r8) (:)]  Rural 2 m height surface air temperature (Kelvin) 
    !    pes%t_ref2m_max_inst          Output: [real(r8) (:)]  instantaneous daily max of average 2 m height surface air temp (K)
    !    pes%t_ref2m_min_inst          Output: [real(r8) (:)]  instantaneous daily min of average 2 m height surface air temp (K)
    !    pes%t_ref2m_max               Output: [real(r8) (:)]  daily maximum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_min               Output: [real(r8) (:)]  daily minimum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_max_u             Output: [real(r8) (:)]  Urban daily maximum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_max_r             Output: [real(r8) (:)]  Rural daily maximum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_min_u             Output: [real(r8) (:)]  Urban daily minimum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_min_r             Output: [real(r8) (:)]  Rural daily minimum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_max_inst_u        Output: [real(r8) (:)]  Urban instantaneous daily max of average 2 m height surface air temp (K)
    !    pes%t_ref2m_max_inst_r        Output: [real(r8) (:)]  Rural instantaneous daily max of average 2 m height surface air temp (K)
    !    pes%t_ref2m_min_inst_u        Output: [real(r8) (:)]  Urban instantaneous daily min of average 2 m height surface air temp (K)
    !    pes%t_ref2m_min_inst_r        Output: [real(r8) (:)]  Rural instantaneous daily min of average 2 m height surface air temp (K)
    !    pes%t10                       Output: [real(r8) (:)]  10-day running mean of the 2 m temperature (K)    
    !    pes%a10tmin                   Output: [real(r8) (:)]  10-day running mean of min 2-m temperature        
    !    pes%a5tmin                    Output: [real(r8) (:)]  5-day running mean of min 2-m temperature         
    !    pdgvs%t_mo                    Output: [real(r8) (:)]  30-day average temperature (Kelvin)               
    !    pdgvs%t_mo_min                Output: [real(r8) (:)]  annual min of t_mo (Kelvin)                       
    !    pdgvs%prec365                 Output: [real(r8) (:)]  365-day running mean of tot. precipitation        
    !    pdgvs%agddtw                  Output: [real(r8) (:)]  accumulated growing degree days above twmax       
    !    pdgvs%agdd                    Output: [real(r8) (:)]  accumulated growing degree days above 5           
    !    dgv_pftcon%twmax              Output: [real(r8) (:)]  upper limit of temperature of the warmest month   
    !    pps%prec60                    Output: [real(r8) (:)]  60-day running mean of tot. precipitation         
    !    pps%prec10                    Output: [real(r8) (:)]  10-day running mean of tot. precipitation         
    !    pps%gdd0                      Output: [real(r8) (:)]  growing degree-days base 0C'                      
    !    pps%gdd8                      Output: [real(r8) (:)]  growing degree-days base 8C from planting         
    !    pps%gdd10                     Output: [real(r8) (:)]  growing degree-days base 10C from planting        
    !    pps%gddplant                  Output: [real(r8) (:)]  growing degree-days from planting                 
    !    pps%gddtsoi                   Output: [real(r8) (:)]  growing degree-days from planting (top two soil layers)
    !    pvs%t_veg24                   Output: [real(r8) (:)]  24hr average vegetation temperature (K)           
    !    pvs%t_veg240                  Output: [real(r8) (:)]  240hr average vegetation temperature (Kelvin)     
    !    pvs%fsd24                     Output: [real(r8) (:)]  24hr average of direct beam radiation             
    !    pvs%fsd240                    Output: [real(r8) (:)]  240hr average of direct beam radiation            
    !    pvs%fsi24                     Output: [real(r8) (:)]  24hr average of diffuse beam radiation            
    !    pvs%fsi240                    Output: [real(r8) (:)]  240hr average of diffuse beam radiation           
    !    pvs%fsun24                    Output: [real(r8) (:)]  24hr average of sunlit fraction of canopy         
    !    pvs%fsun240                   Output: [real(r8) (:)]  240hr average of sunlit fraction of canopy        
    !    pvs%elai_p                    Output: [real(r8) (:)]  leaf area index average over timestep             
    !    pes%t_veg 	               Output: [real(r8) (:)]  pft vegetation temperature (Kelvin)               
    !    pps%fsun 	               Output: [real(r8) (:)]  sunlit fraction of canopy                         
    !    pps%elai 	               Output: [real(r8) (:)]  one-sided leaf area index with burying by snow    
    
    ! Determine calendar information

    dtime = get_step_size()
    nstep = get_nstep()
    call get_curr_date (year, month, day, secs)

    ! Don't do any accumulation if nstep is zero
    ! (only applies to coupled or cam mode)

    if (nstep == 0) return

    ! NOTE: currently only single level pft fields are used below
    ! Variables are declared above that should make it easy to incorporate
    ! multi-level or single-level fields of any subgrid type

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(bounds%begp:bounds%endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Accumulate and extract TREFAV - hourly average 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV', pes%t_ref2m, nstep)
    call extract_accum_field ('TREFAV', rbufslp, nstep)
    end_cd = is_end_curr_day()
    do p = bounds%begp,bounds%endp
       if (rbufslp(p) /= spval) then
          pes%t_ref2m_max_inst(p) = max(rbufslp(p), pes%t_ref2m_max_inst(p))
          pes%t_ref2m_min_inst(p) = min(rbufslp(p), pes%t_ref2m_min_inst(p))
       endif
       if (end_cd) then
          pes%t_ref2m_max(p) = pes%t_ref2m_max_inst(p)
          pes%t_ref2m_min(p) = pes%t_ref2m_min_inst(p)
          pes%t_ref2m_max_inst(p) = -spval
          pes%t_ref2m_min_inst(p) =  spval
       else if (secs == int(dtime)) then
          pes%t_ref2m_max(p) = spval
          pes%t_ref2m_min(p) = spval
       endif
    end do

    ! Accumulate and extract TREFAV_U - hourly average urban 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_U', pes%t_ref2m_u, nstep)
    call extract_accum_field ('TREFAV_U', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       if (rbufslp(p) /= spval) then
          pes%t_ref2m_max_inst_u(p) = max(rbufslp(p), pes%t_ref2m_max_inst_u(p))
          pes%t_ref2m_min_inst_u(p) = min(rbufslp(p), pes%t_ref2m_min_inst_u(p))
       endif
       if (end_cd) then
         if (lun%urbpoi(l)) then
          pes%t_ref2m_max_u(p) = pes%t_ref2m_max_inst_u(p)
          pes%t_ref2m_min_u(p) = pes%t_ref2m_min_inst_u(p)
          pes%t_ref2m_max_inst_u(p) = -spval
          pes%t_ref2m_min_inst_u(p) =  spval
         end if
       else if (secs == int(dtime)) then
          pes%t_ref2m_max_u(p) = spval
          pes%t_ref2m_min_u(p) = spval
       endif
    end do

    ! Accumulate and extract TREFAV_R - hourly average rural 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called

    call update_accum_field  ('TREFAV_R', pes%t_ref2m_r, nstep)
    call extract_accum_field ('TREFAV_R', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       if (rbufslp(p) /= spval) then
          pes%t_ref2m_max_inst_r(p) = max(rbufslp(p), pes%t_ref2m_max_inst_r(p))
          pes%t_ref2m_min_inst_r(p) = min(rbufslp(p), pes%t_ref2m_min_inst_r(p))
       endif
       if (end_cd) then
         if (.not.(lun%ifspecial(l))) then
          pes%t_ref2m_max_r(p) = pes%t_ref2m_max_inst_r(p)
          pes%t_ref2m_min_r(p) = pes%t_ref2m_min_inst_r(p)
          pes%t_ref2m_max_inst_r(p) = -spval
          pes%t_ref2m_min_inst_r(p) =  spval
         end if
       else if (secs == int(dtime)) then
          pes%t_ref2m_max_r(p) = spval
          pes%t_ref2m_min_r(p) = spval
       endif
    end do

    ! Accumulate and extract T_VEG24 & T_VEG240 
    do p = bounds%begp,bounds%endp
       rbufslp(p) = pes%t_veg(p)
    end do
    call update_accum_field  ('T_VEG24' , rbufslp, nstep)
    call extract_accum_field ('T_VEG24' , pvs%t_veg24, nstep)
    call update_accum_field  ('T_VEG240', rbufslp, nstep)
    call extract_accum_field ('T_VEG240', pvs%t_veg240, nstep)

    ! Accumulate and extract forc_solad24 & clm_a2l%forc_solad240 
    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)
       rbufslp(p) = clm_a2l%forc_solad(g,1)
    end do
    call update_accum_field  ('FSD240', rbufslp, nstep)
    call extract_accum_field ('FSD240', pvs%fsd240, nstep)
    call update_accum_field  ('FSD24' , rbufslp, nstep)
    call extract_accum_field ('FSD24' , pvs%fsd24, nstep)

    ! Accumulate and extract forc_solai24 & forc_solai240 
    do p = bounds%begp,bounds%endp
       g = pft%gridcell(p)
       rbufslp(p) = clm_a2l%forc_solai(g,1)
    end do
    call update_accum_field  ('FSI24' , rbufslp, nstep)
    call extract_accum_field ('FSI24' , pvs%fsi24, nstep)
    call update_accum_field  ('FSI240', rbufslp, nstep)
    call extract_accum_field ('FSI240', pvs%fsi240, nstep)

    ! Accumulate and extract fsun24 & fsun240   
    do p = bounds%begp,bounds%endp
       rbufslp(p) = pps%fsun(p)
    end do
    call update_accum_field  ('FSUN24' , rbufslp, nstep)
    call extract_accum_field ('FSUN24' , pvs%fsun24, nstep)
    call update_accum_field  ('FSUN240', rbufslp, nstep)
    call extract_accum_field ('FSUN240', pvs%fsun240, nstep)

    ! Accumulate and extract elai_p 
    do p = bounds%begp,bounds%endp
       rbufslp(p) = pps%elai(p)
    end do
    call update_accum_field  ('LAIP', rbufslp, nstep)
    call extract_accum_field ('LAIP', pvs%elai_p, nstep)

    ! Accumulate and extract T10
    !(acumulates TSA as 10-day running mean)

    call update_accum_field  ('T10', pes%t_ref2m, nstep)
    call extract_accum_field ('T10', pes%t10, nstep)

    if (use_cndv) then
       ! Accumulate and extract TDA
       ! (accumulates TBOT as 30-day average)
       ! Also determine t_mo_min
       
       do p = bounds%begp,bounds%endp
          c = pft%column(p)
          rbufslp(p) = a2l_downscaled_col%forc_t(c)
       end do
       call update_accum_field  ('TDA', rbufslp, nstep)
       call extract_accum_field ('TDA', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pdgvs%t_mo(p) = rbufslp(p)
          pdgvs%t_mo_min(p) = min(pdgvs%t_mo_min(p), rbufslp(p))
       end do

       ! Accumulate and extract PREC365
       ! (accumulates total precipitation as 365-day running mean)

       do p = bounds%begp,bounds%endp
          c = pft%column(p)
          rbufslp(p) = a2l_downscaled_col%forc_rain(c) + a2l_downscaled_col%forc_snow(c)
       end do
       call update_accum_field  ('PREC365', rbufslp, nstep)
       call extract_accum_field ('PREC365', pdgvs%prec365, nstep)

       ! Accumulate growing degree days based on 10-day running mean temperature.
       ! The trigger to reset the accumulated values to zero is -99999.

       ! Accumulate and extract AGDDTW (gdd base twmax, which is 23 deg C
       ! for boreal woody pfts)

       do p = bounds%begp,bounds%endp
          rbufslp(p) = max(0._r8, (pes%t10(p) - SHR_CONST_TKFRZ - dgv_pftcon%twmax(ndllf_dcd_brl_tree)) &
               * dtime/SHR_CONST_CDAY)
          if (month==1 .and. day==1 .and. secs==int(dtime)) rbufslp(p) = accumResetVal
       end do
       call update_accum_field  ('AGDDTW', rbufslp, nstep)
       call extract_accum_field ('AGDDTW', pdgvs%agddtw, nstep)

       ! Accumulate and extract AGDD

       do p = bounds%begp,bounds%endp
          rbufslp(p) = max(0.0_r8, (pes%t_ref2m(p) - (SHR_CONST_TKFRZ + 5.0_r8)) &
               * dtime/SHR_CONST_CDAY)
          !
          ! Fix (for bug 1858) from Sam Levis to reset the annual AGDD variable
          ! 
          if (month==1 .and. day==1 .and. secs==int(dtime)) rbufslp(p) = accumResetVal
       end do
       call update_accum_field  ('AGDD', rbufslp, nstep)
       call extract_accum_field ('AGDD', pdgvs%agdd, nstep)
    end if
    
    do p = bounds%begp,bounds%endp
       c = pft%column(p)
       rbufslp(p) = a2l_downscaled_col%forc_rain(c) + a2l_downscaled_col%forc_snow(c)
    end do
    call update_accum_field  ('PREC60', rbufslp, nstep)
    call extract_accum_field ('PREC60', pps%prec60, nstep)

    ! Accumulate and extract PREC10
    ! (accumulates total precipitation as 10-day running mean)
     do p = bounds%begp,bounds%endp
       c = pft%column(p)
       rbufslp(p) = a2l_downscaled_col%forc_rain(c) + a2l_downscaled_col%forc_snow(c)
    end do
    call update_accum_field  ('PREC10', rbufslp, nstep)
    call extract_accum_field ('PREC10', pps%prec10, nstep)
 
    if ( crop_prog )then
       ! Accumulate and extract TDM10

       do p = bounds%begp,bounds%endp
          rbufslp(p) = min(pes%t_ref2m_min(p),pes%t_ref2m_min_inst(p)) !slevis: ok choice?
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                                     !'min_inst' not initialized?
       call update_accum_field  ('TDM10', rbufslp, nstep)
       call extract_accum_field ('TDM10', pes%a10tmin, nstep)

       ! Accumulate and extract TDM5

       do p = bounds%begp,bounds%endp
          rbufslp(p) = min(pes%t_ref2m_min(p),pes%t_ref2m_min_inst(p)) !slevis: ok choice?
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                         !'min_inst' not initialized?
       call update_accum_field  ('TDM5', rbufslp, nstep)
       call extract_accum_field ('TDM5', pes%a5tmin, nstep)

       ! Accumulate and extract GDD0

       do p = bounds%begp,bounds%endp
          itypveg = pft%itype(p)
          g = pft%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(26._r8, pes%t_ref2m(p)-SHR_CONST_TKFRZ)) &
                  * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD0', rbufslp, nstep)
       call extract_accum_field ('GDD0', pps%gdd0, nstep)

       ! Accumulate and extract GDD8

       do p = bounds%begp,bounds%endp
          itypveg = pft%itype(p)
          g = pft%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                  pes%t_ref2m(p)-(SHR_CONST_TKFRZ + 8._r8))) &
                  * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD8', rbufslp, nstep)
       call extract_accum_field ('GDD8', pps%gdd8, nstep)

       ! Accumulate and extract GDD10

       do p = bounds%begp,bounds%endp
          itypveg = pft%itype(p)
          g = pft%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                  pes%t_ref2m(p)-(SHR_CONST_TKFRZ + 10._r8))) &
                  * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       call update_accum_field  ('GDD10', rbufslp, nstep)
       call extract_accum_field ('GDD10', pps%gdd10, nstep)

       ! Accumulate and extract GDDPLANT

       do p = bounds%begp,bounds%endp
          if (pps%croplive(p)) then ! relative to planting date
             itypveg = pft%itype(p)
             rbufslp(p) = max(0._r8, min(mxtmp(itypveg), &
                  pes%t_ref2m(p)-(SHR_CONST_TKFRZ + baset(itypveg)))) &
                          * dtime/SHR_CONST_CDAY
             if (itypveg == nwcereal .or. itypveg == nwcerealirrig) rbufslp(p) = rbufslp(p)*pps%vf(p)
          else
             rbufslp(p) = accumResetVal
          end if
       end do
       call update_accum_field  ('GDDPLANT', rbufslp, nstep)
       call extract_accum_field ('GDDPLANT', pps%gddplant, nstep)

       ! Accumulate and extract GDDTSOI
       ! In agroibis this variable is calculated
       ! to 0.05 m, so here we use the top two soil layers
   
       do p = bounds%begp,bounds%endp
          if (pps%croplive(p)) then ! relative to planting date
             itypveg = pft%itype(p)
             c = pft%column(p)
             rbufslp(p) = max(0._r8, min(mxtmp(itypveg), &
                  ((ces%t_soisno(c,1)*cps%dz(c,1) + &
                  ces%t_soisno(c,2)*cps%dz(c,2))/(cps%dz(c,1)+cps%dz(c,2))) - &
                  (SHR_CONST_TKFRZ + baset(itypveg)))) * dtime/SHR_CONST_CDAY
             if (itypveg == nwcereal .or. itypveg == nwcerealirrig) rbufslp(p) = rbufslp(p)*pps%vf(p)
          else
             rbufslp(p) = accumResetVal
          end if
       end do
       call update_accum_field  ('GDDTSOI', rbufslp, nstep)
       call extract_accum_field ('GDDTSOI', pps%gddtsoi, nstep)

    end if

    ! Deallocate dynamic memory

    deallocate(rbufslp)

  end subroutine updateAccFlds

  !-----------------------------------------------------------------------
  subroutine initAccClmtype(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize clmtype variables that are associated with
    ! time accumulated fields. This routine is called in an initial run
    ! at nstep=0 for cam and csm mode.
    ! This routine is also always called for a restart run and
    ! therefore must be called after the restart file is read in
    ! and the accumulated fields are obtained.
    !
    ! !USES:
    use clmtype
    use accumulMod  , only : extract_accum_field
    use clm_time_manager, only : get_nstep
    use clm_varctl  , only : nsrest, nsrStartup
    use clm_varcon  , only : spval
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: nstep        ! time step
    integer :: ier          ! error status
    real(r8), pointer :: rbufslp(:)  ! temporary
    character(len=32) :: subname = 'initAccClmtype'  ! subroutine name
    !-----------------------------------------------------------------------

    !    pes%t_ref2m_max_inst      Output: [real(r8) (:)]  instantaneous daily max of average 2 m height surface air temp (K)
    !    pes%t_ref2m_min_inst      Output: [real(r8) (:)]  instantaneous daily min of average 2 m height surface air temp (K)
    !    pes%t_ref2m_max           Output: [real(r8) (:)]  daily maximum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_min           Output: [real(r8) (:)]  daily minimum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_max_inst_u    Output: [real(r8) (:)]  Urban instantaneous daily max of average 2 m height surface air temp (K)
    !    pes%t_ref2m_max_inst_r    Output: [real(r8) (:)]  Rural instantaneous daily max of average 2 m height surface air temp (K)
    !    pes%t_ref2m_min_inst_u    Output: [real(r8) (:)]  Urban instantaneous daily min of average 2 m height surface air temp (K)
    !    pes%t_ref2m_min_inst_r    Output: [real(r8) (:)]  Rural instantaneous daily min of average 2 m height surface air temp (K)
    !    pes%t_ref2m_max_u         Output: [real(r8) (:)]  Urban daily maximum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_max_r         Output: [real(r8) (:)]  Rural daily maximum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_min_u         Output: [real(r8) (:)]  Urban daily minimum of average 2 m height surface air temperature (K)
    !    pes%t_ref2m_min_r         Output: [real(r8) (:)]  Rural daily minimum of average 2 m height surface air temperature (K)
    !    pes%t10                   Output: [real(r8) (:)]  10-day running mean of the 2 m temperature (K)    
    !    pes%a10tmin               Output: [real(r8) (:)]  10-day running mean of min 2-m temperature        
    !    pes%a5tmin                Output: [real(r8) (:)]  5-day running mean of min 2-m temperature         
    !    pdgvs%t_mo                Output: [real(r8) (:)]  30-day average temperature (Kelvin)               
    !    pdgvs%prec365             Output: [real(r8) (:)]  365-day running mean of tot. precipitation        
    !    pdgvs%agddtw              Output: [real(r8) (:)]  accumulated growing degree days above twmax       
    !    pdgvs%agdd                Output: [real(r8) (:)]  accumulated growing degree days above 5           
    !    pps%prec60                Output: [real(r8) (:)]  60-day running mean of tot. precipitation         
    !    pps%prec10                Output: [real(r8) (:)]  10-day running mean of tot. precipitation         
    !    pps%gdd0                  Output: [real(r8) (:)]  growing degree-days base 0C'                      
    !    pps%gdd8                  Output: [real(r8) (:)]  growing degree-days base 8C from planting         
    !    pps%gdd10                 Output: [real(r8) (:)]  growing degree-days base 10C from planting        
    !    pps%gddplant              Output: [real(r8) (:)]  growing degree-days from planting                 
    !    pps%gddtsoi               Output: [real(r8) (:)]  growing degree-days from planting (top two soil layers)
    !    pvs%t_veg24               Output: [real(r8) (:)]  24hr average vegetation temperature (K)           
    !    pvs%t_veg240              Output: [real(r8) (:)]  240hr average vegetation temperature (Kelvin)     
    !    pvs%fsd24                 Output: [real(r8) (:)]  24hr average of direct beam radiation             
    !    pvs%fsd240                Output: [real(r8) (:)]  240hr average of direct beam radiation            
    !    pvs%fsi24                 Output: [real(r8) (:)]  24hr average of diffuse beam radiation            
    !    pvs%fsi240                Output: [real(r8) (:)]  240hr average of diffuse beam radiation           
    !    pvs%fsun24                Output: [real(r8) (:)]  24hr average of sunlit fraction of canopy         
    !    pvs%fsun240               Output: [real(r8) (:)]  240hr average of sunlit fraction of canopy        
    !    pvs%elai_p                Output: [real(r8) (:)]  leaf area index average over timestep             
    
    ! Determine time step

    nstep = get_nstep()

    ! Initialize 2m ref temperature max and min values

    if (nsrest == nsrStartup) then ! Why not restart&branch? These vars are not in clmr.
       do p = bounds%begp,bounds%endp
          pes%t_ref2m_max(p) = spval
          pes%t_ref2m_min(p) = spval
          pes%t_ref2m_max_inst(p) = -spval
          pes%t_ref2m_min_inst(p) =  spval
          pes%t_ref2m_max_u(p) = spval
          pes%t_ref2m_max_r(p) = spval
          pes%t_ref2m_min_u(p) = spval
          pes%t_ref2m_min_r(p) = spval
          pes%t_ref2m_max_inst_u(p) = -spval
          pes%t_ref2m_max_inst_r(p) = -spval
          pes%t_ref2m_min_inst_u(p) =  spval
          pes%t_ref2m_min_inst_r(p) =  spval
       end do
    end if

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(bounds%begp:bounds%endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'extract_accum_hist allocation error for rbufslp in '//subname
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Initialize clmtype variables that are to be time accumulated

    call extract_accum_field ('T_VEG24', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pvs%t_veg24(p) = rbufslp(p)
    end do

    call extract_accum_field ('T_VEG240', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pvs%t_veg240(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSD24', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pvs%fsd24(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSD240', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pvs%fsd240(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSI24', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pvs%fsi24(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSI240', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pvs%fsi240(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSUN24', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pvs%fsun24(p) = rbufslp(p)
    end do

    call extract_accum_field ('FSUN240', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pvs%fsun240(p) = rbufslp(p)
    end do

    call extract_accum_field ('LAIP', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pvs%elai_p(p) = rbufslp(p)
    end do

    if ( crop_prog )then

       call extract_accum_field ('GDD0', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pps%gdd0(p) = rbufslp(p)
       end do
   
       call extract_accum_field ('GDD8', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pps%gdd8(p) = rbufslp(p)
       end do

       call extract_accum_field ('GDD10', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pps%gdd10(p) = rbufslp(p)
       end do

       call extract_accum_field ('GDDPLANT', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pps%gddplant(p) = rbufslp(p)
       end do

       call extract_accum_field ('GDDTSOI', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pps%gddtsoi(p) = rbufslp(p)
       end do

       call extract_accum_field ('TDM10', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pes%a10tmin(p) = rbufslp(p)
       end do

       call extract_accum_field ('TDM5', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pes%a5tmin(p) = rbufslp(p)
       end do

    end if

    call extract_accum_field ('T10', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pes%t10(p) = rbufslp(p)
    end do

    if (use_cndv) then
       call extract_accum_field ('TDA', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pdgvs%t_mo(p) = rbufslp(p)
       end do

       call extract_accum_field ('PREC365', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pdgvs%prec365(p) = rbufslp(p)
       end do

       call extract_accum_field ('AGDDTW', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pdgvs%agddtw(p) = rbufslp(p)
       end do

       call extract_accum_field ('AGDD', rbufslp, nstep)
       do p = bounds%begp,bounds%endp
          pdgvs%agdd(p) = rbufslp(p)
       end do
    end if

    call extract_accum_field ('PREC60', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pps%prec60(p) = rbufslp(p)
    end do

    call extract_accum_field ('PREC10', rbufslp, nstep)
    do p = bounds%begp,bounds%endp
       pps%prec10(p) = rbufslp(p)
    end do

    deallocate(rbufslp)

  end subroutine initAccClmtype

end module accFldsMod
