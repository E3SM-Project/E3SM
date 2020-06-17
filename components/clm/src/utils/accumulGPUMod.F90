module accumulGPUMod

  !module that contains the class methods needed to update accumulators
  use accumulMod       , only : accum_type, accumResetVal, naccflds
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use atm2lndType      , only : atm2lnd_type
  use TopounitDataType , only : topounit_atmospheric_flux
  use decompMod        , only : bounds_type
  use CropType         , only : crop_type
  use CNStateType      , only : cnstate_type
  use CanopyStateType  , only : canopystate_type
  use VegetationDataType, only : vegetation_energy_state
  use VegetationType        , only : veg_pp
  use LandunitType      , only : lun_pp
  use GridcellType      , only : grc_pp
  use clm_varpar        , only : crop_prog

  type accum_field_gpu_type
     character(len=  8), pointer :: name    => null()  !field name
     character(len=128), pointer :: desc    => null()  !field description
     character(len=  8), pointer :: units   => null()  !field units
     character(len=  8), pointer :: acctype => null()  !accumulation type: ["timeavg","runmean","runaccum"]
     character(len=  8), pointer :: type1d  => null()  !subgrid type: ["gridcell","topounit","landunit","column" or "pft"]
     character(len=  8), pointer :: type2d  => null()  !type2d ('','levsoi','numrad',..etc. )
     integer , pointer           :: beg1d   => null() !subgrid type beginning index
     integer , pointer           :: end1d   => null() !subgrid type ending index
     integer , pointer           :: num1d   => null() !total subgrid points
     integer , pointer           :: numlev  => null() !number of vertical levels in field
     real(r8), pointer           :: initval => null() !initial value of accumulated field
     real(r8), pointer           :: val(:,:) => null()         !accumulated field
     integer , pointer           :: period  => null()  !field accumulation period (in model time steps)
  end type accum_field_gpu_type

  type(accum_field_gpu_type), public, allocatable :: accum_field_gpu(:)
  !$acc declare create(accum_field_gpu)

  public :: init_accum_gpu

  interface extract_accum_field_gpu
     module procedure extract_accum_field_sl ! Extract current val of single-level accumulator field
     module procedure extract_accum_field_ml ! Extract current val of multi-level accumulator field
  end interface

  interface update_accum_field_gpu           ! Updates the current value of an accumulator field
     module procedure update_accum_field_sl  ! Update single-level accumulator field
     module procedure update_accum_field_ml  ! Update multi-level accumulator field
  end interface

  contains
    subroutine init_accum_gpu()
      use accumulMod, only : accum, naccflds
      implicit none

      integer :: nf
      integer :: i, j, lb, ub

      print *, "allocating ",naccflds,"for gpu accum"
      allocate(accum_field_gpu(naccflds))

      do nf = 1, naccflds
        allocate(accum_field_gpu(nf)%name    ); accum_field_gpu(nf)%name   = accum(nf)%name
        allocate(accum_field_gpu(nf)%desc    ); accum_field_gpu(nf)%desc   = accum(nf)%desc
        allocate(accum_field_gpu(nf)%units   ); accum_field_gpu(nf)%units  = accum(nf)%units
        allocate(accum_field_gpu(nf)%acctype ); accum_field_gpu(nf)%acctype= accum(nf)%acctype
        allocate(accum_field_gpu(nf)%type1d  ); accum_field_gpu(nf)%type1d = accum(nf)%type1d
        allocate(accum_field_gpu(nf)%type2d  ); accum_field_gpu(nf)%type2d = accum(nf)%type2d
        allocate(accum_field_gpu(nf)%beg1d   ); accum_field_gpu(nf)%beg1d  = accum(nf)%beg1d
        allocate(accum_field_gpu(nf)%end1d   ); accum_field_gpu(nf)%end1d  = accum(nf)%end1d
        allocate(accum_field_gpu(nf)%num1d   ); accum_field_gpu(nf)%num1d  = accum(nf)%num1d
        allocate(accum_field_gpu(nf)%numlev  ); accum_field_gpu(nf)%numlev = accum(nf)%numlev
        allocate(accum_field_gpu(nf)%initval ); accum_field_gpu(nf)%initval= accum(nf)%initval
        allocate(accum_field_gpu(nf)%period  ); accum_field_gpu(nf)%period = accum(nf)%period
        lb = accum_field_gpu(nf)%beg1d
        ub = accum_field_gpu(nf)%end1d
        allocate(accum_field_gpu(nf)%val(lb:ub ,accum(nf)%numlev ) )
        accum_field_gpu(nf)%val(:,:) = accum(nf)%val(:,:)
        print *, nf, accum_field_gpu(nf)%name
      end do

      !$acc enter data copyin(naccflds,accum_field_gpu)

    end subroutine init_accum_gpu

    !------------------------------------------------------------------------
    subroutine update_accum_field_sl (nf, name, beg1d, end1d, field, nstep)
      !
      ! !DESCRIPTION:
      ! Accumulate single level field over specified time interval.
      ! The appropriate field is accumulated in the array [accval].
      !$acc routine seq
      ! !ARGUMENTS:
      implicit none
      integer, value,intent(in)  :: nf
      character(len=8), intent(in) :: name     !field name
      integer,value, intent(in) :: beg1d
      integer,value, intent(in) :: end1d
      real(r8), intent(inout) :: field(beg1d : ) !field values for current time step
      integer ,value, intent(in) :: nstep            !time step index
      !
      ! !LOCAL VARIABLES:
      integer :: i,k              !indices
      integer :: accper              !temporary accumulation period
      !integer :: beg1d,end1d             !subgrid beginning,ending indices
      !------------------------------------------------------------------------
      ! error check
      ! beg1d = accum_field_gpu(nf)%beg1d
      ! end1d = accum_field_gpu(nf)%end1d

      ! reset accumulated field value if necessary and  update
      ! accumulation field
      ! running mean never reset

      if (accum_field_gpu(nf)%acctype == 'timeavg') then

         !time average field reset every accumulation period
         !normalize at end of accumulation period

         if ((mod(nstep,accum_field_gpu(nf)%period) == 1 .or. accum_field_gpu(nf)%period == 1) .and. (nstep /= 0))then
          do k = beg1d, end1d
             accum_field_gpu(nf)%val(k,1) = 0._r8
          end do
         end if
         do k = beg1d, end1d
              accum_field_gpu(nf)%val(k,1) =  accum_field_gpu(nf)%val(k,1) + field(k)
         end do
         if (mod(nstep,accum_field_gpu(nf)%period) == 0) then
           do k = beg1d, end1d
              accum_field_gpu(nf)%val(k,1) = accum_field_gpu(nf)%val(k,1) / accum_field_gpu(nf)%period
           end do
         endif

      else if (accum_field_gpu(nf)%acctype == 'runmean') then

         !running mean - reset accumulation period until greater than nstep

         accper = min(nstep,accum_field_gpu(nf)%period)
         do k = beg1d, end1d
           accum_field_gpu(nf)%val(k,1) = ((accper-1)*accum_field_gpu(nf)%val(k,1) + field(k)) / accper
         end do
      else if (accum_field_gpu(nf)%acctype == 'runaccum') then

         !running accumulation field reset at trigger -99999

         do k = beg1d,end1d
            if (nint(field(k)) == -99999) then
               accum_field_gpu(nf)%val(k,1) = 0._r8
            end if
            accum_field_gpu(nf)%val(k,1) = min(max(accum_field_gpu(nf)%val(k,1) + field(k), 0._r8), 99999._r8)

         end do

      end if

    end subroutine update_accum_field_sl

    !------------------------------------------------------------------------
    subroutine update_accum_field_ml (nf, name, beg, end,field, nstep)
      !
      ! !DESCRIPTION:
      ! Accumulate multi level field over specified time interval.
      !$acc routine seq
      ! !ARGUMENTS:
      implicit none
      integer, value, intent(in)   :: nf
      character(len=8), intent(in) :: name       !field name
      integer,value, intent(in) :: beg
      integer,value, intent(in) :: end
      real(r8),  intent(inout)  :: field(beg: ,:) !field values for current time step
      integer , value , intent(in) :: nstep              !time step index
      !
      ! !LOCAL VARIABLES:
      integer :: i,j,k            !indices
      integer :: accper              !temporary accumulation period
      integer :: numlev              !number of vertical levels
      !------------------------------------------------------------------------

      numlev = accum_field_gpu(nf)%numlev
      ! beg = accum_field_gpu(nf)%beg1d
      ! end = accum_field_gpu(nf)%end1d

      ! reset accumulated field value if necessary and  update
      ! accumulation field
      ! running mean never reset

      if (accum_field_gpu(nf)%acctype == 'timeavg') then

         !time average field reset every accumulation period
         !normalize at end of accumulation period

         if ((mod(nstep,accum_field_gpu(nf)%period) == 1 .or. accum_field_gpu(nf)%period == 1) .and. (nstep /= 0))then
            accum_field_gpu(nf)%val(beg:end,1:numlev) = 0._r8
         endif
         accum_field_gpu(nf)%val(beg:end,1:numlev) =  accum_field_gpu(nf)%val(beg:end,1:numlev) + field(beg:end,1:numlev)
         if (mod(nstep,accum_field_gpu(nf)%period) == 0) then
            accum_field_gpu(nf)%val(beg:end,1:numlev) = accum_field_gpu(nf)%val(beg:end,1:numlev) / accum_field_gpu(nf)%period
         endif

      else if (accum_field_gpu(nf)%acctype == 'runmean') then

         !running mean - reset accumulation period until greater than nstep

         accper = min (nstep,accum_field_gpu(nf)%period)
         accum_field_gpu(nf)%val(beg:end,1:numlev) = &
              ((accper-1)*accum_field_gpu(nf)%val(beg:end,1:numlev) + field(beg:end,1:numlev)) / accper

      else if (accum_field_gpu(nf)%acctype == 'runaccum') then

         !running accumulation field reset at trigger -99999

         do j = 1,numlev
            do k = beg,end
               if (nint(field(k,j)) == -99999) then
                  accum_field_gpu(nf)%val(k,j) = 0._r8
               end if
            end do
         end do
         accum_field_gpu(nf)%val(beg:end,1:numlev) = &
              min(max(accum_field_gpu(nf)%val(beg:end,1:numlev) + field(beg:end,1:numlev), 0._r8), 99999._r8)

      end if

    end subroutine update_accum_field_ml


    !------------------------------------------------------------------------
    subroutine extract_accum_field_sl (nf, name,beg1d,end1d, field, nstep)
      !
      ! !DESCRIPTION:
      ! Extract single-level accumulated field.
      ! This routine extracts the field values from the multi-level
      ! accumulation field. It extracts the current value except if
      ! the field type is a time average. In this case, an absurd value
      ! is assigned to  indicate the time average is not yet valid.
      !$acc routine seq
      ! !USES:
      use clm_varcon, only : spval, ispval
      !
      ! !ARGUMENTS:
      implicit none
      integer,value, intent(in)  :: nf
      character(len=8), intent(in) :: name     !field name
      integer,value, intent(in) :: beg1d
      integer,value, intent(in) :: end1d
      real(r8),  intent(inout)  :: field(beg1d : ) !field values for current time step
      integer ,value, intent(in) :: nstep            !timestep index
      !
      ! !LOCAL VARIABLES:
      integer :: i,k         !indices
      !------------------------------------------------------------------------
      ! extract field

      if (accum_field_gpu(nf)%acctype == 'timeavg' .and. &
           mod(nstep,accum_field_gpu(nf)%period) /= 0) then
         do k = beg1d,end1d
            field(k) = spval  !assign absurd value when avg not ready
         end do
      else
         do k = beg1d,end1d
            field(k) = accum_field_gpu(nf)%val(k,1)
         end do
      end if

    end subroutine extract_accum_field_sl

    !------------------------------------------------------------------------
    subroutine extract_accum_field_ml (nf, name,beg,end, field, nstep)
      !
      ! !DESCRIPTION:
      ! Extract mutli-level accumulated field.
      ! This routine extracts the field values from the multi-level
      ! accumulation field. It extracts the current value except if
      ! the field type is a time average. In this case, an absurd value
      ! is assigned to  indicate the time average is not yet valid.
      !$acc routine seq
      ! !USES:
      use clm_varcon, only : spval
      !
      ! !ARGUMENTS:
      implicit none
      integer, value,intent(in)    :: nf
      character(len=8), intent(in) :: name       !field name
      integer,value, intent(in) :: beg
      integer,value, intent(in) :: end
      real(r8),  intent(inout)  :: field(beg :,: ) !field values for current time step
      integer,value, intent(in) :: nstep               !timestep index
      !
      ! !LOCAL VARIABLES:
      integer :: i,j,k           !indices
      integer :: numlev          !number of vertical levels
      !------------------------------------------------------------------------

      numlev = accum_field_gpu(nf)%numlev
      ! beg = accum_field_gpu(nf)%beg1d
      ! end = accum_field_gpu(nf)%end1d

      !extract field

      if (accum_field_gpu(nf)%acctype == 'timeavg' .and. &
           mod(nstep,accum_field_gpu(nf)%period) /= 0) then
         do j = 1,numlev
            do k = beg,end
               field(k,j) = spval  !assign absurd value when avg not ready
            end do
         end do
      else
         do j = 1,numlev
            do k = beg,end
               field(k,j) = accum_field_gpu(nf)%val(k,j)
            end do
         end do
      end if

    end subroutine extract_accum_field_ml

     subroutine atm2lnd_UpdateAccVars (atm2lnd_vars, bounds, nstep)
       !
       ! USES
       !$acc routine seq
       ! !ARGUMENTS:
       type(atm2lnd_type)  , intent(inout) :: atm2lnd_vars
       type(bounds_type)   , intent(in)    :: bounds
       integer         , value ,intent(in) :: nstep          ! timestep number
       !
       ! !LOCAL VARIABLES:
       integer :: g,c,p                     ! indices
       integer :: begp, endp
       character(len=8) :: name
       real(r8) :: rbufslp(bounds%begp:bounds%endp)      ! temporary single level - pft level
       !---------------------------------------------------------------------
       associate( &
         fsd240_patch => atm2lnd_vars%fsd240_patch  , &
         fsd24_patch  => atm2lnd_vars%fsd24_patch   , &
         fsi24_patch  => atm2lnd_vars%fsi24_patch   , &
         fsi240_patch => atm2lnd_vars%fsi240_patch  , &
         prec60_patch => atm2lnd_vars%prec60_patch  , &
         prec10_patch => atm2lnd_vars%prec10_patch   &
         )
       begp = bounds%begp; endp = bounds%endp

       ! Allocate needed dynamic memory for single level pft field
      ! Accumulate and extract forc_solad24 & forc_solad240
       do p = begp,endp
          g = veg_pp%gridcell(p)
          rbufslp(p) = atm2lnd_vars%forc_solad_grc(g,1)
       end do
       name = 'FSD240'
       call update_accum_field_gpu  (2,name,begp,endp, rbufslp(begp:endp), nstep)
       call extract_accum_field_gpu (2,name,begp,endp, fsd240_patch(begp:endp), nstep)
       name = 'FSD24'
       call update_accum_field_gpu  (1, name ,begp,endp, rbufslp(begp:endp), nstep)
       call extract_accum_field_gpu (1, name ,begp,endp, fsd24_patch(begp:endp), nstep)
           ! Accumulate and extract forc_solai24 & forc_solai240
       do p = begp,endp
          g = veg_pp%gridcell(p)
          rbufslp(p) = atm2lnd_vars%forc_solai_grc(g,1)
       end do
       name = 'FSI24'
       call update_accum_field_gpu  (3, name, begp,endp, rbufslp(begp:endp) , nstep)
       call extract_accum_field_gpu (3, name, begp,endp, fsi24_patch(begp:endp) , nstep)
       name = 'FSI240'
       call update_accum_field_gpu  (4, name,begp,endp, rbufslp(begp:endp) , nstep)
       call extract_accum_field_gpu (4, name,begp,endp, fsi240_patch(begp:endp), nstep)
       do p = begp,endp
          c = veg_pp%column(p)
          rbufslp(p) = atm2lnd_vars%forc_rain_downscaled_col(c) + atm2lnd_vars%forc_snow_downscaled_col(c)
       end do

       if (use_cn) then
          name = 'PREC60'
          ! Accumulate and extract PREC60 (accumulates total precipitation as 60-day running mean)
          call update_accum_field_gpu  (6,name,begp,endp, rbufslp(begp:endp), nstep)
          call extract_accum_field_gpu (6,name,begp,endp, prec60_patch(begp:endp), nstep)

          ! Accumulate and extract PREC10 (accumulates total precipitation as 10-day running mean)
          name = 'PREC10'
          call update_accum_field_gpu  (5, name,begp,endp, rbufslp(begp:endp), nstep)
          call extract_accum_field_gpu (5, name,begp,endp, prec10_patch(begp:endp), nstep)
       end if
       if (use_fates) then
         name = 'PREC24'
          !call update_accum_field_gpu  (name, rbufslp, nstep)
          !call extract_accum_field_gpu (name, atm2lnd_vars%prec24_patch, nstep)
          do p = bounds%begp, bounds%endp
             c = veg_pp%column(p)
             g = veg_pp%gridcell(p)
             rbufslp(p) = atm2lnd_vars%forc_wind_grc(g)
          end do
          name = 'WIND24'
          !call update_accum_field_gpu  (NAME, rbufslp, nstep)
          !call extract_accum_field_gpu (NAME, atm2lnd_vars%wind24_patch, nstep)
          do p = bounds%begp, bounds%endp
             c = veg_pp%column(p)
             g = veg_pp%gridcell(p)
             rbufslp(p) = atm2lnd_vars%forc_rh_grc(g)
          end do
          name = 'RH24'
          !call update_accum_field_gpu  (name, rbufslp, nstep)
          !call extract_accum_field_gpu (name, atm2lnd_vars%rh24_patch, nstep)
       end if
       end associate
     end subroutine atm2lnd_UpdateAccVars

     !-----------------------------------------------------------------------
     subroutine update_acc_vars_top_afGPU (top_af, bounds, nstep)
       !
       ! USES
       !$acc routine seq
       ! !ARGUMENTS:
       type(topounit_atmospheric_flux)    :: top_af
       type(bounds_type) , intent(in) :: bounds
       integer,   value  , intent(in) :: nstep  !timestep number
       character(len=8) :: name
       !
       ! !LOCAL VARIABLES:
       integer :: g,t,c,p                   ! indices
       integer :: begt, endt
       real(r8) :: rbufslt(bounds%begt:bounds%endt)       ! temporary single level - topounit level
       !---------------------------------------------------------------------
       associate( &
         prec10d  =>  top_af%prec10d , &
         prec60d  =>  top_af%prec60d  &
         )
       begt = bounds%begt; endt = bounds%endt
       ! Allocate needed dynamic memory for single level topounit field

       ! Accumulate and extract total precip
       do t = begt,endt
          rbufslt(t) = top_af%rain(t) + top_af%snow(t)
       end do
       if (use_cn) then
          ! Accumulate and extract PREC60D (accumulates total precipitation as 60-day running mean)
          name = 'PREC60D'
          call update_accum_field_gpu  (8,NAME,begt,endt, rbufslt(begt:endt), nstep)
          call extract_accum_field_gpu (8,NAME,begt,endt, prec60d(begt:endt), nstep)
          ! Accumulate and extract PREC10D (accumulates total precipitation as 10-day running mean)
          name = 'PREC10D'
          call update_accum_field_gpu  (7,name,begt,endt, rbufslt(begt:endt), nstep)
          call extract_accum_field_gpu (7,name,begt,endt, prec10d(begt:endt), nstep)
       end if
       if (use_fates) then
         name = 'PREC24H'
          !call update_accum_field_gpu  (name, rbufslt, nstep)
          !call extract_accum_field_gpu (name, top_af%prec24h, nstep)
       end if
       end associate
     end subroutine update_acc_vars_top_afGPU

     !----------------------  -------------------------------------------------
  subroutine update_acc_vars_veg_es_GPU(veg_es, bounds, dtime, nstep, &
            end_cd, month,day, secs)
    !$acc routine seq
    ! USES
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    type(vegetation_energy_state)    :: veg_es
    type(bounds_type)      , intent(in) :: bounds
    integer, intent(in) :: dtime                     ! timestep size [seconds]
    integer, intent(in) :: nstep                     ! timestep number
    logical, intent(in) :: end_cd                    ! temporary for is_end_curr_day() value
    integer, intent(in) :: month                     ! month (1, ..., 12) for nstep
    integer, intent(in) :: day                       ! day of month (1, ..., 31) for nstep
    integer, intent(in) :: secs                      ! seconds into current date for nstep
    !
    ! !LOCAL VARIABLES:
    integer :: m,g,l,c,p                 ! indices
    integer :: ier                       ! error status
    character(len=8) :: name
    integer :: begp, endp
    real(r8) :: rbufslp(bounds%begp:bounds%endp)      ! temporary single level - pft level
    !---------------------------------------------------------------------
    associate( &
      t_ref2m  =>  veg_es%t_ref2m , &
      t_a10    =>  veg_es%t_a10   , &
      t_a10min =>  veg_es%t_a10min   , &
      t_a5min  =>  veg_es%t_a5min   , &
      t_veg24  =>  veg_es%t_veg24 , &
      t_veg240 =>  veg_es%t_veg240 , &
      t_ref2m_u => veg_es%t_ref2m_u , &
      t_ref2m_r => veg_es%t_ref2m_r , &
      gdd0      =>  veg_es%gdd0   , &
      gdd8      =>  veg_es%gdd8   , &
      gdd10     =>  veg_es%gdd10   &
      )
    begp = bounds%begp; endp = bounds%endp
    ! Allocate needed dynamic memory for single level pft field
    ! fill the temporary variable
    do p = begp,endp
       rbufslp(p) = veg_es%t_veg(p)
    end do
    name = 'T10'
    call update_accum_field_gpu  (9,name,begp,endp, t_ref2m(begp:endp), nstep)
    call extract_accum_field_gpu (9,name,begp,endp, t_a10(begp:endp)  , nstep)
    name = 'T_VEG24'
    call update_accum_field_gpu  (13,NAME,begp,endp, rbufslp(begp:endp)       , nstep)
    call extract_accum_field_gpu (13,NAME,begp,endp, t_veg24(begp:endp)  , nstep)
    name = 'T_VEG240'
    call update_accum_field_gpu  (14,name,begp,endp, rbufslp       , nstep)
    call extract_accum_field_gpu (14,name,begp,endp, t_veg240(begp:endp) , nstep)
    ! Accumulate and extract TREFAV - hourly average 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called
    name = 'TREFAV'
    call update_accum_field_gpu  (10,name,begp,endp, t_ref2m(begp:endp), nstep)
    call extract_accum_field_gpu (10,name,begp,endp, rbufslp(begp:endp), nstep)
    do p = begp,endp
       if (rbufslp(p) /= spval) then
          veg_es%t_ref2m_max_inst(p) = max(rbufslp(p), veg_es%t_ref2m_max_inst(p))
          veg_es%t_ref2m_min_inst(p) = min(rbufslp(p), veg_es%t_ref2m_min_inst(p))
       endif
       if (end_cd) then
          veg_es%t_ref2m_max(p) = veg_es%t_ref2m_max_inst(p)
          veg_es%t_ref2m_min(p) = veg_es%t_ref2m_min_inst(p)
          veg_es%t_ref2m_max_inst(p) = -spval
          veg_es%t_ref2m_min_inst(p) =  spval
       else if (secs == int(dtime)) then
          veg_es%t_ref2m_max(p) = spval
          veg_es%t_ref2m_min(p) = spval
       endif
    end do
    ! Accumulate and extract TREFAV_U - hourly average urban 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called
    name = 'TREFAV_U'
    call update_accum_field_gpu  (11,NAME,begp,endp, t_ref2m_u(begp:endp), nstep)
    call extract_accum_field_gpu (11,NAME,begp,endp, rbufslp(begp:endp), nstep)
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (rbufslp(p) /= spval) then
          veg_es%t_ref2m_max_inst_u(p) = max(rbufslp(p), veg_es%t_ref2m_max_inst_u(p))
          veg_es%t_ref2m_min_inst_u(p) = min(rbufslp(p), veg_es%t_ref2m_min_inst_u(p))
       endif
       if (end_cd) then
         if (lun_pp%urbpoi(l)) then
          veg_es%t_ref2m_max_u(p) = veg_es%t_ref2m_max_inst_u(p)
          veg_es%t_ref2m_min_u(p) = veg_es%t_ref2m_min_inst_u(p)
          veg_es%t_ref2m_max_inst_u(p) = -spval
          veg_es%t_ref2m_min_inst_u(p) =  spval
         end if
       else if (secs == int(dtime)) then
          veg_es%t_ref2m_max_u(p) = spval
          veg_es%t_ref2m_min_u(p) = spval
       endif
    end do
    ! Accumulate and extract TREFAV_R - hourly average rural 2m air temperature
    ! Used to compute maximum and minimum of hourly averaged 2m reference
    ! temperature over a day. Note that "spval" is returned by the call to
    ! accext if the time step does not correspond to the end of an
    ! accumulation interval. First, initialize the necessary values for
    ! an initial run at the first time step the accumulator is called
    NAME = 'TREFAV_R'
    call update_accum_field_gpu  (12,name,begp,endp, t_ref2m_r(begp:endp), nstep)
    call extract_accum_field_gpu (12,name,begp,endp, rbufslp(begp:endp), nstep)
    do p = begp,endp
       l = veg_pp%landunit(p)
       if (rbufslp(p) /= spval) then
          veg_es%t_ref2m_max_inst_r(p) = max(rbufslp(p), veg_es%t_ref2m_max_inst_r(p))
          veg_es%t_ref2m_min_inst_r(p) = min(rbufslp(p), veg_es%t_ref2m_min_inst_r(p))
       endif
       if (end_cd) then
         if (.not.(lun_pp%ifspecial(l))) then
          veg_es%t_ref2m_max_r(p) = veg_es%t_ref2m_max_inst_r(p)
          veg_es%t_ref2m_min_r(p) = veg_es%t_ref2m_min_inst_r(p)
          veg_es%t_ref2m_max_inst_r(p) = -spval
          veg_es%t_ref2m_min_inst_r(p) =  spval
         end if
       else if (secs == int(dtime)) then
          veg_es%t_ref2m_max_r(p) = spval
          veg_es%t_ref2m_min_r(p) = spval
       endif
    end do
    if ( crop_prog )then
       ! Accumulate and extract TDM10
       do p = begp,endp
          rbufslp(p) = min(veg_es%t_ref2m_min(p),veg_es%t_ref2m_min_inst(p))
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                                     !'min_inst' not initialized?
       name = 'TDM10'
       call update_accum_field_gpu  (5222,name,begp,endp, rbufslp(begp:endp), nstep)
       call extract_accum_field_gpu (5222,name,begp,endp, t_a10min(begp:endp), nstep)
       ! Accumulate and extract TDM5
       do p = begp,endp
          rbufslp(p) = min(veg_es%t_ref2m_min(p),veg_es%t_ref2m_min_inst(p))
          if (rbufslp(p) > 1.e30_r8) rbufslp(p) = SHR_CONST_TKFRZ !and were 'min'&
       end do                                         !'min_inst' not initialized?
       name = 'TDM5'
       call update_accum_field_gpu  (54554,name,begp,endp, rbufslp(begp:endp), nstep)
       call extract_accum_field_gpu (54554,name,begp,endp, t_a5min(begp:endp), nstep)

       ! Accumulate and extract GDD0
       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(26._r8, veg_es%t_ref2m(p)-SHR_CONST_TKFRZ)) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       name = 'GDD0'
       call update_accum_field_gpu  (454,name,begp,endp, rbufslp(begp:endp), nstep)
       call extract_accum_field_gpu (454,name,begp,endp, gdd0(begp:endp), nstep)
       ! Accumulate and extract GDD8
       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                  veg_es%t_ref2m(p)-(SHR_CONST_TKFRZ + 8._r8))) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       name = 'GDD8'
       call update_accum_field_gpu  (44,name,begp,endp, rbufslp(begp:endp), nstep)
       call extract_accum_field_gpu (44,name,begp,endp, gdd8(begp:endp), nstep)
       ! Accumulate and extract GDD10
       do p = begp,endp
          g = veg_pp%gridcell(p)
          if (month==1 .and. day==1 .and. secs==int(dtime)) then
             rbufslp(p) = accumResetVal ! reset gdd
          else if (( month > 3 .and. month < 10 .and. grc_pp%latdeg(g) >= 0._r8) .or. &
                   ((month > 9 .or.  month < 4) .and. grc_pp%latdeg(g) <  0._r8)     ) then
             rbufslp(p) = max(0._r8, min(30._r8, &
                  veg_es%t_ref2m(p)-(SHR_CONST_TKFRZ + 10._r8))) * dtime/SHR_CONST_CDAY
          else
             rbufslp(p) = 0._r8      ! keeps gdd unchanged at other times (eg, through Dec in NH)
          end if
       end do
       name = 'GDD10'
       call update_accum_field_gpu  (42,name,begp,endp, rbufslp(begp:endp), nstep)
       call extract_accum_field_gpu (42,name,begp,endp, gdd10(begp:endp), nstep)
    end if

    end associate
  end subroutine update_acc_vars_veg_es_GPU

  !-----------------------------------------------------------------------
    subroutine CanopyState_UpdateAccVars (canopystate_vars, bounds, dtime,nstep)
      !
      ! USES
      !$acc routine seq
      ! !ARGUMENTS:
      type(canopystate_type)             :: canopystate_vars
      type(bounds_type)      , intent(in) :: bounds
      integer,value, intent(in) :: dtime                     ! timestep size [seconds]
      integer,value, intent(in) :: nstep                     ! timestep number
      !
      ! !LOCAL VARIABLES:
      integer :: g,p                       ! indices
      integer :: ier                       ! error status
      integer :: begp, endp
      character(len=8) :: name
      real(r8) :: rbufslp(bounds%begp:bounds%endp)      ! temporary single level - pft level
      !---------------------------------------------------------------------
      associate( &
        fsun24_patch  => canopystate_vars%fsun24_patch ,&
        fsun240_patch => canopystate_vars%fsun240_patch  ,&
        elai_p_patch  => canopystate_vars%elai_p_patch  &
        )
      begp = bounds%begp; endp = bounds%endp

      ! Accumulate and extract fsun24 & fsun240
      do p = begp,endp
         rbufslp(p) = canopystate_vars%fsun_patch(p)
      end do
      name = 'FSUN24'
      call update_accum_field_gpu  (15,NAME,begp,endp, rbufslp(begp:endp) , nstep)
      call extract_accum_field_gpu (15,NAME,begp,endp, fsun24_patch(begp:endp), nstep)
      name = 'FSUN240'
      call update_accum_field_gpu  (16,name,begp,endp, rbufslp(begp:endp)         , nstep)
      call extract_accum_field_gpu (16,name,begp,endp, fsun240_patch(begp:endp)   , nstep)

      ! Accumulate and extract elai_patch
      do p = begp,endp
         rbufslp(p) = canopystate_vars%elai_patch(p)
      end do
      name = 'LAIP'
      call update_accum_field_gpu  (17,name,begp,endp, rbufslp (begp:endp)  , nstep)
      call extract_accum_field_gpu (17,name,begp,endp, elai_p_patch(begp:endp)  , nstep)
      end associate
    end subroutine CanopyState_UpdateAccVars

    !-----------------------------------------------------------------------
    subroutine crop_vars_UpdateAccVars(crop_vars, bounds, dtime, nstep)
      !
      ! !DESCRIPTION:
      ! Update accumulated variables. Should be called every time step.
      !
      ! Should only be called if crop_prog is true.
      !$acc routine seq
      ! !USES:
      use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
      use pftvarcon        , only : nwcereal, nwcerealirrig, mxtmp, baset
      use VegetationType   , only : veg_pp
      use ColumnType       , only : col_pp
      use VegetationDataType,only : veg_es
      use ColumnDataType   , only : col_es
      !
      ! !ARGUMENTS:
      type(crop_type)       , intent(inout) :: crop_vars
      type(bounds_type)      , intent(in)    :: bounds
      integer,value, intent(in) :: dtime ! timestep size [seconds]
      integer,value, intent(in) :: nstep ! timestep number
      !
      ! !LOCAL VARIABLES:
      integer :: p,c   ! indices
      integer :: ivt   ! vegetation type
      integer :: begp, endp
      character(len=8) :: name
      real(r8) :: rbufslp(bounds%begp:bounds%endp)      ! temporary single level - pft level

      !-----------------------------------------------------------------------
      associate( &
        gddplant_patch  =>  crop_vars%gddplant_patch  , &
        gddtsoi_patch   =>  crop_vars%gddtsoi_patch   &
        )
      begp = bounds%begp; endp = bounds%endp
      ! Accumulate and extract GDDPLANT

      do p = begp,endp
         if (crop_vars%croplive_patch(p)) then ! relative to planting date
            ivt = veg_pp%itype(p)
            rbufslp(p) = max(0._r8, min(mxtmp(ivt), &
                 veg_es%t_ref2m(p)-(SHR_CONST_TKFRZ + baset(ivt)))) &
                 * dtime/SHR_CONST_CDAY
            if (ivt == nwcereal .or. ivt == nwcerealirrig) then
               rbufslp(p) = rbufslp(p)*crop_vars%vf_patch(p)
            end if
         else
            rbufslp(p) = accumResetVal
         end if
      end do
      name = 'GDDPLANT'
      call update_accum_field_gpu  (0,name,begp,endp, rbufslp(begp:endp), nstep)
      call extract_accum_field_gpu (0,name,begp,endp, gddplant_patch(begp:endp), nstep)

      ! Accumulate and extract GDDTSOI
      ! In agroibis crop_vars variable is calculated
      ! to 0.05 m, so here we use the top two soil layers
      do p = begp,endp
         if (crop_vars%croplive_patch(p)) then ! relative to planting date
            ivt = veg_pp%itype(p)
            c   = veg_pp%column(p)
            rbufslp(p) = max(0._r8, min(mxtmp(ivt), &
                 ((col_es%t_soisno(c,1)*col_pp%dz(c,1) + &
                 col_es%t_soisno(c,2)*col_pp%dz(c,2))/(col_pp%dz(c,1)+col_pp%dz(c,2))) - &
                 (SHR_CONST_TKFRZ + baset(ivt)))) * dtime/SHR_CONST_CDAY
            if (ivt == nwcereal .or. ivt == nwcerealirrig) then
               rbufslp(p) = rbufslp(p)*crop_vars%vf_patch(p)
            end if
         else
            rbufslp(p) = accumResetVal
         end if
      end do
      name = 'GDDTSOI'
      call update_accum_field_gpu  (0,name,begp,endp, rbufslp(begp:endp), nstep)
      call extract_accum_field_gpu (0,name,begp,endp, gddtsoi_patch(begp:endp), nstep)

      end associate
    end subroutine crop_vars_UpdateAccVars

    !-----------------------------------------------------------------------
      subroutine CNState_UpdateAccVars (cnstate_vars, bounds, nstep)
        !
        ! USES
        !$acc routine seq
        ! !ARGUMENTS:
        type(cnstate_type)                    :: cnstate_vars
        type(bounds_type)      , intent(in)    :: bounds
        integer, intent(in) :: nstep                     ! timestep number
        !
        ! !LOCAL VARIABLES:
        integer :: m,g,l,c,p                 ! indices
        integer :: begp, endp
        character(len=8) :: name
        real(r8) :: rbufslp(bounds%begp:bounds%endp)      ! temporary single level - pft level
        !---------------------------------------------------------------------
        associate(  &
          cn_scalar_runmean => cnstate_vars%cn_scalar_runmean , &
          cp_scalar_runmean => cnstate_vars%cp_scalar_runmean &
          )
        begp = bounds%begp; endp = bounds%endp

        ! Accumulate and extract
        do p = begp,endp
           rbufslp(p) = cnstate_vars%cn_scalar(p)
        end do
        name = 'nlim_m'
        call update_accum_field_gpu  (18,name,begp,endp, rbufslp(begp:endp) , nstep)
        call extract_accum_field_gpu (18,name,begp,endp, cn_scalar_runmean(begp:endp), nstep)
        do p = begp,endp
           rbufslp(p) = cnstate_vars%cp_scalar(p)
        end do
        name = 'plim_m'
        call update_accum_field_gpu  (19,name,begp,endp, rbufslp(begp:endp)   , nstep)
        call extract_accum_field_gpu (19,name,begp,endp, cp_scalar_runmean(begp:endp)  , nstep)

        end associate
      end subroutine CNState_UpdateAccVars

end module accumulGPUMod
