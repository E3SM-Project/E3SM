module scamMod
!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: scamMod
! 
! !DESCRIPTION: 
! scam specific routines and data
!
! !USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
  use pmgrid,       only: plon,plev,plevp,plat
  use wrap_nf
  use cam_logfile,  only: iulog
  use time_manager, only: get_nstep,timemgr_time_inc,&
                          get_start_date,get_curr_date,&
                          timemgr_init,get_curr_calday,&
                          is_first_step
  use shr_scam_mod, only: shr_scam_GetCloseLatLon
  use constituents, only: readtrace, cnst_get_ind, pcnst, cnst_name
  use string_utils, only: to_lower
  use cam_abortutils,   only: endrun
  use phys_control, only: phys_getopts
!
  implicit none

  private    ! By default all data is public to this module
!
! !PUBLIC INTERFACES:
!
  public scam_clm_default_opts    ! SCAM default run-time options for CLM
  public scam_default_opts        ! SCAM default run-time options 
  public scam_setopts             ! SCAM run-time options 
  public setiopupdate
  public readiopdata         

!
! !PUBLIC MODULE DATA:
!
  real(r8), public ::  pressure_levels(plev)
  real(r8), public ::  scmlat   ! input namelist latitude for scam
  real(r8), public ::  scmlon   ! input namelist longitude for scam

  real(r8), allocatable, public :: scm_dgnum( : ),scm_std( : ),&
                                   scm_num( :), scm_div(:,:)

  integer, parameter :: num_switches = 20
  integer, parameter :: max_path_len = 128

  logical, public ::  single_column         ! Using IOP file or not
  logical, public ::  use_iop               ! Using IOP file or not
  logical, public ::  use_analysis
  logical, public ::  use_saveinit
  logical, public ::  use_pert_init         ! perturb initial values
  logical, public ::  use_pert_frc          ! perturb forcing 
  logical, public ::  scm_diurnal_avg       ! If using diurnal averaging or not
  logical, public ::  scm_crm_mode          ! column radiation mode
  logical, public ::  use_userdata
  logical, public ::  isrestart             ! If this is a restart step or not
  logical, public ::  switch(num_switches)  ! Logical flag settings from GUI
  logical, public ::  l_uvphys              ! If true, update u/v after TPHYS
  logical, public ::  l_uvadvect            ! If true, T, U & V will be passed to SLT
  logical, public ::  l_conv                ! use flux divergence terms for T and q?     
  logical, public ::  l_divtr               ! use flux divergence terms for constituents?
  logical, public ::  l_diag                ! do we want available diagnostics?

  integer, public ::  error_code            ! Error code from netCDF reads
  integer, public ::  initTimeIdx
  integer, public ::  seedval
  integer :: closelatidx,closelonidx,latid,lonid,levid,timeid
  real(r8):: closelat,closelon

  character*(max_path_len), public ::  modelfile
  character*(max_path_len), public ::  analysisfile
  character*(max_path_len), public ::  sicfile
  character*(max_path_len), public ::  userfile
  character*(max_path_len), public ::  sstfile
  character*(max_path_len), public ::  lsmpftfile
  character*(max_path_len), public ::  pressfile
  character*(max_path_len), public ::  topofile
  character*(max_path_len), public ::  ozonefile
  character*(max_path_len), public ::  iopfile
  character*(max_path_len), public ::  absemsfile
  character*(max_path_len), public ::  aermassfile
  character*(max_path_len), public ::  aeropticsfile
  character*(max_path_len), public ::  timeinvfile
  character*(max_path_len), public ::  lsmsurffile
  character*(max_path_len), public ::  lsminifile

  real(r8), public ::  fixmascam
  real(r8), public ::  betacam
  real(r8), public ::  alphacam(pcnst)
  real(r8), public ::  dqfxcam(plon,plev,pcnst)

  real(r8), public ::      divq3d(plev,pcnst)  ! 3D q advection
  real(r8), public ::      divt3d(plev)        ! 3D T advection
  real(r8), public ::      vertdivq(plev,pcnst)! vertical q advection
  real(r8), public ::      vertdivt(plev)      ! vertical T advection
  real(r8), public ::      ptend               ! surface pressure tendency
  real(r8), public ::      qdiff(plev)         ! model minus observed humidity
  real(r8), public ::      qobs(plev)          ! actual W.V. Mixing ratio
  real(r8), public ::      cldliqobs(plev)     ! actual W.V. Mixing ratio
  real(r8), public ::      cldiceobs(plev)     ! actual W.V. Mixing ratio
  real(r8), public ::      numliqobs(plev)     ! actual 
  real(r8), public ::      numiceobs(plev)     ! actual 
  real(r8), public ::      precobs(1)          ! observed precipitation 
  real(r8), public ::      lhflxobs(1)         ! observed surface latent heat flux 
  real(r8), public ::      shflxobs(1)         ! observed surface sensible heat flux
  real(r8), public ::      q1obs(plev)         ! observed apparent heat source
  real(r8), public ::      q2obs(plev)         ! observed apparent heat sink
  real(r8), public ::      tdiff(plev)         ! model minus observed temp 
  real(r8), public ::      tground(1)          ! ground temperature
  real(r8), public ::      tobs(plev)          ! actual temperature
  real(r8), public ::      psobs
  real(r8), public ::      tsair(1)            ! air temperature at the surface
  real(r8), public ::      udiff(plev)         ! model minus observed uwind
  real(r8), public ::      uobs(plev)          ! actual u wind
  real(r8), public ::      vdiff(plev)         ! model minus observed vwind
  real(r8), public ::      vobs(plev)          ! actual v wind
  real(r8), public ::      cldobs(plev)        ! observed cld
  real(r8), public ::      clwpobs(plev)       ! observed clwp
  real(r8), public ::      aldirobs(1)         ! observed aldir
  real(r8), public ::      aldifobs(1)         ! observed aldif
  real(r8), public ::      asdirobs(1)         ! observed asdir
  real(r8), public ::      asdifobs(1)         ! observed asdif

  real(r8), public ::      wfld(plev)          ! Vertical motion (slt)
  real(r8), public ::      wfldh(plevp)        ! Vertical motion (slt)
  real(r8), public ::      divq(plev,pcnst)    ! Divergence of moisture
  real(r8), public ::      divt(plev)          ! Divergence of temperature
  real(r8), public ::      divu(plev)          ! Horiz Divergence of E/W
  real(r8), public ::      divv(plev)          ! Horiz Divergence of N/S
                                               ! mo_drydep algorithm
					       
  real(r8), public ::  scm_relaxation_low      ! lowest level to apply relaxation
  real(r8), public ::  scm_relaxation_high     ! highest level to apply relaxation					       
					       
  real(r8), public, pointer :: loniop(:)
  real(r8), public, pointer :: latiop(:)
!
  integer, public ::     iopTimeIdx            ! index into iop dataset
  integer, public ::     steplength            ! Length of time-step
  integer, public ::     base_date             ! Date in (yyyymmdd) of start time
  integer, public ::     base_secs             ! Time of day of start time (sec)

  logical*4, public ::  doiopupdate   ! do we need to read next iop timepoint
  logical*4, public ::  have_divq     ! dataset contains divq 
  logical*4, public ::  have_divt     ! dataset contains divt
  logical*4, public ::  have_divq3d   ! dataset contains divq3d 
  logical*4, public ::  have_vertdivt ! dataset contains vertdivt
  logical*4, public ::  have_vertdivq ! dataset contains vertdivq 
  logical*4, public ::  have_divt3d   ! dataset contains divt3d
  logical*4, public ::  have_divu     ! dataset contains divu
  logical*4, public ::  have_divv     ! dataset contains divv 
  logical*4, public ::  have_omega    ! dataset contains omega
  logical*4, public ::  have_phis     ! dataset contains phis
  logical*4, public ::  have_ptend    ! dataset contains ptend
  logical*4, public ::  have_ps       ! dataset contains ps
  logical*4, public ::  have_q        ! dataset contains q
  logical*4, public ::  have_q1       ! dataset contains Q1
  logical*4, public ::  have_q2       ! dataset contains Q2
  logical*4, public ::  have_prec     ! dataset contains prec 
  logical*4, public ::  have_lhflx    ! dataset contains lhflx 
  logical*4, public ::  have_shflx    ! dataset contains shflx
  logical*4, public ::  have_t        ! dataset contains t
  logical*4, public ::  have_tg       ! dataset contains tg
  logical*4, public ::  have_tsair    ! dataset contains tsair
  logical*4, public ::  have_u        ! dataset contains u 
  logical*4, public ::  have_v        ! dataset contains v 
  logical*4, public ::  have_cld      ! dataset contains cld
  logical*4, public ::  have_cldliq   ! dataset contains cldliq
  logical*4, public ::  have_cldice   ! dataset contains cldice
  logical*4, public ::  have_numliq   ! dataset contains numliq
  logical*4, public ::  have_numice   ! dataset contains numice
  logical*4, public ::  have_clwp     ! dataset contains clwp
  logical*4, public ::  have_aldir    ! dataset contains aldir
  logical*4, public ::  have_aldif    ! dataset contains aldif
  logical*4, public ::  have_asdir    ! dataset contains asdir
  logical*4, public ::  have_asdif    ! dataset contains asdif
  logical*4, public ::  scm_iop_srf_prop   ! use the specified surface properties
  logical*4, public ::  scm_relaxation! use relaxation
  logical*4, public ::  scm_observed_aero ! use observed aerosols in SCM file
  logical*4, public ::  swrad_off     ! turn off SW radiation (assume night)
  logical*4, public ::  lwrad_off     ! turn off LW radiation
  logical*4, public ::  precip_off    ! turn off precipitation processes
  logical*4, public ::  use_camiop    ! use cam generated forcing 
  logical*4, public ::  use_3dfrc     ! use 3d forcing

  character(len=200), public ::  scm_clubb_iop_name   ! IOP name for CLUBB

!=======================================================================
  contains
!=======================================================================

!
!-----------------------------------------------------------------------
!


subroutine scam_default_opts( scmlat_out,scmlon_out,iopfile_out, &
	single_column_out,scm_iop_srf_prop_out, scm_relaxation_out, &
	scm_relaxation_low_out, scm_relaxation_high_out, &
        scm_diurnal_avg_out, scm_crm_mode_out, scm_observed_aero_out, &
	swrad_off_out, lwrad_off_out, precip_off_out, scm_clubb_iop_name_out)
!-----------------------------------------------------------------------
   real(r8), intent(out), optional :: scmlat_out,scmlon_out
   character*(max_path_len), intent(out), optional ::  iopfile_out
   logical, intent(out), optional ::  single_column_out
   logical, intent(out), optional ::  scm_iop_srf_prop_out
   logical, intent(out), optional ::  scm_relaxation_out
   logical, intent(out), optional ::  scm_diurnal_avg_out
   logical, intent(out), optional ::  scm_crm_mode_out
   logical, intent(out), optional ::  scm_observed_aero_out
   logical, intent(out), optional ::  swrad_off_out
   logical, intent(out), optional ::  lwrad_off_out
   logical, intent(out), optional ::  precip_off_out
   real(r8), intent(out), optional ::  scm_relaxation_low_out
   real(r8), intent(out), optional ::  scm_relaxation_high_out   
   character(len=*), intent(out), optional ::  scm_clubb_iop_name_out

   if ( present(scmlat_out) )           scmlat_out     = -999._r8
   if ( present(scmlon_out) )           scmlon_out     = -999._r8
   if ( present(iopfile_out) )          iopfile_out    = ''
   if ( present(single_column_out) )    single_column_out  = .false.
   if ( present(scm_iop_srf_prop_out) )scm_iop_srf_prop_out  = .false.
   if ( present(scm_relaxation_out) )   scm_relaxation_out  = .false.
   if ( present(scm_relaxation_low_out) ) scm_relaxation_low_out = 1050.0_r8
   if ( present(scm_relaxation_high_out) ) scm_relaxation_high_out = 0.e3   
   if ( present(scm_diurnal_avg_out) )  scm_diurnal_avg_out = .false.
   if ( present(scm_crm_mode_out) )     scm_crm_mode_out  = .false.
   if ( present(scm_observed_aero_out)) scm_observed_aero_out = .false.
   if ( present(swrad_off_out))         swrad_off_out = .false.
   if ( present(lwrad_off_out))         lwrad_off_out = .false.
   if ( present(precip_off_out))        precip_off_out = .false.
   if ( present(scm_clubb_iop_name_out) ) scm_clubb_iop_name_out  = ' '

end subroutine scam_default_opts

subroutine scam_setopts( scmlat_in, scmlon_in,iopfile_in,single_column_in, &
                         scm_iop_srf_prop_in, scm_relaxation_in, &
			 scm_relaxation_low_in, scm_relaxation_high_in, &
                         scm_diurnal_avg_in, scm_crm_mode_in, scm_observed_aero_in, &
			 swrad_off_in, lwrad_off_in, precip_off_in, scm_clubb_iop_name_in)
!-----------------------------------------------------------------------
  real(r8), intent(in), optional       :: scmlon_in, scmlat_in
  character*(max_path_len), intent(in), optional :: iopfile_in
  logical, intent(in), optional        :: single_column_in
  logical, intent(in), optional        :: scm_iop_srf_prop_in
  logical, intent(in), optional        :: scm_relaxation_in
  logical, intent(in), optional        :: scm_diurnal_avg_in
  logical, intent(in), optional        :: scm_crm_mode_in
  logical, intent(in), optional        :: scm_observed_aero_in
  logical, intent(in), optional        :: swrad_off_in
  logical, intent(in), optional        :: lwrad_off_in
  logical, intent(in), optional        :: precip_off_in
  character(len=*), intent(in), optional :: scm_clubb_iop_name_in
  real(r8), intent(in), optional       :: scm_relaxation_low_in
  real(r8), intent(in), optional       :: scm_relaxation_high_in  
  integer ncid,latdimid,londimid,latsiz,lonsiz,latid,lonid,ret,i
  integer latidx,lonidx
  real(r8) ioplat,ioplon
  
  if (present (single_column_in ) ) then 
     single_column=single_column_in
  endif

  if (present (scm_iop_srf_prop_in)) then
     scm_iop_srf_prop=scm_iop_srf_prop_in
  endif
  
  if (present (scm_relaxation_in)) then
     scm_relaxation=scm_relaxation_in
  endif
  
  if (present (scm_relaxation_low_in)) then
     scm_relaxation_low=scm_relaxation_low_in
  endif  
  
  if (present (scm_relaxation_high_in)) then
     scm_relaxation_high=scm_relaxation_high_in
  endif   
  
  if (present (scm_diurnal_avg_in)) then
     scm_diurnal_avg=scm_diurnal_avg_in
  endif
  
  if (present (scm_crm_mode_in)) then
     scm_crm_mode=scm_crm_mode_in
  endif

  if (present (scm_observed_aero_in)) then
     scm_observed_aero=scm_observed_aero_in
  endif

  if (present (swrad_off_in)) then
     swrad_off=swrad_off_in
  endif

  if (present (lwrad_off_in)) then
     lwrad_off=lwrad_off_in
  endif
  
  if (present (precip_off_in)) then
     precip_off=precip_off_in
  endif

  if (present (scm_clubb_iop_name_in)) then
     scm_clubb_iop_name=scm_clubb_iop_name_in
  endif

  if (present (iopfile_in)) then
     iopfile=trim(iopfile_in)
  endif

  if( single_column) then
     
     if (plon /= 1 .or. plat /=1 ) then 
        call endrun('SCAM_SETOPTS: must compile model for SCAM mode when namelist parameter single_column is .true.')
     endif
     
     if (present (iopfile_in)) then
        iopfile=trim(iopfile_in)
        if (iopfile.ne."") then 
           use_iop = .true.
        else
           call endrun('SCAM_SETOPTS: must specify IOP file for single column mode')
        endif
        call wrap_open (iopfile, NF90_NOWRITE, ncid)
        call wrap_inq_dimid( ncid, 'lon', londimid   )
        call wrap_inq_dimid( ncid, 'lat', latdimid   )
        call wrap_inq_dimlen( ncid, londimid, lonsiz   )
        call wrap_inq_dimlen( ncid, latdimid, latsiz   )
        call wrap_inq_varid( ncid, 'lon', lonid   )
        call wrap_inq_varid( ncid, 'lat', latid   )
        if ( nf90_inquire_attribute( ncid, NF90_GLOBAL, 'CAM_GENERATED_FORCING', attnum=i ).EQ. NF90_NOERR ) then
           use_camiop = .true.
        else
           use_camiop = .false.
        endif

        if (present (scmlat_in) .and. present (scmlon_in) )then
           scmlat=scmlat_in
           scmlon=scmlon_in
           if( scmlat .lt. -90._r8 .or. scmlat .gt. 90._r8) then
              call endrun('SCAM_SETOPTS: SCMLAT must be between -90. and 90. degrees.')
           elseif( scmlon .lt. 0._r8 .or. scmlon .gt. 360._r8) then
              call endrun('SCAM_SETOPTS: SCMLON must be between 0. and 360. degrees.')
           else
              if (latsiz==1 .and. lonsiz==1) then
                 ret = nf90_get_var(ncid, lonid, ioplon)
                 if (ret/=NF90_NOERR) then
                    call endrun('SCAM_SETOPTS: error reading longitude variable from iopfile')
                 end if
                 ret = nf90_get_var(ncid, latid, ioplat)
                 if (ret/=NF90_NOERR) then
                    call endrun('SCAM_SETOPTS: error reading latitude variable from iopfile')
                 end if
!!$                 if (ioplon-scmlon.gt.5.) then
!!$                    write(iulog,*)'WARNING: SCMLON/SCMLAT specified in namelist is different'
!!$                    write(iulog,*)'from the IOP file lat,lon by more than 5 degrees'
!!$                    write(iulog,*)'Using specified SCMLAT and SCMLON for all boundary data'
!!$                 endif
                 call shr_scam_GetCloseLatLon(ncid,scmlat,scmlon,ioplat,ioplon,latidx,lonidx)
                 if (ioplon.lt.0) ioplon=ioplon+360._r8
                 scmlat=ioplat
                 scmlon=ioplon
                 write(iulog,*)'For CAM Generated IOP using closest dataset lat and lon'
              else
                 if (use_camiop) then
                    call shr_scam_GetCloseLatLon(ncid,scmlat,scmlon,ioplat,ioplon,latidx,lonidx)
                    scmlat=ioplat
                    scmlon=ioplon
                    write(iulog,*)'For CAM Generated IOP using closest dataset lat and lon'
                 endif
              endif
           endif
        else   
           call endrun('namelist variables SCMLAT and SCMLON must be specified for single column mode')
        endif
     endif
!!jt fix this for crm
!!jt   if(scm_crm_modes) then
!!jt      iyear_AD_out     = (base_date-mod(base_date,10000))/10000 ! year AD to calculate the orbital parameters for.
!!jt   else
!!jt      iyear_AD_out     = 1950
!!jt   end if

  else
     if (plon ==1 .and. plat ==1) then 
        call endrun('SCAM_SETOPTS: single_column namelist option must be set to true when running in single column mode')
     endif
  endif


end subroutine scam_setopts
!
!-----------------------------------------------------------------------
!

subroutine scam_clm_default_opts( pftfile_out, srffile_out, inifile_out )
!-----------------------------------------------------------------------
   character(len=*), intent(out) :: pftfile_out
   character(len=*), intent(out) :: srffile_out
   character(len=*), intent(out) :: inifile_out

   pftfile_out = lsmpftfile
   inifile_out = lsminifile
   srffile_out = lsmsurffile
end subroutine scam_clm_default_opts

subroutine setiopupdate

!-----------------------------------------------------------------------
!   
! Open and read netCDF file to extract time information
!
!---------------------------Code history--------------------------------
!
! Written by John Truesdale    August, 1996
! 
!-----------------------------------------------------------------------
  implicit none
#if ( defined RS6000 )
  implicit automatic (a-z)
#endif

!------------------------------Locals-----------------------------------

   integer NCID,i
   integer tsec_varID, time_dimID
   integer, allocatable :: tsec(:)
   integer  ntime
   integer bdate, bdate_varID
   integer STATUS
   integer next_date, next_sec, last_date, last_sec
   integer :: ncsec,ncdate                      ! current time of day,date
   integer :: yr, mon, day                      ! year, month, and day component
   integer :: start_ymd,start_tod
   logical :: doiter
   save tsec, ntime, bdate
   save last_date, last_sec
!------------------------------------------------------------------------------

   if ( get_nstep() .eq. 0 ) then
!     
!     Open  IOP dataset
!     
      STATUS = NF90_OPEN( iopfile, NF90_NOWRITE, NCID )
!     
!     Read time (tsec) variable 
!     
      STATUS = NF90_INQ_VARID( NCID, 'tsec', tsec_varID )
      if ( STATUS .NE. NF90_NOERR ) write(iulog,*)'ERROR - setiopupdate.F:', &
         'Cant get variable ID for tsec'

      STATUS = NF90_INQ_VARID( NCID, 'bdate', bdate_varID )
      if ( STATUS .NE. NF90_NOERR ) then
         STATUS = NF90_INQ_VARID( NCID, 'basedate', bdate_varID )
         if ( STATUS .NE. NF90_NOERR )         &
            write(iulog,*)'ERROR - setiopupdate.F:Cant get variable ID for bdate'
      endif

      STATUS = NF90_INQ_DIMID( NCID, 'time', time_dimID )
      if ( STATUS .NE. NF90_NOERR )  then
         STATUS = NF90_INQ_DIMID( NCID, 'tsec', time_dimID )
         if ( STATUS .NE. NF90_NOERR )  then
            write(iulog,* )'ERROR - setiopupdate.F:Could not find variable dim ID for time'
            STATUS = NF90_CLOSE ( NCID )
            return
         end if
      end if

      if ( STATUS .NE. NF90_NOERR )  &
         write(iulog,*)'ERROR - setiopupdate.F:Cant get variable dim ID for time'

      STATUS = NF90_INQUIRE_DIMENSION( NCID, time_dimID, len=ntime )
      if ( STATUS .NE. NF90_NOERR )then
         write(iulog,*)'ERROR - setiopupdate.F:Cant get time dimlen'
      endif

      if (.not.allocated(tsec)) allocate(tsec(ntime))

      STATUS = NF90_GET_VAR( NCID, tsec_varID, tsec )
      if ( STATUS .NE. NF90_NOERR )then
         write(iulog,*)'ERROR - setiopupdate.F:Cant get variable tsec'
      endif
      STATUS = NF90_GET_VAR( NCID, bdate_varID, bdate )
      if ( STATUS .NE. NF90_NOERR )then
         write(iulog,*)'ERROR - setiopupdate.F:Cant get variable bdate'
      endif
!     Close the netCDF file
      STATUS = NF90_CLOSE( NCID )
!     
!     determine the last date in the iop dataset
!     
      call timemgr_time_inc(bdate, 0, last_date, last_sec, inc_s=tsec(ntime))
!     
!     set the iop dataset index
!    
      iopTimeIdx=0
      do i=1,ntime           ! set the first ioptimeidx
         call timemgr_time_inc(bdate, 0, next_date, next_sec, inc_s=tsec(i))
         call get_start_date(yr,mon,day,start_tod)
         start_ymd = yr*10000 + mon*100 + day

         if ( start_ymd .gt. next_date .or. (start_ymd .eq. next_date &
            .and. start_tod .ge. next_sec)) then
            iopTimeIdx = i
         endif
      enddo

      call get_curr_date(yr,mon,day,ncsec)
      ncdate=yr*10000 + mon*100 + day

      if (iopTimeIdx == 0.or.iopTimeIdx .ge. ntime) then
         call timemgr_time_inc(bdate, 0, next_date, next_sec, inc_s=tsec(1))
         write(iulog,*) 'Error::setiopupdate: Current model time does not fall within IOP period'
         write(iulog,*) ' Current CAM Date is ',ncdate,' and ',ncsec,' seconds'
         write(iulog,*) ' IOP start is        ',next_date,' and ',next_sec,'seconds'
         write(iulog,*) ' IOP end is          ',last_date,' and ',last_sec,'seconds'
         call endrun
      endif

      doiopupdate = .true.

!------------------------------------------------------------------------------
!     Check if iop data needs to be updated and set doiopupdate accordingly
!------------------------------------------------------------------------------
   else                      ! endstep > 1

!      call timemgr_time_inc(bdate, 0, next_date, next_sec,
!      inc_s=tsec(iopTimeIdx+1))

! call a second time
!      call timemgr_time_inc(bdate, 0, next_date2, next_sec2,
!      inc_s2=tsec(iopTimeIdx+2))

      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day

      doiopupdate = .false.
      iopTimeIdx = iopTimeIdx
      i=0
      doiter=.true.
      do while(doiter)
        call timemgr_time_inc(bdate, 0, next_date, next_sec,inc_s=tsec(iopTimeIdx+i+1))
        if (ncdate .gt. next_date .or. (ncdate .eq. next_date &
          .and. ncsec .ge. next_sec)) then

          doiopupdate=.true.
          i=i+1
          iopTimeIdx=iopTimeIdx+1
        else
          doiter=.false.
        endif
      enddo

      if (doiopupdate) then

          write(iulog,*) 'iopTimeIdx =', iopTimeIdx
          write(iulog,*) 'nstep = ',get_nstep()
          write(iulog,*) 'ncdate=',ncdate,' ncsec=',ncsec
          write(iulog,*) 'next_date=',next_date,' next_sec=',next_sec
          write(iulog,*)'******* do iop update'
      endif

   endif                     ! if (endstep .eq. 0 )
!
!     make sure we're
!     not going past end of iop data
!
   if ( ncdate .gt. last_date .or. (ncdate .eq. last_date &
      .and. ncsec .gt. last_sec))  then
      if ( .not. use_userdata ) then
         write(iulog,*)'ERROR - setiopupdate.c:Reached the end of the time varient dataset'
         stop
      else
         doiopupdate = .false.
      end if
   endif

#if DEBUG > 1
   write(iulog,*)'iop time index = ' , ioptimeidx
#endif

   return

end subroutine setiopupdate

  subroutine readiopdata(iop_update_surface,hyam,hybm)

!-----------------------------------------------------------------------
!     
!     Open and read netCDF file containing initial IOP  conditions
!     
!---------------------------Code history--------------------------------
!     
!     Written by J.  Truesdale    August, 1996, revised January, 1998
!     
!-----------------------------------------------------------------------
        use comsrf
        use ppgrid, only: begchunk, endchunk
        use getinterpnetcdfdata
        use shr_sys_mod,      only: shr_sys_flush
        use error_messages, only : handle_ncerr
        use netcdf
        use shr_const_mod, only : SHR_CONST_PI
!-----------------------------------------------------------------------
   implicit none
#if ( defined RS6000 )
   implicit automatic ( a-z )
#endif
!------------------------------Locals-----------------------------------
!     
   logical, intent(in) :: iop_update_surface
   integer NCID, status
   integer time_dimID, lev_dimID,lev_varID,mod_dimID,&
           mod_varID,sps_varID,sps_dimID
   integer tsec_varID, bdate_varID,varid
   integer i,j, ie
   integer nlev, nmod, nsps
   integer total_levs

   integer bdate, ntime, thelev
   integer, allocatable :: tsec(:)
   integer k, m
   integer icldliq,icldice
   integer inumliq,inumice

   logical have_srf              ! value at surface is available
   logical fill_ends             ! 
   logical have_cnst(pcnst)
   real(r8), allocatable :: dplevs( : )
   integer, allocatable :: dmods( : )
   real(r8), intent(in) :: hyam(plev), hybm(plev)
   real(r8) dummy
   real(r8) lat,xlat
   real(r8) srf(1)                  ! value at surface
   real(r8) pmid(plev)  ! pressure at model levels (time n)
   real(r8) pint(plevp) ! pressure at model interfaces (n  )
   real(r8) pdel(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) weight
   real(r8) tmpdata(1)
   real(r8) coldata(plev)
   real(r8) ps_surf, thelat, thelon, the_clat
   integer strt4(4),cnt4(4)
   character(len=16) :: lowername
   real(r8), parameter :: rad2deg = 180.0_r8/SHR_CONST_PI

!   type(dyn_export_t), intent(inout) :: dyn_out

   fill_ends= .false.

!   t1f = Timelevel%n0
!     
!     Open IOP dataset
!     
  call handle_ncerr( nf90_open (iopfile, 0, ncid),&
       'readiopdata.F90', __LINE__)

!
!     if the dataset is a CAM generated dataset set use_camiop to true
!       CAM IOP datasets have a global attribute called CAM_GENERATED_IOP      
!
   if ( nf90_inquire_attribute( ncid, NF90_GLOBAL, 'CAM_GENERATED_FORCING',attnum=i ).EQ. NF90_NOERR ) then
      use_camiop = .true.
   else
      use_camiop = .false.
   endif

!=====================================================================
!     
!     Read time variables


   status = nf90_inq_dimid (ncid, 'time', time_dimID )
   if (status /= NF90_NOERR) then
      status = nf90_inq_dimid (ncid, 'tsec', time_dimID )
      if (status /= NF90_NOERR) then
         write(iulog,* )'ERROR - readiopdata.F:Could not find dimension ID for time/tsec'
         status = NF90_CLOSE ( ncid )
         call endrun
      end if
   end if

   call handle_ncerr( nf90_inquire_dimension( ncid, time_dimID, len=ntime ),&
         'readiopdata.F90', __LINE__)

   allocate(tsec(ntime))

   status = nf90_inq_varid (ncid, 'tsec', tsec_varID )
   call handle_ncerr( nf90_get_var (ncid, tsec_varID, tsec),&
           'readiopdata.F90', __LINE__)

   status = nf90_inq_varid (ncid, 'nbdate', bdate_varID )
   if (status /= NF90_NOERR) then
      status = nf90_inq_varid (ncid, 'bdate', bdate_varID )
      if (status /= NF90_NOERR) then
         write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for bdate'
         status = NF90_CLOSE ( ncid )
         call endrun
      end if
   end if
   call handle_ncerr( nf90_get_var (ncid, bdate_varID, bdate),&
        'readiopdata.F90', __LINE__)

!     
!======================================================
!     read level data
!     
   status = NF90_INQ_DIMID( ncid, 'lev', lev_dimID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable dim ID  for lev'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_inquire_dimension( ncid, lev_dimID, len=nlev ),&
         'readiopdata.F90', __LINE__)

   allocate(dplevs(nlev+1))

   status = NF90_INQ_VARID( ncid, 'lev', lev_varID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for lev'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, lev_varID, dplevs(:nlev)),&
                    'readiopdata.F90', __LINE__)

! =====================================================
!     read observed aersol data

 if(scm_observed_aero .and. .not. iop_update_surface) then
   status = NF90_INQ_DIMID( ncid, 'mod', mod_dimID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable dim ID  for lev'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_inquire_dimension( ncid, mod_dimID, len=nmod ),&
         'readiopdata.F90', __LINE__)

   status = NF90_INQ_DIMID( ncid, 'sps', sps_dimID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable dim ID  for sps'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_inquire_dimension( ncid, sps_dimID, len=nsps ),&
         'readiopdata.F90', __LINE__)

   if (.not.allocated(dmods)) then
      allocate(dmods(nmod))
      dmods=-999
   end if
   if (.not.allocated(scm_num)) then
      allocate(scm_num(nmod))
      scm_num= 1.0e30_R8
   end if
   if (.not.allocated(scm_dgnum)) then
      allocate(scm_dgnum(nmod))
      scm_dgnum= 1.0e30_R8
   end if
   if (.not.allocated(scm_std)) then
      allocate(scm_std(nmod))
      scm_std= 1.0e30_R8
   end if
   if (.not.allocated(scm_div)) then
      allocate(scm_div(nmod,nsps))
      scm_div= 1.0e30_R8
   end if

  status = NF90_INQ_VARID( ncid, 'mod', mod_varID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for mode'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, mod_varID, dmods(:nmod)),&
                    'readiopdata.F90', __LINE__)

status = NF90_INQ_VARID( ncid, 'scm_num', mod_varID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for mode'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, mod_varID, scm_num(:nmod)),&
                    'readiopdata.F90', __LINE__)

status = NF90_INQ_VARID( ncid, 'scm_diam', mod_varID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for mode'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, mod_varID, scm_dgnum(:nmod)),&
                    'readiopdata.F90', __LINE__)

status = NF90_INQ_VARID( ncid, 'scm_std', mod_varID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for mode'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, mod_varID, scm_std(:nmod)),&
                    'readiopdata.F90', __LINE__)

status = NF90_INQ_VARID( ncid, 'scm_accum_div', sps_varID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for mode'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, sps_varID, scm_div(1,:nsps)),&
                    'readiopdata.F90', __LINE__)


status = NF90_INQ_VARID( ncid, 'scm_aitken_div', sps_varID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for mode'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, sps_varID, scm_div(2,:nsps)),&
                    'readiopdata.F90', __LINE__)

status = NF90_INQ_VARID( ncid, 'scm_coarse_div', sps_varID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for mode'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, sps_varID, scm_div(3,:nsps)),&
                    'readiopdata.F90', __LINE__)

endif !scm_observed_aero 
!======================================================================
!
!CAM generated forcing already has pressure on millibars
!
   if (.not. use_camiop) then
!
!     convert pressure to millibars ( lev is expressed in pascals in iop
!     datasets )
!
      do i=1,nlev
         dplevs( i ) = dplevs( i )/100._r8
      end do
   endif

   call shr_scam_GetCloseLatLon(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)

   lonid = 0
   latid = 0
   levid = 0
   timeid = 0

   call wrap_inq_dimid(ncid, 'lat', latid)
   call wrap_inq_dimid(ncid, 'lon', lonid)
   call wrap_inq_dimid(ncid, 'lev', levid)
   call wrap_inq_dimid(ncid, 'time', timeid)

   strt4(1) = closelonidx
   strt4(2) = closelatidx
   strt4(3) = iopTimeIdx
   strt4(4) = 1
   cnt4(1)  = 1
   cnt4(2)  = 1
   cnt4(3)  = 1
   cnt4(4)  = 1

   if (.not. iop_update_surface) then

     status = nf90_inq_varid( ncid, 'Ps', varid   )
     if ( status .ne. nf90_noerr ) then
       have_ps = .false.
       write(iulog,*)'Could not find variable Ps'
       if ( .not. use_userdata ) then
         status = NF90_CLOSE( ncid )
         return
       else
         if ( get_nstep() .eq. 0 ) write(iulog,*) 'Using value from Analysis Dataset'
       endif
     else
       !+ PAB, check the time levels for all variables
       status = nf90_get_var(ncid, varid, psobs, strt4)
       have_ps = .true.
     endif

!  for reproducing CAM output don't do interpolation.
!  the most expedient way of doing this is to set      
!  the dataset pressure levels to the current
!  scam model levels
        
     if ( use_camiop ) then
       do i = 1, plev
         dplevs( i ) = 1000.0_r8 * hyam( i ) + psobs * hybm( i ) / 100.0_r8
       end do
     endif

!     add the surface pressure to the pressure level data, so that
!     surface boundary condition will be set properly,
!     making sure that it is the highest pressure in the array.
!

     total_levs = nlev+1
     dplevs(nlev+1) = psobs/100 ! ps is expressed in pascals
     do i= nlev, 1, -1
       if ( dplevs(i) .GT. psobs/100.0_r8) then
         total_levs = i
         dplevs(i) = psobs/100
       end if
     end do
     if (.not. use_camiop ) then
       nlev = total_levs
     endif
     if ( nlev .eq. 1 ) then
       write(iulog,*) 'Error - Readiopdata.F: Ps too low!'
       return
     endif

!=====================================================================


     status =  nf90_inq_varid( ncid, 'Tsair', varid   )
     if ( status .ne. nf90_noerr ) then
       have_tsair = .false.
     else
       call wrap_get_vara_realx (ncid,varid,strt4,cnt4,tsair)
       have_tsair = .true.
     endif

!
!      read in Tobs  For cam generated iop readin small t to avoid confusion
!      with capital T defined in cam
!

!!!!!!!force fill_end to be .true in getinterpncdata () for temperature !!!!!!!!
     if ( use_camiop ) then
       call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,'t', have_tsair, &
          tsair(1), .true. , scm_crm_mode, &
          dplevs, nlev, psobs, hyam, hybm, tobs, status )
     else
       call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,'T', have_tsair, &
          tsair(1), .true. , scm_crm_mode, &
          dplevs, nlev, psobs, hyam, hybm, tobs, status )
     endif
     if ( status .ne. nf90_noerr ) then
       have_t = .false.
       write(iulog,*)'Could not find variable T'
       if ( .not. use_userdata ) then
         status = NF90_CLOSE( ncid )
         return
       else
         write(iulog,*) 'Using value from Analysis Dataset'
       endif
!     
!     set T3 to Tobs on first time step
!     
     else
       have_t = .true.
     endif

     status = nf90_inq_varid( ncid, 'qsrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       have_srf = .true.
     endif
!
!!!!!!!force fill_end to be .true in getinterpncdata () for humidity!!!!!!!!
     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'q', have_srf, &
       srf(1), .true., scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, qobs, status )
     if ( status .ne. nf90_noerr ) then
       have_q = .false.
       write(iulog,*)'Could not find variable q'
       if ( .not. use_userdata ) then
         status = nf90_close( ncid )
         return
       else
         write(iulog,*) 'Using values from Analysis Dataset'
       endif
     else
       have_q=.true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'cld', .false., &
       dummy, fill_ends, scm_crm_mode, dplevs, nlev,psobs, hyam, hybm, cldobs, status )
     if ( status .ne. nf90_noerr ) then
       have_cld = .false.
     else
       have_cld = .true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'clwp', .false., &
       dummy, fill_ends, scm_crm_mode, dplevs, nlev,psobs, hyam, hybm, clwpobs, status )
     if ( status .ne. nf90_noerr ) then
       have_clwp = .false.
     else
       have_clwp = .true.
     endif

!
!       read divq (horizontal advection)
!      
     status = nf90_inq_varid( ncid, 'divqsrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       have_srf = .true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
        'divq', have_srf, srf(1), fill_ends, scm_crm_mode, &
        dplevs, nlev,psobs, hyam, hybm, divq(:,1), status )
     if ( status .ne. nf90_noerr ) then
       have_divq = .false.
     else
       have_divq = .true.
     endif

!
!     read vertdivq if available
!
     status = nf90_inq_varid( ncid, 'vertdivqsrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       have_srf = .true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'vertdivq', &
        have_srf, srf(1), fill_ends, scm_crm_mode, &
        dplevs, nlev,psobs, hyam, hybm, vertdivq(:,1), status )
     if ( status .ne. nf90_noerr ) then
       have_vertdivq = .false.
     else
       have_vertdivq = .true.
     endif

     status = nf90_inq_varid( ncid, 'vertdivqsrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       have_srf = .true.
     endif

!
!   add calls to get dynamics tendencies for all prognostic consts
!
     do m = 1, pcnst

       call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_dten', &
         have_srf, srf(1), fill_ends, scm_crm_mode, &
         dplevs, nlev,psobs, hyam, hybm, divq3d(:,m), status )
       if ( status .ne. nf90_noerr ) then
         have_cnst(m) = .false.
         divq3d(1:,m)=0._r8
       else
         have_cnst(m) = .true.
       endif

       call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_dqfx', &
         have_srf, srf(1), fill_ends, scm_crm_mode, &
         dplevs, nlev,psobs, hyam, hybm, coldata, status )
       if ( STATUS .NE. NF90_NOERR ) then
         dqfxcam=0._r8
       else
         dqfxcam(1,:,m)=coldata(:)
       endif

       call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_alph', &
          have_srf, srf(1), fill_ends, scm_crm_mode, &
         dplevs, nlev,psobs, hyam, hybm, tmpdata, status )
       if ( status .ne. nf90_noerr ) then
!         have_cnst(m) = .false.
         alphacam(m)=0._r8
       else
         alphacam(m)=tmpdata(1)
!        have_cnst(m) = .true.
       endif

     end do


     call cnst_get_ind('NUMLIQ', inumliq, abort=.false.)
     if ( inumliq > 0 ) then
       have_srf = .false.
       call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'NUMLIQ', &
         have_srf, srf(1), fill_ends, scm_crm_mode, &
         dplevs, nlev,psobs, hyam, hybm, numliqobs, status )
       if ( status .ne. nf90_noerr ) then
         have_numliq = .false.
       else
         have_numliq = .true.
       endif
     end if

     call cnst_get_ind('CLDLIQ', icldliq)

     have_srf = .false.
     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'CLDLIQ', &
       have_srf, srf(1), fill_ends, scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, cldliqobs, status )
     if ( status .ne. nf90_noerr ) then
       have_cldliq = .false.
     else
       have_cldliq = .true.
     endif

     call cnst_get_ind('CLDICE', icldice)

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'CLDICE', &
       have_srf, srf(1), fill_ends, scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, cldiceobs, status )
     if ( status .ne. nf90_noerr ) then
       have_cldice = .false.
     else
       have_cldice = .true.
     endif

     call cnst_get_ind('NUMICE', inumice, abort=.false.)
     if ( inumice > 0 ) then
        have_srf = .false.

        call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'NUMICE', &
          have_srf, srf(1), fill_ends, scm_crm_mode, &
          dplevs, nlev,psobs, hyam, hybm, numiceobs, status )
       if ( status .ne. nf90_noerr ) then
         have_numice = .false.
       else
         have_numice = .true.
       endif
     end if

!
!       read divu (optional field)
!      
     status = nf90_inq_varid( ncid, 'divusrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       have_srf = .true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divu', &
       have_srf, srf(1), fill_ends, scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, divu, status )
     if ( status .ne. nf90_noerr ) then
       have_divu = .false.
     else
       have_divu = .true.
     endif
!
!       read divv (optional field)
!      
     status = nf90_inq_varid( ncid, 'divvsrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       have_srf = .true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divv', &
       have_srf, srf(1), fill_ends, scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, divv, status )
     if ( status .ne. nf90_noerr ) then
       have_divv = .false.
     else
       have_divv = .true.
     endif
!
!       read divt (optional field)
!      
     status = nf90_inq_varid( ncid, 'divtsrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       have_srf = .true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
       'divT', have_srf, srf(1), fill_ends, scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, divt, status )
     if ( status .ne. nf90_noerr ) then
       have_divt = .false.
     else
       have_divt = .true.
     endif

!
!     read vertdivt if available
!
     status = nf90_inq_varid( ncid, 'vertdivTsrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       have_srf = .true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'vertdivT', &
       have_srf, srf(1), fill_ends, scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, vertdivt, status )
     if ( status .ne. nf90_noerr ) then
       have_vertdivt = .false.
     else
       have_vertdivt = .true.
     endif
!
!       read divt3d (combined vertical/horizontal advection)
!      (optional field)

     status = nf90_inq_varid( ncid, 'divT3dsrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       have_srf = .true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divT3d', &
       have_srf, srf(1), fill_ends, scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, divt3d, status )
     if ( status .ne. nf90_noerr ) then
       have_divt3d = .false.
     else
       have_divt3d = .true.
     endif

     status = nf90_inq_varid( ncid, 'Ptend', varid   )
     if ( status .ne. nf90_noerr ) then
       have_ptend = .false.
       write(iulog,*)'Could not find variable Ptend. Setting to zero'
       ptend = 0.0_r8
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       have_ptend = .true.
       ptend= srf(1)
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
       'omega', .true., ptend, fill_ends, scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, wfld, status )
     if ( status .ne. nf90_noerr ) then
       have_omega = .false.
       write(iulog,*)'Could not find variable omega'
       if ( .not. use_userdata ) then
         status = nf90_close( ncid )
         return
       else
         write(iulog,*) 'Using value from Analysis Dataset'
       endif
     else
       have_omega = .true.
     endif

     call plevs0(1    ,plon   ,plev    ,psobs   ,pint,pmid ,pdel)
     call shr_sys_flush( iulog )
!
! Build interface vector for the specified omega profile
! (weighted average in pressure of specified level values)
!
     wfldh(1) = 0.0_r8

     do k=2,plev
       weight = (pint(k) - pmid(k-1))/(pmid(k) - pmid(k-1))
       wfldh(k) = (1.0_r8 - weight)*wfld(k-1) + weight*wfld(k)
     end do

     wfldh(plevp) = 0.0_r8


     status = nf90_inq_varid( ncid, 'usrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       call wrap_get_vara_realx (ncid,varid,strt4,cnt4,srf)
       have_srf = .true.
     endif

!!!!!!!force fill_end to be .true in getinterpncdata () for u-wind !!!!!!!!
     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
       'u', have_srf, srf(1), .true. , scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, uobs, status )
     if ( status .ne. nf90_noerr ) then
       have_u = .false.
     else
       have_u = .true.
     endif

     status = nf90_inq_varid( ncid, 'vsrf', varid   )
     if ( status .ne. nf90_noerr ) then
       have_srf = .false.
     else
       call wrap_get_vara_realx (ncid,varid,strt4,cnt4,srf)
       have_srf = .true.
     endif

!!!!!!!force fill_end to be .true in getinterpncdata () for v-wind !!!!!!!!
     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
       'v', have_srf, srf(1), .true. , scm_crm_mode, &
       dplevs, nlev,psobs, hyam, hybm, vobs, status )
     if ( status .ne. nf90_noerr ) then
       have_v = .false.
     else
       have_v = .true.
     endif
     call shr_sys_flush( iulog )

     status = nf90_inq_varid( ncid, 'Prec', varid   )
     if ( status .ne. nf90_noerr ) then
       have_prec = .false.
     else
       call wrap_get_vara_realx (ncid,varid,strt4,cnt4,precobs)
       have_prec = .true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'Q1', &
       .false., dummy, fill_ends, scm_crm_mode, & ! datasets don't contain Q1 at surface
       dplevs, nlev,psobs, hyam, hybm, q1obs, status )
     if ( status .ne. nf90_noerr ) then
       have_q1 = .false.
     else
       have_q1 = .true.
     endif

     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'Q2', &
        .false., dummy, fill_ends, scm_crm_mode, & ! datasets don't contain Q2 at surface
        dplevs, nlev,psobs, hyam, hybm, q1obs, status )
     if ( status .ne. nf90_noerr ) then
       have_q2 = .false.
     else
       have_q2 = .true.
     endif

!  Test for BOTH 'lhflx' and 'lh' without overwriting 'have_lhflx'.  
!  Analagous changes made for the surface heat flux

     call shr_sys_flush( iulog )

!
!     fill in 3d forcing variables if we have both horizontal
!     and vertical components, but not the 3d
!
     if ( .not. have_cnst(1) .and. have_divq .and. have_vertdivq ) then
       do k=1,plev
         do m=1,pcnst
           divq3d(k,m) = divq(k,m) + vertdivq(k,m)
         enddo
       enddo
       have_divq3d = .true.
     endif

     if ( .not. have_divt3d .and. have_divt .and. have_vertdivt ) then
       write(iulog,*) 'Don''t have divt3d - using divt and vertdivt'
       do k=1,plev
         divt3d(k) = divt(k) + vertdivt(k)
        enddo
        have_divt3d = .true.
     endif
!
!     make sure that use_3dfrc flag is set to true if we only have
!     3d forcing available
!
!   if ( .not. have_divt .or. .not. have_divq ) then
     if (have_divt3d .or. have_divq3d) then
       use_3dfrc = .true.
     endif
     call shr_sys_flush( iulog )

!   status =  nf90_inq_varid( ncid, 'CLAT', varid   )
!   if ( status .eq. nf90_noerr ) then
!      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,clat)
!      clat_p(1)=clat(1)
!      latdeg(1) = clat(1)*45._r8/atan(1._r8)
!   endif

     status =  nf90_inq_varid( ncid, 'beta', varid   )
     if ( status .ne. nf90_noerr ) then
       betacam = 0._r8
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       betacam=srf(1)
     endif

     status =  nf90_inq_varid( ncid, 'fixmas', varid   )
     if ( status .ne. nf90_noerr ) then
       fixmascam=1.0_r8
     else
       status = nf90_get_var(ncid, varid, srf(1), strt4)
       fixmascam=srf(1)
     endif
   
   else ! if read in surface information
   
     status = nf90_inq_varid( ncid, 'Tg', varid   )
     if (status .ne. nf90_noerr) then
       write(iulog,*)'Could not find variable Tg'
       if ( have_tsair ) then
         write(iulog,*) 'Using Tsair'
         tground = tsair     ! use surface value from T field
       else
         write(iulog,*) 'Using T at lowest level'
         tground = tobs(plev) 
       endif
     else
       call wrap_get_vara_realx (ncid,varid,strt4,cnt4,tground)
       have_tg = .true.
     endif 
   
     status = nf90_inq_varid( ncid, 'lhflx', varid   )
     if ( status .ne. nf90_noerr ) then
       status = nf90_inq_varid( ncid, 'lh', varid   )
       if ( status .ne. nf90_noerr ) then
         have_lhflx = .false.
       else
         call wrap_get_vara_realx (ncid,varid,strt4,cnt4,lhflxobs)
         have_lhflx = .true.
       endif
     else
       call wrap_get_vara_realx (ncid,varid,strt4,cnt4,lhflxobs)
       have_lhflx = .true.
     endif

     status = nf90_inq_varid( ncid, 'shflx', varid   )
     if ( status .ne. nf90_noerr ) then
       status = nf90_inq_varid( ncid, 'sh', varid   )
       if ( status .ne. nf90_noerr ) then
         have_shflx = .false.
       else
         call wrap_get_vara_realx (ncid,varid,strt4,cnt4,shflxobs)
         have_shflx = .true.
       endif
     else
       call wrap_get_vara_realx (ncid,varid,strt4,cnt4,shflxobs)
       have_shflx = .true.
     endif
     
   endif

   call shr_sys_flush( iulog )

   status = nf90_close( ncid )
   call shr_sys_flush( iulog )
   status = nf90_close( ncid )
   call shr_sys_flush( iulog )

   deallocate(dplevs,tsec)

   return
end subroutine readiopdata


end module scamMod
