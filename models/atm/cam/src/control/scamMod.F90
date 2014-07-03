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
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_scam_mod, only: shr_scam_getCloseLatLon
  use constituents, only: pcnst
  use pmgrid,       only: plon,plev,plevp,plat
  use wrap_nf
  use cam_logfile,  only: iulog
!
  implicit none

  private    ! By default all data is public to this module
!
! !PUBLIC INTERFACES:
!
  public scam_clm_default_opts    ! SCAM default run-time options for CLM
  public scam_default_opts        ! SCAM default run-time options 
  public scam_setopts             ! SCAM run-time options 

!
! !PUBLIC MODULE DATA:
!
  real(r8), public ::  pressure_levels(plev)
  real(r8), public ::  scmlat   ! input namelist latitude for scam
  real(r8), public ::  scmlon   ! input namelist longitude for scam


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
        scm_diurnal_avg_out, scm_crm_mode_out, scm_clubb_iop_name_out)
!-----------------------------------------------------------------------
   real(r8), intent(out), optional :: scmlat_out,scmlon_out
   character*(max_path_len), intent(out), optional ::  iopfile_out
   logical, intent(out), optional ::  single_column_out
   logical, intent(out), optional ::  scm_iop_srf_prop_out
   logical, intent(out), optional ::  scm_relaxation_out
   logical, intent(out), optional ::  scm_diurnal_avg_out
   logical, intent(out), optional ::  scm_crm_mode_out
   character(len=*), intent(out), optional ::  scm_clubb_iop_name_out

   if ( present(scmlat_out) )           scmlat_out     = -999._r8
   if ( present(scmlon_out) )           scmlon_out     = -999._r8
   if ( present(iopfile_out) )          iopfile_out    = ''
   if ( present(single_column_out) )    single_column_out  = .false.
   if ( present(scm_iop_srf_prop_out) )scm_iop_srf_prop_out  = .false.
   if ( present(scm_relaxation_out) )   scm_relaxation_out  = .false.
   if ( present(scm_diurnal_avg_out) )  scm_diurnal_avg_out = .false.
   if ( present(scm_crm_mode_out) )     scm_crm_mode_out  = .false.
   if ( present(scm_clubb_iop_name_out) ) scm_clubb_iop_name_out  = ' '

end subroutine scam_default_opts

subroutine scam_setopts( scmlat_in, scmlon_in,iopfile_in,single_column_in, &
                         scm_iop_srf_prop_in, scm_relaxation_in, &
                         scm_diurnal_avg_in,scm_crm_mode_in, scm_clubb_iop_name_in)
!-----------------------------------------------------------------------
  real(r8), intent(in), optional       :: scmlon_in, scmlat_in
  character*(max_path_len), intent(in), optional :: iopfile_in
  logical, intent(in), optional        :: single_column_in
  logical, intent(in), optional        :: scm_iop_srf_prop_in
  logical, intent(in), optional        :: scm_relaxation_in
  logical, intent(in), optional        :: scm_diurnal_avg_in
  logical, intent(in), optional        :: scm_crm_mode_in
  character(len=*), intent(in), optional :: scm_clubb_iop_name_in
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
  
  if (present (scm_diurnal_avg_in)) then
     scm_diurnal_avg=scm_diurnal_avg_in
  endif
  
  if (present (scm_crm_mode_in)) then
     scm_crm_mode=scm_crm_mode_in
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
                 if (ioplon.lt.0) ioplon=ioplon+360._r8
!!$                 if (ioplon-scmlon.gt.5.) then
!!$                    write(iulog,*)'WARNING: SCMLON/SCMLAT specified in namelist is different'
!!$                    write(iulog,*)'from the IOP file lat,lon by more than 5 degrees'
!!$                    write(iulog,*)'Using specified SCMLAT and SCMLON for all boundary data'
!!$                 endif
                 call shr_scam_GetCloseLatLon(ncid,scmlat,scmlon,ioplat,ioplon,latidx,lonidx)
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

!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!

end module scamMod
