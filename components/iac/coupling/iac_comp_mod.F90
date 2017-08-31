#define DEBUG
Module iac_comp_mod
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: iac_comp_mod
!
!  Interface of the integrated assessment component in CCSM
!
! !DESCRIPTION:
!
! !USES:
  use gcam_comp_mod
  use glm_comp_mod
  use iac2gcam_mod
  use gcam2glm_mod
  use gcam2emisfile_mod
  use glm2iac_mod
  use iac_fields_mod
  use shr_cal_mod
  use shr_file_mod
  use shr_sys_mod
  use shr_kind_mod, only : r8 => shr_kind_r8,r4 => shr_kind_r4
  use netcdf

  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: iac_init_mod               ! clm initialization
  public :: iac_run_mod                ! clm run phase
  public :: iac_final_mod              ! clm finalization/cleanup

! !PUBLIC DATA MEMBERS: None

! !REVISION HISTORY:
! Author: T Craig


! !PRIVATE DATA MEMBERS:

  real*8, pointer :: gcami(:,:)
  real*8, pointer :: gcamiold(:,:)
  real*8, pointer :: gcamiinterp(:,:,:)
  real*8, pointer :: gcamiptr(:,:)
  real*8, pointer :: gcamo(:,:)
  real*8, pointer :: gcamoemis(:,:)
  real*8, pointer :: co2gcam2005base(:,:,:)
  real*8, pointer :: glmi(:,:)
  real*8, pointer :: glmi_wh(:)
  real*8, pointer :: glmo(:,:)
  integer,save :: iulog
  integer,save :: iacymd_fudge = -1

  character(len=512) :: clmC_bfn_dir = 'unknown'
  character(len=512) :: clm2gcam_mapfile = 'unknown'
  character(len=512) :: iac_base_clmfile = 'unknown'
  character(len=512) :: gcam2glm_basecrop = 'unknown'
  character(len=512) :: gcam2glm_basepast = 'unknown'
  character(len=512) :: gcam2glm_baseothr = 'unknown'
  character(len=512) :: gcam2glm_aezmap = 'unknown'
  character(len=512) :: gcam2glm_basebiomass = 'unknown'
  character(len=512) :: gcam2emisfile_co2base2000 = 'unknown'
  character(len=512) :: gcam2emisfile_co2shipbase2000 = 'unknown'
  character(len=512) :: gcam2emisfile_grid720x360 = 'unknown'
  character(len=512) :: gcam2emisfile_grid288x192 = 'unknown'
  character(len=512) :: gcam2emisfile_lut720x360map = 'unknown'
  character(len=512) :: gcam2emisfile_downscaleinfo = 'unknown'
  character(len=512) :: gcam2emisfile_rcp45allsteps = 'unknown'
  character(len=*),parameter :: iac_restfile = 'iac_restart.'
  character(len=*),parameter :: iac_rpointer = 'rpointer.iac'
  character(256) :: rstfilename
  
  logical :: fast_oneway_iac_coupling = .false.
  logical :: clm_iac_carbon_scaling = .true.
  logical :: sneakermode = .false.
  logical :: co2flux_coupling = .false.
  logical :: npp_hr_on = .false.
  logical :: initial_run = .true.
  integer :: clm_nx, clm_ny = -1
  integer :: long_gcam_timestep

  namelist /iacnml/   &
       fast_oneway_iac_coupling,sneakermode,npp_hr_on,initial_run, &
       clm2gcam_mapfile, iac_base_clmfile, clm_nx, clm_ny, clmC_bfn_dir, &
       gcam2glm_basecrop,gcam2glm_basepast,gcam2glm_baseothr,gcam2glm_aezmap, &
       gcam2glm_basebiomass,clm_iac_carbon_scaling,co2flux_coupling,gcam2emisfile_co2base2000, &
       gcam2emisfile_grid720x360,gcam2emisfile_grid288x192,gcam2emisfile_co2shipbase2000, &
       gcam2emisfile_lut720x360map,gcam2emisfile_downscaleinfo,gcam2emisfile_rcp45allsteps
!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_init_mod

! !INTERFACE:
  subroutine iac_init_mod( EClock, cdata, iaci, iaco)

! !DESCRIPTION:
! Initialize interface for iac

! !USES:
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real*8, pointer :: iaci(:,:)
    real*8, pointer :: iaco(:,:)

! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='(iac_init_mod)'
    integer :: ymdtest,iacymd,iactod,run_up_date,ncid
    integer :: gcam_data_size
    integer :: glm_data_size
    integer :: iac_data_size
    integer :: nunit, ier,varid
    integer :: iacymd_hold
    logical :: lexist
    character(len=128) :: casename

! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

  nunit = shr_file_getUnit()
  open(nunit,file="iac_in",status="old",action="read")
  read(nunit, iacnml, iostat=ier)
  if (ier /= 0) then
     write(iulog,*)'error: iacnml namelist input resulted in error code ',ier
     call shr_sys_abort(subname//' ERROR: iacnml error')
  endif
  close(nunit)
  call shr_file_freeUnit(nunit)

  iulog = cdata%i(iac_cdatai_logunit)
  cdata%l(iac_cdatal_fastiac) = fast_oneway_iac_coupling
  cdata%l(iac_cdatal_sneakermode) = sneakermode
  cdata%l(iac_cdatal_co2flux_coupling) = co2flux_coupling
  cdata%l(iac_cdatal_npphr)   = npp_hr_on
  cdata%l(iac_cdatal_nocarbonscale)   = clm_iac_carbon_scaling
  cdata%l(iac_cdatal_initrun) = initial_run
  cdata%c(iac_cdatac_clmcbfndir) = trim(clmC_bfn_dir)
  cdata%c(iac_cdatac_gcam2glm_basecrop) = trim(gcam2glm_basecrop)
  cdata%c(iac_cdatac_gcam2glm_basepast) = trim(gcam2glm_basepast)
  cdata%c(iac_cdatac_gcam2glm_baseothr) = trim(gcam2glm_baseothr)
  cdata%c(iac_cdatac_gcam2glm_aezmap)   = trim(gcam2glm_aezmap)
  cdata%c(iac_cdatac_gcam2glm_basebiomass)  = trim(gcam2glm_basebiomass)
  cdata%c(iac_cdatac_clm2gcam) = trim(clm2gcam_mapfile)
  cdata%c(iac_cdatac_ibclmfile) = trim(iac_base_clmfile)
  cdata%c(iac_cdatac_gcam2emisfile_co2base2000)=trim(gcam2emisfile_co2base2000)
  cdata%c(iac_cdatac_gcam2emisfile_co2shipbase2000)=trim(gcam2emisfile_co2shipbase2000)
  cdata%c(iac_cdatac_gcam2emisfile_grid720x360)=trim(gcam2emisfile_grid720x360)
  cdata%c(iac_cdatac_gcam2emisfile_grid288x192)=trim(gcam2emisfile_grid288x192)
  cdata%c(iac_cdatac_gcam2emisfile_lut720x360map)=trim(gcam2emisfile_lut720x360map)
  cdata%c(iac_cdatac_gcam2emisfile_downscaleinfo)=trim(gcam2emisfile_downscaleinfo)
  cdata%c(iac_cdatac_gcam2emisfile_rcp45allsteps)=trim(gcam2emisfile_rcp45allsteps)

  cdata%i(iac_cdatai_iac_nx)  = clm_nx
  cdata%i(iac_cdatai_iac_ny)  = clm_ny
  cdata%i(iac_cdatai_iac_size) = clm_nx * clm_ny
  casename = trim(cdata%c(iac_cdatac_casename))

  write(iulog,*) subname,' iacnml settings:'
  write(iulog,*) subname,' clm_nx    = ',cdata%i(iac_cdatai_iac_nx)
  write(iulog,*) subname,' clm_ny    = ',cdata%i(iac_cdatai_iac_ny)
  write(iulog,*) subname,' fastiac   = ',cdata%l(iac_cdatal_fastiac)
  write(iulog,*) subname,' sneakermode= ',cdata%l(iac_cdatal_sneakermode)
  write(iulog,*) subname,' co2flux_coupling= ',cdata%l(iac_cdatal_co2flux_coupling)
  write(iulog,*) subname,' clm->iac carbon scaling = ',cdata%l(iac_cdatal_nocarbonscale)
  write(iulog,*) subname,' npphr     = ',cdata%l(iac_cdatal_npphr)
  write(iulog,*) subname,' initrun   = ',cdata%l(iac_cdatal_initrun)
  write(iulog,*) subname,' clm2gcam  = ',trim(cdata%c(iac_cdatac_clm2gcam))
  write(iulog,*) subname,' clmC_bfn_dir = ',trim(cdata%c(iac_cdatac_clmcbfndir))
  write(iulog,*) subname,' ibclmfile = ',trim(cdata%c(iac_cdatac_ibclmfile))
  write(iulog,*) subname,' gcam2glm_basecrop = ',cdata%c(iac_cdatac_gcam2glm_basecrop)
  write(iulog,*) subname,' gcam2glm_basepast = ',cdata%c(iac_cdatac_gcam2glm_basepast)
  write(iulog,*) subname,' gcam2glm_baseothr = ',cdata%c(iac_cdatac_gcam2glm_baseothr)
  write(iulog,*) subname,' gcam2emisfile_co2base2000 = ',cdata%c(iac_cdatac_gcam2emisfile_co2base2000)
  write(iulog,*) subname,' gcam2emisfile_co2shipbase2000 = ',cdata%c(iac_cdatac_gcam2emisfile_co2shipbase2000)
  write(iulog,*) subname,' gcam2emisfile_grid720x360 = ',cdata%c(iac_cdatac_gcam2emisfile_grid720x360)
  write(iulog,*) subname,' gcam2emisfile_grid288x192 = ',cdata%c(iac_cdatac_gcam2emisfile_grid288x192)
  write(iulog,*) subname,' gcam2emisfile_lut720x360map = ',cdata%c(iac_cdatac_gcam2emisfile_lut720x360map)
  write(iulog,*) subname,' gcam2emisfile_downscaleinfo = ',cdata%c(iac_cdatac_gcam2emisfile_downscaleinfo)
  write(iulog,*) subname,' gcam2emisfile_rcp45allsteps = ',cdata%c(iac_cdatac_gcam2emisfile_rcp45allsteps)
  write(iulog,*) subname,' gcam2glm_aezmap = ',cdata%c(iac_cdatac_gcam2glm_aezmap)
  write(iulog,*) subname,' gcam2glm_basebiomass = ',cdata%c(iac_cdatac_gcam2glm_basebiomass)

  call iac_fields_init

  iac_data_size = cdata%i(iac_cdatai_iac_size)
  allocate(iaci(iac_iaci_nflds,iac_data_size))
  allocate(iaco(iac_iaco_nflds,iac_data_size))
  iaci = iac_spval
  iaco = iac_spval
  cdata%l(iac_cdatal_iac_present) = .true.
  cdata%l(iac_cdatal_iac_prognostic) = .true.

  if (sneakermode) then
     long_gcam_timestep=15
  else
     long_gcam_timestep=5
  end if


  !--- initialize models ---
  call shr_sys_flush(iulog)
  call gcam_init_mod(EClock, cdata, gcami, gcamo, gcamoemis  )

  !--- Initialize gcam2emisfile and run up to current date

  if (co2flux_coupling) call gcam2emisfile_init_mod( EClock, cdata, gcamoemis)

  !
  ! allocate arrays to interpolate the gcami data
  !

  allocate(gcamiinterp(iac_gcami_nflds,cdata%i(iac_cdatai_gcami_size),long_gcam_timestep))
  allocate(gcamiold(iac_gcami_nflds,cdata%i(iac_cdatai_gcami_size)))

  call shr_sys_flush(iulog)
  call glm_init_mod (EClock, cdata, glmi , glmi_wh, glmo )

  !--- initialize couplers ---
  call shr_sys_flush(iulog)
  call iac2gcam_init_mod(EClock, cdata, iaci, gcami)

  call shr_sys_flush(iulog)
  call gcam2glm_init_mod(EClock, cdata, gcamo, glmi, glmi_wh )

  call shr_sys_flush(iulog)
  call glm2iac_init_mod (EClock, cdata, glmo,  iaco)

  gcam_data_size = size(gcamo,2)
  glm_data_size = size(glmi,2)
  write(iulog,*) trim(subname),' case name  = ',trim(casename)
  write(iulog,*) trim(subname),' iac size   = ',iac_data_size
  write(iulog,*) trim(subname),' gcamo size = ',gcam_data_size
  write(iulog,*) trim(subname),' glm size   = ',glm_data_size

  !--- Run GCAM up to current date and initialize gcam2glm with starting data ---

  iacymd = EClock(iac_EClock_ymd)
  iactod = EClock(iac_Eclock_tod)
  iacymd_hold=iacymd

  if (initial_run) then
     gcamiold=iac_spval
     call write_iac_restart(iac_rpointer,iac_restfile,iacymd,'gcami',gcamiold)
  else
     ! read restart and set gcamiold
     inquire(file=trim(iac_rpointer),exist=lexist)
     if (lexist) then
        
#ifdef DEBUG
        write(iulog,*) subname,' read_restart rpointer ',trim(iac_rpointer)
#endif
        nunit = shr_file_getunit()
        open(nunit,file=trim(iac_rpointer),form='formatted')
        read(nunit,'(a)') rstfilename
        close(nunit)
        call shr_file_freeunit(nunit)
        
#ifdef DEBUG
        write(iulog,*) subname,' read_restart file ',trim(rstfilename)
#endif
        
        inquire(file=trim(rstfilename),exist=lexist)
        if (.not.lexist) then
           write(iulog,*) subname,' ERROR: missing file ',trim(rstfilename)
           call shr_sys_abort(subname//' ERROR: missing file')
        endif
        
        call iac_ncerr(nf90_open(rstfilename,nf90_nowrite,ncid),subname,__LINE__)
        
        ! Need to save off gcami carbon scalar data
        
        call iac_ncerr(nf90_inq_varid(ncid,'gcami',varid),subname,__LINE__)
        call iac_ncerr(nf90_get_var(ncid,varid,gcamiold),subname,__LINE__)
        call iac_ncerr(nf90_inq_varid(ncid,'ymd',varid),subname,__LINE__)
        call iac_ncerr(nf90_get_var(ncid,varid,iacymd_fudge),subname,__LINE__)
        
        if (fast_oneway_iac_coupling) then 
           iacymd=iacymd_fudge
        end if

        call iac_ncerr(nf90_close(ncid),subname//':close',__LINE__)

     else
        write(iulog,*) subname,' read_restart rpointer NOT found ',trim(iac_rpointer)
        call shr_sys_abort(subname//' ERROR: missing file')
     end if ! rpointer exist
  end if ! initial_run

  run_up_date=iacymd

  if (.not.initial_run) then
     run_up_date=iacymd+long_gcam_timestep*10000
  end if

  write(iulog,*)'running gcam up to date ',run_up_date
  ! during runup turn off gcam restart writes for efficiency
  cdata%l(iac_cdatal_write_rest) = .false.

  do ymdtest=19750101,run_up_date,50000
     EClock(iac_EClock_ymd) = ymdtest
     call gcam_run_mod(EClock, cdata, gcami, gcamo, gcamoemis)
!    Running emissions 
     if (ymdtest.ge.20000101.and.co2flux_coupling) &
     call gcam2emisfile_run_mod( EClock, cdata, gcamoemis)
  enddo

  ! after runup turn gcam restart back on 
  cdata%l(iac_cdatal_write_rest) = .true.

  if (fast_oneway_iac_coupling) then
     iacymd=iacymd_hold
     EClock(iac_EClock_ymd) = iacymd
  end if

  end subroutine iac_init_mod

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_run_mod

! !INTERFACE:
  subroutine iac_run_mod( EClock, cdata, iaci, iaco)

! !DESCRIPTION:
! Run interface for iac

! !USES:
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real*8, pointer :: iaci(:,:)
    real*8, pointer :: iaco(:,:)

! !LOCAL VARIABLES:
    integer :: i,j,ij,h,k,iu,nx,ny,ifld,ji,jr,num
    integer :: ymdtest,iacymd,iactod,iacymd_orig,iactod_orig,iacymd_hold
    integer :: yyyy,mm,dd,yyyym1,gcamiinterpyr(15),dum1,dum2,iyr
    character(len=128) :: hfile,vname,string,filename
    integer :: ncid,dimid,varid,nmode,n,ierr
    integer :: dimidiac(2),dimidgcamo(3),dimidgcami(2),dimidglm(2),dimidglmiwh(2)
    integer :: start3(3),count3(3),start(1),count(1),datearr(12),yeararr(12),mtharr(12)
    integer :: ncidco2,londimid,latdimid,timedimid,numco2lat,numco2lon,numtime,startyr,endyr,ind,yr,mth,startday,endday,calday
    real*8, allocatable :: array3(:,:,:)
    real*8, allocatable :: array2(:,:)
    real*8, allocatable :: arin(:)
    real*8  :: fact1,fact2
    character(len=*),parameter :: subname='(iac_run_mod)'
    character(len=128) :: casename
    logical :: fast_oneway_iac_coupling
    integer :: glmyear=14990101
    integer :: day,naez,nreg,ncrop,nfld,rungcam2emisfile
    logical :: datagcam

! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------
  datagcam=.false.

  iulog = cdata%i(iac_cdatai_logunit)
  fast_oneway_iac_coupling = cdata%l(iac_cdatal_fastiac)
  casename = trim(cdata%c(iac_cdatac_casename))

  if (fast_oneway_iac_coupling) then
     ! save original ymd/tod value
     iacymd_orig = EClock(iac_EClock_ymd)
     iactod_orig = EClock(iac_Eclock_tod)

#ifdef DEBUG
       write(iulog,*) trim(subname),'current model ymd =',iacymd_orig,' tod =',iactod_orig
#endif
     !--- reset iac ymd Clock
     if (iacymd_fudge < 0) then
        iacymd_fudge = 20050101
     else
        iacymd_fudge = iacymd_fudge + 100
        if (mod(iacymd_fudge,10000) > 1231) iacymd_fudge = (iacymd_fudge/10000+1)*10000 + 0101
     endif

     EClock(iac_Eclock_ymd) = iacymd_fudge
     EClock(iac_Eclock_tod) = EClock(iac_Eclock_dt)
  endif

  ! iac time
  iacymd = EClock(iac_EClock_ymd)
  iactod = EClock(iac_Eclock_tod)
  call shr_cal_date2ymd(iacymd,yyyy,mm,dd)

  if (yyyy.eq.2096.and.mm.eq.01) &
     call shr_sys_abort(subname//' ERROR: IESM is not designed to operate past 2095')


  ! compute "alarms" 0 = off, 1 = on
  EClock(iac_EClock_Agcam) = 0
  EClock(iac_EClock_Aglm)  = 0
  EClock(iac_EClock_AclmC) = 0
  EClock(iac_EClock_Agcamsetden) = 0
  rungcam2emisfile=0
  if (iactod == EClock(iac_EClock_dt)) then   ! first timestep of day
     if (dd==1) &
          EClock(iac_EClock_AclmC) = 1   ! every month
     if (dd==1 .and. mm==1) then
        EClock(iac_EClock_Aglm)  = 1   ! every year
        EClock(iac_EClock_Agcam) = 1   ! every year
     end if
     if(dd==1 .and. mm==1 .and. yyyy>2005 .and. mod(yyyy-2005,long_gcam_timestep)==0) &
          EClock(iac_EClock_Agcamsetden) = 1   ! every long_gcam_timestep year
     if(dd==1 .and. mm==1 .and. yyyy>=2005 .and. mod(yyyy-2005,long_gcam_timestep)==0) &
          rungcam2emisfile = 1   ! every long_gcam_timestep year
  endif

#ifdef DEBUG
     write(iulog,*) trim(subname),'current model date1 ',iacymd,iactod
     write(iulog,*) trim(subname),'current model date2 ',yyyy,mm,dd
     write(iulog,*) trim(subname),'current model alarm ', &
          EClock(iac_Eclock_Agcam),EClock(iac_Eclock_Aglm),EClock(iac_Eclock_AclmC)
#endif

  if (EClock(iac_EClock_AclmC) == 1) then
#ifdef DEBUG
      write(iulog,*) trim(subname),'calling iac2gcam_run',EClock(iac_EClock_ymd),EClock(iac_EClock_tod)
#endif
     call shr_sys_flush(iulog)
     call iac2gcam_run_mod(EClock,cdata,iaci,gcami)
     call iac_diag(' gcami: ',gcami)
  endif

  if (EClock(iac_EClock_Agcamsetden) == 1) then
     ! First time through gcamiold is set to iac_spval 
     ! initialize gcamiold
     if (maxval(gcamiold).lt.0.9*iac_spval) then
        gcamiold=1.
#ifdef DEBUG
        write(iulog,*) trim(subname),'initializing gcamiold to 1.'
        write(iulog,*) trim(subname),'gcamiold max/min/sum=',maxval(gcamiold),minval(gcamiold),sum(gcamiold)
#endif
     else
        ! if clm sets a previously good point to 0 (pft disappears)
        ! hold scaling constant.
        where (gcami.eq.0. .and. gcamiold.gt.0.)
           gcami=gcamiold
        end where
     end if
     
     ! we have averaged clm scalers for each gcam timestep (5 years) 
     ! and will use consecutive averages to interpolate all interrim years thus smoothing
     ! the clm carbon scaler data. Initial scalars for first interpolation are 1.
     
     iacymd_hold = EClock(iac_Eclock_ymd)
     do i = 1,long_gcam_timestep
        EClock(iac_Eclock_ymd) = EClock(iac_Eclock_ymd)+10000
        fact1=real(long_gcam_timestep-i)/real(long_gcam_timestep)
        fact2=real(i)/real(long_gcam_timestep)
#ifdef DEBUG
        write(iulog,*) trim(subname),'setting gcam density for gcam year',EClock(iac_EClock_ymd),EClock(iac_EClock_tod),fact1
#endif
        call shr_cal_date2ymd(EClock(iac_EClock_ymd),gcamiinterpyr(i),dum1,dum2)
        ! only interp gcami for non missing data, trap 0s and missing
        write(iulog,*) trim(subname),'figuring gcamiinterp gcami max/min=',maxval(gcami),minval(gcami),sum(gcami),fact1,fact2
        write(iulog,*) trim(subname),'gcamiold max/min/sum=',maxval(gcamiold),minval(gcamiold),sum(gcamiold)
        write(iulog,*) trim(subname),'gcamiinterp sum old=',sum(gcamiinterp)
        where (gcami.gt.0.)
           gcamiinterp(:,:,i)=fact1*gcamiold(:,:)+fact2*gcami(:,:)
        elsewhere
           gcamiinterp(:,:,i)=iac_spval
        end where
        if (clm_iac_carbon_scaling) then
           write(iulog,*) trim(subname),'gcamiinterp sum new=',sum(gcamiinterp)
	   gcamiptr=>gcamiinterp(:,:,i)
           call gcam_setdensity_mod(EClock, cdata, gcamiptr)
        end if
     end do
     gcamiold=gcami
     write(iulog,*) trim(subname),'new gcamiold max/min/sum=',maxval(gcamiold),minval(gcamiold),sum(gcamiold)
     EClock(iac_Eclock_ymd) = iacymd_hold
  end if
  !
  ! write the iac restart pointer that contains the current date and gcamiold values for restart
  ! the date is read off of this file to allow the fast mode restarts at the fudged date.
  ! write every timestep for fastmode or every time gcamiold is changed
  if (EClock(iac_EClock_Agcamsetden) == 1 .or. fast_oneway_iac_coupling) &
       call write_iac_restart(iac_rpointer,iac_restfile,EClock(iac_Eclock_ymd),'gcami',gcamiold)
     
  if (EClock(iac_EClock_Agcam) == 1) then
     iacymd_hold = EClock(iac_Eclock_ymd)
     ! Simulate the long_gcam_timestep by calling gcam_run_mod one or more times
     ! advancing the model year by 5 years (gcams real timestep) each time
     EClock(iac_Eclock_ymd) = EClock(iac_Eclock_ymd) + 50000
#ifdef DEBUG
     write(iulog,*) trim(subname),'calling gcam_run_mod ',EClock(iac_EClock_ymd),EClock(iac_EClock_tod)
#endif
     call gcam_run_mod(EClock, cdata, gcami, gcamo, gcamoemis)
     call iac_diag(' gcamo: ',gcamo)
     
     if (sneakermode.and.EClock(iac_EClock_Agcamsetden)) then 
        do i = 1,long_gcam_timestep/5 - 1
           ! Simulate the long_gcam_timestep by calling gcam_run_mod one or more times
           ! advancing the model year by 5 years (gcams real timestep) each time
           EClock(iac_Eclock_ymd) = EClock(iac_Eclock_ymd) + 50000
#ifdef DEBUG
           write(iulog,*) trim(subname),'calling gcam_run_mod ',EClock(iac_EClock_ymd),EClock(iac_EClock_tod)
#endif
           call gcam_run_mod(EClock, cdata, gcami, gcamo, gcamoemis)
           call iac_diag(' gcamo: ',gcamo)
        end do
     end if
     EClock(iac_Eclock_ymd) = iacymd_hold
  end if
!
! After GCAM has completed full 5/15 year timestep write CO2 emissions to CO2 flux file
!     
  if (rungcam2emisfile == 1.and.co2flux_coupling) then
     iacymd_hold = EClock(iac_Eclock_ymd)
     EClock(iac_Eclock_ymd) = EClock(iac_Eclock_ymd) + 50000
     call gcam2emisfile_run_mod( EClock, cdata, gcamoemis)
     EClock(iac_Eclock_ymd) = iacymd_hold
  end if

  if (EClock(iac_EClock_Aglm) == 1) then
     iacymd_hold = EClock(iac_Eclock_ymd)
     EClock(iac_Eclock_ymd) = EClock(iac_Eclock_ymd) + 10000
#ifdef DEBUG
      write(iulog,*) trim(subname),'calling gcam2glm_run',EClock(iac_EClock_ymd),EClock(iac_EClock_tod)
#endif
     call gcam2glm_run_mod(EClock, cdata, gcamo, glmi, glmi_wh)
     call iac_diag(' glmi: ',glmi)
     write(iulog,'(2a,i3,2f13.6)') '(iac_diag)'//' ',' glmi_wh:',1, &
           minval(glmi_wh(:)),maxval(glmi_wh(:))
#ifdef DEBUG
      write(iulog,*) trim(subname),'calling glm_run',EClock(iac_EClock_ymd),EClock(iac_EClock_tod)
#endif
     call glm_run_mod (EClock, cdata, glmi, glmi_wh, glmo )
     EClock(iac_Eclock_ymd) = iacymd_hold
  endif

  ! Need to reverse glmo order, north->south, east->west to south->north, east->west
  ! use array2 as temporary

  allocate(array2(size(glmo,dim=1),size(glmo,dim=2)))
  array2 = glmo
  do j = 1,iac_glm_ny
  do i = 1,iac_glm_nx
    ji = (j-1)*iac_glm_nx + i
    jr = (iac_glm_ny - j)*iac_glm_nx + i
    do k = 1,iac_glmo_nflds
      glmo(k,ji) = array2(k,jr)
    enddo
  enddo
  enddo
  deallocate(array2)


  !--- tcx end fill

  if (EClock(iac_EClock_Aglm) == 1) then
     iacymd_hold = EClock(iac_Eclock_ymd)
     EClock(iac_Eclock_ymd) = EClock(iac_Eclock_ymd) + 10000
#ifdef DEBUG
        write(iulog,*) trim(subname),'calling glm2iac_run',EClock(iac_EClock_ymd),EClock(iac_EClock_tod)
#endif
     call iac_diag(' glmo: ',glmo)
     call glm2iac_run_mod(EClock, cdata, glmo, iaco)
     EClock(iac_Eclock_ymd) = iacymd_hold

!-------- history file ---------------------------------------------------
  write(hfile,'(a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)') trim(casename)//'.iac.hi.',yyyy,'-',mm,'-',dd,'-',iactod,'.nc'
  write(iulog,*) trim(subname),' writing history file ',trim(hfile)

  nmode = ior(NF90_CLOBBER,NF90_64BIT_OFFSET)
  call iac_ncerr(nf90_create(trim(hfile),nmode,ncid),subname//'create',__LINE__)
  call iac_ncerr(nf90_put_att(ncid,NF90_GLOBAL,'missing_value',iac_spval),subname//'putatt_missval',__LINE__)

  call iac_ncerr(nf90_def_dim(ncid,'iac_nx' ,cdata%i(iac_cdatai_iac_nx) ,dimid),subname//'defdim_iacnx',__LINE__)
  dimidiac(1) = dimid
  call iac_ncerr(nf90_def_dim(ncid,'iac_ny' ,cdata%i(iac_cdatai_iac_ny) ,dimidiac(2)),subname//'defdim_iacny',__LINE__)
  call iac_ncerr(nf90_def_dim(ncid,'gcami_nreg',cdata%i(iac_cdatai_gcam_nreg),dimidgcami(2)),subname//'defdim_gcaminreg',__LINE__)
  call iac_ncerr(nf90_def_dim(ncid,'gcami_naez',cdata%i(iac_cdatai_gcam_naez),dimidgcami(1)),subname//'defdim_gcaminaez',__LINE__)
  call iac_ncerr(nf90_def_dim(ncid,'gcamo_ntime',cdata%i(iac_cdatai_gcamo_ntime),dimidgcamo(3)),subname//'defdim_gcamontime',__LINE__)
  call iac_ncerr(nf90_def_dim(ncid,'gcamo_nreg',cdata%i(iac_cdatai_gcam_nreg),dimidgcamo(2)),subname//'defdim_gcamonreg',__LINE__)
  call iac_ncerr(nf90_def_dim(ncid,'gcamo_naez',cdata%i(iac_cdatai_gcam_naez),dimidgcamo(1)),subname//'defdim_gcamonaez',__LINE__)
  call iac_ncerr(nf90_def_dim(ncid,'glm_nx' ,cdata%i(iac_cdatai_glm_nx) ,dimidglm(1)),subname//'defdim_glmnx',__LINE__)
  call iac_ncerr(nf90_def_dim(ncid,'glm_ny' ,cdata%i(iac_cdatai_glm_ny) ,dimidglm(2)),subname//'defdim_glmny',__LINE__)
  call iac_ncerr(nf90_def_dim(ncid,'glmiwh_nreg' ,cdata%i(iac_cdatai_gcam_nreg) ,dimidglmiwh(2)),subname//'defdim_glminreg',__LINE__)
  call iac_ncerr(nf90_def_dim(ncid,'glmiwh_naez' ,cdata%i(iac_cdatai_gcam_naez) ,dimidglmiwh(1)),subname//'defdim_glminaez',__LINE__)

!  do n = 1,size(iaci,dim=1)
!     write(vname,'(a,i2.2)') 'iaci',n
!     ierr = nf90_def_var(ncid,vname,NF90_DOUBLE,dimidiac,varid)
!     call iac_ncerr(ierr,'defvar_'//trim(vname))
!  enddo
!  do n = 1,size(iaco,dim=1)
!     write(vname,'(a,i2.2)') 'iaco',n
!     ierr = nf90_def_var(ncid,vname,NF90_DOUBLE,dimidiac,varid)
!     call iac_ncerr(ierr,'defvar_'//trim(vname))
!  enddo

  if (EClock(iac_EClock_Agcamsetden) == 1 .and. clm_iac_carbon_scaling ) then
  do i=1,long_gcam_timestep
     num=0
     do n = 1,size(gcamiinterp,dim=1)
        do ncrop = 1,iac_gcam_ncrops
           num=num+1
           write(vname,'(a,i2.2,a,i4.4)') 'gcami_scalar',num,'_'//trim(iac_gcami_fld_names(n))//'_'//trim(iac_gcami_crop_names(ncrop))//'_',gcamiinterpyr(i)
           call iac_ncerr(nf90_def_var(ncid,vname,NF90_DOUBLE,dimidgcami,varid),subname//':defvar_'//trim(vname),__LINE__)
        enddo
     enddo
  enddo
  endif

  do n = 1,size(gcamo,dim=1)
     write(vname,'(a,i2.2,a)') 'gcamo',n,'_'//trim(iac_gcamo_fld_names(n))
     call iac_ncerr(nf90_def_var(ncid,vname,NF90_DOUBLE,dimidgcamo,varid),subname//':defvar_'//trim(vname),__LINE__)
  enddo

  do n = 1,size(glmi,dim=1)
     write(vname,'(a,i2.2,a)') 'glmi',n,'_'//trim(iac_glmi_fld_names(n))
     call iac_ncerr(nf90_def_var(ncid,vname,NF90_DOUBLE,dimidglm,varid),subname//':defvar_'//trim(vname),__LINE__)
  enddo

  write(vname,'(a,i2.2,a)') 'glmi_wh',1,'_woodharvest'
  call iac_ncerr(nf90_def_var(ncid,vname,NF90_DOUBLE,dimidglmiwh,varid),subname//':defvar_'//trim(vname),__LINE__)

  do n = 1,size(glmo,dim=1)
     write(vname,'(a,i2.2,a)') 'glmo',n,'_'//trim(iac_glmo_fld_names(n))
     call iac_ncerr(nf90_def_var(ncid,vname,NF90_DOUBLE,dimidglm,varid),subname//':defvar_'//trim(vname),__LINE__)
  enddo

  call iac_ncerr(nf90_enddef(ncid),subname//':enddef',__LINE__)

  if (EClock(iac_EClock_Agcamsetden) == 1 .and. clm_iac_carbon_scaling ) then
     iacymd_hold = EClock(iac_Eclock_ymd)
     allocate(array3(iac_gcam_ncrops,cdata%i(iac_cdatai_gcam_naez),cdata%i(iac_cdatai_gcam_nreg)))
     do iyr=1,long_gcam_timestep
        num=0
        do n = 1,size(gcamiinterp,dim=1)
           ij = 0
           do j = 1,cdata%i(iac_cdatai_gcam_nreg)
              do i = 1,cdata%i(iac_cdatai_gcam_naez)
                 do ncrop = 1,iac_gcam_ncrops
                    ij = ij+1
                    array3(ncrop,i,j) = gcamiinterp(n,ij,iyr)
                 enddo
              enddo
           enddo
           
           do ncrop = 1,iac_gcam_ncrops
              num=num+1
              write(vname,'(a,i2.2,a,i4.4)') 'gcami_scalar',num,'_'//trim(iac_gcami_fld_names(n))//'_'//trim(iac_gcami_crop_names(ncrop))//'_',gcamiinterpyr(iyr)
              call iac_ncerr(nf90_inq_varid(ncid,vname,varid),subname//':inqvar_'//trim(vname),__LINE__)
              call iac_ncerr(nf90_put_var(ncid,varid,array3(ncrop,:,:)),subname//':putvar_'//trim(vname),__LINE__)
           enddo
        enddo
     enddo
     EClock(iac_Eclock_ymd) = iacymd_hold
     deallocate(array3)
  end if

  allocate(array3(cdata%i(iac_cdatai_gcam_naez),cdata%i(iac_cdatai_gcam_nreg),cdata%i(iac_cdatai_gcamo_ntime)))
  do n = 1,size(gcamo,dim=1)
     ij = 0
     do k = 1,cdata%i(iac_cdatai_gcamo_ntime)
        do j = 1,cdata%i(iac_cdatai_gcam_nreg)
           do i = 1,cdata%i(iac_cdatai_gcam_naez)
              ij = ij+1
              array3(i,j,k) = gcamo(n,ij)
           enddo
        enddo
     enddo
     write(vname,'(a,i2.2,a)') 'gcamo',n,'_'//trim(iac_gcamo_fld_names(n))
     call iac_ncerr(nf90_inq_varid(ncid,vname,varid),subname//':inqvar_'//trim(vname),__LINE__)
     call iac_ncerr(nf90_put_var(ncid,varid,array3),subname//':putvar_'//trim(vname),__LINE__)
  enddo
  deallocate(array3)

  allocate(array2(cdata%i(iac_cdatai_glm_nx),cdata%i(iac_cdatai_glm_ny)))
  do n = 1,size(glmi,dim=1)
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
     do i = 1,cdata%i(iac_cdatai_glm_nx)
        ij = ij+1
        array2(i,j) = glmi(n,ij)
     enddo
     enddo
     write(vname,'(a,i2.2,a)') 'glmi',n,'_'//trim(iac_glmi_fld_names(n))
     call iac_ncerr(nf90_inq_varid(ncid,vname,varid),subname//':inqvar_'//trim(vname),__LINE__)
     call iac_ncerr(nf90_put_var(ncid,varid,array2),subname//':putvar_'//trim(vname),__LINE__)
  enddo

  do n = 1,size(glmo,dim=1)
     ij = 0
     do j = 1,cdata%i(iac_cdatai_glm_ny)
     do i = 1,cdata%i(iac_cdatai_glm_nx)
        ij = ij+1
        array2(i,j) = glmo(n,ij)
     enddo
     enddo
     write(vname,'(a,i2.2,a)') 'glmo',n,'_'//trim(iac_glmo_fld_names(n))
     call iac_ncerr(nf90_inq_varid(ncid,vname,varid),subname//':inqvar_'//trim(vname),__LINE__)
     call iac_ncerr(nf90_put_var(ncid,varid,array2),subname//':putvar_'//trim(vname),__LINE__)
  enddo
  deallocate(array2)

  allocate(array2(cdata%i(iac_cdatai_gcam_naez),cdata%i(iac_cdatai_gcam_nreg)))
  ij = 0
  do j = 1,cdata%i(iac_cdatai_gcam_nreg)
     do i = 1,cdata%i(iac_cdatai_gcam_naez)
        ij = ij+1
        array2(i,j) = glmi_wh(ij)
     enddo
  enddo
  write(vname,'(a,i2.2,a)') 'glmi_wh',1,'_woodharvest'
  call iac_ncerr(nf90_inq_varid(ncid,vname,varid),subname//':inqvar_'//trim(vname),__LINE__)
  call iac_ncerr(nf90_put_var(ncid,varid,array2),subname//':putvar_'//trim(vname),__LINE__)
  deallocate(array2)

  call iac_ncerr(nf90_close(ncid),subname//':close',__LINE__)

!-------- end history file ------------
  
  endif ! Aglm = 1

  if (fast_oneway_iac_coupling) then
     EClock(iac_EClock_ymd) = iacymd_orig
     EClock(iac_Eclock_tod) = iactod_orig
  endif

  end subroutine iac_run_mod

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_final_mod

! !INTERFACE:
  subroutine iac_final_mod( )

! !DESCRIPTION:
! Finalize iac model

!------------------------------------------------------------------------------

   implicit none
! !ARGUMENTS:

! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='(iac_final_mod)'

! !REVISION HISTORY:
! Author: T Craig

!EOP
!---------------------------------------------------------------------------

  !  Cleanup GCAM
  call gcam_final_mod()
  call glm_final_mod()
  deallocate(gcamiinterp)
  deallocate(gcamiold)
  end subroutine iac_final_mod

!====================================================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_diag

! !INTERFACE:
  subroutine iac_diag(string,array)

! !DESCRIPTION:
! iac array diagnostic

!------------------------------------------------------------------------------

   implicit none
! !ARGUMENTS:
    character(len=*) :: string
    real*8, pointer :: array(:,:)

! !LOCAL VARIABLES:
    integer :: nj,j
    character(len=*),parameter :: subname='(iac_diag)'

! !REVISION HISTORY:
! Author: T Craig

!EOP
!---------------------------------------------------------------------------

     nj = size(array,dim=1)
     do j = 1,nj
        write(iulog,'(2a,i3,2f13.6)') trim(subname)//' ',trim(string),j, &
              minval(array(j,:)),maxval(array(j,:))
     enddo

  end subroutine iac_diag

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: iac_ncerr

! !INTERFACE:
  subroutine iac_ncerr( ret, mes, line )

! !DESCRIPTION:
! iac netcdf error diagnostic

!------------------------------------------------------------------------------

   implicit none
! !ARGUMENTS:
    integer, intent(in) ::&
         ret                 ! return code from netCDF library routine
    character(len=*), intent(in) ::&
         mes                 ! message to be printed if error detected
    integer, intent(in), optional :: line

    !-----------------------------------------------------------------------

    if ( ret .ne. NF90_NOERR ) then
       if(present(line)) then
          write(6,*) mes, line
       else   
          write(6,*) mes
       end if
       write(6,*) nf90_strerror( ret )
       !jt        call endrun ('HANDLE_NCERR')
       call abort
    endif

    return
  end subroutine iac_ncerr

!====================================================================================
!---------------------------------------------------------------------------
!BOP
  
! !IROUTINE: write_iac_restart

! !INTERFACE:
  subroutine write_iac_restart(rpointfile,rstfile,ymd,varname,var)
    
    ! !DESCRIPTION:
    ! write restart for iac model
    ! !USES:
    use shr_cal_mod
    
    implicit none
    
    ! !ARGUMENTS:
    integer :: ymd
    real(r8) :: var(:,:)
    character(len=*) :: varname
    character(len=*) :: rpointfile
    character(len=*) :: rstfile
    
    ! !LOCAL VARIABLES:
    integer :: year,mon,day,iun,status,dimid2(2),ncid,varid,varidymd,dimid,nmode
    character(len=*),parameter :: subname='(write_iac_restart)'
    character(256) :: filename
    
    ! !REVISION HISTORY:
    ! Author: T Craig
    
    !EOP
    
    !---------------------------------------------------------------------------
    
    
    call shr_cal_date2ymd(ymd,year,mon,day)
    write(filename,'(a,i4.4,a,i2.2,a)') trim(rstfile)//'r.',year,'.nc'
    
    iun = shr_file_getunit()
    open(iun,file=trim(rpointfile),form='formatted')
    write(iun,'(a)') trim(filename)
    close(iun)
    call shr_file_freeunit(iun)
    
    write(iulog,*) subname,' write_restart rpointer ',trim(rpointfile)
    write(iulog,*) subname,' write_restart file     ',trim(filename)
    
    nmode = ior(NF90_CLOBBER,NF90_64BIT_OFFSET)

    call iac_ncerr(nf90_create(filename,nmode,ncid),subname,__LINE__)
    call iac_ncerr(nf90_def_dim(ncid,'dim1',size(var,dim=1),dimid2(1)),subname,__LINE__)
    call iac_ncerr(nf90_def_dim(ncid,'dim2',size(var,dim=2),dimid2(2)),subname,__LINE__)
    call iac_ncerr(nf90_def_var(ncid,varname,NF90_DOUBLE,dimid2,varid),subname,__LINE__)
    call iac_ncerr(nf90_def_var(ncid,'ymd',NF90_INT,varidymd),subname,__LINE__)
    
    call iac_ncerr(nf90_enddef(ncid),subname,__LINE__)
    
    call iac_ncerr(nf90_inq_varid(ncid,varname,varid),subname,__LINE__)
    call iac_ncerr(nf90_put_var(ncid,varid,var),subname,__LINE__)
    call iac_ncerr(nf90_put_var(ncid,varidymd,ymd),subname,__LINE__)
    call iac_ncerr(nf90_close(ncid),subname,__LINE__)
    
  end subroutine write_iac_restart

!====================================================================================

end module iac_comp_mod
