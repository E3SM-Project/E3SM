module mpaso_DA_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: mpaso_DA_mod
! ! interfacce:
      
! !DESCRIPTION:
!   !This is the online DA program for the Model for Predication Across Scales Ocean Model (MPASO).
!   !REVISION HISTORY:

!   !USES:
    use mct_mod
    use esmf
    use seq_timemgr_mod
    use perf_mod
    use shr_file_mod
    use shr_cal_mod,       only : shr_cal_date2ymd
    use shr_sys_mod,  only : shr_sys_flush
!    use shr_kind_mod,      only : r8=> shr_kind_r8                     
!   !mpaso
    use mpaso_cpl_indices
    use mpaso_mct_vars

    use mpas_framework
    use mpas_derived_types
    use mpas_pool_routines
    use mpas_io
    use mpas_io_streams
    use mpas_stream_manager
    use mpas_kind_types
    use mpas_io_units
    use mpas_timekeeping
    use mpas_bootstrapping
    use mpas_dmpar
    use mpas_constants
    use mpas_log
    
!   ! DA module
    use seq_comm_mct,  only : logunit
    use seq_comm_mct,  only : num_inst_driver
    use seq_comm_mct,  only : DAOCNID,DAOCNALLID
    use seq_comm_mct,  only : seq_DA_comm_iamin
    use seq_comm_mct,  only : seq_DA_comm_iamroot
    use seq_comm_mct,  only : seq_DA_comm_iam
    use seq_comm_mct,  only : seq_DA_comm_mpicom
    use seq_comm_mct,  only : seq_DA_comm_inst
    use seq_comm_mct,  only : seq_DA_comm_npes

    use seq_comm_mct,  only : OCNID
    use seq_comm_mct,  only : seq_comm_iamin
    use seq_comm_mct,  only : seq_comm_iamroot
    use seq_comm_mct,  only : seq_comm_iam
    use seq_comm_mct,  only : seq_comm_mpicom

    use ECDA_mpaso_state_mod
    use ECDA_mpaso_obs_mod

#ifdef _MPI
#ifndef NOMPIMOD
    use mpi
#endif
#endif

    implicit none

#ifdef _MPI
#ifdef NOMPIMOD
    include 'mpif.h'
#endif
#endif
    integer, parameter :: MPI_INTEGERKIND = MPI_INTEGER
    integer, parameter :: MPI_2INTEGERKIND = MPI_2INTEGER

#ifdef SINGLE_PRECISION
    integer, parameter :: MPI_REALKIND = MPI_REAL
    integer, parameter :: MPI_2REALKIND = MPI_2REAL
#else
    integer, parameter :: MPI_REALKIND = MPI_DOUBLE_PRECISION
    integer, parameter :: MPI_2REALKIND = MPI_2DOUBLE_PRECISION
#endif

!
!   !PUBLIC MEMBER FUNCTIONS:

    public :: mpaso_DA_init
    public :: mpaso_DA_run
!    public :: state_update_IAU
!    public :: ocean_ave4da

!   Public variables !
!
    
!    integer,public             ::  my_inst, ensemble_size
!    integer :: allCellsSolve, DACellsSolve,DACellsSolve_halo
!    integer :: da_nts,da_nte, da_nbs,da_nbe   ! [da_nbs,da-nbe]=[da_nts da_nte]+halo
!EOP


!   !Private by default
    private

    type (domain_type), pointer:: domain 

!    variable for Data assimilation, move to ECDA_mpaso_state_mod
!    integer :: mpicom_ocn_local           ! mpaso model communicator
!    integer :: mpicom_da_all,mpicom_da  !comunicator
!    integer :: iam_da_all, iam_da       ! comunicater ID
!    logical :: iamin_da_all,iamin_da
!    logical :: iamroot_da_all,iamroot_da ! root
    
    integer :: iam_local                ! local comunicater

    integer :: init_obs_YR               ! ex:2016
    integer :: init_obs_MO              ! ex: 1~12
    integer :: init_obs_DY             ! ex:1~31
    integer :: init_obs_HR              ! ex:0~23
    integer :: init_obs_doy            

    type(ESMF_Time) :: obs_init_Time                  ! use ESMF to manage observation time, can it be public? @Yun
!  DA namelist
   character(len = 10) :: init_obs_YMDH = '1980010112'   ! start time for DA ex:2010010103
   integer :: DA_interv = 1                    ! By default the obervation record time interval = da_interv, record 1 on the time of init_obs_YMDH, record 2 at time init_obs_YMDH+DA_interv
   character(len =2) :: DA_interv_unit = 'dy'  ! it can be hour:'hr' or day:'dy', 'mo',month. for hours, it should be full divided by day: 1,2,3,4,6,8,12,24 for DA_interv
   integer :: obs_interv  = 1               ! @YL: the observation interval, may small than da_interv (for sensitivty test),  use the same unit as da_interv for now.
   logical :: update_states(n_OCN_state)= .False.    ! define the states to be updated by observation, default is no update
   real(kind=RKIND) :: relax_inf_rate=0.5                      ! relaxiation to prior weight on prior, for RTPS,RTPP 
   !- increment analysis updates (IAU) to prevent the model crash by abrupt modifying model state. (analysis increment add into model within multiple couple interval)
   !   need know the couple interval and ocean integrate interval. @YL to be develop
   integer :: IAU_cpl_cycles = 0    ! 0, default, no IAU applied

   namelist /ocn_DA_nml/init_obs_YMDH, DA_interv, DA_interv_unit,obs_interv, &
                           update_states,relax_inf_rate,&
                          IAU_cpl_cycles

contains
!==============================================================================
!BOP
!   ! get the basic grid for current PE for DA
!   !ROUTINE: mpaso_DA_init
!   !INTERFACE:

subroutine mpaso_DA_init(mpicom_ocn,domain_ptr)
!BOC
      implicit none
      integer,intent(in) :: mpicom_ocn
      type (domain_type), pointer:: domain_ptr 
!  begin
      character(*), parameter :: subName =   '(mpaso_DA_init) '
      character(*), parameter :: nml_inpfile = 'ocean_da.namelist'
      type(ESMF_CalKind_Flag) :: esmf_calkind
      type(ESMF_Calendar) :: esmf_cal
      integer :: init_obs_sod                    ! seconds of the day.
      integer :: ierr,iunit,rc
      logical :: exists
      integer :: ymd, tod, ihour, iminute, isecond
      integer :: iyear, imonth, iday, curr_ymd, curr_tod
      integer :: shrloglev, shrlogunit
      real (kind=RKIND) :: dt, current_wallclock_time
      integer :: ivar

!      call seq_cdata_setptrs(cdata_o, ID=OCNID, mpicom=mpicom_o, &
!         gsMap=gsMap_o, dom=dom_o, infodata=infodata)

      mpicom_ocn_local=mpicom_ocn
      call MPI_comm_rank(mpicom_ocn_local,iam_local,ierr) 
      domain => domain_ptr    ! acess the domain at current PE
      if (num_inst_driver <= 1) then
          write(logunit,*) trim(subname), 'single ocean copy, no ensmeble DA, OI to be developed '
          return
      endif
      ensemble_size = num_inst_driver         !       multi-driver for DA ensemble

!      ocnDALogUnit = shr_file_getUnit()      
      ! open the logfile for DA
!      if(iam_local==0) then
!        open( ocnDALogUnit, file = './mpaso_da_log')
!      endif

      call get_OCN_DA_comm()                
      call mpi_barrier(mpicom_da_all, ierr)   !       make sure all the ocn pes reach to here 

!-      read in the name list for ocean da
        inquire(file=trim(nml_inpfile), exist=exists)
        if (exists) then
             iunit = shr_file_getUnit()
             open(unit = iunit, file = trim(nml_inpfile), action = 'read')
             rewind(unit = iunit)
             read(unit = iunit, nml = ocn_DA_nml, iostat = ierr)
             close(unit = iunit)
             call shr_file_freeUnit( iunit )
             if (ierr /= 0) write(logunit,*) trim(subname),'ERROR reading namelist ocn_DA_nml from ',trim(nml_inpfile)
        end if

!-     convert the observation initial time to ESMF format
        read(init_obs_YMDH(1:4),'(i)') init_obs_YR
        read(init_obs_YMDH(5:6),'(i)') init_obs_MO
        read(init_obs_YMDH(7:8),'(i)') init_obs_DY
        read(init_obs_YMDH(9:10),'(i)') init_obs_HR
        esmf_calkind = ESMF_CALKIND_NOLEAP   ! @YL, need input from domain%clock
        !esmf_calkind = ESMF_CALKIND_GREGORIAN   ! default for
        esmf_cal = ESMF_CalendarCreate( name="mpaso_obs_calendar", calkindflag=esmf_calkind, rc=rc )
        init_obs_sod= init_obs_HR * 3600
        call ESMF_TimeSet( obs_init_Time, yy=init_obs_YR, mm=init_obs_MO, dd=init_obs_DY, s=init_obs_sod, calendar=esmf_cal, rc=rc )
        call ESMF_TimeGet( obs_init_Time, dayOfYear=init_obs_doy, rc=rc )
       ! write(logunit,*) trim(subname),'obs initial time day of year',init_obs_doy

!-       initial state, observation         
        call ocn_DAstate_init() ! initial state for DA
        call ocn_obs_type_register() ! register assimilating observations type,need after ocn_DAstate_init
        call Define_OCN_DA_grids() 

        if(iamroot_da_all) write(logunit,*) trim(subname),'NUM_OCN_obs_type:',NUM_OCN_obs_type
        if (NUM_OCN_obs_type==0) return   ! no valid observation type. @YL

!-     set ocn state for DA
        ! all the updated variables are required

       if(iamroot_da_all) write(logunit,*) trim(subname),'update states:',update_states
        do ivar=1, n_OCN_state
!           OCN_DAstate_data(ivar)%averaged = average_obs
           if(update_states(ivar)) then
              OCN_DAstate_data(ivar)%r2p_fac = relax_inf_rate
              OCN_DAstate_data(ivar)%required = .true.
              OCN_DAstate_data(ivar)%update = .true.
           endif
        end do
        call ocn_DAstate_register()
!-    test read in obs
!        call read_ocn_obs4da(1,exists)

!-    No IAU applied ,set IAU weight as 1 @YL, for IAU need to calculate based on coupling timesteps 
       steps_IAU=1

       call mpi_barrier(mpicom_da_all, ierr)   !       make sure all the ocn pes reach to here
       if(iamroot_da_all) write(logunit,*) trim(subname), 'finsh OCN DA initial'

    end subroutine mpaso_DA_init
!==============================================================================
!BOP
!   !ROUTINE: mpaso_DA_run
!   !INTERFACE:

    subroutine mpaso_DA_run(MPASclock)
!BOC

      implicit none
      type (MPAS_Clock_type), intent(in) :: MPASclock

!  begin
      character(*), parameter :: subName =   '(mpaso_DA_run) '
      logical :: isdatime
      integer :: shrloglev, shrlogunit
      integer :: ierr,i,j,k,n
!     observation related
      integer :: obs_rec      ! the file record for observations @YL
      logical ::  obs_exists 

      type (MPAS_Time_Type) :: currTime
      type (domain_type), pointer :: domain_ptr
      character(len=StrKIND) :: timeStamp, streamName, WCstring
      real (kind=RKIND) :: pertubation
      

      if (num_inst_driver <= 1) return

      if (NUM_OCN_obs_type==0) return

!       IAU update in multiple coupling steps, when it reach to 0, stop IAU.
        IAU_update= IAU_update -1

      !if reach to assimilation time @YL, maybe use the alarm way
      call is_ocn_da_time(isdatime,obs_rec,MPASclock)
      if( .not. isdatime)  return

      call mpi_barrier(mpicom_da_all, ierr)   !       make sure all the ocn pes reach to here

!     read in observations

!      obs_rec = 1
      call read_ocn_obs4da(obs_rec,obs_exists)

      call collect_mpaso_state4da()

      call mpi_barrier(mpicom_da_all, ierr)   !       make sure all the ocn pes reach to here

      call ecda_eakf_core()    !do the assimilation

      !=    inflation after DA
      call mpas_relax_inflation()     ! relaxiation  inflation  with RPTS

      !=    calculate increment for IAU
      call ecda_calculate_IAU()
     
      != update state
      IAU_update = IAU_cpl_cycles ! set IAU update coupling cycles
!      write(logunit,*)trim(subname),'IAU_update',IAU_update,'steps_IAU',steps_IAU
!    if IAU applied, update all the increment on model states
      if(steps_IAU == 1) then
          IAU_update=1
          call mpaso_update_IAU()  ! no IAU selection, Update the all increment ocean
          IAU_update = 0           ! make sure no additional update happens.
       endif 

!=   release DA related variable
      call  ocn_DAstate_deallocate()
      call observation_deallocate()

      !if(average_obs) call average_state_data_initial()  !reinitial the average variables
      call mpi_barrier(mpicom_da_all, ierr)   !       make sure all the ocn pes reach to here
     if(iamroot_da_all) write(logunit,*)trim(subname),'Finished '
      flush(logunit)
      
      return

    end subroutine mpaso_DA_run

!==============================================================================
!BOP
!   !ROUTINE: ecda_eakf_core
!   ! the main subroutine of ECDA, use EAKF
!   ! 
!   ! loop through all assimilated obs data
!   !     step 1  observation oprerator followed by calculate the ensemble increment at observation space
!   !     Step 2 loop through all updated states to calculate the ensemble increment at model space.

    subroutine ecda_eakf_core()

       use sEnKF_mod,     only : obs_inflation
       use sEnKF_mod,     only : obs_increment,update_from_obs_inc
       use sEnKF_mod,     only : cal_lcwt_latlon
!BOC
       implicit none
       include 'mpif.h'
       character(*), parameter :: subName =   '(ecda_eakf_core) '
       integer                 :: ID
       integer                 :: i,j,k, iobs,ivar,ierr
       integer                 :: icell_edge
       real(kind=RKIND)        :: fobs(ensemble_size),obs_inc(ensemble_size)
       real(kind=RKIND)        :: obs,obs_var
       real(kind=RKIND)        :: fstate(ensemble_size),state_inc(ensemble_size)
       real(kind=RKIND)        :: rlon,rlat
       real(kind=RKIND)        :: loc_wt0,loc_wt,vertical_distance, horiz_distance

       do ID=1,NUM_OCN_obs_type
          !write(logunit,*) trim(subname),trim(OCN_obs_data(ID)%obs_name),OCN_obs_data(ID)%assimilate
          if(.not. OCN_obs_data(ID)%assimilate .or. OCN_obs_data(ID)%nobs <=0) cycle
          !write(logunit,*) trim(subname),trim(OCN_obs_data(ID)%obs_name),'for obs of', OCN_obs_data(ID)%nobs, &
          !               'update:',OCN_obs_data(ID)%update_states
           do iobs=1, OCN_obs_data(ID)%nobs
               if (abs(OCN_obs_data(ID)%values(iobs))>999.0) cycle  ! an ad hoc value
               !  @YL skip ice cover region
               !if(OCN_obs_data(ID)%obs_type == OCN_TEMP_index .and. &
               !           OCN_obs_data(ID)%values(iobs) <= 5.0  ) cycle
               call  mpaso_obs_operator(ID,iobs,fobs,obs,obs_var, icell_edge)
               !if(iobs==10 .and. iamroot_da ) &
               !       write(logunit,*) trim(subname),'iobs,icell,obs, fobs,location',iobs,icell_edge,obs,fobs,&
               !       OCN_obs_data(ID)%lon(iobs),OCN_obs_data(ID)%lat(iobs),DA_lonCell(icell_edge),DA_latCell(icell_edge)
               !--eakf step 1
               call obs_increment(fobs, obs,ensemble_size, obs_var, obs_inc)
               !write(logunit,*) trim(subname),'iobs,da inc',iobs,obs_inc
               !--eakf step two
               do ivar=1,n_OCN_state
                   if(.not. OCN_obs_data(ID)%update_states(ivar)) cycle
                   do i=da_Cell_nbs,da_Cell_nbe
                   !do i=icell_edge,icell_edge
                      ! the horizontal distance between obs and state for localization
                       call cal_lcwt_latlon(DA_lonCell(i),DA_latCell(i),OCN_obs_data(ID)%lon(iobs),OCN_obs_data(ID)%lat(iobs), &
                               OCN_obs_data(ID)%RR_xy, loc_wt0)

                       if (loc_wt0 < 0.05) cycle 
                       if (OCN_DAstate_data(ivar)%state_dim==1) then
                          fstate=OCN_DAstate_data(ivar)%posterior1d(i,:)
                          call update_from_obs_inc(fobs,obs_inc, fstate, ensemble_size, state_inc, loc_wt0)
                          OCN_DAstate_data(ivar)%posterior1d(i,:)=OCN_DAstate_data(ivar)%posterior1d(i,:)+state_inc
                       endif    ! state_dim==1

                       if (OCN_DAstate_data(ivar)%state_dim==2) then
                           do k=1,DA_maxLevelCell(i)       !-1       !@YL, skip the bottom    
                                vertical_distance=abs(OCN_obs_data(ID)%depth(iobs)-DA_depth(k))
                                if (vertical_distance > 3*OCN_obs_data(ID)%RR_z) CYCLE  ! exp(-3) ~= 0.05, weight will be small

                                loc_wt=exp(-1.0*vertical_distance/OCN_obs_data(ID)%RR_z)*loc_wt0
                                !if(iobs==10 .and. iamroot_da) write(logunit,*) trim(subname),'dpth/wt',k,DA_depth(k),loc_wt0,loc_wt  
                             
                                if (loc_wt < 0.05) cycle 
                                fstate=OCN_DAstate_data(ivar)%posterior2d(i,k,:)
                                call update_from_obs_inc(fobs,obs_inc, fstate, ensemble_size, state_inc, loc_wt)
                                OCN_DAstate_data(ivar)%posterior2d(i,k,:)=OCN_DAstate_data(ivar)%posterior2d(i,k,:)+state_inc
                            end do
                      endif    ! state_dim ==2
                   enddo  ! i
               enddo  ! ivar
           enddo   ! iobs

           !write(logunit,*) trim(subname),'increment', state_inc
       end do !ID

       call mpi_barrier(mpicom_da_all, ierr)   !       make sure all the ocn pes reach to here
       if(iamroot_da_all) write(logunit,*)trim(subname),'finished'
       flush(logunit)
    end subroutine ecda_eakf_core

!==============================================================================
!BOP
!   !ROUTINE: mpaso_update_IAU
!   !INTERFACE:

    subroutine mpaso_update_IAU()
!BOC
      implicit none
!  begin
      character(*), parameter :: subName =   '(mpaso_update_IAU) '
      logical :: isdatime
      integer :: shrloglev, shrlogunit
      integer :: ierr,i,j,k,n,ivar
!     temporal for testing communication
      real (kind=RKIND), dimension(:,:),allocatable :: state_da_2D, state_da_2D_all
      real (kind=RKIND), dimension(:,:,:),allocatable :: state_da_3D

      type (block_type), pointer :: block_ptr    
      integer, pointer :: nCellsSolve,nVertLevels
      integer, pointer :: index_temperature, index_salinity
!     define array      
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: SSHValue
      real (kind=RKIND), dimension(:,:), pointer :: TValue
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      integer :: idx
      type (field1DReal), pointer :: SSHfield
      type (field2DReal), pointer :: normalVelocity

      type (mpas_pool_type), pointer :: meshPool, statePool, diagnosticsPool, forcingPool, averagePool, scratchPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_iterator_type) :: groupItr

!     to add the analysis increment to the state using IAU, One should consider for multi steps of model integration
      if(IAU_update <= 0) RETURN   ! not updat time
      
      n=0
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
          call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
          call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature) !  OCN_TEMP_index
          call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)       ! OCN_SALT_index
  !        call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity,1)           ! OCN_VEL_index
          call mpas_pool_get_array(statePool, 'ssh', SSHValue,1)                            ! OCN_SSH_index

          if(iamroot_da_all) write(logunit,*)trim(subname),'before update sst', activeTracers(index_temperature,1,1), &
                             OCN_DAstate_data(1)%IAU_increm2D(n+1,1)
          do i = 1, nCellsSolve
             n = n + 1
             do ivar=1,n_OCN_state
                if(.not. OCN_DAstate_data(ivar)%update) cycle
                select case (ivar)
                    case (OCN_TEMP_index)
                            idx=index_temperature
                            do k=1,maxLevelCell(i)  
                               activeTracers(idx,k,i)=activeTracers(idx,k,i) + OCN_DAstate_data(ivar)%IAU_increm2D(n,k)
                            end do
                    case (OCN_SALT_index)
                            idx=index_salinity
                            do k=1,maxLevelCell(i)  
                               activeTracers(idx,k,i)=activeTracers(idx,k,i) + OCN_DAstate_data(ivar)%IAU_increm2D(n,k)
                            end do
                    case (OCN_SSH_index)
                            SSHValue(i)=SSHValue(i) +OCN_DAstate_data(ivar)%IAU_increm1D(n)

                    case (OCN_VEL_index)
                            write(logunit,*) trim(subname),'OCN_VELOCITY to be develop'

                    CASE DEFAULT
                        write(logunit,*) trim(subname),'DA state not valide:', ivar
                end select
             end do ! ivar,n_OCN_state
          end do  ! i, nCellsSolve

         if(iamroot_da_all) write(logunit,*)trim(subname),'after update sst', activeTracers(index_temperature,1,1)

          block_ptr => block_ptr % next
      enddo
      call mpi_barrier(mpicom_ocn_local, ierr)

!   exchange halo field 
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
         ! for regular state
         if (OCN_DAstate_data(OCN_SSH_index)%update) then 
             call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
             call mpas_pool_get_field(statePool, 'ssh', SSHfield)
             call mpas_dmpar_exch_halo_field(SSHfield)
          end if
         ! for tracers : temperature and salinity
         if (OCN_DAstate_data(OCN_TEMP_index)%update .or. OCN_DAstate_data(OCN_SALT_index)%update) then
              call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
              call mpas_pool_begin_iteration(tracersPool)
               do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
                  if (groupItr%memberType == MPAS_POOL_FIELD) then
                      call mpas_dmpar_field_halo_exch(domain,  groupItr%memberName, timeLevel=1)
                  end if
               end do
         end if

         block_ptr => block_ptr % next
      enddo 
     
    end subroutine mpaso_update_IAU

!==============================================================================
!BOP
!   !ROUTINE: ecda_calculate_IAU
!   ! calculate IAU and redistribute to each model PEs

!   !INTERFACE:
    subroutine  ecda_calculate_IAU()
!EOP
!BOC
        ! --- Local variables ---
        character(*), parameter :: subName =   '(ecda_calculate_IAU) '
        integer :: i, ierr, j, k, lu, lv, ien, iz,ix,iy, ivar
        real(kind=RKIND), dimension(:), allocatable :: increment1d
        real(kind=RKIND), dimension(:,:), allocatable :: increment2d
        integer ::  da_nbs,da_nbe


!EOC
        call mpi_barrier(mpicom_da_all, ierr)
        da_nbs=da_Cell_nbs
        da_nbe=da_Cell_nbe

        allocate(increment1d(da_nbs:da_nbe), &
                increment2d(da_nbs:da_nbe,DA_nvert), stat=ierr)


        do ivar=1,n_OCN_state   ! loop through all  variables
           !write(logunit,*) trim(subname),trim(OCN_DAstate_data(ivar)%name),OCN_DAstate_data(ivar)%update
           if(.not.  OCN_DAstate_data(ivar)%update) cycle


           if (OCN_DAstate_data(ivar)%state_dim==1) then
               increment1d=0.0
               do i=da_nbs,da_nbe
                  OCN_DAstate_data(ivar)%posterior1d(i,my_inst)=max(OCN_DAstate_data(ivar)%posterior1d(i,my_inst),OCN_DAstate_data(ivar)%mn_value)
                  OCN_DAstate_data(ivar)%posterior1d(i,my_inst)=min(OCN_DAstate_data(ivar)%posterior1d(i,my_inst),OCN_DAstate_data(ivar)%mx_value)
                  increment1d(i)=OCN_DAstate_data(ivar)%posterior1d(i,my_inst)-OCN_DAstate_data(ivar)%prior1d(i,my_inst)
                  increment1d(i)=sign(min(Mx_da_increm(ivar),abs( increment1d(i))),increment1d(i))
               enddo
               if(iamroot_da_all) write(logunit,*) trim(subname),trim(OCN_DAstate_data(ivar)%name),' increment max',maxval(abs(increment1d)), maxloc(abs(increment1d))
               OCN_DAstate_data(ivar)%IAU_increm1D=0.0
               ! halo issue 
               ! call  ensem_async_IAU2d(increment1d,OCN_DAstate_data(ivar)%IAU_increm1d)

               OCN_DAstate_data(ivar)%IAU_increm1D=increment1d /steps_IAU
           endif

!           if (OCN_DAstate_data(ivar)%state_dim==1) then
!              increment1d=OCN_DAstate_data(ivar)%posterior1d-OCN_DAstate_data(ivar)%prior1d
!               do i=da_nbs,da_nbe
!                  do ien=1,ensemble_size
!                           increment1d(i,ien)=sign(min(Mx_da_increm(ivar),abs( increment1d(i,ien))),increment1d(i,ien))
!                  enddo
!               enddo
!               if(iamroot_da_all) write(logunit,*) trim(subname),trim(OCN_DAstate_data(ivar)%name),' increment max',maxval(abs(increment1d))
!           !   call   ensem_async_IAU2d(increment1d,OCN_DAstate_data(ivar)%IAU_increm1d)
!               OCN_DAstate_data(ivar)%IAU_increm1D=increment1d(:,my_inst) /steps_IAU
!           endif

           if (OCN_DAstate_data(ivar)%state_dim==2) then
               increment2d=0.0
               do i=da_nbs,da_nbe
                   do k=1,DA_maxLevelCell(i)
                      OCN_DAstate_data(ivar)%posterior2d(i,k,my_inst)=max(OCN_DAstate_data(ivar)%posterior2d(i,k,my_inst),OCN_DAstate_data(ivar)%mn_value)
                      OCN_DAstate_data(ivar)%posterior2d(i,k,my_inst)=min(OCN_DAstate_data(ivar)%posterior2d(i,k,my_inst),OCN_DAstate_data(ivar)%mx_value)
                      increment2d(i,k)=OCN_DAstate_data(ivar)%posterior2d(i,k,my_inst)-OCN_DAstate_data(ivar)%prior2d(i,k,my_inst)
                      increment2d(i,k)=sign(min(Mx_da_increm(ivar),abs(increment2d(i,k))),increment2d(i,k))
                   enddo
               enddo
              if(iamroot_da_all) write(logunit,*) trim(subname),trim(OCN_DAstate_data(ivar)%name),' increment max/loc',maxval(abs(increment2d)),maxloc(abs(increment2d))
              OCN_DAstate_data(ivar)%IAU_increm2D=0.0
              !for halo @YL
              !call ensem_async_IAU3d(increment3d,OCN_DAstate_data(ivar)%IAU_increm)
              OCN_DAstate_data(ivar)%IAU_increm2D=increment2d /steps_IAU
           endif

!           if (OCN_DAstate_data(ivar)%state_dim==2) then
!               increment2d=OCN_DAstate_data(ivar)%posterior2d-OCN_DAstate_data(ivar)%prior2d
!               do i=da_nbs,da_nbe
!                   do ien=1,ensemble_size
!                      do k=1,DA_maxLevelCell(i)
!                           increment2d(i,k,ien)=sign(min(Mx_da_increm(ivar),abs(increment2d(i,k,ien))),increment2d(i,k,ien))
!                      enddo
!                   enddo
!               enddo
!              if(iamroot_da_all) write(logunit,*) trim(subname),trim(OCN_DAstate_data(ivar)%name),' increment max/loc',maxval(abs(increment2d)),maxloc(abs(increment2d))
!           !    call ensem_async_IAU3d(increment3d,OCN_DAstate_data(ivar)%IAU_increm)
!                OCN_DAstate_data(ivar)%IAU_increm2D=increment2d(:,:,my_inst) /steps_IAU
!           endif

        end do !ivar
        deallocate(increment2d,increment1d)
!        if(iamroot_da_all) write(logunit,*) trim(subname),'finished'

!EOC
    end subroutine ecda_calculate_IAU
!==============================================================================



!==============================================================================
!BOP
!   !ROUTINE: get_OCN_DA_comm
!   !DESCRIPTION:
!       set the DA grids for each PE, 
    subroutine get_OCN_DA_comm()
!EOP
!BOC
        implicit none
        character(*), parameter :: subName =   '(get_OCN_DA_comm) '
        integer :: npes_da_all,npes_da  ! for diagnostic
        integer :: ierr
!EOC
        ! --- Local variables ---

!-      get current comm for DA
        mpicom_da_all = seq_DA_comm_mpicom (DAOCNALLID)
        iam_da_all = seq_DA_comm_iam (DAOCNALLID)
        iamin_da_all = seq_DA_comm_iamin (DAOCNALLID)
        iamroot_da_all = seq_DA_comm_iamroot (DAOCNALLID)
        call mpi_comm_size(mpicom_da_all,npes_da_all,ierr)

        my_inst=seq_DA_comm_inst(DAOCNALLID)  ! use for DA ensemble ID

        mpicom_da = seq_DA_comm_mpicom (DAOCNID)
        iam_da = seq_DA_comm_iam (DAOCNID)
        iamin_da = seq_DA_comm_iamin (DAOCNID)
        iamroot_da = seq_DA_comm_iamroot (DAOCNID)
        call mpi_comm_size(mpicom_da,npes_da,ierr)

!       iamin_da  and iamin_da_all should be consistant

        if(iamroot_da_all) write(logunit,*)trim(subname),'npes :',npes_da_all,npes_da, &
                    'my peID:',iam_da_all,iam_da,'root:',iamroot_da_all,iamroot_da, &
                    'my_inst',my_inst
        flush(logunit)
     end subroutine get_OCN_DA_comm
!=============================================================================
!BOP
!   !ROUTINE: Define_OCN_DA_grids
!   !DESCRIPTION:
!       set the DA grids for each PE, 
    subroutine Define_OCN_DA_grids()
!EOP
!BOC
       implicit none
!EOC
      ! --- Local variables ---
      character(*), parameter :: subName =   '(Define_OCN_DA_grids) '
      integer, pointer :: nCellsSolve,nVertLevels,nEdgesSolve
      type (block_type), pointer :: block_ptr
      type (mpas_pool_type), pointer :: meshPool,verticalMeshPool
      real (kind=RKIND), dimension(:), pointer :: lonCell,latCell, lonEdge, latEdge, refZMid
      integer, dimension(:), pointer :: maxLevelCell,indexToCellID, indexToEdgeID
      integer :: i,ierr,icell,iedge
      real (kind=RKIND), dimension(:),allocatable :: local_lonCell, local_latCell, local_lonEdge, local_latEdge
      integer, dimension(:), allocatable:: local_maxLevelCell
      integer, dimension(:), allocatable::  local_indexToCellID,local_indexToEdgeID

      block_ptr => domain % blocklist

      !=get the verital layer depth: z-coordinator, [-5,-15,...]
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh',    meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh',    verticalMeshPool)
      call mpas_pool_get_array(verticalMeshPool, 'refZMid',refZMid)
      DA_nvert=nVertLevels   ! vertical layers
      allocate(DA_depth(DA_nvert),stat=ierr)    
      do i=1,DA_nvert
             DA_depth(i)=refZMid(i)
      enddo
      if(iamroot_da_all) write(logunit,*) trim(subname),'ocn depth',DA_depth(1:10)
      
!==     derive  DA Cell/Edge from model Cell/Edge at current PE  

      != get size of Cell and edge of current PE for allocate variables
      allCellsSolve=0
      allEdgesSolve=0
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
          
         !write(logunit,*) trim(subname),'blockID', block_ptr %blockID, block_ptr %localBlockID
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh',    meshPool)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

         allCellsSolve=allCellsSolve + nCellsSolve
         allEdgesSolve=allEdgesSolve + nEdgesSolve

         block_ptr => block_ptr % next
      end do
      allocate(local_lonCell(allCellsSolve),local_latCell(allCellsSolve),local_maxLevelCell(allCellsSolve),stat=ierr)
      allocate(local_lonEdge(allEdgesSolve),local_latEdge(allEdgesSolve),stat=ierr)
      !allocate(local_indexToCellID(allCellsSolve),local_indexToEdgeID(allEdgesSolve),stat=ierr)
      !write(logunit,*) trim(subname), 'Cells/Edge on current block/domian',nCellsSolve, nEdgesSolve,allCellsSolve,allEdgesSolve
      !write(logunit,*) trim(subname),'test local_lonCell',local_lonCell(1:4)
      !write(logunit,*) trim(subname),'test local_latCell',local_latCell(1:4)
      
      !=get the lontitude, latitude for each Cell, Edge
      icell = 0
      iedge = 0
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh',    meshPool)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_array(meshPool, 'lonCell',lonCell)
         call mpas_pool_get_array(meshPool, 'latCell',latCell)
         call mpas_pool_get_array(meshPool, 'lonEdge',lonEdge)
         call mpas_pool_get_array(meshPool, 'latEdge',latEdge)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
 
         !the output indexToCellID values are not as expected 
         !call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
         !call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
         !write(logunit,*) trim(subname), 'index cell/edge',indexToCellID(1),indexToEdgeID(1),maxLevelCell(3)

         !=local cell lon/lat for T/S,...
         do i = 1, nCellsSolve
             icell = icell + 1
             local_lonCell(icell)=lonCell(i)
             local_latCell(icell)=latCell(i)
             local_maxLevelCell(icell)=maxLevelCell(i)
             !local_indexToCellID(icell)=indexToCellID(i)

             !write(logunit,*) trim(subname), 'lon, i,icell',i,icell,lonCell(i),local_lonCell(icell)
             !write(logunit,*) trim(subname), 'lat, i,icell',i,icell,latCell(i),local_latCell(icell)
             !write(logunit,*) trim(subname), 'maxLevelCell, i,icell',i,icell,maxLevelCell(i),local_maxLevelCell(icell)
         end do
        !=local edge lon/lat for velocity DA
         do i = 1, nEdgesSolve
             iedge = iedge + 1
             local_lonEdge(iedge)=lonEdge(i)
             local_latEdge(iedge)=latEdge(i)
             !local_indexToEdgeID(iedge)=indexToEdgeID(i)
         end do

         block_ptr => block_ptr % next
      end do

      != DA local Cells/Edges for current PE   
      DACellsSolve = ceiling(allCellsSolve *1.0 /ensemble_size)
      DAEdgesSolve = ceiling(allEdgesSolve *1.0 /ensemble_size)
     ! Cells
      da_Cell_nts=(my_inst-1)*DACellsSolve+1
      da_Cell_nte=my_inst * DACellsSolve
      da_Cell_nte=min(da_cell_nte,allCellsSolve)
     ! Edges
      da_Edge_nts=(my_inst-1)*DAEdgesSolve+1
      da_Edge_nte=my_inst * DAEdgesSolve
      da_Edge_nte=min(da_Edge_nte,allEdgesSolve)

     ! @YL, no local decomposation for DA, redundant, use for EAKF first version
      da_Cell_nts=1
      da_Cell_nte=allCellsSolve
      da_Edge_nts=1
      da_Edge_nte=allEdgesSolve

      !  set halo as zero for now, Good for LETKF, to be develop for EAKF @YL
      da_Cell_nbs=da_Cell_nts
      da_Cell_nbe=da_Cell_nte
      da_Edge_nbs=da_Edge_nts
      da_Edge_nbe=da_Edge_nte

      allocate(DA_lonCell(da_Cell_nbs:da_Cell_nbe),DA_latCell(da_Cell_nbs:da_Cell_nbe),stat=ierr)
      allocate(DA_maxLevelCell(da_Cell_nbs:da_Cell_nbe),stat=ierr)
      allocate(DA_lonEdge(da_Edge_nbs:da_Edge_nbe),DA_latEdge(da_Edge_nbs:da_Edge_nbe),stat=ierr)

!      allocate(DA_indexToCellID(da_Cell_nbs:da_Cell_nbe),DA_indexToEdgeID(da_Edge_nbs:da_Edge_nbe),stat=ierr)
           
!   Set the DA local lon/lat
      DA_lonCell(da_Cell_nbs:da_Cell_nbe)=local_lonCell(da_Cell_nbs:da_Cell_nbe)
      DA_latCell(da_Cell_nbs:da_Cell_nbe)=local_latCell(da_Cell_nbs:da_Cell_nbe)
      DA_maxLevelCell(da_Cell_nbs:da_Cell_nbe)=local_maxLevelCell(da_Cell_nbs:da_Cell_nbe)
!      DA_indexToCellID(da_Cell_nbs:da_Cell_nbe)=local_indexToCellID(da_Cell_nbs:da_Cell_nbe)

      DA_lonEdge(da_Edge_nbs:da_Edge_nbe)=local_lonEdge(da_Edge_nbs:da_Edge_nbe)
      DA_latEdge(da_Edge_nbs:da_Edge_nbe)=local_latEdge(da_Edge_nbs:da_Edge_nbe)
!      DA_indexToEdgeID(da_Edge_nbs:da_Edge_nbe)=local_indexToEdgeID(da_Edge_nbs:da_Edge_nbe)

!      write(logunit,*) trim(subname), 'Cells on current PE',allCellsSolve,DACellsSolve
!      write(logunit,*) trim(subname), 'DA local grid',my_inst,da_Cell_nbs,da_Cell_nbe,da_Edge_nbs,da_Edge_nbe
      if(iam_local .eq. 0) write(logunit,*) trim(subname), 'DA inst,cell/edge min/max lon/lat',my_inst, &
                        minval(DA_lonCell),maxval( DA_lonCell),minval(DA_latCell),maxval( DA_latCell)
!      if(iamroot_da) write(logunit,*) trim(subname), 'local DA lon/lat', DA_lonCell(1:4),DA_latCell(1:4)
!      if(iamroot_da) write(logunit,*) trim(subname), 'local PE lon/lat', lonCell(1:4),latCell(1:4)

!      write(logunit,*) trim(subname), 'DA index cell/edge',DA_indexToCellID(da_Cell_nbs),DA_indexToEdgeID(da_Edge_nbs)
      
      
      deallocate(local_lonCell,local_latCell,local_lonEdge,local_latEdge)
          
      
     end subroutine Define_OCN_DA_grids


!==============================================================================
!BOP
!   !ROUTINE: is_ocn_da_time
!
!   !DESCRIPTION:
!       judge if mode reach to DA time for ocean and return the observation record
!
!   !REVISION HISTORY:
!       Feb 5, 2019 - Y. Liu <liu6@tamu.edu> - initial release
!   !INTERFACE:
    subroutine is_ocn_da_time(isdatime,obs_rec,o_clock)
!EOP
!BOC
      implicit none
      Logical, intent(out) :: isdatime
      integer, intent(out) :: obs_rec
      type (MPAS_Clock_type), intent(in) :: o_clock
!EOC
       !local variables
      character(*), parameter :: subName =   '(is_ocn_da_time) '
      type (MPAS_Time_Type) :: currTime
      integer :: err_tmp, err
      integer :: ymd, tod, ihour, iminute, isecond
      integer :: iyear, imonth, iday, curr_ymd, curr_tod
      integer :: year,month,day,hour,minute,second, doy
      isdatime = .FALSE.
      obs_rec=0   ! the record for grids data

      err = 0
      currTime = mpas_get_clock_time(o_clock, MPAS_NOW, err_tmp)
      err = ior(err,err_tmp)

      !if(currTime < obs_init_Time) return       ! model time before observation initial time, @YL not working

      call mpas_get_time(curr_time=currTime, YYYY=year, MM=month, DD=day,DoY=doy, H=hour, M=minute, S=second, ierr=err_tmp)
      err = ior(err,err_tmp)
      
      ! DA only in whole hours now
      select case (DA_interv_unit)
          case ('mo')
            if(hour /=0 .or. second /= 0 .or. minute /= 0) return
             imonth = (year-init_obs_YR)*12 + (month-init_obs_MO)
             if ( day /=init_obs_DY  .or. mod(imonth,DA_interv) /=0 ) RETURN     ! the DA interval in month, do DA at the middle of the month (ex: 1/16,2/16,3/16),00z
             imonth = (year-init_obs_YR)*12 + (month-init_obs_MO)
             if( mod(imonth,obs_interv) == 0) obs_rec=imonth/obs_interv+1                    ! add 1 is the first record are the obs initial time

          case ('dy')
             if(second /= 0 .or. minute /= 0) return
             iday= (year-init_obs_YR)*365+(doy-init_obs_doy)
             if ( hour /= init_obs_HR .or. mod(iday,DA_interv) /=0 ) RETURN     ! the DA interval in day, do DA at the middle of the day depend on input, ex: 12z  @YL
             !obs_rec=idays/DA_interv+1                    ! add 1 is the first record are the obs initial time
             if( mod(iday,obs_interv) ==0)obs_rec=iday/obs_interv+1                    ! add 1 is the first record are the obs initial time

          case ('hr')                                    ! the 24/DA_interv is whole number
             if(second /= 0 .or. minute /= 0) return
             if(mod(hour-init_obs_HR,DA_interv) /=0 ) RETURN     ! do DA at begining of the hour @YL
             iday= (year-init_obs_YR)*365+(doy-init_obs_doy)   ! ignore the leap year, revisit @Yun ???
             ihour =iday*24+hour-init_obs_HR        ! hours from current model time to observation initial time.
             if( mod(ihour,obs_interv) ==0) obs_rec= ihour/obs_interv+1

          CASE DEFAULT
             write(logunit,*) trim(subname),'data assimilation time interval unit error, ask for hr,dy,mo, but get ', DA_interv_unit
             return

        end select
     
       ! some profile data save in the file of XX.YMDH 
        write(YMDH,'(I4.4,I2.2,I2.2,I2.2)')year,month,day,hour
        if(iamroot_da_all) write(logunit,*) trim(subname),'YMDH/obs_rec', YMDH,obs_rec
        if(obs_rec>0) isdatime = .TRUE.    ! here the obs_rec is for all the grid obs
        RETURN
  

      ! @YL assimilate on 12:00:00
   !   if(hour == 12 .and. iminute == 0 .and. isecond == 0) isdatime = .TRUE.

     end subroutine is_ocn_da_time

!==============================================================================
!BOP
!   !ROUTINE: read_ocn_obs4da
!   ! read oceanic observations into OCN_obs_data
    subroutine read_ocn_obs4da(obs_rec,obs_exists)
!BOC
       implicit none
       integer, intent(in) :: obs_rec
       logical,intent(out) :: obs_exists

       character(*), parameter :: subName =   '(read_ocn_obs4da) '
       integer :: i,j,k, ID
       logical :: exists

       obs_exists=.false.
       do ID=1,NUM_OCN_obs_type
          if (.not. OCN_obs_data(ID)%assimilate) cycle
          SELECT CASE (OCN_obs_data(ID)%grid_type)
            CASE(O_grid_model)
               call read_mpaso_obs_profile(ID,obs_rec,exists)
            CASE(O_profile)
               call read_mpaso_obs_profile(ID,obs_rec,exists)
            CASE DEFAULT
               write(logunit,*) trim(subname), 'to be developed'
          end SELECT
          if(exists)obs_exists=.true.
       enddo

    end subroutine read_ocn_obs4da
!==============================================================================

!BOP
!   !ROUTINE: read_mpaso_obs_profile
!   ! read observation in profile ( in records)
!   ! the time dimension is on the filename

    subroutine read_mpaso_obs_profile(ID,obs_rec,exists)
!EOP
     implicit none
     integer,intent(in) :: ID
     integer,intent(in) :: obs_rec

     logical,intent(out) :: exists
!BOC
     character(*), parameter :: subName = '(read_mpaso_obs_profile)'

     type (MPAS_IO_Handle_type) :: inputFile
     type (MPAS_Stream_type) :: observationStream
     type (mpas_io_context_type), pointer :: iocontext_ptr
     type (field1DReal) :: OB_lon,OB_lat,OB_depth,OB_obs, OB_err
     type (field1DInteger) :: OB_index 
     
     real (kind=RKIND) :: mxlon,mnlon,mxlat,mnlat   
     integer :: nDIMOBS, nobs,iob, iErr,i,j,k
     integer :: blockID
     integer, dimension(:),allocatable :: obs_index 
   
     !      check if the observation exist at current time YMDH
     if(OCN_obs_data(ID)%grid_type == O_profile) then 
           inquire(file=trim(OCN_obs_data(ID)%filename)//YMDH//'.nc', exist=exists)
     else
           inquire(file=trim(OCN_obs_data(ID)%filename), exist=exists)

     endif
     if (.not. exists) return         !no observation file exist

     ! use for decomp OBS
     blockID = domain % blocklist %blockID    ! here I assume only one block for each PE. @YL
     
     !character (len=StrKIND), pointer ::   config_ocean_observation_file, &
     !                                      config_ocean_observation_lat_varname, &
     !                                      config_ocean_observation_nlat_dimname, &
     !                                      config_ocean_observation_lon_varname, &
     !                                      config_ocean_observation_nlon_dimname, &
     !                                      config_ocean_observation_depth_varname, &
     !                                      config_ocean_observation_ndepth_dimname, &
     !                                      config_ocean_observation_obs_varname, &
     !                                      config_ocean_observation_nobs_dimname, &
     !                                      config_ocean_observation_err_varname, &
     !                                      config_ocean_observation_nerr_dimname

     !  the  MPAS_readStream(observationStream, obs_rec, iErr) can not identify the correct obs_rec, it always read tghe first
     !  record. SO need fugure out the issues @YL
       ! get observation dimension
       iocontext_ptr  =>  domain % iocontext
       if(OCN_obs_data(ID)%grid_type == O_profile) then 
           inputFile = MPAS_io_open(trim(OCN_obs_data(ID)%filename)//YMDH//'.nc', MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)
       else
           inputFile = MPAS_io_open(trim(OCN_obs_data(ID)%filename), MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)
       endif
       call MPAS_io_inq_dim(inputFile, trim(OCN_obs_data(ID)%Dm_name), nDIMOBS, iErr)
       call MPAS_io_close(inputFile, iErr)

       if(OCN_obs_data(ID)%grid_type == O_profile) then 
            call MPAS_createStream(observationStream, domain % iocontext, trim(OCN_obs_data(ID)%filename)//YMDH//'.nc', MPAS_IO_NETCDF, &
                              MPAS_IO_READ, ierr=iErr)
       else
            call MPAS_createStream(observationStream, domain % iocontext, trim(OCN_obs_data(ID)%filename), MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)
       endif
       ! the same dimesion name for, lon,lat,depth,obs,obs error

      ! Setup OB_lon,OB_lat,OB_depth,OB_obs, OB_err to be read in 
       OB_lat % fieldName = trim(OCN_obs_data(ID)%latname)
       OB_lat % dimSizes(1) = nDIMOBS
       OB_lat % dimNames(1) = trim(OCN_obs_data(ID)%Dm_name)
       OB_lat % isVarArray = .false.
       OB_lat % isPersistent = .true.
       OB_lat % isActive = .true.
       OB_lat % hasTimeDimension = .false.
       OB_lat % block => domain % blocklist
       allocate(OB_lat % attLists(1))
       allocate(OB_lat % array(nDIMOBS))
       call MPAS_streamAddField(observationStream, OB_Lat, iErr)

       OB_lon % fieldName = trim(OCN_obs_data(ID)%lonname)
       OB_lon % dimSizes(1) = nDIMOBS
       OB_lon % dimNames(1) = trim(OCN_obs_data(ID)%Dm_name)
       OB_lon % isVarArray = .false.
       OB_lon % isPersistent = .true.
       OB_lon % isActive = .true.
       OB_lon % hasTimeDimension = .false.
       OB_lon % block => domain % blocklist
       allocate(OB_lon % attLists(1))
       allocate(OB_lon % array(nDIMOBS))
       call MPAS_streamAddField(observationStream, OB_Lon, iErr)

       OB_depth % fieldName = trim(OCN_obs_data(ID)%vertical_name)
       OB_depth % dimSizes(1) = nDIMOBS
       OB_depth % dimNames(1) = trim(OCN_obs_data(ID)%Dm_name)
       OB_depth % isVarArray = .false.
       OB_depth % isPersistent = .true.
       OB_depth % isActive = .true.
       OB_depth % hasTimeDimension = .false.
       OB_depth % block => domain % blocklist
       allocate(OB_depth % attLists(1))
       allocate(OB_depth % array(nDIMOBS))
       call MPAS_streamAddField(observationStream, OB_depth, iErr)

       OB_obs % fieldName = trim(OCN_obs_data(ID)%varname)
       OB_obs % dimSizes(1) = nDIMOBS
       OB_obs % dimNames(1) = trim(OCN_obs_data(ID)%Dm_name)
       OB_obs % isVarArray = .false.
       OB_obs % isPersistent = .true.
       OB_obs % isActive = .true.
       OB_obs % hasTimeDimension = .true.
       OB_obs % block => domain % blocklist
       allocate(OB_obs % attLists(1))
       allocate(OB_obs % array(nDIMOBS))
       call MPAS_streamAddField(observationStream, OB_obs, iErr)

       if(OCN_obs_data(ID)%err_var0 <=0.0) then
          OB_err % fieldName = trim(OCN_obs_data(ID)%uncert_name)
          OB_err % dimSizes(1) = nDIMOBS
          OB_err % dimNames(1) = trim(OCN_obs_data(ID)%Dm_name)
          OB_err % isVarArray = .false.
          OB_err % isPersistent = .true.
          OB_err % isActive = .true.
          OB_err % hasTimeDimension = .true.
          OB_err % block => domain % blocklist
          allocate(OB_err % attLists(1))
          allocate(OB_err % array(nDIMOBS))
          call MPAS_streamAddField(observationStream, OB_err, iErr)
       end if

     ! use for decomp OBS
       OB_index % fieldName = 'indexToblock'
       OB_index % dimSizes(1) = nDIMOBS
       OB_index % dimNames(1) = trim(OCN_obs_data(ID)%Dm_name)
       OB_index % isVarArray = .false.
       OB_index % isPersistent = .true.
       OB_index % isActive = .true.
       OB_index % hasTimeDimension = .false.
       OB_index % block => domain % blocklist
       allocate(OB_index % attLists(1))
       allocate(OB_index % array(nDIMOBS))
       call MPAS_streamAddField(observationStream, OB_index, iErr)

       ! Add  OB_lon,OB_lat,OB_obs 

       ! Read stream
       if(OCN_obs_data(ID)%grid_type == O_profile) then 
             call MPAS_readStream(observationStream, 1, iErr)
       else
             call MPAS_readStream(observationStream, obs_rec, iErr)
       endif
       ! Close stream
       call MPAS_closeStream(observationStream)      ! If need destroy the Strem for next observations? @YL

       ! output to check the value
        if(iamroot_da_all) then
           write(logunit,*) trim(subname),'obs input begin',OB_lat%array(1),OB_lon%array(1),OB_obs%array(1) 
           write(logunit,*) trim(subname),'obs input end',OB_lat%array(nDIMOBS),OB_lon%array(nDIMOBS),OB_obs%array(nDIMOBS) 
        end if
!      collect local obs 
       ! can be move to grid setting @YL
       ! periodic boundary to be take care @YL
!        mxlon=max(maxval(DA_lonCell),maxval(DA_lonEdge))
!        mnlon=min(minval(DA_lonCell),minval(DA_lonEdge))
!        mxlat=max(maxval(DA_latCell),maxval(DA_latEdge))
!        mnlat=min(minval(DA_latCell),minval(DA_latEdge))
!        write(logunit,*) trim(subname),'local DA grid lon,lat range',mnlon,mxlon,mxlat,mnlat,'@ curent PE:',iam_da_all,iam_da
        
        ! get local obs number
        allocate(obs_index(nDIMOBS),stat=ierr)
        nobs=0
        do i=1,nDIMOBS
!           if(OB_lat%array(i)<mnlat .or. OB_lat%array(i)>mxlat .or. OB_lon%array(i)<mnlon .or. OB_lon%array(i)> mxlon) CYCLE
            if(OB_index%array(i) .NE. blockID) CYCLE
            nobs=nobs+1
            obs_index(nobs)=i
        enddo
        !-     set the observation size  based on observation dimension 
        OCN_obs_data(ID)%nobs=nobs
       ! write(logunit,*) trim(subname),'nobs:',nobs,'@ curent PE:',iam_da_all,iam_da
        if(nobs==0)then
           write(logunit,*) trim(subname),'finish read obs:',trim(OCN_obs_data(ID)%obs_name),OCN_obs_data(ID)%nobs
            return
        endif

        allocate (OCN_obs_data(ID)%values(nobs), stat=ierr)
        allocate (OCN_obs_data(ID)%lon(nobs), stat=ierr)
        allocate (OCN_obs_data(ID)%lat(nobs), stat=ierr)
        allocate (OCN_obs_data(ID)%depth(nobs), stat=ierr)
        allocate (OCN_obs_data(ID)%err_var(nobs), stat=ierr)

        !save the local observations     
        do iob=1,nobs
            OCN_obs_data(ID)%lon(iob) = OB_lon%array(obs_index(iob))
            OCN_obs_data(ID)%lat(iob) = OB_lat%array(obs_index(iob)) 
            OCN_obs_data(ID)%depth(iob) = OB_depth%array(obs_index(iob))
            OCN_obs_data(ID)%values(iob)= OB_obs%array(obs_index(iob))
           ! OCN_obs_data(ID)%err_var(iob)=OB_err%array(obs_index(iob))
       end do

       if(OCN_obs_data(ID)%err_var0 <=0.0) then
            do iob=1,nobs
              OCN_obs_data(ID)%err_var(iob)=OB_err%array(obs_index(iob))
            end do
       else 
            OCN_obs_data(ID)%err_var= OCN_obs_data(ID)%err_var0
       endif 
       
      deallocate(obs_index) 
        !
!     write(logunit,*) trim(subname),'lon range: obs/state',maxval(OCN_obs_data(ID)%lon),minval(OCN_obs_data(ID)%lon),&
!              maxval(DA_lonCell),minval(DA_lonCell)
!     write(logunit,*) trim(subname),'lat range: obs/state',maxval(OCN_obs_data(ID)%lat),minval(OCN_obs_data(ID)%lat),&
!              maxval(DA_latCell),minval(DA_latCell)

!     write(logunit,*) trim(subname),'finished:',trim(OCN_obs_data(ID)%obs_name),OCN_obs_data(ID)%nobs, &
!               '@ curent PE:',iam_da_all,iam_da

       call mpi_barrier(mpicom_da_all, ierr)
       !

    end subroutine read_mpaso_obs_profile


!==============================================================================
!BOP
!   !ROUTINE: collect_mpaso_state4da
!   ! collect DA required ocean state from each model member to creat ensmeble array.
    subroutine collect_mpaso_state4da()

!BOC
       implicit none
       character(*), parameter :: subName =   '(collect_mpaso_state4da) '
       integer :: i,j,k, ivar,ierr

       call mpaso_DAstate_allocate()
       do ivar = 1, n_OCN_state
           if(.not. OCN_DAstate_data(ivar)%required) CYCLE
           !write(logunit,*) trim(subname),ivar, OCN_DAstate_data(ivar)%name,OCN_DAstate_data(ivar)%state_dim
           select case (ivar)
                case(OCN_TEMP_index)
                        call ensem_async_prior_tracers(ivar,'index_temperature')
                case(OCN_SALT_index)
                        call ensem_async_prior_tracers(ivar,'index_salinity')
                case(OCN_VEL_index)
                        call ensem_async_prior_states(ivar,'normalVelocity')
                case(OCN_SSH_index)
                        call ensem_async_prior_states(ivar,'ssh')
                CASE DEFAULT
                       write(logunit,*)trim(subname),'More states to develop'
           end select
           
           if(iamroot_da_all .and. OCN_DAstate_data(ivar)%state_dim == 1) &
           write(logunit,*)trim(subname),trim(OCN_DAstate_data(ivar)%name), &
                   'max/min',maxval(OCN_DAstate_data(ivar)%prior1d),minval(OCN_DAstate_data(ivar)%prior1d)
           if(iamroot_da_all .and. OCN_DAstate_data(ivar)%state_dim == 2) &
           write(logunit,*)trim(subname),trim(OCN_DAstate_data(ivar)%name), &
                   'max/min',maxval(OCN_DAstate_data(ivar)%prior2d),minval(OCN_DAstate_data(ivar)%prior2d)
       end do
    end subroutine collect_mpaso_state4da
!==============================================================================
!BOP
!   !ROUTINE: ensem_async_prior_tracers
!   ! collect DA required ocean tracers from each model member to creat ensmeble array.
!   ! mpaso treat temperature and salinity as tracers 
!   ! All tracers are 2-D, DA will be 3D
    subroutine ensem_async_prior_tracers(VID, idx_name)
!BOC
       implicit none
       integer,intent(in) :: VID
       character(*),intent(in) :: idx_name
    
!BOC
       character(*), parameter :: subName =   '(ensem_async_prior_tracers) '
      
!     temporal for testing communication
      real (kind=RKIND), dimension(:,:,:),allocatable :: state_da,state_da_all

!      type (mpas_pool_type), pointer ::  diagnosticsPool, forcingPool, averagePool, scratchPool
      type (mpas_pool_type), pointer :: meshPool, statePool, tracersPool

      type (block_type), pointer :: block_ptr    
      integer, pointer :: nCellsSolve,nVertLevels
      integer, pointer :: idx
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: SSHValue
      real (kind=RKIND), dimension(:,:), pointer :: TValue
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      integer :: da_nbs,da_nbe
      integer:: i,j,k,n
      integer:: ierr

 !    allocate local  array, tracers are 2D, so ensemble are 3D, @YL no halo used
      allocate (state_da(allCellsSolve,DA_nvert,ensemble_size),state_da_all(allCellsSolve,DA_nvert,ensemble_size),stat=ierr)
      state_da=0.0
      state_da_all=0.0

 !    collect state from mesh
      n = 0
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
         call mpas_pool_get_dimension(tracersPool, trim(idx_name), idx)

         !if (iamroot_da) write(logunit,*) 'vertlevels',nVertLevels, maxLevelCell(1), maxLevelCell(nCellsSolve)
         do i = 1, nCellsSolve
             n = n + 1
             do k=1,maxLevelCell(i)
                  state_da(n,k,my_inst)=activeTracers(idx,k,i)
             end do
         end do

          block_ptr => block_ptr % next
      end do       

     !gather ensemble with DA local comm
      call MPI_ALLREDUCE(state_da,state_da_all,size(state_da),MPI_REALKIND,MPI_SUM,mpicom_da, ierr)
     ! keep the local DA ensemble for current PE
      OCN_DAstate_data(VID)%prior2d(da_Cell_nbs:da_Cell_nbe,:,:)=state_da_all(da_Cell_nbs:da_Cell_nbe,:,:)
      OCN_DAstate_data(VID)%posterior2d=OCN_DAstate_data(VID)%prior2d

!      if (iamroot_da) write(logunit,*) &
!            trim(subname), 'PEIDs,array',my_inst, iam_local, state_da(da_Cell_nbs,1,my_inst),state_da_all(da_Cell_nbs,1,my_inst), &
!                       OCN_DAstate_data(VID)%prior2d(da_Cell_nbs,1,:)
      deallocate(state_da,state_da_all)
      call mpi_barrier(mpicom_da_all, ierr)

    end subroutine ensem_async_prior_tracers

!==============================================================================
!BOP
!   !ROUTINE: ensem_async_prior_states
!   ! collect DA required ocean tracers from each model member to creat ensmeble array.
!   ! mpaso treat temperature and salinity as tracers 
!   ! All tracers are 2-D, DA will be 3D
    subroutine ensem_async_prior_states(VID, var_name)
!BOC
       implicit none
       integer,intent(in) :: VID
       character(*),intent(in) :: var_name
    
!BOC
       character(*), parameter :: subName =   '(ensem_async_prior_states) '
      
!     temporal for testing communication
      real (kind=RKIND), dimension(:,:,:),allocatable :: state_da,state_da_all
      real (kind=RKIND), dimension(:,:),allocatable :: state_1D_da,state_1D_da_all

!      type (mpas_pool_type), pointer ::  diagnosticsPool, forcingPool, averagePool, scratchPool
      type (mpas_pool_type), pointer :: meshPool, statePool, tracersPool

      type (block_type), pointer :: block_ptr    
      integer, pointer :: nCellsSolve,nVertLevels
      integer, pointer :: idx
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: state1DValue
      real (kind=RKIND), dimension(:,:), pointer :: state2DValue

      integer :: da_nbs,da_nbe, allSolve
      integer:: i,j,k,n,ierr

      ! there are three type of grids: CEll, Edge, Vertex. Start from CELL
      allSolve=allCellsSolve
      da_nbs=da_Cell_nbs
      da_nbe=da_Cell_nbe

 !    allocate local  array, tracers are 1/2D, so ensemble are 2/3D, @YL no halo used
     if(OCN_DAstate_data(VID)%state_dim ==2) then 
          allocate (state_da(allSolve,DA_nvert,ensemble_size),state_da_all(allSolve,DA_nvert,ensemble_size),stat=ierr)
          state_da=0.0
          state_da_all=0.0
      end if 
      if(OCN_DAstate_data(VID)%state_dim ==1) then
          allocate (state_1D_da(allSolve,ensemble_size),state_1D_da_all(allSolve,ensemble_size),stat=ierr)
          state_1D_da=0.0
          state_1D_da_all=0.0
      end if 

 !    collect state from mesh
      n = 0
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)      !@YL
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)    !@YL
         if(OCN_DAstate_data(VID)%state_dim ==1)  call mpas_pool_get_array(statePool, trim(var_name), state1DValue,1)
         if(OCN_DAstate_data(VID)%state_dim ==2)  call mpas_pool_get_array(statePool, trim(var_name), state2DValue,1)

         !if (iamroot_da) write(logunit,*) 'vertlevels',nVertLevels, maxLevelCell(1), maxLevelCell(nCellsSolve)

         do i = 1, nCellsSolve
             n = n + 1
             if(OCN_DAstate_data(VID)%state_dim ==1) state_1D_da(n,my_inst)=state1DValue(i)
             
             if(OCN_DAstate_data(VID)%state_dim ==2)then
                do k=1,maxLevelCell(i)
                   state_da(n,k,my_inst)=state2DValue(k,i)
                enddo
             endif
         end do
          block_ptr => block_ptr % next
      end do       

     !gather ensemble with DA local comm
      if (OCN_DAstate_data(VID)%state_dim ==2) then
          call MPI_ALLREDUCE(state_da,state_da_all,size(state_da),MPI_REALKIND,MPI_SUM,mpicom_da, ierr)
          ! keep the local DA ensemble for current PE
          OCN_DAstate_data(VID)%prior2d(da_nbs:da_nbe,:,:)=state_da_all(da_nbs:da_nbe,:,:)
          OCN_DAstate_data(VID)%posterior2d=OCN_DAstate_data(VID)%prior2d
          ! if (iamroot_da) write(logunit,*) &
          !       trim(subname), 'PEIDs,array',my_inst, iam_local, state_da(da_nbs,1,my_inst),state_da_all(da_nbs,1,my_inst), &
          !              OCN_DAstate_data(VID)%prior2d(da_nbs,1,:)
          deallocate(state_da,state_da_all)
      endif

      if (OCN_DAstate_data(VID)%state_dim ==1) then
          call MPI_ALLREDUCE(state_1D_da,state_1D_da_all,size(state_1D_da),MPI_REALKIND,MPI_SUM,mpicom_da, ierr)
          ! keep the local DA ensemble for current PE
          OCN_DAstate_data(VID)%prior1d(da_nbs:da_nbe,:)=state_1D_da_all(da_nbs:da_nbe,:)
          OCN_DAstate_data(VID)%posterior1d=OCN_DAstate_data(VID)%prior1d
          !if (iamroot_da) write(logunit,*) &
          !    trim(subname), 'PEIDs,array',my_inst, iam_local, state_1D_da(da_nbs,my_inst),state_1D_da_all(da_nbs,my_inst), &
          !             OCN_DAstate_data(VID)%prior1d(da_nbs,:)
          deallocate( state_1D_da,state_1D_da_all)
      endif
       call mpi_barrier(mpicom_da_all, ierr)

    end subroutine ensem_async_prior_states

!==============================================================================




end module mpaso_DA_mod
