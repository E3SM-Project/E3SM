module glc_comp_mct

! !USES:

  use shr_sys_mod
  use shr_kind_mod     , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, &
                               CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod     , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel, &
                               shr_file_setlogunit, shr_file_setloglevel, shr_file_setio, &
                               shr_file_freeunit
  use mct_mod
  use esmf

  use seq_flds_mod
  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod

  use glc_cpl_indices
  use glc_constants,       only: verbose, stdout, stderr, nml_in, &
                                 radius,  radian, tkfrz,  glc_nec
  use glc_errormod,        only: glc_success
  use glc_InitMod,         only: glc_initialize
  use glc_RunMod,          only: glc_run
  use glc_FinalMod,        only: glc_final
  use glc_io,              only: glc_io_write_restart, glc_io_write_history
  use glc_communicate,     only: init_communicate, my_task, master_task
  use glc_time_management, only: iyear,imonth,iday,ihour,iminute,isecond,runtype
  use glc_global_fields,   only: ice_sheet
  use glc_global_grid,     only: glc_grid, glc_landmask, glc_landfrac

! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: glc_init_mct
  public :: glc_run_mct
  public :: glc_final_mct

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

  !--- stdin input stuff ---
  character(CS) :: str                  ! cpp  defined model name
  
  !--- other ---
  integer(IN)   :: errorcode            ! glc error code
  
  ! my_task_local and master_task_local are needed for some checks that are done before
  ! init_communicate is called (although, it's possible that init_communicate could be
  ! moved to earlier to prevent the need for these copies)
  integer(IN)           :: my_task_local        ! my task in mpi communicator mpicom 
  integer(IN),parameter :: master_task_local=0  ! task number of master task
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: glc_init_mct
!
! !DESCRIPTION:
!     initialize glc model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine glc_init_mct( EClock, cdata, x2g, g2x, NLFilename )

! !USES:

    use glc_ensemble, only : set_inst_vars, write_inst_vars, get_inst_name
    use glc_files   , only : set_filenames, ionml_filename

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)         , intent(in)    :: EClock
    type(seq_cdata)          , intent(inout) :: cdata
    type(mct_aVect)          , intent(inout) :: x2g, g2x
    character(len=*), optional  , intent(in) :: NLFilename ! Namelist filename

!EOP

    !--- local variables ---
    integer(IN)              :: ierr        ! error code
    integer(IN)              :: i,j,n,nxg,nyg
    integer(IN)              :: COMPID
    integer(IN)              :: mpicom
    integer(IN)              :: lsize
    type(mct_gsMap), pointer :: gsMap
    type(mct_gGrid), pointer :: dom
    type(seq_infodata_type), pointer :: infodata   ! Input init object
    integer(IN)              :: shrlogunit, shrloglev  
    character(CL)            :: starttype
    character(CS)            :: myModelName

    !--- formats ---
    character(*), parameter :: F00   = "('(glc_init_mct) ',8a)"
    character(*), parameter :: F01   = "('(glc_init_mct) ',a,8i8)"
    character(*), parameter :: F02   = "('(glc_init_mct) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(glc_init_mct) ',a,i8,a)"
    character(*), parameter :: F90   = "('(glc_init_mct) ',73('='))"
    character(*), parameter :: F91   = "('(glc_init_mct) ',73('-'))"
    character(*), parameter :: subName = "(glc_init_mct) "
!EOP
!-------------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    ! Determine attribute vector indices
    !----------------------------------------------------------------------------
    
    call glc_cpl_indices_set()

    !----------------------------------------------------------------------------
    ! Set cdata pointers
    !----------------------------------------------------------------------------

    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsMap, dom=dom, infodata=infodata)

    call mpi_comm_rank(mpicom, my_task_local, ierr)

    !---------------------------------------------------------------------------
    ! set variables that depend on ensemble index
    !---------------------------------------------------------------------------

    call set_inst_vars(COMPID)
    call get_inst_name(myModelName)
    call set_filenames()

    !---------------------------------------------------------------------------
    ! use infodata to determine type of run
    !---------------------------------------------------------------------------

    call seq_infodata_GetData( infodata, start_type=starttype)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       write(*,*) 'glc_comp_mct ERROR: unknown starttype'
       call shr_sys_abort()
    end if

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    !--- open log file ---
    if (my_task_local == master_task_local) then
       stdout = shr_file_getUnit()
       call shr_file_setIO(ionml_filename,stdout)
    else
       stdout = 6
    endif
    stderr = stdout
    nml_in = shr_file_getUnit()

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout)

    errorCode = glc_Success
    if (verbose .and. my_task_local == master_task_local) then
       write(stdout,F00) ' Starting'
       write(stdout,*) subname, 'COMPID: ', COMPID
       call write_inst_vars
       call shr_sys_flush(stdout)
    endif
    call init_communicate(mpicom)
    call glc_initialize(errorCode)
    if (verbose .and. my_task == master_task) then
       write(stdout,F01) ' GLC Initial Date ',iyear,imonth,iday,ihour,iminute,isecond
       write(stdout,F01) ' Initialize Done', errorCode
       call shr_sys_flush(stdout)
    endif

    nxg = glc_grid%nx
    nyg = glc_grid%ny
    lsize = nxg*nyg

    ! Initialize MCT gsmap

    call glc_SetgsMap_mct(mpicom, COMPID, gsMap)

    ! Initialize MCT domain

    call glc_domain_mct(lsize,gsMap,dom)

    ! Set flags in infodata

    call seq_infodata_PutData(infodata, glc_present=.true., &
       glc_prognostic = .true., glc_nx=nxg, glc_ny=nyg)

    ! Initialize MCT attribute vectors

    call mct_aVect_init(g2x, rList=seq_flds_g2x_fields, lsize=lsize)
    call mct_aVect_zero(g2x)

    call mct_aVect_init(x2g, rList=seq_flds_x2g_fields, lsize=lsize)
    call mct_aVect_zero(x2g)

   if (my_task == master_task) then
      write(stdout,F91) 
      write(stdout,F00) trim(myModelName),': start of main integration loop'
      write(stdout,F91) 
   end if

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(stdout)

end subroutine glc_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: glc_run_mct
!
! !DESCRIPTION:
!     run method for glc model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine glc_run_mct( EClock, cdata, x2g, g2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2g        ! driver -> glc
   type(mct_aVect)             ,intent(inout) :: g2x        ! glc    -> driver
   
!EOP
   !--- local ---
   integer(IN)   :: cesmYMD           ! cesm model date
   integer(IN)   :: cesmTOD           ! cesm model sec
   integer(IN)   :: glcYMD            ! glc model date
   integer(IN)   :: glcTOD            ! glc model sec 
   integer(IN)   :: n                 ! index
   integer(IN)   :: nf                ! fields loop index
   integer(IN)   :: ki                ! index of ifrac
   integer(IN)   :: lsize             ! size of AttrVect
   integer(IN)   :: nxg,nyg           ! global grid size
   real(R8)      :: lat               ! latitude
   real(R8)      :: lon               ! longitude
   integer(IN)   :: shrlogunit, shrloglev  
   logical       :: stop_alarm        ! is it time to stop
   logical       :: hist_alarm        ! is it time to write a history file
   logical       :: rest_alarm        ! is it time to write a restart
   logical       :: done              ! time loop logical
   integer           :: num 
   character(len= 2) :: cnum
   character(len=64) :: name

   character(*), parameter :: F00   = "('(glc_run_mct) ',8a)"
   character(*), parameter :: F01   = "('(glc_run_mct) ',a,8i8)"
   character(*), parameter :: F04   = "('(glc_run_mct) ',2a,2i8,'s')"
   character(*), parameter :: subName = "(glc_run_mct) "
!-------------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout)

    lsize = mct_avect_lsize(x2g)
    nxg = glc_grid%nx
    nyg = glc_grid%ny

    ! Set internal time info
 
    errorCode = glc_Success
    call seq_timemgr_EClockGetData(EClock,curr_ymd=cesmYMD, curr_tod=cesmTOD)
    stop_alarm = seq_timemgr_StopAlarmIsOn( EClock )

    glcYMD = iyear*10000 + imonth*100 + iday
    glcTOD = ihour*3600 + iminute*60 + isecond
    done = .false.
    if (glcYMD == cesmYMD .and. glcTOD == cesmTOD) done = .true.
    if (verbose .and. my_task == master_task) then
       write(stdout,F01) ' Run Starting ',glcYMD,glcTOD
       call shr_sys_flush(stdout)
    endif

    ! Unpack

    do num = 1,glc_nec
       call glc_import_mct(x2g, num, &
            index_x2g_Ss_tsrf(num), index_x2g_Ss_topo(num), index_x2g_Fgss_qice(num)) 
    end do

    ! Run 

    do while (.not. done) 
       if (glcYMD > cesmYMD .or. (glcYMD == cesmYMD .and. glcTOD > cesmTOD)) then
          write(stdout,*) subname,' ERROR overshot coupling time ',glcYMD,glcTOD,cesmYMD,cesmTOD
          call shr_sys_abort('glc error overshot time')
       endif

       call glc_run

       glcYMD = iyear*10000 + imonth*100 + iday
       glcTOD = ihour*3600 + iminute*60 + isecond
       if (glcYMD == cesmYMD .and. glcTOD == cesmTOD) done = .true.
       if (verbose .and. my_task == master_task) then
          write(stdout,F01) ' GLC  Date ',glcYMD,glcTOD
       endif
    enddo

    if (verbose .and. my_task == master_task) then
       write(stdout,F01) ' Run Done',glcYMD,glcTOD
       call shr_sys_flush(stdout)
    endif
    
    ! Pack

    do num = 1,glc_nec
       call glc_export_mct(g2x, num, &
            index_g2x_Sg_frac(num), index_g2x_Sg_topo(num)  ,&
            index_g2x_Fsgg_rofi(num), index_g2x_Fsgg_rofl(num),index_g2x_Fsgg_hflx(num))
    end do
    
    ! Log output for model date

    if (my_task == master_task) then
       call seq_timemgr_EClockGetData(EClock,curr_ymd=cesmYMD, curr_tod=cesmTOD)
       write(stdout,F01) ' CESM Date ', cesmYMD,cesmTOD
       glcYMD = iyear*10000 + imonth*100 + iday
       glcTOD = ihour*3600 + iminute*60 + isecond
       write(stdout,F01) ' GLC  Date ',glcYMD,glcTOD
       call shr_sys_flush(stdout)
    end if

    ! If time to write history, do so

    hist_alarm = seq_timemgr_HistoryAlarmIsOn( EClock )
    if (hist_alarm) then
       ! TODO loop over instances
       call glc_io_write_history(ice_sheet%instances(1), EClock)
    endif

    ! If time to write restart, do so

    rest_alarm = seq_timemgr_RestartAlarmIsOn( EClock )
    if (rest_alarm) then
       ! TODO loop over instances
       call glc_io_write_restart(ice_sheet%instances(1), EClock)
    endif

    ! Reset shr logging to original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(stdout)
    
end subroutine glc_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: glc_final_mct
!
! !DESCRIPTION:
!     finalize method for glc model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine glc_final_mct( EClock, cdata, x2d, d2x)

! !USES:

    use glc_ensemble, only : get_inst_name

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        
   type(mct_aVect)             ,intent(inout) :: d2x        

!EOP

   integer(IN)             :: shrlogunit, shrloglev  
   character(CS)           :: myModelName


   !--- formats ---
   character(*), parameter :: F00   = "('(glc_final_mct) ',8a)"
   character(*), parameter :: F01   = "('(glc_final_mct) ',a,8i8)"
   character(*), parameter :: F91   = "('(glc_final_mct) ',73('-'))"
   character(*), parameter :: subName = "(glc_final_mct) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   ! Reset shr logging to my log file
   call shr_file_getLogUnit (shrlogunit)
   call shr_file_getLogLevel(shrloglev)
   call shr_file_setLogUnit (stdout)

   call get_inst_name(myModelName)

   if (my_task == master_task) then
      write(stdout,F91) 
      write(stdout,F00) trim(myModelName),': end of main integration loop'
      write(stdout,F91) 
   end if
      
   errorCode = glc_Success

   call glc_final(errorCode)

   if (verbose .and. my_task == master_task) then
      write(stdout,F01) ' Done',errorCode
      call shr_sys_flush(stdout)
   endif

   ! Reset shr logging to original values

   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)
   call shr_sys_flush(stdout)

 end subroutine glc_final_mct

!=================================================================================
  subroutine glc_import_mct(x2g,ndx,&
       index_tsrf,index_topo,index_qice)

    !-------------------------------------------------------------------
    use glc_global_fields, only: tsfc, topo, qsmb       ! from coupler

    type(mct_aVect),intent(inout) :: x2g
    integer(IN), intent(in) :: ndx                      ! elevation class
    integer(IN), intent(in) :: index_tsrf
    integer(IN), intent(in) :: index_topo
    integer(IN), intent(in) :: index_qice

    integer(IN) :: j,jj,i,g,nxg,nyg,n
    character(*), parameter :: subName = "(glc_import_mct) "
    !-------------------------------------------------------------------

    nxg = glc_grid%nx
    nyg = glc_grid%ny
    do j = 1, nyg           ! S to N
       jj = nyg - j + 1     ! reverse j index for glint grid (N to S)
       do i = 1, nxg
          g = (j-1)*nxg + i   ! global index (W to E, S to N)
          tsfc(i,jj,ndx) = x2g%rAttr(index_tsrf,g) - tkfrz
          topo(i,jj,ndx) = x2g%rAttr(index_topo,g)
          qsmb(i,jj,ndx) = x2g%rAttr(index_qice,g)
       enddo
    enddo

    if (verbose .and. my_task==master_task) then
       write(stdout,*) ' '
       write(stdout,*) subname,' x2g tsrf ',ndx,minval(x2g%rAttr(index_tsrf,:)),maxval(x2g%rAttr(index_tsrf,:))
       write(stdout,*) subname,' x2g topo ',ndx,minval(x2g%rAttr(index_topo,:)),maxval(x2g%rAttr(index_topo,:))
       write(stdout,*) subname,' x2g qice ',ndx,minval(x2g%rAttr(index_qice,:)),maxval(x2g%rAttr(index_qice,:))
       call shr_sys_flush(stdout)
    endif

  end subroutine glc_import_mct

!=================================================================================

  subroutine glc_export_mct(g2x,ndx,&
       index_frac,index_topo,&
       index_rofi,index_rofl,index_hflx)

    !-------------------------------------------------------------------
    use glc_global_fields, only: gfrac, gtopo, grofi, grofl, ghflx   ! to coupler

    type(mct_aVect),intent(inout) :: g2x
    integer(IN), intent(in) :: ndx
    integer(IN), intent(in) :: index_frac
    integer(IN), intent(in) :: index_topo
    integer(IN), intent(in) :: index_rofi
    integer(IN), intent(in) :: index_rofl
    integer(IN), intent(in) :: index_hflx

    integer(IN) :: j,jj,i,g,nxg,nyg,n
    character(*), parameter :: subName = "(glc_export_mct) "
    !-------------------------------------------------------------------

    nxg = glc_grid%nx
    nyg = glc_grid%ny
    do j = 1, nyg           ! S to N
       jj = nyg - j + 1     ! reverse j index for glint grid (N to S)
       do i = 1, nxg
          g = (j-1)*nxg + i ! global index (W to E, S to N)
          g2x%rAttr(index_frac,g) = gfrac(i,jj,ndx)
          g2x%rAttr(index_topo,g) = gtopo(i,jj,ndx)
          g2x%rAttr(index_rofi,g) = grofi(i,jj,ndx)
          g2x%rAttr(index_rofl,g) = grofl(i,jj,ndx)
          g2x%rAttr(index_hflx,g) = ghflx(i,jj,ndx)
       enddo
    enddo

    if (verbose .and. my_task==master_task) then
       write(stdout,*) subname,' g2x frac ',ndx,minval(g2x%rAttr(index_frac,:)),maxval(g2x%rAttr(index_frac,:))
       write(stdout,*) subname,' g2x topo ',ndx,minval(g2x%rAttr(index_topo,:)),maxval(g2x%rAttr(index_topo,:))
       write(stdout,*) subname,' g2x rofi ',ndx,minval(g2x%rAttr(index_rofi,:)),maxval(g2x%rAttr(index_rofi,:))
       write(stdout,*) subname,' g2x rofl ',ndx,minval(g2x%rAttr(index_rofl,:)),maxval(g2x%rAttr(index_rofl,:))
       write(stdout,*) subname,' g2x hflx ',ndx,minval(g2x%rAttr(index_hflx,:)),maxval(g2x%rAttr(index_hflx,:))
       call shr_sys_flush(stdout)
    endif

  end subroutine glc_export_mct

!=================================================================================

  subroutine glc_SetgsMap_mct( mpicom_g, GLCID, gsMap_g )

    integer        , intent(in)  :: mpicom_g
    integer        , intent(in)  :: GLCID
    type(mct_gsMap), intent(out) :: gsMap_g

    ! Local Variables

    integer,allocatable :: gindex(:)
    integer :: i, j, n, nxg, nyg
    integer :: lsize
    integer :: ier

    !--- formats ---
    character(*), parameter :: F02   = "('(glc_SetgsMap_mct) ',a,4es13.6)"
    character(*), parameter :: subName = "(glc_SetgsMap_mct) "
    !-------------------------------------------------------------------

    ! Note that the following assumes that the master task is responsible for all points

    if (my_task == master_task) then

       nxg = glc_grid%nx
       nyg = glc_grid%ny
       lsize = nxg*nyg

       ! Initialize MCT global seg map (the simple method used here only works because the
       ! master task is responsible for all points)

       allocate(gindex(lsize))
       do j = 1,nyg
          do i = 1,nxg
             n = (j-1)*nxg + i
             gindex(n) = n
          enddo
       enddo

    else
       
       nxg = 0
       nyg = 0
       lsize = 0
       allocate(gindex(lsize))

    end if

    call mct_gsMap_init( gsMap_g, gindex, mpicom_g, GLCID, lsize )

    deallocate(gindex)

  end subroutine glc_SetgsMap_mct

!===============================================================================

  subroutine glc_domain_mct( lsize, gsMap_g, dom_g )

    !-------------------------------------------------------------------
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(inout) :: gsMap_g
    type(mct_ggrid), intent(out)   :: dom_g      

    ! Local Variables

    integer :: g,i,j,n,nxg,nyg    ! index
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    character(*), parameter :: subName = "(glc_domain_mct) "
    !-------------------------------------------------------------------

    nxg = glc_grid%nx
    nyg = glc_grid%ny

    ! Initialize mct domain type

    call mct_gGrid_init( GGrid=dom_g, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )

    ! Initialize attribute vector with special value

    allocate(data(lsize))
    dom_g%data%rAttr(:,:) = -9999.0_R8
    dom_g%data%iAttr(:,:) = -9999
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_g,"mask" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_g,"frac" ,data,lsize) 

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT

    call mct_gsMap_orderedPoints(gsMap_g, my_task, idata)
    call mct_gGrid_importIAttr(dom_g,'GlobGridNum',idata,lsize)

    ! Fill in correct values for domain components
    ! lat/lon in degrees, area in radians^2, real-valued mask and frac
    do j = 1,nyg
    do i = 1,nxg
       n = (j-1)*nxg + i
       data(n) = glc_grid%lons(i)   ! Note: degrees, not radians
    end do
    end do
    call mct_gGrid_importRattr(dom_g,"lon",data,lsize) 

    do j = 1,nyg
    do i = 1,nxg
       n = (j-1)*nxg + i
       data(n) = glc_grid%lats(j)   ! Note: degrees, not radians
    end do
    end do
    call mct_gGrid_importRattr(dom_g,"lat",data,lsize) 

    do j = 1,nyg
    do i = 1,nxg
       n = (j-1)*nxg + i
       data(n) = glc_grid%box_areas(i,j)/(radius*radius)
    end do
    end do
    call mct_gGrid_importRattr(dom_g,"area",data,lsize) 

    ! Note: glc_landmask and glc_landfrac are read from file when grid is initialized
    do j = 1,nyg
    do i = 1,nxg
       n = (j-1)*nxg + i
       data(n) = real(glc_landmask(i,j), r8)  ! data is r8, glc_landmask is i4
    end do
    end do

    call mct_gGrid_importRattr(dom_g,"mask",data,lsize) 

    do j = 1,nyg
    do i = 1,nxg
       n = (j-1)*nxg + i
       data(n) = glc_landfrac(i,j)
    end do
    end do
    call mct_gGrid_importRattr(dom_g,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

    if (verbose .and. my_task==master_task) then
       i = mct_aVect_nIattr(dom_g%data)
       do n = 1,i
          write(stdout,*) subname,' dom_g ',n,minval(dom_g%data%iAttr(n,:)),maxval(dom_g%data%iAttr(n,:))
       enddo
       i = mct_aVect_nRattr(dom_g%data)
       do n = 1,i
          write(stdout,*) subname,' dom_g ',n,minval(dom_g%data%rAttr(n,:)),maxval(dom_g%data%rAttr(n,:))
       enddo
       call shr_sys_flush(stdout)
    endif

  end subroutine glc_domain_mct
    
!===============================================================================

end module glc_comp_mct
