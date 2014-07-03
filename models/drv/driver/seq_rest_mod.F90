!===============================================================================
! SVN $Id: seq_rest_mod.F90 45286 2013-03-26 18:17:04Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/trunk_tags/drvseq4_2_33/driver/seq_rest_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_rest_mod -- cpl7 restart reading/writing routines
!
! !DESCRIPTION:
!
!    Reads & writes cpl7 restart files
!
! !REMARKS:
!
!    aVect, domain, and fraction info accessed via seq_avdata_mod
!    to avoid excessively long routine arg lists.
!
! !REVISION HISTORY:
!     2009-Sep-25 - B. Kauffman - move from cpl7 main program into rest module
!     2007-mmm-dd - T. Craig - initial restart functionality
!
! !INTERFACE: ------------------------------------------------------------------

module seq_rest_mod

! !USES:

   use shr_kind_mod,      only: R8 => SHR_KIND_R8, IN => SHR_KIND_IN
   use shr_kind_mod,      only: CL => SHR_KIND_CL, CS => SHR_KIND_CS
   use shr_sys_mod,       only: shr_sys_abort, shr_sys_flush
   use shr_mpi_mod,       only: shr_mpi_bcast
   use shr_cal_mod,       only: shr_cal_date2ymd
   use shr_file_mod,      only: shr_file_getunit, shr_file_freeunit
   use mct_mod
   use ESMF

   use seq_avdata_mod    ! drv aVects & associated domain, fract, cdata (public read/write data)
   use seq_diag_mct      ! diagnostic routines                          (public read/write data)
   use seq_comm_mct      ! Sets mpi communicators, logunit and loglevel
   use seq_cdata_mod     ! "cdata" type & methods (domain + decomp + infodata in one datatype)
   use seq_infodata_mod  ! "infodata" gathers various control flags into one datatype
   use seq_timemgr_mod   ! clock & alarm routines
   use seq_io_mod        ! lower level io routines

   implicit none

   private

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS

   public :: seq_rest_read   ! read  cpl7 restart data
   public :: seq_rest_write  ! write cpl7 restart data

! !PUBLIC DATA MEMBERS:

   ! no public data

!EOP

   !----------------------------------------------------------------------------
   ! local data
   !----------------------------------------------------------------------------

   logical     :: iamin_CPLID            ! pe associated with CPLID
   integer(IN) :: mpicom_GLOID           ! MPI global communicator
   integer(IN) :: mpicom_CPLID           ! MPI cpl communicator

   integer(IN) :: nthreads_GLOID         ! OMP global number of threads
   integer(IN) :: nthreads_CPLID         ! OMP cpl number of threads
   logical     :: drv_threading          ! driver threading control

   logical     :: atm_present            ! .true.  => atm is present
   logical     :: lnd_present            ! .true.  => land is present
   logical     :: ice_present            ! .true.  => ice is present
   logical     :: ocn_present            ! .true.  => ocn is present
   logical     :: rof_present            ! .true.  => land runoff is present
   logical     :: rof_prognostic         ! .true.  => rof comp expects input
   logical     :: glc_present            ! .true.  => glc is present
   logical     :: sno_present            ! .true.  => land sno is present
   logical     :: wav_present            ! .true.  => wav is present

   logical     :: atm_prognostic         ! .true.  => atm comp expects input
   logical     :: lnd_prognostic         ! .true.  => lnd comp expects input
   logical     :: ice_prognostic         ! .true.  => ice comp expects input
   logical     :: ocn_prognostic         ! .true.  => ocn comp expects input
   logical     :: ocnrof_prognostic      ! .true.  => ocn comp expects runoff input
   logical     :: glc_prognostic         ! .true.  => glc comp expects input
   logical     :: sno_prognostic         ! .true.  => sno comp expects input
   logical     :: wav_prognostic         ! .true.  => wav comp expects input

   integer(IN) :: info_debug = 0         ! local info_debug level

!===============================================================================
contains
!===============================================================================

subroutine seq_rest_read(rest_file)

   implicit none

   character(*),intent(in) :: rest_file  ! restart file path/name

   integer(IN)          :: n,n1,n2,n3
   real(r8),allocatable :: ds(:)         ! for reshaping diag data for restart file
   real(r8),allocatable :: ns(:)         ! for reshaping diag data for restart file
   character(CS)        :: string
   integer(IN)          :: ierr          ! MPI error return
   type(mct_gsMap),pointer :: gsmap
   character(len=*),parameter :: subname = "(seq_rest_read) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! get required infodata
   !----------------------------------------------------------------------------
   iamin_CPLID  = seq_comm_iamin(CPLID)
   call seq_comm_setptrs(GLOID,mpicom=mpicom_GLOID,nthreads=nthreads_GLOID)
   call seq_comm_setptrs(CPLID,mpicom=mpicom_CPLID,nthreads=nthreads_CPLID)
   call seq_infodata_getData(infodata,drv_threading=drv_threading)
   call seq_infodata_getData(infodata, &
        atm_present=atm_present, &
        lnd_present=lnd_present, &
        rof_present=rof_present, &
        ice_present=ice_present, &
        ocn_present=ocn_present, &
        glc_present=glc_present, &
        wav_present=wav_present, &
        sno_present=sno_present  )
   call seq_infodata_getData(infodata, &
        atm_prognostic=atm_prognostic, &
        lnd_prognostic=lnd_prognostic, &
        ice_prognostic=ice_prognostic, &
        ocn_prognostic=ocn_prognostic, &
        rof_prognostic=rof_prognostic, &
        ocnrof_prognostic=ocnrof_prognostic, &
        glc_prognostic=glc_prognostic, &
        wav_prognostic=wav_prognostic, &
        sno_prognostic=sno_prognostic  )

   if (iamin_CPLID) then
      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
      if (atm_present) then
         call seq_cdata_setptrs(cdata_ax,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,fractions_ax,'fractions_ax')
         call seq_io_read(rest_file,gsmap,a2x_ax,'a2x_ax')
      endif
      if (lnd_present) then
         call seq_cdata_setptrs(cdata_lx,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,fractions_lx,'fractions_lx')
      endif
      if (ocn_present) then
         call seq_cdata_setptrs(cdata_ox,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,fractions_ox,'fractions_ox')
         call seq_io_read(rest_file,gsmap,o2x_ox,'o2x_ox')
         call seq_io_read(rest_file,gsmap,x2oacc_ox,'x2oacc_ox')
         call seq_io_read(rest_file,gsmap,xao_ox,'xao_ox')
         call seq_io_read(rest_file,      x2oacc_ox_cnt,'x2oacc_ox_cnt')
         call seq_cdata_setptrs(cdata_ax,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,xao_ax,'xao_ax')
      endif
      if (ice_present) then
         call seq_cdata_setptrs(cdata_ix,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,fractions_ix,'fractions_ix')
         call seq_io_read(rest_file,gsmap,i2x_ix,'i2x_ix')
      endif
      if (rof_present) then
         call seq_cdata_setptrs(cdata_rx,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,fractions_rx,'fractions_rx')
         call seq_io_read(rest_file,gsmap,r2x_rx,'r2x_rx')
      endif
      if (lnd_present .and. rof_present .and. rof_prognostic) then
         call seq_cdata_setptrs(cdata_lx,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,x2racc_lx,'x2racc_lx')
         call seq_io_read(rest_file,      x2racc_lx_cnt ,'x2racc_lx_cnt')
      end if
      if (rof_present .and. ocnrof_prognostic) then
         call seq_cdata_setptrs(cdata_rx,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,r2xacc_rx,'r2xacc_rx')
         call seq_io_read(rest_file,      r2xacc_rx_cnt,'r2xacc_rx_cnt')
      endif
      if (glc_present) then
         call seq_cdata_setptrs(cdata_gx,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,fractions_gx,'fractions_gx')
      endif
      if (sno_present) then
         call seq_cdata_setptrs(cdata_sx,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,x2s_sx,'x2s_sx')
      endif
      if (wav_present) then
         call seq_cdata_setptrs(cdata_wx,gsmap=gsmap)
         call seq_io_read(rest_file,gsmap,fractions_wx,'fractions_wx')
         call seq_io_read(rest_file,gsmap,w2x_wx,'w2x_wx')
      endif

      n = size(budg_dataG)
      allocate(ds(n),ns(n))
      call seq_io_read(rest_file,ds,'budg_dataG')
      call seq_io_read(rest_file,ns,'budg_ns')

      n = 0
      do n1 = 1,size(budg_dataG,dim=1)
      do n2 = 1,size(budg_dataG,dim=2)
      do n3 = 1,size(budg_dataG,dim=3)
         n = n + 1
         budg_dataG(n1,n2,n3) = ds(n)
         budg_ns   (n1,n2,n3) = ns(n)
      enddo
      enddo
      enddo
!     call shr_mpi_bcast(budg_dataG,cpl_io_root) ! not necessary, io lib does bcast

      deallocate(ds,ns)

      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

   endif

end subroutine seq_rest_read

!===============================================================================

subroutine seq_rest_write(EClock_d,seq_SyncClock)

   implicit none

   type(ESMF_Clock)      ,intent(in)    :: EClock_d      ! driver clock
   type(seq_timemgr_type),intent(inout) :: seq_SyncClock ! contains ptr to driver clock

   integer(IN)   :: n,n1,n2,n3,fk
   integer(IN)   :: curr_ymd         ! Current date YYYYMMDD
   integer(IN)   :: curr_tod         ! Current time-of-day (s)
   integer(IN)   :: yy,mm,dd         ! year, month, day
   character(CL) :: case_name        ! case name
   character(CL) :: cvar             ! char variable
   integer(IN)   :: ivar             ! integer variable
   real(r8)      :: rvar             ! real variable
   logical       :: whead,wdata      ! flags header/data writing
   logical       :: cdf64            ! true => create netCDF with 64 bit addressing
   logical       :: cplroot          ! root pe on cpl id
   integer(IN)   :: iun              ! unit number
   character(CL) :: rest_file        ! Local path to restart filename
   integer(IN)   :: ierr             ! MPI error return
   type(mct_gsMap),pointer :: gsmap

   real(r8),allocatable :: ds(:)     ! for reshaping diag data for restart file
   real(r8),allocatable :: ns(:)     ! for reshaping diag data for restart file
   character(len=*),parameter :: subname = "(seq_rest_write) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! get required infodata
   !----------------------------------------------------------------------------
   iamin_CPLID  = seq_comm_iamin(CPLID)
   call seq_comm_setptrs(GLOID,mpicom=mpicom_GLOID,nthreads=nthreads_GLOID)
   call seq_comm_setptrs(CPLID,mpicom=mpicom_CPLID,nthreads=nthreads_CPLID)
   call seq_comm_setptrs(CPLID,iamroot=cplroot)
   call seq_infodata_getData(infodata,drv_threading=drv_threading)
   call seq_infodata_getData(infodata, &
        atm_present=atm_present, &
        lnd_present=lnd_present, &
        rof_present=rof_present, &
        ice_present=ice_present, &
        ocn_present=ocn_present, &
        glc_present=glc_present, &
        wav_present=wav_present, &
        sno_present=sno_present  )
   call seq_infodata_getData(infodata, &
        atm_prognostic=atm_prognostic, &
        lnd_prognostic=lnd_prognostic, &
        ice_prognostic=ice_prognostic, &
        rof_prognostic=rof_prognostic, &
        ocn_prognostic=ocn_prognostic, &
        ocnrof_prognostic=ocnrof_prognostic, &
        glc_prognostic=glc_prognostic, &
        wav_prognostic=wav_prognostic, &
        sno_prognostic=sno_prognostic  )
   call seq_infodata_getData(infodata, cpl_cdf64=cdf64 )

   ! Write out infodata and time manager data to restart file

   call seq_infodata_GetData( infodata, case_name=case_name)
   call seq_timemgr_EClockGetData( EClock_d, curr_ymd=curr_ymd, curr_tod=curr_tod)
   call shr_cal_date2ymd(curr_ymd,yy,mm,dd)
   write(rest_file,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
      trim(case_name), '.cpl.r.', yy,'-',mm,'-',dd,'-',curr_tod,'.nc'

   ! Write driver data to restart file

   if (iamin_CPLID) then

      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

      ! copy budg_dataG into 1d array
      n = size(budg_dataG)
      allocate(ds(n),ns(n))
      call shr_mpi_bcast(budg_dataG,mpicom_CPLID) ! pio requires data on all pe's?

      n = 0
      do n1 = 1,size(budg_dataG,dim=1)
      do n2 = 1,size(budg_dataG,dim=2)
      do n3 = 1,size(budg_dataG,dim=3)
         n = n + 1
         ds(n) = budg_dataG(n1,n2,n3)
         ns(n) = budg_ns(n1,n2,n3)
      enddo
      enddo
      enddo

      if (cplroot) then
         iun = shr_file_getUnit()
         call seq_infodata_GetData(infodata,restart_pfile=cvar)
         if (loglevel > 0) write(logunit,"(3A)") subname," write rpointer file ", &
            trim(cvar)
         open(iun, file=cvar, form='FORMATTED')
         write(iun,'(a)') rest_file
         close(iun)
         call shr_file_freeUnit( iun )
      endif

      call shr_mpi_bcast(rest_file,mpicom_CPLID)
      call seq_io_wopen(rest_file,clobber=.true.,cdf64=cdf64)

      ! loop twice (for perf), first time write header, second time write data
      do fk = 1,2
         if (fk == 1) then
            whead = .true.
            wdata = .false.
!!            call seq_io_redef(rest_file)
         elseif (fk == 2) then
            whead = .false.
            wdata = .true.
            call seq_io_enddef(rest_file)
         else
            call shr_sys_abort('driver_write_rstart fk illegal')
         end if
         call seq_infodata_GetData(infodata,nextsw_cday=rvar)
         call seq_io_write(rest_file,rvar,'seq_infodata_nextsw_cday',whead=whead,wdata=wdata)
         call seq_infodata_GetData(infodata,precip_fact=rvar)
         call seq_io_write(rest_file,rvar,'seq_infodata_precip_fact',whead=whead,wdata=wdata)
         call seq_infodata_GetData(infodata,case_name=cvar)
         call seq_io_write(rest_file,trim(cvar),'seq_infodata_case_name',whead=whead,wdata=wdata)

         call seq_timemgr_EClockGetData( EClock_d, start_ymd=ivar)
         call seq_io_write(rest_file,ivar,'seq_timemgr_start_ymd',whead=whead,wdata=wdata)
         call seq_timemgr_EClockGetData( EClock_d, start_tod=ivar)
         call seq_io_write(rest_file,ivar,'seq_timemgr_start_tod',whead=whead,wdata=wdata)
         call seq_timemgr_EClockGetData( EClock_d, ref_ymd=ivar)
         call seq_io_write(rest_file,ivar,'seq_timemgr_ref_ymd'  ,whead=whead,wdata=wdata)
         call seq_timemgr_EClockGetData( EClock_d, ref_tod=ivar)
         call seq_io_write(rest_file,ivar,'seq_timemgr_ref_tod'  ,whead=whead,wdata=wdata)
         call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ivar)
         call seq_io_write(rest_file,ivar,'seq_timemgr_curr_ymd' ,whead=whead,wdata=wdata)
         call seq_timemgr_EClockGetData( EClock_d, curr_tod=ivar)
         call seq_io_write(rest_file,ivar,'seq_timemgr_curr_tod' ,whead=whead,wdata=wdata)

         call seq_io_write(rest_file,ds,'budg_dataG',whead=whead,wdata=wdata)
         call seq_io_write(rest_file,ns,'budg_ns',whead=whead,wdata=wdata)

         if (atm_present) then
            call seq_cdata_setptrs(cdata_ax,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,fractions_ax,'fractions_ax',whead=whead,wdata=wdata)
            call seq_io_write(rest_file,gsmap,a2x_ax,'a2x_ax',whead=whead,wdata=wdata)
         endif
         if (lnd_present) then
            call seq_cdata_setptrs(cdata_lx,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,fractions_lx,'fractions_lx',whead=whead,wdata=wdata)
         endif
         if (ocn_present) then
            call seq_cdata_setptrs(cdata_ox,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,fractions_ox,'fractions_ox',whead=whead,wdata=wdata)
            call seq_io_write(rest_file,gsmap,o2x_ox,'o2x_ox',whead=whead,wdata=wdata)
            call seq_io_write(rest_file,gsmap,x2oacc_ox,'x2oacc_ox',whead=whead,wdata=wdata)
            call seq_io_write(rest_file,      x2oacc_ox_cnt,'x2oacc_ox_cnt',whead=whead,wdata=wdata)
            call seq_io_write(rest_file,gsmap,xao_ox,'xao_ox',whead=whead,wdata=wdata)
            call seq_cdata_setptrs(cdata_ax,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,xao_ax,'xao_ax',whead=whead,wdata=wdata)
         endif
         if (ice_present) then
            call seq_cdata_setptrs(cdata_ix,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,fractions_ix,'fractions_ix',whead=whead,wdata=wdata)
            call seq_io_write(rest_file,gsmap,i2x_ix,'i2x_ix',whead=whead,wdata=wdata)
         endif
         if (rof_present) then
            call seq_cdata_setptrs(cdata_rx,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,fractions_rx,'fractions_rx',whead=whead,wdata=wdata)
            call seq_io_write(rest_file,gsmap,r2x_rx,'r2x_rx',whead=whead,wdata=wdata)
         endif
         if (lnd_present .and. rof_present .and. rof_prognostic) then
            call seq_cdata_setptrs(cdata_lx,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,x2racc_lx,'x2racc_lx'    ,whead=whead,wdata=wdata)
            call seq_io_write(rest_file,      x2racc_lx_cnt ,'x2racc_lx_cnt',whead=whead,wdata=wdata)
         end if
         if (rof_present .and. ocnrof_prognostic) then
            call seq_cdata_setptrs(cdata_rx,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,r2xacc_rx,'r2xacc_rx',whead=whead,wdata=wdata)
            call seq_io_write(rest_file,      r2xacc_rx_cnt,'r2xacc_rx_cnt',whead=whead,wdata=wdata)
         endif
         if (glc_present) then
            call seq_cdata_setptrs(cdata_gx,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,fractions_gx,'fractions_gx',whead=whead,wdata=wdata)
         endif
         if (sno_present) then
            call seq_cdata_setptrs(cdata_sx,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,x2s_sx,'x2s_sx',whead=whead,wdata=wdata)
         endif
         if (wav_present) then
            call seq_cdata_setptrs(cdata_wx,gsmap=gsmap)
            call seq_io_write(rest_file,gsmap,fractions_wx,'fractions_wx',whead=whead,wdata=wdata)
            call seq_io_write(rest_file,gsmap,w2x_wx,'w2x_wx',whead=whead,wdata=wdata)
         endif
      enddo

      call seq_io_close(rest_file)
      deallocate(ds,ns)

      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   endif
end subroutine seq_rest_write

!===============================================================================

end module seq_rest_mod
