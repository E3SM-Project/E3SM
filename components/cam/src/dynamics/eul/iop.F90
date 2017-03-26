
module iop
!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: iop
! 
! !DESCRIPTION: 
! iop specific routines
!
! !USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_scam_mod, only: shr_scam_GetCloseLatLon
  use constituents, only: readtrace, cnst_get_ind, pcnst, cnst_name
  use string_utils, only: to_lower
  use pmgrid
  use prognostics
  use time_manager, only: timemgr_init, get_curr_date, get_curr_calday,&
                          get_nstep,get_start_date,timemgr_time_inc
  use cam_abortutils,   only: endrun
  use scamMod
  use wrap_nf
  use cam_logfile,  only: iulog
  use phys_control, only: phys_getopts
  use eul_control_mod,only: eul_nsplit
!
! !PUBLIC TYPES:
  implicit none

  private

  real(r8), allocatable, target :: dqfx3sav(:,:,:,:)       
  real(r8), allocatable, target :: t2sav(:,:,:)       
  real(r8), allocatable, target :: divq3dsav(:,:,:,:)
  real(r8), allocatable, target :: divt3dsav(:,:,:)       
  real(r8), allocatable, target :: betasav(:)
  integer :: closelatidx,closelonidx,latid,lonid,levid,timeid

  real(r8):: closelat,closelon

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_iop_fields
  public :: readiopdata         ! read iop boundary data
  public :: setiopupdate        ! find index in iopboundary data for current time
!  public :: scam_use_iop_srf
! !PUBLIC DATA:
  public betasav, &
         dqfx3sav, divq3dsav, divt3dsav,t2sav

!
! !REVISION HISTORY:
! Created by John Truesdale
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
!----------------------------------------------------------------------- 

contains
   subroutine init_iop_fields()
!------------------------------------------------------------------------------
! Coupler for converting dynamics output variables into physics input variables
! also writes dynamics variables (on physics grid) to history file
!------------------------------------------------------------------------------
   implicit none

!-----------------------------------------------------------------------
   if (eul_nsplit>1) then
      call endrun('iop module cannot be used with eul_nsplit>1')
   endif
	        
   if(.not.allocated(betasav)) then
      allocate (betasav(beglat:endlat))
      betasav(:)=0._r8
   endif

   if(.not.allocated(dqfx3sav)) then
      allocate (dqfx3sav(plon,plev,pcnst,beglat:endlat))
      dqfx3sav(:,:,:,:)=0._r8
   endif
   if(.not.allocated(divq3dsav)) then
      allocate (divq3dsav(plon,plev,pcnst,beglat:endlat))
      divq3dsav(:,:,:,:)=0._r8
   endif
   if(.not.allocated(divt3dsav)) then
      allocate (divt3dsav(plon,plev,beglat:endlat))
      divt3dsav(:,:,:)=0._r8
   endif
   if(.not.allocated(t2sav)) then
      allocate (t2sav(plon,plev,beglat:endlat))  ! temp tendency
      t2sav(:,:,:)=0._r8
   endif
  end subroutine init_iop_fields

subroutine readiopdata( )


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
        use ppgrid,           only: begchunk, endchunk
	use phys_grid,        only: clat_p
	use commap
        use string_utils,     only: to_lower
        use getinterpnetcdfdata
        use shr_sys_mod,      only: shr_sys_flush
	use hycoef, only : hyam, hybm
        use eul_control_mod
        use error_messages, only : handle_ncerr
        use netcdf
!-----------------------------------------------------------------------
   implicit none
#if ( defined RS6000 )
   implicit automatic ( a-z )
#endif
!------------------------------Locals-----------------------------------
!     
   integer NCID, status
   integer time_dimID, lev_dimID,lev_varID,mod_dimID,&
           mod_varID,sps_varID,sps_dimID
   integer tsec_varID, bdate_varID,varid
   integer i,j
   integer nlev, nmod, nsps
   integer total_levs

   integer bdate, ntime
   integer, allocatable :: tsec(:)
   integer k, m
   integer icldliq,icldice
   integer inumliq,inumice

   logical have_srf              ! value at surface is available
   logical fill_ends             ! 
   logical have_cnst(pcnst)
   real(r8), allocatable :: dplevs( : ), aitken( :)
   integer, allocatable :: dmods( : ), dsps( : )
   real(r8) dummy
   real(r8) lat,xlat
   real(r8) srf(1)                  ! value at surface
   real(r8) pmid(plev)  ! pressure at model levels (time n)
   real(r8) pint(plevp) ! pressure at model interfaces (n  )
   real(r8) pdel(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) weight
   real(r8) tmpdata(1)
   real(r8) coldata(plev)
   integer strt4(4),cnt4(4)
   character(len=16) :: lowername

   fill_ends= .false.

!     
!     Open IOP dataset
!     
  call handle_ncerr( nf90_open (iopfile, 0, ncid),&
       'readiopdata.F90', __LINE__)

!
!     if the dataset is a CAM generated dataset set use_camiop to true
!       CAM IOP datasets have a global attribute called CAM_GENERATED_IOP      
!
   if ( nf90_inquire_attribute( ncid, NF90_GLOBAL, 'CAM_GENERATED_FORCING', attnum=i ).EQ. NF90_NOERR ) then
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
 
 if(scm_observed_aero) then
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
   if (.not.allocated(aitken)) then
      allocate(aitken(nmod))
      aitken= 1.0e30_R8
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
   if (.not.allocated(dsps)) then
      allocate(dsps(nsps))
      dsps=-999
   end if
   
  status = NF90_INQ_VARID( ncid, 'mod', mod_varID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for mode'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, mod_varID, dmods(:nmod)),&
                    'readiopdata.F90', __LINE__)
  
status = NF90_INQ_VARID( ncid, 'sps', sps_varID )
   if ( status .ne. nf90_noerr ) then
      write(iulog,* )'ERROR - readiopdata.F:Could not find variable ID for mode'
      status = NF90_CLOSE ( ncid )
      return
   end if

   call handle_ncerr( nf90_get_var (ncid, sps_varID, dsps(:nsps)),&
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
    
   if (.not.allocated(scm_div)) then
      allocate(scm_div(nmod,nsps))
      scm_div= 1.0e30_R8
   end if
   ! allocate(scm_aitken_div(nsps))
   ! allocate(scm_coarse_div(nsps))

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
!     convert pressure to millibars ( lev is expressed in pascals in iop datasets )
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
      status = nf90_get_var(ncid, varid, ps(1,1,n3), strt4)
      have_ps = .true.
   endif


!  for reproducing CAM output don't do interpolation.
!  the most expedient way of doing this is to set      
!  the dataset pressure levels to the current
!  scam model levels
	
   if ( use_camiop ) then
      do i = 1, plev
         dplevs( i ) = 1000.0_r8 * hyam( i ) + ps(1,1,n3) * hybm( i ) / 100.0_r8
      end do
   endif

!     add the surface pressure to the pressure level data, so that
!     surface boundary condition will be set properly,
!     making sure that it is the highest pressure in the array.
!

   total_levs = nlev+1
   dplevs(nlev+1) = ps(1,1,n3)/100 ! ps is expressed in pascals
   do i= nlev, 1, -1
      if ( dplevs(i) .GT. ps(1,1,n3)/100.0_r8) then
         total_levs = i
         dplevs(i) = ps(1,1,n3)/100
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

   if ( get_nstep() .le. 1  ) then
     tobs(:)=t3(1,:,1,1)
   else
     tobs(:)= t3(1,:,1,n3)
   endif

   if ( use_camiop ) then
     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,'t', have_tsair, &
          tsair(1), fill_ends, &
          dplevs, nlev,ps(1,1,n3),tobs, status )
   else
     call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,'T', have_tsair, &
          tsair(1), fill_ends, &
          dplevs, nlev,ps(1,1,n3), tobs, status )
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
      if (.not.use_camiop .and. get_nstep() .eq. 0) then
         do i=1, PLEV
            t3(1,i,1,n3)=tobs(i)  !     set t to tobs at first time step
         end do
      endif
   endif

   status = nf90_inq_varid( ncid, 'Tg', varid   )
   if (status .ne. nf90_noerr) then
      write(iulog,*)'Could not find variable Tg'
      if ( have_tsair ) then
         write(iulog,*) 'Using Tsair'
         tground = tsair     ! use surface value from T field
      else
         write(iulog,*) 'Using T at lowest level'
         tground = t3(1,plev,1,n3)
      endif
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,tground)
      have_tg = .true.
   endif

   status = nf90_inq_varid( ncid, 'qsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif
!
   if ( get_nstep() .le. 1  ) then
     qobs(:)= q3(1,:,1,1,1)
   else
     qobs(:)= q3(1,:,1,1,n3)
   endif

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'q', have_srf, &
      srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), qobs, status )
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
      if (.not. use_camiop .and. get_nstep() .eq. 0) then
         do i=1, PLEV
            q3(1,i,1,1,n3)=qobs(i)
         end do
         have_q = .true.
      endif
   endif

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'cld', .false., &
      dummy, fill_ends, dplevs, nlev,ps(1,1,n3), cldobs, status )
   if ( status .ne. nf90_noerr ) then
      have_cld = .false.
   else
      have_cld = .true.
   endif
   
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx,  'clwp', .false., &
      dummy, fill_ends, dplevs, nlev,ps(1,1,n3), clwpobs, status )
   if ( status .ne. nf90_noerr ) then
      have_clwp = .false.
   else
      have_clwp = .true.
   endif

!
!	read divq (horizontal advection)
!      
   status = nf90_inq_varid( ncid, 'divqsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
        'divq', have_srf, srf(1), fill_ends, &
        dplevs, nlev,ps(1,1,n3), divq(:,1), status )
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
        have_srf, srf(1), fill_ends, &
        dplevs, nlev,ps(1,1,n3), vertdivq(:,1), status )
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
      have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), divq3d(:,m), status )
   if ( status .ne. nf90_noerr ) then
      have_cnst(m) = .false.
      divq3d(1:,m)=0._r8
   else
      have_cnst(m) = .true.
   endif

      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_dqfx', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), coldata, status )
       if ( STATUS .NE. NF90_NOERR ) then
         dqfxcam=0._r8
      else
         dqfxcam(1,:,m)=coldata(:)
      endif

      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, trim(cnst_name(m))//'_alph', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), tmpdata, status )
      if ( status .ne. nf90_noerr ) then
!         have_cnst(m) = .false.
         alphacam(m)=0._r8
      else
          alphacam(m)=tmpdata(1)
!         have_cnst(m) = .true.
      endif

   end do


   call cnst_get_ind('NUMLIQ', inumliq, abort=.false.)
   if ( inumliq > 0 ) then
      have_srf = .false.
      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'NUMLIQ', &
         have_srf, srf(1), fill_ends, &
         dplevs, nlev,ps(1,1,n3), numliqobs, status )
      if ( status .ne. nf90_noerr ) then
         have_numliq = .false.
      else
         have_numliq = .true.
         if ( get_nstep() .eq. 0) then
            do i=1, PLEV
               q3(1,i,inumliq,1,n3)=numliqobs(i)
            end do
         endif
     endif
   end if

   call cnst_get_ind('CLDLIQ', icldliq)

   have_srf = .false.
   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'CLDLIQ', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), cldliqobs, status )
   if ( status .ne. nf90_noerr ) then
      have_cldliq = .false.
   else
      have_cldliq = .true.
      if ( get_nstep() .eq. 0) then
         do i=1, PLEV
            q3(1,i,icldliq,1,n3)=cldliqobs(i)
         end do
      endif
   endif

   call cnst_get_ind('CLDICE', icldice)

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'CLDICE', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), cldiceobs, status )
   if ( status .ne. nf90_noerr ) then
      have_cldice = .false.
   else
      have_cldice = .true.
      if ( get_nstep() .eq. 0) then
         do i=1, PLEV
            q3(1,i,icldice,1,n3)=cldiceobs(i)
         end do
      endif
   endif

   call cnst_get_ind('NUMICE', inumice, abort=.false.)
   if ( inumice > 0 ) then
      have_srf = .false.

      call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'NUMICE', &
         have_srf, srf(1), fill_ends, &
         dplevs, nlev,ps(1,1,n3), numiceobs, status )
      if ( status .ne. nf90_noerr ) then
         have_numice = .false.
      else
         have_numice = .true.
         if ( get_nstep() .eq. 0) then
            do i=1, PLEV
               q3(1,i,inumice,1,n3)=numiceobs(i)
            end do
         endif
      endif
   end if

!
!	read divu (optional field)
!      
   status = nf90_inq_varid( ncid, 'divusrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divu', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), divu, status )
   if ( status .ne. nf90_noerr ) then
      have_divu = .false.
   else
      have_divu = .true.
   endif
!
!	read divv (optional field)
!      
   status = nf90_inq_varid( ncid, 'divvsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divv', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), divv, status )
   if ( status .ne. nf90_noerr ) then
      have_divv = .false.
   else
      have_divv = .true.
   endif
!
!	read divt (optional field)
!      
   status = nf90_inq_varid( ncid, 'divtsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      'divT', have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), divt, status )
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
      have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), vertdivt, status )
   if ( status .ne. nf90_noerr ) then
      have_vertdivt = .false.
   else
      have_vertdivt = .true.
   endif
!
!	read divt3d (combined vertical/horizontal advection)
!      (optional field)

   status = nf90_inq_varid( ncid, 'divT3dsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      status = nf90_get_var(ncid, varid, srf(1), strt4)
      have_srf = .true.
   endif

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'divT3d', &
      have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), divt3d, status )
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
      'omega', .true., ptend, fill_ends, &
      dplevs, nlev,ps(1,1,n3), wfld, status )
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
   call plevs0(1    ,plon   ,plev    ,ps(1,1,n3)   ,pint,pmid ,pdel)
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

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      'u', have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), uobs, status )
   if ( status .ne. nf90_noerr ) then
      have_u = .false.
   else
      have_u = .true.
      if (.not. use_camiop .and. get_nstep() .eq. 0 ) then
         do i=1, PLEV
            u3(1,i,1,n3) = uobs(i)  !     set u to uobs at first time step
         end do
      endif
   endif

   status = nf90_inq_varid( ncid, 'vsrf', varid   )
   if ( status .ne. nf90_noerr ) then
      have_srf = .false.
   else
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,srf)
      have_srf = .true.
   endif

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, &
      'v', have_srf, srf(1), fill_ends, &
      dplevs, nlev,ps(1,1,n3), vobs, status )
   if ( status .ne. nf90_noerr ) then
      have_v = .false.
   else
      have_v = .true.
      if (.not. use_camiop .and. get_nstep() .eq. 0 ) then
         do i=1, PLEV
            v3(1,i,1,n3) = vobs(i)  !     set u to uobs at first time step
         end do
      endif
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
      .false., dummy, fill_ends, & ! datasets don't contain Q1 at surface
      dplevs, nlev,ps(1,1,n3), q1obs, status )
   if ( status .ne. nf90_noerr ) then
      have_q1 = .false.
   else
      have_q1 = .true.
   endif

   call getinterpncdata( ncid, scmlat, scmlon, ioptimeidx, 'Q2', &
      .false., dummy, fill_ends, & ! datasets don't contain Q2 at surface
      dplevs, nlev,ps(1,1,n3), q1obs, status )
   if ( status .ne. nf90_noerr ) then
      have_q2 = .false.
   else
      have_q2 = .true.
   endif

!  Test for BOTH 'lhflx' and 'lh' without overwriting 'have_lhflx'.  
!  Analagous changes made for the surface heat flux

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
   if ( .not. have_divt .or. .not. have_divq ) then
      use_3dfrc = .true.
   endif
   call shr_sys_flush( iulog )

   status =  nf90_inq_varid( ncid, 'CLAT', varid   )
   if ( status .eq. nf90_noerr ) then
      call wrap_get_vara_realx (ncid,varid,strt4,cnt4,clat)
      clat_p(1)=clat(1)
      latdeg(1) = clat(1)*45._r8/atan(1._r8)
   endif

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

   call shr_sys_flush( iulog )

   status = nf90_close( ncid )
   call shr_sys_flush( iulog )

   deallocate(dplevs,tsec)

   return
end subroutine readiopdata

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
         write(iulog,*) ' IOP start is        ',next_date,' and ',next_sec,' seconds'
         write(iulog,*) ' IOP end is          ',last_date,' and ',last_sec,' seconds'
         call endrun
      endif

      doiopupdate = .true.

!------------------------------------------------------------------------------
!     Check if iop data needs to be updated and set doiopupdate accordingly
!------------------------------------------------------------------------------
   else                      ! endstep > 1

      call timemgr_time_inc(bdate, 0, next_date, next_sec, inc_s=tsec(iopTimeIdx+1))

      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day

      if ( ncdate .gt. next_date .or. (ncdate .eq. next_date &
         .and. ncsec .ge. next_sec)) then
         iopTimeIdx = iopTimeIdx + 1
         doiopupdate = .true.
#if DEBUG > 2
         write(iulog,*) 'nstep = ',get_nstep()
         write(iulog,*) 'ncdate=',ncdate,' ncsec=',ncsec
         write(iulog,*) 'next_date=',next_date,' next_sec=',next_sec
         write(iulog,*)'******* do iop update'
#endif 
      else
         doiopupdate = .false.
      end if
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

end module iop

