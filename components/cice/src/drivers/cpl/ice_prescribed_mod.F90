!===================================================================
!BOP 
!
! !MODULE: ice_prescribed_mod - Prescribed Ice Model
!
! !DESCRIPTION:
!
! The prescribed ice model reads in ice concentration data from a netCDF
! file.  Ice thickness, temperature, the ice temperature profile are
! prescribed.  Air/ice fluxes are computed to get surface temperature,
! Ice/ocean fluxes are set to zero, and ice dynamics are not calculated.
! Regridding and data cycling capabilities are included.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_prescribed_mod.F90 40 2006-12-01 19:09:30Z eclare $
!
! 2010-May-15 - Tony Craig and Mariana Vertenstein - updated to latest streams
! 2006-Aug-22 - D. Bailey, E. Hunke, modified to fit with CICE
! 2005-May-19 - J. Schramm - first version
! 2005-Apr-19 - B. Kauffman, J. Schramm, M. Vertenstein, NCAR - design
!
! !INTERFACE: ----------------------------------------------------------
 
module ice_prescribed_mod

! !USES:

   use shr_strdata_mod
   use shr_dmodel_mod
   use shr_string_mod
   use shr_ncread_mod
   use shr_sys_mod
   use shr_mct_mod
   use mct_mod
   use pio

   use ice_broadcast
   use ice_communicate, only : my_task, master_task, MPI_COMM_ICE
   use ice_kinds_mod
   use ice_fileunits
   use ice_exit,        only : abort_ice
   use ice_domain_size, only : nx_global, ny_global, ncat, nilyr, nslyr, max_blocks
   use ice_constants
   use ice_blocks,     only : nx_block, ny_block
   use ice_domain,     only : nblocks, distrb_info, blocks_ice
   use ice_grid,       only : TLAT,TLON,hm,tmask
   use ice_calendar,   only : idate, sec, calendar_type
   use ice_itd,        only : ilyr1, slyr1, hin_max
   use ice_read_write

   implicit none
   save

   private ! except


! !PUBLIC TYPES:

! !PUBLIC MEMBER FUNCTIONS:

   public :: ice_prescribed_init      ! initialize input data stream
   public :: ice_prescribed_run       ! get time slices and time interp
   public :: ice_prescribed_phys      ! set prescribed ice state and fluxes

! !PUBLIC DATA MEMBERS:

   logical(kind=log_kind), public :: prescribed_ice      ! true if prescribed ice

!EOP

   integer(SHR_KIND_IN),parameter :: nFilesMaximum = 400 ! max number of files
   integer(kind=int_kind)         :: stream_year_first   ! first year in stream to use
   integer(kind=int_kind)         :: stream_year_last    ! last year in stream to use
   integer(kind=int_kind)         :: model_year_align    ! align stream_year_first 
                                                         ! with this model year

   character(len=char_len_long)   :: stream_fldVarName
   character(len=char_len_long)   :: stream_fldFileName(nFilesMaximum)
   character(len=char_len_long)   :: stream_domTvarName
   character(len=char_len_long)   :: stream_domXvarName
   character(len=char_len_long)   :: stream_domYvarName
   character(len=char_len_long)   :: stream_domAreaName
   character(len=char_len_long)   :: stream_domMaskName
   character(len=char_len_long)   :: stream_domFileName
   character(len=char_len_long)   :: stream_mapread
   logical(kind=log_kind)         :: prescribed_ice_fill        ! true if data fill required

   type(shr_strdata_type)       :: sdat         ! prescribed data stream
   character(len=char_len_long) :: fldList      ! list of fields in data stream
   real(kind=dbl_kind)          :: ice_cov(nx_block,ny_block,max_blocks) ! ice cover 

    real (kind=dbl_kind), parameter :: &
       cp_sno = 0.0_dbl_kind & ! specific heat of snow                (J/kg/K)
    ,  rLfi = Lfresh*rhoi & ! latent heat of fusion ice               (J/m^3)
    ,  rLfs = Lfresh*rhos & ! latent heat of fusion snow              (J/m^3)
    ,  rLvi = Lvap*rhoi   & ! latent heat of vapor*rhoice             (J/m^3)
    ,  rLvs = Lvap*rhos   & ! latent heat of vapor*rhosno             (J/m^3)
    ,  rcpi = cp_ice*rhoi & ! heat capacity of fresh ice              (J/m^3)
    ,  rcps = cp_sno*rhos & ! heat capacity of snow                   (J/m^3)
    ,  rcpidepressT = rcpi*depressT & ! param for finding T(z) from q (J/m^3)
    ,  rLfidepressT = rLfi*depressT ! param for heat capacity   (J deg/m^3)
         ! heat capacity of sea ice, rhoi*C=rcpi+rLfidepressT*salinity/T^2

!=======================================================================
contains
!===============================================================================
!BOP
!
! !IROUTINE: ice_prescribed_init -  prescribed ice initialization
!
! !INTERFACE: 
 subroutine ice_prescribed_init(compid, gsmap, dom)
   use shr_pio_mod, only : shr_pio_getiotype, shr_pio_getiosys
! !DESCRIPTION:
!    Prescribed ice initialization - needed to 
!    work with new shr_strdata module derived type 
!
! !REVISION HISTORY:
!    2009-Oct-12 - M. Vertenstein
!
! !INPUT/OUTPUT PARAMETERS:
!
   implicit none
   include 'mpif.h'
   integer(kind=int_kind), intent(in) :: compid
   type(mct_gsMap) :: gsmap
   type(mct_gGrid) :: dom

!EOP
   !----- Local ------
   integer(kind=int_kind) :: nml_error ! namelist i/o error flag
   integer(kind=int_kind) :: n, nFile, ierr   
   character(len=8)       :: fillalgo
   character(*),parameter :: subName = "('ice_prescribed_init2')"
   character(*),parameter :: F00 = "('(ice_prescribed_init2) ',4a)"

   namelist /ice_prescribed_nml/  &
        prescribed_ice,      &
        model_year_align,    &
	stream_year_first ,  &
        stream_year_last  ,  &
        stream_fldVarName ,  &
        stream_fldFileName,  &
        stream_domTvarName,  &
        stream_domXvarName,  &
        stream_domYvarName,  &
        stream_domAreaName,  &
        stream_domMaskName,  &
        stream_domFileName,  &
        stream_mapread,      &
        prescribed_ice_fill

   ! default values for namelist
   prescribed_ice         = .false.          ! if true, prescribe ice
   stream_year_first      = 1                ! first year in  pice stream to use
   stream_year_last       = 1                ! last  year in  pice stream to use
   model_year_align       = 1                ! align stream_year_first with this model year
   stream_fldVarName      = 'ice_cov'
   stream_fldFileName(:)  = ' '
   stream_domTvarName     = 'time'
   stream_domXvarName     = 'lon'
   stream_domYvarName     = 'lat'
   stream_domAreaName     = 'area'
   stream_domMaskName     = 'mask'
   stream_domFileName     = ' '
   stream_mapread         = 'NOT_SET'
   prescribed_ice_fill    = .false.          ! true if pice data fill required

   ! read from input file
   call get_fileunit(nu_nml)
   if (my_task == master_task) then
      open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nu_nml, nml=ice_prescribed_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nu_nml)
   endif
   call release_fileunit(nu_nml)
   call broadcast_scalar(nml_error,master_task)
   if (nml_error /= 0) then
      call abort_ice ('ice: Namelist read error in ice_prescribed_mod')
   endif

   call broadcast_scalar(prescribed_ice,master_task)

   ! *** If not prescribed ice then return ***
   if (.not. prescribed_ice) RETURN

   call broadcast_scalar(model_year_align,master_task)
   call broadcast_scalar(stream_year_first,master_task)
   call broadcast_scalar(stream_year_last,master_task)
   call broadcast_scalar(stream_fldVarName,master_task)
   call broadcast_scalar(stream_domTvarName,master_task)
   call broadcast_scalar(stream_domXvarName,master_task)
   call broadcast_scalar(stream_domYvarName,master_task)
   call broadcast_scalar(stream_domAreaName,master_task)
   call broadcast_scalar(stream_domMaskName,master_task)
   call broadcast_scalar(stream_domFileName,master_task)
   call broadcast_scalar(stream_mapread,master_task)
   call broadcast_scalar(prescribed_ice_fill,master_task)
   call mpi_bcast(stream_fldFileName, len(stream_fldFileName(1))*NFilesMaximum, &
        MPI_CHARACTER, 0, MPI_COMM_ICE, ierr)

   nFile = 0
   do n=1,nFilesMaximum
      if (stream_fldFileName(n) /= ' ') nFile = nFile + 1
   end do

   ! Read shr_strdata_nml namelist
   if (prescribed_ice_fill) then
      fillalgo='nn'
   else
      fillalgo='none'
   endif

   if (my_task == master_task) then
      write(nu_diag,*) ' '
      write(nu_diag,*) 'This is the prescribed ice coverage option.'
      write(nu_diag,*) '  stream_year_first  = ',stream_year_first  
      write(nu_diag,*) '  stream_year_last   = ',stream_year_last   
      write(nu_diag,*) '  model_year_align   = ',model_year_align   
      write(nu_diag,*) '  stream_fldVarName  = ',trim(stream_fldVarName)
      do n = 1,nFile
         write(nu_diag,*) '  stream_fldFileName = ',trim(stream_fldFileName(n)),n
      end do
      write(nu_diag,*) '  stream_domTvarName = ',trim(stream_domTvarName)
      write(nu_diag,*) '  stream_domXvarName = ',trim(stream_domXvarName)
      write(nu_diag,*) '  stream_domYvarName = ',trim(stream_domYvarName)
      write(nu_diag,*) '  stream_domFileName = ',trim(stream_domFileName)
      write(nu_diag,*) '  stream_mapread     = ',trim(stream_mapread)
      write(nu_diag,*) '  stream_fillalgo    = ',trim(fillalgo)
      write(nu_diag,*) ' '
   endif

   call shr_strdata_create(sdat,name="prescribed_ice", &
        mpicom=MPI_COMM_ICE, compid=compid, &
        gsmap=gsmap, ggrid=dom,          &
        nxg=nx_global,nyg=ny_global,     &
        yearFirst=stream_year_first,     &
        yearLast=stream_year_last,       &
        yearAlign=model_year_align,      &
        offset=0,                        &
        domFilePath='',                  &
        domFileName=trim(stream_domFileName), &
        domTvarName=stream_domTvarName,  &
        domXvarName=stream_domXvarName,  &
        domYvarName=stream_domYvarName,  &
        domAreaName=stream_domAreaName,  &
        domMaskName=stream_domMaskName,  &
        filePath='',                     &
        filename=stream_fldFileName(1:nFile), &
        fldListFile=stream_fldVarName,   &
        fldListModel=stream_fldVarName,  &
        pio_subsystem=shr_pio_getiosys(inst_name), &
        pio_iotype=shr_pio_getiotype(inst_name),   &
        fillalgo=trim(fillalgo),       &
        calendar=trim(calendar_type),  &
        mapread=trim(stream_mapread))

   if (my_task == master_task) then
      call shr_strdata_print(sdat,'SPRESICE data')
   endif

   !-----------------------------------------------------------------
   ! For one ice category, set hin_max(1) to something big
   !-----------------------------------------------------------------
   if (ncat == 1) then
      hin_max(1) = 999._dbl_kind
   end if
end subroutine ice_prescribed_init
  
!=======================================================================
!BOP ===================================================================
!
! !IROUTINE: ice_prescribed_run -- Update ice coverage
!
! !DESCRIPTION:
!
!  Finds two time slices bounding current model time, remaps if necessary
!
! !REVISION HISTORY:
!     2005-May-19 - J. Schramm - first version
!     2009-Oct-15 - M. Vertenstein - update to new data model changes
!
! !INTERFACE: -----------------------------------------------------------

subroutine ice_prescribed_run(mDateIn, secIn)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(kind=int_kind), intent(in) :: mDateIn  ! Current model date (yyyymmdd)
   integer(kind=int_kind), intent(in) :: secIn    ! Elapsed seconds on model date

!EOP

   integer(kind=int_kind) :: i,j,n,iblk       ! loop indices and counter
   integer(kind=int_kind) :: ilo,ihi,jlo,jhi  ! beginning and end of physical domain
   type (block)           :: this_block
   real(kind=dbl_kind)    :: aice_max         ! maximun ice concentration
   logical, save          :: first_time = .true.
   character(*),parameter :: subName = "('ice_prescribed_run')"
   character(*),parameter :: F00 = "('(ice_prescribed_run) ',a,2g20.13)"
 
   !------------------------------------------------------------------------
   ! Interpolate to new ice coverage
   !------------------------------------------------------------------------

   call shr_strdata_advance(sdat,mDateIn,SecIn,MPI_COMM_ICE,'cice_pice')
   
   ice_cov(:,:,:) = c0  ! This initializes ghost cells as well 

   n=0
   do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      
      do j = jlo, jhi
      do i = ilo, ihi
         n = n+1
         ice_cov(i,j,iblk) = sdat%avs(1)%rAttr(1,n)
      end do
      end do
   end do

   !--------------------------------------------------------------------
   ! Check to see that ice concentration is in fraction, not percent
   !--------------------------------------------------------------------
   if (first_time) then
      aice_max = maxval(ice_cov)

      if (aice_max > c10) then
         write(nu_diag,F00) "ERROR: Ice conc data must be in fraction, aice_max= ",&
              aice_max
         call abort_ice(subName)
      end if
      first_time = .false.
   end if

   !-----------------------------------------------------------------
   ! Set prescribed ice state and fluxes
   !-----------------------------------------------------------------

   call ice_prescribed_phys()

end subroutine ice_prescribed_run

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ice_prescribed_phys -- set prescribed ice state and fluxes
!
! !DESCRIPTION:
!
! Set prescribed ice state using input ice concentration;
! set surface ice temperature to atmospheric value; use
! linear temperature gradient in ice to ocean temperature.
!
! !REVISION HISTORY:
!     2005-May-23 - J. Schramm - Updated with data models
!     2004-July   - J. Schramm - Modified to allow variable snow cover
!     2001-May    - B. P. Briegleb - Original version
!
! !INTERFACE: ------------------------------------------------------------------
 
subroutine ice_prescribed_phys

! !USES:
 
   use ice_flux
!  use ice_grid, only : bound
   use ice_state
   use ice_itd, only  : aggregate
   use ice_dyn_evp

   implicit none
 
! !INPUT/OUTPUT PARAMETERS:
 
!EOP

   !----- Local ------
   integer(kind=int_kind) :: layer    ! level index
   integer(kind=int_kind) :: nc       ! ice category index
   integer(kind=int_kind) :: i,j,k    ! longitude, latitude and level indices
   integer(kind=int_kind) :: iblk

   real(kind=dbl_kind) :: slope     ! diff in underlying ocean tmp and ice surface tmp
   real(kind=dbl_kind) :: Ti        ! ice level temperature
   real(kind=dbl_kind) :: Tmlt      ! ice level melt temperature
   real(kind=dbl_kind) :: qin_save(nilyr) 
   real(kind=dbl_kind) :: qsn_save(nslyr)
   real(kind=dbl_kind) :: hi        ! ice prescribed (hemispheric) ice thickness
   real(kind=dbl_kind) :: hs        ! snow thickness
   real(kind=dbl_kind) :: zn        ! normalized ice thickness
   real(kind=dbl_kind) :: salin(nilyr)  ! salinity (ppt) 

   real(kind=dbl_kind), parameter :: nsal    = 0.407_dbl_kind
   real(kind=dbl_kind), parameter :: msal    = 0.573_dbl_kind
   real(kind=dbl_kind), parameter :: saltmax = 3.2_dbl_kind   ! max salinity at ice base (ppm)

   !-----------------------------------------------------------------
   ! Initialize ice state
   !-----------------------------------------------------------------

   ! TODO  - can we now get rid of the following???

   !  aicen(:,:,:,:) = c0
   !  vicen(:,:,:,:) = c0
   !  eicen(:,:,:,:) = c0
   
   !  do nc=1,ncat
   !     trcrn(:,:,nt_Tsfc,nc,:) = Tf(:,:,:)
   !  enddo
   
   !-----------------------------------------------------------------
   ! Set ice cover over land to zero, not sure if this should be
   ! be done earier, before time/spatial interp??????
   !-----------------------------------------------------------------
   do iblk = 1,nblocks
   do j = 1,ny_block
   do i = 1,nx_block
      if (tmask(i,j,iblk)) then
         if (ice_cov(i,j,iblk) .lt. eps04) ice_cov(i,j,iblk) = c0
         if (ice_cov(i,j,iblk) .gt. c1)    ice_cov(i,j,iblk) = c1
      else
         ice_cov(i,j,iblk) = c0
      end if
   enddo
   enddo
   enddo

   do iblk = 1,nblocks
   do j = 1,ny_block
   do i = 1,nx_block

      if (tmask(i,j,iblk)) then   ! Over ocean points

         !--------------------------------------------------------------
         ! Place ice where ice concentration > .0001
         !--------------------------------------------------------------
         if (ice_cov(i,j,iblk) >= eps04) then

            hi = 0.0_dbl_kind
            !----------------------------------------------------------
            ! Set ice thickness in each hemisphere
            !----------------------------------------------------------
            if(TLAT(i,j,iblk)*rad_to_deg > 40.0_dbl_kind) then
              hi  = 2.0_dbl_kind
            else if(TLAT(i,j,iblk)*rad_to_deg < -40.0_dbl_kind) then
              hi  = 1.0_dbl_kind
            end if

            !----------------------------------------------------------
            ! All ice in appropriate thickness category
            !----------------------------------------------------------
            do nc = 1,ncat

              if(hin_max(nc-1) < hi .and. hi < hin_max(nc)) then

                  if (aicen(i,j,nc,iblk) > c0) then
                     hs = vsnon(i,j,nc,iblk) / aicen(i,j,nc,iblk)
                  else
                     hs = c0
                  endif

                  qin_save(:) = c0
                  qsn_save(:) = c0

                  if (vicen(i,j,nc,iblk) > c0) then
                     do k=1,nilyr
                        qin_save(k) = eicen(i,j,ilyr1(nc)+k-1,iblk)         &
                                    * real(nilyr,kind=dbl_kind)             &
                                    / Vicen(i,j,nc,iblk)
                     enddo
                  endif

                  if (vsnon(i,j,nc,iblk) > c0) then
                     do k=1,nslyr
                        qsn_save(k) = esnon(i,j,slyr1(nc)+k-1,iblk)         &
                                    * real(nslyr,kind=dbl_kind)             &
                                    / vsnon(i,j,nc,iblk)
                     enddo
                  endif

                  aicen(i,j,nc,iblk) = ice_cov(i,j,iblk)
                  vicen(i,j,nc,iblk) = hi*aicen(i,j,nc,iblk) 
                  vsnon(i,j,nc,iblk) = hs*aicen(i,j,nc,iblk) 

                  ! remember enthalpy profile to compute energy
                  do k=1,nilyr
                     eicen(i,j,ilyr1(nc)+k-1,iblk)                          &
                        = qin_save(k) * vicen(i,j,nc,iblk)                  &
                        / real(nilyr,kind=dbl_kind)
                  enddo

                  do k=1,nslyr
                     esnon(i,j,slyr1(nc)+k-1,iblk)                          &
                        = qsn_save(k) * vsnon(i,j,nc,iblk)                  &
                        / real(nslyr,kind=dbl_kind)
                  enddo

                  !---------------------------------------------------------
                  ! make linear temp profile and compute enthalpy
                  !---------------------------------------------------------

                  if (abs(eicen(i,j,ilyr1(nc),iblk)) < puny) then

                  if (aice(i,j,iblk) < puny) &
                     trcrn(i,j,nt_Tsfc,nc,iblk) = Tf(i,j,iblk)

                  slope = Tf(i,j,iblk) - trcrn(i,j,nt_Tsfc,nc,iblk)
                  do k = 1, nilyr
                     zn = (real(k,kind=dbl_kind)-p5) / real(nilyr,kind=dbl_kind)
                     Ti = trcrn(i,j,nt_Tsfc,nc,iblk) + slope*zn
                     salin(k) = (saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
                     Tmlt = -salin(k)*depressT
                     eicen(i,j,ilyr1(nc)+k-1,iblk) =                        &
                       -(rhoi * (cp_ice*(Tmlt-Ti) &
                       + Lfresh*(c1-Tmlt/Ti) - cp_ocn*Tmlt)) &
                       * vicen(i,j,nc,iblk)/real(nilyr,kind=dbl_kind)
                  enddo

                  do k=1,nslyr
                     esnon(i,j,slyr1(nc)+k-1,iblk) =                       &
                        -rhos*(Lfresh - cp_ice*trcrn(i,j,nt_Tsfc,nc,iblk)) &
                         *vsnon(i,j,nc,iblk)
                  enddo

                  endif  ! aice < puny
               end if    ! hin_max
            enddo        ! ncat
         else
            trcrn(i,j,nt_Tsfc,:,iblk) = Tf(i,j,iblk)
            aicen(i,j,:,iblk) = c0
            vicen(i,j,:,iblk) = c0
            vsnon(i,j,:,iblk) = c0
            esnon(i,j,:,iblk) = c0
            eicen(i,j,:,iblk) = c0
         end if          ! ice_cov >= eps04
      end if             ! tmask
   enddo                 ! i
   enddo                 ! j

   !--------------------------------------------------------------------
   ! compute aggregate ice state and open water area
   !--------------------------------------------------------------------
   call aggregate (nx_block, ny_block,                      &
                   aicen(:,:,:,iblk),  trcrn(:,:,:,:,iblk), &
                   vicen(:,:,:,iblk),  vsnon(:,:,:,iblk),   &
                   eicen(:,:,:,iblk),  esnon(:,:,:,iblk),   &
                   aice(:,:,iblk),     trcr(:,:,:,iblk),    &
                   vice(:,:,iblk),     vsno(:,:,iblk),      &
                   eice(:,:,iblk),     esno(:,:,iblk),      &
                   aice0(:,:,iblk),    tmask(:,:,iblk),     &
                   ntrcr,              trcr_depend) 

   enddo                 ! iblk

   do iblk = 1, nblocks
   do j = 1, ny_block
     do i = 1, nx_block
       aice_init(i,j,iblk) = aice(i,j,iblk)
     enddo
   enddo
   enddo

   !--------------------------------------------------------------------
   ! set non-computed fluxes, ice velocities, ice-ocn stresses to zero
   !--------------------------------------------------------------------

   frzmlt    (:,:,:) = c0
   uvel      (:,:,:) = c0
   vvel      (:,:,:) = c0
   strocnxT  (:,:,:) = c0
   strocnyT  (:,:,:) = c0

   !-----------------------------------------------------------------
   ! other atm and ocn fluxes
   !-----------------------------------------------------------------
   call init_flux_atm
   call init_flux_ocn

end subroutine ice_prescribed_phys

!==============================================================================

end module ice_prescribed_mod

!==============================================================================
