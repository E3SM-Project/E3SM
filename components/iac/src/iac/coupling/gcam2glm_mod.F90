#define DEBUG
Module gcam2glm_mod
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: gcam2glm_mod
!
!  Interface of the integrated assessment component in CCSM
!
! !DESCRIPTION:
!
! !USES:

  use iac_fields_mod
  use shr_file_mod, only: shr_file_getunit, shr_file_freeunit
  use shr_cal_mod
  use shr_sys_mod
  use shr_kind_mod, only : r8 => shr_kind_r8,r4 => shr_kind_r4
  use netcdf
  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: gcam2glm_init_mod               ! clm initialization
  public :: gcam2glm_run_mod                ! clm run phase
  public :: gcam2glm_final_mod              ! clm finalization/cleanup
  public :: handle_err
  public :: fround
  private :: D_mrgrnk

! !PUBLIC DATA MEMBERS: 
    integer, parameter :: kdp = selected_real_kind(15)
    real(r8), dimension(:, :), allocatable :: hydeGCROP2005,&
         hydeGPAST2005,   &
         hydeGWH2005,     &
         hydeGOTHR2005,   &
         cellarea,        &
         cellarea_forest, &
         cellarea_nonforest,&
         glm_crop_ann,    &
         glm_past_ann,    &
         glm_othr_ann,    &
         aez_regions,     &
         aez_zones,       &
         fnfforest,       &
         fnfnonforest,    &
         pot_veg,         &
         pot_veg_rev,         &
         crop_area,       &
         avail_land0,     &
         avail_landA,     &
         gcam_past,       &
         gcam_forest_area,&
         gcam_crop,       &
         gcam_farea

    real(r8), dimension(:), allocatable :: unmet_neg_past,unmet_neg_crop,unmet_farea, &
         cumsum_sorted_farea

    real(r8), dimension(:, :), allocatable,save :: gcam_wh, &
         pctland_in2005,rAEZ_sites,rAEZ_sites_rev,sortsitesup,sortsitesdn

    integer, dimension(:,:), allocatable      ::  raezs
    real(r8), dimension(:, :, :), allocatable,save :: glm_crop, &
         glm_past

    integer,save :: year1,year2

    real(r8), dimension(:), allocatable :: datearr, gcam_wh_ann, glm_wh_ann,lon,lat
    real(r8), dimension(:), pointer :: array1d
    real(r8)  :: forested_past, forested_past_percent
    real(r8)  :: nonforested_past, nonforested_past_percent
    real(r8)  :: forested_crop, forested_crop_percent
    real(r8)  :: nonforested_crop, nonforested_crop_percent

    real(r4)  :: miss_val = 1.0e36

    integer :: n,np1,nflds,gcamsize
    integer, dimension(3) :: start3,count3
    integer, dimension(2) :: start2,count2
    integer                               :: ncid,tmp(1), &
                                               lonDimID, latDimId, timeDimId, &
                                               numLons, numLats, numTimes,    &
                                               status,GCROPVarId,timevarid,varid
    integer, dimension(nf90_max_var_dims) :: dimIDs
    character(len=*),parameter :: gcam2glm_restfile = 'gcam2glm_restart.'
    character(len=*),parameter :: gcam2glm_rpointer = 'rpointer.gcam2glm'
    character(256) :: filename

! !REVISION HISTORY:
! Author: T Craig


! !PRIVATE DATA MEMBERS:

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam2glm_init_mod

! !INTERFACE:
  subroutine gcam2glm_init_mod( EClock, cdata, gcamo, glmi, glmi_wh)

! !DESCRIPTION:
! Initialize interface for glm

! !USES:
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real(r8), pointer :: gcamo(:,:)
    real(r8), pointer :: glmi(:,:)
    real(r8), pointer :: glmi_wh(:)




! !LOCAL VARIABLES:
    logical :: restart_run,lexist
    logical :: initial_run
    integer :: iu,iun,tmpyears(2)
    character(len=*),parameter :: subname='(gcam2glm_init_mod)'
    character(len=512) :: gcam2glm_basecrop
    character(len=512) :: gcam2glm_basepast
    character(len=512) :: gcam2glm_baseothr
    character(len=512) :: gcam2glm_aezmap
    character(len=512) :: gcam2glm_basebiomass

! !REVISION HISTORY:
! Author: T Craig
! Author: JET          ! rewrite of matlat preprocessing script new_grids_matchforest4.m

!EOP
!-----------------------------------------------------------------------

    iu  = cdata%i(iac_cdatai_logunit)
#ifdef DEBUG
     write(iu,*) subname,' starting subroutine '
#endif
    restart_run  = cdata%l(iac_cdatal_rest)
    gcamsize = cdata%i(iac_cdatai_gcam_naez)*cdata%i(iac_cdatai_gcam_nreg)
    initial_run = cdata%l(iac_cdatal_initrun)
    gcam2glm_basecrop = trim(cdata%c(iac_cdatac_gcam2glm_basecrop))
    gcam2glm_basepast = trim(cdata%c(iac_cdatac_gcam2glm_basepast))
    gcam2glm_baseothr = trim(cdata%c(iac_cdatac_gcam2glm_baseothr))
    gcam2glm_aezmap   = trim(cdata%c(iac_cdatac_gcam2glm_aezmap))
    gcam2glm_basebiomass = trim(cdata%c(iac_cdatac_gcam2glm_basebiomass))

! initialize two level time indexes 

    n=1
    np1=2

! Variable to keep track of when to calculate the next GCAM time step.
! Initialize this with the year passed in via EClock ... 2005 for now

    status= nf90_open(gcam2glm_basecrop,nf90_nowrite,ncid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inq_varid(ncid, "cropland", GCROPVarId)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inquire_variable(ncid, GCROPVarId, dimids = dimIDs)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inquire_dimension(ncid, dimIDs(1), len = numLons)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inquire_dimension(ncid, dimIDs(2), len = numLats)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inquire_dimension(ncid, dimIDs(3), len = numTimes)
    if(status /= nf90_NoErr) call handle_err(status)

    allocate(lon(numLons))
    allocate(lat(numLats))
    
    status = nf90_inq_varid(ncid, "lon", varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid,varid,lon)
    if(status /= nf90_NoErr) call handle_err(status)
    
    status = nf90_inq_varid(ncid, "lat", varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid,varid,lat)
    if(status /= nf90_NoErr) call handle_err(status)
    
    allocate(hydeGCROP2005(numLons, numLats))
    allocate(hydeGPAST2005(numLons, numLats))
    allocate(hydeGOTHR2005(numLons, numLats))
    allocate(hydeGWH2005(numLons, numLats))
    allocate(cellarea(numLons, numLats))
    allocate(cellarea_forest(numLons, numLats))
    allocate(cellarea_nonforest(numLons, numLats))
    allocate(glm_crop_ann(numLons, numLats))
    allocate(glm_past_ann(numLons, numLats))
    allocate(glm_othr_ann(numLons, numLats))
    allocate(cumsum_sorted_farea(numLons*numLats))
    allocate(glm_wh_ann(gcamsize))
    allocate(aez_regions(numLons, numLats))
    allocate(aez_zones(numLons, numLats))
    allocate(fnfforest(numLons, numLats))
    allocate(fnfnonforest(numLons, numLats))
    allocate(pot_veg(numLons, numLats))
    allocate(pot_veg_rev(numLats, numLons))
    allocate(crop_area(numLons, numLats))
    allocate(pctland_in2005(numLons, numLats))
    allocate(sortsitesup(numLons, numLats))
    allocate(sortsitesdn(numLons, numLats))
    allocate(rAEZ_sites(numLons, numLats))
    allocate(rAEZ_sites_rev(numLats, numLons))
    allocate(raezs(numLons, numLats))
    allocate(datearr(numTimes))
    allocate(glm_crop(numLons, numLats, 2))
    allocate(glm_past(numLons, numLats, 2))
    allocate(gcam_crop(gcamsize, 2))
    allocate(gcam_wh(gcamsize, 2))
    allocate(gcam_past(gcamsize, 2))
    allocate(gcam_forest_area(gcamsize, 2))
    allocate(unmet_neg_past(gcamsize))
    allocate(unmet_neg_crop(gcamsize))
    allocate(unmet_farea(gcamsize))
    allocate(avail_land0(numLons, numLats))
    allocate(avail_landA(numLons, numLats))

    glm_crop=iac_spval
    glm_past=iac_spval
    avail_land0=iac_spval
    avail_landA=iac_spval
    hydeGCROP2005=iac_spval
    hydeGPAST2005=iac_spval
    hydeGOTHR2005=iac_spval
    hydeGWH2005=iac_spval
    cellarea=iac_spval
    cellarea_forest=iac_spval
    cellarea_nonforest=iac_spval
    pctland_in2005=iac_spval
    datearr=iac_spval


    status = nf90_inq_varid(ncid,'time',timeVarId)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid,timeVarId,datearr)
    if(status /= nf90_NoErr) call handle_err(status)
    start3(1)=1
    count3(1)=numLons
    start3(2)=1
    count3(2)=numLats
    tmp=MAXLOC(datearr,mask=datearr.EQ.2005)
    start3(3)=tmp(1)
    count3(3)=1
    start2=1
    count2(1)=numLons
    count2(2)=numLats
! read in hyde data

    status = nf90_inq_varid(ncid,'cropland',varid)
    status = nf90_get_var(ncid,varid,hydeGCROP2005,start3,count3)
    status = nf90_close(ncid)
    status = nf90_open(gcam2glm_basepast,nf90_nowrite,ncid)
    status = nf90_inq_varid(ncid,'pasture',varid)
    status = nf90_get_var(ncid,varid,hydeGPAST2005,start3,count3)
    status = nf90_close(ncid)

    status = nf90_open(gcam2glm_baseothr,nf90_nowrite,ncid)
    status = nf90_inq_varid(ncid,'primary',varid)
    status = nf90_get_var(ncid,varid,hydeGOTHR2005,start3,count3)
    status = nf90_inq_varid(ncid,'cell_area',varid)
    status = nf90_get_var(ncid,varid,cellarea,start3,count3)
    status = nf90_close(ncid)
    cellarea=cellarea/1.e6
    status = nf90_open(gcam2glm_aezmap,nf90_nowrite,ncid)
    status = nf90_inq_varid(ncid,'aez_regions',varid)
    status = nf90_get_var(ncid,varid,aez_regions,start2,count2)
    status = nf90_inq_varid(ncid,'aez_zones',varid)
    status = nf90_get_var(ncid,varid,aez_zones,start2,count2)
    status = nf90_close(ncid)

    status = nf90_open(gcam2glm_basebiomass,nf90_nowrite,ncid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inq_varid(ncid,'biomass',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid,varid,pot_veg,start2,count2)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_close(ncid)
    if(status /= nf90_NoErr) call handle_err(status)
    
    cellarea_nonforest(:,:)=0.
    cellarea_forest(:,:)=0.
    pot_veg=pot_veg*0.75
    pctland_in2005(:,:) = 0.    
    pctland_in2005=hydeGCROP2005+hydeGPAST2005+hydeGOTHR2005

    if (.not.restart_run) then
       glm_crop(:,:,1)=hydeGCROP2005;
       glm_past(:,:,1)=hydeGPAST2005;
    else
       ! read restart and set crop and past
       
       inquire(file=trim(gcam2glm_rpointer),exist=lexist)
       if (lexist) then
          
#ifdef DEBUG
           write(iu,*) subname,' read_restart rpointer ',trim(gcam2glm_rpointer)
#endif
          
          iun = shr_file_getunit()
          open(iun,file=trim(gcam2glm_rpointer),form='formatted')
          read(iun,'(a)') filename
          close(iun)
          call shr_file_freeunit(iun)
          
#ifdef DEBUG
           write(iu,*) subname,' read_restart file ',trim(filename)
#endif
          
          inquire(file=trim(filename),exist=lexist)
          if (.not.lexist) then
             write(iu,*) subname,' ERROR: missing file ',trim(filename)
             call shr_sys_abort(subname//' ERROR: missing file')
          endif
          
          status= nf90_open(filename,nf90_nowrite,ncid)
          if(status /= nf90_NoErr) call handle_err(status)

	  ! Need to set cdata year1 and year 2 for restart

          status = nf90_inq_varid(ncid,'gcam_years',varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid,varid,tmpyears)
          if(status /= nf90_NoErr) call handle_err(status)
	  cdata%i(iac_cdatai_gcam_yr1)=tmpyears(1)
	  cdata%i(iac_cdatai_gcam_yr2)=tmpyears(2)
          
          status = nf90_inq_varid(ncid,'glm_crop',varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid,varid,glm_crop)
          if(status /= nf90_NoErr) call handle_err(status)
          
          status = nf90_inq_varid(ncid,'glm_past',varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid,varid,glm_past)
          if(status /= nf90_NoErr) call handle_err(status)
       else
          write(iu,*) subname,' read_restart rpointer NOT found ',trim(gcam2glm_rpointer)
          call shr_sys_abort(subname//' ERROR: missing file')
       end if ! rpointer exist
    end if ! restart run
  end subroutine gcam2glm_init_mod

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam2glm_run_mod

! !INTERFACE:
  subroutine gcam2glm_run_mod( EClock, cdata, gcamo, glmi, glmi_wh)

! !DESCRIPTION:
! Run interface for glm

! !USES:
    implicit none

! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real(r8), pointer :: gcamo(:,:)
    real(r8), pointer :: glmi_wh(:)
    real(r8), pointer :: glmi(:,:)

! !PARAMETERS:

    character(len=*),parameter :: subname='(gcam2glm_run_mod)'

    real(r8), parameter :: crop_forest_abandon_percent = 0.9
    real(r8), parameter :: past_forest_abandon_percent = 0.9

! !LOCAL VARIABLES:

    character*4 :: yearc
    character(256) :: filename
    integer :: i,j,ij,r,i1,j1,aez,ind,h,z
    integer :: iu,iun,iyr
    integer :: ymd, tod, dt,naez,nreg,ii,ntime,year,mon,day
    logical :: restart_now,gcam_alarm
    real(r8)  :: crop_d,past_d,crop_neg,crop_pos,past_neg,past_pos,farea_d
    real(r8)  :: gcam_crop_tmp(2,18,14),gcam_past_tmp(2,18,14),gcam_forest_area_tmp(2,18,14)
    real(r8)  :: fact1,fact2,eclockyr,fact1yrm1, fact2yrm1,delyr,eclockyrm1
    real(r8)  :: tmp0
    integer,save :: ncid,varid,dimid,dimid3(3)
    integer :: totraezs
    integer, allocatable  ::indxdn(:),indxup(:),sortlatsup(:),sortlonsup(:),sortlatsdn(:),sortlonsdn(:),indxa(:),indxadn(:),v1u(:),v2u(:),v1d(:),v2d(:)
    real(r8), allocatable   :: tmparr(:)
    integer :: sortind(1),max_aez_ind(1),regional_unmet_reassign
    real(r8)  :: crop_pos_f,crop_pos_nf,pot_forest_to_crop,pot_forest_from_crop,reassign_crop, &
               crop_neg_f,crop_neg_nf,crop_nfarea,crop_farea, &
               past_pos_f,past_pos_nf,pot_forest_to_past,pot_forest_from_past,reassign_past, &
               past_neg_f,past_neg_nf,past_nfarea,past_farea,GLM_nfarea,GLM_farea,regional_farea_needed, &
               crop_after_decrease,past_after_decrease,crop_before_decrease,past_before_decrease, &
               total_ag_decrease,crop_decrease_ratio,past_decrease_ratio,ag_area_avail,crop_area_avail, &
               final_area_needed,f_diff,reassign_ag_at_max_aez_ind,sumavail_land0,sumavail_landA

   real(r8), allocatable    ::  avail_farea(:),avail_nfarea(:),avail_ag_farea(:),reassign_ag(:), &
                              unmet_aez_farea(:),cumsum_sorted_reassign_ag(:),unmet_regional_farea(:)
   integer ntimes,nntimes,zz
! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------
    regional_unmet_reassign=1
    iu  = cdata%i(iac_cdatai_logunit)
    ymd = EClock(iac_EClock_ymd)
    tod = EClock(iac_EClock_tod)
    dt  = EClock(iac_EClock_dt)
    gcam_alarm=(EClock(iac_EClock_Agcam)==1)
#ifdef DEBUG
    write(iu,*) trim(subname),' date= ',ymd,tod
#endif

    year1=cdata%i(iac_cdatai_gcam_yr1)
    year2=cdata%i(iac_cdatai_gcam_yr2)
    
    naez=cdata%i(iac_cdatai_gcam_naez)
    nreg=cdata%i(iac_cdatai_gcam_nreg)
    ntime=cdata%i(iac_cdatai_gcamo_ntime)

    allocate(indxup(numLons*numLats))
    allocate(indxdn(numLons*numLats))
    allocate(sortlatsup(numLons*numLats))
    allocate(sortlatsdn(numLons*numLats))
    allocate(sortlonsup(numLons*numLats))
    allocate(sortlonsdn(numLons*numLats))
    allocate(tmparr(numLons*numLats))
    allocate(indxa(naez))
    allocate(indxadn(naez))
    allocate(avail_farea(naez))
    allocate(avail_nfarea(naez))
    allocate(avail_ag_farea(naez))
    allocate(reassign_ag(naez))
    allocate(unmet_aez_farea(naez))
    allocate(cumsum_sorted_reassign_ag(naez))
    allocate(unmet_regional_farea(nreg))

    avail_farea=0.
    avail_nfarea=0.
    avail_ag_farea=0.
    reassign_ag=0.
    unmet_aez_farea=0.
    cumsum_sorted_reassign_ag=0.
    unmet_regional_farea=0.

    unmet_farea=0.

    if (gcam_alarm) then
    where (  aez_regions  >= 1 .and.   aez_regions <= 14)
       glm_crop(:,:,np1)=0
       glm_past(:,:,np1)=0
    elsewhere
       glm_crop(:,:,np1)=hydeGCROP2005;
       glm_past(:,:,np1)=hydeGPAST2005;
    end where
#ifdef DEBUG
    write(6,*) 'sum 0 glm_crop='
    write(6,fmt="(1ES25.15)") sum(glm_crop(:,:,np1))
    write(6,*) 'sum 0 glm_past='
    write(6,fmt="(1ES25.15)") sum(glm_past(:,:,np1))
#endif     
    ! Unpack gcamo field
    ij = 0
    do j = n,np1
       do i = 1,nreg
          do h = 1,naez
             ij = ij + 1
             gcam_crop((i-1)*naez+h,j) = gcamo(iac_gcamo_crop,ij)
             gcam_past((i-1)*naez+h,j) = gcamo(iac_gcamo_pasture,ij)
             gcam_wh((i-1)*naez+h,j) = gcamo(iac_gcamo_woodharv,ij)
             gcam_forest_area((i-1)*naez+h,j) = gcamo(iac_gcamo_forest,ij)
          enddo
       enddo
    enddo

#ifdef DEBUG
    write(iu,*) subname,' year1 year2 ',year1,year2
#endif
    call shr_sys_flush(iu)
    ind=0
    do r = 1,nreg
       do aez = 1,naez
          ind=ind+1
          crop_d = (gcam_crop(ind,np1)-gcam_crop(ind,n))*1000
          past_d = (gcam_past(ind,np1)-gcam_past(ind,n))*1000
          farea_d = (gcam_forest_area(ind,np1)-gcam_forest_area(ind,n))*1000

          ! convert area changes into positive and negative changes
          
          if (crop_d <= 0) then
             crop_neg = -1.*crop_d
             crop_pos = 0.
          else
             crop_neg = 0.
             crop_pos = crop_d
          end if
          if (past_d <= 0) then
             past_neg = -1.*past_d
             past_pos = 0.
          else
             past_neg = 0.
             past_pos = past_d
          end if
          
          ! compute the potentially forested and potentially non-forested
          ! area currently occupied by cropland or pasture  in each
          ! region/AEZ. Also compute the *actual* forested and non-forested area
          ! currently available in each region/AEZ        
          raezs=0
          where ( aez_regions .EQ. r .and. aez_zones .EQ. aez)
             raezs=1
             rAEZ_sites=1.
          elsewhere
             rAEZ_sites=0.
          end where

          cellarea_forest(:,:)=0.
          cellarea_nonforest(:,:)=0.
          fnfnonforest=0.
          fnfforest=0.
          where ( pot_veg > 1.)
             cellarea_forest(:,:)=cellarea(:,:)
             fnfforest(:,:)=1.
          elsewhere
             fnfnonforest(:,:)=1.
             cellarea_nonforest(:,:)=cellarea(:,:)
          end where
	  cellarea_nonforest=cellarea_nonforest*rAEZ_sites
	  cellarea_forest=cellarea_forest*rAEZ_sites
          fnfforest=fnfforest*rAEZ_sites
          fnfnonforest=fnfnonforest*rAEZ_sites

          crop_farea =sum(glm_crop(:,:,n)*cellarea_forest)
          crop_nfarea =sum(glm_crop(:,:,n)*cellarea_nonforest)
          past_farea =sum(glm_past(:,:,n)*cellarea_forest)
          past_nfarea =sum(glm_past(:,:,n)*cellarea_nonforest)
          GLM_nfarea = sum((pctland_in2005 - glm_crop(:,:,n) - glm_past(:,:,n))*cellarea_nonforest)
          GLM_farea = sum((pctland_in2005 - glm_crop(:,:,n) - glm_past(:,:,n))*cellarea_forest)
          totraezs=sum(raezs)
          if (allocated(v1u)) deallocate(v1u)
          if (allocated(v2u)) deallocate(v2u)
          if (allocated(v1d)) deallocate(v1d)
          if (allocated(v2d)) deallocate(v2d)
          allocate(v1u(totraezs))
          allocate(v2u(totraezs))
          allocate(v1d(totraezs))
          allocate(v2d(totraezs))
          indxup(:)=(/(i,i=1,numlons*numlats)/)
          indxdn=indxup
          pot_veg_rev=transpose(pot_veg)
          rAEZ_sites_rev=transpose(rAEZ_sites)
          call D_mrgrnk(pot_veg_rev*rAEZ_sites_rev,indxup,numlons*numlats)
          call D_mrgrnk(pot_veg_rev*rAEZ_sites_rev*-1.,indxdn,numlons*numlats)
          !	  
          !  The sortxxxup and sortxxxdn arrays are only good 1:totraezs these arrays are also based on arrays lat first (360,720)
          !  sorting to match original matlab scripts.  Thats why I transposed the arrays above.  Also didn't use qsort because it
          !  isn't a stable sort.
          !
          sortlonsdn=(indxdn-1)/numlats+1
          sortlatsdn=mod(indxdn-1,numlats)+1
          sortlonsup=(indxup-1)/numlats+1
          sortlatsup=mod(indxup-1,numlats)+1
          v1u=sortlonsup(numlons*numlats-totraezs+1:numlons*numlats)
          v2u=sortlatsup(numlons*numlats-totraezs+1:numlons*numlats)
          v1d=(sortlonsdn(:totraezs))
          v2d=(sortlatsdn(:totraezs))

          if ((past_neg >= (past_farea + past_nfarea)).and.(past_neg>0)) then
             unmet_neg_past(ind) = past_neg - (past_farea + past_nfarea)
             past_neg = past_farea + past_nfarea
          end if
          if ((crop_neg >= (crop_farea + crop_nfarea)).and.(crop_neg>0)) then
             unmet_neg_crop(ind) = crop_neg - (crop_farea + crop_nfarea)
             crop_neg = crop_farea + crop_nfarea
          end if
          
          ! compute the potential forest that could be created from cropland
          ! and pasture abandonment
          if (crop_neg <= crop_farea) then
             pot_forest_from_crop = crop_neg
          else
             pot_forest_from_crop = crop_farea
          end if
          
          if (past_neg <= past_farea) then
             pot_forest_from_past = past_neg
          else
             pot_forest_from_past = past_farea
          end if
          
          ! compute the potential forest that could be lost from cropland and
          ! pasture expansion
          if ((crop_pos>0).and.(crop_pos > GLM_nfarea)) then
             pot_forest_to_crop = crop_pos - GLM_nfarea
             ! QUESTION: what about using abandoned pasture if needed?
             ! (MODIFY IN FUTURE?)
          else
             pot_forest_to_crop = 0
          end if
          
          if ((past_pos>0).and.(past_pos > (GLM_nfarea - (crop_pos - pot_forest_to_crop)))) then
             pot_forest_to_past = past_pos - (GLM_nfarea - (crop_pos - pot_forest_to_crop))
             ! QUESTION: what about using abandoned cropland if needed?
             ! (MODIFY IN FUTURE)
          else
             pot_forest_to_past = 0
          end if
          
          ! if GCAM forest area not expanding, compute the crop and pasture
          ! changes that need to come from forest and non-forest
          if (farea_d <= 0) then
             ! farea not expanding
             ! MODIFY THIS IN FUTURE? SO THAT WE DON'T "OVERSHOOT" FAREA_D
             crop_pos_f = pot_forest_to_crop
             crop_pos_nf = crop_pos - pot_forest_to_crop
             crop_neg_f = pot_forest_from_crop
             crop_neg_nf = (crop_neg - pot_forest_from_crop)/(crop_nfarea+1e-12)
             past_pos_f = pot_forest_to_past
             past_pos_nf = past_pos - pot_forest_to_past
             past_neg_f = pot_forest_from_past
             past_neg_nf = (past_neg - pot_forest_from_past)/(past_nfarea+1e-12)
             
          elseif ((pot_forest_from_crop + pot_forest_from_past - pot_forest_to_crop - pot_forest_to_past) >= farea_d) then
             ! compute the areas of crop and past needed to move off forest
             
             ! enough forested area available just by abandoning
             ! cropland and pasture on pot. forest
             ! compute forest percentages etc
             
             ! MODIFY THIS IN FUTURE? ... currently doing max forest abaondon ... might need less
             ! WHAT ABOUT CROPLAND ABANDONMENT THAT ALLOWS PASTURE
             ! EXPANSION?
             crop_pos_f = pot_forest_to_crop
             crop_pos_nf = crop_pos - pot_forest_to_crop
             crop_neg_f = pot_forest_from_crop
             crop_neg_nf = (crop_neg - pot_forest_from_crop)/(crop_nfarea+1e-12)
             past_pos_f = pot_forest_to_past
             past_pos_nf = past_pos - pot_forest_to_past
             past_neg_f = pot_forest_from_past
             past_neg_nf = (past_neg - pot_forest_from_past)/(past_nfarea+1e-12)
             
          else
             ! in addition to abandoning crop and past on pot. forest,
             ! also reassign some crop and past to non-forest
             ! compute forest percentages and also reassignment
             ! percentages
             f_diff = farea_d - (pot_forest_from_crop + pot_forest_from_past - pot_forest_to_crop - pot_forest_to_past)
             if (((GLM_nfarea - (crop_pos - pot_forest_to_crop) -(past_pos - pot_forest_to_past))>=f_diff) .and. ((crop_farea+past_farea - pot_forest_from_crop - pot_forest_from_past)>=f_diff)) then
                reassign_crop = (crop_farea - pot_forest_from_crop)/(crop_farea + past_farea - pot_forest_from_crop - pot_forest_from_past)*f_diff
                reassign_past = (past_farea - pot_forest_from_past)/(crop_farea + past_farea - pot_forest_from_crop - pot_forest_from_past)*f_diff
                ! what if availability is mostly in crop or past?
             else
                reassign_crop = min(crop_farea - pot_forest_from_crop, GLM_nfarea - (crop_pos - pot_forest_to_crop) - (past_pos - pot_forest_to_past))
                reassign_past = min(past_farea - pot_forest_from_past, GLM_nfarea - (past_pos - pot_forest_to_past) - (crop_pos - pot_forest_to_crop) - reassign_crop)
                unmet_farea(ind) = f_diff - (reassign_crop + reassign_past)
             end if
             
             crop_pos_f = pot_forest_to_crop
             crop_pos_nf = crop_pos - pot_forest_to_crop + reassign_crop
             crop_neg_f = (pot_forest_from_crop + reassign_crop)
             crop_neg_nf = (crop_neg - pot_forest_from_crop)/(crop_nfarea+1e-12)
             past_pos_f = pot_forest_to_past
             past_pos_nf = past_pos - pot_forest_to_past + reassign_past
             past_neg_f = (pot_forest_from_past + reassign_past)
             past_neg_nf = (past_neg - pot_forest_from_past)/(past_nfarea+1e-12)
             
          end if
          
          ! apply the cropland and pature changes (both forested and non-forested)
          ! to individual gridcells in each region/AEZ
          
          ! crop_neg_nf and crop_neg_f
	  where (raez_sites > 0) 
             glm_crop(:,:,np1) = glm_crop(:,:,n)   -  glm_crop(:,:,n)*fnfnonforest*crop_neg_nf
          endwhere
   !jt             [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'descend')
          if (crop_neg_f>0) then 
             cumsum_sorted_farea=0.
             call cumsum(glm_crop(:,:,n)*cellarea_forest,v1d,v2d,cumsum_sorted_farea(:totraezs),totraezs)
             sortind=MINLOC(cumsum_sorted_farea(:totraezs),mask=cumsum_sorted_farea(:totraezs)>crop_neg_f)
             if (sortind(1)==0) then
                if ( abs(crop_neg_f-cumsum_sorted_farea(totraezs)) .lt. 1e-10) then 
                   crop_neg_f=cumsum_sorted_farea(totraezs)
                   sortind=MINLOC(cumsum_sorted_farea(:totraezs),mask=cumsum_sorted_farea(:totraezs)>=crop_neg_f)
                else
                   sortind(1) = totraezs
                end if
             end if
             if (sortind(1)>1) then
                sortsitesdn=0
                do i=1,sortind(1)-1
                   sortsitesdn(v1d(i),v2d(i))=1
                end do
                where(sortsitesdn>0)
                   glm_crop(:,:,np1)= glm_crop(:,:,np1) - glm_crop(:,:,n) * fnfforest
                end where
                final_area_needed =  crop_neg_f - cumsum_sorted_farea(sortind(1)-1)
             else
                final_area_needed =  crop_neg_f
             end if
             glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) = &
                  glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) - &
                  min(final_area_needed/cellarea(v1d(sortind(1)),v2d(sortind(1))), &
                  glm_crop(v1d(sortind(1)),v2d(sortind(1)),n)*fnfforest(v1d(sortind(1)),v2d(sortind(1))))
          end if
          
          ! past_neg_f and past_neg_nf
          
	  where(raez_sites > 0) 
             glm_past(:,:,np1) = glm_past(:,:,n) - glm_past(:,:,n)*fnfnonforest*past_neg_nf
          end where
          if (past_neg_f>0) then
             !jt             [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'descend')
             cumsum_sorted_farea=0.
             call cumsum(glm_past(:,:,n)*cellarea_forest(:,:),v1d,v2d,cumsum_sorted_farea(:totraezs),totraezs)
             sortind=MINLOC(cumsum_sorted_farea(:totraezs),mask=cumsum_sorted_farea(:totraezs)>past_neg_f)
             if (sortind(1)==0) then
                if ( abs(past_neg_f-cumsum_sorted_farea(totraezs)) .lt. 1e-10) then 
                   past_neg_f=cumsum_sorted_farea(totraezs)
                   sortind=MINLOC(cumsum_sorted_farea(:totraezs),mask=cumsum_sorted_farea(:totraezs)>=past_neg_f)
                else
                   sortind(1) = totraezs
                end if
             end if
             if (sortind(1)>1) then
                sortsitesdn=0
                do i=1,sortind(1)-1
                   sortsitesdn(v1d(i),v2d(i))=1
                end do
                where(sortsitesdn>0)
                   glm_past(:,:,np1) = glm_past(:,:,np1) - glm_past(:,:,n) * fnfforest(:,:)
                end where
                final_area_needed =  past_neg_f - cumsum_sorted_farea(sortind(1)-1)
             else
                final_area_needed =  past_neg_f
             end if
             glm_past(v1d(sortind(1)),v2d(sortind(1)),np1) = &
                  glm_past(v1d(sortind(1)),v2d(sortind(1)),np1) - &
                  min(final_area_needed/cellarea(v1d(sortind(1)),v2d(sortind(1))), &
                  glm_past(v1d(sortind(1)),v2d(sortind(1)),n)*fnfforest(v1d(sortind(1)),v2d(sortind(1))))
          end if
          
          ! crop_pos_nf
          avail_land0 = 0
          where ( aez_regions .EQ. r .and. aez_zones .EQ. aez.and.glm_crop(:,:,np1)>0.)
             avail_land0=(pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
          endwhere
          sumavail_land0=sum(avail_land0)
          
          if ( sumavail_land0 >= crop_pos_nf .or. abs(crop_pos_nf)<=1e-6) then
             if (abs(crop_pos_nf)<=1e-6) then
                crop_pos_nf=0.
             end if
#ifdef DEBUG
             write(6,*)'cropland increase on non-forested land - land available'
#endif
             where(raez_sites > 0) 
                glm_crop(:,:,np1) = glm_crop(:,:,np1) + (avail_land0/(sumavail_land0+1e-12)*crop_pos_nf*fnfnonforest)/cellarea
             end where
          else
             avail_landA = 0.
             where ( aez_regions .EQ. r .and. aez_zones .EQ. aez)
                avail_landA=(pctland_in2005 - glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
             end where
             sumavail_landA=sum(avail_landA)
             if (crop_pos_nf - sumavail_landA < 1e-6) then
#ifdef DEBUG
                write(6,*)'cropland increase - land available'
#endif 
                where(raez_sites > 0) 
                   glm_crop(:,:,np1) = glm_crop(:,:,np1) + avail_land0*fnfnonforest/cellarea
                   glm_crop(:,:,np1) = glm_crop(:,:,np1) + &
                        ((avail_landA-avail_land0)/(sumavail_landA-sumavail_land0+1e-12)*(crop_pos_nf-sumavail_land0)*fnfnonforest)/cellarea
                end where
             else
                write(6,*)'crop increase on non-forest - land not available'
                call abort
             end if
          end if
        
          ! past_pos_nf
          avail_land0 = 0
          where ( aez_regions .EQ. r .and. aez_zones .EQ. aez.and.glm_past(:,:,np1)>0.)
             avail_land0=(pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
          end where
          sumavail_land0=sum(avail_land0)
          if ( sumavail_land0 >= past_pos_nf .or. abs(past_pos_nf)<=1e-6) then
             if (abs(past_pos_nf)<=1e-6) then
                past_pos_nf=0
             end if
#ifdef DEBUG
             write(6,*)'pasture increase on non-forested land - land available'
#endif
             where(raez_sites > 0) 
                glm_past(:,:,np1) = glm_past(:,:,np1) + (avail_land0/(sumavail_land0+1e-12)*past_pos_nf*fnfnonforest)/cellarea
             end where
          else
             avail_landA = 0
             where ( aez_regions .EQ. r .and. aez_zones .EQ. aez)
                avail_landA=(pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
             end where
             sumavail_landA=sum(avail_landA)
             if ( sumavail_landA >= past_pos_nf) then
#ifdef DEBUG
                write(6,*)'pasture increase on non-forested land - land available'
#endif
                where(raez_sites > 0) 
                   glm_past(:,:,np1) = glm_past(:,:,np1) + avail_land0*fnfnonforest/cellarea
                   glm_past(:,:,np1) = glm_past(:,:,np1) + &
                        ((avail_landA-avail_land0)/(sumavail_landA -sumavail_land0+1e-12) * &
                        (past_pos_nf-sumavail_land0) * fnfnonforest)/cellarea
                end where
             else
#ifdef DEBUG
                write(6,*)'pasture increase on non-forest - land NOT available, reducing pasture increase to accomodate'
#endif
                past_pos_nf = sumavail_landA
                where(raez_sites > 0) 
                   glm_past(:,:,np1) = glm_past(:,:,np1) + (avail_land0*fnfnonforest/cellarea)
                   glm_past(:,:,np1) = glm_past(:,:,np1) + ((avail_landA - avail_land0) / &
                        (sumavail_landA - sumavail_land0 + 1e-12)*(past_pos_nf-sumavail_land0)*fnfnonforest)/cellarea
                end where
             end if
          end if
        
          ! crop_pos_f
          avail_land0 = 0
          where ( aez_regions .EQ. r .and. aez_zones .EQ. aez.and.glm_crop(:,:,np1)>0.)
             avail_land0=(pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
          end where
          sumavail_land0=sum(avail_land0)
          if (abs(crop_pos_f)>1e-6) then
             if (sumavail_land0>=crop_pos_f) then
#ifdef DEBUG
                write(6,*)'cropland increase on forested land - land available'
#endif
                !jt              [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
                cumsum_sorted_farea=0.
                call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totraezs),totraezs)
                sortind=MINLOC(cumsum_sorted_farea(:totraezs),mask=cumsum_sorted_farea(:totraezs)>crop_pos_f)
                if (sortind(1)==0) then
                   sortind(1) = totraezs
                end if
                if (sortind(1)>1) then
                   sortsitesup=0
                   do i=1,sortind(1)-1
                      sortsitesup(v1u(i),v2u(i))=1
                   end do
                   where(sortsitesup>0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + avail_land0 * fnfforest / cellarea
                   end where
                   final_area_needed =  crop_pos_f - cumsum_sorted_farea(sortind(1)-1)
                else
                   final_area_needed =  crop_pos_f
                end if
                glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                     glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) + &
                     min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))), &
                     avail_land0(v1u(sortind(1)),v2u(sortind(1))) * &
                     fnfforest(v1u(sortind(1)),v2u(sortind(1)))/cellarea(v1u(sortind(1)),v2u(sortind(1))))
             else
                avail_landA = 0.
                where ( aez_regions .EQ. r .and. aez_zones .EQ. aez)
                   avail_landA=(pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
                end where
                sumavail_landA=sum(avail_landA)
                if (sumavail_landA >=crop_pos_f) then
                   !jt  [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
                   cumsum_sorted_farea=0.
                   call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totraezs),totraezs)
                   where(raez_sites > 0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + avail_land0*fnfforest/cellarea
                   end where
#ifdef DEBUG
                   write(6,*)'cropland increase on forested land - land available'
#endif                   
                   !jt  [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
                   call cumsum(avail_landA(:,:)-avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totraezs),totraezs)
                   sortind = MINLOC(cumsum_sorted_farea(:totraezs),mask=cumsum_sorted_farea(:totraezs)>(crop_pos_f-sum(avail_land0(:,:),mask=raez_sites>0)))
                   if (sortind(1)==0) then
                      sortind(1) = totraezs
                   end if
                   if (sortind(1)>1) then
                      sortsitesup=0
                      do i=1,sortind(1)-1
                         sortsitesup(v1u(i),v2u(i))=1
                      end do
                      where(sortsitesup>0)
                         glm_crop(:,:,np1) = glm_crop(:,:,np1) + (avail_landA - avail_land0) * fnfforest / cellarea
                      end where
                      final_area_needed =  crop_pos_f - sumavail_land0-cumsum_sorted_farea(sortind(1)-1)
                   else
                      final_area_needed =  crop_pos_f - sumavail_land0
                   end if
                   glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                        glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) +&
                        min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))), &
                        (avail_landA(v1u(sortind(1)),v2u(sortind(1)))) * &
                        fnfforest(v1u(sortind(1)),v2u(sortind(1))) / &
                        cellarea(v1u(sortind(1)),v2u(sortind(1))))
                else
                   write(6,*)'crop increase on forest - land not available'
                   call abort
                end if
             end if
             
          end if
        
          ! past_pos_f
          avail_land0 = 0.
          where ( aez_regions .EQ. r .and. aez_zones .EQ. aez.and.glm_past(:,:,np1)>0.)
             avail_land0=(pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
          end where
          sumavail_land0=sum(avail_land0)
          if (abs(past_pos_f)>1e-6) then 
             if (sumavail_land0>=past_pos_f) then
#ifdef DEBUG
                write(6,*)'pasture increase on forest - land available'
#endif                
                !jt [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
                cumsum_sorted_farea=0.
                call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totraezs),totraezs)
                sortind = MINLOC(cumsum_sorted_farea(:totraezs),mask=cumsum_sorted_farea(:totraezs)>past_pos_f)
                if (sortind(1)==0) then
                   sortind(1) = totraezs
                end if
                if (sortind(1)>1) then 
                   sortsitesup=0
                   do i=1,sortind(1)-1
                      sortsitesup(v1u(i),v2u(i))=1
                   end do
                   where(sortsitesup>0)
                      glm_past(:,:,np1) = glm_past(:,:,np1) + avail_land0 * fnfforest / cellarea
                   end where
                   final_area_needed =  past_pos_f - cumsum_sorted_farea(sortind(1)-1)
                else
                   final_area_needed =  past_pos_f
                end if
                glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                     glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) + &
                     min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))) , &
                     avail_land0(v1u(sortind(1)),v2u(sortind(1))) * &
                     fnfforest(v1u(sortind(1)),v2u(sortind(1))) / &
                     cellarea(v1u(sortind(1)),v2u(sortind(1))))
             else
                avail_landA = 0.
                where ( aez_regions .EQ. r .and. aez_zones .EQ. aez)
                   avail_landA=(pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
                end where
                sumavail_landA=sum(avail_landA)
                if (sumavail_landA >= past_pos_f) then 
#ifdef DEBUG
                   write(6,*)'pasture increase on forest - land available'
#endif
                   !jt [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
                   cumsum_sorted_farea=0.
                   call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totraezs),totraezs)
                   where(raez_sites > 0)
                      glm_past(:,:,np1) = glm_past(:,:,np1) + avail_land0 * fnfforest / cellarea
                   end where
                   
                   !jt [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
                   
                   call cumsum(avail_landA(:,:)-avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totraezs),totraezs)
                   sortind = MINLOC(cumsum_sorted_farea(:totraezs),mask=cumsum_sorted_farea(:totraezs)>(past_pos_f-sum(avail_land0(:,:),mask=raez_sites>0)))
                   if (sortind(1)==0) then
                      sortind(1) = totraezs
                   end if
                   if (sortind(1)>1) then
                      sortsitesup=0
                      do i=1,sortind(1)-1
                         sortsitesup(v1u(i),v2u(i))=1
                      end do
                      where(sortsitesup>0)
                         glm_past(:,:,np1) = glm_past(:,:,np1)+(avail_landA-avail_land0)*fnfforest/cellarea
                      end where
                      final_area_needed =  past_pos_f - sumavail_land0-cumsum_sorted_farea(sortind(1)-1)
                   else
                      final_area_needed =  past_pos_f - sumavail_land0
                   end if
                   glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                        glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) + &
                        min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))), &
                        (avail_landA(v1u(sortind(1)),v2u(sortind(1))) - &
                        avail_landA(v1u(sortind(1)),v2u(sortind(1)))) * &
                        fnfforest(v1u(sortind(1)),v2u(sortind(1))) / &
                        cellarea(v1u(sortind(1)),v2u(sortind(1))))
                else
                   write(6,*)'pasture increase on forest - land not available'
                   call abort
                end if
             end if
             
          end if
       end do ! end n loop
    end do ! end n loop

    if (regional_unmet_reassign==1) then
       do r = 1,nreg
          regional_farea_needed = sum(unmet_farea(((r-1)*18+1):r*18))
          do aez = 1,naez
             raezs=0
             where ( aez_regions .EQ. r .and. aez_zones .EQ. aez)
                raezs=1
                rAEZ_sites=1.
             elsewhere
                rAEZ_sites=0.
             end where
             totraezs=sum(raezs)
             
             
             cellarea_forest(:,:)=0.
             cellarea_nonforest(:,:)=0.
             fnfnonforest=0.
             fnfforest=0.
             where ( pot_veg > 1.)
                cellarea_forest(:,:)=cellarea(:,:)
                fnfforest(:,:)=1.
             elsewhere
                fnfnonforest(:,:)=1.
                cellarea_nonforest(:,:)=cellarea(:,:)
             end where
             cellarea_nonforest=cellarea_nonforest*rAEZ_sites
             cellarea_forest=cellarea_forest*rAEZ_sites
             fnfforest=fnfforest*rAEZ_sites
             fnfnonforest=fnfnonforest*rAEZ_sites
             
             avail_farea(aez) = sum((pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest)
             avail_nfarea(aez) = sum((pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest)
             avail_ag_farea(aez) = sum((glm_crop(:,:,np1)+glm_past(:,:,np1))*cellarea_forest)
             reassign_ag(aez) = min(avail_ag_farea(aez), avail_nfarea(aez), regional_farea_needed)
             unmet_aez_farea(aez) = regional_farea_needed - reassign_ag(aez)
          end do
          
          !jt [sorted_reassign_ag,sort_aez] = sort(reassign_ag,'descend')
          indxa=(/(i,i=1,naez)/)
          indxadn=(/(i,i=1,naez)/)
          call D_mrgrnk(reassign_ag,indxa,naez)
          call D_mrgrnk(reassign_ag*-1.,indxadn,naez)
          cumsum_sorted_reassign_ag(1)=reassign_ag(indxadn(1))
          do i=2,naez
             cumsum_sorted_reassign_ag(i)=cumsum_sorted_reassign_ag(i-1)+reassign_ag(indxadn(i))
          end do
          sortind = MINLOC(cumsum_sorted_reassign_ag(:),mask=cumsum_sorted_reassign_ag(:) >= regional_farea_needed)
          if (sortind(1).eq.0) then
             sortind(1) = naez
             reassign_ag_at_max_aez_ind = reassign_ag(indxadn(sortind(1)))
             unmet_regional_farea(r) = regional_farea_needed - sum(reassign_ag(indxadn))
          elseif (sortind(1)>1) then 
             reassign_ag_at_max_aez_ind = regional_farea_needed - cumsum_sorted_reassign_ag(sortind(1)-1)
             unmet_regional_farea(r) = 0
          else
             reassign_ag_at_max_aez_ind = regional_farea_needed
             unmet_regional_farea(r) = 0
          end if
          
          do zz=1,sortind(1)
              z=indxadn(zz)
             if (reassign_ag(z)>0) then
                crop_before_decrease = sum(glm_crop(:,:,np1)*cellarea)
                past_before_decrease = sum(glm_past(:,:,np1)*cellarea)
                raezs=0
                where ( aez_regions .EQ. r .and. aez_zones .EQ. z)
                   raezs=1
                   rAEZ_sites=1.
                elsewhere
                   rAEZ_sites=0.
                end where
                totraezs=sum(raezs)
                
                
                cellarea_forest(:,:)=0.
                cellarea_nonforest(:,:)=0.
                fnfnonforest=0.
                fnfforest=0.
                where ( pot_veg > 1.)
                   cellarea_forest(:,:)=cellarea(:,:)
                   fnfforest(:,:)=1.
                elsewhere
                   fnfnonforest(:,:)=1.
                   cellarea_nonforest(:,:)=cellarea(:,:)
                end where
                cellarea_nonforest=cellarea_nonforest*rAEZ_sites
                cellarea_forest=cellarea_forest*rAEZ_sites
                fnfforest=fnfforest*rAEZ_sites
                fnfnonforest=fnfnonforest*rAEZ_sites
                indxup(:)=(/(i,i=1,numlons*numlats)/)
                indxdn=indxup
                pot_veg_rev=transpose(pot_veg)
                rAEZ_sites_rev=transpose(rAEZ_sites)
                call D_mrgrnk(pot_veg_rev*rAEZ_sites_rev, indxup,numlons*numlats)
                call D_mrgrnk(pot_veg_rev*rAEZ_sites_rev*-1., indxdn,numlons*numlats)
                
                !jt  The sortxxxup and sortxxxdn arrays are only good 1:totraezs these arrays are also based on row major
                !jt  sorting to match original matlab scripts.  We can get rid of this after validation.
                
                sortlonsdn=(indxdn-1)/numlats+1
                sortlatsdn=mod(indxdn-1,numlats)+1
                sortlonsup=(indxup-1)/numlats+1
                sortlatsup=mod(indxup-1,numlats)+1
                
                if (allocated(v1u)) deallocate(v1u)
                if (allocated(v2u)) deallocate(v2u)
                if (allocated(v1d)) deallocate(v1d)
                if (allocated(v2d)) deallocate(v2d)
                allocate(v1u(totraezs))
                allocate(v2u(totraezs))
                allocate(v1d(totraezs))
                allocate(v2d(totraezs))
                v1u=sortlonsup(numlons*numlats-totraezs+1:numlons*numlats)
                v2u=sortlatsup(numlons*numlats-totraezs+1:numlons*numlats)
                v1d=(sortlonsdn(:totraezs))
                v2d=(sortlatsdn(:totraezs))
                
                !jt [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'descend')
                cumsum_sorted_farea=0.
                call cumsum((glm_crop(:,:,np1)+glm_past(:,:,np1))*cellarea_forest,v1d,v2d,cumsum_sorted_farea(:totraezs),totraezs)
                sortind = MINLOC(cumsum_sorted_farea(:totraezs),cumsum_sorted_farea(:totraezs) > reassign_ag(z))
                if (sortind(1)==0) then
                      sortind(1) = totraezs
                end if
                if (sortind(1)>1) then
                   sortsitesdn=0
                   do i=1,sortind(1)-1
                      sortsitesdn(v1d(i),v2d(i))=1
                   end do
                   where(sortsitesdn>0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) - glm_crop(:,:,np1) * fnfforest
                   end where

                   where(sortsitesdn>0)
                      glm_past(:,:,np1) = glm_past(:,:,np1) - glm_past(:,:,np1) * fnfforest
                   end where

                   final_area_needed =  reassign_ag(z) - cumsum_sorted_farea(sortind(1)-1)
                else
                   final_area_needed =  reassign_ag(z)
                end if
                ag_area_avail = (glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) + glm_past(v1d(sortind(1)),v2d(sortind(1)),np1))*cellarea(v1d(sortind(1)),v2d(sortind(1)))
                crop_area_avail = glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1)*cellarea(v1d(sortind(1)),v2d(sortind(1)))
                glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) = &
                     glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) - &
                     min(crop_area_avail/(ag_area_avail+1e-12)*final_area_needed / &
                     cellarea(v1d(sortind(1)),v2d(sortind(1))), &
                     crop_area_avail/cellarea(v1d(sortind(1)),v2d(sortind(1))))

                glm_past(v1d(sortind(1)),v2d(sortind(1)),np1) = &
                     glm_past(v1d(sortind(1)),v2d(sortind(1)),np1) - &
                     (final_area_needed/cellarea(v1d(sortind(1)),v2d(sortind(1))) - &
                     min(crop_area_avail/(ag_area_avail+1e-12)*final_area_needed / &
                     cellarea(v1d(sortind(1)),v2d(sortind(1))),crop_area_avail / &
                     cellarea(v1d(sortind(1)),v2d(sortind(1)))))
                
                crop_after_decrease = sum(glm_crop(:,:,np1)*cellarea)
                past_after_decrease = sum(glm_past(:,:,np1)*cellarea)
                total_ag_decrease = crop_before_decrease+past_before_decrease - crop_after_decrease - past_after_decrease
                crop_decrease_ratio = (crop_before_decrease - crop_after_decrease)/(total_ag_decrease+1e-12)
                past_decrease_ratio = 1-crop_decrease_ratio
                
                avail_land0 = 0.
                
                where ( aez_regions .EQ. r .and. aez_zones .EQ. z.and.(glm_crop(:,:,np1)>0..or.glm_past(:,:,np1)>0.))
                   avail_land0=(pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
                end where
                sumavail_land0=sum(avail_land0)
                if (sumavail_land0 >=reassign_ag(z)) then 
#ifdef DEBUG
                   write(6,*)'crop and pasture increase on nonforest - land available'
#endif
                   where(raez_sites > 0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + crop_decrease_ratio*avail_land0 / &
                           (sumavail_land0+1e-12) * reassign_ag(z)*fnfnonforest/cellarea
                      glm_past(:,:,np1) = glm_past(:,:,np1) + past_decrease_ratio*avail_land0 / &
                           (sumavail_land0+1e-12) * reassign_ag(z)*fnfnonforest/cellarea
                   end where
                else
                   avail_landA = 0.
                   where ( aez_regions .EQ. r .and. aez_zones .EQ. z)
                      avail_landA=(pctland_in2005-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
                   end where
                   sumavail_landA=sum(avail_landA)
                   if (reassign_ag(z) > sumavail_landA) then 
                      reassign_ag(z) = sum(avail_landA)
                      call abort
                   end if
#ifdef DEBUG
                   write(6,*)'cropland and pasture increase on non-forest - land available'
#endif
                   where(raez_sites > 0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + crop_decrease_ratio*(avail_land0*fnfnonforest)/cellarea
                      glm_past(:,:,np1) = glm_past(:,:,np1) + past_decrease_ratio*(avail_land0*fnfnonforest)/cellarea
                   end where
                      
                   where(raez_sites > 0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + crop_decrease_ratio*((avail_landA-avail_land0) / &
                           (sumavail_landA-sumavail_land0+1e-12)*(reassign_ag(z)-sumavail_land0)*fnfnonforest)/cellarea
                      glm_past(:,:,np1) = glm_past(:,:,np1) + past_decrease_ratio*((avail_landA-avail_land0)/ &
                           (sumavail_landA-sumavail_land0+1e-12)*(reassign_ag(z)-sumavail_land0)*fnfnonforest)/cellarea
                   end where
                end if
             end if
          end do ! end z loop
       end do !  end r loop
    end if   ! ! if regional_unmet_reassign
 end if

#ifdef DEBUG
    write(6,*) 'sum final gcrop/gpast n'
    write(6,fmt="(1ES25.15)") sum(glm_crop(:,:,n))
    write(6,fmt="(1ES25.15)") sum(glm_past(:,:,n))
    write(6,*) 'sum final gcrop/gpast np1'
    write(6,fmt="(1ES25.15)") sum(glm_crop(:,:,np1))
    write(6,fmt="(1ES25.15)") sum(glm_past(:,:,np1))
#endif 
!
! Interpolate data to years needed by GLM
!
! glm calculates year+1 using 
! gcam constructed states at year+1
! and harvest transitions at year
!
! since eclock is already advanced 
! by 1 year for gcam2glm and glm we have 
!
! glm calculates eclock year using 
! gcam constructed states at eclock year
! and harvest transitions at eclock year minus 1

    eclockyr=ymd/10000
    eclockyrm1=eclockyr-1
    delyr= year2-year1
    fact1yrm1=(year2-eclockyrm1)/delyr
    fact2yrm1=(eclockyrm1-year1)/delyr
    fact1=(year2-eclockyr)/delyr
    fact2=(eclockyr-year1)/delyr

#ifdef DEBUG
     write(iu,*)'crop interpolation factors fact1,fact2,year1,year2=',fact1,fact2,year1,year2
#endif
! use eclock year year for interpolating crop past and othr
    glm_crop_ann(:,:)=glm_crop(:,:,n)*fact1+glm_crop(:,:,np1)*fact2
    glm_past_ann(:,:)=glm_past(:,:,n)*fact1+glm_past(:,:,np1)*fact2
    glm_othr_ann(:,:)=pctland_in2005-glm_past_ann-glm_crop_ann

! use previous year for fractions fact1 and fact2 for woodharvest
    glm_wh_ann(:)=gcam_wh(:,n)*fact1yrm1+gcam_wh(:,np1)*fact2yrm1

!    do j=1,360
!       do i=1,720
!          glm_crop_ann(i,j)=fround(glm_crop_ann(i,j),6)
!          glm_past_ann(i,j)=fround(glm_past_ann(i,j),6)
!          glm_othr_ann(i,j)=fround(glm_othr_ann(i,j),6)
!       end do
!    end do


!  Conversion factor of 1.3 to account for a slash fraction 
!  of 30% for non-harvested carbon lost as a result of wood harvest

    glm_wh_ann(:)=glm_wh_ann(:)*1.3

    write(yearc,fmt="(I4)") ymd/10000
    open (unit=55,file='gcrop'//yearc//'.txt',action="write",status="unknown")
    do j=1,360
       write(55,fmt="(720F9.6)") glm_crop_ann(:,j)
    end do
    close(55)

    open (unit=55,file='gpast'//yearc//'.txt',action="write",status="unknown")
    do j=1,360
       write(55,fmt="(720F9.6)") glm_past_ann(:,j)
    end do
    close(55)

    open (unit=55,file='gothr'//yearc//'.txt',action="write",status="unknown")
    do j=1,360
       write(55,fmt="(720F9.6)") glm_othr_ann(:,j)
    end do
    close(55)

    ij=0
    do j=1,numLats
       do i=1,numLons
          ij=ij+1
          glmi(iac_glmi_cropland,ij)=glm_crop_ann(i,j)
          glmi(iac_glmi_pasture,ij)=glm_past_ann(i,j)
          glmi(iac_glmi_natveg,ij)=glm_othr_ann(i,j)
       end do
    end do
    glmi_wh(:)=glm_wh_ann(:)

!
! If we are at the boundary year of the gcam data. Need to move np1 info
! into timelevel n to calculate the next set of years.
!
    if (eclockyr==year2) then 
       glm_crop(:,:,n)=glm_crop(:,:,np1)
       glm_past(:,:,n)=glm_past(:,:,np1)
    end if

    ! lets write a restart every time this routine is called

    call shr_cal_date2ymd(ymd,year,mon,day)
    write(filename,'(a,i4.4,a,i2.2,a)') trim(gcam2glm_restfile)//'r.',year,'.nc'

    iun = shr_file_getunit()
    open(iun,file=trim(gcam2glm_rpointer),form='formatted')
    write(iun,'(a)') trim(filename)
    close(iun)
    call shr_file_freeunit(iun)

    write(iu,*) subname,' write_restart rpointer ',trim(gcam2glm_rpointer)
    write(iu,*) subname,' write_restart file     ',trim(filename)

    status= nf90_create(filename,nf90_clobber,ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_dim(ncid,'lon',numlons,dimid3(1))
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_def_dim(ncid,'lat',numlats,dimid3(2))
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_def_dim(ncid,'time',2,dimid3(3))
    if(status /= nf90_NoErr) call handle_err(status)

    dimid = dimid3(1)
    status = nf90_def_var(ncid,'lon',NF90_DOUBLE,dimid,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"units","degrees_east")
    if(status /= nf90_NoErr) call handle_err(status)

    dimid = dimid3(2)
    status = nf90_def_var(ncid,'lat',NF90_DOUBLE,dimid,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"units","degrees_north")
    if(status /= nf90_NoErr) call handle_err(status)

    dimid = dimid3(3)
    status = nf90_def_var(ncid,'gcam_years',NF90_INT,dimid,varid)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_var(ncid,'glm_crop',NF90_DOUBLE,dimid3,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"fraction","frac")
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"missing_value",miss_val)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_var(ncid,'glm_past',NF90_DOUBLE,dimid3,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"fraction","frac")
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"missing_value",miss_val)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_enddef(ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'lon',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,lon)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'lat',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,lat)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'gcam_years',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,(/cdata%i(iac_cdatai_gcam_yr1),cdata%i(iac_cdatai_gcam_yr2)/))
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'glm_crop',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,glm_crop)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'glm_past',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,glm_past)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_close(ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    deallocate(indxup)
    deallocate(indxdn)
    deallocate(sortlatsup)
    deallocate(sortlatsdn)
    deallocate(sortlonsup)
    deallocate(sortlonsdn)
    deallocate(tmparr)
    deallocate(indxa)
    deallocate(indxadn)
    deallocate(avail_farea)
    deallocate(avail_nfarea)
    deallocate(avail_ag_farea)
    deallocate(reassign_ag)
    deallocate(unmet_aez_farea)
    deallocate(cumsum_sorted_reassign_ag)
    deallocate(unmet_regional_farea)

  end subroutine gcam2glm_run_mod
  
  
!---------------------------------------------------------------------------
!BOP
  
! !IROUTINE: gcam2glm_final_mod

! !INTERFACE:
  subroutine gcam2glm_final_mod( )

! !DESCRIPTION:
! Finalize glm model
! !USES:
    implicit none

! !ARGUMENTS:

! !LOCAL VARIABLES:
    integer :: iu
    character(len=*),parameter :: subname='(gcam2glm_final_mod)'

! !REVISION HISTORY:
! Author: T Craig

!EOP

!---------------------------------------------------------------------------

!    iu  = nint(cdata(iac_cdata_logunit))
!    write(iu,*) trim(subname)

    deallocate(hydeGCROP2005)
    deallocate(hydeGPAST2005)
    deallocate(hydeGOTHR2005)
    deallocate(hydeGWH2005)
    deallocate(cellarea)
    deallocate(cellarea_forest)
    deallocate(cellarea_nonforest)
    deallocate(glm_crop_ann)
    deallocate(glm_past_ann)
    deallocate(glm_othr_ann)
    deallocate(cumsum_sorted_farea)
    deallocate(glm_wh_ann)
    deallocate(aez_regions)
    deallocate(aez_zones)
    deallocate(fnfforest)
    deallocate(fnfnonforest)
    deallocate(pot_veg)
    deallocate(pot_veg_rev)
    deallocate(crop_area)
    deallocate(pctland_in2005)
    deallocate(sortsitesup)
    deallocate(sortsitesdn)
    deallocate(rAEZ_sites)
    deallocate(rAEZ_sites_rev)
    deallocate(raezs)
    deallocate(datearr)
    deallocate(glm_crop)
    deallocate(glm_past)
    deallocate(gcam_crop)
    deallocate(gcam_wh)
    deallocate(gcam_past)
    deallocate(gcam_forest_area)
    deallocate(unmet_neg_past)
    deallocate(unmet_neg_crop)
    deallocate(unmet_farea)
    deallocate(avail_land0)
    deallocate(avail_landA)
    deallocate(lon)
    deallocate(lat)

  end subroutine gcam2glm_final_mod
!====================================================================================
real(r8) function fround(n, d)
real(r8) n
integer d
  fround= floor(n * 10.**d + .5) / 10.**d
end function fround
!====================================================================================
 subroutine handle_err (status)
    implicit none
    integer, intent (in) :: status
    print *, nf90_strerror(status)
    stop
 end subroutine handle_err
 subroutine cumsum (arrayin,sortcol,sortrow,cumsumvecout,n)
   implicit none
   real(r8), intent (in) :: arrayin(:,:)
   integer, intent (in) :: sortcol(:),sortrow(:)
   integer, intent (in) :: n
   real(r8), intent (out):: cumsumvecout(:)
   integer             ::j
   cumsumvecout(1)=arrayin(sortcol(1),sortrow(1))
   do j=2,n
      cumsumvecout(j)=cumsumvecout(j-1)+arrayin(sortcol(j),sortrow(j))
   end do
 end subroutine cumsum
!====================================================================================
Subroutine D_mrgrnk (XDONT, IRNGT,N)
!                                                                                                                                                      
!   MRGRNK = Merge-sort ranking of an array                                                                                                            
!   For performance reasons, the first 2 passes are taken                                                                                              
!   out of the standard loop, and use dedicated coding.                                                                                                
!                                                                                                                                                      
!   The routine is part of ORDERPACK 2.0 -- Unconditional, Unique, and Partial Ranking,                                                                
!   Sorting, and Permutation Downloadable Fortran 90 source code                                                                                       
!   Author: Michel Olagnon                                                                                                                             
!   http://www.fortran-2000.com/rank/                                                                                                                  
!   Users can freely download ORDERPACK 2.0 from this site.                                                                                            
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (N), Intent (In) :: XDONT
      Integer, Dimension (N), Intent (Out) :: IRNGT
      Integer, Intent (In) :: N
! __________________________________________________________
      Real (kind=kdp) :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine D_mrgrnk


!====================================================================================
end module gcam2glm_mod

